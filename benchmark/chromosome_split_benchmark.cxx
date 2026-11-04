#include <benchmark/benchmark.h>
#include "benchmark_config.h"
#include "benchmark_utils.h"
#include "ramcore/SamToNTuple.h"
#include "generate_sam_benchmark.h"
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <vector>
#include <string>
#include <thread>

static void CleanupFiles(const std::string &pattern)
{
   for (const auto &entry : std::filesystem::directory_iterator(".")) {
      if (entry.path().filename().string().find(pattern) != std::string::npos) {
         std::remove(entry.path().c_str());
      }
   }
}

static size_t GetTotalFileSize(const std::string &pattern)
{
   size_t total = 0;
   for (const auto &entry : std::filesystem::directory_iterator(".")) {
      if (entry.path().filename().string().find(pattern) != std::string::npos) {
         total += std::filesystem::file_size(entry.path());
      }
   }
   return total;
}

static std::vector<std::string> GetChromosomes(const std::string &sam_file)
{
   std::vector<std::string> chroms;
   std::ifstream f(sam_file);
   std::string line;

   while (std::getline(f, line) && line[0] == '@') {
      if (line.find("@SQ\tSN:") == 0) {
         size_t start = 7;
         size_t end = line.find('\t', start);
         chroms.push_back(line.substr(start, end - start));
      }
   }
   return chroms;
}

static void BM_SamtoolsSplit(benchmark::State &state)
{
   int num_reads = state.range(0);
   std::string sam_file = "bench_st_" + std::to_string(num_reads) + ".sam";

   GenerateSAMFile(sam_file, num_reads);
   auto chromosomes = GetChromosomes(sam_file);

   for (auto _ : state) {
      std::string bam_file = "bench_st_tmp.bam";
      std::string sorted_bam = "bench_st_sorted.bam";

      system(("samtools view -bS " + sam_file + " -o " + bam_file + " 2>/dev/null").c_str());

      system(("samtools sort " + bam_file + " -o " + sorted_bam + " 2>/dev/null").c_str());

      system(("samtools index " + sorted_bam + " 2>/dev/null").c_str());

      for (const auto &chr : chromosomes) {
         std::string cmd = "samtools view -b " + sorted_bam + " " + chr + " > bench_st_" + chr + ".bam 2>/dev/null";
         system(cmd.c_str());
      }

      state.counters["size_MB"] = GetTotalFileSize("bench_st_chr") / (1024.0 * 1024.0);

      CleanupFiles("bench_st_chr");
      std::remove(bam_file.c_str());
      std::remove(sorted_bam.c_str());
      std::remove((sorted_bam + ".bai").c_str());
   }

   std::remove(sam_file.c_str());
   state.counters["reads/s"] = benchmark::Counter(num_reads, benchmark::Counter::kIsRate);
}

static void BM_SamtoolsSplitThreaded(benchmark::State &state)
{
   int num_reads = state.range(0);
   int num_threads = state.range(1);
   std::string sam_file = "bench_st_mt_" + std::to_string(num_reads) + ".sam";

   GenerateSAMFile(sam_file, num_reads);
   auto chromosomes = GetChromosomes(sam_file);

   for (auto _ : state) {
      std::string bam_file = "bench_st_mt_tmp.bam";
      std::string sorted_bam = "bench_st_mt_sorted.bam";

      std::string cmd =
         "samtools view -@ " + std::to_string(num_threads) + " -bS " + sam_file + " -o " + bam_file + " 2>/dev/null";
      system(cmd.c_str());

      cmd = "samtools sort -@ " + std::to_string(num_threads) + " -m 1G " + bam_file + " -o " + sorted_bam +
            " 2>/dev/null";
      system(cmd.c_str());

      cmd = "samtools index -@ " + std::to_string(num_threads) + " " + sorted_bam + " 2>/dev/null";
      system(cmd.c_str());

      std::vector<std::thread> threads;
      for (const auto &chr : chromosomes) {
         threads.emplace_back([&sorted_bam, &chr]() {
            std::string cmd =
               "samtools view -@ 2 -b " + sorted_bam + " " + chr + " > bench_st_mt_" + chr + ".bam 2>/dev/null";
            system(cmd.c_str());
         });

         if (threads.size() >= num_threads) {
            for (auto &t : threads) {
               t.join();
            }
            threads.clear();
         }
      }

      for (auto &t : threads) {
         t.join();
      }

      state.counters["size_MB"] = GetTotalFileSize("bench_st_mt_chr") / (1024.0 * 1024.0);
      state.counters["threads"] = num_threads;

      CleanupFiles("bench_st_mt_chr");
      std::remove(bam_file.c_str());
      std::remove(sorted_bam.c_str());
      std::remove((sorted_bam + ".bai").c_str());
   }

   std::remove(sam_file.c_str());
   state.counters["reads/s"] = benchmark::Counter(num_reads, benchmark::Counter::kIsRate);
}

static void BM_ChromosomeSplit(benchmark::State &state, const BenchmarkOptions &opts)
{
   int num_reads = state.range(0);
   std::string sam_file = "bench_split_par_" + std::to_string(num_reads) + ".sam";

   GenerateSAMFile(sam_file, num_reads);

   for (auto _ : state) {
      {
         const StdoutSilencer silence;
         samtoramntuple_split_by_chromosome(sam_file.c_str(), /*output_prefix=*/"bench_split_par_out", opts.compression,
                                            opts.quality_policy);
      }

      state.counters["size_MB"] = GetTotalFileSize("bench_split_par_out_") / (1024.0 * 1024.0);
      CleanupFiles("bench_split_par_out_");
   }

   std::remove(sam_file.c_str());
   state.counters["reads/s"] = benchmark::Counter(num_reads, benchmark::Counter::kIsRate);
}

int main(int argc, char **argv)
{
   BenchmarkOptions opts;
   opts.reads = {100000, 500000, 1000000};
   opts.threads = {2, 4};
   if (!ParseBenchmarkFlags(&argc, argv, kFlagCompression | kFlagQuality | kFlagReads | kFlagThreads, opts))
      return 1;

   auto *samtools = ::benchmark::RegisterBenchmark("BM_SamtoolsSplit", BM_SamtoolsSplit);
   auto *samtools_threaded = ::benchmark::RegisterBenchmark("BM_SamtoolsSplitThreaded", BM_SamtoolsSplitThreaded);
   auto *split = ::benchmark::RegisterBenchmark("BM_ChromosomeSplit", BM_ChromosomeSplit, opts);
   for (const int reads : opts.reads) {
      samtools->Arg(reads);
      for (const int threads : opts.threads)
         samtools_threaded->Args({reads, threads});
      split->Arg(reads);
   }
   for (auto *bench : {samtools, samtools_threaded, split})
      bench->Unit(benchmark::kMillisecond);

   ::benchmark::Initialize(&argc, argv);
   if (::benchmark::ReportUnrecognizedArguments(argc, argv))
      return 1;
   ::benchmark::RunSpecifiedBenchmarks();
   ::benchmark::Shutdown();
   return 0;
}
