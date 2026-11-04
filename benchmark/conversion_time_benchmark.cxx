#include <benchmark/benchmark.h>
#include "benchmark_config.h"
#include "benchmark_utils.h"
#include "generate_sam_benchmark.h"
#include "ramcore/SamToNTuple.h"
#include <filesystem>
#include <iostream>
#include <cstdio>

static void BM_SamToRNTuple(benchmark::State &state, const BenchmarkOptions &opts)
{
   int num_reads = state.range(0);
   const int threads = static_cast<int>(state.range(1));
   std::string sam_file = "rntuple_test.sam";
   std::string rntuple_file = "rntuple_output.root";

   GenerateSAMFile(sam_file, num_reads);

   for (auto _ : state) {
      bool converted = false;
      {
         const StdoutSilencer silence;
         converted =
            samtoramntuple(sam_file.c_str(), rntuple_file.c_str(), opts.compression, opts.quality_policy, threads);
      }
      if (!converted) {
         state.SkipWithError("samtoramntuple failed");
         break;
      }

      if (std::filesystem::exists(rntuple_file)) {
         auto file_size = std::filesystem::file_size(rntuple_file);
         state.counters["file_size_mb"] = file_size / (1024.0 * 1024.0);
      }

      std::remove(rntuple_file.c_str());
   }

   std::remove(sam_file.c_str());
   state.counters["reads_per_second"] = benchmark::Counter(num_reads, benchmark::Counter::kIsRate);
}

int main(int argc, char **argv)
{
   BenchmarkOptions opts;
   opts.reads = {1000, 10000, 100000};
   opts.threads = {1};
   if (!ParseBenchmarkFlags(&argc, argv, kFlagCompression | kFlagQuality | kFlagReads | kFlagThreads, opts))
      return 1;

   std::cout << "Conversion Time Benchmark" << std::endl;
   std::cout << std::endl;

   auto *bench = ::benchmark::RegisterBenchmark("RNTuple_Conversion", BM_SamToRNTuple, opts);
   bench->ArgNames({"reads", "threads"});
   for (const int reads : opts.reads) {
      for (const int threads : opts.threads)
         bench->Args({reads, threads});
   }
   // CPU time counts the main thread only, which with -threads mostly waits:
   // rates and iteration counts go by wall time instead.
   bench->UseRealTime()->Unit(benchmark::kMillisecond);

   ::benchmark::Initialize(&argc, argv);
   if (::benchmark::ReportUnrecognizedArguments(argc, argv))
      return 1;
   ::benchmark::RunSpecifiedBenchmarks();
   ::benchmark::Shutdown();
   return 0;
}
