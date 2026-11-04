#include <benchmark/benchmark.h>
#include "benchmark_config.h"
#include "benchmark_utils.h"
#include "generate_sam_benchmark.h"
#include "ramcore/SamToNTuple.h"
#include <filesystem>
#include <iostream>
#include <cstdio>

static void BM_SamToRam(benchmark::State &state, const BenchmarkOptions &opts)
{
   int num_reads = state.range(0);

   std::string sam_file = "benchmark_test.sam";
   std::string rntuple_file = "benchmark_rntuple.root";

   GenerateSAMFile(sam_file, num_reads);

   for (auto _ : state) {
      bool converted = false;
      {
         const StdoutSilencer silence;
         converted = samtoramntuple(sam_file.c_str(), rntuple_file.c_str(), opts.compression, opts.quality_policy);
      }
      if (!converted) {
         state.SkipWithError("samtoramntuple failed");
         break;
      }

      if (std::filesystem::exists(rntuple_file)) {
         auto rntuple_size = std::filesystem::file_size(rntuple_file);
         state.counters["rntuple_size_mb"] = rntuple_size / (1024.0 * 1024.0);
      }
   }

   std::remove(sam_file.c_str());
   std::remove(rntuple_file.c_str());

   state.counters["reads"] = num_reads;
   state.counters["reads_per_second"] = benchmark::Counter(num_reads, benchmark::Counter::kIsRate);
}

int main(int argc, char **argv)
{
   BenchmarkOptions opts;
   opts.reads = {1000, 10000, 100000};
   if (!ParseBenchmarkFlags(&argc, argv, kFlagCompression | kFlagQuality | kFlagReads, opts))
      return 1;

   std::cout << "SAM to RAM Conversion Benchmark" << std::endl;
   std::cout << std::endl;

   auto *bench = ::benchmark::RegisterBenchmark("BM_SamToRam", BM_SamToRam, opts);
   for (const int reads : opts.reads)
      bench->Arg(reads);
   bench->Unit(benchmark::kMillisecond);

   ::benchmark::Initialize(&argc, argv);
   if (::benchmark::ReportUnrecognizedArguments(argc, argv))
      return 1;
   ::benchmark::RunSpecifiedBenchmarks();
   ::benchmark::Shutdown();
   return 0;
}
