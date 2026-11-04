#include <benchmark/benchmark.h>
#include "benchmark_config.h"
#include "benchmark_utils.h"
#include "ramcore/SamToNTuple.h"
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>
#include <system_error>

static void
BM_SamToRNTuple(benchmark::State &state, const SamInput &input, const std::string &output, const BenchmarkOptions &opts)
{
   const int threads = static_cast<int>(state.range(0));
   PrepareSamInput(input);
   const int64_t records = input.reads > 0 ? input.reads : CountSamRecords(input.path);

   for ([[maybe_unused]] const auto &_ : state) {
      bool converted = false;
      {
         const StdoutSilencer silence;
         converted = samtoramntuple(input.path.c_str(), output.c_str(), opts.compression, opts.quality_policy, threads);
      }
      if (!converted) {
         state.SkipWithError("samtoramntuple failed on " + input.path);
         break;
      }
   }

   std::error_code ec;
   const auto ram_size = std::filesystem::file_size(output, ec);
   if (!ec)
      state.counters["file_size_mb"] = static_cast<double>(ram_size) / (1024.0 * 1024.0);
   std::remove(output.c_str());

   const auto sam_size = std::filesystem::file_size(input.path, ec);
   if (!ec)
      state.SetBytesProcessed(state.iterations() * static_cast<int64_t>(sam_size));
   state.counters["reads"] = static_cast<double>(records);
   state.counters["reads_per_second"] = benchmark::Counter(static_cast<double>(records), benchmark::Counter::kIsRate);
}

int main(int argc, char **argv)
{
   BenchmarkOptions opts;
   opts.reads = {1000, 10000, 100000};
   opts.threads = {1};
   if (!ParseBenchmarkFlags(
          &argc, argv, kFlagCompression | kFlagQuality | kFlagReads | kFlagThreads | kFlagSam | kFlagWorkDir, opts))
      return 1;
   // Before any input is prepared, so that a mistyped flag or --help returns at once.
   ::benchmark::Initialize(&argc, argv);
   if (::benchmark::ReportUnrecognizedArguments(argc, argv))
      return 1;
   const WorkDir dir(opts.work_dir);

   std::cout << "Conversion Time Benchmark" << std::endl;
   std::cout << std::endl;

   const std::string output = dir.File("conversion_output.ram");
   for (const SamInput &input : SamInputs(opts, dir, /*sorted=*/false)) {
      auto *bench =
         ::benchmark::RegisterBenchmark("RNTuple_Conversion/" + input.name, BM_SamToRNTuple, input, output, opts);
      bench->ArgName("threads");
      for (const int threads : opts.threads)
         bench->Arg(threads);
      // CPU time counts the main thread only, which with -threads mostly waits:
      // rates and iteration counts go by wall time instead.
      bench->UseRealTime()->Unit(benchmark::kMillisecond);
   }

   ::benchmark::RunSpecifiedBenchmarks();
   ::benchmark::Shutdown();
   return 0;
}
