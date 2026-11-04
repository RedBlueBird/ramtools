#include <benchmark/benchmark.h>
#include "bam_ops.h"
#include "benchmark_config.h"
#include "benchmark_utils.h"
#include "ramcore/SamToNTuple.h"
#include <TROOT.h>
#include <cstdint>
#include <string>

// Both sides start from the SAM and end with one file per reference. BAM needs
// a sorted, indexed BAM to extract references from, so its time includes making
// one, as the samtools pipeline did (view, sort, index, view per reference).

static double MiB(uint64_t bytes)
{
   return static_cast<double>(bytes) / (1024.0 * 1024.0);
}

static void SetRates(benchmark::State &state, const SamInput &input)
{
   const int64_t records = input.reads > 0 ? input.reads : CountSamRecords(input.path);
   state.counters["reads/s"] = benchmark::Counter(static_cast<double>(records), benchmark::Counter::kIsRate);
}

// ROOT's implicit threads compress the split's pages, as samtoramntuple -split
// -threads N enables them; 1 disables them.
static void SetImplicitMT(int threads)
{
   if (ROOT::IsImplicitMTEnabled())
      ROOT::DisableImplicitMT();
   if (threads > 1)
      ROOT::EnableImplicitMT(threads);
}

static void BM_BamSplit(benchmark::State &state, const SamInput &input, const std::string &dir)
{
   PrepareSamInput(input);
   const int threads = static_cast<int>(state.range(0));
   const std::string sorted_bam = dir + "/bam_sorted.bam";
   const std::string prefix = dir + "/bam_split_";

   for ([[maybe_unused]] const auto &_ : state) {
      if (!SamToSortedBam(input.path, sorted_bam, threads) || !SplitBam(sorted_bam, prefix, threads)) {
         state.SkipWithError("the BAM split failed on " + input.path);
         break;
      }

      state.counters["size_MB"] = MiB(TotalSizeWithPrefix(dir, "bam_split_"));
      RemoveWithPrefix(dir, "bam_split_");
      RemoveWithPrefix(dir, "bam_sorted.bam");
   }
   SetRates(state, input);
}

static void
BM_ChromosomeSplit(benchmark::State &state, const SamInput &input, const std::string &dir, const BenchmarkOptions &opts)
{
   PrepareSamInput(input);
   const std::string prefix = dir + "/ram_split";
   SetImplicitMT(static_cast<int>(state.range(0)));

   for ([[maybe_unused]] const auto &_ : state) {
      {
         const StdoutSilencer silence;
         samtoramntuple_split_by_chromosome(input.path.c_str(), prefix.c_str(), opts.compression, opts.quality_policy);
      }

      state.counters["size_MB"] = MiB(TotalSizeWithPrefix(dir, "ram_split_"));
      RemoveWithPrefix(dir, "ram_split_");
   }
   SetImplicitMT(/*threads=*/1);
   SetRates(state, input);
}

int main(int argc, char **argv)
{
   BenchmarkOptions opts;
   opts.reads = {100000, 500000, 1000000};
   opts.threads = {1, 2, 4};
   if (!ParseBenchmarkFlags(
          &argc, argv, kFlagCompression | kFlagQuality | kFlagReads | kFlagThreads | kFlagSam | kFlagWorkDir, opts))
      return 1;
   // Before any input is prepared, so that a mistyped flag or --help returns at once.
   ::benchmark::Initialize(&argc, argv);
   if (::benchmark::ReportUnrecognizedArguments(argc, argv))
      return 1;
   const WorkDir dir(opts.work_dir);

   // Synthetic inputs are sorted, as aligned reads are before they are split;
   // the BAM side sorts an unsorted --sam itself. The work happens on threads
   // the CPU time of this one does not see: every rate goes by wall time.
   for (const SamInput &input : SamInputs(opts, dir, /*sorted=*/true)) {
      auto *bam = ::benchmark::RegisterBenchmark("BM_BamSplit/" + input.name, BM_BamSplit, input, dir.Path());
      auto *ram = ::benchmark::RegisterBenchmark("BM_ChromosomeSplit/" + input.name, BM_ChromosomeSplit, input,
                                                 dir.Path(), opts);
      for (auto *bench : {bam, ram}) {
         bench->ArgName("threads");
         for (const int threads : opts.threads)
            bench->Arg(threads);
         bench->UseRealTime()->Unit(benchmark::kMillisecond);
      }
   }

   ::benchmark::RunSpecifiedBenchmarks();
   ::benchmark::Shutdown();
   return 0;
}
