#include <benchmark/benchmark.h>
#include "benchmark_config.h"
#include "benchmark_utils.h"
#include "ramcore/RAMNTupleView.h"
#include <string>
#include <cstddef>
#include <cstdint>

// The file and the regions come from flags, or from the environment as before,
// so the numbers are reproducible on any machine. A hard-coded path is not a
// benchmark result.
//
//   --ram=FILE      or RAMTOOLS_BENCH_RNTUPLE   RAM file written by samtoramntuple (required)
//   --regions=R,... or RAMTOOLS_BENCH_REGIONS   regions to query (optional)

static void BM_RegionQuery(benchmark::State &state, const std::string &ram, const std::string &region)
{
   if (ram.empty()) {
      state.SkipWithError("pass --ram=FILE, or set RAMTOOLS_BENCH_RNTUPLE, to a file written by samtoramntuple");
      return;
   }

   int64_t total_reads_processed = 0;
   std::int64_t reads_in_this_run = 0;
   const RAMNTupleViewOpts opts{};

   for (auto _ : state) {
      {
         const StdoutSilencer silence;
         reads_in_this_run = ramntupleview(ram.c_str(), region.c_str(), opts);
      }

      total_reads_processed += reads_in_this_run;
   }

   state.SetItemsProcessed(total_reads_processed);
   state.SetLabel(region + ": " + std::to_string(reads_in_this_run) + " reads");
}

int main(int argc, char **argv)
{
   BenchmarkOptions opts;
   opts.regions = {"chr1:1000000-1001000", "chr1:1000000-2000000", "chr1:1-50000000", "chr21:1-48129895"};
   if (!ParseBenchmarkFlags(&argc, argv, kFlagRam | kFlagRegions, opts))
      return 1;

   // One case per region, numbered in the order given.
   for (std::size_t i = 0; i < opts.regions.size(); ++i) {
      ::benchmark::RegisterBenchmark("RegionQueryFixture/RNTuple/" + std::to_string(i), BM_RegionQuery, opts.ram,
                                     opts.regions[i])
         ->Unit(benchmark::kSecond);
   }

   ::benchmark::Initialize(&argc, argv);
   if (::benchmark::ReportUnrecognizedArguments(argc, argv))
      return 1;
   ::benchmark::RunSpecifiedBenchmarks();
   ::benchmark::Shutdown();
   return 0;
}
