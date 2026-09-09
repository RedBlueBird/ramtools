#include "bam_ops.h"
#include "benchmark_config.h"
#include "benchmark_utils.h"
#include <RtypesCore.h>
#include <benchmark/benchmark.h>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

Long64_t ramview(const char *file, const char *query, bool cache = true, bool perfstats = false,
                 const char *perfstatsfilename = "perf.root");
Long64_t ramntupleview(const char *file, const char *query, bool cache = true, bool perfstats = false,
                       const char *perfstatsfilename = "perf.root");

// A scale ladder of regions, from single-position to whole-chromosome, matching
// GenerateSAMFile()'s synthetic "chrN" naming and hg19 lengths. Used when no
// --regions-file is given. Indices into this list (or the --regions-file contents)
// are what --regions selects (default 0,3,6,9 keeps the historical set).
static const std::vector<std::string> kSyntheticRegions = {"chr1:1000000-1001000",
                                                           "chr2:5000000-5010000",
                                                           "chrX:100000-150000",
                                                           "chr1:1000000-2000000",
                                                           "chr5:10000000-15000000",
                                                           "chr10:50000000-60000000",
                                                           "chr1:1-50000000",
                                                           "chr2:1-100000000",
                                                           "chr7:50000000-150000000",
                                                           "chr21:1-48129895",
                                                           "chrM:1-16571",
                                                           "chrY:2600000-2700000",
                                                           "GL000227.1:1-100000",
                                                           "chr1:1-1000",
                                                           "chr1:249250621-249250621",
                                                           "chr22:51304566-51304566",
                                                           "chr17:41196312-41277500",
                                                           "chr13:32889611-32973805"};

// Time one region query, whatever the backend. `query` returns the number of records
// matched, and is called once per iteration.
//
// Every backend opens its file inside the timed loop. That is not an oversight: ramview()
// and ramntupleview() are whole-tool entry points that reopen and re-read the index on
// each call, so the BAM side opens its samFile and loads the .bai per iteration too.
// Hoisting only the BAM open would make its column a steady-state latency while the ROOT
// columns stayed cold-open latency, and the three would no longer be the same measurement.
template <typename QueryFn>
static void RunQuery(benchmark::State &state, QueryFn &&query, int region_idx)
{
   int64_t total_reads_processed = 0;
   int64_t reads_in_this_run = 0;

   {
      // Suppress the per-call chatter (stopwatch + "Found N records") for the whole run;
      // constructed outside the timed loop so it adds no measured overhead.
      benchutil::ScopedStdoutSuppressor quiet;
      for ([[maybe_unused]] auto _ : state) {
         reads_in_this_run = query();
         total_reads_processed += reads_in_this_run;
      }
   }

   if (reads_in_this_run < 0) {
      state.SkipWithError("query failed (missing file or index?)");
      return;
   }
   state.SetItemsProcessed(total_reads_processed);
   state.counters["region_idx"] = region_idx;
   state.SetLabel(std::to_string(reads_in_this_run) + " reads");
}

// Reads "name:start-end" regions from `path`, one per line. Blank lines and lines
// starting with '#' are ignored. Use this to query a real dataset, since
// kSyntheticRegions assumes GenerateSAMFile()'s "chrN" naming and hg19 lengths.
static std::vector<std::string> LoadRegionsFromFile(const std::string &path)
{
   std::vector<std::string> regions;
   std::ifstream in(path);
   if (!in) {
      std::cerr << "region_query_benchmark: could not open --regions-file=" << path << "\n";
      return regions;
   }
   std::string line;
   while (std::getline(in, line)) {
      const auto first = line.find_first_not_of(" \t\r\n");
      if (first == std::string::npos || line[first] == '#')
         continue;
      const auto last = line.find_last_not_of(" \t\r\n");
      regions.push_back(line.substr(first, last - first + 1));
   }
   return regions;
}

// {0,3,6,9} is a curated subset of kSyntheticRegions' 18 entries; it means nothing
// for a user-supplied --regions-file, so default to every region in that case.
static std::vector<int> SelectedRegions(const benchutil::BenchmarkConfig &cfg, size_t region_count)
{
   if (cfg.allRegions || (!cfg.regionsFile.empty() && cfg.regions.empty())) {
      std::vector<int> v(region_count);
      std::iota(v.begin(), v.end(), 0);
      return v;
   }
   if (!cfg.regions.empty())
      return cfg.regions;
   return {0, 3, 6, 9};
}

int main(int argc, char **argv)
{
   benchutil::BenchmarkConfig cfg = benchutil::BenchmarkConfig::FromArgs(&argc, argv);

   // Resolve (or convert once) the files the queries run against. EnsureBam() returns an
   // empty string if no sorted, indexed BAM could be produced, in which case the BAM rows
   // are simply not registered rather than failing the whole run.
   const std::string ttree_root = cfg.EnsureTTreeRoot();
   const std::string rntuple_root = cfg.EnsureRNTupleRoot();
   const std::string bam_file = cfg.EnsureBam();
   if (bam_file.empty())
      std::cerr << "region_query_benchmark: no usable BAM - skipping the BAM rows\n";

   // The ROOT views seek via the sparse index and stop at the first record past the
   // region end, so they only see every overlapping read when the input is coordinate
   // sorted. GenerateSAMFile() emits SO:unsorted, so warn rather than let the resulting
   // short read counts be compared against BAM's complete ones.
   if (!benchutil::FileIsCoordinateSorted(cfg.EnsureSam())) {
      std::cerr << "region_query_benchmark: WARNING - input is not coordinate sorted.\n"
                << "  The TTree and RNTuple queries stop scanning at the first record past the\n"
                << "  region end, so they will under-report matches and are NOT comparable with\n"
                << "  the BAM rows (which are sorted and indexed). Pass a coordinate-sorted\n"
                << "  --sam=PATH for a meaningful comparison.\n";
   }

   const std::vector<std::string> regions =
      cfg.regionsFile.empty() ? kSyntheticRegions : LoadRegionsFromFile(cfg.regionsFile);
   if (regions.empty()) {
      std::cerr << "region_query_benchmark: no regions to run (check --regions-file)\n";
      return 1;
   }

   for (int idx : SelectedRegions(cfg, regions.size())) {
      const std::string &region = regions[idx % regions.size()];

      benchmark::RegisterBenchmark("RegionQuery/TTree/r" + std::to_string(idx), [ttree_root, region,
                                                                                 idx](benchmark::State &state) {
         RunQuery(
            state,
            [&] {
               return static_cast<int64_t>(ramview(ttree_root.c_str(), region.c_str(), /*cache=*/true,
                                                   /*perfstats=*/false, /*perfstatsfilename=*/"perf.root"));
            },
            idx);
      })->Unit(benchmark::kSecond);

      benchmark::RegisterBenchmark("RegionQuery/RNTuple/r" + std::to_string(idx), [rntuple_root, region,
                                                                                   idx](benchmark::State &state) {
         RunQuery(
            state,
            [&] {
               return static_cast<int64_t>(ramntupleview(rntuple_root.c_str(), region.c_str(), /*cache=*/true,
                                                         /*perfstats=*/false, /*perfstatsfilename=*/"perf.root"));
            },
            idx);
      })->Unit(benchmark::kSecond);

      if (!bam_file.empty()) {
         benchmark::RegisterBenchmark("RegionQuery/BAM/r" + std::to_string(idx), [bam_file, region,
                                                                                  idx](benchmark::State &state) {
            RunQuery(
               state,
               [&] {
                  const benchutil::BamReader reader(bam_file);
                  return reader.Ok() ? reader.Query(region) : -1;
               },
               idx);
         })->Unit(benchmark::kSecond);
      }
   }

   benchmark::Initialize(&argc, argv);
   benchmark::RunSpecifiedBenchmarks();
   benchmark::Shutdown();
   return 0;
}
