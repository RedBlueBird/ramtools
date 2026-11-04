#include <benchmark/benchmark.h>
#include "benchmark_config.h"
#include "benchmark_utils.h"
#include "ramcore/RAMNTupleView.h"
#include "ramcore/SamToNTuple.h"
#include "rntuple/RAMFile.h"
#include "rntuple/RAMNTupleRecord.h"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// The file and the regions come from flags, or from the environment as before,
// so the numbers are reproducible on any machine. A hard-coded path is not a
// benchmark result.
//
//   --ram=FILE      or RAMTOOLS_BENCH_RNTUPLE   RAM file written by samtoramntuple
//   --sam=FILE                                  SAM file, converted once before the queries
//   --reads=N                                   synthetic sorted records, when neither is given
//   --regions=R,... or RAMTOOLS_BENCH_REGIONS   regions; by default sizes from 1 kb to a
//                                               whole reference, see DefaultRegions()

static void BM_RegionQuery(benchmark::State &state, const std::string &ram, const std::string &region)
{
   int64_t total_reads_processed = 0;
   std::int64_t reads_in_this_run = 0;
   const RAMNTupleViewOpts opts{};

   for ([[maybe_unused]] const auto &_ : state) {
      {
         const StdoutSilencer silence;
         reads_in_this_run = ramntupleview(ram.c_str(), region.c_str(), opts);
      }

      total_reads_processed += reads_in_this_run;
   }

   state.SetItemsProcessed(total_reads_processed);
   state.SetLabel(region + ": " + std::to_string(reads_in_this_run) + " reads");
}

// The LN of every @SQ line of the SAM header, by reference name.
static std::unordered_map<std::string, int64_t> ReferenceLengths(const RAMFileMetadata &meta)
{
   std::unordered_map<std::string, int64_t> lengths;
   for (const std::string &line : meta.header_lines) {
      if (line.rfind("@SQ\t", 0) != 0)
         continue;
      std::istringstream fields(line);
      std::string field;
      std::string sn;
      int64_t ln = 0;
      while (std::getline(fields, field, '\t')) {
         if (field.rfind("SN:", 0) == 0)
            sn = field.substr(3);
         else if (field.rfind("LN:", 0) == 0)
            ln = std::strtoll(field.substr(3).c_str(), nullptr, 10);
      }
      lengths[sn] = ln;
   }
   return lengths;
}

// Regions of 1 kb, 1 Mb and 10 Mb around the middle record of the longest
// reference that has records, then that whole reference. Taken from the file
// itself, they fit any reference naming and hold records wherever the reads
// are, gaps and centromeres included. In an unsorted file a reference's rows
// are not together, so the regions centre on the middle of the reference.
static std::vector<std::string> DefaultRegions(RAMFile &file)
{
   const RAMFileMetadata &meta = file.GetMetadata();
   const RAMRefCounts &counts = meta.counts;
   const auto lengths = ReferenceLengths(meta);

   std::string name;
   int64_t length = -1;
   uint64_t middle_row = 0;
   uint64_t first_row = 0; // In a sorted file a reference's rows follow those of the references before it.
   for (std::size_t refid = 0; refid < counts.mapped.size() && refid < meta.rname_refs.size(); ++refid) {
      const uint64_t records = counts.mapped[refid] + counts.unmapped[refid];
      const auto known = lengths.find(meta.rname_refs[refid]);
      const int64_t ln = known == lengths.end() ? 0 : known->second;
      if (records > 0 && ln > length) {
         name = meta.rname_refs[refid];
         length = ln;
         middle_row = first_row + records / 2;
      }
      first_row += records;
   }
   if (name.empty())
      return {};

   int64_t centre = length / 2;
   if (meta.coordinate_sorted) {
      auto pos = file.GetReader().GetView<int32_t>("record.pos");
      centre = static_cast<int64_t>(pos(middle_row)) + 1; // pos is 0-based
   }

   std::vector<std::string> regions;
   for (const int64_t size : {int64_t{1000}, int64_t{1000000}, int64_t{10000000}}) {
      if (size >= length)
         break;
      const int64_t start = std::clamp(centre - size / 2, int64_t{1}, length - size + 1);
      regions.push_back(name + ":" + std::to_string(start) + "-" + std::to_string(start + size - 1));
   }
   regions.push_back(name);
   return regions;
}

int main(int argc, char **argv)
{
   BenchmarkOptions opts;
   opts.reads = {1000000};
   if (!ParseBenchmarkFlags(&argc, argv, kFlagRam | kFlagRegions | kFlagSam | kFlagReads | kFlagWorkDir, opts))
      return 1;
   // Before any input is prepared, so that a mistyped flag or --help returns at once.
   ::benchmark::Initialize(&argc, argv);
   if (::benchmark::ReportUnrecognizedArguments(argc, argv))
      return 1;
   if (opts.reads.size() != 1) {
      std::cerr << "region_query_benchmark queries one file: give one --reads count\n";
      return 1;
   }
   const WorkDir dir(opts.work_dir);

   // Without a RAM file, convert the SAM once, as samtoramntuple would by default.
   std::string ram = opts.ram;
   if (ram.empty()) {
      const SamInput input = SamInputs(opts, dir, /*sorted=*/true).front();
      PrepareSamInput(input);
      ram = dir.File("region_query.ram");
      bool converted = false;
      {
         const StdoutSilencer silence;
         converted = samtoramntuple(input.path.c_str(), ram.c_str(), /*compression_algorithm=*/505,
                                    /*quality_policy=*/RAMNTupleRecord::kPhred33);
      }
      if (!converted) {
         std::cerr << "region_query_benchmark: could not convert " << input.path << "\n";
         return 1;
      }
   }

   std::vector<std::string> regions = opts.regions;
   {
      const auto file = RAMFile::Open(ram);
      if (!file)
         return 1;
      if (!file->GetMetadata().coordinate_sorted)
         std::cerr << "region_query_benchmark: " << ram << " is not coordinate sorted, so every query reads the "
                   << "whole file\n";
      if (regions.empty())
         regions = DefaultRegions(*file);
   }
   if (regions.empty()) {
      std::cerr << "region_query_benchmark: " << ram << " has no placed records; give --regions\n";
      return 1;
   }

   // One case per region, numbered in the order given.
   for (std::size_t i = 0; i < regions.size(); ++i) {
      ::benchmark::RegisterBenchmark("RegionQueryFixture/RNTuple/" + std::to_string(i), BM_RegionQuery, ram, regions[i])
         ->Unit(benchmark::kSecond);
   }

   ::benchmark::RunSpecifiedBenchmarks();
   ::benchmark::Shutdown();
   return 0;
}
