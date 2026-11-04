#include <benchmark/benchmark.h>
#include "bam_ops.h"
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
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

// Records overlapping each region, counted in each format: RAM with
// ramntupleview(), BAM and CRAM through their index, and SAM, which has none,
// by reading it all. Every count opens its file, and the index, for the call,
// so all columns are the same measurement. There is no whole-file row: RAM
// counts a whole file from its METADATA without reading a record, which says
// nothing against a format that decodes them; the whole-reference region is
// the scan that compares.
//
// The inputs come from flags, or from the environment as before, so the numbers
// are reproducible on any machine. A hard-coded path is not a benchmark result.
//
//   --ram=FILE      or RAMTOOLS_BENCH_RNTUPLE   RAM file written by samtoramntuple; RAM only
//   --sam=FILE                                  SAM file, converted once to each format
//   --reads=N                                   synthetic sorted records, when neither is given
//   --ref=FASTA                                 the reference, for the CRAM column
//   --formats=F,...                             the columns: ram, bam, cram, sam
//   --regions=R,... or RAMTOOLS_BENCH_REGIONS   regions; by default sizes from 1 kb to a
//                                               whole reference, see DefaultRegions()

// Counts the records of a region in one format; negative on failure.
using CountFn = std::function<int64_t(const std::string &region)>;

static void BM_Count(benchmark::State &state, const CountFn &count, const std::string &region, double file_mb)
{
   int64_t total_reads_processed = 0;
   int64_t reads_in_this_run = 0;

   for ([[maybe_unused]] const auto &_ : state) {
      reads_in_this_run = count(region);
      if (reads_in_this_run < 0) {
         state.SkipWithError("the count failed");
         break;
      }
      total_reads_processed += reads_in_this_run;
   }

   state.SetItemsProcessed(total_reads_processed);
   state.counters["file_mb"] = file_mb;
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

// A format's column: its name in the benchmark names, its file, and its count.
struct Column {
   std::string name;
   std::string path;
   CountFn count;
};

int main(int argc, char **argv)
{
   BenchmarkOptions opts;
   opts.reads = {1000000};
   if (!ParseBenchmarkFlags(
          &argc, argv, kFlagRam | kFlagRegions | kFlagSam | kFlagReads | kFlagWorkDir | kFlagRef | kFlagFormats, opts))
      return 1;
   // Before any input is prepared, so that a mistyped flag or --help returns at once.
   ::benchmark::Initialize(&argc, argv);
   if (::benchmark::ReportUnrecognizedArguments(argc, argv))
      return 1;
   if (opts.reads.size() != 1) {
      std::cerr << "region_query_benchmark queries one file: give one --reads count\n";
      return 1;
   }
   // A RAM file alone gives nothing to make the other formats from.
   if (!opts.ram.empty() && !opts.formats.empty() &&
       (WantsFormat(opts, "bam") || WantsFormat(opts, "cram") || WantsFormat(opts, "sam"))) {
      std::cerr << "region_query_benchmark: --ram measures RAM only; give --sam for the other formats\n";
      return 1;
   }
   const WorkDir dir(opts.work_dir);
   const int threads = std::max(1, static_cast<int>(std::thread::hardware_concurrency()));

   // Each format's file, made once from the SAM before any timing, as its own
   // tool makes it: samtoramntuple by default, a sorted and indexed BAM and CRAM.
   std::string ram = opts.ram;
   std::vector<Column> columns;
   if (ram.empty()) {
      const SamInput input = SamInputs(opts, dir, /*sorted=*/true).front();
      PrepareSamInput(input);
      ram = dir.File(input.stem + ".ram");
      bool converted = false;
      {
         const StdoutSilencer silence;
         converted = samtoramntuple(input.path.c_str(), ram.c_str(), /*compression_algorithm=*/505,
                                    /*quality_policy=*/RAMNTupleRecord::kPhred33, threads);
      }
      if (!converted) {
         std::cerr << "region_query_benchmark: could not convert " << input.path << "\n";
         return 1;
      }
      const std::vector<std::pair<std::string, std::string>> indexed = {{"bam", "BAM"}, {"cram", "CRAM"}};
      for (const auto &[format, name] : indexed) {
         if (!WantsFormat(opts, format) || (format == "cram" && opts.ref.empty()))
            continue;
         const std::string path = dir.File(input.stem + "." + format);
         if (!ConvertSorted(input.path, path, threads, opts.ref))
            return 1;
         const std::string ref = opts.ref;
         columns.push_back({name, path, [path, ref](const std::string &r) { return CountRegion(path, r, ref); }});
      }
      if (WantsFormat(opts, "sam")) {
         const std::string path = input.path;
         columns.push_back({"SAM", path, [path](const std::string &r) { return ScanRegion(path, r); }});
      }
   }
   if (WantsFormat(opts, "ram")) {
      columns.insert(columns.begin(),
                     Column{"RAM", ram, [ram](const std::string &r) {
                               const StdoutSilencer silence;
                               return static_cast<int64_t>(ramntupleview(ram.c_str(), r.c_str(), RAMNTupleViewOpts{}));
                            }});
   }

   std::vector<std::string> regions = opts.regions;
   {
      const auto file = RAMFile::Open(ram);
      if (!file)
         return 1;
      if (!file->GetMetadata().coordinate_sorted)
         std::cerr << "region_query_benchmark: " << ram << " is not coordinate sorted, so every RAM query reads "
                   << "the whole file\n";
      if (regions.empty())
         regions = DefaultRegions(*file);
   }
   if (regions.empty()) {
      std::cerr << "region_query_benchmark: " << ram << " has no placed records; give --regions\n";
      return 1;
   }

   // Region by region, each format's row after the other's.
   for (std::size_t i = 0; i < regions.size(); ++i) {
      for (const Column &column : columns)
         ::benchmark::RegisterBenchmark("RegionQuery/r" + std::to_string(i) + "/" + column.name, BM_Count, column.count,
                                        regions[i], FileMiB(column.path))
            ->Unit(benchmark::kMillisecond);
   }

   ::benchmark::RunSpecifiedBenchmarks();
   ::benchmark::Shutdown();
   return 0;
}
