#include "bam_ops.h"
#include "benchmark_config.h"
#include "benchmark_utils.h"
#include "ramcore/RecordBlocks.h"
#include "rntuple/RAMFile.h"
#include "rntuple/RAMNTupleRecord.h"
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleTypes.hxx>
#include <ROOT/RNTupleView.hxx>
#include <algorithm>
#include <benchmark/benchmark.h>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

// Scans of every record of a file, each reading what its format must read to
// answer a question. RAM reads only the columns asked for: FLAG alone, the four
// columns `ramstats flagstat` reads, or FLAG with a typed tag column. A BAM or
// SAM record is always decoded whole; a CRAM can decode only some fields, so it
// gets a row with flagstat's fields too. The full scans decode every field of
// every format, RAM's names and qualities blocks included.
//
// Every row adds up FLAG over the file, and its label shows the sum and the
// records read: rows of the same input must agree. Every scan opens its file
// for the call, as the region queries do.

// RAM, reading FLAG and \p columns of "record" and nothing else.
static void BM_RamColumns(benchmark::State &state, const std::string &ram, const std::vector<std::string> &columns)
{
   uint64_t flag_sum = 0;
   uint64_t rows = 0;
   for ([[maybe_unused]] const auto &_ : state) {
      const auto file = RAMFile::Open(ram);
      if (!file) {
         state.SkipWithError("cannot open " + ram);
         break;
      }
      auto &reader = file->GetReader();
      auto flag = reader.GetView<uint16_t>("record.flag");
      std::vector<ROOT::RNTupleView<void>> others;
      others.reserve(columns.size());
      for (const std::string &column : columns)
         others.push_back(reader.GetView<void>("record." + column));

      rows = reader.GetNEntries();
      flag_sum = 0;
      for (ROOT::NTupleSize_t row = 0; row < rows; ++row) {
         flag_sum += flag(row);
         for (auto &view : others)
            view(row);
      }
      benchmark::DoNotOptimize(flag_sum);
   }
   state.SetItemsProcessed(static_cast<int64_t>(state.iterations() * rows));
   state.counters["columns"] = static_cast<double>(1 + columns.size());
   state.counters["file_mb"] = FileMiB(ram);
   state.SetLabel(std::to_string(rows) + " records, FLAG sum " + std::to_string(flag_sum));
}

// RAM, decoding every column of "record", whatever columns the file has, and
// every name and quality from their blocks: what a whole record costs.
static void BM_RamFull(benchmark::State &state, const std::string &ram)
{
   uint64_t flag_sum = 0;
   uint64_t rows = 0;
   for ([[maybe_unused]] const auto &_ : state) {
      const auto file = RAMFile::Open(ram);
      if (!file) {
         state.SkipWithError("cannot open " + ram);
         break;
      }
      auto &reader = file->GetReader();
      // QualityBlockReader takes QUAL from the record when it is not in a block,
      // so those two columns are read into a record; FLAG is summed.
      RAMNTupleRecord rec;
      auto flag = reader.GetView("record.flag", &rec.flag);
      auto qual = reader.GetView("record.qual", &rec.qual);
      auto compression_flags = reader.GetView("record.compression_flags", &rec.compression_flags);
      std::vector<ROOT::RNTupleView<void>> others;
      const auto &desc = reader.GetDescriptor();
      for (const auto &field : desc.GetFieldIterable(desc.FindFieldId("record"))) {
         const std::string &name = field.GetFieldName();
         if (name != "flag" && name != "qual" && name != "compression_flags")
            others.push_back(reader.GetView<void>("record." + name));
      }
      NameBlockReader names(*file);
      QualityBlockReader quals(*file);

      rows = reader.GetNEntries();
      flag_sum = 0;
      size_t bytes = 0;
      for (ROOT::NTupleSize_t row = 0; row < rows; ++row) {
         flag(row);
         qual(row);
         compression_flags(row);
         for (auto &view : others)
            view(row);
         bytes += names.Get(row).size() + quals.Get(rec, row).size();
         flag_sum += rec.flag;
      }
      benchmark::DoNotOptimize(bytes);
      benchmark::DoNotOptimize(flag_sum);
   }
   state.SetItemsProcessed(static_cast<int64_t>(state.iterations() * rows));
   state.counters["file_mb"] = FileMiB(ram);
   state.SetLabel(std::to_string(rows) + " records, FLAG sum " + std::to_string(flag_sum));
}

// SAM, BAM or CRAM, through htslib.
static void BM_HtsScan(benchmark::State &state, const std::string &path, const std::string &ref, bool flagstat_fields)
{
   ScanTotals totals;
   for ([[maybe_unused]] const auto &_ : state) {
      if (!ScanRecords(path, ref, flagstat_fields, totals)) {
         state.SkipWithError("cannot scan " + path);
         break;
      }
   }
   state.SetItemsProcessed(state.iterations() * totals.records);
   state.counters["file_mb"] = FileMiB(path);
   state.SetLabel(std::to_string(totals.records) + " records, FLAG sum " + std::to_string(totals.flag_sum));
}

int main(int argc, char **argv)
{
   BenchmarkOptions opts;
   opts.reads = {1000000};
   if (!ParseBenchmarkFlags(&argc, argv, kFlagRam | kFlagSam | kFlagReads | kFlagWorkDir | kFlagRef | kFlagFormats,
                            opts))
      return 1;
   // Before any input is prepared, so that a mistyped flag or --help returns at once.
   ::benchmark::Initialize(&argc, argv);
   if (::benchmark::ReportUnrecognizedArguments(argc, argv))
      return 1;
   if (opts.reads.size() != 1) {
      std::cerr << "columnar_read_benchmark scans one file per format: give one --reads count\n";
      return 1;
   }
   if (!opts.ram.empty() && !opts.formats.empty() &&
       (WantsFormat(opts, "bam") || WantsFormat(opts, "cram") || WantsFormat(opts, "sam"))) {
      std::cerr << "columnar_read_benchmark: --ram measures RAM only; give --sam for the other formats\n";
      return 1;
   }
   const WorkDir dir(opts.work_dir);
   const int threads = std::max(1, static_cast<int>(std::thread::hardware_concurrency()));

   // Each format's file, made once from the SAM before any timing.
   std::string ram = opts.ram;
   std::string bam;
   std::string cram;
   std::string sam;
   if (ram.empty()) {
      const SamInput input = SamInputs(opts, dir, /*sorted=*/false).front();
      PrepareSamInput(input);
      ram = dir.File(input.stem + ".ram");
      if (!ConvertToRam(input.path, ram, threads))
         return 1;
      if (WantsFormat(opts, "bam")) {
         bam = dir.File(input.stem + ".bam");
         if (!ConvertAlignments(input.path, bam, threads))
            return 1;
      }
      if (WantsFormat(opts, "cram") && !opts.ref.empty()) {
         cram = dir.File(input.stem + ".cram");
         if (!ConvertAlignments(input.path, cram, threads, opts.ref))
            return 1;
      }
      if (WantsFormat(opts, "sam"))
         sam = input.path;
   }

   const std::vector<std::string> flagstat = {"refid", "refnext", "mapq"};
   if (WantsFormat(opts, "ram")) {
      ::benchmark::RegisterBenchmark("ColumnScan/RAM/flag", BM_RamColumns, ram, std::vector<std::string>{})
         ->Unit(benchmark::kMillisecond);
      ::benchmark::RegisterBenchmark("ColumnScan/RAM/flagstat", BM_RamColumns, ram, flagstat)
         ->Unit(benchmark::kMillisecond);
      ::benchmark::RegisterBenchmark("ColumnScan/RAM/flag+tag_nm", BM_RamColumns, ram,
                                     std::vector<std::string>{"tag_nm"})
         ->Unit(benchmark::kMillisecond);
   }
   if (!cram.empty())
      ::benchmark::RegisterBenchmark("ColumnScan/CRAM/flagstat", BM_HtsScan, cram, opts.ref, /*flagstat_fields=*/true)
         ->Unit(benchmark::kMillisecond);
   if (WantsFormat(opts, "ram"))
      ::benchmark::RegisterBenchmark("FullScan/RAM", BM_RamFull, ram)->Unit(benchmark::kMillisecond);
   if (!bam.empty())
      ::benchmark::RegisterBenchmark("FullScan/BAM", BM_HtsScan, bam, std::string(), /*flagstat_fields=*/false)
         ->Unit(benchmark::kMillisecond);
   if (!cram.empty())
      ::benchmark::RegisterBenchmark("FullScan/CRAM", BM_HtsScan, cram, opts.ref, /*flagstat_fields=*/false)
         ->Unit(benchmark::kMillisecond);
   if (!sam.empty())
      ::benchmark::RegisterBenchmark("FullScan/SAM", BM_HtsScan, sam, std::string(), /*flagstat_fields=*/false)
         ->Unit(benchmark::kMillisecond);

   ::benchmark::RunSpecifiedBenchmarks();
   ::benchmark::Shutdown();
   return 0;
}
