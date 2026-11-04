#include <benchmark/benchmark.h>
#include "bam_ops.h"
#include "benchmark_config.h"
#include "benchmark_utils.h"
#include "ramcore/BamtoNTuple.h"
#include "ramcore/SamToNTuple.h"
#include "rntuple/RAMNTupleRecord.h"
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>
#include <system_error>
#include <vector>

// Conversion of each input to each format, timed from the input file to the
// complete output: SAM to RAM, BAM and CRAM, as samtoramntuple, `samtools view
// -b` and `samtools view -C -T ref` write them, and BAM to RAM, as
// bamtoramntuple does. Every row reports its output's size and the SAM's, so
// the formats compare by size as well as by time.

// The counters every conversion reports, or an error when \p written, the
// records in the output, is not the input's count: a time for a conversion
// that lost records is not a result. The rates are of the input's records and
// bytes, per second of wall time.
static void Report(benchmark::State &state, const SamInput &input, const std::string &source, const std::string &output,
                   int64_t written)
{
   const int64_t records = input.reads > 0 ? input.reads : CountSamRecords(input.path);
   if (written != records) {
      state.SkipWithError(output + " holds " + std::to_string(written) + " records, " + input.path + " " +
                          std::to_string(records));
      return;
   }
   state.counters["file_size_mb"] = FileMiB(output);
   state.counters["sam_size_mb"] = FileMiB(input.path);
   state.counters["reads"] = static_cast<double>(records);
   state.counters["reads_per_second"] = benchmark::Counter(static_cast<double>(records), benchmark::Counter::kIsRate);
   std::error_code ec;
   const auto source_size = std::filesystem::file_size(source, ec);
   if (!ec)
      state.SetBytesProcessed(state.iterations() * static_cast<int64_t>(source_size));
}

static void
BM_SamToRam(benchmark::State &state, const SamInput &input, const std::string &output, const BenchmarkOptions &opts)
{
   const int threads = static_cast<int>(state.range(0));
   PrepareSamInput(input);

   const PeakMemory memory;
   for ([[maybe_unused]] const auto &_ : state) {
      bool converted = false;
      {
         const StdoutSilencer silence;
         converted = samtoramntuple(input.path.c_str(), output.c_str(), opts.compression, opts.quality_policy, threads);
      }
      if (!converted) {
         state.SkipWithError("samtoramntuple failed on " + input.path);
         return;
      }
   }
   memory.Report(state);
   Report(state, input, input.path, output, RamRecords(output));
   std::remove(output.c_str());
}

// SAM to BAM or CRAM, by the output's extension.
static void
BM_SamToHts(benchmark::State &state, const SamInput &input, const std::string &output, const BenchmarkOptions &opts)
{
   const int threads = static_cast<int>(state.range(0));
   PrepareSamInput(input);

   const PeakMemory memory;
   for ([[maybe_unused]] const auto &_ : state) {
      if (!ConvertAlignments(input.path, output, threads, opts.ref)) {
         state.SkipWithError("the conversion of " + input.path + " to " + output + " failed");
         return;
      }
   }
   memory.Report(state);
   ScanTotals written;
   if (!ScanRecords(output, opts.ref, /*flagstat_fields=*/false, written))
      written.records = -1;
   Report(state, input, input.path, output, written.records);
   std::remove(output.c_str());
}

// BAM to RAM. The BAM is made from the SAM before the timing starts, once for
// all the thread counts.
static void BM_BamToRam(benchmark::State &state, const SamInput &input, const std::string &bam,
                        const std::string &output, const BenchmarkOptions &opts)
{
   PrepareSamInput(input);
   if (!MakeOnce(bam, [&] { return ConvertAlignments(input.path, bam, /*threads=*/1); })) {
      state.SkipWithError("cannot make a BAM of " + input.path);
      return;
   }
   // The BAM is read on htslib's threads whatever this is; ROOT's implicit
   // threads compress the pages, as bamtoramntuple -threads N enables them.
   SetImplicitMT(static_cast<int>(state.range(0)));

   const PeakMemory memory;
   for ([[maybe_unused]] const auto &_ : state) {
      std::remove(output.c_str());
      {
         const StdoutSilencer silence;
         bamtoramntuple(bam.c_str(), output.c_str(), /*split=*/false, /*cache=*/true, opts.compression,
                        opts.quality_policy);
      }
      // bamtoramntuple() returns nothing: a failure shows as no output.
      if (FileMiB(output) == 0.0) {
         SetImplicitMT(/*threads=*/1);
         state.SkipWithError("bamtoramntuple failed on " + bam);
         return;
      }
   }
   memory.Report(state);
   SetImplicitMT(/*threads=*/1);
   Report(state, input, bam, output, RamRecords(output));
   std::remove(output.c_str());
}

int main(int argc, char **argv)
{
   BenchmarkOptions opts;
   opts.reads = {1000, 10000, 100000};
   opts.threads = {1};
   if (!ParseBenchmarkFlags(&argc, argv,
                            kFlagCompression | kFlagQuality | kFlagReads | kFlagThreads | kFlagSam | kFlagWorkDir |
                               kFlagRef | kFlagFormats,
                            opts))
      return 1;
   // Before any input is prepared, so that a mistyped flag or --help returns at once.
   ::benchmark::Initialize(&argc, argv);
   if (::benchmark::ReportUnrecognizedArguments(argc, argv))
      return 1;
   StartROOT();
   const WorkDir dir(opts.work_dir);

   std::cout << "Conversion Time Benchmark" << std::endl;
   std::cout << std::endl;
   if (opts.quality_policy != RAMNTupleRecord::kPhred33 && (WantsFormat(opts, "bam") || WantsFormat(opts, "cram")))
      std::cout << "Note: BAM and CRAM keep qualities losslessly, so only --quality=lossless compares the RAM "
                   "rows with theirs like for like.\n\n";

   const bool ram = WantsFormat(opts, "ram");
   const bool bam = WantsFormat(opts, "bam");
   const bool cram = WantsFormat(opts, "cram") && !opts.ref.empty();
   for (const SamInput &input : SamInputs(opts, dir, /*sorted=*/false)) {
      std::vector<benchmark::internal::Benchmark *> benches;
      if (ram)
         benches.push_back(::benchmark::RegisterBenchmark("SamToRam/" + input.name, BM_SamToRam, input,
                                                          dir.File(input.stem + "_out.ram"), opts));
      if (bam)
         benches.push_back(::benchmark::RegisterBenchmark("SamToBam/" + input.name, BM_SamToHts, input,
                                                          dir.File(input.stem + "_out.bam"), opts));
      if (cram)
         benches.push_back(::benchmark::RegisterBenchmark("SamToCram/" + input.name, BM_SamToHts, input,
                                                          dir.File(input.stem + "_out.cram"), opts));
      if (ram && bam)
         benches.push_back(::benchmark::RegisterBenchmark("BamToRam/" + input.name, BM_BamToRam, input,
                                                          dir.File(input.stem + ".bam"),
                                                          dir.File(input.stem + "_from_bam.ram"), opts));
      for (auto *bench : benches) {
         bench->ArgName("threads");
         for (const int threads : opts.threads)
            bench->Arg(threads);
         // CPU time counts the main thread only, which with threads mostly
         // waits: rates and iteration counts go by wall time instead.
         bench->UseRealTime()->Unit(benchmark::kMillisecond);
      }
   }

   ::benchmark::RunSpecifiedBenchmarks();
   ::benchmark::Shutdown();
   return 0;
}
