// Chromosome split: RNTuple's parallel per-chromosome writer vs the same operation on BAM.
//
// Both sides start from the same SAM and finish with one file per populated chromosome,
// so the comparison is end-to-end and symmetric:
//
//   RNTuple  samtoramntuple_split_by_chromosome(): parse SAM -> group by chromosome ->
//            sort each chromosome -> parallel write, one .root per chromosome.
//   BAM      BamConvert(): parse SAM -> sort -> write BAM -> build .bai; then BamSplit():
//            indexed extraction into one .bam per chromosome across `threads` workers.
//
// The BAM route writes an intermediate sorted+indexed BAM that the RNTuple route does not
// need. That is a real cost of the format's workflow rather than an artefact of the
// benchmark -- it is what a samtools user actually pays -- so it is inside the timed
// region. Only the per-chromosome outputs are counted in size_MB; the intermediate is
// excluded and deleted.
//
// This replaces an earlier version that shelled out to `samtools view/sort/index` via
// system(). Running in-process takes fork/exec out of every iteration, puts the sort on
// the same footing as the RNTuple side's internal sort, and drops the suite's dependency
// on a samtools binary being installed.

#include "bam_ops.h"
#include "benchmark_config.h"
#include "benchmark_utils.h"
#include "generate_sam_benchmark.h"
#include "ramcore/SamToNTuple.h"
#include <benchmark/benchmark.h>
#include <cstdio>
#include <iostream>
#include <string>

namespace {

// Empty => generate synthetic data per benchmark arg; non-empty => split this real SAM.
std::string g_realSam;      // NOLINT(cppcoreguidelines-avoid-non-const-global-variables)
int g_compression = 505;    // NOLINT(cppcoreguidelines-avoid-non-const-global-variables)
unsigned int g_quality = 0; // NOLINT(cppcoreguidelines-avoid-non-const-global-variables)

constexpr double kBytesPerMB = 1024.0 * 1024.0;

// Resolve the SAM to operate on: the real dataset, or a freshly generated synthetic file.
// Sets `generated` so the caller knows whether to delete it afterwards.
std::string PrepareSam(int num_reads, const std::string &gen_name, bool &generated)
{
   if (!g_realSam.empty()) {
      generated = false;
      return g_realSam;
   }
   generated = true;
   GenerateSAMFile(gen_name, num_reads);
   return gen_name;
}

void BM_RNTupleSplit(benchmark::State &state)
{
   const int num_reads = static_cast<int>(state.range(0));
   const int num_threads = static_cast<int>(state.range(1));
   bool generated = false;
   const std::string sam_file =
      PrepareSam(num_reads, "split_rntuple_in_" + std::to_string(num_reads) + ".sam", generated);
   const std::string prefix = "split_rntuple_out";

   for ([[maybe_unused]] auto _ : state) {
      {
         benchutil::ScopedStdoutSuppressor quiet(/*suppress_stderr=*/true);
         // NOTE: the compression argument is currently ignored by the implementation,
         // which hardcodes ZSTD-1 for the split path; it is passed for correctness of the
         // call site, not because it takes effect.
         samtoramntuple_split_by_chromosome(sam_file.c_str(), prefix.c_str(), g_compression, g_quality, num_threads);
      }

      state.counters["size_MB"] = static_cast<double>(benchutil::TotalSizeByPrefix(prefix + "_")) / kBytesPerMB;
      state.counters["threads"] = num_threads;
      benchutil::CleanupByPrefix(prefix + "_");
   }

   if (generated)
      std::remove(sam_file.c_str());
   if (g_realSam.empty())
      state.counters["reads/s"] = benchmark::Counter(num_reads, benchmark::Counter::kIsRate);
}

void BM_BamSplit(benchmark::State &state)
{
   const int num_reads = static_cast<int>(state.range(0));
   const int num_threads = static_cast<int>(state.range(1));
   bool generated = false;
   const std::string sam_file = PrepareSam(num_reads, "split_bam_in_" + std::to_string(num_reads) + ".sam", generated);
   const std::string intermediate = "split_bam_tmp.bam";
   const std::string prefix = "split_bam_out";

   for ([[maybe_unused]] auto _ : state) {
      const auto converted = benchutil::BamConvert(sam_file, intermediate, num_threads);
      if (!converted.ok || !converted.indexed) {
         state.SkipWithError("could not build a sorted, indexed BAM to split");
         break;
      }
      const std::size_t bytes = benchutil::BamSplit(intermediate, prefix, num_threads);

      state.counters["size_MB"] = static_cast<double>(bytes) / kBytesPerMB;
      state.counters["threads"] = num_threads;
      benchutil::CleanupByPrefix(prefix + "_");
      std::remove(intermediate.c_str());
      std::remove((intermediate + ".bai").c_str());
   }

   if (generated)
      std::remove(sam_file.c_str());
   if (g_realSam.empty())
      state.counters["reads/s"] = benchmark::Counter(num_reads, benchmark::Counter::kIsRate);
}

} // namespace

int main(int argc, char **argv)
{
   benchutil::BenchmarkConfig cfg = benchutil::BenchmarkConfig::FromArgs(&argc, argv);
   g_realSam = cfg.sam;
   g_compression = cfg.compression;
   g_quality = cfg.quality;

   if (g_quality != 0)
      std::cout << "NOTE: --quality=" << g_quality << " is lossy; BAM has no lossy mode, so the BAM rows\n"
                << "      are not a like-for-like comparison at this setting.\n\n";

   if (!g_realSam.empty()) {
      // Real dataset: vary only the thread count (read count is fixed by the file).
      benchmark::RegisterBenchmark("ChromosomeSplit/RNTuple/real", BM_RNTupleSplit)
         ->Args({0, 2})
         ->Args({0, 4})
         ->Unit(benchmark::kMillisecond);
      benchmark::RegisterBenchmark("ChromosomeSplit/BAM/real", BM_BamSplit)
         ->Args({0, 2})
         ->Args({0, 4})
         ->Unit(benchmark::kMillisecond);
   } else {
      benchmark::RegisterBenchmark("ChromosomeSplit/RNTuple", BM_RNTupleSplit)
         ->Args({100000, 2})
         ->Args({100000, 4})
         ->Args({500000, 2})
         ->Args({500000, 4})
         ->Args({1000000, 2})
         ->Args({1000000, 4})
         ->Unit(benchmark::kMillisecond);
      benchmark::RegisterBenchmark("ChromosomeSplit/BAM", BM_BamSplit)
         ->Args({100000, 2})
         ->Args({100000, 4})
         ->Args({500000, 2})
         ->Args({500000, 4})
         ->Args({1000000, 2})
         ->Args({1000000, 4})
         ->Unit(benchmark::kMillisecond);
   }

   benchmark::Initialize(&argc, argv);
   benchmark::RunSpecifiedBenchmarks();
   benchmark::Shutdown();
   return 0;
}
