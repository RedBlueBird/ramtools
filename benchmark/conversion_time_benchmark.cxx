#include "bam_ops.h"
#include "benchmark_config.h"
#include "benchmark_utils.h"
#include "generate_sam_benchmark.h"
#include "ramcore/SamToNTuple.h"
#include "ramcore/SamToTTree.h"
#include <benchmark/benchmark.h>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>

// Time SAM->TTree conversion of an on-disk SAM, recording the output size.
static void TimeTTreeConversion(benchmark::State &state, const std::string &sam_file)
{
   const std::string out = "conv_ttree_out.root";
   for ([[maybe_unused]] auto _ : state) {
      {
         benchutil::ScopedStdoutSuppressor quiet;
         samtoram(sam_file.c_str(), out.c_str(), /*index=*/true, /*split=*/true, /*cache=*/true,
                  /*compression_algorithm=*/1, /*quality_policy=*/0);
      }
      if (std::filesystem::exists(out))
         state.counters["file_size_mb"] = static_cast<double>(std::filesystem::file_size(out)) / (1024.0 * 1024.0);
      std::remove(out.c_str());
   }
}

// Time SAM->RNTuple conversion of an on-disk SAM, recording the output size.
static void
TimeRNTupleConversion(benchmark::State &state, const std::string &sam_file, int compression, unsigned int quality)
{
   const std::string out = "conv_rntuple_out.root";
   for ([[maybe_unused]] auto _ : state) {
      {
         benchutil::ScopedStdoutSuppressor quiet;
         samtoramntuple(sam_file.c_str(), out.c_str(), /*index=*/true, /*split=*/true, /*cache=*/true, compression,
                        quality);
      }
      if (std::filesystem::exists(out))
         state.counters["file_size_mb"] = static_cast<double>(std::filesystem::file_size(out)) / (1024.0 * 1024.0);
      std::remove(out.c_str());
   }
}

// Time SAM->BAM conversion of an on-disk SAM, recording the output size.
//
// "Conversion" here means the same thing it means for the ROOT backends: produce a file
// that can answer a region query. For RNTuple that is one pass with index=true; for BAM
// it is the write plus the .bai (and, for an input that is not already coordinate sorted
// -- which the synthetic generator's is not -- the sort that indexing requires). The
// `sorted_in_loop` counter records whether that sort was paid, so the two run modes are
// not silently conflated.
static void TimeBamConversion(benchmark::State &state, const std::string &sam_file, int threads)
{
   const std::string out = "conv_bam_out.bam";
   for ([[maybe_unused]] auto _ : state) {
      const auto res = benchutil::BamConvert(sam_file, out, threads);
      if (!res.ok) {
         state.SkipWithError("BAM conversion failed");
         break;
      }
      state.counters["file_size_mb"] = static_cast<double>(res.bytes) / (1024.0 * 1024.0);
      state.counters["indexed"] = res.indexed ? 1 : 0;
      state.counters["sorted_in_loop"] = res.sorted ? 1 : 0;
      std::remove(out.c_str());
      std::remove((out + ".bai").c_str());
   }
}

// Generated-data sweep: build a synthetic SAM of state.range(0) reads, then time conversion.
static void BM_TTreeGenerated(benchmark::State &state)
{
   const int num_reads = static_cast<int>(state.range(0));
   const std::string sam = "conv_gen_ttree_" + std::to_string(num_reads) + ".sam";
   GenerateSAMFile(sam, num_reads);
   TimeTTreeConversion(state, sam);
   std::remove(sam.c_str());
   state.counters["reads_per_second"] = benchmark::Counter(num_reads, benchmark::Counter::kIsRate);
}

static void BM_RNTupleGenerated(benchmark::State &state, int compression, unsigned int quality)
{
   const int num_reads = static_cast<int>(state.range(0));
   const std::string sam = "conv_gen_rntuple_" + std::to_string(num_reads) + ".sam";
   GenerateSAMFile(sam, num_reads);
   TimeRNTupleConversion(state, sam, compression, quality);
   std::remove(sam.c_str());
   state.counters["reads_per_second"] = benchmark::Counter(num_reads, benchmark::Counter::kIsRate);
}

static void BM_BamGenerated(benchmark::State &state, int threads)
{
   const int num_reads = static_cast<int>(state.range(0));
   const std::string sam = "conv_gen_bam_" + std::to_string(num_reads) + ".sam";
   GenerateSAMFile(sam, num_reads);
   TimeBamConversion(state, sam, threads);
   std::remove(sam.c_str());
   state.counters["reads_per_second"] = benchmark::Counter(num_reads, benchmark::Counter::kIsRate);
}

int main(int argc, char **argv)
{
   benchutil::BenchmarkConfig cfg = benchutil::BenchmarkConfig::FromArgs(&argc, argv);

   std::cout << "Individual Conversion Time Benchmark\n";
   std::cout << "====================================\n";
   std::cout << "Measuring TTree, RNTuple and BAM conversion times separately\n";
   if (cfg.quality != 0)
      std::cout << "NOTE: --quality=" << cfg.quality << " is lossy; BAM has no lossy mode, so the BAM row\n"
                << "      is not a like-for-like comparison at this setting.\n";
   std::cout << '\n';

   const int compression = cfg.compression;
   const unsigned int quality = cfg.quality;
   const int threads = cfg.threads;

   if (cfg.HasRealDataset()) {
      const std::string sam = cfg.sam;
      benchmark::RegisterBenchmark("TTree_Conversion/real", [sam](benchmark::State &state) {
         TimeTTreeConversion(state, sam);
      })->Unit(benchmark::kMillisecond);
      benchmark::RegisterBenchmark("RNTuple_Conversion/real", [sam, compression, quality](benchmark::State &state) {
         TimeRNTupleConversion(state, sam, compression, quality);
      })->Unit(benchmark::kMillisecond);
      benchmark::RegisterBenchmark("BAM_Conversion/real", [sam, threads](benchmark::State &state) {
         TimeBamConversion(state, sam, threads);
      })->Unit(benchmark::kMillisecond);
   } else {
      benchmark::RegisterBenchmark("TTree_Conversion", BM_TTreeGenerated)
         ->Arg(1000)
         ->Arg(10000)
         ->Arg(100000)
         ->Unit(benchmark::kMillisecond);
      benchmark::RegisterBenchmark(
         "RNTuple_Conversion",
         [compression, quality](benchmark::State &state) { BM_RNTupleGenerated(state, compression, quality); })
         ->Arg(1000)
         ->Arg(10000)
         ->Arg(100000)
         ->Unit(benchmark::kMillisecond);
      benchmark::RegisterBenchmark("BAM_Conversion",
                                   [threads](benchmark::State &state) { BM_BamGenerated(state, threads); })
         ->Arg(1000)
         ->Arg(10000)
         ->Arg(100000)
         ->Unit(benchmark::kMillisecond);
   }

   benchmark::Initialize(&argc, argv);
   benchmark::RunSpecifiedBenchmarks();
   benchmark::Shutdown();
   return 0;
}
