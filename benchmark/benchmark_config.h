#ifndef BENCHMARK_CONFIG_H
#define BENCHMARK_CONFIG_H

#include "rntuple/RAMNTupleRecord.h"

#include <cstdint>
#include <string>
#include <vector>

/// The command-line flags a benchmark binary accepts, besides Google
/// Benchmark's own. Combine them with | to say which flags a binary takes.
enum BenchmarkFlag : uint16_t {
   kFlagCompression = 1U << 0, ///< --compression=N
   kFlagQuality = 1U << 1,     ///< --quality=lossless|fqzcomp|illumina|drop
   kFlagReads = 1U << 2,       ///< --reads=N[,N...]
   kFlagThreads = 1U << 3,     ///< --threads=N[,N...]
   kFlagRam = 1U << 4,         ///< --ram=FILE, or $RAMTOOLS_BENCH_RNTUPLE
   kFlagRegions = 1U << 5,     ///< --regions=REGION[,REGION...], or $RAMTOOLS_BENCH_REGIONS
   kFlagSam = 1U << 6,         ///< --sam=FILE, instead of synthetic records
   kFlagWorkDir = 1U << 7,     ///< --work-dir=DIR
   kFlagRef = 1U << 8,         ///< --ref=FASTA, which CRAM is encoded against
   kFlagFormats = 1U << 9,     ///< --formats=ram,bam,cram,sam
};

/// What a benchmark runs on. The members a binary does not take a flag for
/// keep the defaults it set before calling ParseBenchmarkFlags().
struct BenchmarkOptions {
   int compression = 505; ///< ROOT compression code, algorithm*100+level
   /// The converters' quality policy, as samtoramntuple's flags set it.
   uint32_t quality_policy = RAMNTupleRecord::kPhred33;
   std::vector<int> reads;   ///< Synthetic record counts, one benchmark each
   std::vector<int> threads; ///< Thread counts, one benchmark each
   std::string ram;          ///< A RAM file written by samtoramntuple
   std::vector<std::string> regions;
   std::string sam;      ///< A SAM file to run on; empty for synthetic records, one input per count in reads
   std::string work_dir; ///< Where inputs and outputs are written; empty for a temporary directory
   std::string ref;      ///< The reference FASTA for CRAM; empty for no CRAM
   /// The formats to measure, each of "ram", "bam", "cram" and "sam" at most
   /// once; empty for every format the inputs allow.
   std::vector<std::string> formats;
};

/// Whether \p format is among \p opts.formats, or \p opts.formats is empty.
bool WantsFormat(const BenchmarkOptions &opts, const std::string &format);

/// Takes the flags in \p accepted out of argv, lowering \p argc, and sets
/// \p opts from them; the rest of argv is left for benchmark::Initialize().
/// A flag that is not accepted stays in argv, so Google Benchmark reports it as
/// unrecognized. With --help, prints the accepted flags and leaves --help in
/// argv for Google Benchmark to print its own. False, with the reason printed,
/// when a value is malformed, a --sam or --ram file does not exist, or inputs
/// that exclude each other are given: --sam with --reads, --ram with either.
bool ParseBenchmarkFlags(int *argc, char **argv, unsigned accepted, BenchmarkOptions &opts);

#endif // BENCHMARK_CONFIG_H
