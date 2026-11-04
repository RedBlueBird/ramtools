#ifndef BENCHMARK_BAM_OPS_H
#define BENCHMARK_BAM_OPS_H

// BAM operations for the benchmarks, done with htslib in this process rather
// than by starting samtools: no process start-up is timed, threads show in the
// wall time as they do for the RAM side, and no samtools binary is needed.
// BAM keeps qualities losslessly, so it holds the same data as a RAM file only
// with --quality=lossless.

#include <string>

/// Writes \p bam, coordinate sorted, and its .bai index from the SAM or BAM
/// \p input, on \p threads threads. An input that is not coordinate sorted is
/// sorted in memory, which is refused beyond kMaxSortRecords. False, with the
/// reason printed, on failure.
bool SamToSortedBam(const std::string &input, const std::string &bam, int threads);

/// The most records SamToSortedBam() sorts in memory.
constexpr long kMaxSortRecords = 10000000;

/// Writes one BAM per reference that has records, <prefix><name>.bam, from the
/// sorted, indexed \p bam, with \p threads threads extracting references side by
/// side. Records with no reference are left out, as the RAM split leaves them
/// out. False, with the reason printed, on failure.
bool SplitBam(const std::string &bam, const std::string &prefix, int threads);

#endif // BENCHMARK_BAM_OPS_H
