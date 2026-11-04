#ifndef BENCHMARK_BAM_OPS_H
#define BENCHMARK_BAM_OPS_H

// BAM and CRAM operations for the benchmarks, done with htslib in this process
// rather than by starting samtools: no process start-up is timed, threads show
// in the wall time as they do for the RAM side, and no samtools binary is
// needed. BAM and CRAM keep qualities losslessly, so they hold the same data as
// a RAM file only with --quality=lossless.
//
// The output format follows the file name, as with samtools: .cram is CRAM,
// anything else BAM. A CRAM is encoded and decoded against \p ref, a FASTA,
// which every function touching a CRAM requires: without it htslib would look
// the reference up over the network.

#include <cstdint>
#include <string>

/// Writes \p output from the SAM, BAM or CRAM \p input as it comes, unsorted and
/// without an index, as `samtools view -b` (or `-C -T ref`) does, on \p threads
/// threads. False, with the reason printed, on failure.
bool ConvertAlignments(const std::string &input, const std::string &output, int threads, const std::string &ref = "");

/// The same, written in coordinate order and indexed (.bai, or .crai for a
/// CRAM). An input that is not coordinate sorted is sorted in memory, which is
/// refused beyond kMaxSortRecords.
bool ConvertSorted(const std::string &input, const std::string &output, int threads, const std::string &ref = "");

/// The most records ConvertSorted() sorts in memory.
constexpr long kMaxSortRecords = 10000000;

/// Writes one BAM per reference that has records, <prefix><name>.bam, from the
/// sorted, indexed \p bam, with \p threads threads extracting references side by
/// side. Records with no reference are left out, as the RAM split leaves them
/// out. False, with the reason printed, on failure.
bool SplitBam(const std::string &bam, const std::string &prefix, int threads);

/// The records overlapping \p region in an indexed BAM or CRAM, counted as
/// `samtools view -c` counts them. The file and its index are opened for every
/// call, as ramntupleview() opens its file. -1, with the reason printed, on
/// failure.
int64_t CountRegion(const std::string &path, const std::string &region, const std::string &ref = "");

/// The same count by reading the whole file, for SAM, which has no index: a
/// record overlaps when it is on the region's reference and its CIGAR span meets
/// the region, the rule the index uses.
int64_t ScanRegion(const std::string &path, const std::string &region);

/// What a scan of every record adds up, for comparing scans across formats.
struct ScanTotals {
   int64_t records = 0;
   uint64_t flag_sum = 0; ///< The sum of every record's FLAG
};

/// Reads every record of a SAM, BAM or CRAM and adds up \p totals. With
/// \p flagstat_fields a CRAM decodes only FLAG, RNAME, RNEXT and MAPQ, what
/// `samtools flagstat` reads; a BAM or a SAM record is always decoded whole.
/// False, with the reason printed, on failure.
bool ScanRecords(const std::string &path, const std::string &ref, bool flagstat_fields, ScanTotals &totals);

#endif // BENCHMARK_BAM_OPS_H
