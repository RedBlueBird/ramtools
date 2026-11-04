#ifndef GENERATE_SAM_BENCHMARK_H
#define GENERATE_SAM_BENCHMARK_H

#include <string>

/// Writes num_reads random 36 bp reads over the 24 hg19 references, in random
/// order (SO:unsorted), from a random seed.
void GenerateSAMFile(const std::string &filename, int num_reads);

/// Writes num_reads random 36 bp reads over the same references in coordinate
/// order (SO:coordinate). The same seed gives the same file, so benchmarks on it
/// can be repeated, and region queries on it take the sorted path.
void GenerateSortedSAMFile(const std::string &filename, int num_reads, unsigned seed = 1);

#endif
