#include "generate_sam_benchmark.h"
#include <fstream>
#include <random>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cstddef>
#include <string_view>
#include <utility>

namespace {

const std::vector<std::pair<std::string, int>> &Chromosomes()
{
   static const std::vector<std::pair<std::string, int>> chromosomes = {
      {"chr1", 249250621},  {"chr2", 243199373},  {"chr3", 198022430},  {"chr4", 191154276},  {"chr5", 180915260},
      {"chr6", 171115067},  {"chr7", 159138663},  {"chr8", 146364022},  {"chr9", 141213431},  {"chr10", 135534747},
      {"chr11", 135006516}, {"chr12", 133851895}, {"chr13", 115169878}, {"chr14", 107349540}, {"chr15", 102531392},
      {"chr16", 90354753},  {"chr17", 81195210},  {"chr18", 78077248},  {"chr19", 59128983},  {"chr20", 63025520},
      {"chr21", 48129895},  {"chr22", 51304566},  {"chrM", 16571},      {"chrX", 155270560}};
   return chromosomes;
}

constexpr int kReadLength = 36;
constexpr std::string_view kBases = "ACGT";

void WriteHeader(std::ostream &out, const char *sort_order)
{
   out << "@HD\tVN:1.6\tSO:" << sort_order << "\n";
   for (const auto &[chrom, length] : Chromosomes()) {
      out << "@SQ\tSN:" << chrom << "\tLN:" << length << "\n";
   }
}

// One 36 bp record at chrom:position; every other field is drawn from rng.
void WriteRecord(std::ostream &out, std::mt19937 &rng, const std::string &chrom, int position)
{
   std::uniform_int_distribution<size_t> base_dist(0, kBases.size() - 1);
   std::uniform_int_distribution<> flag_dist(0, 1);
   std::uniform_int_distribution<> qual_dist(33, 73);
   std::uniform_int_distribution<> fc_dist(10000, 99999);
   std::uniform_int_distribution<> lane_dist(1, 8);
   std::uniform_int_distribution<> tile_dist(1, 300);
   std::uniform_int_distribution<> coord_dist(1, 1000);
   std::uniform_int_distribution<> mismatch_dist(0, 9);

   out << "SOLEXA-1GA-2_2_FC" << fc_dist(rng) << ":" << lane_dist(rng) << ":" << tile_dist(rng) << ":"
       << coord_dist(rng) << ":" << coord_dist(rng) << "\t";

   out << (flag_dist(rng) * 16) << "\t" << chrom << "\t" << position << "\t" << 25 << "\t";

   out << kReadLength << "M\t";
   out << "*\t0\t0\t";

   std::string sequence;
   sequence.reserve(kReadLength);
   for (int j = 0; j < kReadLength; ++j) {
      sequence += kBases[base_dist(rng)];
   }
   out << sequence << "\t";

   for (int j = 0; j < kReadLength; ++j) {
      out << static_cast<char>(qual_dist(rng));
   }

   int nm = mismatch_dist(rng) < 7 ? 0 : 1;
   out << "\tNM:i:" << nm;
   out << "\tX" << nm << ":i:1";
   out << "\tMD:Z:" << kReadLength;
   out << "\n";
}

// A position a whole read fits in on a reference of this length.
int DrawPosition(std::mt19937 &rng, int chrom_length)
{
   std::uniform_int_distribution<> pos_dist(1, std::max(1, chrom_length - kReadLength));
   return pos_dist(rng);
}

} // namespace

void GenerateSAMFile(const std::string &filename, int num_reads)
{
   std::random_device rd;
   std::mt19937 rng(rd());

   const auto &chromosomes = Chromosomes();
   std::ofstream out(filename);
   WriteHeader(out, /*sort_order=*/"unsorted");

   std::uniform_int_distribution<size_t> chrom_dist(0, chromosomes.size() - 1);
   std::uniform_int_distribution<> chr1_pos_dist(1, 1000000);

   for (int i = 0; i < num_reads; ++i) {
      if (i < 10) {
         WriteRecord(out, rng, "chr1", chr1_pos_dist(rng));
      } else {
         const auto &[chrom, length] = chromosomes[chrom_dist(rng)];
         WriteRecord(out, rng, chrom, DrawPosition(rng, length));
      }
   }
   out.close();
}

void GenerateSortedSAMFile(const std::string &filename, int num_reads, unsigned seed)
{
   std::mt19937 rng(seed);

   // Draw every placement first, then write the records in (reference, position)
   // order: 8 bytes a record in memory instead of the whole record.
   const auto &chromosomes = Chromosomes();
   std::uniform_int_distribution<size_t> chrom_dist(0, chromosomes.size() - 1);
   std::vector<std::pair<size_t, int>> placements;
   placements.reserve(num_reads);
   for (int i = 0; i < num_reads; ++i) {
      const size_t chrom = chrom_dist(rng);
      placements.emplace_back(chrom, DrawPosition(rng, chromosomes[chrom].second));
   }
   std::sort(placements.begin(), placements.end());

   std::ofstream out(filename);
   WriteHeader(out, /*sort_order=*/"coordinate");
   for (const auto &[chrom, position] : placements) {
      WriteRecord(out, rng, chromosomes[chrom].first, position);
   }
}
