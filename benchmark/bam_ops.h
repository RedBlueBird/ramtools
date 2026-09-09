#ifndef BENCHMARK_BAM_OPS_H
#define BENCHMARK_BAM_OPS_H

// In-process BAM operations for the benchmark suite.
//
// htslib is already a hard dependency of ramcore (linked PUBLIC by the top-level
// CMakeLists), so every benchmark binary picks it up transitively -- including this
// header needs no extra CMake wiring.
//
// These call htslib directly rather than shelling out to `samtools`, which puts BAM on
// the same footing as the ROOT backends: no fork/exec inside a timed loop, and index
// loading can be hoisted out of the measured region. It also drops the suite's runtime
// dependency on a `samtools` binary being installed.
//
// Fidelity: BAM has no lossy mode. An RNTuple written with kIlluminaBinning or kDrop
// does not hold the same data as the BAM it is being compared against, so pass
// --quality=0 for any number meant to be read as a format comparison.

#include "benchmark_utils.h"

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

namespace benchutil {

// Flags a region query skips, matching RAMNTupleView's filtering: unmapped (0x4),
// secondary (0x100), supplementary (0x800). Without this the htslib iterator returns
// records the ROOT views drop, and the two record counts stop being comparable.
constexpr uint16_t kQuerySkipFlags = 0x4 | 0x100 | 0x800;

// Refuse to sort in memory beyond this many records. The synthetic ladder tops out at
// 1M reads (~200 MB of bam1_t), which is fine; a large unsorted real dataset is not, and
// should be pre-sorted and passed via --bam instead.
constexpr std::size_t kMaxInMemorySortRecords = 5000000;

namespace detail {

struct HtsFileDeleter {
   void operator()(samFile *f) const noexcept
   {
      if (f != nullptr)
         sam_close(f);
   }
};
struct HtsHdrDeleter {
   void operator()(sam_hdr_t *h) const noexcept
   {
      if (h != nullptr)
         sam_hdr_destroy(h);
   }
};
struct HtsIdxDeleter {
   void operator()(hts_idx_t *i) const noexcept
   {
      if (i != nullptr)
         hts_idx_destroy(i);
   }
};
struct HtsItrDeleter {
   void operator()(hts_itr_t *i) const noexcept
   {
      if (i != nullptr)
         hts_itr_destroy(i);
   }
};
struct Bam1Deleter {
   void operator()(bam1_t *b) const noexcept
   {
      if (b != nullptr)
         bam_destroy1(b);
   }
};

using HtsFilePtr = std::unique_ptr<samFile, HtsFileDeleter>;
using HtsHdrPtr = std::unique_ptr<sam_hdr_t, HtsHdrDeleter>;
using HtsIdxPtr = std::unique_ptr<hts_idx_t, HtsIdxDeleter>;
using HtsItrPtr = std::unique_ptr<hts_itr_t, HtsItrDeleter>;
using Bam1Ptr = std::unique_ptr<bam1_t, Bam1Deleter>;

// True when the header's @HD line declares SO:coordinate. Indexing -- and therefore any
// region query or indexed split -- requires it. GenerateSAMFile() writes SO:unsorted, so
// synthetic runs take the false branch and get sorted first.
inline bool IsCoordinateSorted(sam_hdr_t *hdr)
{
   const char *text = sam_hdr_str(hdr);
   if (text == nullptr)
      return false;
   const std::string_view all(text);
   if (all.rfind("@HD", 0) != 0)
      return false;
   return all.substr(0, all.find('\n')).find("SO:coordinate") != std::string_view::npos;
}

// Rewrite the @HD SO: value, so a file we sorted ourselves advertises it and htslib will
// index it. Returns false if the header has no @HD line to update.
inline bool MarkCoordinateSorted(sam_hdr_t *hdr)
{
   if (sam_hdr_update_line(hdr, "HD", nullptr, nullptr, "SO", "coordinate", nullptr) == 0)
      return true;
   return sam_hdr_add_line(hdr, "HD", "VN", "1.6", "SO", "coordinate", nullptr) == 0;
}

} // namespace detail

// True when `path` (SAM or BAM -- htslib reads both) declares @HD SO:coordinate.
//
// Worth checking before a region query: ramview() and ramntupleview() seek to an
// index-derived start row and then stop scanning at the first record past the region's
// end, which only reaches every overlapping record if the file is coordinate sorted.
// GenerateSAMFile() writes SO:unsorted, so on synthetic data those two backends
// under-report and their counts are not comparable with BAM's.
inline bool FileIsCoordinateSorted(const std::string &path)
{
   detail::HtsFilePtr in(sam_open(path.c_str(), "r"));
   if (!in)
      return false;
   detail::HtsHdrPtr hdr(sam_hdr_read(in.get()));
   return hdr && detail::IsCoordinateSorted(hdr.get());
}

// Outcome of BamConvert: total bytes written (BAM plus its .bai), whether an index was
// produced, and whether the input had to be sorted on the way through.
struct BamConvertResult {
   std::size_t bytes = 0;
   bool indexed = false;
   bool sorted = false; // true => a sort happened inside this call
   bool ok = false;
};

// SAM/BAM -> a *queryable* BAM: coordinate sorted, written, and indexed.
//
// An already-sorted input streams straight through and pays nothing extra. An unsorted
// one (the synthetic generator) is sorted in memory first -- that cost is inside the
// call, and is reported via BamConvertResult::sorted, because a BAM that is not
// coordinate sorted cannot be indexed and so cannot answer a region query at all. The
// returned byte count includes the .bai, making it comparable with an RNTuple's
// self-contained file.
inline BamConvertResult BamConvert(const std::string &input, const std::string &outputBam, int threads)
{
   using namespace detail;
   BamConvertResult result;

   HtsFilePtr in(sam_open(input.c_str(), "r"));
   if (!in)
      return result;
   HtsHdrPtr hdr(sam_hdr_read(in.get()));
   if (!hdr)
      return result;
   if (threads > 1)
      hts_set_threads(in.get(), threads);

   const bool preSorted = IsCoordinateSorted(hdr.get());

   // Records are buffered only when a sort is actually required; the sorted path streams.
   std::vector<Bam1Ptr> buffered;
   {
      HtsFilePtr out(sam_open(outputBam.c_str(), "wb"));
      if (!out)
         return result;
      if (threads > 1)
         hts_set_threads(out.get(), threads);

      if (preSorted) {
         if (sam_hdr_write(out.get(), hdr.get()) < 0)
            return result;
         Bam1Ptr rec(bam_init1());
         while (sam_read1(in.get(), hdr.get(), rec.get()) >= 0) {
            if (sam_write1(out.get(), hdr.get(), rec.get()) < 0)
               return result;
         }
      } else {
         Bam1Ptr rec(bam_init1());
         while (sam_read1(in.get(), hdr.get(), rec.get()) >= 0) {
            if (buffered.size() >= kMaxInMemorySortRecords)
               return result; // too large to sort here: pre-sort and pass --bam
            buffered.emplace_back(bam_dup1(rec.get()));
         }
         std::stable_sort(buffered.begin(), buffered.end(), [](const Bam1Ptr &a, const Bam1Ptr &b) {
            // Unplaced reads (tid < 0) sort last, as in a coordinate-sorted BAM.
            const auto atid = static_cast<uint32_t>(a->core.tid);
            const auto btid = static_cast<uint32_t>(b->core.tid);
            if (atid != btid)
               return atid < btid;
            return a->core.pos < b->core.pos;
         });
         result.sorted = true;
         MarkCoordinateSorted(hdr.get());
         if (sam_hdr_write(out.get(), hdr.get()) < 0)
            return result;
         for (const auto &b : buffered) {
            if (sam_write1(out.get(), hdr.get(), b.get()) < 0)
               return result;
         }
      }
   } // out closed (and flushed) here -- required before the file can be indexed

   result.bytes = FileSize(outputBam);
   if (sam_index_build3(outputBam.c_str(), nullptr, /*min_shift=*/0, threads) == 0) {
      result.indexed = true;
      result.bytes += FileSize(outputBam + ".bai");
   }
   result.ok = true;
   return result;
}

// An open BAM plus its index, so that opening and index loading stay outside the timed
// region and only Query() is measured -- the same treatment the ROOT views get.
class BamReader {
public:
   explicit BamReader(const std::string &path) : m_file(sam_open(path.c_str(), "r"))
   {
      if (!m_file)
         return;
      m_hdr.reset(sam_hdr_read(m_file.get()));
      if (!m_hdr)
         return;
      m_idx.reset(sam_index_load(m_file.get(), path.c_str()));
      m_rec.reset(bam_init1());
   }

   [[nodiscard]] bool Ok() const { return m_file && m_hdr && m_idx && m_rec; }

   // Records overlapping `region` ("name:start-end"), after the kQuerySkipFlags filter.
   // Returns -1 if the reader is unusable, 0 for a region naming an unknown reference
   // (which is what the ROOT views report for the same input).
   [[nodiscard]] int64_t Query(const std::string &region) const
   {
      if (!Ok())
         return -1;
      detail::HtsItrPtr itr(sam_itr_querys(m_idx.get(), m_hdr.get(), region.c_str()));
      if (!itr)
         return 0;
      int64_t matched = 0;
      while (sam_itr_next(m_file.get(), itr.get(), m_rec.get()) >= 0) {
         if ((m_rec->core.flag & kQuerySkipFlags) == 0)
            ++matched;
      }
      return matched;
   }

private:
   detail::HtsFilePtr m_file;
   detail::HtsHdrPtr m_hdr;
   detail::HtsIdxPtr m_idx;
   detail::Bam1Ptr m_rec;
};

// Split a coordinate-sorted, indexed BAM into one file per populated chromosome, named
// "<prefix>_<chr>.bam". Only references that actually carry records get a file, matching
// samtoramntuple_split_by_chromosome, which writes a .root only for chromosomes it saw.
// `threads` worker threads each own their own handles and take a strided slice of the
// reference list -- the analogue of the RNTuple side's parallel per-chromosome writers.
// Returns the total bytes written (0 on failure).
inline std::size_t BamSplit(const std::string &sortedBam, const std::string &prefix, int threads)
{
   using namespace detail;

   std::vector<std::string> refs;
   {
      HtsFilePtr in(sam_open(sortedBam.c_str(), "r"));
      if (!in)
         return 0;
      HtsHdrPtr hdr(sam_hdr_read(in.get()));
      if (!hdr)
         return 0;
      HtsIdxPtr idx(sam_index_load(in.get(), sortedBam.c_str()));
      if (!idx)
         return 0;
      const int nref = sam_hdr_nref(hdr.get());
      for (int tid = 0; tid < nref; ++tid) {
         uint64_t mapped = 0;
         uint64_t unmapped = 0;
         if (hts_idx_get_stat(idx.get(), tid, &mapped, &unmapped) == 0 && (mapped + unmapped) > 0)
            refs.emplace_back(sam_hdr_tid2name(hdr.get(), tid));
      }
   }
   if (refs.empty())
      return 0;

   const auto workers = static_cast<std::size_t>(std::max(1, std::min<int>(threads, static_cast<int>(refs.size()))));
   std::vector<std::thread> pool;
   pool.reserve(workers);
   for (std::size_t w = 0; w < workers; ++w) {
      pool.emplace_back([&refs, &sortedBam, &prefix, w, workers]() {
         HtsFilePtr in(sam_open(sortedBam.c_str(), "r"));
         if (!in)
            return;
         HtsHdrPtr hdr(sam_hdr_read(in.get()));
         if (!hdr)
            return;
         HtsIdxPtr idx(sam_index_load(in.get(), sortedBam.c_str()));
         if (!idx)
            return;
         Bam1Ptr rec(bam_init1());
         for (std::size_t i = w; i < refs.size(); i += workers) {
            const std::string &chr = refs[i];
            HtsFilePtr out(sam_open((prefix + "_" + chr + ".bam").c_str(), "wb"));
            if (!out)
               continue;
            if (sam_hdr_write(out.get(), hdr.get()) < 0)
               continue;
            HtsItrPtr itr(sam_itr_querys(idx.get(), hdr.get(), chr.c_str()));
            if (!itr)
               continue;
            while (sam_itr_next(in.get(), itr.get(), rec.get()) >= 0) {
               if (sam_write1(out.get(), hdr.get(), rec.get()) < 0)
                  break;
            }
         }
      });
   }
   for (auto &t : pool)
      t.join();

   return TotalSizeByPrefix(prefix + "_");
}

} // namespace benchutil

#endif
