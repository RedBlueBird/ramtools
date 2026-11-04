#include "bam_ops.h"

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <memory>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

namespace {

struct SamFileCloser {
   void operator()(samFile *f) const { sam_close(f); }
};
struct HeaderDestroyer {
   void operator()(sam_hdr_t *h) const { sam_hdr_destroy(h); }
};
struct IndexDestroyer {
   void operator()(hts_idx_t *i) const { hts_idx_destroy(i); }
};
struct IteratorDestroyer {
   void operator()(hts_itr_t *i) const { hts_itr_destroy(i); }
};
struct RecordDestroyer {
   void operator()(bam1_t *b) const { bam_destroy1(b); }
};

using SamFile = std::unique_ptr<samFile, SamFileCloser>;
using Header = std::unique_ptr<sam_hdr_t, HeaderDestroyer>;
using Index = std::unique_ptr<hts_idx_t, IndexDestroyer>;
using Iterator = std::unique_ptr<hts_itr_t, IteratorDestroyer>;
using Record = std::unique_ptr<bam1_t, RecordDestroyer>;

bool Fail(const std::string &what)
{
   std::cerr << "bam_ops: " << what << "\n";
   return false;
}

bool HeaderSaysCoordinateSorted(sam_hdr_t *hdr)
{
   const std::string_view text = sam_hdr_str(hdr) ? sam_hdr_str(hdr) : "";
   if (text.rfind("@HD", 0) != 0)
      return false;
   return text.substr(0, text.find('\n')).find("\tSO:coordinate") != std::string_view::npos;
}

// An open BAM with its header and index, one per thread that reads it.
struct IndexedBam {
   SamFile file;
   Header hdr;
   Index idx;

   explicit IndexedBam(const std::string &path) : file(sam_open(path.c_str(), "r"))
   {
      if (file)
         hdr.reset(sam_hdr_read(file.get()));
      if (hdr)
         idx.reset(sam_index_load(file.get(), path.c_str()));
   }
   [[nodiscard]] bool Ok() const { return file && hdr && idx; }
};

} // namespace

bool SamToSortedBam(const std::string &input, const std::string &bam, int threads)
{
   SamFile in(sam_open(input.c_str(), "r"));
   if (!in)
      return Fail("cannot open " + input);
   Header hdr(sam_hdr_read(in.get()));
   if (!hdr)
      return Fail("cannot read the header of " + input);
   if (threads > 1)
      hts_set_threads(in.get(), threads);

   // A sorted input streams through; any other is read into memory and sorted,
   // since only a coordinate-sorted BAM can be indexed.
   const bool sorted = HeaderSaysCoordinateSorted(hdr.get());
   std::vector<Record> records;
   Record rec(bam_init1());
   if (!sorted) {
      int read = 0;
      while ((read = sam_read1(in.get(), hdr.get(), rec.get())) >= 0) {
         if (static_cast<long>(records.size()) >= kMaxSortRecords)
            return Fail(input + " is not coordinate sorted and too large to sort here; sort it first");
         records.emplace_back(bam_dup1(rec.get()));
      }
      if (read < -1)
         return Fail("cannot read " + input);
      // Records with no reference (tid -1) go last, as in a sorted BAM; equal
      // positions keep their input order.
      std::stable_sort(records.begin(), records.end(), [](const Record &a, const Record &b) {
         const auto a_tid = static_cast<uint32_t>(a->core.tid);
         const auto b_tid = static_cast<uint32_t>(b->core.tid);
         return a_tid != b_tid ? a_tid < b_tid : a->core.pos < b->core.pos;
      });
      // A header may have no @HD line to update.
      if (sam_hdr_update_hd(hdr.get(), "SO", "coordinate") < 0 &&
          sam_hdr_add_lines(hdr.get(), "@HD\tVN:1.6\tSO:coordinate", /*len=*/0) < 0)
         return Fail("cannot mark the header of " + bam + " as sorted");
   }

   {
      SamFile out(sam_open(bam.c_str(), "wb"));
      if (!out)
         return Fail("cannot create " + bam);
      if (threads > 1)
         hts_set_threads(out.get(), threads);
      if (sam_hdr_write(out.get(), hdr.get()) < 0)
         return Fail("cannot write " + bam);
      if (sorted) {
         int read = 0;
         while ((read = sam_read1(in.get(), hdr.get(), rec.get())) >= 0) {
            if (sam_write1(out.get(), hdr.get(), rec.get()) < 0)
               return Fail("cannot write " + bam);
         }
         if (read < -1)
            return Fail("cannot read " + input);
      }
      for (const Record &r : records) {
         if (sam_write1(out.get(), hdr.get(), r.get()) < 0)
            return Fail("cannot write " + bam);
      }
   } // Closed, and so complete, before it is indexed.

   if (sam_index_build3(bam.c_str(), nullptr, /*min_shift=*/0, threads) < 0)
      return Fail("cannot index " + bam);
   return true;
}

bool SplitBam(const std::string &bam, const std::string &prefix, int threads)
{
   std::vector<int> tids;
   {
      const IndexedBam in(bam);
      if (!in.Ok())
         return Fail("cannot open " + bam + " with its index");
      for (int tid = 0; tid < sam_hdr_nref(in.hdr.get()); ++tid) {
         uint64_t mapped = 0;
         uint64_t unmapped = 0;
         if (hts_idx_get_stat(in.idx.get(), tid, &mapped, &unmapped) == 0 && mapped + unmapped > 0)
            tids.push_back(tid);
      }
   }

   // Each worker opens its own handles and takes every workers-th reference.
   const auto workers = static_cast<size_t>(std::clamp(threads, 1, std::max(1, static_cast<int>(tids.size()))));
   std::atomic<bool> ok = true;
   auto extract = [&](size_t first) {
      const IndexedBam in(bam);
      Record rec(bam_init1());
      if (!in.Ok() || !rec) {
         ok = false;
         return;
      }
      for (size_t i = first; i < tids.size() && ok; i += workers) {
         const std::string name = sam_hdr_tid2name(in.hdr.get(), tids[i]);
         const std::string path = prefix + name + ".bam";
         SamFile out(sam_open(path.c_str(), "wb"));
         Iterator itr(sam_itr_queryi(in.idx.get(), tids[i], 0, HTS_POS_MAX));
         if (!out || !itr || sam_hdr_write(out.get(), in.hdr.get()) < 0) {
            ok = false;
            return;
         }
         int read = 0;
         while ((read = sam_itr_next(in.file.get(), itr.get(), rec.get())) >= 0) {
            if (sam_write1(out.get(), in.hdr.get(), rec.get()) < 0) {
               ok = false;
               return;
            }
         }
         if (read < -1)
            ok = false;
      }
   };
   std::vector<std::thread> pool;
   for (size_t w = 1; w < workers; ++w)
      pool.emplace_back(extract, w);
   extract(0);
   for (auto &t : pool)
      t.join();
   return ok ? true : Fail("cannot split " + bam);
}
