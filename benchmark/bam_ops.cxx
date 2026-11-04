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

int64_t FailCount(const std::string &what)
{
   Fail(what);
   return -1;
}

bool IsCram(const std::string &path)
{
   const std::string_view ext = ".cram";
   return path.size() >= ext.size() && path.compare(path.size() - ext.size(), ext.size(), ext) == 0;
}

// Opens a file for reading. A CRAM is decoded against ref, and without one it
// is refused rather than have htslib fetch the reference over the network.
SamFile OpenForReading(const std::string &path, const std::string &ref, int threads)
{
   SamFile file(sam_open(path.c_str(), "r"));
   if (!file) {
      Fail("cannot open " + path);
      return nullptr;
   }
   if (hts_get_format(file.get())->format == cram) {
      if (ref.empty() || hts_set_fai_filename(file.get(), ref.c_str()) < 0) {
         Fail("reading the CRAM " + path + " needs its reference: give --ref");
         return nullptr;
      }
   }
   if (threads > 1)
      hts_set_threads(file.get(), threads);
   return file;
}

SamFile OpenForWriting(const std::string &path, const std::string &ref, int threads)
{
   const bool cram = IsCram(path);
   if (cram && ref.empty()) {
      Fail("writing the CRAM " + path + " needs a reference: give --ref");
      return nullptr;
   }
   SamFile file(sam_open(path.c_str(), cram ? "wc" : "wb"));
   if (!file) {
      Fail("cannot create " + path);
      return nullptr;
   }
   if (cram && hts_set_fai_filename(file.get(), ref.c_str()) < 0) {
      Fail("cannot use the reference " + ref);
      return nullptr;
   }
   if (threads > 1)
      hts_set_threads(file.get(), threads);
   return file;
}

bool HeaderSaysCoordinateSorted(sam_hdr_t *hdr)
{
   const std::string_view text = sam_hdr_str(hdr) ? sam_hdr_str(hdr) : "";
   if (text.rfind("@HD", 0) != 0)
      return false;
   return text.substr(0, text.find('\n')).find("\tSO:coordinate") != std::string_view::npos;
}

// The records of in, read to the end, in coordinate order. Records with no
// reference (tid -1) go last, as in a sorted file; equal positions keep their
// input order.
bool ReadSorted(samFile *in, sam_hdr_t *hdr, const std::string &input, std::vector<Record> &records)
{
   Record rec(bam_init1());
   int read = 0;
   while ((read = sam_read1(in, hdr, rec.get())) >= 0) {
      if (static_cast<long>(records.size()) >= kMaxSortRecords)
         return Fail(input + " is not coordinate sorted and too large to sort here; sort it first");
      records.emplace_back(bam_dup1(rec.get()));
   }
   if (read < -1)
      return Fail("cannot read " + input);
   std::stable_sort(records.begin(), records.end(), [](const Record &a, const Record &b) {
      const auto a_tid = static_cast<uint32_t>(a->core.tid);
      const auto b_tid = static_cast<uint32_t>(b->core.tid);
      return a_tid != b_tid ? a_tid < b_tid : a->core.pos < b->core.pos;
   });
   // A header may have no @HD line to update.
   if (sam_hdr_update_hd(hdr, "SO", "coordinate") < 0 &&
       sam_hdr_add_lines(hdr, "@HD\tVN:1.6\tSO:coordinate", /*len=*/0) < 0)
      return Fail("cannot mark the header of " + input + " as sorted");
   return true;
}

bool Convert(const std::string &input, const std::string &output, int threads, const std::string &ref, bool sort)
{
   const SamFile in = OpenForReading(input, ref, threads);
   if (!in)
      return false;
   const Header hdr(sam_hdr_read(in.get()));
   if (!hdr)
      return Fail("cannot read the header of " + input);

   // A sorted input streams through, as does any input when no order is asked
   // for; else every record is read into memory and sorted first.
   const bool stream = !sort || HeaderSaysCoordinateSorted(hdr.get());
   std::vector<Record> records;
   if (!stream && !ReadSorted(in.get(), hdr.get(), input, records))
      return false;

   {
      const SamFile out = OpenForWriting(output, ref, threads);
      if (!out)
         return false;
      if (sam_hdr_write(out.get(), hdr.get()) < 0)
         return Fail("cannot write " + output);
      if (stream) {
         const Record rec(bam_init1());
         int read = 0;
         while ((read = sam_read1(in.get(), hdr.get(), rec.get())) >= 0) {
            if (sam_write1(out.get(), hdr.get(), rec.get()) < 0)
               return Fail("cannot write " + output);
         }
         if (read < -1)
            return Fail("cannot read " + input);
      }
      for (const Record &r : records) {
         if (sam_write1(out.get(), hdr.get(), r.get()) < 0)
            return Fail("cannot write " + output);
      }
   } // Closed, and so complete, before it is indexed.

   if (sort && sam_index_build3(output.c_str(), nullptr, /*min_shift=*/0, threads) < 0)
      return Fail("cannot index " + output);
   return true;
}

// An open BAM or CRAM with its header and index, one per thread that reads it.
struct IndexedFile {
   SamFile file;
   Header hdr;
   Index idx;

   IndexedFile(const std::string &path, const std::string &ref) : file(OpenForReading(path, ref, /*threads=*/1))
   {
      if (file)
         hdr.reset(sam_hdr_read(file.get()));
      if (hdr)
         idx.reset(sam_index_load(file.get(), path.c_str()));
   }
   [[nodiscard]] bool Ok() const { return file && hdr && idx; }
};

} // namespace

bool ConvertAlignments(const std::string &input, const std::string &output, int threads, const std::string &ref)
{
   return Convert(input, output, threads, ref, /*sort=*/false);
}

bool ConvertSorted(const std::string &input, const std::string &output, int threads, const std::string &ref)
{
   return Convert(input, output, threads, ref, /*sort=*/true);
}

bool SplitBam(const std::string &bam, const std::string &prefix, int threads)
{
   std::vector<int> tids;
   {
      const IndexedFile in(bam, /*ref=*/"");
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
      const IndexedFile in(bam, /*ref=*/"");
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

int64_t CountRegion(const std::string &path, const std::string &region, const std::string &ref)
{
   const IndexedFile in(path, ref);
   if (!in.Ok())
      return FailCount("cannot open " + path + " with its index");
   const Iterator itr(sam_itr_querys(in.idx.get(), in.hdr.get(), region.c_str()));
   if (!itr)
      return FailCount("cannot query " + region + " in " + path);
   const Record rec(bam_init1());
   int64_t count = 0;
   int read = 0;
   while ((read = sam_itr_next(in.file.get(), itr.get(), rec.get())) >= 0)
      count++;
   return read < -1 ? FailCount("cannot read " + path) : count;
}

int64_t ScanRegion(const std::string &path, const std::string &region)
{
   const SamFile in = OpenForReading(path, /*ref=*/"", /*threads=*/1);
   if (!in)
      return -1;
   const Header hdr(sam_hdr_read(in.get()));
   if (!hdr)
      return FailCount("cannot read the header of " + path);
   int tid = -1;
   hts_pos_t beg = 0;
   hts_pos_t end = 0;
   if (!sam_parse_region(hdr.get(), region.c_str(), &tid, &beg, &end, /*flags=*/0) || tid < 0)
      return FailCount("cannot parse the region " + region + " for " + path);

   const Record rec(bam_init1());
   int64_t count = 0;
   int read = 0;
   while ((read = sam_read1(in.get(), hdr.get(), rec.get())) >= 0) {
      if (rec->core.tid == tid && rec->core.pos < end && bam_endpos(rec.get()) > beg)
         count++;
   }
   return read < -1 ? FailCount("cannot read " + path) : count;
}

bool ScanRecords(const std::string &path, const std::string &ref, bool flagstat_fields, ScanTotals &totals)
{
   const SamFile in = OpenForReading(path, ref, /*threads=*/1);
   if (!in)
      return false;
   if (flagstat_fields && hts_get_format(in.get())->format == cram &&
       // NOLINTNEXTLINE(cppcoreguidelines-pro-type-vararg): htslib's option setter
       hts_set_opt(in.get(), CRAM_OPT_REQUIRED_FIELDS, SAM_FLAG | SAM_RNAME | SAM_RNEXT | SAM_MAPQ) < 0)
      return Fail("cannot limit the fields decoded from " + path);
   const Header hdr(sam_hdr_read(in.get()));
   if (!hdr)
      return Fail("cannot read the header of " + path);
   const Record rec(bam_init1());
   totals = ScanTotals{};
   int read = 0;
   while ((read = sam_read1(in.get(), hdr.get(), rec.get())) >= 0) {
      totals.records++;
      if (rec->core.tid >= 0)
         totals.placed++;
      totals.flag_sum += rec->core.flag;
   }
   return read < -1 ? Fail("cannot read " + path) : true;
}
