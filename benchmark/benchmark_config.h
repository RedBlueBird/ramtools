#ifndef BENCHMARK_CONFIG_H
#define BENCHMARK_CONFIG_H

// Shared command-line configuration for the RAMtools benchmark suite.
//
// Lets every benchmark binary run on a user-supplied dataset instead of only on
// synthetic data or hardcoded paths. Recognized flags (all optional):
//
//   --sam=PATH            SAM file to benchmark (default: generate synthetic)
//   --bam=PATH            pre-built BAM (default: derive from --sam).
//   --ttree-root=PATH     pre-converted TTree .root (default: derive from --sam)
//   --rntuple-root=PATH   pre-converted RNTuple .root (default: derive from --sam)
//   --compression=N       ROOT compression setting for RNTuple (default 505 = ZSTD-5)
//   --quality=N           quality policy bitmask (default 0)
//   --threads=N           worker threads for split benchmarks (default 4)
//   --reads=N             synthetic read count when no --sam given (default 100000)
//   --regions=all|i,j,k   region indices for region_query (default: 0,3,6,9)
//   --regions-file=PATH   region_query: load "name:start-end" regions from this file
//                         (one per line, blank lines and lines starting with '#'
//                         ignored) instead of the built-in synthetic-data region
//                         list. Use this when querying a real dataset.
//
// Parse with BenchmarkConfig::FromArgs(&argc, argv) BEFORE benchmark::Initialize so
// these flags are stripped and only --benchmark_* flags reach Google Benchmark.
// When no dataset path is supplied, EnsureSam()/EnsureTTreeRoot()/EnsureRNTupleRoot()/
// EnsureBam() fall back to GenerateSAMFile() and convert on demand, cleaning up temp
// files on destruction.
//
// Note on --quality: it selects an RNTuple quality policy (0 = lossless Phred33). BAM
// has no lossy mode, so only --quality=0 gives an RNTuple-vs-BAM number that compares
// like with like.

#include "bam_ops.h"
#include "benchmark_utils.h"
#include "generate_sam_benchmark.h"
#include "ramcore/SamToNTuple.h"
#include "ramcore/SamToTTree.h"

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace benchutil {

struct BenchmarkConfig {
public:
   // NOLINTBEGIN(cppcoreguidelines-non-private-member-variables-in-classes)
   std::string sam;         // empty => generate synthetic
   std::string bam;         // empty => derive by converting sam
   std::string ttreeRoot;   // empty => derive by converting sam
   std::string rntupleRoot; // empty => derive by converting sam
   int compression = 505;   // ROOT compression setting for RNTuple (ZSTD-5)
   unsigned int quality = 0;
   int threads = 4;
   int reads = 100000;       // synthetic read count for fallback
   std::vector<int> regions; // explicit region indices for region_query
   bool allRegions = false;  // --regions=all
   std::string regionsFile;  // --regions-file=PATH: real "name:start-end" regions
   // NOLINTEND(cppcoreguidelines-non-private-member-variables-in-classes)

   // Parse and strip recognized "--key=value" flags from argv (compacting it in place
   // and lowering *argc), leaving --benchmark_* and other flags for the caller.
   static BenchmarkConfig FromArgs(int *argc, char **argv)
   {
      BenchmarkConfig cfg;
      int out = 1;
      for (int i = 1; i < *argc; ++i) {
         if (!cfg.Consume(argv[i])) // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
            argv[out++] = argv[i];  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      }
      *argc = out;
      cfg.RequireInputsExist();
      return cfg;
   }

   // Abort on a dataset path that does not exist.
   //
   // Without this a mistyped -- or relative-to-the-wrong-directory -- path still "runs":
   // samtoram()/samtoramntuple() accept the missing input, write an essentially empty
   // file, and the suite reports plausible microsecond timings over a couple of KB of
   // nothing. Failing loudly here is the only thing that distinguishes that from a real
   // result. Note that the binaries run from the build tree, so relative paths resolve
   // against it; scripts/run_benchmarks.sh rewrites them to absolute for you.
   void RequireInputsExist() const
   {
      bool missing = false;
      auto check = [&missing](const char *flag, const std::string &path) {
         if (path.empty() || std::filesystem::exists(path))
            return;
         std::cerr << "benchmark_config: " << flag << "=" << path << " does not exist\n";
         missing = true;
      };
      check("--sam", sam);
      check("--bam", bam);
      check("--ttree-root", ttreeRoot);
      check("--rntuple-root", rntupleRoot);
      check("--regions-file", regionsFile);
      if (missing) {
         std::cerr << "  (working directory: " << std::filesystem::current_path().string() << ")\n";
         std::exit(1);
      }
   }

   // Path to a SAM file: the user-supplied one, or a freshly generated synthetic file.
   const std::string &EnsureSam()
   {
      if (!sam.empty())
         return sam;
      if (m_genSam.empty()) {
         m_genSam = "bench_gen_" + std::to_string(reads) + ".sam";
         GenerateSAMFile(m_genSam, reads);
         m_created.push_back(m_genSam);
      }
      return m_genSam;
   }

   // Path to a queryable BAM: the user-supplied one, or one converted from the SAM.
   //
   // The result is coordinate sorted with a .bai beside it, since an unsorted BAM cannot
   // be indexed and so cannot answer a region query. GenerateSAMFile() emits SO:unsorted,
   // so the synthetic path sorts (in memory) here -- outside any timed loop.
   const std::string &EnsureBam()
   {
      if (!bam.empty())
         return bam;
      if (m_genBam.empty()) {
         const std::string out = "bench_gen.bam";
         const std::string &s = EnsureSam();
         const auto res = benchutil::BamConvert(s, out, threads);
         // Register for cleanup even on failure: a partial file may still be on disk.
         m_created.push_back(out);
         m_created.push_back(out + ".bai");
         if (!res.ok || !res.indexed) {
            std::cerr << "benchmark_config: could not build a sorted, indexed BAM from " << s
                      << " (pass --bam=PATH with a coordinate-sorted, indexed file)\n";
            return m_genBam; // still empty => caller skips its BAM cases
         }
         m_genBam = out;
      }
      return m_genBam;
   }

   // Path to an RNTuple .root: the user-supplied one, or one converted from the SAM.
   const std::string &EnsureRNTupleRoot()
   {
      if (!rntupleRoot.empty())
         return rntupleRoot;
      if (m_genRntuple.empty()) {
         m_genRntuple = "bench_gen_rntuple.root";
         const std::string &s = EnsureSam();
         ScopedStdoutSuppressor quiet(/*suppress_stderr=*/true);
         samtoramntuple(s.c_str(), m_genRntuple.c_str(), /*index=*/true, /*split=*/true, /*cache=*/true, compression,
                        quality);
         m_created.push_back(m_genRntuple);
      }
      return m_genRntuple;
   }

   // Path to a TTree .root: the user-supplied one, or one converted from the SAM.
   const std::string &EnsureTTreeRoot()
   {
      if (!ttreeRoot.empty())
         return ttreeRoot;
      if (m_genTtree.empty()) {
         m_genTtree = "bench_gen_ttree.root";
         const std::string &s = EnsureSam();
         ScopedStdoutSuppressor quiet(/*suppress_stderr=*/true);
         samtoram(s.c_str(), m_genTtree.c_str(), /*index=*/true, /*split=*/true, /*cache=*/true,
                  /*compression_algorithm=*/1, quality);
         m_created.push_back(m_genTtree);
      }
      return m_genTtree;
   }

   // True when the user supplied a real SAM (vs. falling back to synthetic data).
   [[nodiscard]] bool HasRealDataset() const { return !sam.empty(); }

   void Cleanup()
   {
      for (const auto &f : m_created)
         std::remove(f.c_str());
      m_created.clear();
   }

   BenchmarkConfig() = default;
   BenchmarkConfig(const BenchmarkConfig &) = delete;
   BenchmarkConfig &operator=(const BenchmarkConfig &) = delete;
   BenchmarkConfig(BenchmarkConfig &&) = default;
   BenchmarkConfig &operator=(BenchmarkConfig &&) = default;
   ~BenchmarkConfig() { Cleanup(); }

private:
   std::vector<std::string> m_created;
   std::string m_genSam;
   std::string m_genBam;
   std::string m_genRntuple;
   std::string m_genTtree;

   static bool Match(const std::string &arg, const char *key, std::string &val)
   {
      const std::string prefix = std::string("--") + key + "=";
      if (arg.rfind(prefix, 0) == 0) {
         val = arg.substr(prefix.size());
         return true;
      }
      return false;
   }

   bool Consume(const std::string &arg)
   {
      std::string v;
      if (Match(arg, /*key=*/"sam", v)) {
         sam = v;
         return true;
      }
      if (Match(arg, /*key=*/"bam", v)) {
         bam = v;
         return true;
      }
      if (Match(arg, /*key=*/"ttree-root", v)) {
         ttreeRoot = v;
         return true;
      }
      if (Match(arg, /*key=*/"rntuple-root", v)) {
         rntupleRoot = v;
         return true;
      }
      if (Match(arg, /*key=*/"compression", v)) {
         compression = std::atoi(v.c_str());
         return true;
      }
      if (Match(arg, /*key=*/"quality", v)) {
         quality = static_cast<unsigned int>(std::strtoul(v.c_str(), nullptr, 10));
         return true;
      }
      if (Match(arg, /*key=*/"threads", v)) {
         threads = std::atoi(v.c_str());
         return true;
      }
      if (Match(arg, /*key=*/"reads", v)) {
         reads = std::atoi(v.c_str());
         return true;
      }
      if (Match(arg, /*key=*/"regions", v)) {
         ParseRegions(v);
         return true;
      }
      if (Match(arg, /*key=*/"regions-file", v)) {
         regionsFile = v;
         return true;
      }
      return false;
   }

   void ParseRegions(const std::string &v)
   {
      regions.clear();
      allRegions = false;
      if (v == "all") {
         allRegions = true;
         return;
      }
      std::size_t start = 0;
      while (start <= v.size()) {
         const std::size_t comma = v.find(',', start);
         const std::string tok = v.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
         if (!tok.empty())
            regions.push_back(std::atoi(tok.c_str()));
         if (comma == std::string::npos)
            break;
         start = comma + 1;
      }
   }
};

} // namespace benchutil

#endif
