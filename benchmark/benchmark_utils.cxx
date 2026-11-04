#include "benchmark_utils.h"

#include "benchmark_config.h"
#include "generate_sam_benchmark.h"

#include <TROOT.h>

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <system_error>
#include <vector>

#include <fcntl.h>
#include <unistd.h>

StdoutSilencer::StdoutSilencer()
{
   std::cout.flush();
   std::fflush(stdout);
   const int null = open("/dev/null", O_WRONLY | O_CLOEXEC); // NOLINT(cppcoreguidelines-pro-type-vararg)
   if (null < 0)
      return;
   m_saved = dup(STDOUT_FILENO);
   if (m_saved >= 0 && dup2(null, STDOUT_FILENO) < 0) {
      close(m_saved);
      m_saved = -1;
   }
   close(null);
}

StdoutSilencer::~StdoutSilencer()
{
   if (m_saved < 0)
      return;
   std::cout.flush();
   std::fflush(stdout);
   dup2(m_saved, STDOUT_FILENO);
   close(m_saved);
}

WorkDir::WorkDir(const std::string &path)
{
   std::error_code ec;
   if (!path.empty()) {
      std::filesystem::create_directories(path, ec);
      if (ec || !std::filesystem::is_directory(path)) {
         std::cerr << "cannot use --work-dir " << path << ": " << ec.message() << "\n";
         std::exit(1); // NOLINT(concurrency-mt-unsafe): before any benchmark runs
      }
      m_path = path;
      return;
   }
   // create_directory() is false when the name is taken, so a directory that
   // already exists is never used, nor later removed.
   const std::filesystem::path temp = std::filesystem::temp_directory_path(ec);
   std::random_device random;
   for (int attempt = 0; !ec && attempt < 100; ++attempt) {
      const std::filesystem::path candidate = temp / ("ramtools-bench-" + std::to_string(random()));
      if (std::filesystem::create_directory(candidate, ec)) {
         m_path = candidate.string();
         m_owned = true;
         return;
      }
   }
   std::cerr << "cannot create a temporary directory for the benchmark files\n";
   std::exit(1); // NOLINT(concurrency-mt-unsafe): before any benchmark runs
}

WorkDir::~WorkDir()
{
   if (!m_owned)
      return;
   std::error_code ec;
   std::filesystem::remove_all(m_path, ec);
}

std::string WorkDir::File(const std::string &name) const
{
   return (std::filesystem::path(m_path) / name).string();
}

std::vector<SamInput> SamInputs(const BenchmarkOptions &opts, const WorkDir &dir, bool sorted)
{
   if (!opts.sam.empty())
      return {SamInput{"sam", "input", opts.sam, 0, false}};
   std::vector<SamInput> inputs;
   for (const int reads : opts.reads) {
      const std::string stem = "synthetic_" + std::to_string(reads) + (sorted ? "_sorted" : "");
      inputs.push_back(SamInput{"reads:" + std::to_string(reads), stem, dir.File(stem + ".sam"), reads, sorted});
   }
   return inputs;
}

void PrepareSamInput(const SamInput &input)
{
   if (input.reads == 0 || std::filesystem::exists(input.path))
      return;
   if (input.sorted)
      GenerateSortedSAMFile(input.path, input.reads);
   else
      GenerateSAMFile(input.path, input.reads);
}

bool MakeOnce(const std::string &path, const std::function<bool()> &make)
{
   static std::map<std::string, bool> made;
   const auto known = made.find(path);
   if (known != made.end())
      return known->second;
   const bool ok = make();
   made[path] = ok;
   return ok;
}

int64_t CountSamRecords(const std::string &path)
{
   static std::map<std::string, int64_t> counted;
   const auto known = counted.find(path);
   if (known != counted.end())
      return known->second;
   std::ifstream in(path);
   int64_t records = 0;
   std::string line;
   while (std::getline(in, line)) {
      if (!line.empty() && line[0] != '@')
         records++;
   }
   counted[path] = records;
   return records;
}

double FileMiB(const std::string &path)
{
   std::error_code ec;
   const auto size = std::filesystem::file_size(path, ec);
   return ec ? 0.0 : static_cast<double>(size) / (1024.0 * 1024.0);
}

void SetImplicitMT(int threads)
{
   if (ROOT::IsImplicitMTEnabled())
      ROOT::DisableImplicitMT();
   if (threads > 1)
      ROOT::EnableImplicitMT(threads);
}

namespace {

// Collected before the caller acts on them, so that removing files does not
// change the directory while it is being listed.
std::vector<std::filesystem::path> FilesWithPrefix(const std::string &dir, const std::string &prefix)
{
   std::vector<std::filesystem::path> files;
   std::error_code ec;
   for (const auto &entry : std::filesystem::directory_iterator(dir, ec)) {
      if (entry.is_regular_file(ec) && entry.path().filename().string().rfind(prefix, 0) == 0)
         files.push_back(entry.path());
   }
   return files;
}

} // namespace

uint64_t TotalSizeWithPrefix(const std::string &dir, const std::string &prefix)
{
   uint64_t total = 0;
   for (const auto &file : FilesWithPrefix(dir, prefix)) {
      std::error_code ec;
      const auto size = std::filesystem::file_size(file, ec);
      if (!ec)
         total += size;
   }
   return total;
}

void RemoveWithPrefix(const std::string &dir, const std::string &prefix)
{
   for (const auto &file : FilesWithPrefix(dir, prefix)) {
      std::error_code ec;
      std::filesystem::remove(file, ec);
   }
}
