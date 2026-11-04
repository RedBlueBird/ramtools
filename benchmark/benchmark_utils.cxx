#include "benchmark_utils.h"

#include "benchmark_config.h"
#include "generate_sam_benchmark.h"
#include "ramcore/SamToNTuple.h"
#include "rntuple/RAMFile.h"
#include "rntuple/RAMNTupleRecord.h"

#include <TROOT.h>
#include <benchmark/benchmark.h>

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
#ifdef __GLIBC__
#include <malloc.h>
#endif

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

namespace {

// A field of /proc/self/status in MiB, such as VmHWM (the peak) or VmRSS; -1
// when there is none.
double ProcStatusMiB(const std::string &field)
{
   std::ifstream status("/proc/self/status");
   std::string line;
   while (std::getline(status, line)) {
      if (line.rfind(field + ":", 0) == 0)
         return std::strtod(line.substr(field.size() + 1).c_str(), nullptr) / 1024.0; // kB
   }
   return -1;
}

} // namespace

PeakMemory::PeakMemory()
{
#ifdef __GLIBC__
   // Hand memory that earlier work freed back to the kernel first, so that this
   // stretch starts from what is in use, not from what the allocator kept.
   malloc_trim(0);
#endif
   std::ofstream clear_refs("/proc/self/clear_refs");
   if (!(clear_refs << "5" << std::flush))
      return;
   m_startMiB = ProcStatusMiB("VmRSS");
}

void PeakMemory::Report(benchmark::State &state) const
{
   const double peak = ProcStatusMiB("VmHWM");
   if (m_startMiB < 0 || peak < 0)
      return;
   state.counters["peak_mib"] = peak;
   state.counters["start_mib"] = m_startMiB;
}

void StartROOT()
{
   const double before = ProcStatusMiB("VmRSS");
   gROOT->GetVersion();
   const double after = ProcStatusMiB("VmRSS");
   if (before >= 0 && after >= 0)
      benchmark::AddCustomContext("root_startup_mib", std::to_string(after - before));
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

bool ConvertToRam(const std::string &sam, const std::string &ram, int threads)
{
   bool converted = false;
   {
      const StdoutSilencer silence;
      converted = samtoramntuple(sam.c_str(), ram.c_str(), /*compression_algorithm=*/505,
                                 /*quality_policy=*/RAMNTupleRecord::kPhred33, threads);
   }
   if (!converted)
      std::cerr << "cannot convert " << sam << " to " << ram << "\n";
   return converted;
}

int64_t RamRecords(const std::string &path)
{
   const auto file = RAMFile::Open(path);
   return file ? static_cast<int64_t>(file->GetReader().GetNEntries()) : -1;
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

// Collected before the caller acts on them, so that removing files does not
// change the directory while it is being listed.
std::vector<std::string> FilesWithPrefix(const std::string &dir, const std::string &prefix)
{
   std::vector<std::string> files;
   std::error_code ec;
   for (const auto &entry : std::filesystem::directory_iterator(dir, ec)) {
      if (entry.is_regular_file(ec) && entry.path().filename().string().rfind(prefix, 0) == 0)
         files.push_back(entry.path().string());
   }
   return files;
}

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
