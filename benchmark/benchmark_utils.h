#ifndef BENCHMARK_UTILS_H
#define BENCHMARK_UTILS_H

#include <cstdio>
#include <filesystem>
#include <stdio.h> // NOLINT(modernize-deprecated-headers) - fileno is a POSIX extension, not in <cstdio>
#include <string>
#include <system_error>
#include <utility>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
constexpr const char *NULL_DEVICE = "NUL";
#define BENCH_DUP _dup
#define BENCH_DUP2 _dup2
#define BENCH_FILENO _fileno
#define BENCH_OPEN _open
#define BENCH_CLOSE _close
#define BENCH_WRONLY _O_WRONLY
#else
#include <fcntl.h>
#include <unistd.h>
constexpr const char *NULL_DEVICE = "/dev/null";
#define BENCH_DUP dup
#define BENCH_DUP2 dup2
#define BENCH_FILENO fileno
#define BENCH_OPEN open
#define BENCH_CLOSE close
#define BENCH_WRONLY O_WRONLY
#endif

namespace benchutil {

// Redirects stdout (and optionally stderr) to NULL_DEVICE for the lifetime of the
// object, restoring the original streams on destruction. Unlike the previous
// freopen("/dev/tty") approach, this saves/restores the real file descriptors with
// dup/dup2, so it works when there is no controlling terminal (CI, redirected runs).
class ScopedStdoutSuppressor {
public:
   explicit ScopedStdoutSuppressor(bool suppress_stderr = false)
      : m_savedStdout(BENCH_DUP(BENCH_FILENO(stdout))),
        m_devnull(BENCH_OPEN(NULL_DEVICE, BENCH_WRONLY)) // NOLINT(cppcoreguidelines-pro-type-vararg)
   {
      std::fflush(stdout);
      if (m_devnull >= 0)
         BENCH_DUP2(m_devnull, BENCH_FILENO(stdout));
      if (suppress_stderr && m_devnull >= 0) {
         std::fflush(stderr);
         m_savedStderr = BENCH_DUP(BENCH_FILENO(stderr));
         BENCH_DUP2(m_devnull, BENCH_FILENO(stderr));
      }
   }

   ~ScopedStdoutSuppressor()
   {
      std::fflush(stdout);
      if (m_savedStdout >= 0) {
         BENCH_DUP2(m_savedStdout, BENCH_FILENO(stdout));
         BENCH_CLOSE(m_savedStdout);
      }
      if (m_savedStderr >= 0) {
         std::fflush(stderr);
         BENCH_DUP2(m_savedStderr, BENCH_FILENO(stderr));
         BENCH_CLOSE(m_savedStderr);
      }
      if (m_devnull >= 0)
         BENCH_CLOSE(m_devnull);
   }

   ScopedStdoutSuppressor(const ScopedStdoutSuppressor &) = delete;
   ScopedStdoutSuppressor &operator=(const ScopedStdoutSuppressor &) = delete;
   ScopedStdoutSuppressor(ScopedStdoutSuppressor &&) = delete;
   ScopedStdoutSuppressor &operator=(ScopedStdoutSuppressor &&) = delete;

private:
   int m_savedStdout = -1;
   int m_savedStderr = -1;
   int m_devnull = -1;
};

// Size of one file in bytes, or 0 if it does not exist.
inline std::size_t FileSize(const std::string &path)
{
   std::error_code ec;
   const auto size = std::filesystem::file_size(path, ec);
   return ec ? 0 : static_cast<std::size_t>(size);
}

namespace detail {

// Split "dir/stem" into the directory to scan and the filename prefix to match. A bare
// prefix scans the current directory, which is where the benchmarks write their output.
inline std::pair<std::filesystem::path, std::string> SplitPrefix(const std::string &prefix)
{
   const std::filesystem::path p(prefix);
   return {p.has_parent_path() ? p.parent_path() : std::filesystem::path("."), p.filename().string()};
}

} // namespace detail

// Sum the sizes of every file whose name *starts with* `prefix`.
//
// Prefer this to GetTotalFileSize: anchoring at the start of the filename stops one
// backend's output from being counted under another's pattern, and it does not depend on
// chromosome names happening to begin with "chr" -- real datasets often use bare `1`,
// `21`, `X` (GRCh37), which a "..._chr" pattern misses entirely, silently reporting 0 MB.
// Sidecar index files (.bai, .csi) share the prefix and are therefore included, so a BAM
// split's total is comparable with an RNTuple's self-contained .root.
inline std::size_t TotalSizeByPrefix(const std::string &prefix)
{
   const auto [dir, stem] = detail::SplitPrefix(prefix);
   std::error_code ec;
   std::size_t total = 0;
   for (const auto &entry : std::filesystem::directory_iterator(dir, ec)) {
      if (entry.is_regular_file(ec) && entry.path().filename().string().rfind(stem, 0) == 0)
         total += FileSize(entry.path().string());
   }
   return total;
}

// Remove every file whose name starts with `prefix` (the counterpart to TotalSizeByPrefix).
inline void CleanupByPrefix(const std::string &prefix)
{
   const auto [dir, stem] = detail::SplitPrefix(prefix);
   std::error_code ec;
   for (const auto &entry : std::filesystem::directory_iterator(dir, ec)) {
      if (entry.is_regular_file(ec) && entry.path().filename().string().rfind(stem, 0) == 0)
         std::remove(entry.path().string().c_str());
   }
}

// Sum the sizes of all files in the current directory whose name contains `pattern`.
inline std::size_t GetTotalFileSize(const std::string &pattern)
{
   std::size_t total = 0;
   for (const auto &entry : std::filesystem::directory_iterator(".")) {
      if (entry.path().filename().string().find(pattern) != std::string::npos)
         total += std::filesystem::file_size(entry.path());
   }
   return total;
}

// Remove all files in the current directory whose name contains `pattern`.
inline void CleanupFiles(const std::string &pattern)
{
   for (const auto &entry : std::filesystem::directory_iterator(".")) {
      if (entry.path().filename().string().find(pattern) != std::string::npos)
         std::remove(entry.path().string().c_str());
   }
}

} // namespace benchutil

#endif
