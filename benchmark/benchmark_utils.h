#ifndef BENCHMARK_UTILS_H
#define BENCHMARK_UTILS_H

#include <cstdint>
#include <string>
#include <vector>

struct BenchmarkOptions;

/// Sends standard output to /dev/null while it lives, so that what a converter
/// or a query prints per call is not timed as terminal output. It redirects the
/// file descriptor, which catches printf, std::cout and ROOT's own printing
/// alike, and puts the original back on destruction, also when output is piped.
/// Standard error is left alone, so failures still show.
class StdoutSilencer {
public:
   StdoutSilencer();
   ~StdoutSilencer();
   StdoutSilencer(const StdoutSilencer &) = delete;
   StdoutSilencer &operator=(const StdoutSilencer &) = delete;
   StdoutSilencer(StdoutSilencer &&) = delete;
   StdoutSilencer &operator=(StdoutSilencer &&) = delete;

private:
   int m_saved = -1; ///< A duplicate of the original stdout; -1 when nothing was redirected
};

/// The directory a benchmark writes its inputs and outputs in. Given a path, it
/// uses that directory, creating it if needed, and leaves it in place. Given
/// none, it makes a new directory under the system's temporary directory and
/// removes it with everything in it on destruction.
class WorkDir {
public:
   /// Exits with a message when the directory cannot be created.
   explicit WorkDir(const std::string &path);
   ~WorkDir();
   WorkDir(const WorkDir &) = delete;
   WorkDir &operator=(const WorkDir &) = delete;
   WorkDir(WorkDir &&) = delete;
   WorkDir &operator=(WorkDir &&) = delete;

   [[nodiscard]] const std::string &Path() const { return m_path; }
   /// The path of a file called \p name in the directory.
   [[nodiscard]] std::string File(const std::string &name) const;

private:
   std::string m_path;
   bool m_owned = false; ///< Made here, so removed here
};

/// A SAM file to run a benchmark on: the one given with --sam, or synthetic
/// records generated into the work directory on first use.
struct SamInput {
   std::string name; ///< "reads:N" for synthetic records, "sam" for --sam
   std::string path;
   int reads = 0; ///< Synthetic records to generate; 0 for --sam
   bool sorted = false;
};

/// --sam when it is given, else one synthetic input per --reads count, written
/// in coordinate order when \p sorted. Nothing is generated yet.
std::vector<SamInput> SamInputs(const BenchmarkOptions &opts, const WorkDir &dir, bool sorted);

/// Generates \p input if it is synthetic and not there yet, outside any timing.
void PrepareSamInput(const SamInput &input);

/// The records in a SAM file, header lines not counted. Counted once per file.
int64_t CountSamRecords(const std::string &path);

/// The sum of the sizes of the files in \p dir whose names start with \p prefix.
uint64_t TotalSizeWithPrefix(const std::string &dir, const std::string &prefix);

/// Removes the files in \p dir whose names start with \p prefix.
void RemoveWithPrefix(const std::string &dir, const std::string &prefix);

#endif // BENCHMARK_UTILS_H
