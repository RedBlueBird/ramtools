#ifndef BENCHMARK_UTILS_H
#define BENCHMARK_UTILS_H

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

struct BenchmarkOptions;

namespace benchmark {
class State;
} // namespace benchmark

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

/// The peak resident memory of this process over a stretch of a benchmark. On
/// construction it resets the kernel's record of the peak (Linux's
/// /proc/self/clear_refs), so the peak it reports is only what came after,
/// every thread of the process included. Memory left from earlier benchmarks
/// stays resident, so Report() gives the resident memory at the start as well;
/// the peak above it is what the stretch cost, a little less when it reuses
/// memory an earlier benchmark freed: run one benchmark per process
/// (--benchmark_filter) to compare memory exactly. Where the kernel keeps no
/// such record, nothing is reported.
class PeakMemory {
public:
   PeakMemory();

   /// Sets the counters peak_mib, the peak since construction, and start_mib,
   /// the resident memory at construction, both in MiB. Call it right after the
   /// timed loop, before any check of the results allocates.
   void Report(benchmark::State &state) const;

private:
   double m_startMiB = -1; ///< -1 when the peak cannot be reset or read
};

/// Starts ROOT before any benchmark runs, so that no benchmark's time or memory
/// takes in its one-time start-up, about 230 MiB, most of it the interpreter.
/// Records the memory the start-up took in the run's context, as
/// root_startup_mib, where the kernel reports it. Call it after
/// benchmark::Initialize().
void StartROOT();

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
   std::string stem; ///< File name, without extension, for what is made from it in the work directory
   std::string path;
   int reads = 0; ///< Synthetic records to generate; 0 for --sam
   bool sorted = false;
};

/// --sam when it is given, else one synthetic input per --reads count, written
/// in coordinate order when \p sorted. Nothing is generated yet.
std::vector<SamInput> SamInputs(const BenchmarkOptions &opts, const WorkDir &dir, bool sorted);

/// Generates \p input if it is synthetic and not there yet, outside any timing.
void PrepareSamInput(const SamInput &input);

/// Runs \p make the first time \p path is asked for in this process, and
/// remembers whether it succeeded. A file made from an input is so never taken
/// from an earlier run in a kept --work-dir, which may have had another input.
bool MakeOnce(const std::string &path, const std::function<bool()> &make);

/// Converts \p sam to the RAM file \p ram as samtoramntuple does by default
/// (ZSTD level 5, lossless qualities), on \p threads threads, with its output
/// silenced. False, with the reason printed, on failure.
bool ConvertToRam(const std::string &sam, const std::string &ram, int threads);

/// The records in a RAM file, from its METADATA; -1, with the reason printed,
/// when it cannot be opened.
int64_t RamRecords(const std::string &path);

/// The records in a SAM file, header lines not counted. Counted once per file.
int64_t CountSamRecords(const std::string &path);

/// The size of a file in MiB; 0 when it does not exist.
double FileMiB(const std::string &path);

/// Turns ROOT's implicit multithreading on with \p threads threads, as the
/// converters' -threads does for the steps that compress pages on it, or off
/// for 1.
void SetImplicitMT(int threads);

/// The files in \p dir whose names start with \p prefix.
std::vector<std::string> FilesWithPrefix(const std::string &dir, const std::string &prefix);

/// The sum of the sizes of the files in \p dir whose names start with \p prefix.
uint64_t TotalSizeWithPrefix(const std::string &dir, const std::string &prefix);

/// Removes the files in \p dir whose names start with \p prefix.
void RemoveWithPrefix(const std::string &dir, const std::string &prefix);

#endif // BENCHMARK_UTILS_H
