#ifndef BENCHMARK_UTILS_H
#define BENCHMARK_UTILS_H

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

#endif // BENCHMARK_UTILS_H
