#include "benchmark_config.h"

#include "rntuple/RAMNTupleRecord.h"

#include <array>
#include <cctype>
#include <climits>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace {

struct FlagSpec {
   BenchmarkFlag flag;
   const char *name;
   const char *value;
   const char *help;
};

constexpr std::array<FlagSpec, 8> kFlagSpecs{{
   {kFlagCompression, "--compression", "N", "ROOT compression code, algorithm*100+level"},
   {kFlagQuality, "--quality", "MODE", "lossless (rANS), fqzcomp, illumina or drop"},
   {kFlagReads, "--reads", "N[,N...]", "synthetic records, one benchmark per count"},
   {kFlagThreads, "--threads", "N[,N...]", "threads, one benchmark per count"},
   {kFlagRam, "--ram", "FILE", "RAM file written by samtoramntuple, or $RAMTOOLS_BENCH_RNTUPLE"},
   {kFlagRegions, "--regions", "R[,R...]",
    "regions to query, or $RAMTOOLS_BENCH_REGIONS (default: 1 kb to a whole reference)"},
   {kFlagSam, "--sam", "FILE", "SAM file to run on, instead of synthetic records (not with --reads)"},
   {kFlagWorkDir, "--work-dir", "DIR", "where inputs and outputs go, kept afterwards (default: a temporary directory)"},
}};

// ROOT compression code: algorithm*100+level, with algorithm 1 ZLIB, 2 LZMA,
// 4 LZ4 or 5 ZSTD and level 1 to 9; 0 means no compression. The same rule as
// the converters' -compression.
bool ParseCompression(const std::string &text, int &code)
{
   char *end = nullptr;
   const long value = std::strtol(text.c_str(), &end, 10);
   if (text.empty() || *end != '\0' || value < 0)
      return false;
   const long algorithm = value / 100;
   const long level = value % 100;
   const bool known = algorithm == 1 || algorithm == 2 || algorithm == 4 || algorithm == 5;
   if (value != 0 && (!known || level < 1 || level > 9))
      return false;
   code = static_cast<int>(value);
   return true;
}

// The quality policies samtoramntuple's flags give: none, -fqzcomp, -illumina
// and -dropqual.
bool ParseQuality(const std::string &text, uint32_t &policy)
{
   if (text == "lossless")
      policy = RAMNTupleRecord::kPhred33;
   else if (text == "fqzcomp")
      policy = RAMNTupleRecord::kPhred33 | RAMNTupleRecord::kFqzcomp;
   else if (text == "illumina")
      policy = RAMNTupleRecord::kIlluminaBinning;
   else if (text == "drop")
      policy = RAMNTupleRecord::kDrop;
   else
      return false;
   return true;
}

// The comma-separated items of text. A trailing comma gives an empty last
// item, so that "1000," is an error rather than 1000.
std::vector<std::string> SplitList(const std::string &text)
{
   std::vector<std::string> items;
   std::stringstream in(text);
   std::string item;
   while (std::getline(in, item, ','))
      items.push_back(item);
   if (!text.empty() && text.back() == ',')
      items.emplace_back();
   return items;
}

// A positive decimal integer and nothing else, so "3x" is an error rather
// than 3.
bool ParsePositive(const std::string &text, int &number)
{
   char *end = nullptr;
   const long value = std::strtol(text.c_str(), &end, 10);
   if (text.empty() || !std::isdigit(static_cast<unsigned char>(text[0])) || *end != '\0' || value < 1 ||
       value > INT_MAX)
      return false;
   number = static_cast<int>(value);
   return true;
}

bool ParsePositiveList(const std::string &text, std::vector<int> &numbers)
{
   std::vector<int> parsed;
   for (const std::string &item : SplitList(text)) {
      int number = 0;
      if (!ParsePositive(item, number))
         return false;
      parsed.push_back(number);
   }
   if (parsed.empty())
      return false;
   numbers = std::move(parsed);
   return true;
}

bool ParseRegions(const std::string &text, std::vector<std::string> &regions)
{
   std::vector<std::string> parsed = SplitList(text);
   if (parsed.empty())
      return false;
   for (const std::string &region : parsed) {
      if (region.empty())
         return false;
   }
   regions = std::move(parsed);
   return true;
}

// A file that exists, so that a mistyped path fails here instead of timing a
// converter or a query that has nothing to read.
bool ParseFile(const std::string &text, std::string &path)
{
   std::error_code ec;
   if (text.empty() || !std::filesystem::is_regular_file(text, ec))
      return false;
   path = text;
   return true;
}

bool ParseValue(BenchmarkFlag flag, const std::string &value, BenchmarkOptions &opts)
{
   switch (flag) {
   case kFlagCompression: return ParseCompression(value, opts.compression);
   case kFlagQuality: return ParseQuality(value, opts.quality_policy);
   case kFlagReads: return ParsePositiveList(value, opts.reads);
   case kFlagThreads: return ParsePositiveList(value, opts.threads);
   case kFlagRam: return ParseFile(value, opts.ram);
   case kFlagRegions: return ParseRegions(value, opts.regions);
   case kFlagSam: return ParseFile(value, opts.sam);
   case kFlagWorkDir:
      if (value.empty())
         return false;
      opts.work_dir = value;
      return true;
   }
   return false;
}

// What else than its form makes a value of this flag invalid.
const char *InvalidBecause(BenchmarkFlag flag)
{
   return (flag == kFlagRam || flag == kFlagSam) ? ": no such file" : "";
}

std::string Join(const std::vector<int> &numbers)
{
   std::string text;
   for (const int number : numbers) {
      if (!text.empty())
         text += ',';
      text += std::to_string(number);
   }
   return text;
}

std::string Join(const std::vector<std::string> &items)
{
   std::string text;
   for (const std::string &item : items) {
      if (!text.empty())
         text += ',';
      text += item;
   }
   return text;
}

std::string QualityName(uint32_t policy)
{
   if (policy & RAMNTupleRecord::kDrop)
      return "drop";
   if (policy & RAMNTupleRecord::kIlluminaBinning)
      return "illumina";
   if (policy & RAMNTupleRecord::kFqzcomp)
      return "fqzcomp";
   return "lossless";
}

// The value a flag takes when it is not given, as this binary set it.
std::string DefaultOf(BenchmarkFlag flag, const BenchmarkOptions &opts)
{
   switch (flag) {
   case kFlagCompression: return std::to_string(opts.compression);
   case kFlagQuality: return QualityName(opts.quality_policy);
   case kFlagReads: return Join(opts.reads);
   case kFlagThreads: return Join(opts.threads);
   case kFlagRam: return opts.ram;
   case kFlagRegions: return Join(opts.regions);
   case kFlagSam: return opts.sam;
   case kFlagWorkDir: return opts.work_dir;
   }
   return "";
}

void PrintUsage(unsigned accepted, const BenchmarkOptions &opts)
{
   std::cout << "Benchmark options:\n";
   for (const FlagSpec &spec : kFlagSpecs) {
      if (!(accepted & spec.flag))
         continue;
      std::cout << "  " << spec.name << "=" << spec.value << "\n      " << spec.help;
      const std::string fallback = DefaultOf(spec.flag, opts);
      if (!fallback.empty())
         std::cout << " (default " << fallback << ")";
      std::cout << "\n";
   }
   std::cout << "\n";
}

// The value of an environment variable, or "" when it is not set.
std::string Env(const char *name)
{
   const char *value = std::getenv(name); // NOLINT(concurrency-mt-unsafe): read before any thread starts
   return value ? value : "";
}

} // namespace

bool ParseBenchmarkFlags(int *argc, char **argv, unsigned accepted, BenchmarkOptions &opts)
{
   const BenchmarkOptions defaults = opts;
   unsigned given = 0;
   int kept = 1;
   for (int i = 1; i < *argc; ++i) {
      const std::string arg = argv[i]; // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): main()'s argv
      if (arg == "--help")
         PrintUsage(accepted, defaults);

      const FlagSpec *match = nullptr;
      for (const FlagSpec &spec : kFlagSpecs) {
         const std::string name = spec.name;
         if ((accepted & spec.flag) && arg.compare(0, name.size(), name) == 0 &&
             (arg.size() == name.size() || arg[name.size()] == '=')) {
            match = &spec;
            break;
         }
      }
      if (!match) {
         argv[kept++] = argv[i]; // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): main()'s argv
         continue;
      }

      const std::string name = match->name;
      if (arg.size() == name.size()) {
         std::cerr << name << " needs a value: " << name << "=" << match->value << "\n";
         return false;
      }
      const std::string value = arg.substr(name.size() + 1);
      if (!ParseValue(match->flag, value, opts)) {
         std::cerr << "invalid " << name << " value '" << value << "'" << InvalidBecause(match->flag) << "\n";
         return false;
      }
      given |= match->flag;
   }
   if ((given & kFlagSam) && (given & kFlagReads)) {
      std::cerr << "--sam and --reads exclude each other: --reads sizes the synthetic input\n";
      return false;
   }
   if ((given & kFlagRam) && (given & (kFlagSam | kFlagReads))) {
      std::cerr << "--ram excludes --sam and --reads: the benchmark runs on the RAM file as it is\n";
      return false;
   }
   argv[kept] = nullptr; // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): main()'s argv
   *argc = kept;

   // The variables the region benchmark read before it took flags; a flag wins,
   // and so does an input given another way.
   struct EnvFallback {
      BenchmarkFlag flag;
      const char *variable;
   };
   for (const EnvFallback fallback :
        {EnvFallback{kFlagRam, "RAMTOOLS_BENCH_RNTUPLE"}, EnvFallback{kFlagRegions, "RAMTOOLS_BENCH_REGIONS"}}) {
      if (!(accepted & fallback.flag) || (given & fallback.flag))
         continue;
      if (fallback.flag == kFlagRam && (given & (kFlagSam | kFlagReads)))
         continue;
      const std::string value = Env(fallback.variable);
      if (!value.empty() && !ParseValue(fallback.flag, value, opts)) {
         std::cerr << "invalid " << fallback.variable << " value '" << value << "'" << InvalidBecause(fallback.flag)
                   << "\n";
         return false;
      }
   }
   return true;
}
