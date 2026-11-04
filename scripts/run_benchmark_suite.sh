#!/usr/bin/env bash
#
# Runs every benchmark binary of a build and writes one JSON result per binary.
#
#   scripts/run_benchmark_suite.sh [BUILD_DIR] [-- FLAG...]
#
#   scripts/run_benchmark_suite.sh build
#   scripts/run_benchmark_suite.sh build -- --sam=reads.sam --ref=ref.fa --threads=1,4
#
# Each FLAG goes to the binaries that take it: --benchmark_* flags to all of
# them, the suite's own flags (--sam, --threads, ...) to those whose --help
# lists them, so one command line serves the whole suite. A flag no binary
# takes is an error, as it is for a single binary.
#
# The results go to BUILD_DIR/benchmark/results, or to $RAMTOOLS_BENCH_RESULTS.
# Exits with 1 when a binary fails or a benchmark reports an error, such as a
# count that disagrees with htslib's.

set -euo pipefail

build_dir=build
if [[ $# -gt 0 && $1 != -- ]]; then
   build_dir=$1
   shift
fi
if [[ ${1:-} == -- ]]; then
   shift
fi
flags=("$@")

binaries=(conversion_time_benchmark region_query_benchmark chromosome_split_benchmark columnar_read_benchmark)
results=${RAMTOOLS_BENCH_RESULTS:-$build_dir/benchmark/results}

# The suite's flags each binary takes, as its --help lists them.
declare -A takes
for binary in "${binaries[@]}"; do
   path=$build_dir/benchmark/$binary
   if [[ ! -x $path ]]; then
      echo "run_benchmark_suite: $path is not built; build with -DRAMTOOLS_BUILD_BENCHMARKS=ON" >&2
      exit 1
   fi
   takes[$binary]=" $("$path" --help 2>/dev/null | sed -n 's/^  \(--[a-z-]*\)=.*/\1/p' | tr '\n' ' ')"
done

for flag in ${flags[@]+"${flags[@]}"}; do
   name=${flag%%=*}
   [[ $name == --benchmark_* ]] && continue
   taken=false
   for binary in "${binaries[@]}"; do
      [[ ${takes[$binary]} == *" $name "* ]] && taken=true
   done
   if [[ $taken == false ]]; then
      echo "run_benchmark_suite: no benchmark takes $name" >&2
      exit 1
   fi
done

mkdir -p "$results"
failed=()
for binary in "${binaries[@]}"; do
   args=()
   for flag in ${flags[@]+"${flags[@]}"}; do
      name=${flag%%=*}
      if [[ $name == --benchmark_* || ${takes[$binary]} == *" $name "* ]]; then
         args+=("$flag")
      fi
   done
   out=$results/$binary.json
   echo "=== $binary ==="
   if ! "$build_dir/benchmark/$binary" ${args[@]+"${args[@]}"} --benchmark_out="$out" --benchmark_out_format=json; then
      failed+=("$binary (exit status)")
   elif grep -q '"error_occurred": true' "$out"; then
      failed+=("$binary (a benchmark reported an error)")
   fi
   echo
done

echo "Results in $results"
if [[ ${#failed[@]} -gt 0 ]]; then
   printf 'run_benchmark_suite: failed: %s\n' "${failed[@]}" >&2
   exit 1
fi
