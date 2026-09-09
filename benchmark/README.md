# RAMtools Benchmark Suite

Performance benchmarks for the RAMtools SAM → ROOT (TTree / RNTuple) pipeline, built on
[Google Benchmark](https://github.com/google/benchmark). Every benchmark runs on
synthetic data out of the box and can be pointed at a real dataset from the command line.

Three of them — conversion time, region query, chromosome split — also measure **BAM**
alongside the two ROOT backends, via in-process htslib (see
[The BAM backend](#the-bam-backend)).

> These benchmarks measure **performance** (time, file size, throughput). Correctness is
> covered separately by the gtest suite in [`../test/`](../test). Memory profiling is not
> yet wired up (planned follow-up).

## Benchmarks

| Binary | What it measures |
|--------|------------------|
| `sam_to_ram_benchmark` | SAM→TTree vs SAM→RNTuple conversion, output sizes, compression ratio |
| `conversion_time_benchmark` | SAM→TTree, SAM→RNTuple and SAM→BAM conversion time separately, output size |
| `region_query_benchmark` | Region-query latency, TTree vs RNTuple vs BAM, across a scale ladder of regions |
| `chromosome_split_benchmark` | Per-chromosome split: RNTuple parallel writer vs BAM indexed extraction |
| `columnar_read_benchmark` | Single-column (`flag`/`mapq`) scan vs full-record read — RNTuple's columnar advantage |

## Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release   # RAMTOOLS_BUILD_BENCHMARKS=ON by default
cmake --build build
```

Binaries land in `build/benchmark/`.

## Command-line flags

Parsed by `benchmark_config.h`; stripped before Google Benchmark sees argv.

| Flag | Default | Meaning |
|------|---------|---------|
| `--sam=PATH` | *(generate)* | SAM file to benchmark |
| `--bam=PATH` | *(derive from `--sam`)* | Pre-built BAM. Must be **coordinate sorted** with a `.bai` beside it |
| `--ttree-root=PATH` | *(derive from `--sam`)* | Pre-converted TTree `.root` |
| `--rntuple-root=PATH` | *(derive from `--sam`)* | Pre-converted RNTuple `.root` |
| `--compression=N` | `505` (ZSTD-5) | ROOT compression setting for RNTuple |
| `--quality=N` | `0` | Quality policy bitmask |
| `--threads=N` | `4` | Worker threads (split benchmarks) |
| `--reads=N` | `100000` | Synthetic read count when no `--sam` given |
| `--regions=all\|i,j,k` | `0,3,6,9` | Region indices for `region_query` |
| `--regions-file=PATH` | *(built-in list)* | `region_query`: real `name:start-end` regions, one per line (`#` comments allowed), instead of the built-in synthetic-data list |

When no `--sam`/`.root`/`--bam` is supplied, temp synthetic files are generated and
cleaned up automatically.

> `--quality` selects an **RNTuple** quality policy (`0` = lossless Phred33). BAM has no
> lossy mode, so only `--quality=0` produces an RNTuple-vs-BAM number that compares like
> with like. The binaries print a warning at any other setting.

## Run

### Everything, with JSON output (one command)

```bash
cmake --build build --target benchmark
```

Runs every benchmark and writes one `*.json` result file per binary into
`build/benchmark/`.

### A single benchmark

```bash
./build/benchmark/region_query_benchmark
./build/benchmark/columnar_read_benchmark --benchmark_out=cols.json --benchmark_out_format=json
```

All standard Google Benchmark flags work (`--benchmark_filter=`, `--benchmark_repetitions=`,
`--benchmark_out=`, `--benchmark_out_format=json|csv`, …).

### On a real dataset

By default the benchmarks generate a synthetic SAM (`GenerateSAMFile`). Pass `--sam=` (and
optionally pre-converted `.root` files) to use your own data:

```bash
# Convert + query a real SAM (region_query / columnar convert once, then read):
./build/benchmark/region_query_benchmark   --sam=/data/HG002.sam --regions-file=/data/HG002_regions.txt
./build/benchmark/conversion_time_benchmark --sam=/data/HG002.sam
./build/benchmark/columnar_read_benchmark   --sam=/data/HG002.sam

# Reuse already-converted files (skip conversion):
./build/benchmark/region_query_benchmark \
    --ttree-root=/data/HG002_ttree.root --rntuple-root=/data/HG002_rntuple.root \
    --bam=/data/HG002.bam \
    --regions-file=/data/HG002_regions.txt
```

`--bam` must point at a coordinate-sorted BAM with its `.bai` alongside; without it the
BAM is derived from `--sam` on startup (sorting first if needed), which is correct but
slower. Build one with `samtools sort` + `samtools index`, or reuse an existing sorted
BAM directly.

> **`region_query_benchmark` and real data:** the built-in region list (selected by
> `--regions=`) is tuned to `GenerateSAMFile()`'s synthetic output — `"chr1"`-style
> names at hg19 lengths. A real SAM/BAM (e.g. GIAB, which uses bare `1`/`21`/`X`
> GRCh37 naming) won't match those names and will silently return 0 reads for every
> region. Pass `--regions-file=PATH` with your own `name:start-end` lines matching
> the dataset's actual `@SQ` names instead — e.g.:
> ```
> # HG002_regions.txt
> 21:20000000-20001000
> 21:20000000-21000000
> 21:1-48129895
> ```
> With `--regions-file` set (and no explicit `--regions=`), every line in the file
> is run — the `0,3,6,9` default only makes sense against the built-in 18-entry list.

To forward dataset flags to *all* binaries at once and collect JSON:

```bash
scripts/run_benchmarks.sh build -- --sam=data/HG002.sam --regions-file=data/regions.txt
```

> **Paths and the working directory.** The binaries are run from `build/benchmark/` so
> that JSON and temp files land in the build tree. `run_benchmarks.sh` rewrites the
> path-valued flags (`--sam`, `--bam`, `--ttree-root`, `--rntuple-root`,
> `--regions-file`) to absolute paths for you, so relative paths on *its* command line
> work as typed. Invoking a binary directly from the repo root does **not** get that
> rewrite — pass an absolute path, or `cd` to the build directory first.
>
> A path that does not exist is now a hard error: every binary checks its dataset flags
> up front and exits rather than running. This used to fail silently — `samtoram()` and
> `samtoramntuple()` accept a missing input, write an essentially empty file, and the
> suite would report convincing sub-millisecond timings over a couple of KB of nothing.

## The BAM backend

`conversion_time_benchmark`, `region_query_benchmark` and `chromosome_split_benchmark`
each measure BAM next to the ROOT backends. BAM is driven through **in-process htslib**
(`bam_ops.h`), not by shelling out to `samtools`: that keeps `fork`/`exec` out of every
timed iteration, and lets file opening and index loading be placed on exactly the same
side of the timer as the ROOT backends' equivalents. htslib is already a `PUBLIC`
dependency of `ramcore`, so no extra build wiring is involved and `samtools` no longer
needs to be installed to run the suite.

What each benchmark actually times:

| Benchmark | RNTuple / TTree | BAM |
|---|---|---|
| `conversion_time` | one conversion pass, `index=true` | write BAM + build `.bai` (+ sort, if the input is not already coordinate sorted) |
| `region_query` | `ramview()` / `ramntupleview()`, which reopen the file per call | `sam_open` + `sam_index_load` + iterate, **also per call** |
| `chromosome_split` | `samtoramntuple_split_by_chromosome()`: parse SAM → group → sort → parallel write | `BamConvert()` then `BamSplit()`: parse SAM → sort → BAM + `.bai` → indexed extraction |

Three choices worth knowing about when reading the numbers:

- **Record counts are aligned.** The htslib iterator is filtered with `0x904` (unmapped,
  secondary, supplementary) so it returns the same set `RAMNTupleView` does. On real
  sorted data RNTuple and BAM agree exactly; TTree returns slightly *more*, because it
  does no flag filtering at all.
- **Sizes include the index.** A BAM's `.bai` is counted alongside the `.bam`, since an
  RNTuple carries its index inside the single `.root`. On small synthetic runs the `.bai`
  dominates — 1 000 reads scattered over 26 hg19-length chromosomes produce a sparse,
  disproportionately large index — so only read the size columns at realistic scale.
- **The BAM split writes an intermediate.** Producing per-chromosome BAMs requires a
  sorted, indexed BAM first, and that cost is inside the timed region because it is what
  the workflow genuinely costs. Only the per-chromosome outputs are counted in `size_MB`.

### Coordinate sorting

`ramview()` and `ramntupleview()` seek to a row derived from the sparse index and then
**stop scanning at the first record past the region's end**. That reaches every
overlapping read only if the file is coordinate sorted. `GenerateSAMFile()` writes
`SO:unsorted`, so on synthetic data both ROOT backends under-report matches, while the
BAM side — which must be sorted before it can be indexed at all — reports the true count.
`region_query_benchmark` prints a warning when its input is unsorted; treat the backend
comparison as meaningful only on a coordinate-sorted `--sam`/`--bam`.

## Reading the results

Each row reports wall time, CPU time, iterations, and benchmark-specific counters:
`file_size_mb`, `ttree_size_mb` / `rntuple_size_mb` / `compression_ratio`,
`reads_per_second`, `size_MB`, `threads`, `region_idx`. The label column shows the read
count for the run (e.g. region matches, dataset rows).

BAM conversion rows add two more: `indexed` (1 when a `.bai` was produced, so the output
is actually queryable) and `sorted_in_loop` (1 when the input was not coordinate sorted
and the sort was therefore paid inside the timed region — true for every synthetic run,
false for an already-sorted real dataset).

### Notebook

The JSON files are visualized by [`io_performance/benchmark_report.ipynb`](io_performance/benchmark_report.ipynb),
a self-contained notebook that loads every `*.json` in `build/benchmark/` and plots each
benchmark (sizes/compression, conversion time, region-query latency, chromosome split,
columnar read), with BAM shown alongside TTree and RNTuple wherever it is measured.
Missing files are skipped rather than erroring.

```bash
python3 -m venv .venv && source .venv/bin/activate      # system python is often PEP-668 managed
pip install -r benchmark/io_performance/requirements.txt
jupyter lab benchmark/io_performance/benchmark_report.ipynb
```

By default the notebook reads `../../build/benchmark`; set `RESULTS_DIR` in the second cell to
analyze a different run.

> **Sourcing ROOT?** `thisroot.sh` sets `JUPYTER_CONFIG_DIR` to ROOT's own
> `$ROOTSYS/etc/notebook`, which is read-only for a system-wide install and makes
> `jupyter lab` abort on first-run config migration (`PermissionError: .../migrated`).
> Launch with a writable config dir instead:
> ```bash
> JUPYTER_CONFIG_DIR="$HOME/.jupyter" jupyter lab benchmark/io_performance/benchmark_report.ipynb
> ```

> The older `ramview_perf.ipynb` / `samtoram_perf.ipynb` notebooks are a superseded workflow
> that parsed a `.perf`/`.log` text format and do not read this JSON output.

## Reproducibility

When reporting numbers, record the environment:

```
Date:
CPU model / cores:
RAM:
Disk (SSD/NVMe/HDD):
OS / kernel:
ROOT version:
htslib version (`pkg-config --modversion htslib`):
Build type: Release
Dataset (name, #reads, source):
Command line:
```

## Notes

- The suite no longer requires a `samtools` binary: `chromosome_split_benchmark` used to
  shell out to `samtools view/sort/index`, and now uses in-process htslib instead.
