# **PrecPack**

<p align="center">
  <img src="assets/precpack-cover.svg" alt="PrecPack: exact precedence-constrained packing" width="86%">
</p>

<p align="center">
  <b>An Efficient Open-Source Exact Solver for Bin Packing with Generalized Precedence Constraints</b>
</p>
<p align="center">
  <a href="https://github.com/Sunkanghong-Wang/PrecPack/actions/workflows/build.yml"><img alt="Build" src="https://github.com/Sunkanghong-Wang/PrecPack/actions/workflows/build.yml/badge.svg"></a>
  <a href="#requirements"><img alt="Language: C++" src="https://img.shields.io/badge/language-C%2B%2B-00599C?logo=cplusplus&amp;logoColor=white"></a>
  <a href="#build"><img alt="Build: CMake" src="https://img.shields.io/badge/build-CMake-064F8C?logo=cmake&amp;logoColor=white"></a>
  <a href="#command-line"><img alt="Problems: SALBP-I, BPP-P, and BPP-GP" src="https://img.shields.io/badge/problems-SALBP--I%20%7C%20BPP--P%20%7C%20BPP--GP-2b9348"></a>
  <a href="#overview"><img alt="Solver: exact" src="https://img.shields.io/badge/solver-exact-f59f00"></a>
  <a href="#build"><img alt="Platforms: macOS, Linux, and Windows" src="https://img.shields.io/badge/platforms-macOS%20%7C%20Linux%20%7C%20Windows-6c63ff"></a>
  <a href="LICENSE"><img alt="License: MIT" src="https://img.shields.io/badge/license-MIT-green"></a>
</p>

<p align="center">
  <a href="#cite">Cite</a> |
  <a href="#overview">Overview</a> |
  <a href="#requirements">Requirements</a> |
  <a href="#quick-start">Quick Start</a> |
  <a href="#benchmark-instances">Instances</a> |
  <a href="#contact">Contact</a> |
  <a href="#license">License</a>
</p>

## Cite

**PrecPack** is continuously maintained in this repository as a solver. A separate INFORMS Journal on Computing (IJOC) companion repository will archive the computational results and the source-code snapshot used in the paper. Once the final identifiers are assigned, please cite both the paper and that archival repository.

The following IJOC-formatted identifiers are deliberate placeholders so that the final digital object identifier (DOI) and repository address can be inserted directly when they become available.

- Paper DOI: `10.1287/ijoc.XXXX.YYYY`
- Repository DOI: `10.1287/ijoc.XXXX.YYYY.cd`

Below is the provisional BibTeX for citing the IJOC archival repository.

```bibtex
@misc{Wang2026PrecPack,
  author    = {Wang, Sunkanghong and You, Zhengzhong Ricky and Baldacci, Roberto and
               Mo, Baichuan and Qin, Hu and Wei, Lijun and Xu, Zhou},
  publisher = {INFORMS Journal on Computing},
  title     = {{PrecPack: An Efficient Open-Source Exact Solver for Bin Packing
                with Generalized Precedence Constraints}},
  year      = {2026},
  doi       = {10.1287/ijoc.XXXX.YYYY.cd},
  url       = {https://github.com/INFORMSJoC/XXXX.YYYY},
  note      = {Available for download at https://github.com/INFORMSJoC/XXXX.YYYY},
}
```

Citation metadata for the continuously maintained **PrecPack** solver is available in [`CITATION.cff`](CITATION.cff), and the complete author list is recorded in [`AUTHORS`](AUTHORS).

## Overview

**PrecPack** is a single, efficient open-source exact solver for the Bin Packing Problem with Generalized Precedence Constraints (BPP-GP), which encompasses the following two classical problems as special cases:

- Simple Assembly Line Balancing Problem type I (SALBP-I);
- Bin Packing Problem with Precedence Constraints (BPP-P).

All three problems use the same branch-bound-and-remember (BBR) search, with state representations and applicable rules determined by their precedence weights.

This repository contains the solver source code, regression tests, benchmark instances, and build/run documentation. Generated results, experiment logs, profiler output, and build artifacts are excluded from version control.

```text
PrecPack/
├── .github/
│   └── workflows/build.yml    # Gurobi-free three-platform tests
├── assets/                     # README artwork
├── code/
│   ├── CMakeLists.txt          # Build definition
│   ├── cmake/                  # Optional Gurobi discovery
│   ├── include/precpack/       # Shared declarations and data types
│   ├── scripts/                # Build/run tools and regression tests
│   └── src/
│       ├── bbr.cpp             # Exact BBR search
│       ├── initial_bounds.cpp  # Preprocessing, bounds, and primal heuristics
│       ├── root_column_generation.cpp  # Position-free root bound
│       ├── bin_packing_bound.cpp       # Ordinary bin-packing bound
│       ├── conflict_bin_packing.cpp    # Conflict-aware bin-packing bound
│       ├── dff.cpp                     # Dual-feasible functions
│       └── ...                         # Input/output, configuration, and orchestration
├── data/
│   ├── README.md               # Data formats and sources
│   ├── instances/              # Shared Otto and Scholl item data
│   └── bpp-gp-graphs/          # Generalized-precedence graph data
├── .gitattributes
├── .gitignore
├── AUTHORS
├── CITATION.cff
├── LICENSE
└── README.md
```

All implementation and developer-facing tooling lives under `code/`. The components are separated by algorithmic responsibility, while **PrecPack** remains one focused solver rather than a plug-in framework.

## Requirements

| Dependency | Requirement | Notes |
| --- | --- | --- |
| Operating system | 64-bit macOS, Linux, or Windows | The supplied build and test workflow covers these platforms. |
| C++ compiler and standard library | C++20 support | Required features include `<bit>` operations and associative-container `contains()`. |
| CMake | 3.20 or newer | Enforced by the build definition. |
| Build tool | A tool supported by the selected CMake generator | For example, Make or Ninja on macOS/Linux, or MSBuild with Visual Studio on Windows. |
| Launcher shell | Bash on macOS/Linux; Command Prompt on Windows | Required only for the supplied shell or batch wrappers; manual CMake commands are also available. |
| Gurobi Optimizer | Optional; the CMake finder requests version 9.1 or newer | A compatible C++ library and a valid license are required to use the Gurobi-dependent components. |

**PrecPack** has no mandatory mathematical-programming-solver dependency. The compiler, CMake, and a compatible build tool suffice for a commercial-solver-free build. Python is not required by the supplied build or run launchers.

### Optional Root-Bound Strengthening (Gurobi)

Gurobi supports root-bound strengthening and reference tests, not the BBR search itself. The root lower bound is computed from column-generation dual values using fixed-point arithmetic, ensuring numerical validity for pruning and optimality decisions. Root strengthening may improve the lower bound but is not guaranteed to reduce solve time. A build without Gurobi retains preprocessing, primal heuristics, complete BBR, and assignment validation.

The CMake option `PRECPACK_GUROBI` selects one of three build modes:

| Value | Behavior |
| --- | --- |
| `AUTO` | Default. Enable Gurobi support when a compatible installation is found; otherwise build without it. |
| `OFF` | Exclude the Gurobi-dependent source files and library dependency. No Gurobi installation or license is needed. |
| `ON` | Require Gurobi and stop configuration if a suitable installation cannot be found. |

The version floor used by CMake is not a claim that every later release and compiler combination has been tested. The manuscript experiments use Gurobi 13.0.2. The installed Gurobi C++ library must be binary-compatible with the compiler and platform used to build **PrecPack**. Gurobi calls use one thread. For deployment without a commercial license, select `OFF`.

To enable Gurobi support, set `GUROBI_HOME` to the platform directory containing its `include`, `lib`, and, on Windows, `bin` directories:

| Platform | Example `GUROBI_HOME` |
| --- | --- |
| macOS | `/Library/gurobi<version>/macos_universal2` |
| Linux x86-64 | `/opt/gurobi<version>/linux64` |
| Windows x64 | `C:\gurobi<version>\win64` |

## Quick Start

### Build

The recommended commercial-solver-free build requires no Gurobi installation or license.

macOS or Linux:

```bash
./code/scripts/build.sh --gurobi off
```

Windows Command Prompt:

```bat
code\scripts\build.bat --gurobi off
```

The macOS, Linux, and Windows build launchers call CMake directly. On macOS, `build.sh` also detects CMake installed with CMake.app, Homebrew, MacPorts, or CLion when it is not on `PATH`. Set `CMAKE_BIN` to an executable path to override automatic detection. None of the build or run launchers requires Python.

Omit `--gurobi off` to use the default `AUTO` detection. To require the optional Gurobi backend, set `GUROBI_HOME` and use `--gurobi on`:

```bash
export GUROBI_HOME=/path/to/gurobi/platform
./code/scripts/build.sh --gurobi on
```

```bat
set "GUROBI_HOME=C:\path\to\gurobi\win64"
code\scripts\build.bat --gurobi on
```

These launchers configure a Release build, compile **PrecPack**, and run the test suite. Every configuration runs the commercial-solver-free regressions for assignment validation, independent brute-force exactness oracles, dual-feasible-function bounds and the ordinary-bin-packing lower-bound routine (BINLB), conflict-aware BINLB optima and cutoff certificates, initialization, time-limit handling, and the public interface. When Gurobi is available, root certificates and reference integer-programming models are tested in addition. Compilation may use multiple jobs; every solver run is single-threaded.

The checked-in GitHub Actions workflow performs the Gurobi-free Release build and test suite on macOS, Linux, and Windows. Optional Gurobi builds remain local because they require a separately licensed installation.

Equivalent manual commands for a Gurobi-free build are:

```bash
cmake -S code -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DPRECPACK_GUROBI=OFF
cmake --build build --config Release --parallel
ctest --test-dir build -C Release --output-on-failure
```

For a required Gurobi build, replace `-DPRECPACK_GUROBI=OFF` with `-DPRECPACK_GUROBI=ON -DGUROBI_ROOT="$GUROBI_HOME"`. Run `precpack --help` after compilation to see whether optional Gurobi root strengthening is enabled in that executable.

With a standard Windows multi-configuration generator, executables are placed under `build\Release`; single-configuration generators place them directly under `build`.

### Command Line

After building, solve one included instance from each supported problem:

macOS or Linux:

```bash
./code/scripts/run_examples.sh
```

Windows Command Prompt:

```bat
code\scripts\run_examples.bat
```

Example outputs are written below `results/examples/`, which is ignored by Git.

The single-instance solver interface is:

```text
precpack --problem TYPE --instance FILE [options]
```

The production command line contains only the problem type, input locations, resource limits, and the output location.

Memory limits use binary units: one mebibyte (MiB) is `2^20` bytes, and one gibibyte (GiB) is `2^30` bytes. Despite the `mb` suffix, `--memory-limit-mb` and the CSV field `memory_limit_mb` both use MiB, not decimal megabytes (MB, `10^6` bytes). Their names are retained for compatibility.

| Option | Meaning | Default |
| --- | --- | --- |
| `--problem` | `salbp-i`, `bpp-p`, or `bpp-gp` | required |
| `--instance` | Input assembly-line-balancing-format `.txt` file | required in single-instance mode |
| `--graph` | Labeled `.graph` file | required for single-instance BPP-GP only |
| `--time-limit` | Wall-clock limit per solve | `300` for single instances and external batches; bundled batches use the schedule below |
| `--memory-limit-mb` | Global memory cap used by BBR accounting, in MiB | `24576` MiB (24 GiB) |
| `--output-dir` | Directory for `<PROBLEM>_Results.csv` and `solutions/` | `results/<problem>` |

The same executable provides a sequential, resume-safe batch mode used by the supplied launchers:

| Batch option | Meaning | Default |
| --- | --- | --- |
| `--batch` | Enable sequential batch selection and execution | off |
| `--input` | Input `.txt` file or directory | complete benchmark set for the problem |
| `--graph-dir` | BPP-GP `.graph` file or directory | both included labeled-graph sets |

When `--time-limit` is omitted, bundled benchmark batches use the same per-instance limits as the paper:

| Problem and bundled benchmark set | Time limit |
| --- | ---: |
| SALBP-I: Scholl; Otto `n_0100`, `n_1000` | 350 seconds |
| SALBP-I: Otto `n_0020`, `n_0050`, `n_0050_permuted` | 1000 seconds |
| SALBP-I: Otto `n_0250`, `n_0500`, `n_0750` | 75 seconds |
| BPP-P: Scholl; Otto `n_0100` | 1000 seconds |
| BPP-P: all other bundled Otto sizes | 75 seconds |
| BPP-GP: both graph sets at every bundled size | 75 seconds |

Batch inputs outside the bundled benchmark directories retain the 300-second product default. Supplying `--time-limit` explicitly overrides the schedule uniformly for every selected instance.

The deterministic seed is fixed to 1, and solver execution is strictly single-threaded. A time- or memory-limited run still writes its validated incumbent and valid lower bound; only `status=OPTIMAL` certifies optimality. Every successful run appends one compact summary row in comma-separated values (CSV) format to `SALBP-I_Results.csv`, `BPP-P_Results.csv`, or `BPP-GP_Results.csv`, according to the selected problem. The corresponding `.sol` file below `solutions/` contains only the bin assignment, avoiding duplicated data in large benchmark sets. Routine runs create no log files; if an instance raises an error, the batch launcher appends its path, resource limit, error category, and diagnostic to `errors.log` in the output directory.

#### Output Files

Each problem-specific results CSV contains one row per run with the following fixed schema:

| Column | Meaning |
| --- | --- |
| `instance_set` | Benchmark-set directory relative to `data/instances` or `data/bpp-gp-graphs`, such as `otto/n_0020` or `separation-03/n_0020` |
| `instance` | Input filename without the `.txt` suffix |
| `n` | Number of items |
| `time_limit_seconds` | Configured global wall-clock limit in seconds |
| `memory_limit_mb` | Configured global BBR memory limit in MiB |
| `status` | Termination status; only `OPTIMAL` certifies optimality |
| `opt` | `1` for `OPTIMAL`; `0` otherwise |
| `lower_bound` | Valid objective lower bound |
| `upper_bound` | Objective value of the saved assignment |
| `time_seconds` | Total wall-clock solution time in seconds |
| `bbr_peak_memory_bytes` | Peak memory tracked by BBR, not whole-process resident memory |
| `bbr_states_created` | Total BBR states created across all BBR phases |

The problem is identified by the problem-specific output directory and CSV filename. Capacity is stored in the input file, the relative gap is `(upper_bound - lower_bound) / upper_bound`, and the assignment path is deterministically `solutions/<instance-set>/<instance>.sol`; these derived fields are not duplicated in the CSV. CSV fields containing commas, quotes, or line breaks use standard double-quote escaping. **PrecPack** refuses to append to an existing problem-specific results CSV with a different header, preventing rows with incompatible schemas from being mixed. Use a separate output directory for each problem and independent configuration.

Before solving, both single-instance and batch modes check the deterministic `.sol` path. A nonempty file is treated as an existing result and the instance is skipped without comparing time, memory, or Gurobi settings. Use a different output directory for an independent configuration or rerun. **PrecPack** locks the problem-specific results CSV itself while writing, so concurrent writers are rejected without creating a standalone `.precpack.lock` file.

The fixed public BBR interface reports three result statuses:

| Status | Meaning |
| --- | --- |
| `OPTIMAL` | The exact proof is complete and the valid lower bound equals the saved upper bound. |
| `TIME_LIMIT` | The global wall-clock limit was reached before the proof completed. |
| `MEMORY_LIMIT` | The global BBR accounted-memory limit was reached before the proof completed. |

Both limited statuses retain a validated feasible assignment and a valid lower bound. They do not certify optimality.

Each `.sol` file contains one line for every bin position from 1 through the saved upper bound. Generalized precedence constraints may require an unused intermediate position, which is written as an empty line after the colon:

```text
Bin 1: 1 4 7
Bin 2:
Bin 3: 2 3 5 6
```

Bin numbers and item numbers are one-based, and item numbers refer to the original input order. Item weights, total bin weights, bounds, status, resource limits, and runtime statistics are omitted because they can be recovered from the input files and the matching CSV row.

Examples from the repository root:

```bash
./build/precpack \
  --problem bpp-p \
  --instance data/instances/otto/n_0020/instance_n=20_1.txt

./build/precpack \
  --problem salbp-i \
  --instance data/instances/scholl/Bowman/Bowman_c20.txt \
  --time-limit 60 \
  --output-dir results/salbp-i

./build/precpack \
  --problem bpp-gp \
  --instance data/instances/otto/n_0020/instance_n=20_1.txt \
  --graph data/bpp-gp-graphs/separation-03/n_0020/instance_n=20_1.graph
```

Use `build\Release\precpack.exe` for a standard Windows build.

### Batch Runs

Batch selection and recovery are implemented by the same compiled `precpack` executable. The platform launchers execute instances sequentially and skip an instance whenever its deterministically named `.sol` file already exists and is nonempty. Resource settings are intentionally not compared during recovery, so independent configurations must use separate output directories.

Use `./code/scripts/run_batch.sh` on macOS or Linux and `code\scripts\run_batch.bat` on Windows Command Prompt. Both launch the same C++ batch implementation with the same arguments and bundled benchmark schedule. Relative `--input`, `--graph-dir`, and `--output-dir` paths are resolved from the directory in which the launcher is called; the default benchmark data remain located from the repository root. **PrecPack** prevents two processes from writing the same output directory concurrently.

Run all bundled instances for one problem by specifying only its problem type. On macOS or Linux:

```bash
./code/scripts/run_batch.sh --problem salbp-i
./code/scripts/run_batch.sh --problem bpp-p
./code/scripts/run_batch.sh --problem bpp-gp
```

On Windows Command Prompt:

```bat
code\scripts\run_batch.bat --problem salbp-i
code\scripts\run_batch.bat --problem bpp-p
code\scripts\run_batch.bat --problem bpp-gp
```

Each command uses the documented per-instance time-limit schedule and writes to `results/<problem>/` by default. The complete BPP-GP batch pairs every bundled Otto base instance with both labeled-graph sets.

Run a small BPP-P directory:

```bash
./code/scripts/run_batch.sh \
  --problem bpp-p \
  --input data/instances/otto/n_0020 \
  --time-limit 60
```

Run the Scholl SALBP-I benchmark set:

```bash
./code/scripts/run_batch.sh \
  --problem salbp-i \
  --input data/instances/scholl \
  --output-dir results/salbp-i-scholl
```

Run one BPP-GP benchmark set and instance size:

```bash
./code/scripts/run_batch.sh \
  --problem bpp-gp \
  --input data/instances/otto \
  --graph-dir data/bpp-gp-graphs/separation-01/n_0020
```

A complete batch can take a long time; use a benchmark subset or an instance-size directory for an initial check.

## Benchmark Instances

**PrecPack** includes the Otto and Scholl benchmark instance sets together with two BPP-GP labeled-graph sets under [`data/`](data/). See [`data/README.md`](data/README.md) for the directory structure, exact instance counts, file formats, pairing rules, sources, and acknowledgements.

## Contact

If you have any questions or suggestions, or if you encounter an issue, please feel free to contact Sunkanghong Wang at [wskh0929@gmail.com](mailto:wskh0929@gmail.com).

## License

**PrecPack** source code is released under the [MIT License](LICENSE).

Copyright (c) 2026 Sunkanghong Wang.

Gurobi is governed by its own license. The included benchmark instances retain the source attribution and acknowledgements given in the [data documentation](data/README.md).
