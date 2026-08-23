# PrecPack

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
  <a href="#highlights"><img alt="Solver: exact" src="https://img.shields.io/badge/solver-exact-f59f00"></a>
  <a href="#build"><img alt="Platforms: macOS, Linux, and Windows" src="https://img.shields.io/badge/platforms-macOS%20%7C%20Linux%20%7C%20Windows-6c63ff"></a>
  <a href="LICENSE"><img alt="License: MIT" src="https://img.shields.io/badge/license-MIT-green"></a>
</p>

<p align="center">
  <a href="#cite">Cite</a> |
  <a href="#overview">Overview</a> |
  <a href="#highlights">Highlights</a> |
  <a href="#requirements">Requirements</a> |
  <a href="#quick-start">Quick Start</a> |
  <a href="#benchmark-instances">Instances</a> |
  <a href="#contact">Contact</a> |
  <a href="#license">License</a>
</p>

## Cite

To cite the contents of this repository, please cite both the paper and the software repository using their respective DOIs.

- Paper DOI: `10.1287/ijoc.XXXX.YYYY`
- Repository DOI: `10.1287/ijoc.XXXX.YYYY.cd`

Below is the BibTeX for citing this snapshot of the repository.

```bibtex
@misc{Wang2026PrecPack,
  author    = {Wang, Sunkanghong and You, Zhengzhong and Baldacci, Roberto and
               Mo, Baichuan and Wei, Lijun and Xu, Zhou},
  publisher = {INFORMS Journal on Computing},
  title     = {{PrecPack: An Efficient Open-Source Exact Solver for Bin Packing
                with Generalized Precedence Constraints}},
  year      = {2026},
  doi       = {10.1287/ijoc.XXXX.YYYY.cd},
  url       = {https://github.com/INFORMSJoC/XXXX.YYYY},
  note      = {Available for download at https://github.com/INFORMSJoC/XXXX.YYYY},
}
```

The provisional software metadata is also available in [`CITATION.cff`](CITATION.cff), and the complete author list is recorded in [`AUTHORS`](AUTHORS).

## Overview

PrecPack is a single, efficient open-source exact solver for the Bin Packing Problem with Generalized Precedence Constraints (BPP-GP), which encompasses the following two classical problems as special cases:

- Simple Assembly Line Balancing Problem type I (SALBP-I);
- Bin Packing Problem with Precedence Constraints (BPP-P).

This repository is intentionally software-first. It contains the solver source code, regression tests, benchmark instances, and build/run documentation. It does not contain computational results, experiment logs, profiler output, or development artifacts.

```text
PrecPack/
├── .github/
│   └── workflows/build.yml    # Gurobi-free three-platform CI
├── assets/                     # README artwork
├── code/
│   ├── CMakeLists.txt          # Build definition
│   ├── cmake/                  # Optional Gurobi discovery
│   ├── include/precpack/       # Shared declarations and data types
│   ├── scripts/                # Build/run tools and regression tests
│   └── src/
│       ├── bbr.cpp             # Serial and shared-memory exact BBR
│       ├── initial_bounds.cpp  # Preprocessing, bounds, and incumbents
│       ├── root_column_generation.cpp  # Position-free root bound
│       ├── bin_indexed_root_bound.cpp  # Bin-indexed root bound
│       ├── bin_packing_bound.cpp       # Exact residual relaxation
│       ├── dff.cpp                     # Dual-feasible functions
│       └── ...                         # CLI, I/O, profiles, and orchestration
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

All implementation and developer-facing tooling lives under `code/`. The components are separated by algorithmic responsibility, while PrecPack remains one focused solver rather than a plug-in framework.

## Highlights

| Goal | PrecPack implementation |
| --- | --- |
| One generalized exact state | BBR stores the assigned-item set together with a nested residual-distance profile. The profile becomes empty for SALBP-I and BPP-P. |
| Complete branching | Precedence-closed maximal-load generation is combined with a unique forced empty-bin transition when no nonempty load is feasible. |
| Proof-aware bounds | Integer capacity, path, machine, closure, DFF, and completed exact-relaxation bounds are ordered by cost around collision-safe state memory. |
| Safe auxiliary work | Preliminary searches and root modules may improve an incumbent or certified bound, but unfinished work cannot lend memory or an uncertified value to the final proof. |
| Commercial-solver-free exactness | PrecPack can be compiled without Gurobi or any other commercial optimizer. The complete BBR search still proves optimality; only Gurobi-dependent root strengthening and small-instance reference checks are omitted. |
| Validated outcomes | Every incumbent is checked against the original instance; limited runs retain a certified lower bound and a validated feasible assignment. |
| Reproducible defaults | One fixed production profile per problem; algorithmic ablation switches are not exposed through the public command line. |
| Optional shared-memory search | The default serial BBR path is preserved; an explicit worker count enables exact shallow-subtree parallelism shared by all three problem profiles. |
| Portable workflow | CMake definitions and convenience scripts support macOS, Linux, and Windows. |
| Documented data | The included benchmark sets have format and source documentation. |

## Requirements

The requirements below are derived from the checked-in source and build files, not from the versions installed on a development machine.

| Dependency | Required or supported version | Basis |
| --- | --- | --- |
| Operating system | 64-bit macOS, Linux, or Windows | Supported by the supplied CMake and launcher workflows. |
| C++ compiler and standard library | Full C++20 support | PrecPack uses C++20 library features, including `<bit>` operations and associative-container `contains()`. |
| CMake | 3.20 or newer | Enforced by `cmake_minimum_required(VERSION 3.20)`. |
| Gurobi Optimizer | Optional; 9.1 through 13.x | Used only for optional root strengthening and small-instance reference checks. The build reads the installed version from `gurobi_c.h`; 9.1.1 and 13.0.2 are the currently tested versions. |

PrecPack has no mandatory mathematical-programming-solver dependency. A C++20 compiler, CMake, and the operating-system thread library are sufficient to build and run the exact solver and its supplied launchers.

### Optional Gurobi acceleration

The CMake option `PRECPACK_GUROBI` selects one of three build modes:

| Value | Behavior |
| --- | --- |
| `AUTO` | Default. Enable Gurobi support when a compatible installation is found; otherwise build the commercial-solver-free exact solver. |
| `OFF` | Exclude all Gurobi-dependent source files and link no Gurobi library. Use this mode for a guaranteed commercial-solver-free binary. |
| `ON` | Require Gurobi and stop configuration with an error if it cannot be found. |

Gurobi 13.0 is not required. PrecPack supports versions 9.1 through 13.x. It conditionally handles the `WORK_LIMIT` status introduced in Gurobi 9.5, and derives the matching core-library name from the version macros in `gurobi_c.h`. The finder also keeps headers and libraries within the same installation root. Versions 9.1.1 and 13.0.2 are tested; intermediate releases rely on the same public API but remain untested on the current development machine. Older Gurobi C++ libraries must be binary-compatible with the compiler used to build PrecPack.

Without Gurobi, PrecPack disables the optional root-strengthening models and the compact-MIP oracle used by regression tests. The public SALBP-I, BPP-P, and BPP-GP interface continues to run the same complete BBR search, return validated incumbents and certified lower bounds under limits, and report `OPTIMAL` only after the BBR proof is complete. Thus Gurobi can affect running time and search statistics, but it is not required for correctness or exactness.

If a Gurobi-enabled binary is run without an accessible license, the optional root module is abandoned safely and PrecPack continues with BBR. For predictable deployment on machines without a commercial license, prefer an `OFF` build.

To enable the optional backend, set `GUROBI_HOME` to the platform directory containing Gurobi's `include`, `lib`, and, on Windows, `bin` directories:

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

The macOS, Linux, and Windows build launchers call CMake directly. They print the CMake configuration, compilation, and test phases, and report a specific error when CMake is unavailable or a phase fails. On macOS, `build.sh` also detects CMake installed with CMake.app, Homebrew, MacPorts, or CLion when it is not on `PATH`. Set `CMAKE_BIN` to an executable path to override automatic detection. None of the build or run launchers requires Python.

Omit `--gurobi off` to use the default `AUTO` detection. To require the optional Gurobi backend, set `GUROBI_HOME` and use `--gurobi on`:

```bash
export GUROBI_HOME=/path/to/gurobi/platform
./code/scripts/build.sh --gurobi on
```

```bat
set "GUROBI_HOME=C:\path\to\gurobi\win64"
code\scripts\build.bat --gurobi on
```

These launchers configure a Release build, compile PrecPack, and run the test suite. Every configuration runs the commercial-solver-free regressions for assignment validation, independent brute-force exactness oracles, DFF and BINLB bounds, initialization, resource-limit termination, the public interface, and serial/parallel agreement. When Gurobi is available, root-bound and compact-MIP oracle tests run in addition. Compilation may use multiple jobs; solver execution remains single-threaded unless `--threads` is explicitly set.

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

The public command line contains only problem semantics, input locations, resource limits, and the output location:

| Option | Meaning | Default |
| --- | --- | --- |
| `--problem` | `salbp-i`, `bpp-p`, or `bpp-gp` | required |
| `--instance` | Input ALB-format `.txt` file | required in single-instance mode |
| `--graph` | Labeled `.graph` file | required for single-instance BPP-GP only |
| `--time-limit` | Global wall-clock limit in seconds | `300` |
| `--memory-limit-mb` | Global memory cap used by BBR accounting | `24576` |
| `--threads` | Exact BBR workers; `-1` uses the available hardware concurrency | `1` |
| `--output-dir` | Directory for `<PROBLEM>_Results.csv` and `solutions/` | `results` for one instance; `results/<problem>` for a batch |

The same executable provides a sequential, resume-safe batch mode used by the supplied launchers:

| Batch option | Meaning | Default |
| --- | --- | --- |
| `--batch` | Enable sequential batch selection and execution | off |
| `--input` | Input `.txt` file or directory | complete benchmark set for the problem |
| `--graph-dir` | BPP-GP `.graph` file or directory | both included labeled-graph sets |
| `--check-only` | Validate the selection and compatible existing output without solving | off |

The deterministic seed is fixed to 1, and the production remembered-state cap is fixed to 60,000,000. In parallel runs this is one global cap, not 60,000,000 states per worker. `--threads -1` uses `std::thread::hardware_concurrency()` and falls back to one worker if the platform cannot report it. A time-, state-, or memory-limited run still writes its validated incumbent and certified lower bound; only `status=OPTIMAL` certifies optimality. Every run appends its metadata and summary statistics to `SALBP-I_Results.csv`, `BPP-P_Results.csv`, or `BPP-GP_Results.csv`, according to the selected problem. The corresponding `.sol` file below `solutions/` contains only the bin assignment, avoiding duplicated data in large benchmark sets.

Every writable batch also appends a flushed event stream to `logs/batch-events.log`. A `START` record is persisted before each solve, followed by `SUCCESS` or `ERROR`; an unmatched final `START` therefore identifies the active instance if the process is killed or crashes before C++ can report an exception. A caught solver, input, Gurobi, or output exception is additionally appended to `logs/<problem>__<instance_key>.failure.log` with the complete resource profile and message. If a later resume solves that instance, the same file retains the failure history and receives a `RECOVERED` record.

#### Output Files

Each problem-specific results CSV contains one row per run with the following fixed schema:

| Column | Meaning |
| --- | --- |
| `instance_key` | File-safe identifier consisting of the input stem and a deterministic fingerprint of the normalized input paths; used to name the solution |
| `problem` | `SALBP-I`, `BPP-P`, or `BPP-GP` |
| `instance_file` | Input `.txt` path |
| `graph_file` | Input `.graph` path; empty for SALBP-I and BPP-P |
| `n` | Number of items |
| `capacity` | Bin capacity or SALBP-I cycle time |
| `status` | Termination status; only `OPTIMAL` certifies optimality |
| `lower_bound` | Certified objective lower bound |
| `upper_bound` | Objective value of the saved assignment |
| `gap` | `(upper_bound - lower_bound) / upper_bound` |
| `time_seconds` | Total wall-clock solution time in seconds |
| `time_limit_seconds` | Configured global wall-clock limit in seconds |
| `threads` | Resolved number of exact BBR workers |
| `state_limit` | Configured global remembered-state limit |
| `memory_limit_mb` | Configured global BBR memory limit in MiB |
| `bbr_peak_memory_bytes` | Peak memory tracked by BBR, not whole-process resident memory |
| `gurobi_enabled` | `1` if the executable was built with optional Gurobi support, otherwise `0` |
| `gurobi_required` | `1` when the run required an accessible Gurobi runtime and prohibited fallback, otherwise `0` |
| `solution_file` | `.sol` path relative to the output directory |

CSV fields containing commas, quotes, or line breaks use standard double-quote escaping. Input paths below the repository root established by a supplied launcher, or below the working directory when the executable is called directly, are recorded as normalized relative paths. PrecPack refuses to append to an existing problem-specific results CSV with a different header, preventing rows with incompatible schemas from being mixed. In single-instance mode, it also refuses a repeated `instance_key` before solving so that an existing CSV row can never refer to an overwritten `.sol`; use a different output directory for an independent rerun, or use batch mode for profile-checked resume behavior.

Each output directory contains a zero-byte `.precpack.lock` coordination file. PrecPack uses an operating-system file lock on it to serialize writers; the lock is released automatically when the process exits, including after an abnormal termination. The file may remain in place and is not solver output.

The fixed public BBR interface reports four termination statuses:

| Status | Meaning |
| --- | --- |
| `OPTIMAL` | The exact proof is complete and the certified lower bound equals the saved upper bound. |
| `TIME_LIMIT` | The global wall-clock limit was reached before the proof completed. |
| `STATE_LIMIT` | The global remembered-state limit was reached before the proof completed. |
| `MEMORY_LIMIT` | The global BBR accounted-memory limit was reached before the proof completed. |

All three limited statuses retain a validated feasible assignment and a certified lower bound. They do not certify optimality.

Each `.sol` file contains one line for every bin or station position from 1 through the saved upper bound. Generalized separations may require an unused intermediate position, which is written as an empty line after the colon:

```text
Bin 1: 1 4 7
Bin 2:
Bin 3: 2 3 5 6
```

Bin numbers and item numbers are one-based, and item numbers refer to the original input order. Weights, loads, bounds, status, resource limits, and runtime statistics are omitted because they can be recovered from the input files and the matching CSV row.

Examples from the repository root:

```bash
./build/precpack \
  --problem bpp-p \
  --instance data/instances/otto/n_0020/instance_n=20_1.txt

./build/precpack \
  --problem salbp-i \
  --instance data/instances/scholl/Bowman/Bowman_c20.txt \
  --time-limit 60 \
  --threads 4 \
  --output-dir results/salbp-i

./build/precpack \
  --problem bpp-gp \
  --instance data/instances/otto/n_0020/instance_n=20_1.txt \
  --graph data/bpp-gp-graphs/separation-03/n_0020/instance_n=20_1.graph
```

Use `build\Release\precpack.exe` for a standard Windows build.

### Batch Runs

Batch selection and recovery are implemented by the same compiled `precpack` executable. The platform launchers execute instances sequentially and skip an instance only when its CSV row matches the requested problem, resource profile, build capability, and solution path, and the referenced `.sol` file exists and is nonempty. A conflicting profile or a CSV row without its solution file is reported as inconsistent instead of being silently mixed into or skipped within the same output directory.

Use `./code/scripts/run_batch.sh` on macOS or Linux and `code\scripts\run_batch.bat` on Windows Command Prompt. Both launch the same C++ batch implementation with the same arguments. Relative `--input`, `--graph-dir`, and `--output-dir` paths are resolved from the directory in which the launcher is called; the default benchmark data remain located from the repository root. PrecPack prevents two processes from writing the same output directory concurrently.

Run a small BPP-P directory:

```bash
./code/scripts/run_batch.sh \
  --problem bpp-p \
  --input data/instances/otto/n_0020 \
  --time-limit 60 \
  --threads 4
```

Run the Scholl SALBP-I benchmark set:

```bash
./code/scripts/run_batch.sh \
  --problem salbp-i \
  --input data/instances/scholl \
  --output-dir results/salbp-i-scholl
```

Run one BPP-GP separation family and size:

```bash
./code/scripts/run_batch.sh \
  --problem bpp-gp \
  --input data/instances/otto \
  --graph-dir data/bpp-gp-graphs/separation-01/n_0020
```

Omitting `--input` selects the complete benchmark set for the requested problem. A complete batch can take a long time; use a benchmark subset or an instance-size directory for an initial check.

#### Parallel Experiments

`run_parallel.sh` and `run_parallel.bat` run the multithreaded configurations for the controlled experiment reported in the paper's “Parallel Scalability” table. They use the complete 100-item Otto benchmark sets for SALBP-I, BPP-P, and both BPP-GP labeled-graph sets. Each 525-instance set is run sequentially with 2, 4, and 8 solver threads, for 12 configurations and 6,300 runs. The corresponding single-thread results are produced by the problem-level batch runs and are not repeated here. The limits are 350 seconds for SALBP-I, 1,000 seconds for BPP-P, and 75 seconds for each BPP-GP graph set; every configuration uses 24,576 MiB of globally accounted BBR memory and the fixed 60,000,000-state cap.

The experiment uses one frozen Gurobi-enabled Release executable so that the BPP-P and BPP-GP configurations retain the reported `price-and-switch` root policy. It therefore requires an accessible Gurobi license, although the general PrecPack solver remains exact without Gurobi. The experiment launcher enables a strict runtime profile: it verifies the license before solving, aborts on any later Gurobi failure instead of falling back, and records `gurobi_required=1`. Build once and do not rebuild while a batch is being resumed.

macOS or Linux:

```bash
./code/scripts/build.sh --gurobi on
./code/scripts/run_parallel.sh
```

Windows Command Prompt:

```bat
code\scripts\build.bat --gurobi on
code\scripts\run_parallel.bat
```

Results are isolated below `results/parallel-experiments/<problem>/threads-<n>/`, preventing one thread configuration from overwriting another configuration's solution files. The launcher executes only one configuration at a time, restricts third-party numerical libraries to one thread, and verifies each 525-instance selection, the Gurobi runtime, problem, resource profile, CSV uniqueness, and solution-file presence before resuming. `--output-dir DIR` selects another output root. `--check-only` validates all 12 selections and compatible existing output without solving or requiring a Gurobi license, allowing the experiment design to be checked in the commercial-solver-free CI build. After all 6,300 runs complete, speedup and parallel efficiency can be computed by joining these CSV rows with the existing single-thread results on `instance_key` and retaining the common optimally solved subset. The total worst-case configured sequential time budget is 656.25 hours, excluding build and launcher overhead.

## Benchmark Instances

The repository includes the following benchmark sets:

| Problem | Benchmark set | Count |
| --- | --- | ---: |
| BPP-P | Seven Otto instance-size sets with 525 instances each, plus the 269-instance Scholl set | 3,944 |
| BPP-GP | Two labeled-graph sets, each paired with the 3,675 Otto base instances | 7,350 |
| SALBP-I | The Otto and Scholl sets, plus 4,725 additional permutations of the Otto 50-item instances | 8,669 |

The benchmark instances are drawn from or based on the following benchmark sets:

- the public [Scholl 1993 SALBP benchmark set](https://assembly-line-balancing.de/salbp/benchmark-data-sets-1993/);
- the public [Otto et al. 2013 SALBP benchmark set](https://assembly-line-balancing.de/salbp/benchmark-data-sets-2013/) accompanying [“Systematic data generation and test design for solution algorithms on the example of SALBPGen for assembly line balancing”](https://doi.org/10.1016/j.ejor.2012.12.029);
- the BPP-P and BPP-GP benchmark sets described by Kramer, Dell'Amico, and Iori in [“A batching-move iterated local search algorithm for the bin packing problem with generalized precedence constraints”](https://doi.org/10.1080/00207543.2017.1341065).

We gratefully thank the authors and maintainers of these benchmark sets for making them available to the research community. See [`data/README.md`](data/README.md) for the mathematical interpretation, directory map, `.txt` and `.graph` formats, file pairing rules, exact counts, and detailed attribution. The benchmark files contain no PrecPack result rows, solutions, logs, or timing records.

## Contact

If you have any questions, suggestions, or encounter an issue, please feel free to contact Sunkanghong Wang at [wskh0929@gmail.com](mailto:wskh0929@gmail.com).

## License

PrecPack source code is released under the [MIT License](LICENSE).

Copyright (c) 2026 Sunkanghong Wang.

Gurobi is governed by its own license. The included benchmark instances retain the source attribution and acknowledgements given in the [data documentation](data/README.md).
