#!/usr/bin/env bash

set -uo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
repository_root="$(cd -- "${script_dir}/../.." && pwd -P)"
code_directory="${repository_root}/code"
build_directory="${repository_root}/build"
gurobi_mode="AUTO"

usage() {
    printf 'Usage: %s [--gurobi auto|on|off]\n' "${0##*/}"
}

while (($# > 0)); do
    case "$1" in
        --gurobi)
            if (($# < 2)); then
                printf 'ERROR: --gurobi requires one of: auto, on, off.\n' >&2
                usage >&2
                exit 2
            fi
            case "$2" in
                [Aa][Uu][Tt][Oo]) gurobi_mode="AUTO" ;;
                [Oo][Nn]) gurobi_mode="ON" ;;
                [Oo][Ff][Ff]) gurobi_mode="OFF" ;;
                *)
                    printf 'ERROR: Invalid --gurobi value: %s\n' "$2" >&2
                    usage >&2
                    exit 2
                    ;;
            esac
            shift 2
            ;;
        -h|--help)
            usage
            printf '\nConfigure, build, and test PrecPack in Release mode.\n'
            printf 'This launcher calls CMake directly and does not require Python.\n'
            exit 0
            ;;
        *)
            printf 'ERROR: Unknown argument: %s\n' "$1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

printf 'PrecPack Release build\n'
printf '  Source       : %s\n' "$code_directory"
printf '  Build        : %s\n' "$build_directory"
printf '  Gurobi mode  : %s\n' "$gurobi_mode"

cmake_command=""
if [[ -n "${CMAKE_BIN:-}" ]]; then
    cmake_command="$CMAKE_BIN"
elif cmake_path="$(command -v cmake 2>/dev/null)"; then
    cmake_command="$cmake_path"
elif [[ "$(uname -s)" == "Darwin" ]]; then
    cmake_candidates=(
        "/Applications/CMake.app/Contents/bin/cmake"
        "/opt/homebrew/bin/cmake"
        "/usr/local/bin/cmake"
        "/opt/local/bin/cmake"
    )
    if [[ "$(uname -m)" == "arm64" ]]; then
        cmake_candidates+=(
            "/Applications/CLion.app/Contents/bin/cmake/mac/aarch64/bin/cmake"
            "/Applications/CLion.app/Contents/bin/cmake/mac/x64/bin/cmake"
        )
    else
        cmake_candidates+=(
            "/Applications/CLion.app/Contents/bin/cmake/mac/x64/bin/cmake"
            "/Applications/CLion.app/Contents/bin/cmake/mac/aarch64/bin/cmake"
        )
    fi
    for candidate in "${cmake_candidates[@]}"; do
        if [[ -x "$candidate" ]]; then
            cmake_command="$candidate"
            break
        fi
    done
fi

if [[ -z "$cmake_command" ]]; then
    printf '\nERROR: CMake was not found.\n' >&2
    printf 'Install CMake 3.20 or newer and make cmake available on PATH,\n' >&2
    printf 'or set CMAKE_BIN to the full path of the CMake executable.\n' >&2
    exit 1
fi
if [[ ! -x "$cmake_command" ]]; then
    printf '\nERROR: The selected CMake executable is not runnable:\n' >&2
    printf '  %s\n' "$cmake_command" >&2
    printf 'Correct CMAKE_BIN or install CMake 3.20 or newer.\n' >&2
    exit 1
fi

printf '  CMake        : %s\n' "$cmake_command"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

gurobi_root="${GUROBI_HOME:-}"
if [[ -n "$gurobi_root" ]]; then
    if [[ "$(uname -s)" == "Darwin" ]]; then
        export DYLD_LIBRARY_PATH="${gurobi_root}/lib${DYLD_LIBRARY_PATH:+:${DYLD_LIBRARY_PATH}}"
    else
        export LD_LIBRARY_PATH="${gurobi_root}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
    fi
fi

print_command() {
    printf '+'
    printf ' %q' "$@"
    printf '\n'
}

run_step() {
    local label="$1"
    shift
    printf '\n%s\n' "$label"
    print_command "$@"
    "$@"
    local status=$?
    if ((status != 0)); then
        printf '\nERROR: %s failed. Review the messages above.\n' "$label" >&2
    fi
    return "$status"
}

configure_command=(
    "$cmake_command"
    -S "$code_directory"
    -B "$build_directory"
    -DCMAKE_BUILD_TYPE=Release
    -DBUILD_TESTING=ON
    -DPRECPACK_NATIVE_OPTIMIZATION=ON
    -DPRECPACK_ENABLE_IPO=ON
    -DPRECPACK_ENABLE_SANITIZERS=OFF
    "-DPRECPACK_GUROBI=${gurobi_mode}"
    "-DGUROBI_ROOT=${gurobi_root}"
)
build_command=(
    "$cmake_command"
    --build "$build_directory"
    --config Release
    --target precpack
    --parallel
)
test_command=(
    "$cmake_command"
    --build "$build_directory"
    --config Release
    --target check
)

run_step "[1/3] Configuring PrecPack..." "${configure_command[@]}" || exit $?
run_step "[2/3] Building PrecPack..." "${build_command[@]}" || exit $?
run_step "[3/3] Running the test suite..." "${test_command[@]}" || exit $?

printf '\nPrecPack was built and tested successfully.\n'
if [[ -x "${build_directory}/precpack" ]]; then
    printf '  Executable: %s\n' "${build_directory}/precpack"
elif [[ -x "${build_directory}/Release/precpack" ]]; then
    printf '  Executable: %s\n' "${build_directory}/Release/precpack"
fi
