#!/usr/bin/env bash

set -euo pipefail

caller_directory="$(pwd -P)"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
repository_root="$(cd -- "${script_dir}/../.." && pwd -P)"
binary="${repository_root}/build/precpack"
if [[ ! -x "$binary" ]]; then
    binary="${repository_root}/build/Release/precpack"
fi
if [[ ! -x "$binary" ]]; then
    printf 'ERROR: PrecPack is not built. Run code/scripts/build.sh first.\n' >&2
    exit 1
fi

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
if [[ -n "${GUROBI_HOME:-}" ]]; then
    if [[ "$(uname -s)" == "Darwin" ]]; then
        export DYLD_LIBRARY_PATH="${GUROBI_HOME}/lib${DYLD_LIBRARY_PATH:+:${DYLD_LIBRARY_PATH}}"
    else
        export LD_LIBRARY_PATH="${GUROBI_HOME}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
    fi
fi

export PRECPACK_REPOSITORY_ROOT="$repository_root"
export PRECPACK_CALLER_DIRECTORY="$caller_directory"
exec "$binary" --batch "$@"
