#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
repository_root="$(cd -- "${script_dir}/../.." && pwd -P)"
caller_directory="$(pwd -P)"
output_root="${repository_root}/results/parallel-experiments"
check_only=false

usage() {
    printf 'Usage: %s [--output-dir DIR] [--check-only]\n' "${0##*/}"
}

while (($# > 0)); do
    case "$1" in
        --output-dir)
            if (($# < 2)) || [[ -z "$2" ]]; then
                printf 'ERROR: --output-dir requires a value.\n' >&2
                exit 2
            fi
            output_root="$2"
            shift 2
            ;;
        --check-only)
            check_only=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            printf 'ERROR: Unknown argument: %s\n' "$1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [[ "$output_root" != /* ]]; then
    output_root="${caller_directory}/${output_root}"
fi

export PRECPACK_REQUIRE_GUROBI_RUNTIME=1
if [[ "$check_only" == false ]]; then
    help_text="$("${script_dir}/run_batch.sh" --help)"
    if [[ "$help_text" != *"Optional Gurobi root strengthening: enabled"* ]]; then
        printf 'ERROR: The controlled experiment requires a Gurobi-enabled build.\n' >&2
        printf 'Rebuild with code/scripts/build.sh --gurobi on.\n' >&2
        exit 1
    fi
fi

item_directory="${repository_root}/data/instances/otto/n_0100"
graph_01="${repository_root}/data/bpp-gp-graphs/separation-01/n_0100"
graph_03="${repository_root}/data/bpp-gp-graphs/separation-03/n_0100"
shopt -s nullglob
item_files=("${item_directory}"/*.txt)
graph_01_files=("${graph_01}"/*.graph)
graph_03_files=("${graph_03}"/*.graph)
shopt -u nullglob
if ((${#item_files[@]} != 525 || ${#graph_01_files[@]} != 525 ||
     ${#graph_03_files[@]} != 525)); then
    printf 'ERROR: Each n=100 experiment data set must contain 525 files.\n' >&2
    exit 1
fi

printf 'PrecPack controlled parallel experiments\n'
printf '  datasets       : 4\n'
printf '  threads        : 2, 4, 8\n'
printf '  configurations : 12\n'
printf '  total runs     : 6300\n'
printf '  output         : %s\n' "$output_root"

run_configuration() {
    local slug="$1"
    local problem="$2"
    local time_limit="$3"
    local threads="$4"
    local graph_directory="${5:-}"
    local arguments=(
        --problem "$problem"
        --input "$item_directory"
        --time-limit "$time_limit"
        --memory-limit-mb 24576
        --threads "$threads"
        --output-dir "${output_root}/${slug}/threads-${threads}"
    )
    if [[ -n "$graph_directory" ]]; then
        arguments+=(--graph-dir "$graph_directory")
    fi
    if [[ "$check_only" == true ]]; then
        arguments+=(--check-only)
    fi
    printf '\n[%s, threads=%s]\n' "$slug" "$threads"
    "${script_dir}/run_batch.sh" "${arguments[@]}"
}

for threads in 2 4 8; do
    run_configuration salbp-i salbp-i 350 "$threads"
    run_configuration bpp-p bpp-p 1000 "$threads"
    run_configuration bpp-gp-separation-01 bpp-gp 75 "$threads" "$graph_01"
    run_configuration bpp-gp-separation-03 bpp-gp 75 "$threads" "$graph_03"
done

if [[ "$check_only" == true ]]; then
    printf '\nParallel experiment validation completed.\n'
else
    printf '\nParallel experiment runs completed.\n'
fi
