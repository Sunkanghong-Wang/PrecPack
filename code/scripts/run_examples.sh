#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
repository_root="$(cd -- "${script_dir}/../.." && pwd -P)"
caller_directory="$(pwd -P)"
output_root="${repository_root}/results/examples"

usage() {
    printf 'Usage: %s [--output-dir DIR]\n' "${0##*/}"
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

base_instance="${repository_root}/data/instances/otto/n_0020/instance_n=20_1.txt"
"${script_dir}/run_batch.sh" \
    --problem bpp-p \
    --input "$base_instance" \
    --time-limit 10 \
    --output-dir "${output_root}/bpp-p"
"${script_dir}/run_batch.sh" \
    --problem salbp-i \
    --input "${repository_root}/data/instances/scholl/Bowman/Bowman_c20.txt" \
    --time-limit 10 \
    --output-dir "${output_root}/salbp-i"
"${script_dir}/run_batch.sh" \
    --problem bpp-gp \
    --input "$base_instance" \
    --graph-dir "${repository_root}/data/bpp-gp-graphs/separation-03/n_0020/instance_n=20_1.graph" \
    --time-limit 10 \
    --output-dir "${output_root}/bpp-gp"
