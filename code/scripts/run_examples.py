#!/usr/bin/env python3
"""Run one bundled instance from each supported problem."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from script_utils import REPOSITORY_ROOT, controlled_environment, find_solver

DATA = REPOSITORY_ROOT / "data" / "instances"
BPP_GP_GRAPHS = REPOSITORY_ROOT / "data" / "bpp-gp-graphs"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="exact BBR workers; -1 uses the available CPU threads",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPOSITORY_ROOT / "results" / "examples",
        help="output directory (default: results/examples)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.threads == 0 or args.threads < -1:
        print("error: threads must be -1 or positive", file=sys.stderr)
        return 2
    try:
        binary = find_solver()
        output_root = args.output_dir.resolve()
        base = DATA / "otto" / "n_0020" / "instance_n=20_1.txt"
        cases = (
            ("bpp-p", base, None),
            (
                "salbp-i",
                DATA
                / "scholl269"
                / "Bowman"
                / "Bowman_c20.txt",
                None,
            ),
            (
                "bpp-gp",
                base,
                BPP_GP_GRAPHS
                / "separation-03"
                / "n_0020"
                / "instance_n=20_1.graph",
            ),
        )
        for problem, instance, graph in cases:
            command = [
                str(binary),
                "--problem",
                problem,
                "--instance",
                str(instance),
                "--time-limit",
                "10",
                "--threads",
                str(args.threads),
                "--output-dir",
                str(output_root / problem),
            ]
            if graph is not None:
                command.extend(("--graph", str(graph)))
            print("+", " ".join(command), flush=True)
            subprocess.run(
                command,
                cwd=REPOSITORY_ROOT,
                env=controlled_environment(),
                check=True,
            )
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
