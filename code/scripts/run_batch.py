#!/usr/bin/env python3
"""Run a resume-safe PrecPack batch, one instance at a time."""

from __future__ import annotations

import argparse
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from script_utils import REPOSITORY_ROOT, controlled_environment, find_solver

DATA = REPOSITORY_ROOT / "data" / "instances"
BPP_GP_GRAPHS = REPOSITORY_ROOT / "data" / "bpp-gp-graphs"
OTTO_BASE_DIRECTORIES = (
    "n_0020",
    "n_0050",
    "n_0100",
    "n_0250",
    "n_0500",
    "n_0750",
    "n_1000",
)


@dataclass(frozen=True)
class Case:
    instance: Path
    graph: Path | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--problem", required=True, choices=("salbp-i", "bpp-p", "bpp-gp")
    )
    parser.add_argument(
        "--input",
        type=Path,
        help="ALB-format .txt file or directory; defaults to the bundled collection",
    )
    parser.add_argument(
        "--graph-dir",
        type=Path,
        help=".graph file or directory for bpp-gp; defaults to both bundled separation sets",
    )
    parser.add_argument("--time-limit", type=float, default=300.0)
    parser.add_argument("--memory-limit-mb", type=int, default=24 * 1024)
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="exact BBR workers; -1 uses the available CPU threads",
    )
    parser.add_argument("--output-dir", type=Path)
    return parser.parse_args()


def files_below(path: Path, suffix: str) -> list[Path]:
    if path.is_file():
        if path.suffix.lower() != suffix.lower():
            raise ValueError(f"expected a {suffix} file: {path}")
        return [path.resolve()]
    if not path.is_dir():
        raise ValueError(f"input does not exist: {path}")
    return sorted(
        file.resolve()
        for file in path.rglob("*")
        if file.is_file() and file.suffix.lower() == suffix.lower()
    )


def collect_cases(
    problem: str,
    input_path: Path | None,
    graph_path: Path | None,
) -> list[Case]:
    if problem == "salbp-i":
        if graph_path is not None:
            raise ValueError("--graph-dir is only valid for bpp-gp")
        roots = (
            (input_path,)
            if input_path is not None
            else (DATA / "otto", DATA / "scholl269")
        )
        return [Case(file) for root in roots for file in files_below(root, ".txt")]
    if problem == "bpp-p":
        if graph_path is not None:
            raise ValueError("--graph-dir is only valid for bpp-gp")
        roots = (
            (input_path,)
            if input_path is not None
            else tuple(DATA / "otto" / name for name in OTTO_BASE_DIRECTORIES)
            + (DATA / "scholl269",)
        )
        return [Case(file) for root in roots for file in files_below(root, ".txt")]

    item_root = input_path or DATA / "otto"
    graph_root = graph_path or BPP_GP_GRAPHS
    graphs = files_below(graph_root, ".graph")
    cases: list[Case] = []
    for graph in graphs:
        n_directory = graph.parent.name
        instance = item_root / n_directory / f"{graph.stem}.txt"
        if not instance.is_file():
            raise ValueError(f"no matching instance file for {graph}: {instance}")
        cases.append(Case(instance.resolve(), graph))
    return cases


def sanitized_filename(value: str) -> str:
    result = "".join(
        (
            character
            if character.isascii() and (character.isalnum() or character in "-_.")
            else "_"
        )
        for character in value
    )
    return result or "instance"


def output_key(problem: str, case: Case) -> str:
    key = case.instance.stem
    if case.graph is not None:
        key += "__" + case.graph.parent.parent.name
    return f"{problem}__{sanitized_filename(key)}.sol"


def main() -> int:
    args = parse_args()
    if (
        args.time_limit <= 0
        or args.memory_limit_mb <= 0
        or args.threads == 0
        or args.threads < -1
    ):
        print(
            "error: time and memory limits must be positive; "
            "threads must be -1 or positive",
            file=sys.stderr,
        )
        return 2
    try:
        binary = find_solver()
        cases = collect_cases(args.problem, args.input, args.graph_dir)
        if not cases:
            raise ValueError("no instances selected")
        output_keys = [output_key(args.problem, case) for case in cases]
        if len(output_keys) != len(set(output_keys)):
            raise ValueError(
                "selected files contain duplicate output names; "
                "run the conflicting directories separately"
            )
        output_dir = (
            args.output_dir or REPOSITORY_ROOT / "results" / args.problem
        )
        failures = 0
        for index, case in enumerate(cases, start=1):
            solution = output_dir / "solutions" / output_key(args.problem, case)
            if solution.is_file() and solution.stat().st_size > 0:
                print(f"[{index}/{len(cases)}] skip {case.instance.name}")
                continue
            command = [
                str(binary),
                "--problem",
                args.problem,
                "--instance",
                str(case.instance),
                "--time-limit",
                str(args.time_limit),
                "--memory-limit-mb",
                str(args.memory_limit_mb),
                "--threads",
                str(args.threads),
                "--output-dir",
                str(output_dir),
            ]
            if case.graph is not None:
                command.extend(("--graph", str(case.graph)))
            print(f"[{index}/{len(cases)}] {case.instance.name}", flush=True)
            completed = subprocess.run(
                command,
                cwd=REPOSITORY_ROOT,
                env=controlled_environment(),
                check=False,
            )
            if completed.returncode != 0:
                failures += 1
        return 1 if failures else 0
    except (OSError, RuntimeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
