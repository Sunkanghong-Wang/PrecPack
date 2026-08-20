#!/usr/bin/env python3
"""Run a resume-safe PrecPack batch, one instance at a time."""

from __future__ import annotations

import argparse
import csv
import os
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
PROBLEM_NAMES = {
    "salbp-i": "SALBP-I",
    "bpp-p": "BPP-P",
    "bpp-gp": "BPP-GP",
}
RESULT_COLUMNS = (
    "instance_key",
    "problem",
    "instance_file",
    "graph_file",
    "n",
    "capacity",
    "status",
    "lower_bound",
    "upper_bound",
    "gap",
    "time_seconds",
    "time_limit_seconds",
    "threads",
    "state_limit",
    "memory_limit_mb",
    "bbr_peak_memory_bytes",
    "gurobi_enabled",
    "solution_file",
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


def collect_bpp_gp_cases(item_path: Path, graph_path: Path) -> list[Case]:
    instances = files_below(item_path, ".txt")
    graphs = files_below(graph_path, ".graph")

    if item_path.is_file():
        instance = instances[0]
        matching_graphs = [graph for graph in graphs if graph.stem == instance.stem]
        if not matching_graphs:
            raise ValueError(
                f"no graph with stem {instance.stem!r} below {graph_path}"
            )
        return [Case(instance, graph) for graph in matching_graphs]

    instances_by_stem: dict[str, list[Path]] = {}
    for instance in instances:
        instances_by_stem.setdefault(instance.stem, []).append(instance)

    item_root = item_path.resolve()
    cases: list[Case] = []
    for graph in graphs:
        candidates = instances_by_stem.get(graph.stem, [])
        if not candidates:
            raise ValueError(
                f"no matching instance file for {graph} below {item_path}"
            )

        direct_candidates = (
            item_root / graph.parent.name / f"{graph.stem}.txt",
            item_root / f"{graph.stem}.txt",
        )
        instance = next(
            (
                candidate.resolve()
                for candidate in direct_candidates
                if candidate.is_file() and candidate.resolve() in candidates
            ),
            None,
        )
        if instance is None:
            same_size = [
                candidate
                for candidate in candidates
                if candidate.parent.name == graph.parent.name
            ]
            if len(same_size) == 1:
                instance = same_size[0]
            elif len(candidates) == 1:
                instance = candidates[0]
            else:
                matches = ", ".join(str(candidate) for candidate in candidates)
                raise ValueError(
                    f"ambiguous instance match for {graph}: {matches}"
                )
        cases.append(Case(instance, graph))
    return cases


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

    return collect_bpp_gp_cases(
        input_path or DATA / "otto", graph_path or BPP_GP_GRAPHS
    )


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


def recorded_path(path: Path) -> str:
    normalized = Path(os.path.normpath(path))
    if not normalized.is_absolute():
        return normalized.as_posix()
    try:
        return normalized.relative_to(REPOSITORY_ROOT).as_posix()
    except ValueError:
        return normalized.as_posix()


def path_hash(instance: Path, graph: Path | None) -> str:
    value = recorded_path(instance).encode("utf-8") + b"\0"
    if graph is not None:
        value += recorded_path(graph).encode("utf-8")
    result = 14_695_981_039_346_656_037
    for byte in value:
        result ^= byte
        result = (result * 1_099_511_628_211) & 0xFFFFFFFFFFFFFFFF
    return f"{result:016x}"


def result_key(case: Case) -> str:
    return (
        f"{sanitized_filename(case.instance.stem)}__"
        f"{path_hash(case.instance, case.graph)}"
    )


def output_key(problem: str, case: Case) -> str:
    return f"{problem}__{result_key(case)}.sol"


def solution_reference(problem: str, case: Case) -> str:
    return (Path("solutions") / output_key(problem, case)).as_posix()


def read_completed_results(output_dir: Path, problem: str) -> set[tuple[str, str]]:
    result_path = output_dir / f"{PROBLEM_NAMES[problem]}_Results.csv"
    if not result_path.exists() or result_path.stat().st_size == 0:
        return set()
    if not result_path.is_file():
        raise ValueError(f"result CSV is not a file: {result_path}")

    with result_path.open(encoding="utf-8", newline="") as input_file:
        reader = csv.DictReader(input_file)
        if tuple(reader.fieldnames or ()) != RESULT_COLUMNS:
            raise ValueError(
                f"existing result CSV has an incompatible header: {result_path}"
            )
        completed: set[tuple[str, str]] = set()
        for line_number, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(
                    f"malformed result CSV row {line_number}: {result_path}"
                )
            completed.add((row["instance_key"], row["solution_file"]))
        return completed


def case_is_complete(
    output_dir: Path,
    problem: str,
    case: Case,
    completed_results: set[tuple[str, str]],
) -> bool:
    reference = solution_reference(problem, case)
    solution = output_dir / reference
    return (
        solution.is_file()
        and solution.stat().st_size > 0
        and (result_key(case), reference) in completed_results
    )


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
        output_dir = args.output_dir or REPOSITORY_ROOT / "results" / args.problem
        if not output_dir.is_absolute():
            output_dir = (Path.cwd() / output_dir).resolve()
        completed_results = read_completed_results(output_dir, args.problem)
        failures = 0
        for index, case in enumerate(cases, start=1):
            if case_is_complete(
                output_dir, args.problem, case, completed_results
            ):
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
