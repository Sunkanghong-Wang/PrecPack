#!/usr/bin/env python3
"""Configure, build, and test PrecPack in Release mode."""

from __future__ import annotations

import argparse
import subprocess
import sys

from script_utils import (
    BUILD_DIRECTORY,
    CODE_ROOT,
    REPOSITORY_ROOT,
    controlled_environment,
    find_cmake,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gurobi",
        choices=("auto", "on", "off"),
        default="auto",
        help=(
            "optional Gurobi support: auto detects it, on requires it, "
            "and off builds the commercial-solver-free exact solver"
        ),
    )
    return parser.parse_args()


def run(command: list[str], environment: dict[str, str]) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=REPOSITORY_ROOT, env=environment, check=True)


def main() -> int:
    try:
        args = parse_args()
        cmake = find_cmake()
        environment = controlled_environment()
        configure = [
            cmake,
            "-S",
            str(CODE_ROOT),
            "-B",
            str(BUILD_DIRECTORY),
            "-DCMAKE_BUILD_TYPE=Release",
            "-DBUILD_TESTING=ON",
            "-DPRECPACK_NATIVE_OPTIMIZATION=ON",
            "-DPRECPACK_ENABLE_IPO=ON",
            "-DPRECPACK_ENABLE_SANITIZERS=OFF",
            f"-DPRECPACK_GUROBI={args.gurobi.upper()}",
        ]
        gurobi_home = environment.get("GUROBI_HOME", "")
        configure.append(f"-DGUROBI_ROOT={gurobi_home}")
        run(configure, environment)
        run(
            [
                cmake,
                "--build",
                str(BUILD_DIRECTORY),
                "--config",
                "Release",
                "--target",
                "precpack",
                "--parallel",
            ],
            environment,
        )
        run(
            [
                cmake,
                "--build",
                str(BUILD_DIRECTORY),
                "--config",
                "Release",
                "--target",
                "check",
            ],
            environment,
        )
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
