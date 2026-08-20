from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
CODE_ROOT = SCRIPT_DIRECTORY.parent
REPOSITORY_ROOT = CODE_ROOT.parent
BUILD_DIRECTORY = REPOSITORY_ROOT / "build"


def find_cmake() -> str:
    configured = os.environ.get("CMAKE_BIN")
    if configured:
        return configured
    discovered = shutil.which("cmake")
    if discovered:
        return discovered
    raise RuntimeError("CMake was not found; install it or set CMAKE_BIN")


def find_solver() -> Path:
    name = "precpack.exe" if os.name == "nt" else "precpack"
    for directory in (BUILD_DIRECTORY / "Release", BUILD_DIRECTORY):
        candidate = directory / name
        if candidate.is_file():
            return candidate
    raise RuntimeError("PrecPack is not built; run code/scripts/build.py first")


def controlled_environment() -> dict[str, str]:
    environment = os.environ.copy()
    for variable in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "BLIS_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ):
        environment[variable] = "1"

    gurobi_home = environment.get("GUROBI_HOME")
    if not gurobi_home:
        return environment

    runtime_directory = Path(gurobi_home) / ("bin" if os.name == "nt" else "lib")
    variable = (
        "PATH"
        if os.name == "nt"
        else "DYLD_LIBRARY_PATH" if sys.platform == "darwin" else "LD_LIBRARY_PATH"
    )
    environment[variable] = (
        str(runtime_directory) + os.pathsep + environment.get(variable, "")
    )
    return environment
