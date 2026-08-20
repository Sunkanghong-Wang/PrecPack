#!/usr/bin/env python3

from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import run_batch


class BatchPairingTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory(
            prefix="precpack-batch-test-"
        )
        self.root = Path(self.temporary_directory.name)

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def make_file(self, path: Path) -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
        return path.resolve()

    def write_result_csv(
        self, output_dir: Path, problem: str, case: run_batch.Case
    ) -> Path:
        result_path = output_dir / f"{run_batch.PROBLEM_NAMES[problem]}_Results.csv"
        result_path.parent.mkdir(parents=True, exist_ok=True)
        row = {column: "" for column in run_batch.RESULT_COLUMNS}
        row["instance_key"] = run_batch.result_key(case)
        row["solution_file"] = run_batch.solution_reference(problem, case)
        with result_path.open("w", encoding="utf-8", newline="") as output_file:
            writer = csv.DictWriter(output_file, fieldnames=run_batch.RESULT_COLUMNS)
            writer.writeheader()
            writer.writerow(row)
        return result_path

    def test_single_instance_and_graph_files(self) -> None:
        instance = self.make_file(self.root / "items" / "case.txt")
        graph = self.make_file(self.root / "graphs" / "case.graph")

        self.assertEqual(
            run_batch.collect_cases("bpp-gp", instance, graph),
            [run_batch.Case(instance, graph)],
        )

    def test_size_directories(self) -> None:
        item_directory = self.root / "items" / "n_0020"
        graph_directory = self.root / "graphs" / "n_0020"
        instance = self.make_file(item_directory / "case.txt")
        graph = self.make_file(graph_directory / "case.graph")

        self.assertEqual(
            run_batch.collect_cases("bpp-gp", item_directory, graph_directory),
            [run_batch.Case(instance, graph)],
        )

    def test_single_instance_with_graph_collection(self) -> None:
        instance = self.make_file(self.root / "items" / "case.txt")
        graph_01 = self.make_file(
            self.root / "graphs" / "separation-01" / "case.graph"
        )
        graph_03 = self.make_file(
            self.root / "graphs" / "separation-03" / "case.graph"
        )
        self.make_file(self.root / "graphs" / "separation-03" / "other.graph")

        self.assertEqual(
            run_batch.collect_cases("bpp-gp", instance, self.root / "graphs"),
            [
                run_batch.Case(instance, graph_01),
                run_batch.Case(instance, graph_03),
            ],
        )

    def test_instance_collection_with_single_graph(self) -> None:
        item_root = self.root / "items"
        instance = self.make_file(item_root / "n_0020" / "case.txt")
        self.make_file(item_root / "n_0020" / "other.txt")
        graph = self.make_file(self.root / "graphs" / "n_0020" / "case.graph")

        self.assertEqual(
            run_batch.collect_cases("bpp-gp", item_root, graph),
            [run_batch.Case(instance, graph)],
        )

    def test_collection_roots(self) -> None:
        item_root = self.root / "items"
        graph_root = self.root / "separation-01"
        instance = self.make_file(item_root / "n_0020" / "case.txt")
        graph = self.make_file(graph_root / "n_0020" / "case.graph")

        self.assertEqual(
            run_batch.collect_cases("bpp-gp", item_root, graph_root),
            [run_batch.Case(instance, graph)],
        )

    def test_default_collection(self) -> None:
        data_root = self.root / "instances"
        graph_root = self.root / "bpp-gp-graphs"
        instance = self.make_file(data_root / "otto" / "n_0020" / "case.txt")
        graph_01 = self.make_file(
            graph_root / "separation-01" / "n_0020" / "case.graph"
        )
        graph_03 = self.make_file(
            graph_root / "separation-03" / "n_0020" / "case.graph"
        )

        with mock.patch.object(run_batch, "DATA", data_root), mock.patch.object(
            run_batch, "BPP_GP_GRAPHS", graph_root
        ):
            cases = run_batch.collect_cases("bpp-gp", None, None)

        self.assertEqual(
            cases,
            [
                run_batch.Case(instance, graph_01),
                run_batch.Case(instance, graph_03),
            ],
        )

    def test_mismatched_single_files_are_rejected(self) -> None:
        instance = self.make_file(self.root / "items" / "item.txt")
        graph = self.make_file(self.root / "graphs" / "other.graph")

        with self.assertRaisesRegex(ValueError, "no graph with stem"):
            run_batch.collect_cases("bpp-gp", instance, graph)

    def test_same_stem_in_different_directories_has_unique_output(self) -> None:
        first = run_batch.Case(self.make_file(self.root / "first" / "case.txt"))
        second = run_batch.Case(self.make_file(self.root / "second" / "case.txt"))

        self.assertNotEqual(run_batch.result_key(first), run_batch.result_key(second))
        self.assertNotEqual(
            run_batch.solution_reference("bpp-p", first),
            run_batch.solution_reference("bpp-p", second),
        )

    def test_same_instance_with_different_graphs_has_unique_output(self) -> None:
        instance = self.make_file(self.root / "items" / "case.txt")
        first = run_batch.Case(
            instance, self.make_file(self.root / "first" / "case.graph")
        )
        second = run_batch.Case(
            instance, self.make_file(self.root / "second" / "case.graph")
        )

        self.assertNotEqual(run_batch.result_key(first), run_batch.result_key(second))

    def test_path_hash_is_stable(self) -> None:
        case = run_batch.Case(Path("external-a/case.txt"))
        self.assertEqual(run_batch.result_key(case), "case__3a95170ceb4e4d41")

    def test_resume_requires_matching_csv_row_and_solution(self) -> None:
        output_dir = self.root / "results"
        instance = self.make_file(self.root / "items" / "case.txt")
        case = run_batch.Case(instance)
        solution = output_dir / run_batch.solution_reference("bpp-p", case)
        solution.parent.mkdir(parents=True)
        solution.write_text("Bin 1: 1\n", encoding="utf-8")

        completed = run_batch.read_completed_results(output_dir, "bpp-p")
        self.assertFalse(
            run_batch.case_is_complete(output_dir, "bpp-p", case, completed)
        )

        self.write_result_csv(output_dir, "bpp-p", case)
        completed = run_batch.read_completed_results(output_dir, "bpp-p")
        self.assertTrue(
            run_batch.case_is_complete(output_dir, "bpp-p", case, completed)
        )

        solution.write_text("", encoding="utf-8")
        self.assertFalse(
            run_batch.case_is_complete(output_dir, "bpp-p", case, completed)
        )

    def test_resume_rejects_incompatible_csv(self) -> None:
        output_dir = self.root / "results"
        result_path = output_dir / "BPP-P_Results.csv"
        result_path.parent.mkdir(parents=True)
        result_path.write_text("old,header\n", encoding="utf-8")

        with self.assertRaisesRegex(ValueError, "incompatible header"):
            run_batch.read_completed_results(output_dir, "bpp-p")


if __name__ == "__main__":
    unittest.main()
