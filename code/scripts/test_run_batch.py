#!/usr/bin/env python3

from __future__ import annotations

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


if __name__ == "__main__":
    unittest.main()
