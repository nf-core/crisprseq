#!/usr/bin/env python3

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
VALIDATE = ROOT / "bin" / "crisprdecode_validate_library.py"
ASSIGN = ROOT / "bin" / "crisprdecode_assign.py"
AGGREGATE = ROOT / "bin" / "crisprdecode_aggregate.py"


def write_fastq(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w") as handle:
        for name, sequence in records:
            handle.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class CrisprdecodePairedGuideTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.work = Path(self.tempdir.name)
        self.library = self.work / "library.tsv"
        self.validated = self.work / "validated.tsv"
        self.r1 = self.work / "sample_R1.fastq"
        self.r2 = self.work / "sample_R2.fastq"
        self.counts = self.work / "sample.crisprdecode.counts.tsv"
        self.summary = self.work / "sample.crisprdecode.assignment_summary.tsv"

        self.library.write_text(
            "construct_id\ttarget_id\tspacer_r1\tspacer_r2\n"
            "construct_a\tGENE_A\tAAAAA\tCCCCC\n"
            "construct_b\tGENE_B\tGGGGG\tTTTTT\n"
            "construct_dup_1\tGENE_DUP\tACGTA\tTGCAT\n"
            "construct_dup_2\tGENE_DUP\tACGTA\tTGCAT\n"
            "construct_zero\tGENE_ZERO\tCATCA\tGTACG\n"
        )

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def run_script(
        self, script: Path, *args: object, check: bool = True
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(script), *map(str, args)],
            text=True,
            capture_output=True,
            check=check,
        )

    def validate_library(self) -> None:
        self.run_script(VALIDATE, "--library", self.library, "--output", self.validated)

    def test_exact_assignment_ambiguity_and_accounting(self) -> None:
        self.validate_library()
        write_fastq(
            self.r1,
            [
                ("read_a_1/1", "GGAAAAA"),
                ("read_a_2/1", "GGAAAAA"),
                ("read_b/1", "GGGGGGG"),
                ("read_ambiguous/1", "GGACGTA"),
                ("read_unassigned/1", "GGTAAAA"),
                ("read_failed/1", "TTAAAAA"),
            ],
        )
        write_fastq(
            self.r2,
            [
                ("read_a_1/2", "GGCCCCC"),
                ("read_a_2/2", "GGCCCCC"),
                ("read_b/2", "GGTTTTT"),
                ("read_ambiguous/2", "GGTGCAT"),
                ("read_unassigned/2", "GGAAAAA"),
                ("read_failed/2", "GGCCCCC"),
            ],
        )
        self.run_script(
            ASSIGN,
            "--r1",
            self.r1,
            "--r2",
            self.r2,
            "--library",
            self.validated,
            "--sample",
            "sample",
            "--r1-anchor",
            "GG",
            "--r2-anchor",
            "GG",
            "--counts",
            self.counts,
            "--summary",
            self.summary,
        )

        counts = {
            row["construct_id"]: int(row["count"]) for row in read_tsv(self.counts)
        }
        self.assertEqual(
            counts,
            {
                "construct_a": 2,
                "construct_b": 1,
                "construct_dup_1": 0,
                "construct_dup_2": 0,
                "construct_zero": 0,
            },
        )
        summary = read_tsv(self.summary)[0]
        self.assertEqual(
            {key: int(value) for key, value in summary.items() if key != "sample"},
            {
                "total_reads": 6,
                "extracted_reads": 5,
                "unique_reads": 3,
                "ambiguous_reads": 1,
                "unassigned_reads": 1,
                "extraction_failed_reads": 1,
            },
        )
        validated = {row["construct_id"]: row for row in read_tsv(self.validated)}
        self.assertEqual(
            validated["construct_dup_1"]["assignability_class"], "ambiguous"
        )
        self.assertEqual(validated["construct_dup_1"]["duplicate_signature_size"], "2")

        matrix = self.work / "count_table.count.txt"
        combined_summary = self.work / "assignment_summary.tsv"
        recovery = self.work / "library_recovery.tsv"
        self.run_script(
            AGGREGATE,
            "--library",
            self.validated,
            "--sample-counts",
            self.counts,
            "--sample-summaries",
            self.summary,
            "--count-matrix",
            matrix,
            "--assignment-summary",
            combined_summary,
            "--library-recovery",
            recovery,
        )
        self.assertEqual(matrix.read_text().splitlines()[0], "sgRNA\tGene\tsample")
        recovery_by_id = {row["construct_id"]: row for row in read_tsv(recovery)}
        self.assertEqual(recovery_by_id["construct_a"]["total_count"], "2")
        self.assertEqual(recovery_by_id["construct_zero"]["zero_count"], "true")

    def test_desynchronised_fastq_fails(self) -> None:
        self.validate_library()
        write_fastq(self.r1, [("read_1/1", "AAAAA")])
        write_fastq(self.r2, [("different_read/2", "CCCCC")])
        result = self.run_script(
            ASSIGN,
            "--r1",
            self.r1,
            "--r2",
            self.r2,
            "--library",
            self.validated,
            "--sample",
            "sample",
            "--counts",
            self.counts,
            "--summary",
            self.summary,
            check=False,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("desynchronised", result.stderr)

    def test_different_fastq_record_counts_fail(self) -> None:
        self.validate_library()
        write_fastq(self.r1, [("read_1/1", "AAAAA"), ("read_2/1", "GGGGG")])
        write_fastq(self.r2, [("read_1/2", "CCCCC")])
        result = self.run_script(
            ASSIGN,
            "--r1",
            self.r1,
            "--r2",
            self.r2,
            "--library",
            self.validated,
            "--sample",
            "sample",
            "--counts",
            self.counts,
            "--summary",
            self.summary,
            check=False,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("different numbers of records", result.stderr)

    def test_negative_offset_fails(self) -> None:
        self.validate_library()
        write_fastq(self.r1, [("read_1/1", "AAAAA")])
        write_fastq(self.r2, [("read_1/2", "CCCCC")])
        result = self.run_script(
            ASSIGN,
            "--r1",
            self.r1,
            "--r2",
            self.r2,
            "--library",
            self.validated,
            "--sample",
            "sample",
            "--r1-offset",
            -1,
            "--counts",
            self.counts,
            "--summary",
            self.summary,
            check=False,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("R1 offset must be non-negative", result.stderr)

    def test_duplicate_construct_id_fails_validation(self) -> None:
        self.library.write_text(
            "construct_id\ttarget_id\tspacer_r1\tspacer_r2\n"
            "duplicate\tGENE_A\tAAAAA\tCCCCC\n"
            "duplicate\tGENE_B\tGGGGG\tTTTTT\n"
        )
        result = self.run_script(
            VALIDATE,
            "--library",
            self.library,
            "--output",
            self.validated,
            check=False,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("construct_id values must be unique", result.stderr)

    def test_nonuniform_spacer_lengths_fail_validation(self) -> None:
        self.library.write_text(
            "construct_id\ttarget_id\tspacer_r1\tspacer_r2\n"
            "construct_a\tGENE_A\tAAAAA\tCCCCC\n"
            "construct_b\tGENE_B\tGGGG\tTTTTT\n"
        )
        result = self.run_script(
            VALIDATE,
            "--library",
            self.library,
            "--output",
            self.validated,
            check=False,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("spacer_r1 must have one uniform length", result.stderr)


if __name__ == "__main__":
    unittest.main()
