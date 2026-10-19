#!/usr/bin/env python3
"""Aggregate per-sample paired-guide counts into crisprseq outputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--library", required=True, type=Path)
    parser.add_argument("--sample-counts", required=True, nargs="+", type=Path)
    parser.add_argument("--sample-summaries", required=True, nargs="+", type=Path)
    parser.add_argument("--count-matrix", required=True, type=Path)
    parser.add_argument("--assignment-summary", required=True, type=Path)
    parser.add_argument("--library-recovery", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    args = parse_args()
    library = read_tsv(args.library)
    if not library:
        raise SystemExit("ERROR: validated construct library is empty")

    counts_by_sample: dict[str, dict[str, int]] = {}
    target_by_construct = {row["construct_id"]: row["target_id"] for row in library}
    for path in args.sample_counts:
        rows = read_tsv(path)
        if not rows:
            raise SystemExit(f"ERROR: sample count file is empty: {path}")
        samples = {row["sample"] for row in rows}
        if len(samples) != 1:
            raise SystemExit(
                f"ERROR: expected one sample in {path}, found {sorted(samples)}"
            )
        sample = samples.pop()
        if sample in counts_by_sample:
            raise SystemExit(f"ERROR: duplicate sample count file for '{sample}'")
        counts_by_sample[sample] = {
            row["construct_id"]: int(row["count"]) for row in rows
        }

    samples = sorted(counts_by_sample)
    with args.count_matrix.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sgRNA", "Gene", *samples])
        for row in library:
            construct_id = row["construct_id"]
            writer.writerow(
                [
                    construct_id,
                    target_by_construct[construct_id],
                    *[
                        counts_by_sample[sample].get(construct_id, 0)
                        for sample in samples
                    ],
                ]
            )

    summaries: list[dict[str, str]] = []
    for path in args.sample_summaries:
        rows = read_tsv(path)
        if len(rows) != 1:
            raise SystemExit(f"ERROR: expected exactly one summary row in {path}")
        summaries.extend(rows)
    summary_fields = [
        "sample",
        "total_reads",
        "extracted_reads",
        "unique_reads",
        "ambiguous_reads",
        "unassigned_reads",
        "extraction_failed_reads",
    ]
    with args.assignment_summary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=summary_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(sorted(summaries, key=lambda row: row["sample"]))

    recovery_fields = [
        "construct_id",
        "target_id",
        "assignability_class",
        "duplicate_signature_size",
        "total_count",
        "zero_count",
    ]
    with args.library_recovery.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=recovery_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in library:
            construct_id = row["construct_id"]
            total_count = sum(
                counts_by_sample[sample].get(construct_id, 0) for sample in samples
            )
            writer.writerow(
                {
                    "construct_id": construct_id,
                    "target_id": row["target_id"],
                    "assignability_class": row["assignability_class"],
                    "duplicate_signature_size": row["duplicate_signature_size"],
                    "total_count": total_count,
                    "zero_count": str(total_count == 0).lower(),
                }
            )


if __name__ == "__main__":
    main()
