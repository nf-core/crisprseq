#!/usr/bin/env python3
"""Validate and normalize a paired-guide construct library."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


REQUIRED_COLUMNS = ("construct_id", "target_id", "spacer_r1", "spacer_r2")
DNA_ALPHABET = frozenset("ACGT")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--library", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def load_library(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("construct library is empty")
        missing = [
            column for column in REQUIRED_COLUMNS if column not in reader.fieldnames
        ]
        if missing:
            raise ValueError(
                "construct library is missing required column(s): " + ", ".join(missing)
            )
        rows = []
        for line_number, raw_row in enumerate(reader, start=2):
            row = {
                column: (raw_row.get(column) or "").strip()
                for column in REQUIRED_COLUMNS
            }
            empty = [column for column, value in row.items() if not value]
            if empty:
                raise ValueError(
                    f"construct library line {line_number} has empty field(s): {', '.join(empty)}"
                )
            for column in ("spacer_r1", "spacer_r2"):
                row[column] = row[column].upper()
                invalid = sorted(set(row[column]) - DNA_ALPHABET)
                if invalid:
                    raise ValueError(
                        f"construct library line {line_number} column {column} contains "
                        f"unsupported base(s): {''.join(invalid)}; only A/C/G/T are supported"
                    )
            rows.append(row)

    if not rows:
        raise ValueError("construct library contains a header but no constructs")

    ids = [row["construct_id"] for row in rows]
    duplicate_ids = sorted(key for key, count in Counter(ids).items() if count > 1)
    if duplicate_ids:
        raise ValueError(
            "construct_id values must be unique: " + ", ".join(duplicate_ids)
        )

    for column in ("spacer_r1", "spacer_r2"):
        lengths = sorted({len(row[column]) for row in rows})
        if len(lengths) != 1:
            raise ValueError(
                f"{column} must have one uniform length for fixed-window extraction; found {lengths}"
            )

    return sorted(rows, key=lambda row: row["construct_id"])


def write_library(rows: list[dict[str, str]], path: Path) -> None:
    signature_sizes = Counter((row["spacer_r1"], row["spacer_r2"]) for row in rows)
    fieldnames = [
        *REQUIRED_COLUMNS,
        "assignability_class",
        "duplicate_signature_size",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            duplicate_size = signature_sizes[(row["spacer_r1"], row["spacer_r2"])]
            writer.writerow(
                {
                    **row,
                    "assignability_class": "unique"
                    if duplicate_size == 1
                    else "ambiguous",
                    "duplicate_signature_size": duplicate_size,
                }
            )


def main() -> None:
    args = parse_args()
    try:
        rows = load_library(args.library)
        write_library(rows, args.output)
    except ValueError as error:
        raise SystemExit(f"ERROR: {error}") from error


if __name__ == "__main__":
    main()
