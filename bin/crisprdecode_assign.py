#!/usr/bin/env python3
"""Assign synchronized paired-end reads to exact paired-guide constructs."""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter, defaultdict
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--library", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--r1-anchor", default="")
    parser.add_argument("--r2-anchor", default="")
    parser.add_argument("--r1-offset", default=0, type=int)
    parser.add_argument("--r2-offset", default=0, type=int)
    parser.add_argument("--reverse-complement-r1", action="store_true")
    parser.add_argument("--reverse-complement-r2", action="store_true")
    parser.add_argument("--counts", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    return parser.parse_args()


@contextmanager
def open_text(path: Path) -> Iterator[TextIO]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as handle:
            yield handle
    else:
        with path.open() as handle:
            yield handle


def fastq_records(handle: TextIO, source: Path) -> Iterator[tuple[str, str]]:
    record_number = 0
    while True:
        header = handle.readline()
        if not header:
            return
        sequence = handle.readline()
        separator = handle.readline()
        quality = handle.readline()
        record_number += 1
        if not sequence or not separator or not quality:
            raise ValueError(f"truncated FASTQ record {record_number} in {source}")
        if not header.startswith("@") or not separator.startswith("+"):
            raise ValueError(f"malformed FASTQ record {record_number} in {source}")
        sequence = sequence.strip().upper()
        if len(sequence) != len(quality.strip()):
            raise ValueError(
                f"sequence/quality length mismatch in record {record_number} of {source}"
            )
        yield header[1:].split()[0].removesuffix("/1").removesuffix("/2"), sequence


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def extract(
    sequence: str,
    *,
    anchor: str,
    offset: int,
    length: int,
    use_reverse_complement: bool,
) -> str | None:
    if use_reverse_complement:
        sequence = reverse_complement(sequence)
    start = offset
    if anchor:
        anchor = anchor.upper()
        anchor_start = sequence.find(anchor)
        if anchor_start < 0:
            return None
        start = anchor_start + len(anchor) + offset
    if start < 0 or start + length > len(sequence):
        return None
    return sequence[start : start + length]


def load_library(
    path: Path,
) -> tuple[list[dict[str, str]], dict[tuple[str, str], list[str]]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError("validated construct library is empty")
    signatures: dict[tuple[str, str], list[str]] = defaultdict(list)
    for row in rows:
        signatures[(row["spacer_r1"], row["spacer_r2"])].append(row["construct_id"])
    return rows, signatures


def write_counts(
    rows: list[dict[str, str]], counts: Counter[str], sample: str, path: Path
) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["construct_id", "target_id", "sample", "count"])
        for row in rows:
            writer.writerow(
                [
                    row["construct_id"],
                    row["target_id"],
                    sample,
                    counts[row["construct_id"]],
                ]
            )


def write_summary(sample: str, metrics: Counter[str], path: Path) -> None:
    extracted = metrics["unique"] + metrics["ambiguous"] + metrics["unassigned"]
    if extracted + metrics["extraction_failed"] != metrics["total"]:
        raise AssertionError("read-accounting invariant failed")
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "sample",
                "total_reads",
                "extracted_reads",
                "unique_reads",
                "ambiguous_reads",
                "unassigned_reads",
                "extraction_failed_reads",
            ]
        )
        writer.writerow(
            [
                sample,
                metrics["total"],
                extracted,
                metrics["unique"],
                metrics["ambiguous"],
                metrics["unassigned"],
                metrics["extraction_failed"],
            ]
        )


def main() -> None:
    args = parse_args()
    try:
        for read, anchor, offset in (
            ("R1", args.r1_anchor, args.r1_offset),
            ("R2", args.r2_anchor, args.r2_offset),
        ):
            if offset < 0:
                raise ValueError(f"{read} offset must be non-negative")
            invalid_anchor_bases = sorted(set(anchor.upper()) - set("ACGT"))
            if invalid_anchor_bases:
                raise ValueError(
                    f"{read} anchor contains unsupported base(s): "
                    f"{''.join(invalid_anchor_bases)}; only A/C/G/T are supported"
                )

        rows, signatures = load_library(args.library)
        r1_length = len(rows[0]["spacer_r1"])
        r2_length = len(rows[0]["spacer_r2"])
        counts: Counter[str] = Counter()
        metrics: Counter[str] = Counter()

        with open_text(args.r1) as r1_handle, open_text(args.r2) as r2_handle:
            r1_records = fastq_records(r1_handle, args.r1)
            r2_records = fastq_records(r2_handle, args.r2)
            while True:
                r1_record = next(r1_records, None)
                r2_record = next(r2_records, None)
                if r1_record is None and r2_record is None:
                    break
                if r1_record is None or r2_record is None:
                    raise ValueError(
                        "paired FASTQ files contain different numbers of records"
                    )
                r1_id, r1_sequence = r1_record
                r2_id, r2_sequence = r2_record
                if r1_id != r2_id:
                    raise ValueError(
                        f"paired FASTQ records are desynchronised: R1 '{r1_id}' != R2 '{r2_id}'"
                    )

                metrics["total"] += 1
                spacer_r1 = extract(
                    r1_sequence,
                    anchor=args.r1_anchor,
                    offset=args.r1_offset,
                    length=r1_length,
                    use_reverse_complement=args.reverse_complement_r1,
                )
                spacer_r2 = extract(
                    r2_sequence,
                    anchor=args.r2_anchor,
                    offset=args.r2_offset,
                    length=r2_length,
                    use_reverse_complement=args.reverse_complement_r2,
                )
                if spacer_r1 is None or spacer_r2 is None:
                    metrics["extraction_failed"] += 1
                    continue

                matched_constructs = signatures.get((spacer_r1, spacer_r2), [])
                if not matched_constructs:
                    metrics["unassigned"] += 1
                elif len(matched_constructs) > 1:
                    metrics["ambiguous"] += 1
                else:
                    metrics["unique"] += 1
                    counts[matched_constructs[0]] += 1

        write_counts(rows, counts, args.sample, args.counts)
        write_summary(args.sample, metrics, args.summary)
    except ValueError as error:
        raise SystemExit(f"ERROR: {error}") from error


if __name__ == "__main__":
    main()
