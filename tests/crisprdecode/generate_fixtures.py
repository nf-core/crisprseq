#!/usr/bin/env python3
"""Regenerate the deterministic CRISPRDecode paired-guide truth fixture."""

import gzip
from pathlib import Path


FIXTURE_DIR = Path(__file__).parent / "fixtures"


def fastq(records: list[tuple[str, str]]) -> str:
    return "".join(
        f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n" for name, sequence in records
    )


def main() -> None:
    FIXTURE_DIR.mkdir(exist_ok=True)
    (FIXTURE_DIR / "library.tsv").write_text(
        "construct_id\ttarget_id\tspacer_r1\tspacer_r2\n"
        "construct_a\tGENE_A\tAAAAA\tCCCCC\n"
        "construct_b\tGENE_B\tGGGGG\tTTTTT\n"
        "construct_dup_1\tGENE_DUP\tACGTA\tTGCAT\n"
        "construct_dup_2\tGENE_DUP\tACGTA\tTGCAT\n"
        "construct_zero\tGENE_ZERO\tCATCA\tGTACG\n"
    )
    r1_text = fastq(
        [
            ("read_a/1", "AAAAA"),
            ("read_b/1", "GGGGG"),
            ("read_ambiguous/1", "ACGTA"),
            ("read_unassigned/1", "TAAAA"),
            ("read_failed/1", "AAA"),
        ]
    )
    r2_text = fastq(
        [
            ("read_a/2", "CCCCC"),
            ("read_b/2", "TTTTT"),
            ("read_ambiguous/2", "TGCAT"),
            ("read_unassigned/2", "AAAAA"),
            ("read_failed/2", "CCCCC"),
        ]
    )
    (FIXTURE_DIR / "sample_R1.fastq").write_text(r1_text)
    (FIXTURE_DIR / "sample_R2.fastq").write_text(r2_text)
    (FIXTURE_DIR / "sample_R1.fastq.gz").write_bytes(
        gzip.compress(r1_text.encode(), mtime=0)
    )
    (FIXTURE_DIR / "sample_R2.fastq.gz").write_bytes(
        gzip.compress(r2_text.encode(), mtime=0)
    )


if __name__ == "__main__":
    main()
