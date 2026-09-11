# CRISPRDecode collaborator test release

This test release covers paired-end, paired-guide decoding without UMI/iBAR processing. It includes the duplicate-TSV-header validation fix and explicit Docker Hub images with the process tools required by Nextflow. It is not a full real-data or production validation of the pipeline.

## Fixed version

Repository: <https://github.com/Yifan1835/crisprseq>  
Test branch: `test/crisprdecode-20260910`  
Release tag (created after CI passes): `crisprdecode-test-20260910`

Use the complete commit SHA shown by the successful **CRISPRDecode collaborator tests** run. Artifacts include `commit.txt`; feedback must include that SHA. Branches may receive fixes, so the branch name alone does not identify the tested version.

## Test on GitHub

1. Open the repository's **Actions** tab and select **CRISPRDecode collaborator tests**. Pushing a `test/crisprdecode-*` branch starts it automatically.
2. Open the run matching the supplied commit. Check all three jobs: **Python validation and assignment**, **Nextflow decoding**, and **Nextflow screening-routing**.
3. Download the three result artifacts. Inspect their logs and `commit.txt`, even if the jobs are green.
4. Return the feedback below. A maintainer with the required repository permissions can rerun a failed run. `Run workflow` is only available when this workflow also exists on the default branch; publishing it on this test branch does not add that button to the default branch.

This workflow uses GitHub-hosted Ubuntu 24.04 runners. It does not require the upstream custom runners, repository secrets, or collaborator data. Nextflow plugins and pinned process containers still require network access. A dependency/download failure is an environment failure until its cause is established.

## Expected outcomes

| Suite             | Expected result                                                                                                      | Scope                                                                                 |
| ----------------- | -------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------- |
| Python            | 8 tests pass, including rejection of all four duplicated required headers and acceptance of a unique metadata column | Python scripts and input validation                                                   |
| Decoding          | 2 nf-tests pass (one stub, one real synthetic decode)                                                                | Three decoding processes; real fixture output below                                   |
| Screening routing | 2 nf-tests pass                                                                                                      | Stub wiring and software-version snapshots for CRISPRDecode and default MAGeCK routes |

The real synthetic fixture has 5 read pairs: 2 unique, 1 ambiguous, 1 unassigned, 1 extraction failure; extracted reads = 4. Counts are 1 for each of `construct_a` and `construct_b`, and 0 for `construct_dup_1`, `construct_dup_2`, and `construct_zero`. Check:

```text
unique_reads + ambiguous_reads + unassigned_reads + extraction_failed_reads = total_reads
2 + 1 + 1 + 1 = 5
```

The count matrix starts with `sgRNA`, `Gene`, then `sample`. Ambiguous signatures must not be arbitrarily assigned. Screening tests use `-stub`: they do not validate real MAGeCK statistics or complete biological data processing. A green run is evidence for these test cases only.

## Run the same tests locally

Use Python 3.12.8 for the unit tests. For the Nextflow tests, install Java 17, Nextflow 25.04.0, nf-test 0.9.3, and Docker. These match the dedicated workflow; inspect its pinned actions for the automated setup.

```bash
git clone https://github.com/Yifan1835/crisprseq.git
cd crisprseq
git checkout --detach crisprdecode-test-20260910
git rev-parse HEAD
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests/crisprdecode -v
NXF_VER=25.04.0 nf-test test subworkflows/local/crisprdecode_paired_guide/tests/main.nf.test --profile=+docker --ci
NXF_VER=25.04.0 nf-test test tests/main_screening_crisprdecode.nf.test --profile=+docker --ci
```

If the release tag is not yet available, use the exact published test commit instead. The fixture provenance and regeneration command are in [README.md](README.md). Do not regenerate fixtures or accept new snapshots merely to make a failing test pass.

## Feedback template

```text
Commit SHA:
Actions run URL (or local OS/Python/Java/Nextflow/nf-test/Docker versions):
Python: PASS / FAIL / NOT RUN
Decoding: PASS / FAIL / NOT RUN
Screening routing: PASS / FAIL / NOT RUN
Expected versus observed counts/QC:
First failing command and error:
Log/artifact location:
Setup minutes / test minutes / review minutes:
Suggested correction:
```

Post feedback in the agreed repository issue/PR, or send the completed template to the project owner. Include logs for failures; never include credentials. This release has not changed the default upstream workflow or expanded the paired-guide feature scope.

## Container and version checks

The three decoding modules use `docker.io/library/python:3.11.4-bookworm@sha256:d7df302a1bcf4db50650da79c174f5d8d973fa4753e0275696644c5bdb477c00`. The explicit registry avoids the pipeline's default `quay.io` prefix. The full image provides `ps` for Nextflow task metrics; the previous slim image did not.

Screening snapshots compare parsed software-version maps: Python 3.11.4 for all three decoding modules, FastQC 0.12.1, and MAGeCK 0.5.9.5 for the default route. Workflow name/version metadata is checked for presence separately, so runtime or pipeline revision metadata does not invalidate software-version comparisons. These expected software values are specified from the pinned modules, not accepted blindly from a failing run.
