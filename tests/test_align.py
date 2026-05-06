"""Tests for varus.align.

The full hisat2 | samtools sort pipe is exercised on the cluster, not here. We
unit-test the log parser and verify that the wrapper checks for the binaries.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from varus import align


def test_parse_hisat2_log_paired(tmp_path: Path):
    # Real-shape HISAT2 paired-end summary; the legacy regex picks up the
    # "aligned concordantly exactly 1 time" line.
    log = tmp_path / "Log.final.out"
    log.write_text(
        "50000 reads; of these:\n"
        "  50000 (100.00%) were paired; of these:\n"
        "    21662 (43.32%) aligned concordantly exactly 1 time\n"
        "    100 (0.20%) aligned concordantly >1 times\n"
        "    1234 (2.47%) aligned discordantly 1 time\n"
    )
    stats = align.parse_hisat2_log(log, batch_size=50000)
    assert stats["num_uniq"] == 21662
    assert abs(stats["uniq_pct"] - 43.324) < 0.01


def test_parse_hisat2_log_single_end(tmp_path: Path):
    log = tmp_path / "Log.final.out"
    log.write_text(
        "50000 reads; of these:\n"
        "  29826 (59.65%) aligned exactly 1 time\n"
        "  900 (1.80%) aligned >1 times\n"
    )
    stats = align.parse_hisat2_log(log, batch_size=50000)
    assert stats["num_uniq"] == 29826


def test_parse_hisat2_log_no_match(tmp_path: Path):
    log = tmp_path / "Log.final.out"
    log.write_text("nothing useful here\n")
    stats = align.parse_hisat2_log(log, batch_size=50000)
    assert stats["num_uniq"] == 0
    assert stats["uniq_pct"] == 0.0


def test_parse_hisat2_log_missing_file(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        align.parse_hisat2_log(tmp_path / "nope.out", batch_size=10)


def test_align_batch_requires_hisat2(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(align.shutil, "which", lambda t: None)
    with pytest.raises(RuntimeError, match="hisat2 not found"):
        align.align_batch_hisat2(
            r1=tmp_path / "r1.fa",
            r2=None,
            index_prefix=tmp_path / "idx",
            batch_dir=tmp_path / "b",
        )


def test_align_batch_requires_samtools(tmp_path: Path, monkeypatch):
    # hisat2 found, samtools missing
    monkeypatch.setattr(align.shutil, "which",
                        lambda t: "/fake/hisat2" if t == "hisat2" else None)
    with pytest.raises(RuntimeError, match="samtools not found"):
        align.align_batch_hisat2(
            r1=tmp_path / "r1.fa",
            r2=None,
            index_prefix=tmp_path / "idx",
            batch_dir=tmp_path / "b",
        )
