"""Tests for varus.download. We mock ``shutil.which`` and ``subprocess.run``
so no SRA download actually happens."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from varus import download


def _stub_run_creates(file_paths: list[Path]):
    """Return a fake subprocess.run that touches the expected output files."""
    def _run(cmd, check):
        for p in file_paths:
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text("")
        return subprocess.CompletedProcess(cmd, 0)
    return _run


def test_batch_dir_layout(tmp_path: Path):
    bdir = download.batch_dir_for(tmp_path, "SRR000001", 0, 49999)
    assert bdir == tmp_path / "batches" / "SRR000001" / "N0X49999"


def test_download_batch_single_end_command(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(download.shutil, "which", lambda _: "/fake/fastq-dump")
    captured = {}

    def fake_run(cmd, check):
        captured["cmd"] = list(cmd)
        bdir = Path(cmd[cmd.index("-O") + 1])
        (bdir / "SRR1.fasta").write_text("")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(download.subprocess, "run", fake_run)

    paths = download.download_batch(
        accession="SRR1", n=0, x=49999, paired=False, outdir=tmp_path,
    )
    assert paths.r2 is None
    assert paths.r1.name == "SRR1.fasta"
    cmd = captured["cmd"]
    assert cmd[:3] == ["fastq-dump", "-N", "0"]
    assert "-X" in cmd and cmd[cmd.index("-X") + 1] == "49999"
    assert "--fasta" in cmd
    assert "--split-files" not in cmd
    assert cmd[-1] == "SRR1"


def test_download_batch_paired_adds_split_files(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(download.shutil, "which", lambda _: "/fake/fastq-dump")

    def fake_run(cmd, check):
        bdir = Path(cmd[cmd.index("-O") + 1])
        (bdir / "SRR2_1.fasta").write_text("")
        (bdir / "SRR2_2.fasta").write_text("")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(download.subprocess, "run", fake_run)

    paths = download.download_batch(
        accession="SRR2", n=100000, x=149999, paired=True, outdir=tmp_path,
    )
    assert paths.r1.name == "SRR2_1.fasta"
    assert paths.r2 is not None and paths.r2.name == "SRR2_2.fasta"


def test_download_batch_paired_three_file_fallback(tmp_path: Path, monkeypatch):
    """Some SRA runs (e.g. SRR097898) yield 3 FASTAs; legacy uses [0] and [-1]."""
    monkeypatch.setattr(download.shutil, "which", lambda _: "/fake/fastq-dump")

    def fake_run(cmd, check):
        bdir = Path(cmd[cmd.index("-O") + 1])
        # Note: no _2.fasta -- the second file has a different name
        (bdir / "SRR3_1.fasta").write_text("")
        (bdir / "SRR3_2.fasta_intermediate").write_text("")  # ignored, not .fasta
        (bdir / "SRR3_3.fasta").write_text("")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(download.subprocess, "run", fake_run)

    paths = download.download_batch(
        accession="SRR3", n=0, x=49999, paired=True, outdir=tmp_path,
    )
    assert paths.r1.name == "SRR3_1.fasta"
    assert paths.r2 is not None and paths.r2.name == "SRR3_3.fasta"


def test_download_batch_retries_then_raises(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(download.shutil, "which", lambda _: "/fake/fastq-dump")
    calls = {"n": 0}

    def always_fail(cmd, check):
        calls["n"] += 1
        raise subprocess.CalledProcessError(returncode=3, cmd=cmd)

    monkeypatch.setattr(download.subprocess, "run", always_fail)

    with pytest.raises(RuntimeError, match="fastq-dump failed"):
        download.download_batch(
            accession="SRRfail", n=0, x=49999, paired=False, outdir=tmp_path,
            retries=2,
        )
    assert calls["n"] == 2


def test_download_batch_missing_tool(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(download.shutil, "which", lambda _: None)
    with pytest.raises(RuntimeError, match="fastq-dump not found"):
        download.download_batch("X", 0, 1, False, tmp_path)


def test_download_full_uses_fasterq_dump(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(download.shutil, "which", lambda _: "/fake/fasterq-dump")
    captured = {}

    def fake_run(cmd, check):
        captured["cmd"] = list(cmd)
        out_idx = cmd.index("-O") + 1
        bdir = Path(cmd[out_idx])
        (bdir / "SRR9.fasta").write_text("")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(download.subprocess, "run", fake_run)

    paths = download.download_full(
        accession="SRR9", paired=False, outdir=tmp_path, threads=8,
    )
    assert paths.r2 is None
    assert "fasterq-dump" in captured["cmd"][0]
    assert "--threads" in captured["cmd"]
    assert captured["cmd"][captured["cmd"].index("--threads") + 1] == "8"
