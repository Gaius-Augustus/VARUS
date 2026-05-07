"""Tests for varus.index.

We don't actually run hisat2-build; we mock subprocess.run and shutil.which to
verify the command construction and error paths.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from varus import index


def test_build_hisat2_index_missing_genome(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        index.build_hisat2_index(
            genome=tmp_path / "nope.fa",
            outdir=tmp_path / "idx",
        )


def test_build_hisat2_index_missing_tool(tmp_path: Path, monkeypatch):
    genome = tmp_path / "g.fa"
    genome.write_text(">chr1\nACGT\n")
    monkeypatch.setattr(index.shutil, "which", lambda _: None)
    with pytest.raises(RuntimeError, match="hisat2-build not found"):
        index.build_hisat2_index(genome=genome, outdir=tmp_path / "idx")


def test_build_hisat2_index_invokes_correct_command(tmp_path: Path, monkeypatch):
    genome = tmp_path / "g.fa"
    genome.write_text(">chr1\nACGT\n")
    outdir = tmp_path / "idx"

    monkeypatch.setattr(index.shutil, "which", lambda _: "/fake/hisat2-build")
    captured = {}

    def fake_run(cmd, check):
        captured["cmd"] = cmd
        captured["check"] = check
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(index.subprocess, "run", fake_run)

    out = index.build_hisat2_index(
        genome=genome, outdir=outdir, threads=8, prefix="myidx"
    )
    assert out == outdir / "myidx"
    assert captured["check"] is True
    assert captured["cmd"] == [
        "hisat2-build", "-p", "8", str(genome), str(outdir / "myidx"),
    ]
    # outdir was created by the function.
    assert outdir.is_dir()


def test_build_minimap2_index_missing_genome(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        index.build_minimap2_index(
            genome=tmp_path / "nope.fa",
            outdir=tmp_path / "idx",
        )


def test_build_minimap2_index_missing_tool(tmp_path: Path, monkeypatch):
    genome = tmp_path / "g.fa"
    genome.write_text(">chr1\nACGT\n")
    monkeypatch.setattr(index.shutil, "which", lambda _: None)
    with pytest.raises(RuntimeError, match="minimap2 not found"):
        index.build_minimap2_index(genome=genome, outdir=tmp_path / "idx")


def test_build_minimap2_index_invokes_correct_command(tmp_path: Path, monkeypatch):
    genome = tmp_path / "g.fa"
    genome.write_text(">chr1\nACGT\n")
    outdir = tmp_path / "idx"

    monkeypatch.setattr(index.shutil, "which", lambda _: "/fake/minimap2")
    captured = {}

    def fake_run(cmd, check):
        captured["cmd"] = cmd
        captured["check"] = check
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(index.subprocess, "run", fake_run)

    out = index.build_minimap2_index(
        genome=genome, outdir=outdir, threads=8, prefix="myidx"
    )
    assert out == outdir / "myidx.mmi"
    assert captured["check"] is True
    assert captured["cmd"] == [
        "minimap2", "-t", "8", "-x", "splice",
        "-d", str(outdir / "myidx.mmi"), str(genome),
    ]
    assert outdir.is_dir()
