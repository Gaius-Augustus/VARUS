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


# ---------------------------------------------------------------------------
# minimap2 long-read path
# ---------------------------------------------------------------------------

def test_align_batch_minimap2_unknown_preset(tmp_path: Path):
    with pytest.raises(ValueError, match="unknown long-read preset"):
        align.align_batch_minimap2(
            reads=tmp_path / "r.fa",
            index=tmp_path / "g.mmi",
            batch_dir=tmp_path / "b",
            preset="bogus",
        )


def test_align_batch_minimap2_requires_minimap2(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(align.shutil, "which", lambda t: None)
    with pytest.raises(RuntimeError, match="minimap2 not found"):
        align.align_batch_minimap2(
            reads=tmp_path / "r.fa",
            index=tmp_path / "g.mmi",
            batch_dir=tmp_path / "b",
        )


def test_align_batch_minimap2_requires_samtools(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(
        align.shutil, "which",
        lambda t: "/fake/minimap2" if t == "minimap2" else None,
    )
    with pytest.raises(RuntimeError, match="samtools not found"):
        align.align_batch_minimap2(
            reads=tmp_path / "r.fa",
            index=tmp_path / "g.mmi",
            batch_dir=tmp_path / "b",
        )


def test_align_batch_minimap2_command_construction(tmp_path: Path, monkeypatch):
    """Capture the cmd lists handed to subprocess.Popen."""
    monkeypatch.setattr(
        align.shutil, "which",
        lambda t: f"/fake/{t}",
    )
    captured: list[list[str]] = []

    class FakeProc:
        def __init__(self, cmd, **kwargs):
            captured.append(cmd)
            self.stdout = kwargs.get("stdin") and None
            # mm2 needs a stdout pipe for the sort to read; fake one
            if "stdout" in kwargs and kwargs["stdout"] is align.subprocess.PIPE:
                # provide a closeable object
                import io
                self.stdout = io.BytesIO()
            self._returncode = 0

        def wait(self):
            return self._returncode

    monkeypatch.setattr(align.subprocess, "Popen", FakeProc)

    junc = tmp_path / "junc.bed"
    junc.write_text("chr1\t0\t10\tj\t1\t+\t0\t10\t0\t2\t1,1\t0,9\n")

    align.align_batch_minimap2(
        reads=tmp_path / "r.fa",
        index=tmp_path / "g.mmi",
        batch_dir=tmp_path / "b",
        threads=8,
        preset="ont",
        junc_bed=junc,
    )

    assert len(captured) == 2  # mm2 then sort
    mm2_cmd, sort_cmd = captured
    assert mm2_cmd[0] == "minimap2"
    assert mm2_cmd[1:3] == ["-t", "8"]
    # ONT preset: -ax splice -uf -k14
    assert "-ax" in mm2_cmd and "splice" in mm2_cmd
    assert "-uf" in mm2_cmd and "-k14" in mm2_cmd
    assert "--junc-bed" in mm2_cmd
    assert str(junc) in mm2_cmd
    # Reads file is the last positional arg.
    assert mm2_cmd[-1] == str(tmp_path / "r.fa")
    assert mm2_cmd[-2] == str(tmp_path / "g.mmi")
    assert sort_cmd[0] == "samtools"
    assert "sort" in sort_cmd


def test_align_batch_minimap2_pacbio_preset(tmp_path: Path, monkeypatch):
    """PacBio preset is just '-ax splice', no -uf / -k14."""
    monkeypatch.setattr(align.shutil, "which", lambda t: f"/fake/{t}")
    captured: list[list[str]] = []

    class FakeProc:
        def __init__(self, cmd, **kwargs):
            captured.append(cmd)
            import io
            self.stdout = io.BytesIO() if kwargs.get("stdout") is align.subprocess.PIPE else None

        def wait(self):
            return 0

    monkeypatch.setattr(align.subprocess, "Popen", FakeProc)

    align.align_batch_minimap2(
        reads=tmp_path / "r.fa",
        index=tmp_path / "g.mmi",
        batch_dir=tmp_path / "b",
        preset="pacbio",
    )

    mm2_cmd = captured[0]
    assert "-uf" not in mm2_cmd
    assert "-k14" not in mm2_cmd
    assert "splice" in mm2_cmd


# count_minimap2_quality is gated on pysam since it scans a real BAM.
try:
    from tests.conftest import requires_pysam
except ImportError:
    requires_pysam = pytest.mark.skip(reason="conftest not found")


@requires_pysam
def test_count_minimap2_quality(tmp_path: Path):
    """Build a tiny BAM and check primary/MAPQ accounting."""
    import pysam

    header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"LN": 1000, "SN": "chr1"}],
    }
    bam_path = tmp_path / "test.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        # 3 primary, MAPQ=60
        for i in range(3):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.flag = 0
            a.reference_id = 0
            a.reference_start = i * 10
            a.mapping_quality = 60
            a.cigar = [(0, 10)]
            a.query_sequence = "A" * 10
            a.query_qualities = pysam.qualitystring_to_array("I" * 10)
            bam.write(a)
        # 1 primary MAPQ=0
        a = pysam.AlignedSegment()
        a.query_name = "r_low"
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 100
        a.mapping_quality = 0
        a.cigar = [(0, 10)]
        a.query_sequence = "A" * 10
        a.query_qualities = pysam.qualitystring_to_array("I" * 10)
        bam.write(a)
        # 1 secondary
        a = pysam.AlignedSegment()
        a.query_name = "r_sec"
        a.flag = 256
        a.reference_id = 0
        a.reference_start = 200
        a.mapping_quality = 60
        a.cigar = [(0, 10)]
        a.query_sequence = "A" * 10
        a.query_qualities = pysam.qualitystring_to_array("I" * 10)
        bam.write(a)
        # 1 unmapped
        a = pysam.AlignedSegment()
        a.query_name = "r_unmap"
        a.flag = 4
        a.reference_id = -1
        bam.write(a)

    pysam.index(str(bam_path))
    stats = align.count_minimap2_quality(bam_path, min_mapq=1)
    # 4 primaries (3 high-MAPQ + 1 MAPQ=0), 3 uniques.
    assert stats["num_uniq"] == 3.0
    assert abs(stats["uniq_pct"] - 75.0) < 1e-6


def test_count_minimap2_quality_missing_file(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        align.count_minimap2_quality(tmp_path / "nope.bam")
