"""Tests for varus.introns.

The pure-Python paths (GFF read/write, IntronCounts.merge) are tested without
pysam. The BAM-walking path needs pysam and is gated on ``requires_pysam``.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from varus import introns
from tests.conftest import requires_pysam


def test_intron_counts_merge():
    a = introns.IntronCounts({("chr1", 100, 200, "."): 3})
    b = introns.IntronCounts({
        ("chr1", 100, 200, "."): 2,
        ("chr2", 50, 80, "+"): 1,
    })
    merged = a.merge(b)
    assert merged.counts == {
        ("chr1", 100, 200, "."): 5,
        ("chr2", 50, 80, "+"): 1,
    }


def test_write_then_read_introns_gff_roundtrip(tmp_path: Path):
    counts = introns.IntronCounts({
        ("chr1", 100, 200, "."): 5,
        ("chr1", 50, 80, "+"): 2,
        ("chr2", 10, 20, "."): 1,
    })
    out = tmp_path / "cumintrons.gff"
    n = introns.write_introns_gff(counts, out)
    assert n == 3

    # Sorted by chrom, start, end, strand
    lines = out.read_text().splitlines()
    assert lines[0].startswith("chr1\tb2h\tintron\t50\t80\t2\t+")
    assert lines[1].startswith("chr1\tb2h\tintron\t100\t200\t5\t.")
    assert lines[2].startswith("chr2\tb2h\tintron\t10\t20\t1\t.")
    # mult attribute prefix matches join_mult_hints output.
    assert "mult=5;pri=4;src=E" in lines[1]

    rt = introns.read_introns_gff(out)
    assert rt.counts == counts.counts


def test_read_introns_gff_skips_non_intron_lines(tmp_path: Path):
    p = tmp_path / "g.gff"
    p.write_text(
        "# comment\n"
        "\n"
        "chr1\tb2h\texon\t1\t100\t.\t+\t.\tfoo=bar\n"
        "chr1\tb2h\tintron\t101\t200\t3\t.\t.\tmult=3;pri=4;src=E\n"
    )
    rt = introns.read_introns_gff(p)
    assert rt.counts == {("chr1", 101, 200, "."): 3}


@requires_pysam
def test_extract_introns_from_bam_simple(tmp_path: Path):
    """Build a tiny BAM with one spliced read and verify the intron coords."""
    import pysam

    bam_path = tmp_path / "test.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 1000}],
    }
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        a = pysam.AlignedSegment(out.header)
        a.query_name = "r1"
        a.query_sequence = "A" * 60
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 99      # 0-based; SAM POS = 100
        a.mapping_quality = 60
        # 30M 200N 30M => intron at SAM pos 100 + 30 = 130, length 200
        # intron_start (1-based) = ref_pos + 1 after 30M consumed = 99+30 + 1 = 130
        # intron_end   (1-based) = 130 + 200 - 1 = 329
        a.cigartuples = [(0, 30), (3, 200), (0, 30)]
        a.query_qualities = pysam.qualitystring_to_array("I" * 60)
        out.write(a)

    pysam.index(str(bam_path))
    counts = introns.extract_introns_from_bam(bam_path)
    assert counts.counts == {("chr1", 130, 329, "."): 1}


@requires_pysam
def test_extract_introns_unmapped_skipped(tmp_path: Path):
    import pysam

    bam_path = tmp_path / "u.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        a = pysam.AlignedSegment(out.header)
        a.query_name = "u1"
        a.query_sequence = "A" * 60
        a.flag = 4  # unmapped
        a.reference_id = -1
        a.reference_start = -1
        out.write(a)

    counts = introns.extract_introns_from_bam(bam_path)
    assert counts.counts == {}
