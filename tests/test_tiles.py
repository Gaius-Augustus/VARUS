"""Tests for varus.tiles. BAM-touching tests require pysam."""

from __future__ import annotations

from pathlib import Path

import pytest

from varus import tiles
from tests.conftest import requires_pysam


def test_count_umrs_rejects_zero_tile_size(tmp_path: Path):
    with pytest.raises(ValueError, match="tile_size"):
        tiles.count_umrs_per_tile(tmp_path / "x.bam", tile_size=0)


def test_count_bam_stats_rejects_zero_tile_size(tmp_path: Path):
    with pytest.raises(ValueError, match="tile_size"):
        tiles.count_bam_stats(tmp_path / "x.bam", tile_size=0)


@requires_pysam
def test_count_umrs_basic(tmp_path: Path):
    """One single-mapper, one multi-mapper, one paired-end pair on same tile."""
    import pysam

    bam_path = tmp_path / "t.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 100_000}],
    }
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        # r1: single alignment in tile 0 (POS=1) -> UMR, +1 to ('chr1', 0)
        r1 = pysam.AlignedSegment(out.header)
        r1.query_name = "r1"
        r1.query_sequence = "A" * 50
        r1.flag = 0
        r1.reference_id = 0
        r1.reference_start = 0  # SAM POS=1, tile = 1//5000 = 0
        r1.mapping_quality = 60
        r1.cigartuples = [(0, 50)]
        r1.query_qualities = pysam.qualitystring_to_array("I" * 50)
        out.write(r1)

        # r2: two alignments to different tiles -> not UMR
        for ref_start in (10, 6_000):  # tile 0 and tile 1
            r = pysam.AlignedSegment(out.header)
            r.query_name = "r2"
            r.query_sequence = "C" * 50
            r.flag = 0 if ref_start == 10 else 256  # secondary
            r.reference_id = 0
            r.reference_start = ref_start
            r.mapping_quality = 60
            r.cigartuples = [(0, 50)]
            r.query_qualities = pysam.qualitystring_to_array("I" * 50)
            out.write(r)

        # r3: paired-end pair both on tile 2 (POS in [10001, 15000]) -> UMR, +1
        for flag, ref_start in [(99, 10_010), (147, 10_200)]:
            r = pysam.AlignedSegment(out.header)
            r.query_name = "r3"
            r.query_sequence = "G" * 50
            r.flag = flag
            r.reference_id = 0
            r.reference_start = ref_start
            r.mapping_quality = 60
            r.cigartuples = [(0, 50)]
            r.query_qualities = pysam.qualitystring_to_array("I" * 50)
            out.write(r)

    counts = tiles.count_umrs_per_tile(bam_path, tile_size=5000)
    assert counts.get(("chr1", 0)) == 1   # r1
    assert counts.get(("chr1", 2)) == 1   # r3 pair
    # r2 split across tiles 0 and 1 -> not counted
    assert ("chr1", 1) not in counts
    # 3 distinct reads, but only 2 are UMR
    assert sum(counts.values()) == 2


@requires_pysam
def test_count_bam_stats_spliced_count(tmp_path: Path):
    """count_bam_stats tracks spliced reads (N-op in CIGAR)."""
    import pysam

    bam_path = tmp_path / "s.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 100_000}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        # unspliced read
        r1 = pysam.AlignedSegment(out.header)
        r1.query_name = "r1"
        r1.query_sequence = "A" * 50
        r1.flag = 0
        r1.reference_id = 0
        r1.reference_start = 0
        r1.mapping_quality = 60
        r1.cigartuples = [(0, 50)]   # 50M – no intron
        r1.query_qualities = pysam.qualitystring_to_array("I" * 50)
        out.write(r1)

        # spliced read with an N-op
        r2 = pysam.AlignedSegment(out.header)
        r2.query_name = "r2"
        r2.query_sequence = "G" * 60
        r2.flag = 0
        r2.reference_id = 0
        r2.reference_start = 0
        r2.mapping_quality = 60
        r2.cigartuples = [(0, 30), (3, 100), (0, 30)]  # 30M 100N 30M
        r2.query_qualities = pysam.qualitystring_to_array("I" * 60)
        out.write(r2)

    stats = tiles.count_bam_stats(bam_path, tile_size=5000)
    assert stats.n_reads == 2
    assert stats.n_spliced == 1   # only r2 has an N-op


@requires_pysam
def test_count_umrs_three_hits_same_tile_excluded(tmp_path: Path):
    """RNAread::UMR rejects reads with > 2 hits to a single tile."""
    import pysam

    bam_path = tmp_path / "t.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 100_000}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        for i, ref_start in enumerate([0, 100, 200]):
            r = pysam.AlignedSegment(out.header)
            r.query_name = "rN"
            r.query_sequence = "A" * 50
            r.flag = 0 if i == 0 else 256
            r.reference_id = 0
            r.reference_start = ref_start
            r.mapping_quality = 60
            r.cigartuples = [(0, 50)]
            r.query_qualities = pysam.qualitystring_to_array("I" * 50)
            out.write(r)

    counts = tiles.count_umrs_per_tile(bam_path, tile_size=5000)
    assert counts == {}
