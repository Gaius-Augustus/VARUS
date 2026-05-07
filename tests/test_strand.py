"""Tests for varus.strand: assign_strand and write_hisat2_splice_sites.

The pyfaidx-dependent assign_strand is gated on requires_pysam for the same
reason as other tests that touch the Linux-only native deps stack. pyfaidx
itself builds everywhere, but in practice the strand test writes a synthetic
FASTA which requires pyfaidx too, so we use the same guard.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from varus.introns import IntronCounts
from varus.strand import (
    DEFAULT_ALLOWED, _rc4, write_hisat2_splice_sites, write_minimap2_junc_bed,
)

try:
    from tests.conftest import requires_pysam
except ImportError:
    requires_pysam = pytest.mark.skip(reason="conftest not found")


# ---------------------------------------------------------------------------
# Pure-Python helpers
# ---------------------------------------------------------------------------

def test_rc4_canonical_motifs():
    """RC of the canonical GT-AG motif pair should reproduce itself on RC strand."""
    # GT-AG on + strand encodes as "gtag"; RC("gtag") should be "ctac"
    # and "ctac" RC should be "gtag" again
    assert _rc4("gtag") == "ctac"
    assert _rc4("ctac") == "gtag"
    assert _rc4("gcag") == "ctgc"
    assert _rc4("atac") == "gtat"


def test_write_hisat2_splice_sites_plus_strand(tmp_path: Path):
    """Plus-strand intron → donor = start-2, acceptor = end."""
    introns = IntronCounts({("chr1", 10, 50, "+"): 3})
    out = tmp_path / "ss.txt"
    n = write_hisat2_splice_sites(introns, out)
    assert n == 1
    line = out.read_text().strip()
    parts = line.split("\t")
    assert parts[0] == "chr1"
    assert int(parts[1]) == 8   # donor = start - 2 = 10 - 2 = 8
    assert int(parts[2]) == 50  # acceptor = end = 50
    assert parts[3] == "+"


def test_write_hisat2_splice_sites_minus_strand(tmp_path: Path):
    introns = IntronCounts({("chr2", 100, 200, "-"): 1})
    out = tmp_path / "ss.txt"
    n = write_hisat2_splice_sites(introns, out)
    assert n == 1
    line = out.read_text().strip()
    parts = line.split("\t")
    assert parts[3] == "-"
    assert int(parts[1]) == 98   # 100 - 2
    assert int(parts[2]) == 200  # end


def test_write_hisat2_splice_sites_skips_dot_strand(tmp_path: Path):
    """Introns with strand '.' are silently skipped."""
    introns = IntronCounts({
        ("chr1", 10, 50, "."): 5,
        ("chr1", 60, 90, "+"): 2,
    })
    out = tmp_path / "ss.txt"
    n = write_hisat2_splice_sites(introns, out)
    assert n == 1
    assert "chr1\t58\t90\t+" in out.read_text()


def test_write_hisat2_splice_sites_skips_donor_below_zero(tmp_path: Path):
    """Intron starting at position 1 (1-based) would give donor=-1; skip it."""
    introns = IntronCounts({("chr1", 1, 10, "+"): 1})
    out = tmp_path / "ss.txt"
    n = write_hisat2_splice_sites(introns, out)
    assert n == 0


def test_write_minimap2_junc_bed_plus_strand(tmp_path: Path):
    """Plus-strand intron → BED12 with two 1bp blocks flanking the intron."""
    introns = IntronCounts({("chr1", 10, 50, "+"): 3})
    out = tmp_path / "junc.bed"
    n = write_minimap2_junc_bed(introns, out)
    assert n == 1
    parts = out.read_text().strip().split("\t")
    assert len(parts) == 12
    assert parts[0] == "chr1"
    assert int(parts[1]) == 8     # bed_start = start - 2 = 8
    assert int(parts[2]) == 51    # bed_end = end + 1 = 51
    assert parts[5] == "+"
    assert int(parts[4]) == 3     # score = multiplicity (clipped at 1000)
    assert int(parts[9]) == 2     # blockCount
    assert parts[10] == "1,1"     # blockSizes
    # blockStarts: first at 0, second at bed_end - bed_start - 1 = 42
    assert parts[11] == "0,42"


def test_write_minimap2_junc_bed_minus_strand(tmp_path: Path):
    introns = IntronCounts({("chr2", 100, 200, "-"): 1})
    out = tmp_path / "junc.bed"
    n = write_minimap2_junc_bed(introns, out)
    assert n == 1
    parts = out.read_text().strip().split("\t")
    assert parts[0] == "chr2"
    assert parts[5] == "-"
    assert int(parts[1]) == 98
    assert int(parts[2]) == 201


def test_write_minimap2_junc_bed_skips_dot_strand(tmp_path: Path):
    introns = IntronCounts({
        ("chr1", 10, 50, "."): 5,
        ("chr1", 60, 90, "+"): 2,
    })
    out = tmp_path / "junc.bed"
    n = write_minimap2_junc_bed(introns, out)
    assert n == 1
    text = out.read_text()
    # Only the + record should appear.
    assert "\t+\t" in text
    assert "\t.\t" not in text


def test_write_minimap2_junc_bed_skips_donor_below_zero(tmp_path: Path):
    introns = IntronCounts({("chr1", 1, 10, "+"): 1})
    out = tmp_path / "junc.bed"
    n = write_minimap2_junc_bed(introns, out)
    assert n == 0


def test_write_minimap2_junc_bed_clips_score_at_1000(tmp_path: Path):
    introns = IntronCounts({("chr1", 10, 50, "+"): 5_000})
    out = tmp_path / "junc.bed"
    write_minimap2_junc_bed(introns, out)
    parts = out.read_text().strip().split("\t")
    assert int(parts[4]) == 1000


# ---------------------------------------------------------------------------
# pyfaidx-dependent: assign_strand
# ---------------------------------------------------------------------------

@requires_pysam
def test_assign_strand_gt_ag(tmp_path: Path):
    """GT-AG canonical intron on + strand."""
    # Build a simple FASTA: 100 Ns then GT...AG at positions 11-50 (1-based)
    # genome[10:12] = "GT"  (0-based = positions 11-12 in 1-based = intron donor)
    # genome[48:50] = "AG"  (0-based = positions 49-50 in 1-based = intron acceptor)
    genome_seq = "N" * 10 + "GT" + "N" * 36 + "AG" + "N" * 100
    fasta = tmp_path / "genome.fa"
    fasta.write_text(f">chr1\n{genome_seq}\n")

    # Intron: start=11, end=50 (1-based, inclusive)
    # donor = genome[start-1 : start+1] = genome[10:12] = "GT"
    # acceptor = genome[end-2 : end] = genome[48:50] = "AG"
    # motif = "gtag" → strand "+"
    introns = IntronCounts({("chr1", 11, 50, "."): 7})
    from varus.strand import assign_strand
    result = assign_strand(introns, fasta)
    assert result.counts == {("chr1", 11, 50, "+"): 7}


@requires_pysam
def test_assign_strand_minus(tmp_path: Path):
    """CT-AC motif → minus strand (RC of GT-AG)."""
    # RC of "gtag" is "ctac"
    # genome[10:12] = "CT", genome[48:50] = "AC"
    genome_seq = "N" * 10 + "CT" + "N" * 36 + "AC" + "N" * 100
    fasta = tmp_path / "genome.fa"
    fasta.write_text(f">chr1\n{genome_seq}\n")

    introns = IntronCounts({("chr1", 11, 50, "."): 4})
    from varus.strand import assign_strand
    result = assign_strand(introns, fasta)
    assert result.counts == {("chr1", 11, 50, "-"): 4}


@requires_pysam
def test_assign_strand_unknown_dropped(tmp_path: Path):
    """Introns whose splice site doesn't match any allowed motif are dropped."""
    # genome[10:12] = "AA", genome[48:50] = "TT"  → motif "aatt" not in allowed
    genome_seq = "N" * 10 + "AA" + "N" * 36 + "TT" + "N" * 100
    fasta = tmp_path / "genome.fa"
    fasta.write_text(f">chr1\n{genome_seq}\n")

    introns = IntronCounts({("chr1", 11, 50, "."): 2})
    from varus.strand import assign_strand
    result = assign_strand(introns, fasta)
    assert result.counts == {}
