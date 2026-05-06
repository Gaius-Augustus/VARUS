"""Assign strand to intron hints using splice-site dinucleotides.

Reimplements ``filterIntronsFindStrand.pl`` (AUGUSTUS/BRAKER scripts).

For each intron at GFF coordinates (chrom, start, end) [1-based inclusive]:
- Donor dinucleotide: genome[start-1 : start+1]  (0-based, 2 chars)
- Acceptor dinucleotide: genome[end-2 : end]      (0-based, 2 chars)
- Concatenate → 4-char motif (lowercase)
- If motif ∈ allowed      → strand '+'
- If RC(motif) ∈ allowed  → strand '-'
- Otherwise               → intron is dropped (matching Perl behaviour)

The default allowed set matches the Perl default: gtag, gcag, atac.

Also provides :func:`write_hisat2_splice_sites` to convert an
:class:`~varus.introns.IntronCounts` into the tab-delimited format expected
by ``hisat2 --known-splicesite-infile``.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import FrozenSet

from varus.introns import IntronCounts, IntronKey

log = logging.getLogger(__name__)

DEFAULT_ALLOWED: FrozenSet[str] = frozenset({"gtag", "gcag", "atac"})

_COMPLEMENT = str.maketrans("acgtACGT", "tgcaTGCA")


def _rc4(motif: str) -> str:
    """Reverse-complement a 4-char splice-site motif."""
    return motif[::-1].translate(_COMPLEMENT).lower()


def assign_strand(
    introns: IntronCounts,
    genome_fasta: Path,
    allowed: FrozenSet[str] = DEFAULT_ALLOWED,
) -> IntronCounts:
    """Return a new IntronCounts with strand assigned; unrecognized introns dropped.

    Uses pyfaidx (optional dep) to fetch dinucleotides from the genome FASTA.
    Requires pyfaidx ≥ 0.7 and the FASTA to be indexed (``samtools faidx``
    or pyfaidx will build the index automatically on first run).
    """
    from pyfaidx import Fasta  # optional on Linux/macOS; not available on Windows

    fa = Fasta(str(genome_fasta), as_raw=True)
    new: dict[IntronKey, int] = {}
    n_kept = n_dropped = 0

    for (chrom, start, end, _), mult in introns.counts.items():
        if chrom not in fa:
            log.warning("Chrom %s absent from genome FASTA; dropping intron", chrom)
            n_dropped += 1
            continue
        # 0-based half-open slices; pyfaidx follows Python convention
        donor = str(fa[chrom][start - 1 : start + 1]).lower()
        acceptor = str(fa[chrom][end - 2 : end]).lower()
        motif = donor + acceptor
        if motif in allowed:
            new[(chrom, start, end, "+")] = mult
            n_kept += 1
        elif _rc4(motif) in allowed:
            new[(chrom, start, end, "-")] = mult
            n_kept += 1
        else:
            n_dropped += 1

    log.info("assign_strand: %d kept, %d dropped", n_kept, n_dropped)
    return IntronCounts(new)


def write_hisat2_splice_sites(introns: IntronCounts, path: Path) -> int:
    """Write a HISAT2 --known-splicesite-infile from stranded introns.

    Each line::

        chrom \\t donor_0based \\t acceptor_0based \\t strand

    where ``donor`` = last base of left exon (0-based) = start - 2,
    and ``acceptor`` = first base of right exon (0-based) = end.

    Introns with strand '.' are silently skipped.
    Returns the number of records written.
    """
    n = 0
    with path.open("w", encoding="utf-8") as f:
        for (chrom, start, end, strand), _ in sorted(introns.counts.items()):
            if strand not in ("+", "-"):
                continue
            donor = start - 2   # 0-based last base before the intron
            acceptor = end      # 0-based first base after the intron
            if donor < 0:
                continue
            f.write(f"{chrom}\t{donor}\t{acceptor}\t{strand}\n")
            n += 1
    log.info("Wrote %d splice sites to %s", n, path)
    return n
