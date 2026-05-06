"""Count uniquely-mapped reads per genome tile from a sorted BAM.

The tile is the proxy for a "transcribed unit" in the VARUS algorithm: the
genome is partitioned into non-overlapping windows of ``tile_size`` bp
(default 5 kb). For each read, all of its alignments are pooled by tile; the
read is counted iff:

* it mapped to exactly one tile, **and**
* that tile saw at most two alignment records for the read (so a paired-end
  pair both landing on the same tile still counts once, but a multi-mapper
  is excluded).

This matches ``RNAread::UMR()`` in the legacy code (see
``legacy/Implementation/src/RNAread.cpp``).

Coordinate convention
---------------------
The legacy reads the SAM POS column directly (1-based). pysam exposes
``reference_start`` as 0-based, so the equivalent tile index is
``(reference_start + 1) // tile_size``.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Tuple

log = logging.getLogger(__name__)

Tile = Tuple[str, int]  # (chromosome, tile_index)


@dataclass
class BAMStats:
    """Summary statistics from a single BAM-pass."""
    umr_counts: Dict[Tile, int]
    n_reads: int           # distinct mapped read names
    n_spliced: int         # reads with at least one N (intron) op


def count_umrs_per_tile(bam_path: Path, tile_size: int) -> Dict[Tile, int]:
    """Return ``{(chrom, tile_idx): umr_count}`` for the given BAM."""
    return count_bam_stats(bam_path, tile_size).umr_counts


def count_bam_stats(bam_path: Path, tile_size: int) -> BAMStats:
    """Single-pass BAM scan: UMR counts per tile + spliced-read count.

    Returns a :class:`BAMStats` with:
    - ``umr_counts``: ``{(chrom, tile_idx): n}`` — uniquely-mapped read count.
    - ``n_reads``: total distinct mapped read names.
    - ``n_spliced``: reads with at least one CIGAR N-op (intron).
    """
    if tile_size <= 0:
        raise ValueError("tile_size must be > 0")

    import pysam  # local import: pysam is an optional install (extras "align")

    per_read: dict[str, dict[Tile, int]] = defaultdict(lambda: defaultdict(int))
    spliced_reads: set[str] = set()

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.reference_name is None:
                continue
            tile = (read.reference_start + 1) // tile_size
            key: Tile = (read.reference_name, tile)
            per_read[read.query_name][key] += 1
            if read.cigartuples and any(op == 3 for op, _ in read.cigartuples):
                spliced_reads.add(read.query_name)

    umr_counts: Dict[Tile, int] = defaultdict(int)
    for read_name, tiles in per_read.items():
        if len(tiles) != 1:
            continue
        ((tile, hits),) = tiles.items()
        if hits <= 2:  # one pair landing on the same tile is allowed
            umr_counts[tile] += 1

    n_reads = len(per_read)
    n_spliced = len(spliced_reads)
    log.info(
        "BAM %s: %d UMRs across %d tiles; %d/%d reads spliced",
        bam_path, sum(umr_counts.values()), len(umr_counts), n_spliced, n_reads,
    )
    return BAMStats(
        umr_counts=dict(umr_counts),
        n_reads=n_reads,
        n_spliced=n_spliced,
    )
