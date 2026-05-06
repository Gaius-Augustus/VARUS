"""Extract spliced-alignment introns from a BAM file.

Replaces ``bam2hints --intronsonly`` (AUGUSTUS) and ``join_mult_hints.pl``
(BRAKER) for the per-batch and cumulative-merge cases.

CIGAR semantics
---------------
The ``N`` operation in the SAM CIGAR string means "skipped from the reference"
and is the standard encoding of an intron. For an alignment starting at
``reference_start`` (0-based pysam) we walk the CIGAR; whenever we hit ``N``
the implied intron coordinates are::

    intron_start (1-based, inclusive) = ref_pos + 1
    intron_end   (1-based, inclusive) = ref_pos + N_length

where ``ref_pos`` is the 0-based reference offset accumulated from prior
``M/D/=/X`` operations.

Strand
------
The legacy chain emits strand ``.`` from ``bam2hints --intronsonly`` and
fills it in afterwards with ``filterIntronsFindStrand.pl`` using the genome
FASTA. We do the same: this function leaves strand as ``.`` so the strand
assignment step (Phase 2.5) stays a separate, testable concern.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, Tuple

log = logging.getLogger(__name__)

# CIGAR ops that consume the reference (BAM_CMATCH, BAM_CDEL, BAM_CEQUAL,
# BAM_CDIFF). N is handled separately.
_REF_CONSUMING = frozenset({0, 2, 7, 8})

# (chromosome, start_1based_inclusive, end_1based_inclusive, strand)
IntronKey = Tuple[str, int, int, str]


@dataclass(frozen=True)
class IntronCounts:
    """Multiplicity-counted introns.

    Attribute name matches GFF semantics; ``items()`` lets the caller iterate
    in insertion order which is useful for downstream sorting tests.
    """
    counts: Dict[IntronKey, int]

    def __len__(self) -> int:
        return len(self.counts)

    def items(self) -> Iterator[Tuple[IntronKey, int]]:
        return iter(self.counts.items())

    def merge(self, other: "IntronCounts") -> "IntronCounts":
        out: Dict[IntronKey, int] = dict(self.counts)
        for k, v in other.counts.items():
            out[k] = out.get(k, 0) + v
        return IntronCounts(out)


def extract_introns_from_bam(bam_path: Path) -> IntronCounts:
    """Walk a BAM file, return per-intron multiplicity.

    Equivalent to ``bam2hints --intronsonly`` followed by
    ``join_mult_hints.pl`` over the resulting GFF.
    """
    import pysam  # extras "align"

    counts: Dict[IntronKey, int] = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.cigartuples is None:
                continue
            chrom = read.reference_name
            if chrom is None:
                continue
            ref_pos = read.reference_start  # 0-based
            for op, length in read.cigartuples:
                if op == 3:  # BAM_CREF_SKIP, the intron 'N' op
                    intron_start = ref_pos + 1
                    intron_end = ref_pos + length
                    key: IntronKey = (chrom, intron_start, intron_end, ".")
                    counts[key] = counts.get(key, 0) + 1
                    ref_pos += length
                elif op in _REF_CONSUMING:
                    ref_pos += length
                # I (1), S (4), H (5), P (6) -- don't consume ref

    log.info("BAM %s: %d distinct introns", bam_path, len(counts))
    return IntronCounts(counts)


def write_introns_gff(introns: IntronCounts, path: Path) -> int:
    """Write a GFF compatible with the legacy ``cumintrons.gff`` format.

    Output lines look like::

        chr1\tb2h\tintron\t101\t200\t15\t.\t.\tmult=15;pri=4;src=E

    Sort order matches what ``join_mult_hints.pl`` expects upstream of itself
    (chrom, start, end, strand) so the result is stable across runs.

    Returns the number of intron records written.
    """
    sorted_keys = sorted(introns.counts.keys(), key=lambda k: (k[0], k[1], k[2], k[3]))
    n = 0
    with path.open("w", encoding="utf-8") as f:
        for k in sorted_keys:
            chrom, start, end, strand = k
            mult = introns.counts[k]
            f.write(
                f"{chrom}\tb2h\tintron\t{start}\t{end}\t{mult}\t{strand}\t.\t"
                f"mult={mult};pri=4;src=E\n"
            )
            n += 1
    return n


def read_introns_gff(path: Path) -> IntronCounts:
    """Round-trip the format written by :func:`write_introns_gff`.

    Useful for resuming a run from disk and for testing.
    """
    counts: Dict[IntronKey, int] = {}
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "intron":
                continue
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attrs = parts[8]
            mult = 1
            for attr in attrs.split(";"):
                if attr.startswith("mult="):
                    try:
                        mult = int(attr[5:])
                    except ValueError:
                        mult = 1
                    break
            key: IntronKey = (chrom, start, end, strand)
            counts[key] = counts.get(key, 0) + mult
    return IntronCounts(counts)
