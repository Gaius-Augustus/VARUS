"""Alignment of one batch (HISAT2 for short reads, minimap2 for long reads).

The HISAT2 path improves on ``legacy/Implementation/src/HISAT_Aligner.cpp``:

* No SAM intermediate. HISAT2's stdout is piped through ``samtools sort -O BAM``
  directly, eliminating one full read+write of the alignments.
* No separate ``samtools sort`` step afterwards.
* ``Log.final.out`` (HISAT2's stderr summary) is captured to a file in the
  batch directory so the legacy quality-parsing logic still works on it.

The minimap2 path follows the same pipe-to-``samtools sort`` pattern; minimap2
emits no end-of-run summary, so the quality gate is computed by scanning the
sorted BAM with pysam (see :func:`count_minimap2_quality`).

The output is ``<batch_dir>/Aligned.out.bam`` (coordinate-sorted), matching the
filename the legacy controller already expects after its own post-conversion.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class AlignmentResult:
    bam: Path
    log: Path


def _require(tool: str) -> None:
    if shutil.which(tool) is None:
        raise RuntimeError(f"{tool} not found on PATH")


def align_batch_hisat2(
    r1: Path,
    r2: Path | None,
    *,
    index_prefix: Path,
    batch_dir: Path,
    threads: int = 4,
    intron_db: Path | None = None,
    hisat2: str = "hisat2",
    samtools: str = "samtools",
) -> AlignmentResult:
    """Align one batch with HISAT2; write a sorted BAM.

    Parameters
    ----------
    r1, r2
        FASTA files. ``r2`` is ``None`` for single-end reads.
    index_prefix
        HISAT2 index prefix (e.g. ``<outdir>/genome/hisatidx``).
    batch_dir
        Output directory; ``Aligned.out.bam`` and ``Log.final.out`` go there.
    intron_db
        Optional path to a known-splice-site file in HISAT2's tab format
        (``--known-splicesite-infile``). If absent, alignment proceeds without
        a splice DB, like the very first batch in the legacy code.
    """
    _require(hisat2)
    _require(samtools)
    batch_dir.mkdir(parents=True, exist_ok=True)
    bam_out = batch_dir / "Aligned.out.bam"
    log_out = batch_dir / "Log.final.out"

    hisat_cmd: list[str] = [
        hisat2,
        "-p", str(threads),
        "-f",
        "-x", str(index_prefix),
    ]
    if r2 is None:
        hisat_cmd += ["-U", str(r1)]
    else:
        hisat_cmd += ["-1", str(r1), "-2", str(r2)]

    if intron_db is not None and intron_db.is_file():
        hisat_cmd += ["--known-splicesite-infile", str(intron_db)]

    sort_cmd = [
        samtools, "sort",
        "-@", str(max(1, threads - 1)),
        "-O", "BAM",
        "-o", str(bam_out),
    ]

    log.info("HISAT2 | samtools sort -> %s", bam_out)
    log.debug("hisat2 cmd: %s", " ".join(hisat_cmd))
    log.debug("samtools cmd: %s", " ".join(sort_cmd))

    with log_out.open("wb") as logf:
        hisat_proc = subprocess.Popen(
            hisat_cmd, stdout=subprocess.PIPE, stderr=logf
        )
        try:
            sort_proc = subprocess.Popen(
                sort_cmd, stdin=hisat_proc.stdout, stdout=subprocess.DEVNULL
            )
            # Allow hisat2 to receive SIGPIPE if samtools dies first.
            assert hisat_proc.stdout is not None
            hisat_proc.stdout.close()
            sort_rc = sort_proc.wait()
        finally:
            hisat_rc = hisat_proc.wait()

    if hisat_rc != 0:
        raise RuntimeError(
            f"hisat2 exited with status {hisat_rc}; see {log_out}"
        )
    if sort_rc != 0:
        raise RuntimeError(f"samtools sort exited with status {sort_rc}")

    return AlignmentResult(bam=bam_out, log=log_out)


LONGREAD_PRESETS: dict[str, list[str]] = {
    "pacbio": ["-ax", "splice"],
    "ont": ["-ax", "splice", "-uf", "-k14"],
}


def align_batch_minimap2(
    reads: Path,
    *,
    index: Path,
    batch_dir: Path,
    threads: int = 4,
    preset: str = "pacbio",
    junc_bed: Path | None = None,
    minimap2: str = "minimap2",
    samtools: str = "samtools",
) -> AlignmentResult:
    """Align one batch of long reads with minimap2; write a sorted BAM.

    Parameters
    ----------
    reads
        FASTA file. Long-read SRA runs are single-end; ``r2`` is never used.
    index
        ``.mmi`` index file from :func:`varus.index.build_minimap2_index`, or
        a genome FASTA (minimap2 will index on the fly).
    batch_dir
        Output directory; ``Aligned.out.bam`` and ``Log.minimap2.err`` go there.
    preset
        ``'pacbio'`` (Iso-Seq / HiFi, ``-ax splice``) or ``'ont'``
        (direct-RNA Nanopore, ``-ax splice -uf -k14``).
    junc_bed
        Optional BED12 of known junctions for minimap2 ``--junc-bed``. If
        absent, alignment proceeds without the hint, like the first batch.
    """
    if preset not in LONGREAD_PRESETS:
        raise ValueError(
            f"unknown long-read preset {preset!r}; "
            f"expected one of {list(LONGREAD_PRESETS)}"
        )
    _require(minimap2)
    _require(samtools)
    batch_dir.mkdir(parents=True, exist_ok=True)
    bam_out = batch_dir / "Aligned.out.bam"
    log_out = batch_dir / "Log.minimap2.err"

    mm2_cmd: list[str] = [
        minimap2,
        "-t", str(threads),
        *LONGREAD_PRESETS[preset],
    ]
    if junc_bed is not None and junc_bed.is_file():
        mm2_cmd += ["--junc-bed", str(junc_bed)]
    mm2_cmd += [str(index), str(reads)]

    sort_cmd = [
        samtools, "sort",
        "-@", str(max(1, threads - 1)),
        "-O", "BAM",
        "-o", str(bam_out),
    ]

    log.info("minimap2 (%s) | samtools sort -> %s", preset, bam_out)
    log.debug("minimap2 cmd: %s", " ".join(mm2_cmd))
    log.debug("samtools cmd: %s", " ".join(sort_cmd))

    with log_out.open("wb") as logf:
        mm2_proc = subprocess.Popen(
            mm2_cmd, stdout=subprocess.PIPE, stderr=logf
        )
        try:
            sort_proc = subprocess.Popen(
                sort_cmd, stdin=mm2_proc.stdout, stdout=subprocess.DEVNULL
            )
            assert mm2_proc.stdout is not None
            mm2_proc.stdout.close()
            sort_rc = sort_proc.wait()
        finally:
            mm2_rc = mm2_proc.wait()

    if mm2_rc != 0:
        raise RuntimeError(
            f"minimap2 exited with status {mm2_rc}; see {log_out}"
        )
    if sort_rc != 0:
        raise RuntimeError(f"samtools sort exited with status {sort_rc}")

    return AlignmentResult(bam=bam_out, log=log_out)


def count_minimap2_quality(
    bam_path: Path,
    *,
    min_mapq: int = 1,
) -> dict[str, float]:
    """Count primary, MAPQ ≥ ``min_mapq`` alignments in a minimap2 BAM.

    Returns the same dict shape as :func:`parse_hisat2_log`:
    ``{'num_uniq': float, 'uniq_pct': float}``.

    Note the denominator difference vs ``parse_hisat2_log``: HISAT2's parser
    divides by ``batch_size`` (input pairs/reads); here we divide by primary
    alignments observed. For long-read SRA spots-are-reads, this is the
    fraction of decoded reads that aligned uniquely — practically equivalent
    for the ``--min-uniq-pct`` quality gate but slightly differently calibrated.
    """
    if not bam_path.is_file():
        raise FileNotFoundError(bam_path)

    import pysam  # local import: pysam is in the optional [align] extra

    n_primary = 0
    n_uniq = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if (
                read.is_unmapped
                or read.is_secondary
                or read.is_supplementary
            ):
                continue
            n_primary += 1
            if read.mapping_quality >= min_mapq:
                n_uniq += 1

    if n_primary == 0:
        return {"num_uniq": 0.0, "uniq_pct": 0.0}
    return {
        "num_uniq": float(n_uniq),
        "uniq_pct": 100.0 * n_uniq / n_primary,
    }


def parse_hisat2_log(log_path: Path, batch_size: int) -> dict[str, float]:
    """Parse HISAT2's ``Log.final.out`` for unique-alignment percentage.

    Mirrors ``HISAT_Aligner::checkQuality``: the first occurrence of either
    "aligned concordantly exactly 1 time" (paired) or "aligned exactly 1 time"
    (single-end) gives the unique alignment count.
    """
    if not log_path.is_file():
        raise FileNotFoundError(log_path)
    num_uniq: int | None = None
    for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        s = line.lstrip()
        if num_uniq is None and (
            "aligned concordantly exactly 1 time" in s
            or "aligned exactly 1 time" in s
        ):
            try:
                num_uniq = int(s.split()[0])
            except (ValueError, IndexError):
                continue
            break
    if num_uniq is None:
        return {"num_uniq": 0, "uniq_pct": 0.0}
    return {
        "num_uniq": float(num_uniq),
        "uniq_pct": 100.0 * num_uniq / batch_size,
    }
