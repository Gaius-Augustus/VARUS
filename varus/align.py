"""HISAT2 alignment of one batch, output coordinate-sorted BAM.

Improvements over ``legacy/Implementation/src/HISAT_Aligner.cpp``:

* No SAM intermediate. HISAT2's stdout is piped through ``samtools sort -O BAM``
  directly, eliminating one full read+write of the alignments.
* No separate ``samtools sort`` step afterwards.
* ``Log.final.out`` (HISAT2's stderr summary) is captured to a file in the
  batch directory so the legacy quality-parsing logic still works on it.

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
