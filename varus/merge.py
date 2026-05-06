"""Merge per-batch BAM files into a single VARUS.bam.

Replaces ``legacy/Implementation/scripts/finalMerge.sh`` and
``mergeAlignments.sh``. We call ``samtools merge`` directly which is simpler
and more portable than the hierarchical shell-script approach in the legacy
code.  The output BAM is coordinate-sorted because every input BAM was already
sorted by ``samtools sort`` in :func:`varus.align.align_batch_hisat2`.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path
from typing import Iterable

log = logging.getLogger(__name__)


def merge_bams(
    bam_files: Iterable[Path],
    out_bam: Path,
    *,
    threads: int = 4,
    samtools: str = "samtools",
    force: bool = True,
) -> Path:
    """Merge sorted BAM files into a single output BAM.

    Parameters
    ----------
    bam_files : Input BAM paths (must be coordinate-sorted).
    out_bam   : Destination path (overwritten when force=True).
    threads   : Threads passed to ``samtools merge -@``.
    samtools  : Executable name or absolute path.
    force     : Pass ``-f`` so an existing output is overwritten.

    Returns the path to the merged BAM.
    """
    if shutil.which(samtools) is None:
        raise RuntimeError(f"{samtools} not found on PATH")

    bam_list = [Path(b) for b in bam_files]
    if not bam_list:
        raise ValueError("merge_bams: no input BAM files provided")

    out_bam.parent.mkdir(parents=True, exist_ok=True)

    cmd = [samtools, "merge", "-@", str(threads)]
    if force:
        cmd.append("-f")
    cmd.append(str(out_bam))
    cmd.extend(str(b) for b in bam_list)

    log.info("samtools merge: %d input BAMs → %s", len(bam_list), out_bam)
    subprocess.run(cmd, check=True)
    return out_bam
