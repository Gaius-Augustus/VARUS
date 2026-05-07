"""Build a HISAT2 or minimap2 index for a genome FASTA.

Replaces the index-building branch of ``legacy/runVARUS.pl``.

For short reads we invoke ``hisat2-build`` directly. The output prefix matches
what the legacy aligner expected (``<outdir>/hisatidx``), so a ``varus run``
task can point at ``--index <outdir>`` and HISAT2 will resolve
``<outdir>/hisatidx.*.ht2``.

For long reads (``--longreads``) we invoke ``minimap2 -d`` with the ``splice``
preset. minimap2 can index on the fly, but pre-building a ``.mmi`` saves the
seed-table construction cost on every batch. The ``splice`` index works for
both PacBio Iso-Seq and ONT direct-RNA; only the alignment-time flags differ.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path

log = logging.getLogger(__name__)


def build_hisat2_index(
    genome: Path,
    outdir: Path,
    threads: int = 4,
    prefix: str = "hisatidx",
) -> Path:
    """Build a HISAT2 index. Returns the index *prefix* path."""
    if not genome.is_file():
        raise FileNotFoundError(f"genome FASTA not found: {genome}")
    if shutil.which("hisat2-build") is None:
        raise RuntimeError("hisat2-build not found on PATH")

    outdir.mkdir(parents=True, exist_ok=True)
    idx_prefix = outdir / prefix

    cmd = [
        "hisat2-build",
        "-p", str(threads),
        str(genome),
        str(idx_prefix),
    ]
    log.info("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)
    log.info("HISAT2 index written with prefix %s", idx_prefix)
    return idx_prefix


def build_minimap2_index(
    genome: Path,
    outdir: Path,
    threads: int = 4,
    prefix: str = "mm2idx",
) -> Path:
    """Build a minimap2 splice index. Returns the ``.mmi`` *file* path.

    Unlike HISAT2, minimap2 produces a single index file rather than a 6-file
    set, so the returned path is the ``.mmi`` itself rather than a stem.
    """
    if not genome.is_file():
        raise FileNotFoundError(f"genome FASTA not found: {genome}")
    if shutil.which("minimap2") is None:
        raise RuntimeError("minimap2 not found on PATH")

    outdir.mkdir(parents=True, exist_ok=True)
    idx_path = outdir / f"{prefix}.mmi"

    cmd = [
        "minimap2",
        "-t", str(threads),
        "-x", "splice",
        "-d", str(idx_path),
        str(genome),
    ]
    log.info("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)
    log.info("minimap2 index written to %s", idx_path)
    return idx_path
