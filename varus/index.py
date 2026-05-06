"""Build a HISAT2 index for a genome FASTA.

Replaces the index-building branch of ``legacy/runVARUS.pl``.

We invoke ``hisat2-build`` directly. The output prefix matches what the legacy
aligner expected (``<outdir>/hisatidx``), so a ``varus run`` task can point at
``--index <outdir>`` and HISAT2 will resolve ``<outdir>/hisatidx.*.ht2``.
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
