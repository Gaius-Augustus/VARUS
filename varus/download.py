"""Download spots from SRA.

The online algorithm requires *spot-range* downloads (read N..X of run R), not
whole-run downloads. ``fasterq-dump`` (sra-toolkit 3.x) does not support range
extraction; it is a full-run multithreaded optimisation. We therefore use:

* ``fastq-dump -N <n> -X <x> --fasta`` for batched range downloads (the proven
  legacy path);
* ``fasterq-dump --threads N --fasta`` for full-run downloads, used only when
  the controller is asked to ``--bootstrap-all``.

If a future sra-toolkit release adds range support to ``fasterq-dump`` we can
swap the batch backend without changing the controller.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

log = logging.getLogger(__name__)

# Number of retry attempts for transient SRA failures. The legacy code retries
# twice (3 attempts total) inside Downloader::getBatch.
DEFAULT_RETRIES = 3


@dataclass(frozen=True)
class BatchPaths:
    """Paths returned by ``download_batch``.

    For paired-end runs both ``r1`` and ``r2`` are populated. For single-end
    only ``r1`` is set; ``r2`` is ``None``.
    """
    r1: Path
    r2: Path | None
    batch_dir: Path

    def as_list(self) -> list[Path]:
        return [self.r1] if self.r2 is None else [self.r1, self.r2]


def batch_dir_for(outdir: Path, accession: str, n: int, x: int) -> Path:
    """Per-batch directory layout, identical to legacy ``Aligner::batchDir``::

        <outdir>/batches/<acc>/N<n>X<x>/
    """
    return outdir / "batches" / accession / f"N{n}X{x}"


def _require(tool: str) -> str:
    path = shutil.which(tool)
    if path is None:
        raise RuntimeError(f"{tool} not found on PATH")
    return path


def download_batch(
    accession: str,
    n: int,
    x: int,
    paired: bool,
    outdir: Path,
    *,
    retries: int = DEFAULT_RETRIES,
    fastq_dump: str = "fastq-dump",
) -> BatchPaths:
    """Download spots [n, x] of an SRA run as FASTA.

    Replicates ``legacy/Implementation/src/Downloader.cpp``'s
    ``shellCommand`` with ``fasta 120`` and ``--split-files`` for paired runs.

    Returns paths to the resulting FASTA file(s) inside the batch directory.
    """
    if shutil.which(fastq_dump) is None:
        raise RuntimeError(f"{fastq_dump} not found on PATH")

    bdir = batch_dir_for(outdir, accession, n, x)
    bdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        fastq_dump,
        "-N", str(n),
        "-X", str(x),
        "-O", str(bdir),
        "--fasta", "120",
    ]
    if paired:
        cmd.append("--split-files")
    cmd.append(accession)

    last_err: subprocess.CalledProcessError | None = None
    for attempt in range(1, retries + 1):
        log.info("fastq-dump %s N=%d X=%d (attempt %d/%d)",
                 accession, n, x, attempt, retries)
        try:
            subprocess.run(cmd, check=True)
            break
        except subprocess.CalledProcessError as e:
            last_err = e
            log.warning("fastq-dump failed (rc=%d) for %s N=%d X=%d",
                        e.returncode, accession, n, x)
    else:
        raise RuntimeError(
            f"fastq-dump failed after {retries} attempts for "
            f"{accession} N={n} X={x}: {last_err}"
        )

    if paired:
        r1 = bdir / f"{accession}_1.fasta"
        r2 = bdir / f"{accession}_2.fasta"
        # SRR097898 has a 3-file split (technical barcode in the middle); legacy
        # code handles that by using files[0] and files[-1]. We mirror that.
        if not r2.is_file():
            fastas = sorted(bdir.glob("*.fasta"))
            if len(fastas) >= 2:
                r1 = fastas[0]
                r2 = fastas[-1]
        if not r1.is_file() or not r2.is_file():
            raise RuntimeError(
                f"paired FASTA files missing in {bdir} for {accession}"
            )
        return BatchPaths(r1=r1, r2=r2, batch_dir=bdir)

    r1 = bdir / f"{accession}.fasta"
    if not r1.is_file():
        raise RuntimeError(f"FASTA missing in {bdir} for {accession}")
    return BatchPaths(r1=r1, r2=None, batch_dir=bdir)


def download_full(
    accession: str,
    paired: bool,
    outdir: Path,
    *,
    threads: int = 4,
    fasterq_dump: str = "fasterq-dump",
    tmpdir: Path | None = None,
) -> BatchPaths:
    """Download an entire SRA run as FASTA via fasterq-dump (multithreaded).

    Used by ``--bootstrap-all`` and for the legacy ``createDice`` workflow.
    fasterq-dump produces FASTQ; we convert via ``--fasta-unsorted``.
    """
    _require(fasterq_dump)
    bdir = outdir / "full" / accession
    bdir.mkdir(parents=True, exist_ok=True)
    tmp = tmpdir or (bdir / "tmp")
    tmp.mkdir(parents=True, exist_ok=True)

    cmd = [
        fasterq_dump,
        "--threads", str(threads),
        "--fasta-unsorted",
        "--skip-technical",
        "-O", str(bdir),
        "-t", str(tmp),
    ]
    if paired:
        cmd.append("--split-files")
    cmd.append(accession)

    log.info("fasterq-dump full run %s threads=%d", accession, threads)
    subprocess.run(cmd, check=True)

    if paired:
        r1 = bdir / f"{accession}_1.fasta"
        r2 = bdir / f"{accession}_2.fasta"
        if not r2.is_file():  # 3-file fallback as above
            fastas = sorted(bdir.glob("*.fasta"))
            if len(fastas) >= 2:
                r1 = fastas[0]
                r2 = fastas[-1]
        return BatchPaths(r1=r1, r2=r2, batch_dir=bdir)

    return BatchPaths(r1=bdir / f"{accession}.fasta", r2=None, batch_dir=bdir)
