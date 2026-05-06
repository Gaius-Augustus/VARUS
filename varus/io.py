"""Write Coverage.csv and RunStatistics.csv output files.

Format matches the legacy ``exportCoverage()`` and ``exportRunStatistics()``
from ``Controller.cpp`` so existing downstream scripts keep working.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Tuple

Tile = Tuple[str, int]


def write_coverage(
    total_obs: Dict[Tile, int],
    path: Path,
) -> None:
    """Write ``Coverage.csv`` (overwritten on each call).

    Columns: ``blockName;totalObservations``

    The blockName follows the legacy ``chrom:tile_idx`` convention so the file
    can be parsed by existing evaluation scripts.
    """
    with path.open("w", encoding="utf-8") as f:
        f.write("blockName;totalObservations\n")
        for (chrom, idx), count in sorted(total_obs.items()):
            f.write(f"{chrom}:{idx};{count}\n")


def write_run_statistics(runs: Iterable, path: Path) -> None:
    """Write ``RunStatistics.csv`` (overwritten on each call).

    Sorted by timesDownloaded descending, then avgUmrPercent descending —
    matching the ``runCount`` comparator in the legacy code.

    Each ``run`` object must expose:
    - ``record.accession``
    - ``times_downloaded``
    - ``avg_umr_pct``
    - ``bad_quality``
    - ``max_batches``
    """
    rows = sorted(
        runs,
        key=lambda r: (-r.times_downloaded, -r.avg_umr_pct),
    )
    with path.open("w", encoding="utf-8") as f:
        f.write("accesion-id;timesDownloaded;avgUmrPercent;badQuality;maxNumOfBatches\n")
        for r in rows:
            f.write(
                f"{r.record.accession};{r.times_downloaded};"
                f"{r.avg_umr_pct:.4f}%;{r.bad_quality};"
                f"{r.max_batches}\n"
            )
