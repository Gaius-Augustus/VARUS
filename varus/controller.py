"""VARUS online run-selection controller.

Implements the algorithm from Stanke et al. (2019) BMC Bioinformatics:
greedy maximization of the expected score gain

    S(c) = Σ_j ln(1 + c_j)

over 5 kb genome tiles j, where c_j is the UMR count in tile j. Each
iteration picks the SRA run whose next batch is expected to increase S the
most, downloads it, aligns it, and updates tile counts.

Key classes
-----------
VARUSConfig : All tunable parameters in one place.
RunState    : Per-run mutable state (observations, sigma, stats).
Controller  : Runs the online loop; owns global state (total_obs, introns).
"""

from __future__ import annotations

import logging
import math
import random
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from varus.align import align_batch_hisat2, parse_hisat2_log
from varus.download import download_batch
from varus.estimator import AdvancedEstimator
from varus.introns import IntronCounts, extract_introns_from_bam, write_introns_gff
from varus.io import write_coverage, write_run_statistics
from varus.merge import merge_bams
from varus.runlist import RunRecord
from varus.strand import assign_strand, write_hisat2_splice_sites
from varus.tiles import BAMStats, count_bam_stats

log = logging.getLogger(__name__)

Tile = Tuple[str, int]


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

@dataclass
class VARUSConfig:
    """All controller parameters in one place."""
    genome: Path
    index_prefix: Path        # hisat2 index prefix (e.g. genome/hisatidx)
    outdir: Path

    batch_size: int = 50_000
    max_batches: int = 1_000
    tile_size: int = 5_000
    min_uniq_pct: float = 5.0  # quality gate; bad-quality if below

    threads: int = 4
    keep_batches: bool = False
    coverage_trace: int = 0    # write Coverage<N>.tsv every N batches; 0=never
    seed: Optional[int] = None
    bootstrap_all: bool = False

    # Estimator hyperparameters
    lambda_: float = 10.0
    pseudo_count: float = 1.0
    cost: float = 0.0          # per-read download cost (default 0 = ignore cost)
    # Stop early when expected profit ≤ 0. Off by default: matches the legacy
    # production pipeline (--profitCondition 0). When on, the check is also
    # skipped while no observations have been collected yet (cold start),
    # so the algorithm always gets at least one batch to bootstrap.
    profit_condition: bool = False


# ---------------------------------------------------------------------------
# Per-run state
# ---------------------------------------------------------------------------

@dataclass
class RunState:
    """Mutable state for a single SRA run."""
    record: RunRecord
    max_batches: int
    sigma: List[int]           # shuffled batch indices

    sigma_idx: int = 0
    times_downloaded: int = 0
    observations: Dict[Tile, int] = field(default_factory=dict)
    p: Dict[Tile, float] = field(default_factory=dict)
    expected_profit: float = 0.0
    avg_umr_pct: float = 0.0
    avg_spliced_pct: float = 0.0
    bad_quality: bool = False

    @classmethod
    def from_record(
        cls, record: RunRecord, batch_size: int, rng: random.Random
    ) -> "RunState":
        """Build a RunState from a RunRecord, initialising the sigma vector."""
        max_batches = max(1, math.ceil(record.total_spots / batch_size))
        sigma = list(range(max_batches))
        if len(sigma) > 1:
            # Shuffle all but the last element: mirrors shuffleExceptLast()
            # in legacy ChromosomeInitializer.cpp so the potentially-short
            # final batch is always downloaded last.
            front, tail = sigma[:-1], sigma[-1]
            rng.shuffle(front)
            sigma = front + [tail]
        return cls(record=record, max_batches=max_batches, sigma=sigma)

    @property
    def is_exhausted(self) -> bool:
        return self.sigma_idx >= len(self.sigma)

    def next_batch_range(self, batch_size: int) -> Tuple[int, int]:
        """Return (n, x) spot range for the next batch index in sigma."""
        k = self.sigma[self.sigma_idx]
        n = k * batch_size
        x = min((k + 1) * batch_size - 1, self.record.total_spots - 1)
        return n, x


# ---------------------------------------------------------------------------
# Controller
# ---------------------------------------------------------------------------

class Controller:
    """Online run-selection loop.

    Usage::

        cfg = VARUSConfig(genome=..., index_prefix=..., outdir=...)
        runs = [RunState.from_record(r, cfg.batch_size, rng) for r in records]
        ctrl = Controller(cfg, runs)
        ctrl.run()
    """

    def __init__(self, config: VARUSConfig, runs: List[RunState]) -> None:
        self.config = config
        self.runs = runs
        self.downloadable: List[RunState] = list(runs)
        self.rng = random.Random(config.seed)

        self.total_obs: Dict[Tile, int] = {}
        self.cumulative_introns: IntronCounts = IntronCounts({})
        self._batch_bams: List[Path] = []

        self.batch_count = 0
        self.total_score = 0.0
        self.total_profit = 0.0
        self.max_profit: float = 1.0  # initialised >0 so continuing() starts True

        # Running average of UMR% and spliced% across all downloaded batches.
        # Priors are 100% to overestimate initially (encourages exploration).
        self._n_stat_runs = 4        # pseudocount (matches legacy numRuns=4)
        self._sum_umr = 100.0 * self._n_stat_runs
        self._sum_spliced = 100.0 * self._n_stat_runs

        self.avg_uniq: float = 100.0
        self.avg_spliced: float = 100.0

        self.estimator = AdvancedEstimator(
            lambda_=config.lambda_, pseudo_count=config.pseudo_count
        )

        self.config.outdir.mkdir(parents=True, exist_ok=True)
        self._splice_db_path = config.outdir / "intronDB.splice_sites"

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

    def run(self) -> None:
        """Execute the online sampling loop."""
        if not self.downloadable:
            log.warning("No downloadable runs; nothing to do.")
            return

        if self.config.bootstrap_all:
            self._bootstrap()

        self._update_downloadable()
        self._estimate_p()
        self._calculate_profit()

        while self._continuing():
            self.batch_count += 1
            log.info(
                "Iteration %d | downloadable=%d | maxProfit=%.4f",
                self.batch_count, len(self.downloadable), self.max_profit,
            )

            run = self._choose_next_run()
            if run is None:
                break

            success = self._process_batch(run)
            if not success:
                run.bad_quality = True

            self.total_score = self._score()
            self.total_profit = (
                self.total_score
                - self.config.cost * self.config.batch_size * self.batch_count
            )

            self._estimate_p()
            self._calculate_profit()
            self._export_stats(run)
            self._update_downloadable()

            if not self.downloadable:
                log.info("No downloadable runs left.")
                break

        self._finalize()

    # ------------------------------------------------------------------
    # Algorithm steps
    # ------------------------------------------------------------------

    def _bootstrap(self) -> None:
        """Download and process one batch from every run (--bootstrap-all)."""
        log.info("Bootstrap: processing one batch from each of %d runs", len(self.runs))
        for run in list(self.runs):
            if run.bad_quality or run.is_exhausted:
                continue
            ok = self._process_batch(run)
            if not ok:
                run.bad_quality = True
        self._estimate_p()

    def _estimate_p(self) -> None:
        """Re-estimate tile probability distributions for all runs."""
        tiles = sorted(self.total_obs.keys())
        if not tiles:
            return
        p_dists = self.estimator.estimate(
            tiles=tiles,
            obs_total=self.total_obs,
            run_obs=[r.observations for r in self.runs],
            times_downloaded=[r.times_downloaded for r in self.runs],
        )
        for run, p in zip(self.runs, p_dists):
            run.p = p

    def _calculate_profit(self) -> None:
        """Compute expectedProfit for every downloadable run; update avg stats."""
        n_stat = 4              # pseudocount (legacy numRuns=4)
        sum_umr = 100.0 * n_stat
        sum_spliced = 100.0 * n_stat

        no_reads_profit: Optional[float] = None

        for run in self.downloadable:
            if run.times_downloaded == 0 and no_reads_profit is not None:
                run.expected_profit = no_reads_profit
            else:
                run.expected_profit = self._profit(run)
                if run.times_downloaded == 0:
                    no_reads_profit = run.expected_profit
                    log.debug(
                        "Prior profit (undownloaded runs): %.4f", no_reads_profit
                    )

            if run.times_downloaded > 0:
                sum_umr += run.avg_umr_pct
                sum_spliced += run.avg_spliced_pct
                n_stat += 1

        self.avg_uniq = sum_umr / n_stat
        self.avg_spliced = sum_spliced / n_stat

        if self.downloadable:
            self.max_profit = max(r.expected_profit for r in self.downloadable)

    def _profit(self, run: RunState) -> float:
        """Expected score gain from downloading one more batch of run r.

        Matches Controller::profit() in the legacy code.
        """
        umr_pct = run.avg_umr_pct if run.times_downloaded > 0 else self.avg_uniq
        spliced_pct = (
            run.avg_spliced_pct if run.times_downloaded > 0 else self.avg_spliced
        )
        effective = (umr_pct + spliced_pct) / 100.0 * self.config.batch_size

        p_dist = run.p
        pr = 0.0
        for tile, prob in p_dist.items():
            x = self.total_obs.get(tile, 0)
            pr += math.log1p(x + prob * effective) - math.log1p(x)
        pr -= self.config.cost * self.config.batch_size
        return pr

    def _choose_next_run(self) -> Optional[RunState]:
        """Return the run with highest expectedProfit; break ties by avg_len.

        Ties are resolved by a weighted random draw proportional to avg_len
        (matching the legacy biasSelect() which prefers longer reads).
        """
        if not self.downloadable:
            return None

        best_profit = max(r.expected_profit for r in self.downloadable)
        self.max_profit = best_profit

        candidates = [
            r for r in self.downloadable if r.expected_profit == best_profit
        ]
        if len(candidates) == 1:
            return candidates[0]

        weights = [r.record.avg_len for r in candidates]
        total_w = sum(weights) or 1.0
        pick = self.rng.random() * total_w
        cumulative = 0.0
        for run, w in zip(candidates, weights):
            cumulative += w
            if pick <= cumulative:
                return run
        return candidates[-1]

    def _continuing(self) -> bool:
        """True while the loop should keep running."""
        if self.config.max_batches > 0 and self.batch_count >= self.config.max_batches:
            log.info("Reached max_batches=%d; stopping.", self.config.max_batches)
            return False
        # Skip the profit check until we actually have observations. Without
        # this, the algorithm cannot bootstrap: every run starts with an empty
        # p distribution, profits are 0, and the loop would stop on iteration 1.
        if (
            self.config.profit_condition
            and self.total_obs
            and self.max_profit <= 0
        ):
            log.info("maxProfit=%.4f ≤ 0; stopping.", self.max_profit)
            return False
        return True

    def _score(self) -> float:
        """S(c) = Σ ln(1 + c_j) over all tiles with observations."""
        return sum(math.log1p(v) for v in self.total_obs.values())

    def _process_batch(self, run: RunState) -> bool:
        """Download, align, and count one batch from run. Return True on success."""
        n, x = run.next_batch_range(self.config.batch_size)
        run.sigma_idx += 1

        # 1. Download
        try:
            paths = download_batch(
                accession=run.record.accession,
                n=n, x=x,
                paired=run.record.paired,
                outdir=self.config.outdir,
            )
        except RuntimeError as e:
            log.warning(
                "Download failed for %s [N=%d X=%d]: %s",
                run.record.accession, n, x, e,
            )
            return False

        # 2. Align (use splice DB if available from a previous batch)
        intron_db = (
            self._splice_db_path
            if self._splice_db_path.is_file()
            else None
        )
        try:
            result = align_batch_hisat2(
                r1=paths.r1,
                r2=paths.r2,
                index_prefix=self.config.index_prefix,
                batch_dir=paths.batch_dir,
                threads=self.config.threads,
                intron_db=intron_db,
            )
        except RuntimeError as e:
            log.warning("Alignment failed for %s: %s", run.record.accession, e)
            if not self.config.keep_batches:
                for p in paths.as_list():
                    p.unlink(missing_ok=True)
            return False

        # 3. Quality gate: reject if too few uniquely mapped reads
        stats = parse_hisat2_log(result.log, batch_size=self.config.batch_size)
        uniq_pct = stats["uniq_pct"]
        if uniq_pct < self.config.min_uniq_pct:
            log.warning(
                "Run %s batch [%d-%d]: uniq_pct=%.1f%% < min=%.1f%%; bad quality",
                run.record.accession, n, x, uniq_pct, self.config.min_uniq_pct,
            )
            if not self.config.keep_batches:
                for p in paths.as_list():
                    p.unlink(missing_ok=True)
            return False

        # 4. Count UMRs and spliced reads
        bam_stats = count_bam_stats(result.bam, self.config.tile_size)

        # 5. Update run and global observations
        for tile, count in bam_stats.umr_counts.items():
            run.observations[tile] = run.observations.get(tile, 0) + count
            self.total_obs[tile] = self.total_obs.get(tile, 0) + count

        spliced_pct = (
            100.0 * bam_stats.n_spliced / bam_stats.n_reads
            if bam_stats.n_reads > 0
            else 0.0
        )
        run.times_downloaded += 1
        nd = run.times_downloaded
        run.avg_umr_pct += (uniq_pct - run.avg_umr_pct) / nd
        run.avg_spliced_pct += (spliced_pct - run.avg_spliced_pct) / nd

        log.info(
            "Run %s batch %d: uniq=%.1f%% spliced=%.1f%% UMRs=%d",
            run.record.accession, nd, uniq_pct, spliced_pct,
            sum(bam_stats.umr_counts.values()),
        )

        # 6. Accumulate introns and rebuild splice-site DB
        batch_introns = extract_introns_from_bam(result.bam)
        self.cumulative_introns = self.cumulative_introns.merge(batch_introns)
        self._rebuild_intron_db()

        # 7. Track BAM for final merge
        self._batch_bams.append(result.bam)

        # 8. Clean up FASTA (keep BAM for final merge)
        if not self.config.keep_batches:
            for p in paths.as_list():
                p.unlink(missing_ok=True)

        return True

    def _rebuild_intron_db(self) -> None:
        """Assign strand to cumulative introns and write HISAT2 splice-site file."""
        if not self.cumulative_introns.counts:
            return
        try:
            stranded = assign_strand(
                self.cumulative_introns, self.config.genome
            )
            write_hisat2_splice_sites(stranded, self._splice_db_path)
        except Exception as e:
            log.warning("Intron DB rebuild failed: %s", e)

    def _update_downloadable(self) -> None:
        """Remove exhausted and bad-quality runs from the downloadable list."""
        before = len(self.downloadable)
        self.downloadable = [
            r for r in self.downloadable
            if not r.bad_quality and not r.is_exhausted
        ]
        removed = before - len(self.downloadable)
        if removed:
            log.info(
                "Removed %d run(s); %d remain downloadable.",
                removed, len(self.downloadable),
            )

    # ------------------------------------------------------------------
    # Output helpers
    # ------------------------------------------------------------------

    def _export_stats(self, last_run: RunState) -> None:
        """Write per-batch coverage trace and run statistics."""
        cov_path = self.config.outdir / "Coverage.csv"
        stats_path = self.config.outdir / "RunStatistics.csv"

        # Coverage trace: every coverage_trace batches, write a snapshot
        if (
            self.config.coverage_trace > 0
            and self.batch_count % self.config.coverage_trace == 0
        ):
            trace = self.config.outdir / f"Coverage{self.batch_count}.tsv"
            write_coverage(self.total_obs, trace)

        # RunStatistics: every batch for the first 10, then every 10
        if self.batch_count < 10 or self.batch_count % 10 == 1:
            write_run_statistics(self.runs, stats_path)

    def _finalize(self) -> None:
        """Merge all batch BAMs, write final Coverage.csv and RunStatistics.csv."""
        log.info("Finalizing: merging %d batch BAMs", len(self._batch_bams))

        write_coverage(self.total_obs, self.config.outdir / "Coverage.csv")
        write_run_statistics(self.runs, self.config.outdir / "RunStatistics.csv")

        if self._batch_bams:
            out_bam = self.config.outdir / "VARUS.bam"
            try:
                merge_bams(
                    self._batch_bams,
                    out_bam,
                    threads=self.config.threads,
                )
                log.info("Final BAM: %s", out_bam)
            except (RuntimeError, ValueError) as e:
                log.error("Final merge failed: %s", e)

            # Delete per-batch BAMs unless the user asked to keep them
            if not self.config.keep_batches:
                for bam in self._batch_bams:
                    bam.unlink(missing_ok=True)

        # Write the cumulative intron GFF alongside the final BAM
        if self.cumulative_introns.counts:
            gff_path = self.config.outdir / "introns.gff"
            write_introns_gff(self.cumulative_introns, gff_path)
            log.info("Cumulative introns: %s", gff_path)


# ---------------------------------------------------------------------------
# Runlist loader (called from CLI)
# ---------------------------------------------------------------------------

def load_runs(
    runlist_path: Path,
    batch_size: int,
    rng: random.Random,
    paired_only: bool = False,
) -> List[RunState]:
    """Parse a Runlist.tsv and return a list of RunState objects."""
    from varus.runlist import RunRecord

    records: List[RunRecord] = []
    with runlist_path.open(encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("@"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            acc = parts[0]
            try:
                spots = int(parts[1])
                bases = int(parts[2])
                avg_len = float(parts[3])
                paired = bool(int(parts[4]))
                colorspace = bool(int(parts[5]))
            except (ValueError, IndexError):
                log.warning("Skipping malformed runlist line: %s", line.rstrip())
                continue
            if colorspace:
                continue
            if paired_only and not paired:
                continue
            records.append(
                RunRecord(
                    accession=acc,
                    total_spots=spots,
                    total_bases=bases,
                    avg_len=avg_len,
                    paired=paired,
                    colorspace=colorspace,
                )
            )

    log.info("Loaded %d runs from %s", len(records), runlist_path)
    return [RunState.from_record(r, batch_size, rng) for r in records]
