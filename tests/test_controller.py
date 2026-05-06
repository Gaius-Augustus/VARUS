"""Tests for varus.controller: RunState, VARUSConfig, Controller logic.

All tests use synthetic data and mock download/align so no actual SRA
downloads or alignments are made.
"""

from __future__ import annotations

import math
import random
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from varus.controller import Controller, RunState, VARUSConfig, load_runs
from varus.runlist import RunRecord


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_record(acc: str = "SRR1", spots: int = 200_000,
                 paired: bool = False) -> RunRecord:
    return RunRecord(
        accession=acc,
        total_spots=spots,
        total_bases=spots * 100,
        avg_len=100.0,
        paired=paired,
        colorspace=False,
    )


def _make_config(tmp_path: Path, **kwargs) -> VARUSConfig:
    defaults = dict(
        genome=tmp_path / "genome.fa",
        index_prefix=tmp_path / "idx/hisatidx",
        outdir=tmp_path / "out",
        batch_size=50_000,
        max_batches=5,
        tile_size=5_000,
        cost=0.0,
        profit_condition=False,  # keep looping for tests unless max_batches hit
    )
    defaults.update(kwargs)
    return VARUSConfig(**defaults)


# ---------------------------------------------------------------------------
# RunState
# ---------------------------------------------------------------------------

def test_runstate_sigma_length():
    rec = _make_record(spots=130_000)
    rng = random.Random(0)
    rs = RunState.from_record(rec, batch_size=50_000, rng=rng)
    assert rs.max_batches == 3        # ceil(130000/50000)
    assert len(rs.sigma) == 3
    assert sorted(rs.sigma) == [0, 1, 2]


def test_runstate_sigma_tail_fixed():
    """Last element in sigma must always be max_batches-1 (shuffleExceptLast)."""
    rng = random.Random(42)
    rec = _make_record(spots=300_000)
    rs = RunState.from_record(rec, batch_size=50_000, rng=rng)
    assert rs.sigma[-1] == rs.max_batches - 1


def test_runstate_next_batch_range():
    rec = _make_record(spots=120_000)
    rng = random.Random(0)
    rs = RunState.from_record(rec, batch_size=50_000, rng=rng)
    # Force sigma[0] = 0 for predictability
    rs.sigma = [0, 1, 2]
    n, x = rs.next_batch_range(batch_size=50_000)
    assert n == 0
    assert x == 49_999


def test_runstate_last_batch_capped():
    """Last batch x should not exceed total_spots-1."""
    rec = _make_record(spots=120_000)
    rng = random.Random(0)
    rs = RunState.from_record(rec, batch_size=50_000, rng=rng)
    rs.sigma = [2]  # batch index 2
    rs.sigma_idx = 0
    n, x = rs.next_batch_range(batch_size=50_000)
    assert n == 100_000
    assert x == 119_999   # capped at total_spots - 1


def test_runstate_is_exhausted():
    rec = _make_record(spots=50_000)
    rng = random.Random(0)
    rs = RunState.from_record(rec, batch_size=50_000, rng=rng)
    assert not rs.is_exhausted
    rs.sigma_idx = rs.max_batches
    assert rs.is_exhausted


# ---------------------------------------------------------------------------
# Controller: pure logic (no I/O)
# ---------------------------------------------------------------------------

def test_controller_score_empty():
    cfg = _make_config(Path("/tmp"), profit_condition=False)
    ctrl = Controller(cfg, [])
    assert ctrl._score() == 0.0


def test_controller_score_nonempty():
    cfg = _make_config(Path("/tmp"), profit_condition=False)
    ctrl = Controller(cfg, [])
    ctrl.total_obs = {("chr1", 0): 1, ("chr1", 1): 2}
    expected = math.log1p(1) + math.log1p(2)
    assert abs(ctrl._score() - expected) < 1e-12


def test_controller_continuing_max_batches(tmp_path: Path):
    cfg = _make_config(tmp_path, max_batches=3, profit_condition=False)
    ctrl = Controller(cfg, [])
    ctrl.batch_count = 2
    assert ctrl._continuing()
    ctrl.batch_count = 3
    assert not ctrl._continuing()


def test_controller_continuing_profit_condition(tmp_path: Path):
    cfg = _make_config(tmp_path, max_batches=0, profit_condition=True)
    ctrl = Controller(cfg, [])
    ctrl.max_profit = 0.5
    assert ctrl._continuing()
    ctrl.max_profit = -0.1
    assert not ctrl._continuing()
    ctrl.max_profit = 0.0
    assert not ctrl._continuing()


def test_controller_choose_next_run_max_profit(tmp_path: Path):
    """Controller should select the run with highest expectedProfit."""
    rng = random.Random(1)
    rec_a = _make_record("SRR_A")
    rec_b = _make_record("SRR_B")
    cfg = _make_config(tmp_path)
    rs_a = RunState.from_record(rec_a, cfg.batch_size, rng)
    rs_b = RunState.from_record(rec_b, cfg.batch_size, rng)
    rs_a.expected_profit = 0.5
    rs_b.expected_profit = 1.2
    ctrl = Controller(cfg, [rs_a, rs_b])
    chosen = ctrl._choose_next_run()
    assert chosen is rs_b


def test_controller_profit_zero_obs(tmp_path: Path):
    """With empty total_obs and no downloads, profit should be ~0 (cost=0)."""
    rng = random.Random(0)
    rec = _make_record()
    cfg = _make_config(tmp_path, cost=0.0)
    rs = RunState.from_record(rec, cfg.batch_size, rng)
    rs.p = {}  # empty p → profit = 0 - 0 = 0
    ctrl = Controller(cfg, [rs])
    assert ctrl._profit(rs) == 0.0


def test_controller_update_downloadable(tmp_path: Path):
    rng = random.Random(0)
    cfg = _make_config(tmp_path)
    rs_ok = RunState.from_record(_make_record("SRR_OK"), cfg.batch_size, rng)
    rs_bad = RunState.from_record(_make_record("SRR_BAD"), cfg.batch_size, rng)
    rs_bad.bad_quality = True
    ctrl = Controller(cfg, [rs_ok, rs_bad])
    ctrl._update_downloadable()
    assert rs_ok in ctrl.downloadable
    assert rs_bad not in ctrl.downloadable


# ---------------------------------------------------------------------------
# load_runs
# ---------------------------------------------------------------------------

def test_load_runs_parses_runlist(tmp_path: Path):
    runlist = tmp_path / "Runlist.tsv"
    runlist.write_text(
        "@Run_acc\ttotal_spots\ttotal_bases\tavg_len\tbool:paired\tcolor_space\n"
        "SRR1\t100000\t10000000\t100.0\t0\t0\n"
        "SRR2\t200000\t20000000\t150.0\t1\t0\n"
        "SRR3\t50000\t5000000\t100.0\t0\t1\n"  # colorspace -> skip
    )
    rng = random.Random(0)
    runs = load_runs(runlist, batch_size=50_000, rng=rng)
    assert len(runs) == 2
    assert runs[0].record.accession == "SRR1"
    assert runs[1].record.accession == "SRR2"


def test_load_runs_skips_colorspace(tmp_path: Path):
    runlist = tmp_path / "Runlist.tsv"
    runlist.write_text(
        "@Run_acc\ttotal_spots\ttotal_bases\tavg_len\tbool:paired\tcolor_space\n"
        "SRR_CS\t50000\t5000000\t100.0\t0\t1\n"
    )
    runs = load_runs(runlist, batch_size=50_000, rng=random.Random(0))
    assert runs == []
