"""Tests for varus.estimator.AdvancedEstimator."""

from __future__ import annotations

import math

import pytest

from varus.estimator import AdvancedEstimator


TILES = [("chr1", 0), ("chr1", 1), ("chr1", 2)]


def test_empty_tiles_returns_empty_dicts():
    est = AdvancedEstimator()
    result = est.estimate(
        tiles=[],
        obs_total={},
        run_obs=[{}, {}],
        times_downloaded=[0, 1],
    )
    assert result == [{}, {}]


def test_undownloaded_run_gets_prior():
    """Runs with timesDownloaded==0 share the same prior distribution."""
    obs_total = {("chr1", 0): 10, ("chr1", 1): 0, ("chr1", 2): 30}
    # run 0 has some observations; run 1 has none (not downloaded yet)
    run_obs = [
        {("chr1", 0): 5, ("chr1", 2): 5},
        {},
    ]
    est = AdvancedEstimator(lambda_=1.0, pseudo_count=1.0)
    result = est.estimate(
        tiles=TILES,
        obs_total=obs_total,
        run_obs=run_obs,
        times_downloaded=[2, 0],
    )
    # Both results must be dicts keyed on TILES
    assert set(result[0].keys()) == set(TILES)
    assert set(result[1].keys()) == set(TILES)

    # Probabilities must sum to 1 (within floating-point tolerance)
    assert abs(sum(result[0].values()) - 1.0) < 1e-9
    assert abs(sum(result[1].values()) - 1.0) < 1e-9


def test_two_undownloaded_runs_share_prior():
    """All undownloaded runs should receive identical distributions."""
    obs_total = {("chr1", 0): 100, ("chr1", 1): 50}
    run_obs = [{}, {}]
    est = AdvancedEstimator()
    result = est.estimate(
        tiles=[("chr1", 0), ("chr1", 1)],
        obs_total=obs_total,
        run_obs=run_obs,
        times_downloaded=[0, 0],
    )
    # Both dicts should be the same object or have the same values
    assert result[0] == result[1]


def test_downloaded_run_higher_mass_on_observed_tile():
    """A run whose observations are concentrated on tile 0 should have
    higher p[tile 0] than the prior."""
    tiles = [("chr1", 0), ("chr1", 1)]
    obs_total = {("chr1", 0): 100, ("chr1", 1): 100}
    run_obs = [
        {("chr1", 0): 90, ("chr1", 1): 10},  # run 0: heavily on tile 0
        {},                                    # run 1: not downloaded
    ]
    est = AdvancedEstimator(lambda_=10.0, pseudo_count=1.0)
    result = est.estimate(
        tiles=tiles,
        obs_total=obs_total,
        run_obs=run_obs,
        times_downloaded=[1, 0],
    )
    # Downloaded run should put more mass on tile 0
    assert result[0][("chr1", 0)] > result[0][("chr1", 1)]
    # Prior (undownloaded run) should be more symmetric but still biased by p_total
    prior = result[1]
    assert abs(sum(prior.values()) - 1.0) < 1e-9


def test_all_zeros_obs_total_gives_uniform_prior():
    """With no observations, prior should be uniform (all tiles equal weight)."""
    tiles = [("chr1", 0), ("chr1", 1), ("chr1", 2)]
    obs_total = {}  # no observations at all
    est = AdvancedEstimator(lambda_=10.0, pseudo_count=1.0)
    result = est.estimate(
        tiles=tiles,
        obs_total=obs_total,
        run_obs=[{}],
        times_downloaded=[0],
    )
    probs = list(result[0].values())
    assert abs(probs[0] - probs[1]) < 1e-9
    assert abs(probs[1] - probs[2]) < 1e-9
    assert abs(sum(probs) - 1.0) < 1e-9


def test_lambda_zero_gives_pseudo_count_uniform():
    """With lambda=0 and uniform pseudo-count, prior is uniform regardless of p_total."""
    tiles = [("chr1", 0), ("chr1", 1)]
    obs_total = {("chr1", 0): 999, ("chr1", 1): 1}  # very skewed
    est = AdvancedEstimator(lambda_=0.0, pseudo_count=1.0)
    result = est.estimate(
        tiles=tiles,
        obs_total=obs_total,
        run_obs=[{}],
        times_downloaded=[0],
    )
    probs = list(result[0].values())
    assert abs(probs[0] - 0.5) < 1e-9
