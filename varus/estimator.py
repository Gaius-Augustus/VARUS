"""AdvancedEstimator: tile probability distribution from paper eq. 3.

Stanke et al. (2019) BMC Bioinformatics, eq. 3:

    p̂_r[j] ∝ c^r_j + λ·T·p̄[j] + a

where:
  c^r_j = UMR count from run r in tile j
  p̄[j]  = c_total[j] / Σ c_total  (pooled UMR fraction)
  T      = number of tiles with at least one total observation
  λ      = smoothing coefficient (default 10.0)
  a      = pseudo-count (default 1.0)

Runs not yet downloaded share a common prior that is computed from the pooled
observations of all downloaded runs.  This matches the ``pRep`` optimization
in ``AdvancedEstimator.cpp``: only one p-vector is computed for undownloaded
runs and the rest point to it.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np

Tile = Tuple[str, int]


class AdvancedEstimator:
    """Compute per-run tile probability distributions."""

    def __init__(self, lambda_: float = 10.0, pseudo_count: float = 1.0) -> None:
        self.lambda_ = lambda_
        self.pseudo_count = pseudo_count

    def estimate(
        self,
        tiles: List[Tile],
        obs_total: Dict[Tile, int],
        run_obs: List[Dict[Tile, int]],
        times_downloaded: List[int],
    ) -> List[Dict[Tile, float]]:
        """Return one {tile: probability} dict per run.

        Parameters
        ----------
        tiles:            Ordered list of all tiles with nonzero pooled count.
        obs_total:        Pooled UMR counts across every run.
        run_obs:          Per-run {tile: count} observations.
        times_downloaded: Download count per run (parallel to run_obs).
        """
        if not tiles:
            return [{} for _ in run_obs]

        T = len(tiles)
        total_arr = np.array(
            [obs_total.get(t, 0) for t in tiles], dtype=np.float64
        )
        total_sum = total_arr.sum()
        p_total = total_arr / total_sum if total_sum > 0 else np.full(T, 1.0 / T)

        # Prior shared by all runs with zero downloads
        raw_prior = self.pseudo_count + self.lambda_ * p_total * T
        prior = raw_prior / raw_prior.sum()
        prior_dict: Dict[Tile, float] = dict(zip(tiles, prior.tolist()))

        results: List[Dict[Tile, float]] = []
        for obs, nd in zip(run_obs, times_downloaded):
            if nd == 0:
                results.append(prior_dict)
                continue
            c = np.array([obs.get(t, 0) for t in tiles], dtype=np.float64)
            raw = c + self.pseudo_count + self.lambda_ * p_total * T
            p = raw / raw.sum()
            results.append(dict(zip(tiles, p.tolist())))
        return results
