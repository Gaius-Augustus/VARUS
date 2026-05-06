"""Shared test fixtures and skip markers.

``requires_pysam`` skips a test if pysam is not importable. pysam needs htslib,
which doesn't build on stock Windows. Tests that build a BAM are gated on this
so the local Windows test pass stays green and the BAM tests still run on the
HPC/Linux CI.
"""

from __future__ import annotations

import importlib.util

import pytest


def _has(mod: str) -> bool:
    return importlib.util.find_spec(mod) is not None


requires_pysam = pytest.mark.skipif(
    not _has("pysam"),
    reason="pysam not installed (extras 'align' on Linux/macOS)",
)
