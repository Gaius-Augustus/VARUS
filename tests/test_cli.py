"""Smoke tests for the argparse interface."""

from __future__ import annotations

import pytest

from varus import cli


def test_parser_builds():
    p = cli.build_parser()
    assert p.prog == "varus"


def test_runlist_parser_accepts_minimal_args():
    args = cli.build_parser().parse_args(
        ["runlist", "Foo bar"]
    )
    assert args.cmd == "runlist"
    assert args.species == "Foo bar"
    assert args.paired_only is False
    assert args.max_runs == 0


def test_index_parser_accepts_minimal_args():
    args = cli.build_parser().parse_args(
        ["index", "genome.fa"]
    )
    assert args.cmd == "index"
    assert str(args.genome) == "genome.fa"
    assert args.threads == 4


def test_run_subcommand_parser_defaults(tmp_path):
    """'varus run' parser is wired up and applies correct defaults."""
    args = cli.build_parser().parse_args(
        [
            "run", "Foo bar", "genome.fa",
            "--runlist", str(tmp_path / "Runlist.tsv"),
            "--index", str(tmp_path / "genome/"),
        ]
    )
    assert args.cmd == "run"
    # Batch size and min-mapq default to None so main() can pick a mode-aware value.
    assert args.batch_size is None
    assert args.min_mapq is None
    assert args.tile_size == 5_000
    assert args.max_batches == 1_000
    assert args.threads == 4
    assert args.keep_batches is False
    assert args.bootstrap_all is False
    assert args.longreads is False


def test_index_parser_longreads_flag(tmp_path):
    """--longreads on the index subcommand toggles minimap2 mode."""
    args = cli.build_parser().parse_args(
        ["index", "genome.fa", "--longreads"]
    )
    assert args.longreads is True
    # Prefix default is None so dispatcher can pick 'mm2idx' / 'hisatidx'.
    assert args.prefix is None


def test_runlist_parser_longreads_flag():
    args = cli.build_parser().parse_args(
        ["runlist", "Foo bar", "--longreads"]
    )
    assert args.longreads is True


def test_run_parser_longreads(tmp_path):
    args = cli.build_parser().parse_args(
        [
            "run", "Foo bar", "genome.fa",
            "--runlist", str(tmp_path / "Runlist.tsv"),
            "--index", str(tmp_path / "genome/mm2idx.mmi"),
            "--longreads",
        ]
    )
    assert args.longreads is True


def test_subcommand_required():
    with pytest.raises(SystemExit):
        cli.build_parser().parse_args([])
