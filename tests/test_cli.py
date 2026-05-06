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
    assert args.batch_size == 50_000
    assert args.tile_size == 5_000
    assert args.max_batches == 1_000
    assert args.threads == 4
    assert args.keep_batches is False
    assert args.bootstrap_all is False


def test_subcommand_required():
    with pytest.raises(SystemExit):
        cli.build_parser().parse_args([])
