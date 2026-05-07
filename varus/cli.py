"""Command-line interface for VARUS.

Subcommands
-----------
runlist : query NCBI SRA for all RNA-seq runs of a species, write Runlist.tsv
index   : build a HISAT2 (default) or minimap2 (``--longreads``) index
run     : execute the online sampling loop (download + align + score)
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path


def _add_runlist(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "runlist",
        help="Fetch SRA RNA-seq run list for a species via NCBI Entrez.",
    )
    p.add_argument("species", help="Binomial species name, e.g. 'Drosophila melanogaster'")
    p.add_argument("--outdir", type=Path, default=Path.cwd(),
                   help="Output directory (default: cwd). Writes Runlist.tsv.")
    p.add_argument("--max-runs", type=int, default=0,
                   help="Limit to first N runs (0 = all available).")
    p.add_argument("--paired-only", action="store_true",
                   help="Keep only paired-end runs.")
    p.add_argument("--longreads", action="store_true",
                   help="Restrict the SRA query to long-read platforms "
                        "(PacBio SMRT, Oxford Nanopore). Without this flag the "
                        "runlist is dominated by Illumina short-read runs.")
    p.add_argument("--email", default=None,
                   help="Contact email for NCBI Entrez (recommended; "
                        "falls back to $NCBI_EMAIL).")
    p.add_argument("--api-key", default=None,
                   help="NCBI API key (optional; falls back to $NCBI_API_KEY).")


def _add_index(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "index",
        help="Build a HISAT2 (default) or minimap2 (--longreads) index.",
    )
    p.add_argument("genome", type=Path, help="Genome FASTA file.")
    p.add_argument("--outdir", type=Path, default=Path("genome"),
                   help="Output directory for the index (default: ./genome/).")
    p.add_argument("--threads", type=int, default=4,
                   help="Threads for the index builder (default: 4).")
    p.add_argument("--prefix", default=None,
                   help="Index file prefix (default: 'hisatidx' for short reads, "
                        "'mm2idx' for --longreads).")
    p.add_argument("--longreads", action="store_true",
                   help="Build a minimap2 splice index instead of HISAT2.")


def _add_run(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "run",
        help="Run the online sampling loop (download + align + score).",
    )
    p.add_argument("species", help="Binomial species name.")
    p.add_argument("genome", type=Path, help="Genome FASTA file.")
    p.add_argument("--runlist", type=Path, required=True, help="Path to Runlist.tsv.")
    p.add_argument("--index", type=Path, required=True,
                   help="HISAT2 index prefix (short reads, e.g. Sp/genome/hisatidx) "
                        "or minimap2 .mmi file (--longreads, e.g. Sp/genome/mm2idx.mmi).")
    p.add_argument("--outdir", type=Path, default=Path.cwd(),
                   help="Output directory.")
    p.add_argument("--batch-size", type=int, default=None,
                   help="Spots per batch (default: 50000 for short reads, "
                        "2000 for --longreads).")
    p.add_argument("--max-batches", type=int, default=1000)
    p.add_argument("--tile-size", type=int, default=5000)
    p.add_argument("--min-uniq-pct", type=float, default=5.0)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--keep-batches", action="store_true",
                   help="Keep per-batch FASTA/BAM files (default: delete after counting).")
    p.add_argument("--coverage-trace", type=int, default=0,
                   help="Write Coverage<N>.tsv every N batches (0 = never; final Coverage.csv "
                        "is always written).")
    p.add_argument("--seed", type=int, default=None,
                   help="Random seed (default: random).")
    p.add_argument("--bootstrap-all", action="store_true",
                   help="Download one batch from every run before starting the online loop "
                        "(equivalent to legacy --loadAllOnce).")
    p.add_argument("--profit-condition", action="store_true",
                   help="Stop early when expected profit ≤ 0. Off by default; matches the "
                        "legacy production setting (--profitCondition 0). The check is "
                        "always skipped on cold start (before any observations).")
    p.add_argument("--pipeline-downloads", action="store_true",
                   help="Overlap round R+1's downloads (network-bound, single-threaded) "
                        "with round R's alignments (CPU-bound, multi-threaded). Adds one "
                        "extra round of staleness to picks; expect 1+T_dl/T_al speedup "
                        "(typically 1.3–1.8×).")
    p.add_argument("--longreads", action="store_true",
                   help="Align with minimap2 instead of HISAT2 (for PacBio Iso-Seq "
                        "or ONT direct-RNA). The platform is auto-detected per run "
                        "from the runlist (PACBIO_SMRT -> '-ax splice'; "
                        "OXFORD_NANOPORE -> '-ax splice -uf -k14'). Implies a "
                        "different splice-DB format and a smaller default --batch-size.")
    p.add_argument("--min-mapq", type=int, default=None,
                   help="MAPQ cutoff for the uniqueness gate (default: 60 for "
                        "short reads, 1 for --longreads).")
    p.add_argument("--advanced", nargs="*", default=[], metavar="KEY=VALUE",
                   help="Advanced overrides, e.g. lambda=10 pseudo-count=1 cost=0.001.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="varus",
        description="VARUS: online sampling of complementary RNA-seq reads from NCBI SRA.",
    )
    parser.add_argument("--log-level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    sub = parser.add_subparsers(dest="cmd", required=True)
    _add_runlist(sub)
    _add_index(sub)
    _add_run(sub)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=args.log_level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    if args.cmd == "runlist":
        from varus.runlist import fetch_runlist
        out = fetch_runlist(
            species=args.species,
            outdir=args.outdir,
            max_runs=args.max_runs,
            paired_only=args.paired_only,
            longreads=args.longreads,
            email=args.email,
            api_key=args.api_key,
        )
        print(out)
        return 0

    if args.cmd == "index":
        if args.longreads:
            from varus.index import build_minimap2_index
            out = build_minimap2_index(
                genome=args.genome,
                outdir=args.outdir,
                threads=args.threads,
                prefix=args.prefix or "mm2idx",
            )
        else:
            from varus.index import build_hisat2_index
            out = build_hisat2_index(
                genome=args.genome,
                outdir=args.outdir,
                threads=args.threads,
                prefix=args.prefix or "hisatidx",
            )
        print(out)
        return 0

    if args.cmd == "run":
        from varus.controller import Controller, VARUSConfig, load_runs
        import random

        # Parse --advanced KEY=VALUE overrides
        advanced: dict[str, str] = {}
        for kv in (args.advanced or []):
            if "=" in kv:
                k, v = kv.split("=", 1)
                advanced[k.strip()] = v.strip()

        # Mode-dependent defaults: long-read SRA runs are smaller and the MAPQ
        # distribution from minimap2 is wider than HISAT2's.
        batch_size = args.batch_size
        if batch_size is None:
            batch_size = 2_000 if args.longreads else 50_000
        min_mapq = args.min_mapq
        if min_mapq is None:
            min_mapq = 1 if args.longreads else 60

        cfg = VARUSConfig(
            genome=args.genome,
            index_prefix=args.index,
            outdir=args.outdir,
            batch_size=batch_size,
            max_batches=args.max_batches,
            tile_size=args.tile_size,
            min_uniq_pct=args.min_uniq_pct,
            threads=args.threads,
            keep_batches=args.keep_batches,
            coverage_trace=args.coverage_trace,
            seed=args.seed,
            bootstrap_all=args.bootstrap_all,
            profit_condition=args.profit_condition,
            pipeline_downloads=args.pipeline_downloads,
            longreads=args.longreads,
            min_mapq=min_mapq,
            lambda_=float(advanced.get("lambda", 10.0)),
            pseudo_count=float(advanced.get("pseudo-count", 1.0)),
            cost=float(advanced.get("cost", 0.0)),
        )

        rng = random.Random(cfg.seed)
        runs = load_runs(args.runlist, cfg.batch_size, rng)
        if not runs:
            raise SystemExit("Runlist is empty or all runs are colorspace / filtered.")
        ctrl = Controller(cfg, runs)
        ctrl.run()
        return 0

    raise SystemExit(f"unknown subcommand: {args.cmd}")


if __name__ == "__main__":
    sys.exit(main())
