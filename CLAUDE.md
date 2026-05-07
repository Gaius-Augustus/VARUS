# Project rules

Bioinformatics group. Unix tools, Snakemake, Nextflow, PyTorch/TensorFlow.
Read this fully before your first action.

## Hard rules (inviolable)

1. Never modify, move, or delete anything under `raw/`, `data/raw/`, or `*_raw/`.
2. Never run `git commit`, `git push`, `reset --hard`, or `rebase`. Propose
   commit messages; the user commits.
3. Never fabricate citations. Delegate to `citation-fetcher`; if it cannot
   verify, say so.
4. Before any real-data run, use `dry-run-first`. Before HPC submission, use
   `hpc-runner`. Cluster time is not free.
5. No compute on HPC login nodes. Mother scheduler processes (Snakemake /
   Nextflow with all rules using SLURM) are the only exception.
6. No destructive shell (`rm -rf`, `find -delete`, `dd`) without explicit ok.

## Routing

- Need a reference → `citation-fetcher` (writes `papers/<citekey>.md` in vault)
- About to run expensive code → `dry-run-first`, then `hpc-runner`
- Need tests → `test-writer`
- Something failing → `debugger` (no HPC until it's fixed locally)
- Changed module behavior → `docs-sync` (also updates `repos/<name>.md` in vault)
- Exploring local code → read `graphify-out/GRAPH_REPORT.md` first (if present); search raw files only after. Graph is built with `/graphify ./scripts --obsidian` (or relevant subfolder — never on `raw/`).
- After structural changes (new file, renamed symbol, changed interface) → if `graphify-out/` exists, run `/graphify --update --obsidian` to keep the graph current.
- Need group context (papers, methods, cross-repo decisions) → `/knowledge lookup`. Vault path in `additionalDirectories`.

## Style

Absolute paths on HPC. Cite with DOI. No emojis. Short answers. End file-
altering replies with a one-line audit list of touched files.
