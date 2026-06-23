# Contributing

Thanks for your interest in VARUS. Bug reports, pull requests, and small
patches are welcome.

## Reporting issues

Open a GitHub issue with:

- a short description of what you ran and what you expected,
- the exact command line (`varus run ...` or the Nextflow invocation),
- the relevant section of the log,
- versions of `python`, `hisat2` / `minimap2`, `samtools`, and `sra-toolkit`
  (`varus --version`, `hisat2 --version`, `samtools --version`,
  `fasterq-dump --version`).

## Development setup

```sh
git clone <REPO_URL>
cd VARUS
pip install -e ".[align,dev]"
pytest
```

`pysam` does not build on stock Windows; tests that build a BAM skip
automatically when `pysam` is not importable (see
[tests/conftest.py](tests/conftest.py)).

## Pull requests

- Keep changes focused -- one PR per behavior change.
- Add or update a test in [tests/](tests/) when you change non-trivial logic.
- Run `pytest` locally before pushing; CI runs the same suite on Linux.
- Match the existing code style (type hints, `from __future__ import
  annotations`, module-level `log = logging.getLogger(__name__)`).
- Do not commit generated artifacts (`*.egg-info/`, `__pycache__/`,
  `.pytest_cache/`). They are gitignored.

## Long runs

Sampling on a real genome takes hours to days and is not appropriate for CI.
For local sanity checks, use a small genome (e.g. *S. cerevisiae* in
[example/](example/)) and `--max-batches 10`.
