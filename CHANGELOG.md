# Changelog

All notable changes to VARUS are documented in this file.

## [2.0.0a0] -- unreleased

Full Python rewrite of the C++/Perl implementation. The online sampling
algorithm (Stanke et al., 2019,
[DOI:10.1186/s12859-019-3182-x](https://doi.org/10.1186/s12859-019-3182-x))
is unchanged.

### Changed vs. v1 (upstream [Gaius-Augustus/VARUS](https://github.com/Gaius-Augustus/VARUS))

| | v1 | v2 |
|---|---|---|
| Language | C++ + Perl + Bash | Python 3.9+ |
| Aligner | STAR or HISAT2 | HISAT2 (short reads), minimap2 (long reads, `--longreads`) |
| Read download | `fastq-dump --fasta` | `fastq-dump` (per-batch ranges) + `fasterq-dump` (full runs) |
| Alignment intermediate | SAM -> samtools sort -> BAM | piped -> coordinate-sorted BAM directly |
| Intron extraction | `bam2hints` (AUGUSTUS) | `pysam` reimplementation |
| Strand assignment | `filterIntronsFindStrand.pl` | `pyfaidx` reimplementation |
| Final merge | hierarchical bash scripts | `samtools merge` |
| Pipeline driver | `runVARUS.pl` + `VARUSparameters.txt` | Nextflow + `varus` Python CLI |
| Per-iteration coverage dump | always (~8 GB for 1000 batches) | off by default, `--coverage-trace N` |
| Per-batch FASTA kept gzipped | yes | deleted by default, `--keep-batches` to retain |
| User-facing parameters | ~25 in a parameters file | ~10 CLI flags + `--advanced KEY=VALUE` |

### Added

- Long-read RNA-seq support (`--longreads`): minimap2 alignment with per-run
  preset selection from SRA platform metadata; BED12 splice-DB feedback;
  uniqueness % from BAM scan instead of HISAT2 log parsing.
- Standalone Nextflow pipeline in [`nextflow/`](nextflow/).
- Entrez retry-with-backoff on HTTP 429 / 5xx in `varus runlist`.

### Removed

- Per-iteration coverage dump is now opt-in (`--coverage-trace N`).
- Per-batch FASTA/BAM are deleted after counting (`--keep-batches` to retain).
