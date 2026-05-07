# VARUS: Drawing Diverse Samples from RNA-Seq Libraries

> **NOTE — v2 Python rewrite in progress** on branch `nextflow_pipeline`. The
> original C++/Perl implementation has been moved to [`legacy/`](legacy/) and
> will be deleted once the new pipeline reaches feature parity. The algorithm
> (Stanke et al., 2019, [DOI](https://doi.org/10.1186/s12859-019-3182-x)) is
> unchanged.

**VARUS** automates the selection and download of a limited number of RNA-seq
reads from NCBI's Sequence Read Archive (SRA) targeting a sufficiently high
coverage for many genes for the purpose of gene-finder training and genome
annotation. Each iteration of the online algorithm

- selects a run to download that is expected to complement previously
  downloaded reads,
- downloads a sample of reads ("batch") with **fasterq-dump**,
- aligns the reads with **HISAT2**,
- evaluates the alignment.

## What changed in v2

| | v1 (legacy) | v2 (this branch) |
|---|---|---|
| Language | C++ + Perl + Bash | Python 3.9+ |
| Aligner | STAR or HISAT2 | HISAT2 (short reads) or minimap2 (long reads, `--longreads`) |
| Read download | `fastq-dump --fasta` | `fastq-dump` (per-batch ranges) + `fasterq-dump` (full runs) |
| Alignment intermediate | SAM → samtools sort → BAM | piped → coordinate-sorted BAM directly |
| Intron extraction | `bam2hints` (AUGUSTUS) | `pysam` reimplementation |
| Strand assignment | `filterIntronsFindStrand.pl` | `pyfaidx` reimplementation |
| Final merge | hierarchical bash scripts | `samtools merge` |
| Pipeline driver | `runVARUS.pl` + `VARUSparameters.txt` | Nextflow + `varus` Python CLI |
| Per-iteration coverage dump | always (~8 GB for 1000 batches) | off by default, `--coverage-trace N` |
| Per-batch FASTA kept gzipped | yes | deleted by default, `--keep-batches` to retain |
| User-facing parameters | ~25 in a parameters file | ~10 CLI flags + `--advanced KEY=VALUE` |

## Installation (v2 Python)

```sh
git clone https://github.com/Gaius-Augustus/VARUS.git
cd VARUS
git checkout nextflow_pipeline
pip install -e ".[align]"      # add ',dev' for the test suite
```

External tools that VARUS shells out to: `hisat2`, `hisat2-build`,
`samtools`, `fasterq-dump` (sra-toolkit). For long-read mode (`--longreads`),
add `minimap2`. Install via your distro or conda.

> Disable the NCBI cache once on the host:
> ```sh
> mkdir -p ~/.ncbi
> echo '/repository/user/cache-disabled = "true"' >> ~/.ncbi/user-settings.mkfg
> ```

## Quick start

Three commands:

```sh
# 1. Query NCBI SRA for all RNA-seq runs of the species
varus runlist "Schizosaccharomyces pombe" --outdir Sp/ --email you@host

# 2. Build a HISAT2 index of the genome
varus index   genome.fa --outdir Sp/genome/ --threads 8

# 3. Run the online sampling loop
varus run     "Schizosaccharomyces pombe" genome.fa \
              --runlist Sp/Runlist.tsv          \
              --index   Sp/genome/hisatidx      \
              --max-batches 1000 --threads 8    \
              --outdir  Sp/
```

Outputs in `Sp/`:

| File | Contents |
|---|---|
| `VARUS.bam` | merged coordinate-sorted alignment of all sampled batches |
| `introns.gff` | cumulative spliced-junction hints, strand-resolved |
| `Coverage.csv` | UMR count per 5 kb tile |
| `RunStatistics.csv` | per-run summary (downloads, UMR%, bad-quality flag) |

### Tuning knobs

| Flag | Default | Notes |
|---|---|---|
| `--batch-size` | 50000 | reads per batch |
| `--max-batches` | 1000 | hard upper bound on download iterations |
| `--tile-size` | 5000 | bp per coverage tile |
| `--min-uniq-pct` | 5.0 | reject batches below this UMR % (low-quality alignment) |
| `--threads` | 4 | alignment threads for HISAT2 / samtools |
| `--seed` | random | random seed for reproducible run order |
| `--bootstrap-all` | off | seed one batch from every run before the greedy loop |
| `--profit-condition` | off | stop early when expected marginal gain ≤ 0 |
| `--pipeline-downloads` | off | overlap round R+1 downloads with round R alignments (1.3–1.8× speedup) |
| `--coverage-trace N` | 0 (off) | snapshot Coverage every N batches |
| `--keep-batches` | off | retain per-batch FASTA/BAM after counting |
| `--advanced KEY=VALUE` | — | estimator hyperparameters: `lambda=10`, `pseudo-count=1`, `cost=0.0` |
| `--longreads` | off | align with minimap2 (long-read RNA-seq); see below |
| `--longread-platform` | `pacbio` | `pacbio` (Iso-Seq, `-ax splice`) or `ont` (direct-RNA, `-ax splice -uf -k14`) |
| `--min-mapq` | 60 / 1 | uniqueness MAPQ cutoff (default 60 short, 1 long) |

### Long-read RNA-seq (`--longreads`)

PacBio Iso-Seq and ONT direct-RNA runs from SRA are aligned with `minimap2 -ax splice`
instead of HISAT2. The same online algorithm runs on top — only the alignment, the
splice-DB feedback format, and a few defaults change.

```sh
# 1. Build a minimap2 splice index instead of HISAT2.
varus index   genome.fa --outdir Sp/genome/ --threads 8 --longreads

# 2. Run with --longreads (and pick the platform). The default --batch-size
#    drops from 50000 to 2000 because long-read SRA runs have far fewer spots.
varus run     "Schizosaccharomyces pombe" genome.fa     \
              --runlist Sp/Runlist.tsv                  \
              --index   Sp/genome/mm2idx.mmi            \
              --longreads --longread-platform ont       \
              --max-batches 1000 --threads 8 --outdir Sp/
```

Differences vs the HISAT2 path:

- `--index` points at the `.mmi` *file* rather than a stem.
- The splice-DB written each round is `intronDB.junc.bed` (BED12 for `minimap2 --junc-bed`)
  instead of `intronDB.splice_sites` (HISAT2 tab format).
- The uniqueness % is computed by scanning the BAM (primary, MAPQ ≥ `--min-mapq`),
  not parsed from a HISAT2-specific log file.

### Nextflow

A standalone Nextflow pipeline lives in [`nextflow/`](nextflow/):

```sh
nextflow run nextflow/main.nf \
  -c nextflow/example.config \
  --species_csv mycsv.csv \
  --outdir results \
  --ncbi_email you@host
```

`mycsv.csv` is a 2-column CSV: `species,genome` (one row per species). The
pipeline runs `VARUS_RUNLIST`, `VARUS_INDEX`, and `VARUS_RUN` in sequence per
species.

#### Nextflow params

| Param | Default | Notes |
|---|---|---|
| `--species_csv` | required | 2-column CSV: `species,genome` |
| `--outdir` | `results` | output root |
| `--ncbi_email` | — | contact email for NCBI Entrez (recommended) |
| `--ncbi_api_key` | — | raises NCBI rate limit to 10 req/s |
| `--varus_max_batches` | 1000 | passed to `varus run --max-batches` |
| `--varus_batch_size` | 50000 | passed to `varus run --batch-size` |
| `--varus_tile_size` | 5000 | passed to `varus run --tile-size` |
| `--varus_min_uniq_pct` | 5.0 | passed to `varus run --min-uniq-pct` |
| `--varus_max_runs` | 0 (all) | passed to `varus runlist --max-runs` |
| `--varus_seed` | 1 | passed to `varus run --seed` |
| `--varus_bootstrap_all` | false | passed to `varus run --bootstrap-all` |
| `--varus_profit_condition` | false | passed to `varus run --profit-condition` |
| `--varus_pipeline_downloads` | false | passed to `varus run --pipeline-downloads` |
| `--varus_index_cpus` | 8 | CPUs for `VARUS_INDEX` |
| `--varus_run_cpus` | 16 | CPUs for `VARUS_RUN` |
| `--longreads` | false | switch to minimap2 + restrict the SRA query to PacBio/ONT |
| `--longread_platform` | `pacbio` | `pacbio` or `ont`, only used with `--longreads` |

`VARUS_RUN` publishes one additional file per species: `runtime.varus.txt`
(`/usr/bin/time -p` wall/user/sys report).

The module can be imported into a larger workflow:

```groovy
include { VARUS_RUNLIST; VARUS_INDEX; VARUS_RUN } from '/path/to/VARUS/nextflow/varus.nf'
```

## Citation

Please cite:
[VARUS: sampling complementary RNA reads from the sequence read archive](https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/s12859-019-3182-x).
2019; *BMC Bioinformatics*, 20:558.

## Legacy

The original C++ implementation, the Perl wrappers, and the bash merge scripts
are preserved verbatim under [`legacy/`](legacy/) for parity testing during the
rewrite. They are not built or installed by default and will be removed in a
future release. See `legacy/README` (the original README, kept) for v1 usage.
