// VARUS v2 Nextflow module — three reusable processes that wrap the Python
// CLI. Designed to be `include`d from a parent workflow:
//
//     include { VARUS_RUNLIST; VARUS_INDEX; VARUS_RUN } from '/path/to/VARUS/nextflow/varus.nf'
//
// Each process expects the `varus` CLI on $PATH (`pip install -e .[align]`)
// plus `hisat2`, `hisat2-build`, `samtools`, and `fastq-dump`.
//
// The three processes feed each other:
//
//     VARUS_RUNLIST -> Runlist.tsv (NCBI Entrez query for the species)
//     VARUS_INDEX   -> HISAT2 index of the genome
//     VARUS_RUN     -> online loop: download SRA batches, align, score tiles
//
// Inputs are passed as a single tuple beginning with `species` and `genome`;
// callers may extend the tuple with arbitrary trailing fields, which are
// forwarded unchanged to make plumbing into larger pipelines easy.

process VARUS_RUNLIST {
    tag { species }
    publishDir { "${params.outdir}/${species.replaceAll(' ', '_')}/varus" }, mode: 'copy', overwrite: true
    cpus 1

    input:
        tuple val(species), path(genome), val(extra)

    output:
        tuple val(species), path(genome), path("Runlist.tsv"), val(extra)

    script:
    def maxRuns  = params.varus_max_runs ?: 0
    def email    = params.ncbi_email ?: ''
    def apiKey   = params.ncbi_api_key ?: ''
    def emailArg = email  ? "--email ${email}"     : ''
    def keyArg   = apiKey ? "--api-key ${apiKey}"  : ''
    """
    set -euo pipefail
    varus runlist '${species}' \\
        --outdir . \\
        --max-runs ${maxRuns} \\
        ${emailArg} ${keyArg}
    test -s Runlist.tsv || { echo "Runlist.tsv is empty for '${species}'" >&2; exit 2; }
    """

    stub:
    """
    printf '@Run_acc\\ttotal_spots\\ttotal_bases\\tavg_len\\tbool:paired\\tcolor_space\\nSRR000001\\t100000\\t10000000\\t100.0\\t0\\t0\\n' > Runlist.tsv
    """
}


process VARUS_INDEX {
    tag { species }
    publishDir { "${params.outdir}/${species.replaceAll(' ', '_')}/varus" }, mode: 'copy', overwrite: true
    cpus { params.varus_index_cpus ?: 8 }

    input:
        tuple val(species), path(genome), path(runlist), val(extra)

    output:
        tuple val(species), path(genome), path(runlist),
              path("genome_index"), val(extra)

    script:
    """
    set -euo pipefail
    mkdir -p genome_index
    varus index ${genome} \\
        --outdir genome_index \\
        --threads ${task.cpus} \\
        --prefix hisatidx
    """

    stub:
    """
    mkdir -p genome_index
    touch genome_index/hisatidx.1.ht2
    """
}


process VARUS_RUN {
    tag { species }
    publishDir { "${params.outdir}/${species.replaceAll(' ', '_')}/varus" }, mode: 'copy', overwrite: true
    cpus { params.varus_run_cpus ?: 16 }

    input:
        tuple val(species), path(genome), path(runlist),
              path(index_dir), val(extra)

    output:
        tuple val(species), path(genome), path("VARUS.bam"), val(extra), emit: bam
        path "introns.gff",       optional: true,                       emit: introns
        path "Coverage.csv",      optional: true,                       emit: coverage
        path "RunStatistics.csv", optional: true,                       emit: stats
        path "runtime.varus.txt",                                       emit: runtime

    script:
    def maxBatches  = params.varus_max_batches ?: 1000
    def batchSize   = params.varus_batch_size  ?: 50000
    def tileSize    = params.varus_tile_size   ?: 5000
    def minUniqPct  = params.varus_min_uniq_pct ?: 5.0
    def seed        = params.varus_seed         ?: 1
    def bootstrap   = params.varus_bootstrap_all ? '--bootstrap-all' : ''
    def profitCond  = params.varus_profit_condition ? '--profit-condition' : ''
    def pipelineDl  = params.varus_pipeline_downloads ? '--pipeline-downloads' : ''
    """
    set -euo pipefail
    /usr/bin/time -p -o runtime.varus.txt \\
      varus run '${species}' ${genome} \\
        --runlist ${runlist} \\
        --index ${index_dir}/hisatidx \\
        --outdir . \\
        --batch-size ${batchSize} \\
        --max-batches ${maxBatches} \\
        --tile-size ${tileSize} \\
        --min-uniq-pct ${minUniqPct} \\
        --threads ${task.cpus} \\
        --seed ${seed} \\
        ${bootstrap} ${profitCond} ${pipelineDl}

    test -s VARUS.bam || { echo "VARUS run produced no BAM" >&2; exit 2; }
    """

    stub:
    """
    touch VARUS.bam runtime.varus.txt introns.gff Coverage.csv RunStatistics.csv
    """
}
