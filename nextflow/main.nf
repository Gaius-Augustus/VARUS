#!/usr/bin/env nextflow
/*
 * Standalone VARUS v2 pipeline.
 *
 *   Input : a CSV with columns species,genome  (paths to per-species genome FASTAs)
 *   Output: per species  →  VARUS.bam, introns.gff, Coverage.csv, RunStatistics.csv
 *
 * Example:
 *   nextflow run nextflow/main.nf \
 *     --species_csv  examples/species.csv \
 *     --outdir       results \
 *     --ncbi_email   you@host
 */

nextflow.enable.dsl = 2

include { VARUS_RUNLIST; VARUS_INDEX; VARUS_RUN } from './varus.nf'


// ---------------------------- params ----------------------------

params.species_csv          = params.species_csv          ?: null
params.outdir               = params.outdir               ?: 'results'

// VARUS run-time hyperparameters (all optional)
params.varus_max_batches    = (params.containsKey('varus_max_batches')   && params.varus_max_batches   != null ? params.varus_max_batches   : 1000) as int
params.varus_batch_size     = (params.containsKey('varus_batch_size')    && params.varus_batch_size    != null ? params.varus_batch_size    : 50000) as int
params.varus_tile_size      = (params.containsKey('varus_tile_size')     && params.varus_tile_size     != null ? params.varus_tile_size     : 5000) as int
params.varus_min_uniq_pct   = (params.containsKey('varus_min_uniq_pct')  && params.varus_min_uniq_pct  != null ? params.varus_min_uniq_pct  : 5.0) as double
params.varus_max_runs       = (params.containsKey('varus_max_runs')      && params.varus_max_runs      != null ? params.varus_max_runs      : 0) as int
params.varus_seed           = (params.containsKey('varus_seed')          && params.varus_seed          != null ? params.varus_seed          : 1) as int
params.varus_bootstrap_all  = (params.containsKey('varus_bootstrap_all') ? params.varus_bootstrap_all : false) as boolean
params.varus_profit_condition = (params.containsKey('varus_profit_condition') ? params.varus_profit_condition : false) as boolean
params.varus_pipeline_downloads = (params.containsKey('varus_pipeline_downloads') ? params.varus_pipeline_downloads : false) as boolean
params.varus_index_cpus     = (params.containsKey('varus_index_cpus')    && params.varus_index_cpus    != null ? params.varus_index_cpus    : 8) as int
params.varus_run_cpus       = (params.containsKey('varus_run_cpus')      && params.varus_run_cpus      != null ? params.varus_run_cpus      : 16) as int

params.ncbi_email           = params.ncbi_email   ?: null
params.ncbi_api_key         = params.ncbi_api_key ?: null


def die(msg) { log.error msg; System.exit(1) }

if (!params.species_csv) die("Missing --species_csv")


// ---------------------------- workflow ----------------------------

workflow {

    ch_input = Channel.fromPath(params.species_csv, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(row.species, file(row.genome), [:]) }
        // (species, genome_path, extra) — `extra` is a placeholder so callers
        // wiring this into a larger workflow can carry context downstream.

    runlist_out = VARUS_RUNLIST(ch_input)
    index_out   = VARUS_INDEX(runlist_out)
    bam_out     = VARUS_RUN(index_out).bam
}
