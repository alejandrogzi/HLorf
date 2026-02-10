#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Copyright (c) 2025 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CHUNKER }      from './modules/chunker/main.nf'
include { CONCAT }       from './modules/concat/main.nf'
include { EMAIL }        from './modules/email/main.nf'

include { PREDICT_ORFS } from './subworkflows/predict_orfs/main.nf'
include { GET_CANDIDATES } from './subworkflows/candidates/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HLORF {
    def ch_regions  = Channel.fromPath(params.regions)
    def ch_sequence = Channel.fromPath(params.sequence)
    def ch_database = Channel.fromPath(params.database)
    def chunkSize   = params.chunk_size ?: 250

    def ch_versions = Channel.empty()

    CHUNKER(
        ch_regions,
        ch_sequence,
        chunkSize,
    )

    CHUNKER.out.chunked_regions
        .flatten()
        .map { region -> [[id: region.baseName], region] }
        .join(
            CHUNKER.out.chunked_sequences
                .flatten()
                .map { fasta -> [[id: fasta.baseName], fasta] }
        )
        .set { ch_pairs }

    ch_versions = ch_versions.mix(CHUNKER.out.versions)

    GET_CANDIDATES(
        ch_pairs,
        ch_database
    )

    PREDICT_ORFS(
        GET_CANDIDATES.out.candidates,
        GET_CANDIDATES.out.counts
    )

    PREDICT_ORFS.out.orfs
        .toList()
        .map { items ->
            def beds = items.collect { meta, bed, tsv -> bed }
            def tsvs = items.collect { meta, bed, tsv -> tsv }
            tuple([id: 'all'], beds, tsvs)
        }
        .set { ch_all }

    ch_versions = ch_versions.mix(GET_CANDIDATES.out.versions)
    ch_versions = ch_versions.mix(PREDICT_ORFS.out.versions)

    CONCAT(
        ch_all
    )

    PREDICT_ORFS.out.counts
    .map { meta, initial, netstart, transaid, ns_td, tai, blast, samba, all, unique, kept ->
        def line = "${meta.id}\t${initial}\t${netstart}\t${transaid}\t${ns_td}\t${tai}\t${blast}\t${samba}\t${all}\t${unique}\t${kept}"
        return line
    }
    .collectFile(
      name: 'counts.tsv',
      newLine: true,
      storeDir: "${params.outdir}/samplesheets"
    )
    .set { ch_counts }

    ch_versions = ch_versions.mix(CONCAT.out.versions)
    ch_pipeline_versions = ch_versions
        .collectFile(
            name:      "HLorf.versions.yml",
            storeDir:  "${params.outdir}/pipeline_info",
            sort:      true,
            keepHeader: false,
            newLine:   true
        )

    emit:
    files = CONCAT.out.files
    counts = ch_counts
    versions = ch_pipeline_versions
}

workflow PIPELINE_COMPLETION {

    take:
    email
    email_on_fail
    plaintext_email
    outdir
    use_mailx
    files
    counts
    ch_versions

    main:

    if (params.sent_email) {
        EMAIL (
            email,
            email_on_fail,
            plaintext_email,
            outdir,
            use_mailx,
            files,
            counts,
            ch_versions
        )
    }

    workflow.onError {
        log.error "ERROR: Pipeline failed!"
        log.error "ERROR: Please check the following error message:\n${workflow.errorMessage}"
        log.error "ERROR: Refer to github issues: https://github.com/alejandrogzi/HLorf/issues"
    }

    workflow.onComplete {
        log.info "\nPipeline completed successfully!"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    HLORF ()

    PIPELINE_COMPLETION (
        params.email_to,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.use_mailx,
        HLORF.out.files,
        HLORF.out.counts,
        HLORF.out.versions
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
