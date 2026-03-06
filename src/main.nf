#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Copyright (c) 2025 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EMAIL }        from './modules/email/main.nf'
include { XORF }         from './subworkflows/xorf/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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
    XORF (
       Channel.fromPath(params.regions).map { it -> [ [id: it.baseName], it ] },
       Channel.fromPath(params.sequence),
       Channel.fromPath(params.database),
       params.outdir,
       params.predict_keep_raw,
       params.chunk_size
    )

    PIPELINE_COMPLETION (
        params.email_to,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.use_mailx,
        XORF.out.files,
        XORF.out.counts,
        XORF.out.versions
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
