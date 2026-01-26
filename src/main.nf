#!/usr/bin/env nextflow

// params.regions = "/home/alejandro/Documents/projects/hiller/champagne/pisco/files/hg38/TOGA/vs_ROSCfam/query_annotation.bed"
params.regions = "/home/alejandro/Documents/projects/HLorf/test.bed"
params.sequence = "/home/alejandro/Documents/projects/hiller/champagne/pisco/files/ROSCfam/fasta/ROSCfam.2bit"
params.database = "/home/alejandro/Documents/projects/HLorf/swissprot_vertebrates.dmnd"

process CHUNKER {
    input:
    path(regions)
    path(sequence)

    output:
    path('tmp') , emit: chunks
    path('tmp/*bed') , emit: chunked_regions
    path('tmp/*fa') , emit: chunked_sequences

    when:
    task.ext.when == null || task.ext.when

    script:
    def upstream = task.ext.upstream ?: 1000
    def downstream = task.ext.downstream ?: 1000
    """
    chunker --regions $regions --sequence $sequence --chunks 250 -u $upstream -d $downstream
    """

    stub:
    """
    touch tmptai --fasta tmp_0.fa --bed tmp_0.bed --outdir tmp_0
    touch tmp/*bed
    touch tmp/*fa
    """
}

process TRANSLATION {
    input:
    tuple val(meta), path(bed), path(sequence)

    output:
    tuple val(meta), path(bed), path(sequence), path("${meta.id}/*result"), emit: predictions

    when:
    task.ext.when == null || task.ext.when

    script:
    def upstream = task.ext.upstream ?: 1000
    def downstream = task.ext.downstream ?: 1000
    """
    tai --fasta $sequence --bed $bed --outdir ${meta.id} -u $upstream -d $downstream
    mv ${meta.id}/tai/*result ${meta.id}/ && rm -rf ${meta.id}/tai
    """

    stub:
    """
    touch ${meta.id}
    touch ${meta.id}/tai
    touch ${meta.id}/tai/*result
    """
}

process RNASAMBA {}

process BLAST {
    input:
    tuple val(meta), path(bed), path(sequence), path(predictions)
    val(database)

    output:
    tuple val(meta), path("${meta.id}/*result"), emit: orfs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def orf_min_len = task.ext.orf_min_len ?: 50
    def orf_min_percent = task.ext.orf_min_percent ?: 0.25
    def upstream = task.ext.upstream ?: 1000
    def downstream = task.ext.downstream ?: 1000
    """
    blast \\
    --fasta $sequence \\
    --bed $bed \\
    --tai $predictions \\
    --outdir ${meta.id} \\
    --orf-min-len $orf_min_len \\
    --orf-min-percent $orf_min_percent \\
    --database $database \\
    --upstream-flank $upstream \\
    --downstream-flank $downstream \\
    $args

    mv ${meta.id}/orf/*result ${meta.id}/ && rm -rf ${meta.id}/orf
    """

    stub:
    """
    touch ${meta.id}
    touch ${meta.id}/orf
    touch ${meta.id}/orf/*
    """
}

workflow {
    ch_regions = Channel.fromPath(params.regions)
    ch_sequence = Channel.fromPath(params.sequence)

    CHUNKER(ch_regions, ch_sequence)

    CHUNKER.out.chunked_regions
        .flatten()
        .map { it -> [[id: it.baseName], it] }
        .join(
            CHUNKER.out.chunked_sequences
                .flatten()
                .map { it -> [[id: it.baseName], it] }
        )
        .set { ch_pairs }

    TRANSLATION(
      ch_pairs
    )

    TRANSLATION.out.predictions.view()

     BLAST(
        TRANSLATION.out.predictions,
        params.database
     )
     .map { meta, prediction -> prediction }
     .collect()
     .set { ch_predictions }

     ch_predictions.view()
}
