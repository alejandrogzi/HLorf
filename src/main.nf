//!/usr/bin/env nextflow

// params.regions = "/home/alejandro/Documents/projects/hiller/champagne/pisco/files/hg38/TOGA/vs_ROSCfam/query_annotation.bed"
params.regions = "/home/alejandro/Documents/projects/HLorf/test.bed"
params.sequence = "/home/alejandro/Documents/projects/hiller/champagne/pisco/files/ROSCfam/fasta/ROSCfam.2bit"
params.database = "/home/alejandro/Documents/projects/HLorf/swissprot_vertebrates.dmnd"

process CHUNKER {
    container 'orf-chunk:latest'
    containerOptions '--pull=never'

    input:
    path(regions)
    path(sequence)
    val(chunk_size)

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
    orf chunk \\
    --regions $regions \\
    --sequence $sequence \\
    --chunks $chunk_size \\
    -u $upstream \\
    -d $downstream
    """

    stub:
    """
    touch tmptai --fasta tmp_0.fa --bed tmp_0.bed --outdir tmp_0
    touch tmp/*bed
    touch tmp/*fa
    """
}

process TRANSLATION {
    container 'orf-tai:latest'
    containerOptions '--pull=never'

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
    orf tai \\
    --fasta $sequence \\
    --bed $bed \\
    --outdir ${meta.id} \\
    -u $upstream \\
    -d $downstream
    
    mv ${meta.id}/tai/*result ${meta.id}/ && rm -rf ${meta.id}/tai
    """

    stub:
    """
    touch ${meta.id}
    touch ${meta.id}/tai
    touch ${meta.id}/tai/*result
    """
}

process RNASAMBA {
    container 'orf-samba:latest'
    containerOptions '--pull=never'

    input:
    tuple val(meta), path(bed), path(sequence)

    output:
    tuple val(meta), path("${meta.id}/*tsv"), emit: samba

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def upstream = task.ext.upstream ?: 1000
    def downstream = task.ext.downstream ?: 1000
    """
    orf samba \\
    --fasta $sequence \\
    --outdir ${meta.id} \\
    --upstream-flank $upstream \\
    --downstream-flank $downstream \\
    $args

    mv ${meta.id}/samba/*tsv ${meta.id}/ && rm -rf ${meta.id}/samba
    rm *strip.fa
    """

    stub:
    """
    touch *strip.fa
    touch ${meta.id}
    touch ${meta.id}/samba
    touch ${meta.id}/samba/*
    """
}

process BLAST {
    container 'orf-blast:latest'
    containerOptions '--pull=never'

    input:
    tuple val(meta), path(bed), path(sequence), path(predictions)
    each path(database)

    output:
    tuple val(meta), path(bed), path("${meta.id}/*result"), emit: blast

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def orf_min_len = task.ext.orf_min_len ?: 50
    def orf_min_percent = task.ext.orf_min_percent ?: 0.25
    def upstream = task.ext.upstream ?: 1000
    def downstream = task.ext.downstream ?: 1000
    """
    orf blast \\
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

process PREDICT {
    container 'orf-predict:latest'
    containerOptions '--pull=never'

    input:
    tuple val(meta), path(bed), path(blast), path(samba)

    output:
    tuple val(meta), path("${meta.id}/*bed"), path("${meta.id}/*tsv"), emit: orfs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def threshold = task.ext.threshold ?: 0.03
    def min_score_max_predictions = task.ext.min_score_max_predictions ?: 0.70
    def max_predictions = task.ext.max_predictions ?: 1
    """
    predict.py \\
    --blast $blast \\
    --samba $samba \\
    --alignments $bed \\
    --outdir ${meta.id} \\
    --prefix ${meta.id} \\
    --threshold $threshold \\
    --min-score-max-predictions $min_score_max_predictions \\
    --max-predictions $max_predictions
    """

    stub:
    """
    touch ${meta.id}
    touch ${meta.id}/*bed
    touch ${meta.id}/*tsv
    """
}

process CONCAT {
    input:
    tuple val(meta), path(beds, stageAs: 'bed_?/*'), path(tsvs, stageAs: 'tsv_?/*')

    output:
    tuple path("*bed"), path("*tsv"), emit: files

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat bed_*/*.bed > all.bed
    cat tsv_*/*.tsv > all.tsv    
    """

    stub:
    """
    touch all.*
    """
}

workflow {
    ch_regions = Channel.fromPath(params.regions)
    ch_sequence = Channel.fromPath(params.sequence)

    CHUNKER(
      ch_regions, 
      ch_sequence,
      250,
    )

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

    RNASAMBA(
        ch_pairs,
    )

    BLAST(
      TRANSLATION.out.predictions,
      Channel.fromPath(params.database)
    )
    .join(RNASAMBA.out.samba)
    .set { ch_candidates }

    PREDICT(
      ch_candidates
    )

    PREDICT.out.orfs
        .toList()
        .map { items ->
            def beds = items.collect { meta, bed, tsv -> bed }
            def tsvs = items.collect { meta, bed, tsv -> tsv }
            tuple([id: 'all'], beds, tsvs)
        }
        .set { ch_all }

    ch_all.view()

    CONCAT(
      ch_all
    )

    CONCAT.out.files.view()
}
