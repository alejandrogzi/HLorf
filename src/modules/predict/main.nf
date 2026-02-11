process PREDICT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'alejandrogzi/orf-predict:latest'

    input:
    tuple val(meta), path(bed), path(blast), path(samba)

    output:
    tuple val(meta), path("${meta.id}/*bed"), path("${meta.id}/*tsv"), emit: orfs
    tuple val(meta), env(BLAST_PREDICTION_COUNT), env(SAMBA_PREDICTION_COUNT), env(ALL_PREDICTED_REGIONS), env(UNIQUE_PREDICTED_REGIONS), env(KEPT_REGIONS), emit: counts
    path "versions.yml", emit: versions

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

    BLAST_PREDICTION_COUNT=\$(wc -l < ${blast})
    SAMBA_PREDICTION_COUNT=\$(wc -l < ${samba})

    ALL_PREDICTED_REGIONS=\$(wc -l ${meta.id}/${meta.id}.predictions.tsv | awk '{print \$1}')
    UNIQUE_PREDICTED_REGIONS=\$(awk '{print \$4}' ${meta.id}/${meta.id}.predictions.tsv | sed -E 's/(_ORF|\\.p[0-9]).*//' | sort | uniq | wc -l)
    KEPT_REGIONS=\$(wc -l ${meta.id}/${meta.id}.predictions.bed | awk '{print \$1}')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        predict.py: \$(predict.py --version 2>&1 | sed 's/.*Version: //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}
    touch ${meta.id}/*bed
    touch ${meta.id}/*tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        predict.py: \$(predict.py --version 2>&1 | sed 's/.*Version: //')
    END_VERSIONS
    """
}
