process TRANSAID {
    tag "$meta.id:$meta.name"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/alejandrogzi/orf-net:latest'

    input:
    tuple val(meta), path(sequence)
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*csv"), emit: transaid
    tuple val(meta), env(PREDICTION_COUNT), emit: count
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    transaid \\
    --input $sequence \\
    --gpu -1 \\
    --output ${meta.id}_transaid \\
    $args

    mv *csv ${meta.id}.${meta.name}.transaid.csv
    PREDICTION_COUNT=\$(wc -l < ${meta.id}.${meta.name}.transaid.csv)
    
    rm *.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transaid: \$(transaid --version 2>&1 | sed 's/.*Version: //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orf-samba: \$(orf --version 2>&1 | sed 's/^.*orf //; s/ .*\$//')
        rnasamba: \$(rnasamba --version 2>&1 | tail -n 1 | sed 's/^rnasamba //')
    END_VERSIONS
    """
}
