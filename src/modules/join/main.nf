process JOIN {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'alejandrogzi/orf-chunk:latest'

    input:
    tuple val(meta), path(netstart)
    tuple val(meta), path(transaid)
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path('net/merged.net'), emit: net
    tuple val(meta), env(PREDICTION_COUNT), emit: count
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    orf net \\
    --bed $bed \\
    --netstart $netstart \\
    --transaid $transaid \\
    --outdir .

    PREDICTION_COUNT=\$(cat net/merged.net | wc -l)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orf-net: \$(orf --version 2>&1 | sed 's/^.*orf //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch net
    touch net/merged.net

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orf-chunk: \$(orf --version 2>&1 | sed 's/^.*orf //; s/ .*\$//')
    END_VERSIONS
    """
}
