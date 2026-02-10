process JOIN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'alejandrogzi/orf-chunk:latest'

    input:
    tuple val(meta), path(bed), path(netstart)
    tuple val(meta), path(transaid)

    output:
    tuple val(meta), path('net/merged.net'), emit: net
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
