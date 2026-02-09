process NETSTART {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'alejandrogzi/orf-net:latest'

    input:
    tuple val(meta), path(bed), path(sequence)

    output:
    tuple val(meta), path(bed), path("${meta.id}*csv"), emit: netstart
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    netstart2 \\
    -in $sequence \\
    -compute_device cpu \\
    -o chordata \\
    -out ${meta.id}
    $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        netstart2: \$(netstart2 --version 2>&1 | sed 's/.*Version: //')
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
