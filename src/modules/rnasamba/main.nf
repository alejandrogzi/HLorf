process RNASAMBA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'alejandrogzi/orf-samba:latest'

    input:
    tuple val(meta), path(bed), path(sequence)

    output:
    tuple val(meta), path("${meta.id}/*tsv"), emit: samba
    tuple val(meta), path("${meta.id}*strip.fa"), emit: fasta
    tuple val(meta), path(bed), emit: bed
    path "versions.yml", emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orf-samba: \$(orf --version 2>&1 | sed 's/^.*orf //; s/ .*\$//')
        rnasamba: \$(rnasamba --version 2>&1 | tail -n 1 | sed 's/^rnasamba //')
    END_VERSIONS
    """

    stub:
    """
    touch *strip.fa
    touch ${meta.id}
    touch ${meta.id}/samba
    touch ${meta.id}/samba/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orf-samba: \$(orf --version 2>&1 | sed 's/^.*orf //; s/ .*\$//')
        rnasamba: \$(rnasamba --version 2>&1 | tail -n 1 | sed 's/^rnasamba //')
    END_VERSIONS
    """
}
