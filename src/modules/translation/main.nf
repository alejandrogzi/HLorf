process TRANSLATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'alejandrogzi/orf-tai:latest'

    input:
    tuple val(meta), path(bed), path(sequence)

    output:
    tuple val(meta), path(bed), path(sequence), path("${meta.id}/*result"), emit: predictions
    path "versions.yml", emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orf-tai: \$(orf --version 2>&1 | sed 's/^.*orf //; s/ .*\$//')
        translationai: 0.0.1
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}
    touch ${meta.id}/tai
    touch ${meta.id}/tai/*result

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orf-tai: \$(orf --version 2>&1 | sed 's/^.*orf //; s/ .*\$//')
        translationai: 0.0.1
    END_VERSIONS
    """
}
