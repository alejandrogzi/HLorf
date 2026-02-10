process BLAST {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'alejandrogzi/orf-blast:latest'

    input:
    tuple val(meta), path(bed), path(sequence), path(predictions)
    tuple val(meta), path(net)
    each path(database)

    output:
    tuple val(meta), path(bed), path("${meta.id}/*result"), emit: blast
    tuple val(meta), env(INITIAL_REGION_COUNT), env(TRANSLATION_COUNT), emit: counts
    path "versions.yml", emit: versions

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
    --net $net \\
    --outdir ${meta.id} \\
    --orf-min-len $orf_min_len \\
    --orf-min-percent $orf_min_percent \\
    --database $database \\
    --upstream-flank $upstream \\
    --downstream-flank $downstream \\
    $args

    mv ${meta.id}/orf/*result ${meta.id}/ && rm -rf ${meta.id}/orf

    INITIAL_REGION_COUNT=\$(wc -l < ${bed})
    TRANSLATION_COUNT=\$(wc -l < ${predictions})

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orf-blast: \$(orf --version 2>&1 | sed 's/^.*orf //; s/ .*\$//')
        diamond: \$(diamond version 2>&1 | sed 's/^.*diamond version //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}
    touch ${meta.id}/orf
    touch ${meta.id}/orf/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orf-blast: \$(orf --version 2>&1 | sed 's/^.*orf //; s/ .*\$//')
        diamond: \$(diamond version 2>&1 | sed 's/^.*diamond version //' )
    END_VERSIONS
    """
}
