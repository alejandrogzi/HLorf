process CONCAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'ubuntu:24.04'

    input:
    tuple val(meta), path(beds, stageAs: 'bed_?/*'), path(tsvs, stageAs: 'tsv_?/*')

    output:
    tuple path("*bed"), path("*tsv"), emit: files
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat bed_*/*.bed > all.bed
    cat tsv_*/*.tsv > all.tsv    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(cat --version 2>&1 | head -n 1 | sed 's/cat (GNU coreutils) //' )
    END_VERSIONS
    """

    stub:
    """
    touch all.*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(cat --version 2>&1 | head -n 1 | sed 's/cat (GNU coreutils) //' )
    END_VERSIONS
    """
}
