process CONCAT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:24.04' :
        'quay.io/biocontainers/python:3.10' }"

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
    awk 'NR==1 || FNR>1' tsv_*/*.tsv > all.tsv

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
