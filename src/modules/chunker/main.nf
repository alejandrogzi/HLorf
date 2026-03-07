process CHUNKER {
    tag "chunk:$meta.id:$meta.name"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/alejandrogzi/orf-chunk:latest'

    input:
    tuple val(meta), path(regions)
    path(sequence)
    val(chunk_size)

    output:
    path('tmp'),          emit: chunks
    tuple val(meta), path('tmp/*bed'),     emit: chunked_regions
    tuple val(meta), path('tmp/*fa'),      emit: chunked_sequences
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def upstream = task.ext.upstream ?: 1000
    def downstream = task.ext.downstream ?: 1000
    """
    orf chunk \\
    --regions $regions \\
    --sequence $sequence \\
    --chunks $chunk_size \\
    -u $upstream \\
    -d $downstream \\
    --ignore-errors

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orf-chunk: \$(orf --version 2>&1 | sed 's/^.*orf //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch tmp
    touch tmp/*bed
    touch tmp/*fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orf-chunk: \$(orf --version 2>&1 | sed 's/^.*orf //; s/ .*\$//')
    END_VERSIONS
    """
}
