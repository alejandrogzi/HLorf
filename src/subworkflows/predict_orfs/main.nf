include { PREDICT }     from '../../modules/predict/main.nf'

workflow PREDICT_ORFS {
    take:
    ch_candidates
    ch_blast_counts

    main:
    ch_versions = Channel.empty()

    PREDICT(ch_candidates)

    ch_blast_counts
    .join(PREDICT.out.counts)
    .set { ch_counts }

    ch_versions = ch_versions.mix(PREDICT.out.versions)

    emit:
    orfs = PREDICT.out.orfs
    counts = ch_counts
    versions = ch_versions
}
