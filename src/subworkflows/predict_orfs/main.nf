include { PREDICT }     from '../../modules/predict/main.nf'

workflow PREDICT_ORFS {
    take:
    ch_candidates
    ch_blast_counts

    main:
    ch_versions = Channel.empty()
    ch_raw = Channel.empty()

    PREDICT(ch_candidates)

    if (params.predict_keep_raw) {
        PREDICT.out.raw.set { ch_raw }
    }

    ch_blast_counts
    .join(PREDICT.out.counts)
    .set { ch_counts }

    ch_versions = ch_versions.mix(PREDICT.out.versions)

    emit:
    orfs = PREDICT.out.orfs
    raw = ch_raw
    counts = ch_counts
    versions = ch_versions
}
