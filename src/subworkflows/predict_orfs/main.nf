include { TRANSLATION } from '../../modules/translation/main.nf'
include { RNASAMBA }    from '../../modules/rnasamba/main.nf'
include { BLAST }       from '../../modules/blast/main.nf'
include { PREDICT }     from '../../modules/predict/main.nf'

workflow PREDICT_ORFS {
    take:
    ch_pairs
    ch_database

    main:
    ch_versions = Channel.empty()

    TRANSLATION(ch_pairs)
    RNASAMBA(ch_pairs)

    BLAST(
        TRANSLATION.out.predictions,
        ch_database
    )
    .blast
    .join(RNASAMBA.out.samba)
    .set { ch_candidates }

    PREDICT(ch_candidates)

    BLAST.out.counts
    .join(PREDICT.out.counts)
    .set { ch_counts }

    ch_versions = ch_versions.mix(TRANSLATION.out.versions)
    ch_versions = ch_versions.mix(RNASAMBA.out.versions)
    ch_versions = ch_versions.mix(BLAST.out.versions)
    ch_versions = ch_versions.mix(PREDICT.out.versions)

    emit:
    orfs = PREDICT.out.orfs
    counts = ch_counts
    versions = ch_versions
}
