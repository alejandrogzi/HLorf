include { TRANSLATION } from '../../modules/translation/main.nf'
include { RNASAMBA }    from '../../modules/rnasamba/main.nf'
include { NETSTART }    from '../../modules/netstart/main.nf'
include { TRANSAID }    from '../../modules/transaid/main.nf'
include { BLAST }       from '../../modules/blast/main.nf'

workflow GET_CANDIDATES {
    take:
    ch_pairs
    ch_database

    main:
    ch_versions = Channel.empty()

    TRANSLATION(ch_pairs)
    RNASAMBA(ch_pairs)
    NETSTART(RNASAMBA.out.fasta)
    TRANSAID(RNASAMBA.out.fasta)

    BLAST(
        TRANSLATION.out.predictions,
        ch_database
    )
    .blast
    .join(RNASAMBA.out.samba)
    .set { ch_candidates }

    ch_versions = ch_versions.mix(TRANSLATION.out.versions)
    ch_versions = ch_versions.mix(RNASAMBA.out.versions)
    ch_versions = ch_versions.mix(BLAST.out.versions)

    emit:
    candidates = ch_candidates
    counts = BLAST.out.counts
    versions = ch_versions
}
