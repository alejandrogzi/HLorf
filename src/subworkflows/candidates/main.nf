include { TRANSLATION } from '../../modules/translation/main.nf'
include { RNASAMBA }    from '../../modules/rnasamba/main.nf'
include { NETSTART }    from '../../modules/netstart/main.nf'
include { TRANSAID }    from '../../modules/transaid/main.nf'
include { BLAST }       from '../../modules/blast/main.nf'
include { JOIN as JOIN_NETS }   from '../../modules/join/main.nf'

workflow GET_CANDIDATES {
    take:
    ch_pairs
    ch_database

    main:
    ch_versions = Channel.empty()

    TRANSLATION(ch_pairs)
    RNASAMBA(ch_pairs)

    NETSTART(
      RNASAMBA.out.fasta,
      RNASAMBA.out.bed
    )

    TRANSAID(
      RNASAMBA.out.fasta,
      RNASAMBA.out.bed
    )

    JOIN_NETS(
        NETSTART.out.netstart,
        TRANSAID.out.transaid,
        RNASAMBA.out.bed
    )

    BLAST(
        TRANSLATION.out.predictions,
        JOIN_NETS.out.net,
        ch_database
    )
    .blast
    .join(RNASAMBA.out.samba)
    .set { ch_candidates }

    BLAST.out.counts
    .join(NETSTART.out.count)
    .join(TRANSAID.out.count)
    .join(JOIN_NETS.out.count)
    .set { ch_counts }

    ch_versions = ch_versions.mix(TRANSLATION.out.versions)
    ch_versions = ch_versions.mix(RNASAMBA.out.versions)
    ch_versions = ch_versions.mix(BLAST.out.versions)

    emit:
    candidates = ch_candidates
    counts = ch_counts
    versions = ch_versions
}
