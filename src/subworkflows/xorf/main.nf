#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Copyright (c) 2025 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CHUNKER }      from '../../modules/chunker/main.nf'
include { CONCAT }       from '../../modules/concat/main.nf'
include { CONCAT as CONCAT_RAW }   from '../../modules/concat/main.nf'

include { PREDICT_ORFS } from '../predict_orfs/main.nf'
include { GET_CANDIDATES } from '../candidates/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow XORF {
    take:
      regions  
      sequence
      database
      output_dir
      predict_keep_raw
      chunk_size

    main:
      def ch_regions  = regions
      def ch_sequence = sequence
      def ch_database = database
      def chunkSize   = chunk_size ?: 20

      def ch_versions = Channel.empty()

      CHUNKER(
          ch_regions,
          ch_sequence,
          chunkSize,
      )

      CHUNKER.out.chunked_regions
          .flatMap { meta, region -> 
              def regions = region instanceof List ? region : [region]
              regions.collect { it ->
                [ meta + [name: it.baseName], it] }
              }
          .join(
              CHUNKER.out.chunked_sequences
                .flatMap { meta, fasta -> 
                    def fas = fasta instanceof List ? fasta : [fasta]
                    fas.collect { it ->
                      [ meta + [name: it.baseName], it] }
                }
          )
          .set { ch_pairs }

      ch_versions = ch_versions.mix(CHUNKER.out.versions)

      GET_CANDIDATES(
          ch_pairs,
          ch_database
      )

      PREDICT_ORFS(
          GET_CANDIDATES.out.candidates,
          GET_CANDIDATES.out.counts
      )

      PREDICT_ORFS.out.orfs
          .toList()
          .map { items ->
              def beds = items.collect { meta, bed, tsv -> bed }
              def tsvs = items.collect { meta, bed, tsv -> tsv }
              tuple([id: 'all'], beds, tsvs)
          }
          .set { ch_all }

      if (predict_keep_raw) {
          PREDICT_ORFS.out.raw
              .toList()
              .map { items ->
                  def beds = items.collect { meta, bed, tsv -> bed }
                  def tsvs = items.collect { meta, bed, tsv -> tsv }
                  tuple([id: 'raw'], beds, tsvs)
              }
              .set { ch_raw }

          CONCAT_RAW(
              ch_raw
          )
      }

      ch_versions = ch_versions.mix(GET_CANDIDATES.out.versions)
      ch_versions = ch_versions.mix(PREDICT_ORFS.out.versions)

      CONCAT(
          ch_all
      )

      PREDICT_ORFS.out.counts
      .map { meta, initial, netstart, transaid, ns_td, tai, blast, samba, all, unique, kept ->
          def line = "${meta.id}\t${initial}\t${netstart}\t${transaid}\t${ns_td}\t${tai}\t${blast}\t${samba}\t${all}\t${unique}\t${kept}"
          return line
      }
      .collectFile(
        name: 'counts.tsv',
        newLine: true,
        storeDir: "${output_dir}/samplesheets"
      )
      .set { ch_counts }

      ch_versions = ch_versions.mix(CONCAT.out.versions)
      ch_pipeline_versions = ch_versions
          .collectFile(
              name:      "HLorf.versions.yml",
              storeDir:  "${output_dir}/XORF_PIPELINE_INFO",
              sort:      true,
              keepHeader: false,
              newLine:   true
          )

    emit:
      files = CONCAT.out.files
      counts = ch_counts
      versions = ch_pipeline_versions
}
