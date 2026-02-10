//! Core module for detecting open reading frames in a query set of reads
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main functions for finding open reading frames (ORFs)
//! in a set of aligned reads.
//!
//! In short, every possible open reading frame (ORF) is detected for every
//! read in the query set. For every potential ORF, learning models and databases
//! are used to determine whether the ORF is a true ORF, a false positive.
//! All the data from each reliable ORF is collected and subjected to another
//! learning model trained with true ORFs and false positives. The process is
//! heavily parallelized to offer fast performance on large datasets.

use std::io::Write;

use crate::{
    cli::SambaArgs,
    utils::{create_fasta, parse_fa},
};

const WEIGHTS: &str =
    "https://github.com/apcamargo/RNAsamba/raw/refs/heads/master/data/full_length_weights.hdf5";

pub fn run_samba(args: SambaArgs) {
    let dir = &args.outdir.join("samba");
    std::fs::create_dir_all(dir)
        .unwrap_or_else(|e| panic!("ERROR: could not create directory -> {e}!"));

    let weights = args.weights.unwrap_or_else(|| {
        let weights = dir.join("full_length_weights.hdf5");

        log::info!("INFO: downloading weights from {WEIGHTS}");
        let cmd = format!(
            "wget {WEIGHTS} -O {weights}",
            WEIGHTS = WEIGHTS,
            weights = weights.display()
        );

        std::process::Command::new("bash")
            .arg("-c")
            .arg(&cmd)
            .status()
            .unwrap_or_else(|e| panic!("ERROR: could not run command -> {cmd} -> {e}!"));

        weights
    });

    // INFO: if flanks activated, we need to correct the fasta
    let fasta = if args.upstream_flank > 0 || args.downstream_flank > 0 {
        let seqs = parse_fa(&args.fasta).unwrap_or_else(|e| {
            panic!("ERROR: failed to parse FASTA file -> {e}");
        });
        let mut writer = create_fasta(&args.fasta, "tmp.strip.fa")
            .unwrap_or_else(|| panic!("ERROR: failed to create tmp file -> {:?}", args.fasta));

        for (header, seq) in seqs.into_iter() {
            if seq.is_empty() {
                println!("WARN: empty sequence for header -> {header}");
                continue;
            }

            let stripped = seq
                .get(args.upstream_flank..seq.len() - args.downstream_flank)
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: failed to get sequence for {:?} from slice -> {:?}",
                        header, seq
                    )
                });

            writer
                .write_all(format!(">{}\n", header).as_bytes())
                .unwrap_or_else(|e| {
                    panic!("ERROR: failed to write to file -> {e}");
                });
            writer.write_all(stripped).unwrap_or_else(|e| {
                panic!("ERROR: failed to write to file -> {e}");
            });
            writer.write_all(b"\n").unwrap_or_else(|e| {
                panic!("ERROR: failed to write to file -> {e}");
            });
        }

        args.fasta.with_extension("tmp.strip.fa")
    } else {
        args.fasta.clone()
    };

    let output = dir.join(
        args.fasta
            .with_extension("tmp.samba.tsv")
            .file_name()
            .unwrap(),
    );

    let cmd = format!(
        "rnasamba classify {output} {input} {weights}",
        output = output.display(),
        input = fasta.display(),
        weights = weights.display()
    );

    log::info!("INFO: running command -> {cmd}");
    let status = std::process::Command::new("bash")
        .arg("-c")
        .arg(&cmd)
        .status()
        .unwrap_or_else(|e| panic!("ERROR: could not run command -> {cmd} -> {e}!"));

    if !status.success() {
        panic!("ERROR: command failed -> {cmd}");
    }
}
