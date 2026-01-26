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

pub const PREDICTIONS: &str = "predictions";
pub const TAI_VENV: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tai/.venv/bin/activate");

pub const SCALE: u64 = 100_000_000_000; // 100Gb

pub const RESULT: &str = "result";
pub const VENV: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/.venv/bin/activate");
pub const ORF_PEP: &str = "orf.pep";

pub const START_CODON: &str = "ATG";
pub const STOP_CODONS: [&str; 3] = ["TAA", "TAG", "TGA"];
pub const STOP_CODONS_BYTES: [&[u8]; 3] = [b"TAA", b"TAG", b"TGA"];

pub const CODON_TABLE: [(&[u8], &str); 64] = [
    // Phenylalanine (F)
    (b"TTT", "F"),
    (b"TTC", "F"),
    // Leucine (L)
    (b"TTA", "L"),
    (b"TTG", "L"),
    (b"CTA", "L"),
    (b"CTC", "L"),
    (b"CTG", "L"),
    (b"CTT", "L"),
    // Isoleucine (I)
    (b"ATT", "I"),
    (b"ATC", "I"),
    (b"ATA", "I"),
    // Methionine (M) - Start codon
    (b"ATG", "M"),
    // Valine (V)
    (b"GTA", "V"),
    (b"GTC", "V"),
    (b"GTG", "V"),
    (b"GTT", "V"),
    // Serine (S)
    (b"TCA", "S"),
    (b"TCC", "S"),
    (b"TCG", "S"),
    (b"TCT", "S"),
    (b"AGT", "S"),
    (b"AGC", "S"),
    // Proline (P)
    (b"CCA", "P"),
    (b"CCC", "P"),
    (b"CCG", "P"),
    (b"CCT", "P"),
    // Threonine (T)
    (b"ACA", "T"),
    (b"ACC", "T"),
    (b"ACG", "T"),
    (b"ACT", "T"),
    // Alanine (A)
    (b"GCA", "A"),
    (b"GCC", "A"),
    (b"GCG", "A"),
    (b"GCT", "A"),
    // Tyrosine (Y)
    (b"TAT", "Y"),
    (b"TAC", "Y"),
    // Stop codons (*)
    (b"TAA", "*"),
    (b"TAG", "*"),
    (b"TGA", "*"),
    // Histidine (H)
    (b"CAT", "H"),
    (b"CAC", "H"),
    // Glutamine (Q)
    (b"CAA", "Q"),
    (b"CAG", "Q"),
    // Asparagine (N)
    (b"AAT", "N"),
    (b"AAC", "N"),
    // Lysine (K)
    (b"AAA", "K"),
    (b"AAG", "K"),
    // Aspartic acid (D)
    (b"GAT", "D"),
    (b"GAC", "D"),
    // Glutamic acid (E)
    (b"GAA", "E"),
    (b"GAG", "E"),
    // Cysteine (C)
    (b"TGT", "C"),
    (b"TGC", "C"),
    // Tryptophan (W)
    (b"TGG", "W"),
    // Arginine (R)
    (b"CGA", "R"),
    (b"CGC", "R"),
    (b"CGG", "R"),
    (b"CGT", "R"),
    (b"AGA", "R"),
    (b"AGG", "R"),
    // Glycine (G)
    (b"GGA", "G"),
    (b"GGC", "G"),
    (b"GGG", "G"),
    (b"GGT", "G"),
];
