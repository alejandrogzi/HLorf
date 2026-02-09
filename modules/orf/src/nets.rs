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

use genepred::GenePred;

use crate::{cli::NetArgs, consts::*, utils::*};

/// Netstart holding struct
pub struct NetNS {
    // information_dict["origin"] = []
    // information_dict["atg_position"] = []
    // information_dict["stop_codon_position"] = []
    // information_dict["peptide_len"] = []
    // information_dict["entry_line"] = []
    // information_dict["strand"] = []
    // information_dict["preds"] = []
}

impl NetNS {}

/// Transaid holding struct
pub struct NetTD {
    // "Sequence_ID,Start,Stop,TIS_Score,TTS_Score,Kozak_Score,"
    // "CAI_Score,GC_Score,Integrated_Score,Protein_Length,"
    // "Passed_Filter,Filter_Reason\n"
}

impl NetTD {}

pub fn run_nets(args: NetArgs) {}
