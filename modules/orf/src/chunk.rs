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

use genepred::{bed::BedFormat, Bed12, GenePred, Gff, Gtf, Reader, ReaderResult, Strand, Writer};
use memchr::memchr;
use memmap2::Mmap;
use rayon::prelude::*;
use twobit::TwoBitFile;

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::{
    fs::{create_dir_all, File},
    io::{BufWriter, Write},
};

use crate::cli::ChunkArgs;

pub fn run_chunk(args: ChunkArgs) {
    let genome = get_sequences(args.sequence);
    let outdir = args.outdir.join("tmp");
    create_dir_all(&outdir).unwrap_or_else(|e| panic!("{}", e));

    match detect_region_format(&args.regions) {
        Some(RegionFormat::Bed) => process_reader::<Bed12>(
            &args.regions,
            args.chunks,
            &outdir,
            &genome,
            args.upstream_flank,
            args.downstream_flank,
            args.ignore_errors,
        ),
        Some(RegionFormat::Gtf) => process_reader::<Gtf>(
            &args.regions,
            args.chunks,
            &outdir,
            &genome,
            args.upstream_flank,
            args.downstream_flank,
            args.ignore_errors,
        ),
        Some(RegionFormat::Gff) => process_reader::<Gff>(
            &args.regions,
            args.chunks,
            &outdir,
            &genome,
            args.upstream_flank,
            args.downstream_flank,
            args.ignore_errors,
        ),
        None => panic!("ERROR: Unsupported file format"),
    }
}

fn process_reader<R>(
    regions: &Path,
    chunks: usize,
    outdir: &Path,
    genome: &HashMap<Vec<u8>, Vec<u8>>,
    upstream_flank: usize,
    downstream_flank: usize,
    ignore_errors: bool,
) where
    R: BedFormat + Into<GenePred> + Send,
{
    Reader::<R>::from_mmap(regions)
        .unwrap_or_else(|e| panic!("{}", e))
        .par_chunks(chunks)
        .unwrap_or_else(|e| panic!("{}", e))
        .for_each(|(idx, chunk)| {
            write_chunk(
                idx,
                chunk,
                genome,
                outdir,
                upstream_flank,
                downstream_flank,
                ignore_errors,
            )
        });
}

fn write_chunk(
    idx: usize,
    chunk: Vec<ReaderResult<GenePred>>,
    genome: &HashMap<Vec<u8>, Vec<u8>>,
    outdir: &Path,
    upstream_flank: usize,
    downstream_flank: usize,
    ignore_errors: bool,
) {
    let tmp = outdir.join(format!("tmp_{}.bed", idx));
    let mut writer = BufWriter::new(File::create(&tmp).unwrap_or_else(|e| panic!("{}", e)));

    let mut f_writer =
        BufWriter::new(File::create(tmp.with_extension("fa")).unwrap_or_else(|e| panic!("{}", e)));

    chunk
        .into_iter()
        .filter_map(|result| match result {
            Ok(record) => Some(record),
            Err(e) => {
                eprintln!("WARN: Failed to process record: {}", e);
                None
            }
        })
        .for_each(|record| {
            let seq = genome.get(&record.chrom).unwrap_or_else(|| {
                panic!(
                    "ERROR: Chromosome {} not found!",
                    std::str::from_utf8(&record.chrom).unwrap()
                )
            });

            let target = extract_seq(
                &record,
                seq,
                upstream_flank,
                downstream_flank,
                ignore_errors,
            );

            if let Some(mut target) = target {
                match &record.strand {
                    Some(Strand::Forward) => {}
                    Some(Strand::Reverse) => {
                        target.reverse();

                        for base in target.iter_mut() {
                            *base = match *base {
                                b'A' => b'T',
                                b'C' => b'G',
                                b'G' => b'C',
                                b'T' => b'A',
                                b'N' => b'N',
                                b'a' => b't',
                                b'c' => b'g',
                                b'g' => b'c',
                                b't' => b'a',
                                b'n' => b'n',
                                _ => panic!("ERROR: Invalid base"),
                            }
                        }
                    }
                    Some(Strand::Unknown) | None => {}
                }

                Writer::<Bed12>::from_record(&record, &mut writer)
                    .unwrap_or_else(|e| panic!("ERROR: Cannot write record to file: {}", e));

                f_writer.write_all(b">").unwrap_or_else(|e| panic!("{}", e));
                f_writer
                    .write_all(record.name().unwrap())
                    .unwrap_or_else(|e| panic!("{}", e));
                f_writer
                    .write_all(b"\n")
                    .unwrap_or_else(|e| panic!("{}", e));
                f_writer
                    .write_all(&target)
                    .unwrap_or_else(|e| panic!("{}", e));
                f_writer
                    .write_all(b"\n")
                    .unwrap_or_else(|e| panic!("{}", e));
            }
        });

    // INFO: flush first so data is written to the file
    writer.flush().unwrap_or_else(|e| panic!("{}", e));
    f_writer.flush().unwrap_or_else(|e| panic!("{}", e));

    // INFO: check actual file size via metadata
    let bed_path = tmp.clone();
    let file_len = std::fs::metadata(&bed_path).map(|m| m.len()).unwrap_or(0);
    if file_len == 0 {
        std::fs::remove_file(&bed_path).unwrap_or_else(|e| panic!("{}", e));
        std::fs::remove_file(bed_path.with_extension("fa")).unwrap_or_else(|e| panic!("{}", e));
    }
}

#[derive(Debug)]
#[allow(dead_code)]
enum RangeError {
    Underflow { feature_coord: usize, flank: usize },
}

impl std::fmt::Display for RangeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RangeError::Underflow {
                feature_coord,
                flank,
            } => write!(
                f,
                "ERROR: Feature coordinate {} is underflowing by {} bases",
                feature_coord, flank
            ),
        }
    }
}

fn slice_range_for_exon(
    exon_idx: usize,
    exon_count: usize,
    exon_start: usize,
    exon_end: usize,
    upstream_flank: usize,
    downstream_flank: usize,
) -> Result<std::ops::Range<usize>, RangeError> {
    let is_first = exon_idx == 0;
    let is_last = exon_idx + 1 == exon_count;

    let start = if is_first {
        exon_start
            .checked_sub(upstream_flank)
            .ok_or(RangeError::Underflow {
                feature_coord: exon_start,
                flank: upstream_flank,
            })?
    } else {
        exon_start
    };

    let end = if is_last {
        exon_end
            .checked_add(downstream_flank)
            .ok_or(RangeError::Underflow {
                feature_coord: exon_end,
                flank: downstream_flank,
            })?
    } else {
        exon_end
    };

    Ok(start..end)
}

fn extend_or_handle_oob(
    record: &GenePred,
    target: &mut Vec<u8>,
    seq: &[u8],
    range: std::ops::Range<usize>,
    exon_idx: usize,
    ignore_errors: bool,
) -> bool {
    if let Some(slice) = seq.get(range.clone()) {
        target.extend_from_slice(slice);
        return true;
    }

    if ignore_errors {
        eprintln!(
            "WARN: out-of-bounds slice for {} exon {}: {:?} (seq_len={})",
            record,
            exon_idx,
            range,
            seq.len()
        );
        false
    } else {
        panic!(
            "ERROR: out-of-bounds slice for {} exon {}: {:?} (seq_len={})",
            record,
            exon_idx,
            range,
            seq.len()
        );
    }
}

fn extract_seq(
    record: &GenePred,
    seq: &[u8],
    upstream_flank: usize,
    downstream_flank: usize,
    ignore_errors: bool,
) -> Option<Vec<u8>> {
    let exons = record.exons();
    let exon_count = exons.len();
    let mut target = Vec::new();

    for (idx, (start, end)) in exons.iter().enumerate() {
        let exon_start = *start as usize;
        let exon_end = *end as usize;

        let range = match slice_range_for_exon(
            idx,
            exon_count,
            exon_start,
            exon_end,
            upstream_flank,
            downstream_flank,
        ) {
            Ok(r) => r,
            Err(err) => {
                if ignore_errors {
                    eprintln!("WARN: {:?} for record {}", err, record);
                    return None;
                } else {
                    panic!("ERROR: {:?} for record {}", err, record);
                }
            }
        };

        if !extend_or_handle_oob(record, &mut target, seq, range, idx, ignore_errors) {
            return None;
        }
    }

    Some(target)
}

#[derive(Clone, Copy)]
enum RegionFormat {
    Bed,
    Gtf,
    Gff,
}

fn detect_region_format(path: &Path) -> Option<RegionFormat> {
    match path.extension().and_then(|ext| ext.to_str()) {
        Some("bed") => Some(RegionFormat::Bed),
        Some("gtf") => Some(RegionFormat::Gtf),
        Some("gff") => Some(RegionFormat::Gff),
        Some("gz") => {
            let stem = path.file_stem()?.to_str()?;
            if stem.ends_with(".bed") {
                Some(RegionFormat::Bed)
            } else if stem.ends_with(".gtf") {
                Some(RegionFormat::Gtf)
            } else if stem.ends_with(".gff") {
                Some(RegionFormat::Gff)
            } else {
                None
            }
        }
        _ => None,
    }
}

pub fn get_sequences(sequence: PathBuf) -> HashMap<Vec<u8>, Vec<u8>> {
    match sequence.extension() {
        Some(ext) => match ext.to_str() {
            Some("2bit") => from_2bit(sequence),
            Some("fa") | Some("gz") => from_fa(sequence),
            _ => panic!("ERROR: Unsupported file format"),
        },
        None => panic!("ERROR: No file extension"),
    }
}

fn from_2bit(twobit: PathBuf) -> HashMap<Vec<u8>, Vec<u8>> {
    let mut genome = TwoBitFile::open_and_read(twobit).expect("ERROR: Cannot open 2bit file");

    let mut sequences = HashMap::new();
    genome.chrom_names().iter().for_each(|chr| {
        let seq = genome
            .read_sequence(chr, ..)
            .unwrap_or_else(|e| panic!("ERROR: {}", e))
            .as_bytes()
            .to_vec();

        sequences.insert(chr.as_bytes().to_vec(), seq);
    });

    sequences
}

pub fn from_fa<F: AsRef<Path>>(f: F) -> HashMap<Vec<u8>, Vec<u8>> {
    let file = File::open(f).unwrap();
    let mmap = unsafe { Mmap::map(&file).unwrap() };
    let data = mmap.as_ref();

    let mut acc = HashMap::new();
    let mut pos = 0;

    while let Some(start) = memchr(b'>', &data[pos..]) {
        let start = pos + start;
        let end = memchr(b'>', &data[start + 1..]).map_or(data.len(), |e| start + 1 + e);
        let entry = &data[start + 1..end];
        let header_end = memchr(b'\n', entry).unwrap();
        let header = &entry[..header_end];
        let record = &entry[header_end + 1..];
        let seq = record
            .iter()
            .filter(|&&b| b != b'\n' && b != b'\r') // Remove newlines and carriage returns
            .cloned()
            .collect::<Vec<u8>>();

        acc.insert(header.to_vec(), seq);
        pos = end;
    }

    acc
}
