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

fn run_chunk(args: ChunkArgs) {
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
        ),
        Some(RegionFormat::Gtf) => process_reader::<Gtf>(
            &args.regions,
            args.chunks,
            &outdir,
            &genome,
            args.upstream_flank,
            args.downstream_flank,
        ),
        Some(RegionFormat::Gff) => process_reader::<Gff>(
            &args.regions,
            args.chunks,
            &outdir,
            &genome,
            args.upstream_flank,
            args.downstream_flank,
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
) where
    R: BedFormat + Into<GenePred> + Send,
{
    Reader::<R>::from_mmap(regions)
        .unwrap_or_else(|e| panic!("{}", e))
        .par_chunks(chunks)
        .unwrap_or_else(|e| panic!("{}", e))
        .for_each(|(idx, chunk)| {
            write_chunk(idx, chunk, genome, outdir, upstream_flank, downstream_flank)
        });
}

fn write_chunk(
    idx: usize,
    chunk: Vec<ReaderResult<GenePred>>,
    genome: &HashMap<Vec<u8>, Vec<u8>>,
    outdir: &Path,
    upstream_flank: usize,
    downstream_flank: usize,
) {
    let tmp = outdir.join(format!("tmp_{}.bed", idx));
    let mut writer = BufWriter::new(File::create(&tmp).unwrap_or_else(|e| panic!("{}", e)));

    let mut f_writer =
        BufWriter::new(File::create(tmp.with_extension("fa")).unwrap_or_else(|e| panic!("{}", e)));

    chunk.into_iter().filter_map(Result::ok).for_each(|record| {
        // record.set_start(record.start() - upstream_flank as u64);
        // record.set_end(record.end() + downstream_flank as u64);

        let seq = genome.get(&record.chrom).unwrap_or_else(|| {
            panic!(
                "ERROR: Chromosome {} not found!",
                std::str::from_utf8(&record.chrom).unwrap()
            )
        });

        let mut target: Vec<u8> = vec![];
        for (idx, feature) in record.exons().iter().enumerate() {
            let feature_start = feature.0 as usize;
            let feature_end = feature.1 as usize;

            if idx == 0 {
                if record.exon_count() == 1 {
                    target.extend_from_slice(
                        &seq[feature_start - upstream_flank..feature_end + downstream_flank],
                    );
                } else {
                    target.extend_from_slice(&seq[feature_start - upstream_flank..feature_start]);
                }
            } else if idx == record.exons().len() - 1 {
                target.extend_from_slice(&seq[feature_start..feature_end + downstream_flank]);
            } else {
                target.extend_from_slice(&seq[feature_start..feature_end]);
            }
        }

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
    });
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

pub fn get_sequences<'a>(sequence: PathBuf) -> HashMap<Vec<u8>, Vec<u8>> {
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
