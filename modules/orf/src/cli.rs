use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(name = "orf", about = "Open reading frame pipeline wrapper", version = env!("CARGO_PKG_VERSION"))]
pub struct Args {
    #[command(subcommand)]
    pub command: Commands,

    #[arg(
        short = 't',
        long = "threads",
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    pub threads: usize,

    #[arg(
        short = 'L',
        long = "level",
        help = "Logging level",
        value_name = "LEVEL",
        default_value_t = log::Level::Info,
    )]
    pub level: log::Level,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    /// Run ORF detection using orfipy + diamond
    Blast(BlastArgs),

    /// Run ORF detection using TranslationAi
    Tai(TaiArgs),

    /// Read and merge TOGA results
    Chunk(ChunkArgs),

    /// Run RNAsamba on a set of reads
    Samba(SambaArgs),
}

#[derive(Debug, Parser)]
pub struct BlastArgs {
    #[arg(
        short = 'f',
        long = "fasta",
        required = true,
        help = "Path to .fa file"
    )]
    pub fasta: PathBuf,

    #[arg(
        short = 'b',
        long = "bed",
        required = true,
        help = "Path to .bed file with candidate regions"
    )]
    pub bed: PathBuf,

    #[arg(
        short = 'o',
        long = "outdir",
        required = false,
        help = "Path to outdir",
        default_value = "."
    )]
    pub outdir: PathBuf,

    #[arg(
        short = 'd',
        long = "database",
        required = true,
        help = "Path to protein database"
    )]
    pub database: PathBuf,

    #[arg(
        short = 't',
        long = "tai",
        required = false,
        help = "Path to translationAi output file"
    )]
    pub tai: Option<PathBuf>,

    #[arg(
        short = 'l',
        long = "orf-min-len",
        required = false,
        help = "Minimum length for ORFs",
        default_value = "50"
    )]
    pub orf_min_len: usize,

    #[arg(
        short = 'p',
        long = "orf-min-percent",
        required = false,
        help = "Minimum percentage of ORF length that is aligned",
        default_value = "0.25"
    )]
    pub orf_min_percent: f32,

    #[arg(
        short = 'u',
        long = "upstream-flank",
        required = false,
        help = "Number of bases upstream of the ORF to include in the context sequence",
        default_value = "0"
    )]
    pub upstream_flank: usize,

    #[arg(
        short = 'd',
        long = "downstream-flank",
        required = false,
        help = "Number of bases downstream of the ORF to include in the context sequence",
        default_value = "0"
    )]
    pub downstream_flank: usize,
}

#[derive(Debug, Parser)]
pub struct TaiArgs {
    #[arg(
        short = 'f',
        long = "fasta",
        required = true,
        help = "Path to .fa/.fa.gz"
    )]
    pub fasta: PathBuf,

    #[arg(
        short = 'b',
        long = "bed",
        required = true,
        help = "Path to .bed file with candidate regions"
    )]
    pub bed: PathBuf,

    #[arg(
        short = 'o',
        long = "outdir",
        required = false,
        help = "Path to outdir",
        default_value = "."
    )]
    pub outdir: PathBuf,

    #[arg(
        short = 'u',
        long = "upstream-flank",
        required = false,
        help = "Number of bases upstream of the TSS to include in the chunk",
        default_value = "0"
    )]
    pub upstream_flank: usize,

    #[arg(
        short = 'd',
        long = "downstream-flank",
        required = false,
        help = "Number of bases downstream of the TSS to include in the chunk",
        default_value = "0"
    )]
    pub downstream_flank: usize,
}

#[derive(Debug, Parser)]
pub struct ChunkArgs {
    #[arg(
        short = 's',
        long = "sequence",
        required = true,
        help = "Path to .fa/.fa.gz/.2bit"
    )]
    pub sequence: PathBuf,

    #[arg(
        short = 'r',
        long = "regions",
        required = true,
        help = "Path to .bed file with candidate regions"
    )]
    pub regions: PathBuf,

    #[arg(
        short = 'o',
        long = "outdir",
        required = false,
        help = "Path to outdir",
        default_value = "."
    )]
    pub outdir: PathBuf,

    #[arg(
        short = 'c',
        long = "chunks",
        required = false,
        help = "Number of chunks to split the bed file into",
        default_value = "10000"
    )]
    pub chunks: usize,

    #[arg(
        short = 'u',
        long = "upstream-flank",
        required = false,
        help = "Number of bases upstream of the TSS to include in the chunk",
        default_value = "0"
    )]
    pub upstream_flank: usize,

    #[arg(
        short = 'd',
        long = "downstream-flank",
        required = false,
        help = "Number of bases downstream of the TSS to include in the chunk",
        default_value = "0"
    )]
    pub downstream_flank: usize,

    #[arg(
        short = 'I',
        long = "ignore-errors",
        help = "Ignore slicing errors coming from out-of-bounds coordinates",
        action = clap::ArgAction::SetTrue
    )]
    pub ignore_errors: bool,
}

#[derive(Debug, Parser)]
pub struct SambaArgs {
    #[arg(
        short = 'f',
        long = "fasta",
        required = true,
        help = "Path to .fa file"
    )]
    pub fasta: PathBuf,

    #[arg(
        short = 'o',
        long = "outdir",
        default_value = ".",
        help = "Output directory"
    )]
    pub outdir: PathBuf,

    #[arg(
        short = 'w',
        long = "weights",
        required = false,
        help = "Path to model weights file"
    )]
    pub weights: Option<PathBuf>,

    #[arg(
        short = 'u',
        long = "upstream-flank",
        required = false,
        help = "Number of bases upstream of the TSS to include in the chunk",
        default_value = "0"
    )]
    pub upstream_flank: usize,

    #[arg(
        short = 'd',
        long = "downstream-flank",
        required = false,
        help = "Number of bases downstream of the TSS to include in the chunk",
        default_value = "0"
    )]
    pub downstream_flank: usize,
}
