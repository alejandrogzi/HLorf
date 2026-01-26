use clap::Parser;
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(name = "orf", about = "Open reading frame pipeline wrapper", version = env!("CARGO_PKG_VERSION"))]
pub struct Args {
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
