use clap::Parser;
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(name = "tai", about = "In-house translationAi caller", version = env!("CARGO_PKG_VERSION"))]
pub struct Args {
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
