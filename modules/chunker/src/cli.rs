use clap::Parser;
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(name = "chunker", about = "In-house .bed/.gtf sequence chunker", version = env!("CARGO_PKG_VERSION"))]
pub struct Args {
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
}
