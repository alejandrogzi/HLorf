use clap::{self, Parser};
use log::info;
use simple_logger::init_with_level;

use orf::{
    blast::run_blast,
    cli::{Args, Commands},
    samba::run_samba,
    tai::run_tai,
};

fn main() {
    let start = std::time::Instant::now();

    let args = Args::parse();

    init_with_level(args.level).unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .unwrap();

    match args.command {
        Commands::Blast(args) => run_blast(args),
        Commands::Tai(args) => run_tai(args),
        Commands::Samba(args) => run_samba(args),
        Commands::Chunk(args) => run_chunk(args),
    }

    let elapsed = start.elapsed();
    info!("Elapsed time: {:.3?}", elapsed);
}
