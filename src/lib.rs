pub mod external_command_checker;
pub mod kmer_coverage;
pub mod genome_pseudoaligner;
pub mod pseudoaligner;
pub mod screen;
pub mod core_genome;
pub mod nucmer_runner;
pub mod nucmer_core_genome_generator;
pub mod ani_clustering;
pub mod pseudoalignment_reference_readers;
pub mod mapping_parameters;

// Modules below shared as direct code copy between coverm and cockatoo
pub mod genomes_and_contigs;
pub mod genome_parsing;

extern crate bio;
#[macro_use]
extern crate log;

extern crate env_logger;
extern crate tempdir;
extern crate tempfile;
extern crate rand;
extern crate debruijn;
extern crate boomphf;
#[macro_use]
extern crate lazy_static;
extern crate rayon;
extern crate failure;
extern crate crossbeam;
extern crate flate2;
extern crate bincode;
#[macro_use]
extern crate serde;
extern crate csv;
extern crate rstar;
extern crate finch;
extern crate seq_io;

use std::io::Read;

fn finish_command_safely(
    mut process: std::process::Child, process_name: &str)
-> std::process::Child {
    let es = process.wait()
        .expect(&format!("Failed to glean exitstatus from failing {} process", process_name));
    debug!("Process {} finished", process_name);
    if !es.success() {
        error!("Error when running {} process.", process_name);
        let mut err = String::new();
        process.stderr.expect(&format!("Failed to grab stderr from failed {} process", process_name))
            .read_to_string(&mut err).expect("Failed to read stderr into string");
        error!("The STDERR was: {:?}", err);
        let mut out = String::new();
        process.stdout.expect(&format!("Failed to grab stdout from failed {} process", process_name))
            .read_to_string(&mut out).expect("Failed to read stdout into string");
        error!("The STDOUT was: {:?}", out);
        error!("Cannot continue after {} failed.", process_name);
        std::process::exit(1);
    }
    return process;
}

fn run_command_safely(
    mut cmd: std::process::Command,
    process_name: &str)
    -> std::process::Child {

    let process = cmd.spawn().expect(&format!("Failed to spawn {}", process_name));
    return finish_command_safely(process, process_name);
}

