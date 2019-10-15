extern crate cockatoo;

extern crate coverm;
//use coverm::CONCATENATED_FASTA_FILE_SEPARATOR;
use coverm::mapping_parameters::MappingParameters;
use cockatoo::external_command_checker; // TODO: check for nucmer
use std::env;
use std::str;
use std::process;
use std::io::Write;

extern crate clap;
use clap::*;

#[macro_use]
extern crate log;
extern crate env_logger;
use log::LevelFilter;
use env_logger::Builder;

extern crate rayon;
use rayon::prelude::*;

fn main(){
    let mut app = build_cli();
    let matches = app.clone().get_matches();
    set_log_level(&matches, false);

    match matches.subcommand_name() {

        Some("genome") => {
            let m = matches.subcommand_matches("genome").unwrap();
            set_log_level(m, true);
            let pseudoalign_params = parse_pseudoaligner_parameters(&m);

            info!("Restoring index ..");
            let core_genome_pseudoaligner = cockatoo::core_genome::restore_index::<cockatoo::pseudoaligner::config::KmerType>(
                m.value_of("index").unwrap());

            info!("Aligning reads and post-processing ..");
            cockatoo::genome_pseudoaligner::core_genome_coverage_pipeline(
                &pseudoalign_params.reads,
                pseudoalign_params.num_threads,
                !m.is_present("no-zeros"),
                &core_genome_pseudoaligner,
            );
        },
        Some("index") => {
            let m = matches.subcommand_matches("index").unwrap();
            set_log_level(m, true);

            let reference = m.value_of("reference").unwrap();
            let output = m.value_of("index").unwrap();
            let num_threads = value_t!(m.value_of("threads"), usize).unwrap();

            rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .build_global()
                .expect("Programming error: rayon initialised multiple times");

            let genomes_and_contigs = if m.is_present("genome-definition") {
                let definition_file = m.value_of("genome-definition").unwrap();
                info!("Reading contig names associated with each genome from {}..", definition_file);
                coverm::read_genome_definition_file(&definition_file)
            } else {
                let genome_fasta_files: Vec<String> = parse_list_of_genome_fasta_files(m);
                info!("Reading contig names for {} genomes ..", genome_fasta_files.len());
                coverm::read_genome_fasta_files(
                    &genome_fasta_files.iter().map(|s| s.as_str()).collect())
            };

            info!("Generating index based on info from individual genome files ..");
            let index = cockatoo::pseudoalignment_reference_readers::generate_debruijn_index_grouping_via_genomes_and_contigs::
            <cockatoo::pseudoaligner::config::KmerType>(
                &genomes_and_contigs,
                reference,
                num_threads
            );

            // For each clade, nucmer against the first genome.
            info!("Reading clade definition ..");
            let clade_definitions_file = m.value_of("clades").unwrap();
            let clades = cockatoo::genome_pseudoaligner::read_clade_definition_file(
                clade_definitions_file);

            // TODO: ProgressBar?
            info!("Calculating core genomes ..");
            let nucmer_core_genomes: Vec<Vec<Vec<cockatoo::core_genome::CoreGenomicRegion>>> = clades
                .par_iter()
                .enumerate()
                .map(
                    |(i, clade_fastas)|
                    cockatoo::nucmer_core_genome_generator::nucmer_core_genomes_from_genome_fasta_files(
                        &clade_fastas.iter().map(|s| &**s).collect::<Vec<&str>>()[..],
                        i as u32)
                ).collect();
            info!("Finished calculating core genomes");

            // Check core genome sizes / report
            cockatoo::genome_pseudoaligner::report_core_genome_sizes(
                &nucmer_core_genomes, &clades);

            // Thread genomes recording the core genome nodes
            // TODO: These data are at least sometimes read in repeatedly, when they
            // maybe should just be cached or something.
            info!("Reading in genome FASTA files to thread graph");
            let dna_strings = cockatoo::genome_pseudoaligner::read_clade_genome_strings(
                &clades);
            info!("Threading DeBruijn graph");
            // TODO: Multithread?
            let core_genome_pseudoaligner = cockatoo::core_genome::generate_core_genome_pseudoaligner(
                &nucmer_core_genomes,
                &dna_strings,
                index,
                genomes_and_contigs,
            );
            debug!("Found index {:?}", core_genome_pseudoaligner);

            info!("Saving index ..");
            cockatoo::core_genome::save_index(
                core_genome_pseudoaligner, &output);
            info!("Saving complete");
        },

        Some("index_stats") => {
            let m = matches.subcommand_matches("index_stats").unwrap();
            set_log_level(m, true);

            let index_path = m.value_of("index").unwrap();
            info!("Reading index {} ..", index_path);
            let index = cockatoo::core_genome::restore_index::<cockatoo::pseudoaligner::config::KmerType>(index_path);
            info!("Finished reading index");

            // Write GFA
            if m.is_present("gfa") { 
                let gfa_filename = m.value_of("gfa").unwrap();
                info!("Writing GFA file {} ..", gfa_filename);
                let mut gfa_writer = std::fs::File::create(gfa_filename).unwrap();
                index
                    .index
                    .dbg
                    .write_gfa(&mut gfa_writer)
                    .expect("Failed to write GFA file");

                // Write CSV data to be loaded into bandage
                let gfa_csv_filename = format!("{}.core_nodes.csv", gfa_filename);
                info!("Writing Bandage-compatible core node annotations to {} ..",
                    gfa_csv_filename);
                let mut csv_writer = std::fs::File::create(gfa_csv_filename).unwrap();
                writeln!(csv_writer, "Node,Clades").unwrap();
                for (node_id, clades) in &index.node_id_to_clade_cores {
                    writeln!(csv_writer, "{},\"{:?}\"", node_id, clades).unwrap();
                }
            }

            // Write info table for each clade
            println!("clade_id\tcore_genome_size");
            for (clade_i, size) in index.core_genome_sizes.iter().enumerate() {
                println!("{}\t{}", clade_i, size);
            }
        },

        Some("cluster") => {
            let m = matches.subcommand_matches("cluster").unwrap();
            set_log_level(m, true);

            let num_threads = value_t!(m.value_of("threads"), usize).unwrap();
            rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .build_global()
                .expect("Programming error: rayon initialised multiple times");

            let genome_fasta_files: Vec<String> = parse_list_of_genome_fasta_files(m);
            let v2: Vec<&str> = genome_fasta_files.iter().map(|s| &**s).collect();
            info!("Clustering {} genomes ..", genome_fasta_files.len());

            let clusters = cockatoo::ani_clustering::minhash_clusterer::minhash_clusters(
                &v2,
                value_t!(m.value_of("ani"), f32).unwrap(),
            );
            info!("Found {} genome clusters", clusters.len());

            for cluster in clusters {
                let rep_index = cluster[0];
                for genome_index in cluster {
                    println!("{}\t{}", v2[rep_index], v2[genome_index]);
                }
            }
            info!("Finished printing genome clusters");
        },
        Some("screen") => {
            let m = matches.subcommand_matches("screen").unwrap();
            set_log_level(m, true);

            // Not really mapping parameters but gets the job done
            // TODO: Remove coverm from this and do our own thing. That'd cut out the rust-htslib too.
            let read_inputs = MappingParameters::generate_from_clap(
                &m, coverm::bam_generator::MappingProgram::BWA_MEM, &None);
            let num_threads = value_t!(m.value_of("threads"), usize).unwrap();

            external_command_checker::check_for_mash();

            cockatoo::screen::screen_and_print_matching_genomes(
                m.value_of("mash-screen-db").unwrap(),
                &read_inputs,
                num_threads,
                parse_percentage(&m, "min-identity"),
            );
        }
        _ => {
            app.print_help().unwrap();
            println!();
        }
    }
}

fn parse_percentage(m: &clap::ArgMatches, parameter: &str) -> f32 {
    match m.is_present(parameter) {
        true => {
            let mut percentage = value_t!(
                m.value_of(parameter), f32).unwrap();
            if percentage >= 1.0 && percentage <= 100.0 {
                percentage = percentage / 100.0;
            } else if percentage < 0.0 || percentage > 100.0 {
                error!("Invalid alignment percentage: '{}'", percentage);
                process::exit(1);
            }
            info!("Using {} {}%", parameter, percentage*100.0);
            percentage
        },
        false => 0.0
    }
}

// fn parse_separator(m: &clap::ArgMatches) -> Option<u8> {
//     let single_genome = m.is_present("single-genome");
//     if single_genome {
//         Some("0".as_bytes()[0])
//     } else if m.is_present("separator") {
//         let separator_str = m.value_of("separator").unwrap().as_bytes();
//         if separator_str.len() != 1 {
//             eprintln!(
//                 "error: Separator can only be a single character, found {} ({}).",
//                 separator_str.len(),
//                 str::from_utf8(separator_str).unwrap());
//             process::exit(1);
//         }
//         Some(separator_str[0])
//     } else if m.is_present("bam-files") || m.is_present("reference") {
//         // Argument parsing enforces that genomes have been specified as FASTA
//         // files.
//         None
//     } else {
//         // Separator is set by CoverM and written into the generated reference
//         // fasta file.
//         Some(CONCATENATED_FASTA_FILE_SEPARATOR.as_bytes()[0])
//     }
// }


fn parse_list_of_genome_fasta_files(m: &clap::ArgMatches) -> Vec<String> {
    match m.is_present("genome-fasta-files") {
        true => {
            m.values_of("genome-fasta-files").unwrap().map(|s| s.to_string()).collect()
        },
        false => {
            if m.is_present("genome-fasta-directory") {
                let dir = m.value_of("genome-fasta-directory").unwrap();
                let paths = std::fs::read_dir(dir).unwrap();
                let mut genome_fasta_files: Vec<String> = vec!();
                let extension = m.value_of("genome-fasta-extension").unwrap();
                for path in paths {
                    let file = path.unwrap().path();
                    match file.extension() {
                        Some(ext) => {
                            if ext == extension {
                                let s = String::from(file.to_string_lossy());
                                genome_fasta_files.push(s);
                            } else {
                                info!(
                                    "Not using directory entry '{}' as a genome FASTA file, as \
                                     it does not end with the extension '{}'",
                                    file.to_str().expect("UTF8 error in filename"),
                                    extension);
                            }
                        },
                        None => {
                            info!("Not using directory entry '{}' as a genome FASTA file",
                                  file.to_str().expect("UTF8 error in filename"));
                        }
                    }
                }
                if genome_fasta_files.len() == 0 {
                    error!("Found 0 genomes from the genome-fasta-directory, cannot continue.");
                    process::exit(1);
                }
                genome_fasta_files // return
            } else {
                error!("Either a separator (-s) or path(s) to genome FASTA files \
                        (with -d or -f) must be given");
                process::exit(1);
            }
        }
    }
}




struct PseudoAlignmentParameters {
    pub num_threads: usize,
    pub reads: Vec<cockatoo::kmer_coverage::PseudoalignmentReadInput>,
}

fn parse_pseudoaligner_parameters(
    m: &clap::ArgMatches)
    -> PseudoAlignmentParameters {

    let num_threads = value_t!(m.value_of("threads"), usize).unwrap();
    let index = m.value_of("index").unwrap();

    let mapping_parameters = MappingParameters::generate_from_clap(
        &m, coverm::bam_generator::MappingProgram::BWA_MEM, &None);

    let mut pseudoalignment_read_input = vec!();
    // TODO: Accept interleaved output here, in genome and in screen
    for readset in mapping_parameters.readsets() {
        pseudoalignment_read_input.push(cockatoo::kmer_coverage::PseudoalignmentReadInput {
            forward_fastq: readset.0.to_string(),
            reverse_fastq: match readset.1 {
                Some(r2) => Some(r2.to_string()),
                None => None
            },
            sample_name: format!("{}/{}", index, readset.0)
        })

    }

    return PseudoAlignmentParameters {
        num_threads: num_threads,
        reads: pseudoalignment_read_input,
    }
}


fn set_log_level(matches: &clap::ArgMatches, is_last: bool) {
    let mut log_level = LevelFilter::Info;
    let mut specified = false;
    if matches.is_present("verbose") {
        specified = true;
        log_level = LevelFilter::Debug;
    }
    if matches.is_present("quiet") {
        specified = true;
        log_level = LevelFilter::Error;
    }
    if specified || is_last {
        let mut builder = Builder::new();
        builder.filter_level(log_level);
        if env::var("RUST_LOG").is_ok() {
            builder.parse_filters(&env::var("RUST_LOG").unwrap());
        }
        if builder.try_init().is_err() {
            panic!("Failed to set log level - has it been specified multiple times?")
        }
    }
    if is_last {
        info!("CoverM version {}", crate_version!());
    }
}

fn build_cli() -> App<'static, 'static> {
    return App::new("coverm")
        .version(crate_version!())
        .author("Ben J. Woodcroft <benjwoodcroft near gmail.com>")
        .about("Mapping coverage analysis for metagenomics")
        .args_from_usage("-v, --verbose       'Print extra debug logging information'
             -q, --quiet         'Unless there is an error, do not print logging information'")
        .global_setting(AppSettings::ArgRequiredElseHelp)
        .subcommand(
            SubCommand::with_name("genome")
                .about("Calculate coverage of genomes")

                .arg(Arg::with_name("index")
                    .short("-d")
                    .long("index")
                    .required(true)
                    .takes_value(true))

                .arg(Arg::with_name("read1")
                     .short("-1")
                     .multiple(true)
                     .takes_value(true)
                     .requires("read2")
                     .required_unless_one(
                         &["coupled","interleaved","single"]))
                .arg(Arg::with_name("read2")
                     .short("-2")
                     .multiple(true)
                     .takes_value(true)
                     .requires("read1")
                     .required_unless_one(
                         &["coupled","interleaved","single"]))
                .arg(Arg::with_name("coupled")
                     .short("-c")
                     .long("coupled")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless_one(
                         &["read1","interleaved","single"]))
                .arg(Arg::with_name("interleaved")
                     .long("interleaved")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless_one(
                         &["read1","coupled","single"]))
                .arg(Arg::with_name("single")
                     .long("single")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless_one(
                         &["read1","coupled","interleaved"])
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("threads")
                     .short("-t")
                     .long("threads")
                     .default_value("1")
                     .takes_value(true))

                .arg(Arg::with_name("no-zeros")
                     .long("no-zeros"))

                .arg(Arg::with_name("verbose")
                     .short("v")
                     .long("verbose"))
                .arg(Arg::with_name("quiet")
                     .short("q")
                     .long("quiet")))

        .subcommand(
            SubCommand::with_name("index")
                .about("Generate a mapping index for a collection of genomes")
                .arg(Arg::with_name("reference")
                     .short("-r")
                     .long("reference")
                     .takes_value(true)
                     .required(true))
                .arg(Arg::with_name("index")
                     .short("-d")
                     .long("index")
                     .takes_value(true)
                     .required(true))
                .arg(Arg::with_name("threads")
                     .short("-t")
                     .long("threads")
                     .default_value("1")
                     .takes_value(true))

                .arg(Arg::with_name("clades")
                     .long("clades")
                     .takes_value(true)
                     .required(true))

                .arg(Arg::with_name("genome-fasta-files")
                     .short("f")
                     .long("genome-fasta-files")
                     .multiple(true)
                     .conflicts_with("genome-fasta-directory")
                     .required_unless_one(
                         &["genome-fasta-directory","separator","genome-definition"])
                     .takes_value(true))
                .arg(Arg::with_name("genome-fasta-directory")
                     .long("genome-fasta-directory")
                     .conflicts_with("separator")
                     .conflicts_with("genome-fasta-files")
                     .required_unless_one(
                         &["genome-fasta-files","separator","genome-definition"])
                     .takes_value(true))
                .arg(Arg::with_name("genome-fasta-extension")
                     .short("x")
                     .long("genome-fasta-extension")
                     // Unsure why, but uncommenting causes test failure (in
                     // genome mode, not sure about here) - clap bug?
                     //.requires("genome-fasta-directory")
                     .default_value("fna")
                     .takes_value(true))                   
                .arg(Arg::with_name("separator")
                     .short("s")
                     .long("separator")
                     .required_unless_one(
                         &["genome-fasta-directory","genome-fasta-files","genome-definition"])
                     .takes_value(true))
                .arg(
                    Arg::with_name("genome-definition")
                        .long("genome-definition")
                        .conflicts_with("separator")
                        .conflicts_with("genome-fasta-files")
                        .conflicts_with("genome-fasta-directory")
                        .required_unless_one(&[
                            "genome-fasta-files",
                            "separator",
                            "genome-fasta-directory",
                        ])
                        .takes_value(true))

                .arg(Arg::with_name("verbose")
                     .short("v")
                     .long("verbose"))
                .arg(Arg::with_name("quiet")
                     .short("q")
                     .long("quiet")))

         .subcommand(
            SubCommand::with_name("index_stats")
                .about("Print information about a generated index")
                .arg(Arg::with_name("index")
                     .short("-d")
                     .long("index")
                     .takes_value(true)
                     .required(true))
                
                .arg(Arg::with_name("gfa")
                    .long("gfa")
                    .takes_value(true))

                .arg(Arg::with_name("verbose")
                     .short("v")
                     .long("verbose"))
                .arg(Arg::with_name("quiet")
                     .short("q")
                     .long("quiet")))

        .subcommand(
            SubCommand::with_name("cluster")
                .about("Cluster FASTA files by average nucleotide identity")
                .arg(Arg::with_name("ani")
                     .long("ani")
                     .takes_value(true)
                     .required(true))
                .arg(Arg::with_name("genome-fasta-files")
                     .short("f")
                     .long("genome-fasta-files")
                     .multiple(true)
                     .conflicts_with("genome-fasta-directory")
                     .conflicts_with("single-genome")
                     .required_unless_one(
                         &["genome-fasta-directory"])
                     .takes_value(true))
                .arg(Arg::with_name("genome-fasta-directory")
                     .long("genome-fasta-directory")
                     .conflicts_with("separator")
                     .conflicts_with("genome-fasta-files")
                     .conflicts_with("single-genome")
                     .required_unless_one(
                         &["genome-fasta-files"])
                     .takes_value(true))
                .arg(Arg::with_name("genome-fasta-extension")
                     .short("x")
                     .long("genome-fasta-extension")
                     // Unsure why, but uncommenting causes test failure (in
                     // genome mode, not sure about here) - clap bug?
                     //.requires("genome-fasta-directory")
                     .default_value("fna")
                     .takes_value(true))

                .arg(Arg::with_name("threads")
                     .short("-t")
                     .long("threads")
                     .default_value("1")
                     .takes_value(true)))

        .subcommand(
            SubCommand::with_name("screen")
                .about("Screen a genome set to determine presence / absence")
                .arg(Arg::with_name("mash-screen-db")
                     .short("-d")
                     .long("mash-screen-db")
                     .required(true)
                     .takes_value(true))
                .arg(Arg::with_name("min-identity")
                     .short("-m")
                     .long("min-identity")
                     .default_value("0.02")
                     .takes_value(true))
                .arg(Arg::with_name("read1")
                     .short("-1")
                     .multiple(true)
                     .takes_value(true)
                     .requires("read2")
                     .required_unless_one(
                         &["coupled","interleaved","single"]))
                .arg(Arg::with_name("read2")
                     .short("-2")
                     .multiple(true)
                     .takes_value(true)
                     .requires("read1")
                     .required_unless_one(
                         &["coupled","interleaved","single"]))
                .arg(Arg::with_name("coupled")
                     .short("-c")
                     .long("coupled")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless_one(
                         &["read1","interleaved","single"]))
                .arg(Arg::with_name("interleaved")
                     .long("interleaved")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless_one(
                         &["read1","coupled","single"]))
                .arg(Arg::with_name("single")
                     .long("single")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless_one(
                         &["read1","coupled","interleaved"]))

                .arg(Arg::with_name("threads")
                     .short("-t")
                     .long("threads")
                     .default_value("1")
                     .takes_value(true))
                .arg(Arg::with_name("verbose")
                     .short("v")
                     .long("verbose"))
                .arg(Arg::with_name("quiet")
                     .short("q")
                     .long("quiet")));
}
