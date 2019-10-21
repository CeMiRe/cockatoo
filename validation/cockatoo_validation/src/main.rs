extern crate seq_io;
use seq_io::fasta::*;
use std::io;
use std::collections::HashSet;

#[macro_use] extern crate log;
extern crate env_logger;
extern crate clap;
use clap::{App, SubCommand};

extern crate rand;
use rand::prelude::*;

fn main() {
    let matches = App::new("cockatoo-validation")
        .subcommand(SubCommand::with_name("chop")
            .about("Remove the last 20% of each contig on stdin"))
        .subcommand(SubCommand::with_name("mutate")
            .about("Mutate 1% of a genome"))
        .get_matches();

    env_logger::init();
    info!("started");

    if let Some(_matches) = matches.subcommand_matches("chop") {
        let mut reader = Reader::new(io::stdin());

        while let Some(record) = reader.next() {
            let record = record.expect("Error reading record");
            println!(">{}",record.id().unwrap());
            let seq = record.seq();
            let bound = (seq.len() as f32 *0.8) as usize;
            
            println!("{}", std::str::from_utf8(&seq[0..bound]).unwrap())
        }
    } else if let Some(_matches) = matches.subcommand_matches("mutate") {
        let mut reader = Reader::new(io::stdin());

        while let Some(record) = reader.next() {
            let record = record.expect("Error reading record");
            println!(">{}",record.id().unwrap());
            let seq = record.seq();

            let mut rng = thread_rng();
            let to_mutate: Vec<usize> = (0..seq.len()).collect();
            let to_mutate_set: HashSet<usize> = HashSet::from(to_mutate.choose_multiple(
                &mut rng, 
                ((seq.len() as f32)*0.01) as usize).map(|x| *x).collect());
            
            let choices = [b'A',b'T',b'G',b'C'];
            for (i,c) in record.seq().iter().enumerate() {
                if to_mutate_set.contains(&i) {
                    let mut new_one = c;
                    while new_one==c {
                        new_one = choices.choose(&mut rng).unwrap();
                    }
                    print!("{}",String::from_utf8(vec![*new_one]).unwrap());
                } else {
                    print!("{}",String::from_utf8(vec![*c]).unwrap());
                }
            }
            println!();
        }
    }
}
