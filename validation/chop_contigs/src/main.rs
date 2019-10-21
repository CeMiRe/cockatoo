extern crate seq_io;
use seq_io::fasta::*;
use std::io;
use std::io::prelude::*;

#[macro_use] extern crate log;
extern crate env_logger;

fn main() {
    env_logger::init();
    error!("error-");
    info!("info-");
    debug!("debug-");
    panic!();

    let mut reader = Reader::new(io::stdin());

    while let Some(record) = reader.next() {
        let record = record.expect("Error reading record");
        println!(">{}",record.id().unwrap());
        let seq = record.seq();
        let bound = (seq.len() as f32 *0.8) as usize;
        
        println!("{}", std::str::from_utf8(&seq[0..bound]).unwrap())
    }
}
