use std::{self, str};
use std::io::{self, Write};
use std::sync::{mpsc, Arc, Mutex};
use std::collections::BTreeMap;
use std::io::Read;

use seq_io::fastq;
use crossbeam;
use debruijn::dna_string::DnaString;
use debruijn::Kmer;
use debruijn_mapping::pseudoaligner::intersect;
use failure::Error;

use crate::core_genome::CoreGenomePseudoaligner;

fn read_coverage_threshold() -> usize {
    debruijn_mapping::config::READ_COVERAGE_THRESHOLD
}

// TODO: Generalise this to take gzipping fastq, fasta etc, and use a faster
// reader.
pub fn process_reads<K: Kmer + Sync + Send, R: Read + Sync + Send>(
    mut forward_or_single_reader: fastq::Reader<R>,
    mut reverse_reader: Option<fastq::Reader<R>>,
    core_genome_pseudoaligner: &CoreGenomePseudoaligner<K>,
    num_threads: usize,
    // Result returned is equivalence class indices mapped to indexes, and
    // counts for each observed index.
) -> Result<(std::collections::BTreeMap<Vec<u32>, usize>, Vec<usize>, Vec<usize>), Error> {
    info!("Done Reading index");
    info!("Starting Multi-threaded Mapping");

    let (tx, rx) = mpsc::sync_channel(num_threads);

    let mapping_pairs = reverse_reader.is_some();
    let atomic_readers = match reverse_reader {
        Some(ref mut s) => Arc::new(Mutex::new(
                (forward_or_single_reader.records(),
                 Some(s.records()))
        )),
        None => Arc::new(Mutex::new(
                (forward_or_single_reader.records(),
                 None)
        ))
    };

    let mut eq_class_indices: BTreeMap<Vec<u32>, usize> = BTreeMap::new();
    let mut eq_class_coverages: Vec<usize> = vec![];
    let mut eq_class_read_counts: Vec<usize> = vec![];

    info!("Spawning {} threads for Mapping.", num_threads);
    // TODO: This scope section has been modified from an older version
    // rust-pseudoaligner. Changes there should potentially be ported here, but
    // that will have to be done manually.
    crossbeam::scope(|scope| {
        for _ in 0..num_threads {
            let tx = tx.clone();
            let readers = Arc::clone(&atomic_readers);

            scope.spawn(move |_| {
                loop {
                    // TODO: Is is slow here to run a function for each read?
                    // Presumably? Fix by benchmarking inline(always)
                    let wrapped_read_data = if mapping_pairs {
                        // If work is available, do that work.
                        match get_next_record_pair(&readers) {
                            Some(recs) => {
                                match recs {
                                    Ok((fwd_record, rev_record)) => {
                                        map_read_pair(&fwd_record, &rev_record, core_genome_pseudoaligner)
                                    },
                                    Err(err) => panic!("Error {:?} in reading paired fastq", err)
                                }
                            },
                            None => {
                                tx.send(None).expect("Could not send finishing sentinal None while mapping paired");
                                break;
                            }
                        }
                    } else {
                        match get_next_first_record(&readers) {
                            Some(result_record) => {
                                match result_record {
                                    Ok(record) => map_single_reads(&record, core_genome_pseudoaligner),
                                    Err(err) => panic!("Error {:?} in reading fastq", err),
                                }
                            },
                            None => {
                                tx.send(None).expect("Could not send finishing sentinal None while mapping single");
                                break;
                            }
                        }
                    };
                    tx.send(Some(wrapped_read_data)).expect("Could not send data!");
                }
            });
        }

        let mut read_counter: usize = 0;
        let mut mapped_read_counter: usize = 0;
        let mut dead_thread_count = 0;

        for eq_class in rx.iter() {
            match eq_class {
                None => {
                    dead_thread_count += 1;
                    if dead_thread_count == num_threads {
                        drop(tx);
                        // can't continue with a flag check
                        // weird Rusty way !
                        // Consume whatever is remaining
                        // Not worrying about counters; hunch is their
                        // should be less
                        for eq_class in rx.iter() {
                            eq_class.map_or((), |eq_class| debug!("eq_class: {:?}", eq_class));
                        }
                        break;
                    }
                }
                Some(read_data) => {
                    //println!("{:?}, {}", read_data, read_data.0);

                    if read_data.0 {
                        mapped_read_counter += 1;

                        let classes = read_data.2;
                        let coverage = read_data.3;

                        trace!("split cov {}, from {:?}", coverage, classes);
                        let mut classes_sorted = classes.clone();
                        classes_sorted.sort();

                        if eq_class_indices.contains_key(&classes_sorted) {
                            let i = *eq_class_indices.get(&classes_sorted).unwrap();
                            eq_class_coverages[i] += coverage;
                            eq_class_read_counts[i] += 1;
                        } else {
                            let index = eq_class_indices.len();
                            eq_class_indices.insert(classes_sorted, index);
                            eq_class_coverages.push(coverage);
                            eq_class_read_counts.push(1);
                        }
                    }

                    read_counter += 1;
                    if read_counter % 1_000_000 == 0 {
                        let frac_mapped = mapped_read_counter as f32 * 100.0 / read_counter as f32;
                        info!(
                            "Done Mapping {} reads w/ Rate: {}",
                            read_counter, frac_mapped
                        );
                        io::stderr().flush().expect("Could not flush stdout");
                    }
                } // end-Some
            } // end-match
        } // end-for

        info!("Found {} reads mapped out of {}", mapped_read_counter, read_counter);
    }).expect("crossbeam result failure"); //end crossbeam

    debug!("Result: {:?}, {:?}, {:?}", eq_class_indices, eq_class_coverages, eq_class_read_counts);

    info!("Done Mapping Reads");
    Ok((eq_class_indices, eq_class_coverages, eq_class_read_counts))
}

fn get_next_first_record<R: io::Read>(
    reader: &Arc<Mutex<(fastq::RecordsIter<R>, Option<fastq::RecordsIter<R>>)>>,
) -> Option<Result<fastq::OwnedRecord, seq_io::fastq::Error>> {
    let mut lock = reader.lock().unwrap();
    lock.0.next()
}

fn get_next_record_pair<R: io::Read>(
    reader: &Arc<Mutex<(fastq::RecordsIter<R>, Option<fastq::RecordsIter<R>>)>>,
) -> Option<Result<(fastq::OwnedRecord, fastq::OwnedRecord), seq_io::fastq::Error>> {
    let mut lock = reader.lock().unwrap();
    let f = lock.0.next();
    let r = match lock.1 {
        Some(ref mut r_reader) => r_reader.next(),
        None => unreachable!()
    };

    if f.is_some() && r.is_some() {
        return match f.unwrap() {
            Ok(fq) => {
                match r.unwrap() {
                    Ok(rq) => Some(Ok((fq,rq))),
                    Err(e) => Some(Err(e))
                }
            },
            Err(e) => Some(Err(e))
        }
    } else if f.is_some() || r.is_some() {
        panic!("Detected a different number of reads in forward and \
                reverse read files, this shouldn't happen.")
    } else {
        // Natural end to both files.
        return None
    }
}

fn add_coverage(c: &Option<(std::vec::Vec<u32>, usize)>)
                -> usize {
    match c {
        Some((_, coverage)) => *coverage,
        None => 0
    }
}


fn map_read_pair<Q: fastq::Record, K: Kmer+Send+Sync>(
    fwd_record: &Q,
    rev_record: &Q,
    core_genome_pseudoaligner: &CoreGenomePseudoaligner<K>)
    -> (bool, String, Vec<u32>, usize) {

    // Read sequences and do the mapping
    let fwd_seq = str::from_utf8(fwd_record.seq()).unwrap();
    trace!("Mapping forward DNA string: {:?}", fwd_seq);
    let fwd_coverages = core_genome_pseudoaligner.map_read(&DnaString::from_dna_string(fwd_seq));

    let rev_seq = str::from_utf8(rev_record.seq()).unwrap();
    trace!("Mapping reverse DNA string: {:?}", rev_seq);
    let rev_coverages = core_genome_pseudoaligner.map_read(&DnaString::from_dna_string(rev_seq));

    // TODO: Check for mismatching sequence names like BWA does

    // Process TODO: The orientations that are reported for the fwd and rev may
    // actually be random if they are mapped to different nodes or sets of
    // nodes, and it isn't clear how to match them together without traversing
    // the graph. Some options here: (1) Take the highest coverage for the fwd
    // and rev separately, and merge those, or (2) Take the highest coverage out
    // of all 4 possibilities and go with that, or (3) Somehow take into account
    // the eq_classes (and maybe even genomes). For now at least that is all too
    // hard, just going with (1).
    trace!("Found forward read coverage {} and reverse read coverage {}",
        add_coverage(&fwd_coverages), add_coverage(&rev_coverages));

    // TODO: Maybe calculating the read_name isn't needed
    let read_name = fwd_record.id().to_owned().expect("UTF-8 parsing error in read name").to_string();

    let total_coverage = add_coverage(&fwd_coverages) + add_coverage(&rev_coverages);
    return if total_coverage > read_coverage_threshold()*2 {
        (
            true,
            read_name,
            match (fwd_coverages, rev_coverages) {
                (Some((mut f_eq_class, _)), Some((r_eq_class, _))) => {
                    // Intersect operates in place
                    intersect(&mut f_eq_class, &r_eq_class);
                    f_eq_class
                },
                (Some((f_eq_class, _)), None) => f_eq_class,
                (None, Some((r_eq_class, _))) => r_eq_class,
                (None, None) => unreachable!()
            },
            total_coverage,
        )
    } else {
        // Don't waste time calculating unused info here
        (false, read_name, vec![], 0)
    };
}




fn map_single_reads<Q: fastq::Record, K: Kmer+Send+Sync>(
    record: &Q,
    core_genome_pseudoaligner: &CoreGenomePseudoaligner<K>)
    -> (bool, String, Vec<u32>, usize)  {

    let read_name = record.id().to_owned().expect("UTF-8 parsing error in read name").to_string();

    let seq = str::from_utf8(record.seq()).unwrap();
    match core_genome_pseudoaligner.map_read(&DnaString::from_dna_string(seq)) {
        Some((eq_class, coverage)) => (true, read_name, eq_class, coverage),
        None => (false, read_name, Vec::new(), 0)
    }
}
