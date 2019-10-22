// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

use std::{self, cmp::Ordering, str};
use std::collections::HashMap;
use std::io::{self, Write};
use std::sync::{mpsc, Arc, Mutex};
use std::collections::BTreeMap;
use std::io::Read;

use seq_io::fastq;
use boomphf::hashmap::NoKeyBoomHashMap;
use crossbeam;
use debruijn::dna_string::DnaString;
use debruijn::filter::EqClassIdType;
use debruijn::graph::DebruijnGraph;
use debruijn::{Dir, Kmer, Mer, Vmer};
use failure::Error;

use pseudoaligner::config::{READ_COVERAGE_THRESHOLD, LEFT_EXTEND_FRACTION};

#[derive(Serialize, Deserialize, Debug)]
pub struct Pseudoaligner<K: Kmer> {
    pub dbg: DebruijnGraph<K, EqClassIdType>,
    pub eq_classes: Vec<Vec<u32>>,
    pub dbg_index: NoKeyBoomHashMap<K, (u32, u32)>,
    pub tx_names: Vec<String>,
    pub tx_gene_mapping: HashMap<String, String>,
}

pub trait PseudoalignmentReadMapper {
    /// Pseudo-align `read_seq` to determine its the equivalence class.
    fn map_read(&self, read_seq: &DnaString) -> Option<(Vec<u32>, usize)>;
}

impl<K: Kmer + Sync + Send> Pseudoaligner<K> {
    pub fn new(
        dbg: DebruijnGraph<K, EqClassIdType>,
        eq_classes: Vec<Vec<u32>>,
        dbg_index: NoKeyBoomHashMap<K, (u32, u32)>,
        tx_names: Vec<String>,
        tx_gene_mapping: HashMap<String, String>
    ) -> Pseudoaligner<K> {
        Pseudoaligner {dbg, eq_classes, dbg_index, tx_names, tx_gene_mapping}
    }
}

impl<K: Kmer + Send + Sync> PseudoalignmentReadMapper for Pseudoaligner<K> {
    fn map_read(&self, read_seq: &DnaString) -> Option<(Vec<u32>, usize)> {
        let read_length = read_seq.len();
        let mut read_coverage: usize = 0;
        let mut colors: Vec<u32> = Vec::new();
        let left_extend_threshold = (LEFT_EXTEND_FRACTION * read_length as f64) as usize;

        let mut kmer_pos: usize = 0;
        let kmer_length = K::k();
        let last_kmer_pos = read_length - kmer_length;

        // Scan the read for the first kmer that exists in the reference
        let find_kmer_match = |kmer_pos: &mut usize| -> Option<(usize, usize)> {
            while *kmer_pos <= last_kmer_pos {
                debug!("Determining kmer at position {}", *kmer_pos);
                let read_kmer = read_seq.get_kmer(*kmer_pos);

                match self.dbg_index.get(&read_kmer) {
                    None => (),
                    Some((nid, offset)) => {
                        let node = self.dbg.get_node(*nid as usize);
                        debug!("kmer hit to node {:?}", node);
                        // Check that this is a real hit and the kmer is
                        // actually in the MPHF.
                        let ref_seq_slice = node.sequence();
                        let ref_kmer: K = ref_seq_slice.get_kmer(*offset as usize);

                        if read_kmer == ref_kmer {
                            return Some((*nid as usize, *offset as usize));
                        }
                    }
                };
                *kmer_pos += 1;
            }

            None
        };

        // extract the first exact matching position of read
        let (mut node_id, mut kmer_offset) =
            // get the first match through mphf
            match find_kmer_match(&mut kmer_pos) {
                None => (None, None),
                Some((nid, offset)) => (Some(nid), Some(offset))
            };

        // check if we can extend back if there were SNP in every kmer query
        if kmer_pos >= left_extend_threshold && node_id.is_some() {
            let mut last_pos = kmer_pos - 1;
            let mut prev_node_id = node_id.unwrap();
            let mut prev_kmer_offset = if kmer_offset.unwrap() > 0 { kmer_offset.unwrap() - 1 } else { 0 };

            loop {
                let node = self.dbg.get_node(prev_node_id);
                //println!("{:?}, {:?}, {:?}, {:?}, {:?}",
                //         node, node.sequence(),
                //         &eq_classes[ *node.data() as usize],
                //         prev_kmer_offset, last_pos);

                // length of remaining read before kmer match
                let skipped_read = last_pos + 1;

                // length of the skipped node sequence before kmer match
                let skipped_ref = prev_kmer_offset + 1;

                // find maximum extention possbile before fork or eof read
                let max_matchable_pos = std::cmp::min(skipped_read, skipped_ref);

                let ref_seq_slice = node.sequence();
                let mut premature_break = false;
                let mut matched_bases = 0;
                let mut seen_snp = 0;
                for idx in 0..max_matchable_pos {
                    let ref_pos = prev_kmer_offset - idx;
                    let read_offset = last_pos - idx;

                    // compare base by base
                    if ref_seq_slice.get(ref_pos) != read_seq.get(read_offset) {
                        if seen_snp > 3 {
                            premature_break = true;
                            break;
                        }

                        // Allowing 2-SNP
                        seen_snp += 1;
                    }

                    matched_bases += 1;
                    read_coverage += 1;
                }

                //break the loop if end of read reached or a premature mismatch
                if last_pos + 1 - matched_bases == 0 || premature_break {
                    break;
                }

                // adjust last position
                last_pos -= matched_bases;

                // If reached here then a fork is found in the reference.
                let exts = node.exts();
                let next_base = read_seq.get(last_pos);
                if exts.has_ext(Dir::Left, next_base) {
                    // found a left extention.
                    let index = exts
                        .get(Dir::Left)
                        .iter()
                        .position(|&x| x == next_base)
                        .unwrap();

                    let edge = node.l_edges()[index];

                    //update the previous node's id
                    prev_node_id = edge.0;
                    let prev_node = self.dbg.get_node(prev_node_id);
                    prev_kmer_offset = prev_node.sequence().len() - kmer_length;

                    // extract colors
                    let color = prev_node.data();
                    colors.push(*color);
                } else {
                    break;
                }
            } // end-loop
        } //end-if

        // forward search
        if kmer_pos <= last_kmer_pos {
            loop {
                let node = self.dbg.get_node(node_id.unwrap());
                //println!("{:?}, {:?}, {:?}, {:?}",
                //         node, node.sequence(),
                //         &eq_classes[ *node.data() as usize],
                //         kmer_offset);
                kmer_pos += kmer_length;
                read_coverage += kmer_length;

                // extract colors
                let color = node.data();
                colors.push(*color);

                // length of remaining read after kmer match
                let remaining_read = read_length - kmer_pos;

                // length of the remaining node sequence after kmer match
                let ref_seq_slice = node.sequence();
                let ref_length = ref_seq_slice.len();
                let ref_offset = kmer_offset.unwrap() + kmer_length;
                let informative_ref = ref_length - ref_offset;

                // find maximum extention possbile before fork or eof read
                let max_matchable_pos = std::cmp::min(remaining_read, informative_ref);

                let mut premature_break = false;
                let mut matched_bases = 0;
                let mut seen_snp = 0;
                for idx in 0..max_matchable_pos {
                    let ref_pos = ref_offset + idx;
                    let read_offset = kmer_pos + idx;

                    // compare base by base
                    if ref_seq_slice.get(ref_pos) != read_seq.get(read_offset) {
                        if seen_snp > 3 {
                            premature_break = true;
                            break;
                        }

                        // Allowing 2-SNP
                        seen_snp += 1;
                    }

                    matched_bases += 1;
                    read_coverage += 1;
                }

                kmer_pos += matched_bases;
                //break the loop if end of read reached or a premature mismatch
                if kmer_pos >= read_length {
                    break;
                }

                // If reached here then a fork is found in the reference.
                let exts = node.exts();
                let next_base = read_seq.get(kmer_pos);

                if !premature_break && exts.has_ext(Dir::Right, next_base) {
                    // found a right extention.
                    let index = exts
                        .get(Dir::Right)
                        .iter()
                        .position(|&x| x == next_base)
                        .unwrap();

                    let edge = node.r_edges()[index];

                    //update the next node's id
                    node_id = Some(edge.0);
                    kmer_offset = Some(0);

                    //adjust for kmer_position
                    kmer_pos -= kmer_length - 1;
                    read_coverage -= kmer_length - 1;
                } else {
                    // can't extend node in dbg extract read using mphf
                    // TODO: might have to check some cases
                    if kmer_pos > last_kmer_pos {
                        // can't search in mphf if no full kmer can be made
                        break;
                    }

                    // get the match through mphf
                    match find_kmer_match(&mut kmer_pos) {
                        None => break,
                        Some((nid, offset)) => {
                            node_id = Some(nid);
                            kmer_offset = Some(offset);
                        }
                    };
                }
            } // end-loop
        } //end-if

        // Take the intersection of the sets
        let colors_len = colors.len();
        if colors_len == 0 {
            if read_coverage != 0 {
                panic!(
                    "Different read coverage {:?} than num of eqclasses {:?}",
                    colors_len, read_coverage
                );
            }

            None
        } else {
            // Intersect the equivalence classes
            let first_color = colors.pop().unwrap();
            let mut eq_class = self.eq_classes[first_color as usize].clone();

            for color in colors {
                intersect(&mut eq_class, &self.eq_classes[color as usize]);
            }

            Some((eq_class, read_coverage))
        }
    }
}

/// Compute the intersection of v1 and v2 inplace on top of v1
/// v1 and v2 must be sorted
pub fn intersect<T: Eq + Ord>(v1: &mut Vec<T>, v2: &[T]) {
    if v1.is_empty() {
        return;
    }

    if v2.is_empty() {
        v1.clear();
    }

    let mut fill_idx1 = 0;
    let mut idx1 = 0;
    let mut idx2 = 0;

    while idx1 < v1.len() && idx2 < v2.len() {
        match v1[idx1].cmp(&v2[idx2]) {
            Ordering::Less => idx1 += 1,
            Ordering::Greater => idx2 += 1,
            Ordering::Equal => {
                v1.swap(fill_idx1, idx1);
                idx1 += 1;
                idx2 += 1;
                fill_idx1 += 1;
            }
        }
    }
    // Can use unreachable!() here because fill_idx1 <= v1.len()
    v1.resize_with(fill_idx1, || { unreachable!() });
}

// TODO: Generalise this to take gzipping fastq, fasta etc, and use a faster
// reader.
pub fn process_reads<K: Kmer + Sync + Send, T: PseudoalignmentReadMapper + Sync + Send, R: Read + Sync + Send>(
    mut forward_or_single_reader: fastq::Reader<R>,
    mut reverse_reader: Option<fastq::Reader<R>>,
    index: &T,
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
                                        map_read_pair(&fwd_record, &rev_record, index)
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
                                    Ok(record) => map_single_reads(&record, index),
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

                        debug!("split cov {}, from {:?}", coverage, classes);
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

fn calculate_fwd_rev<T: PseudoalignmentReadMapper>(
    seq: &str,
    index: &T)
    -> (Option<(Vec<u32>, usize)>, Option<(Vec<u32>, usize)>) {

    let dna_string = DnaString::from_dna_string(seq);
    debug!("Mapping forward DNA string: {:?}", dna_string);
    let fwd_classes = index.map_read(&dna_string);

    debug!("Mapping reverse DNA string: {:?}", dna_string.rc());
    let rev_classes = index.map_read(&dna_string.rc());
    debug!("Found fwd eq_classes {:?} and reverse {:?}", fwd_classes, rev_classes);

    (fwd_classes, rev_classes)
}

fn add_coverage(c: &Option<(std::vec::Vec<u32>, usize)>)
                -> usize {
    match c {
        Some((_, coverage)) => *coverage,
        None => 0
    }
}


fn map_read_pair<T: PseudoalignmentReadMapper, Q: fastq::Record>(
    fwd_record: &Q,
    rev_record: &Q,
    index: &T)
    -> (bool, String, Vec<u32>, usize) {

    // Read sequences and do the mapping
    let fwd_seq = str::from_utf8(fwd_record.seq()).unwrap();
    let (fwd_fwd_classes, fwd_rev_classes) = calculate_fwd_rev(fwd_seq, index);

    let rev_seq = str::from_utf8(rev_record.seq()).unwrap();
    let (rev_fwd_classes, rev_rev_classes) = calculate_fwd_rev(rev_seq, index);

    // TODO: Check for mismatching sequence names like BWA does

    // Process TODO: The orientations that are reported for the fwd and rev may
    // actually be random if they are mapped to different nodes or sets of
    // nodes, and it isn't clear how to match them together without traversing
    // the graph. Some options here: (1) Take the highest coverage for the fwd
    // and rev separately, and merge those, or (2) Take the highest coverage out
    // of all 4 possibilities and go with that, or (3) Somehow take into account
    // the eq_classes (and maybe even genomes). For now at least that is all too
    // hard, just going with (1).
    debug!("Found forward coverage components {} {} and reverse components {} {}",
        add_coverage(&fwd_fwd_classes), add_coverage(&fwd_rev_classes),
        add_coverage(&rev_fwd_classes), add_coverage(&rev_rev_classes));
    let fwd_coverages = if add_coverage(&fwd_fwd_classes) > add_coverage(&fwd_rev_classes) {
        fwd_fwd_classes
    } else {
        fwd_rev_classes
    };
    let rev_coverages = if add_coverage(&rev_fwd_classes) > add_coverage(&rev_rev_classes) {
        rev_fwd_classes
    } else {
        rev_rev_classes
    };

    // TODO: Maybe calculating the read_name isn't needed
    let read_name = fwd_record.id().to_owned().expect("UTF-8 parsing error in read name").to_string();

    let total_coverage = add_coverage(&fwd_coverages) + add_coverage(&rev_coverages);
    return if total_coverage > READ_COVERAGE_THRESHOLD*2 {
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




fn map_single_reads<T: PseudoalignmentReadMapper, Q: fastq::Record>(
    record: &Q,
    index: &T)
    -> (bool, String, Vec<u32>, usize)  {

    let seq = str::from_utf8(record.seq()).unwrap();
    let (fwd_classes, rev_classes) = calculate_fwd_rev(seq, index);
    let read_name = record.id().to_owned().expect("UTF-8 parsing error in read name").to_string();

    return match (fwd_classes, rev_classes) {
        (None, Some((eq_class, coverage))) | (Some((eq_class, coverage)), None) => {
            if coverage >= READ_COVERAGE_THRESHOLD && !eq_class.is_empty() {
                (true, read_name, eq_class, coverage)
            } else {
                (false, read_name, eq_class, coverage)
            }
        },
        (Some((eq_class_fwd, coverage_fwd)),
         Some((eq_class_rev, coverage_rev))) => {
            // Both match, take the read with the highest coverage
            match (eq_class_fwd.is_empty(), eq_class_rev.is_empty()) {
                (false, false) => {
                    if coverage_fwd > coverage_rev {
                        if coverage_fwd > READ_COVERAGE_THRESHOLD {
                            (true, read_name, eq_class_fwd, coverage_fwd)
                        } else {
                            (false, read_name, eq_class_fwd, coverage_fwd)
                        }
                    } else {
                        if coverage_rev > READ_COVERAGE_THRESHOLD {
                            (true, read_name, eq_class_rev, coverage_rev)
                        } else {
                            (false, read_name, eq_class_rev, coverage_rev)
                        }
                    }
                },
                (true, false) => {
                    if coverage_fwd > READ_COVERAGE_THRESHOLD {
                        (true, read_name, eq_class_fwd, coverage_fwd)
                    } else {
                        (false, read_name, eq_class_fwd, coverage_fwd)
                    }
                },
                (false, true) => {
                    if coverage_rev > READ_COVERAGE_THRESHOLD {
                        (true, read_name, eq_class_rev, coverage_rev)
                    } else {
                        (false, read_name, eq_class_rev, coverage_rev)
                    }
                },
                (true, true) => (false, read_name, Vec::new(), 0)
            }
        },
        (None, None) => (false, read_name, Vec::new(), 0)
    }
}
