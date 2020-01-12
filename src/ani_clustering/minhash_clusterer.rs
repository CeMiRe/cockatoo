
use std::collections::{BTreeSet,BTreeMap};
use std::io::{BufReader,Write};
use std::sync::Mutex;
use finch::distance::distance;
use finch::serialization::Sketch;
use rayon::prelude::*;
use finish_command_safely;
use partitions::partition_vec::PartitionVec;

/// Given a list of genomes, return them clustered. If the fastani_threshold is
/// set, use minhash for first pass analysis, then fastani as the actual threshold.
// TODO: Test whether this is a good enough procedure or if there are nasties
// e.g. failure to cluster bad quality genomes.
pub fn minhash_clusters(
    genomes: &[&str],
    minhash_ani: f32,
    n_hashes: usize,
    kmer_length: u8,
    fastani_threshold: Option<f32>,
) -> Vec<Vec<usize>> {

    // Generate sketches for all input files
    let mut filter = finch::filtering::FilterParams { // dummy, no filtering is applied.
        filter_on: None,
        abun_filter: (None, None),
        err_filter: 0f32,
        strand_filter: 0f32,
    };
    info!("Sketching genomes for clustering ..");
    let sketches = finch::mash_files(
        genomes,
        n_hashes,
        n_hashes,
        kmer_length,
        &mut filter,
        true,
        0)
        .expect("Failed to create finch sketches for input genomes");
    info!("Finished sketching genomes for clustering.");

    let distance_threshold: f64 = (100.0 - minhash_ani as f64)/100.0;
    assert!(distance_threshold >= 0.0);
    assert!(distance_threshold < 1.0);

    match fastani_threshold {
        None => {
            // Straight up minhash clustering.

            // Greedily find reps
            let clusters = find_minhash_representatives(&sketches.sketches, distance_threshold);

            // Reassign non-reps based so they are assigned to the nearest
            // representative.
            return find_minhash_memberships(&clusters, &sketches.sketches);
        },
        Some(fastani_threshold) => {
            let minhash_preclusters = partition_sketches(&sketches.sketches, distance_threshold);
            trace!("Found minhash_preclusters: {:?}", minhash_preclusters);

            let all_clusters: Mutex<Vec<Vec<usize>>> = Mutex::new(vec![]);

            for genome_precluster_indices in minhash_preclusters.all_sets() {
                let mut precluster_sketches = vec![];
                let mut precluster_genomes = vec![];
                let mut original_genome_indices: Vec<usize> = vec![];
                for cluster_genome in genome_precluster_indices {
                    precluster_sketches.push(&sketches.sketches[*cluster_genome.1]);
                    precluster_genomes.push(genomes[*cluster_genome.1]);
                    original_genome_indices.push(*cluster_genome.1);
                }
                debug!("Clustering pre-cluster {:?}", original_genome_indices);

                let (clusters, calculated_fastanis) = find_minhash_fastani_representatives(
                    precluster_sketches.as_slice(), precluster_genomes.as_slice(), distance_threshold, fastani_threshold
                );

                let clusters = find_minhash_fastani_memberships(
                    &clusters, precluster_sketches.as_slice(), precluster_genomes.as_slice(), calculated_fastanis, distance_threshold
                );
                // Indices here are within this single linkage cluster only, so
                // here we map them back to their original indices, and then
                // insert back into the original array.
                for cluster in clusters {
                    all_clusters.lock().unwrap().push(
                        cluster.iter().map(|within_index| original_genome_indices[*within_index]).collect()
                    );
                }
            }
            return all_clusters.into_inner().unwrap();
        }
    }

}

/// Create sub-sets by single linkage clustering
fn partition_sketches(
    sketches: &[Sketch],
    distance_threshold: f64) 
    -> PartitionVec<usize> {

    let mut to_return: PartitionVec<usize> = PartitionVec::with_capacity(sketches.len());
    for (i, _) in sketches.iter().enumerate() {
        to_return.push(i);
    }

    for (i, sketch1) in sketches.iter().enumerate() {
        for (j, sketch2) in sketches.iter().enumerate() {
            trace!("Testing precluster between {} and {}", i, j);
            if j >= i {
                continue;
            }
            if distance(&sketch1.hashes, &sketch2.hashes, "", "", true)
                .expect("Failed to calculate distance by sketch comparison")
                .mashDistance
                <= distance_threshold {

                debug!("During preclustering, found a match between genomes {} and {}", i, j);
                to_return.union(i,j)
            }
        }
    }
    return to_return;
}

/// Choose representatives, greedily assigning based on the min_ani threshold.
fn find_minhash_representatives(
    sketches: &[Sketch],
    ani_threshold: f64)
    -> BTreeSet<usize> {

    let mut to_return: BTreeSet<usize> = BTreeSet::new();

    for (i, sketch1) in sketches.iter().enumerate() {
        let mut is_rep = true;
        for j in &to_return {
            let sketch2: &Sketch = &sketches[*j];
            if distance(&sketch1.hashes, &sketch2.hashes, "", "", true)
                .expect("Failed to calculate distance by sketch comparison")
                .mashDistance
                <= ani_threshold {

                is_rep = false;
                break;
            }
        }
        if is_rep {
            to_return.insert(i);
        }
    }
    return to_return;
}

fn find_minhash_fastani_representatives(
    sketches: &[&Sketch],
    genomes: &[&str],
    minhash_ani_threshold: f64,
    fastani_threshold: f32)
    -> (BTreeSet<usize>, BTreeMap<(usize, usize), Option<f32>>) {

    let mut clusters_to_return: BTreeSet<usize> = BTreeSet::new();
    let mut fastani_cache: BTreeMap<(usize, usize), Option<f32>> = BTreeMap::new();

    for (i, sketch1) in sketches.iter().enumerate() {
        // Gather a list of all genomes which pass the minhash threshold.
        let potential_refs: Vec<_> = clusters_to_return.iter().filter( |j| {
            let sketch2: &Sketch = &sketches[**j];
            distance(&sketch1.hashes, &sketch2.hashes, "", "", true)
                .expect("Failed to calculate distance by sketch comparison")
                .mashDistance
                <= minhash_ani_threshold
        }).collect();

        // FastANI all potential reps against the current genome
        let fastanis = calculate_fastani_many_to_one(
            potential_refs.iter().map(|ref_index| genomes[**ref_index]).collect::<Vec<_>>().as_slice(),
            genomes[i]);
        let mut is_rep = true;
        for (potential_ref, fastani) in potential_refs.into_iter().zip(fastanis.iter()) {
            // The reps always have a lower index that the genome under
            // consideration, and the cache key is in ascending order.
            fastani_cache.insert((*potential_ref,i), *fastani);
            match fastani {
                Some(ani) => if *ani >= fastani_threshold {
                    is_rep = false
                },
                None => {}
            }
        }
        if is_rep {
            clusters_to_return.insert(i);
        }
    }
    return (clusters_to_return, fastani_cache);
}

fn calculate_fastani_many_to_one(
    query_genome_paths: &[&str],
    ref_genome_path: &str)
-> Vec<Option<f32>> {
    if query_genome_paths.len() == 0 {
        return vec![];
    }

    let mut query_to_index: BTreeMap<String,usize> = BTreeMap::new();
    for (i, path) in query_genome_paths.iter().enumerate() {
        query_to_index.insert(path.to_string(), i);
    }
    
    let mut cmd = std::process::Command::new("fastANI");
    cmd
        .arg("-o")
        .arg("/dev/stdout")
        .arg("--queryList")
        .arg("/dev/stdin")
        .arg("--ref")
        .arg(&ref_genome_path)
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());
    debug!("Running fastANI command: {:?}", &cmd);
    let mut process = cmd.spawn().expect(&format!("Failed to spawn {}", "fastANI"));
    let stdout = process.stdout.as_mut().unwrap();
    let stdout_reader = BufReader::new(stdout);

    {
        let mut stdin = process.stdin.take().unwrap();
        for query in query_genome_paths {
            writeln!(stdin, "{}", query).expect("Failed to write to FastANI stdin");
        }
        stdin.flush().expect("Failed to flush FastANI stdin");
    }

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(stdout_reader);

    let mut to_return = vec![None; query_genome_paths.len()];

    for record_res in rdr.records() {
        match record_res {
            Ok(record) => {
                assert!(record.len() == 5);
                let ani: f32 = record[2].parse().expect("Failed to convert fastani ANI to float value");
                let query_name = &record[0];
                let query_index: usize = *query_to_index.get(query_name)
                    .expect(&format!("FastANI parse error: Failed to know what query genome {} is", query_name));
                to_return[query_index] = Some(ani);
                debug!("FastANI of {} against {} was {}", ref_genome_path, query_name, ani);
            },
            Err(e) => {
                error!("Error parsing fastani output: {}", e);
                std::process::exit(1);
            }
        }
    }
    finish_command_safely(process, "FastANI");
    return to_return;
}

/// For each genome (sketch) assign it to the closest representative genome:
fn find_minhash_memberships(
    representatives: &BTreeSet<usize>,
    sketches: &[Sketch],
) -> Vec<Vec<usize>> {

    let mut rep_to_index = BTreeMap::new();
    for (i, rep) in representatives.iter().enumerate() {
        rep_to_index.insert(rep, i);
    }

    let mut to_return: Vec<Vec<usize>> = vec![vec![]; representatives.len()];
    for (i, sketch1) in sketches.iter().enumerate() {
        if representatives.contains(&i) {
            to_return[rep_to_index[&i]].push(i);
        } else {
            let mut best_rep_min_ani = None;
            let mut best_rep = None;
            for rep in representatives.iter() {
                let dist = distance(&sketch1.hashes, &sketches[*rep].hashes, "", "", true)
                    .expect("Failed to calculate distance by sketch comparison")
                    .mashDistance;
                if best_rep_min_ani.is_none() || dist < best_rep_min_ani.unwrap() {
                    best_rep = Some(rep);
                    best_rep_min_ani = Some(dist);
                }
            }
            to_return[rep_to_index[best_rep.unwrap()]].push(i);
        }
    }
    return to_return;
}

/// Given a list of representative genomes and a list of genomes, assign each
/// genome to the closest representative.
fn find_minhash_fastani_memberships(
    representatives: &BTreeSet<usize>,
    sketches: &[&Sketch],
    genomes: &[&str],
    calculated_fastanis: BTreeMap<(usize, usize), Option<f32>>,
    minhash_threshold: f64
) -> Vec<Vec<usize>> {
    
    let mut rep_to_index = BTreeMap::new();
    for (i, rep) in representatives.iter().enumerate() {
        rep_to_index.insert(rep, i);
    }

    let to_return: Mutex<Vec<Vec<usize>>> = Mutex::new(vec![vec![]; representatives.len()]);
    let calculated_fastanis_lock = Mutex::new(calculated_fastanis);

    sketches.par_iter().enumerate().for_each( |(i, sketch1)| {
        if representatives.contains(&i) {
            to_return.lock().unwrap()[rep_to_index[&i]].push(i);
        } else {
            let potential_refs: Vec<&usize> = representatives.into_iter().filter(|rep| {
                let (min_i, max_i) = if i < **rep {
                    (i, **rep)
                } else {
                    (**rep, i)
                };
                if calculated_fastanis_lock.lock().unwrap().contains_key(&(min_i,max_i)) {
                    true // FastANI not needed since already cached
                } else {
                    let minhash_dist = distance(&sketch1.hashes, &sketches[**rep].hashes, "", "", true)
                    .expect("Failed to calculate distance by sketch comparison")
                    .mashDistance;
                    minhash_dist < minhash_threshold
                }
            }).collect();

            let new_fastanis = calculate_fastani_many_to_one(
                &potential_refs.iter().map(|ref_i| genomes[**ref_i]).collect::<Vec<_>>(),
                genomes[i],
            );
            for (ref_i, ani_opt) in new_fastanis.iter().enumerate() {
                let (min_i, max_i) = if i < ref_i {
                    (i, ref_i)
                } else {
                    (ref_i, i)
                };
                calculated_fastanis_lock.lock().unwrap().insert((min_i,max_i), *ani_opt);
            };

            // TODO: Is there FastANI bugs here? Donovan/Pierre seem to think so
            let mut best_rep_min_ani = None;
            let mut best_rep = None;
            for rep in representatives.iter() {
                let fastani: Option<f32> = match i < *rep {
                    true => match calculated_fastanis_lock.lock().unwrap().get(&(i, *rep)) {
                        // ani could be known as f32, None for calculated
                        // but below threshold, of not calculated. If not calculated, it isn't nearby
                        Some(ani_opt) => match ani_opt {
                            Some(ani) => Some(*ani),
                            None => None,
                        },
                        None => None
                    },
                    false => match calculated_fastanis_lock.lock().unwrap().get(&(*rep, i)) {
                        Some(ani_opt) => match ani_opt {
                            Some(ani) => Some(*ani),
                            None => None,
                        },
                        None => None
                    }
                };

                if fastani.is_some() && (best_rep_min_ani.is_none() || fastani.unwrap() > best_rep_min_ani.unwrap()) {
                    best_rep = Some(rep);
                    best_rep_min_ani = fastani;
                }
            }
            to_return.lock().unwrap()[rep_to_index[best_rep.unwrap()]].push(i);
        }
    });

    return to_return.into_inner().unwrap();
}

#[cfg(test)]
mod tests {
    use super::*;

    fn init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn test_minhash_hello_world() {
        init();
        let clusters = minhash_clusters(
            &["tests/data/parsnp/1_first_group/73.20120800_S1X.13.fna",
              "tests/data/parsnp/1_first_group/73.20120600_S2D.19.fna",
              "tests/data/parsnp/1_first_group/73.20120700_S3X.12.fna",
              "tests/data/parsnp/1_first_group/73.20110800_S2D.13.fna",
            ],
            95.0,
            1000,
            21,
            None,
        );
        assert_eq!(
            vec![vec![0,1,2,3]],
            clusters
        )
    }

    #[test]
    fn test_minhash_two_clusters() {
        init();
        let clusters = minhash_clusters(
            &["tests/data/parsnp/1_first_group/73.20120800_S1X.13.fna",
              "tests/data/parsnp/1_first_group/73.20120600_S2D.19.fna",
              "tests/data/parsnp/1_first_group/73.20120700_S3X.12.fna",
              "tests/data/parsnp/1_first_group/73.20110800_S2D.13.fna",
            ],
            98.0,
            1000,
            21,
            None,
        );
        assert_eq!(
            vec![vec![0,1,3],vec![2]],
            clusters
        )
    }

    #[test]
    fn test_minhash_fastani_hello_world() {
        init();
        let mut clusters = minhash_clusters(
            &["tests/data/parsnp/1_first_group/73.20120800_S1X.13.fna",
              "tests/data/parsnp/1_first_group/73.20120600_S2D.19.fna",
              "tests/data/parsnp/1_first_group/73.20120700_S3X.12.fna",
              "tests/data/parsnp/1_first_group/73.20110800_S2D.13.fna",
            ],
            95.0,
            1000,
            21,
            Some(95.0),
        );
        for cluster in clusters.iter_mut() { cluster.sort_unstable(); }
        assert_eq!(
            vec![vec![0,1,2,3]],
            clusters
        )
    }

    #[test]
    fn test_minhash_fastani_two_clusters_same_ani() {
        init();
        let mut clusters = minhash_clusters(
            &["tests/data/parsnp/1_first_group/73.20120800_S1X.13.fna",
              "tests/data/parsnp/1_first_group/73.20120600_S2D.19.fna",
              "tests/data/parsnp/1_first_group/73.20120700_S3X.12.fna",
              "tests/data/parsnp/1_first_group/73.20110800_S2D.13.fna",
            ],
            98.0,
            1000,
            21,
            Some(98.0),
        );
        for cluster in clusters.iter_mut() { cluster.sort_unstable(); }
        assert_eq!(
            vec![vec![0,1,3],vec![2]],
            clusters
        )
    }

    #[test]
    fn test_minhash_fastani_two_clusters_low_minhash_ani() {
        init();
        let mut clusters = minhash_clusters(
            &["tests/data/parsnp/1_first_group/73.20120800_S1X.13.fna",
              "tests/data/parsnp/1_first_group/73.20120600_S2D.19.fna",
              "tests/data/parsnp/1_first_group/73.20120700_S3X.12.fna",
              "tests/data/parsnp/1_first_group/73.20110800_S2D.13.fna",
            ],
            90.0,
            1000,
            21,
            Some(98.0),
        );
        for cluster in clusters.iter_mut() { cluster.sort_unstable(); }
        assert_eq!(
            vec![vec![0,1,3],vec![2]],
            clusters
        )
    }
}
