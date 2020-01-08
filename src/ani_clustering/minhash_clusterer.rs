
use std::collections::{BTreeSet,BTreeMap};
use std::io::{BufReader};
use finch::distance::distance;
use finch::serialization::Sketch;
use rayon::prelude::*;

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
            let (clusters, calculated_fastanis) = find_minhash_fastani_representatives(
                &sketches.sketches, genomes, distance_threshold, fastani_threshold
            );

            return find_minhash_fastani_memberships(
                &clusters, &sketches.sketches, genomes, &calculated_fastanis, distance_threshold
            );
        }
    }

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
    sketches: &[Sketch],
    genomes: &[&str],
    minhash_ani_threshold: f64,
    fastani_threshold: f32)
    -> (BTreeSet<usize>, BTreeMap<(usize, usize), f32>) {

    let mut clusters_to_return: BTreeSet<usize> = BTreeSet::new();
    let mut fastani_cache: BTreeMap<(usize, usize), f32> = BTreeMap::new();

    for (i, sketch1) in sketches.iter().enumerate() {
        let mut is_rep = true;
        for j in &clusters_to_return {
            let sketch2: &Sketch = &sketches[*j];
            if distance(&sketch1.hashes, &sketch2.hashes, "", "", true)
                .expect("Failed to calculate distance by sketch comparison")
                .mashDistance
                <= minhash_ani_threshold {

                // Possibly not a represenative. Calculate (and cache) FastANI to check.
                let fastani_result = calculate_fastani(genomes[i], genomes[*j]);

                match fastani_result {
                    None => {
                        warn!(
                            "Potentially unexpectedly genome comparison of \
                            {} and {} showed a minhash hit, but no FastANI hit.", 
                            genomes[i], genomes[*j]);
                        // Insert something so it isn't calculated again later
                        fastani_cache.insert((i,*j), 0.0);
                    },
                    Some(fastani) => {
                        if fastani >= fastani_threshold {
                            fastani_cache.insert((i,*j), fastani);
                            is_rep = false;
                            break;
                        }
                    }
                }
            }
        }
        if is_rep {
            clusters_to_return.insert(i);
        }
    }
    return (clusters_to_return, fastani_cache);
}

fn calculate_fastani(fasta1: &str, fasta2: &str) -> Option<f32> {
    let mut cmd = std::process::Command::new("fastANI");
    cmd
        .arg("-o")
        .arg("/dev/stdout")
        .arg("--query")
        .arg(&fasta1)
        .arg("--ref")
        .arg(&fasta2)
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());
    debug!("Running fastANI command: {:?}", &cmd);
    let mut process = cmd.spawn().expect(&format!("Failed to spawn {}", "fastANI"));
    let stdout = process.stdout.as_mut().unwrap();
    let stdout_reader = BufReader::new(stdout);

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(stdout_reader);

    let mut to_return = None;

    for record_res in rdr.records() {
        match record_res {
            Ok(record) => {
                assert!(record.len() == 5);
                let ani: f32 = record[2].parse().expect("Failed to convert fastani ANI to float value");
                if to_return.is_some() {
                    error!("Unexpectedly found >1 result from fastANI");
                    std::process::exit(1);
                    
                }
                to_return = Some(ani);
            },
            Err(e) => {
                error!("Error parsing fastani output: {}", e);
                std::process::exit(1);
            }
        }
    }
    debug!("FastANI of {} against {} was {:?}", fasta1, fasta2, to_return);
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

fn find_minhash_fastani_memberships(
    representatives: &BTreeSet<usize>,
    sketches: &[Sketch],
    genomes: &[&str],
    calculated_fastanis: &BTreeMap<(usize, usize), f32>,
    minhash_threshold: f64
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
                let minhash_dist = distance(&sketch1.hashes, &sketches[*rep].hashes, "", "", true)
                    .expect("Failed to calculate distance by sketch comparison")
                    .mashDistance;
                if minhash_dist < minhash_threshold {
                    let fastani = match i < *rep {
                        true => match calculated_fastanis.get(&(i, *rep)) {
                            Some(ani) => Some(*ani),
                            None => calculate_fastani(genomes[i], genomes[*rep])
                        },
                        false => match calculated_fastanis.get(&(*rep, i)) {
                            Some(ani) => Some(*ani),
                            None => calculate_fastani(genomes[i], genomes[*rep])
                        }
                    };

                    if best_rep_min_ani.is_none() || fastani > best_rep_min_ani.unwrap() {
                        best_rep = Some(rep);
                        best_rep_min_ani = Some(fastani);
                    }
                }
            }
            to_return[rep_to_index[best_rep.unwrap()]].push(i);
        }
    }
    return to_return;
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
        let clusters = minhash_clusters(
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
        assert_eq!(
            vec![vec![0,1,2,3]],
            clusters
        )
    }

    #[test]
    fn test_minhash_fastani_two_clusters_same_ani() {
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
            Some(98.0),
        );
        assert_eq!(
            vec![vec![0,1,3],vec![2]],
            clusters
        )
    }

    #[test]
    fn test_minhash_fastani_two_clusters_low_minhash_ani() {
        init();
        let clusters = minhash_clusters(
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
        assert_eq!(
            vec![vec![0,1,3],vec![2]],
            clusters
        )
    }
}
