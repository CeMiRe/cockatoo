// use std::path::Path;

// use pseudoaligner::*;
// use debruijn::Kmer;

// use seq_io::fastq;

// use log::Level;

// use pseudoalignment_reference_readers::DebruijnIndex;


// pub fn calculate_contig_kmer_coverage<K: Kmer + Sync + Send>(
//     forward_fastq: &str,
//     reverse_fastq: Option<&str>,
//     num_threads: usize,
//     print_zero_coverage_contigs: bool,
//     index: &DebruijnIndex<K>)
// -> Vec<(usize, f64)>{

//     // Do the mappings
//     let reads = fastq::Reader::from_maybe_gzip_path(Path::new(forward_fastq))
//         .expect(&format!("Failure to read file {}", forward_fastq));
//     let reverse_reads = match reverse_fastq {
//         Some(s) => Some(
//             fastq::Reader::from_maybe_gzip_path(Path::new(s))
//                 .expect(&format!("Failure to read reverse read file {}", s))),
//         None => None
//     };
//     let (eq_class_indices, eq_class_coverages, _eq_class_counts) = pseudoaligner::process_reads::<K, pseudoaligner::Pseudoaligner<K>>(
//         reads, reverse_reads, &index.index, num_threads)
//         .expect("Failure during mapping process");
//     info!("Finished mapping reads!");

//     if log_enabled!(Level::Debug) {
//         for (eq_class, index) in &eq_class_indices {
//             debug!("eq_class_index: {:?}\t{}", eq_class, index);
//         };
//         for (i, eq_class_coverage) in eq_class_coverages.iter().enumerate() {
//             debug!("eq_class_coverage {} {}", i, eq_class_coverage);
//         }
//     }

//     let kmer_coverage_total: f64 = eq_class_coverages.iter().sum::<usize>() as f64;

//     // Make a new vec containing the eq_class memberships
//     let mut eq_classes = vec![None; eq_class_coverages.len()];
//     for (eq_class, index) in &eq_class_indices {
//         eq_classes[*index] = Some(eq_class)
//     }

//     // EM algorithm
//     // Randomly assign random relative abundances to each of the genomes
//     // TODO: remove: for dev, assign each the same abundance
//     info!("Starting EM process ..");
//     let mut contig_to_relative_abundance = vec![1.0; index.seq_lengths.len()];
//     let mut contig_to_read_count;
//     let mut num_covergence_steps: u32 = 0;

//     loop { // loop until converged
//         // E-step: Determine the number of reads we expect come from each contig
//         // given their relative abundance.

//         // for each equivalence class / count pair, we expect for genome A the
//         // abundance of A divided by the sum of abundances of genomes in the
//         // equivalence class.
//         contig_to_read_count = vec![0.0; index.seq_lengths.len()];
//         for (i, coverage) in eq_class_coverages.iter().enumerate() {
//             match eq_classes[i] {
//                 None => unreachable!(),
//                 Some(ref eqs) => {
//                     // Find coverages of each contig to add
//                     let relative_abundances: Vec<f64> = eqs.iter().map(
//                         |eq|
//                         // TODO: Remove the 'as usize' by making eq a usize throughout
//                         contig_to_relative_abundance[*eq as usize]).collect();

//                     // Add coverages divided by the sum of relative abundances
//                     let total_abundance_of_matching_contigs: f64 = relative_abundances.iter().sum();
//                     for (eq, relabund) in eqs.iter().zip(relative_abundances.iter()) {
//                         contig_to_read_count[*eq as usize] += (*coverage as f64) *
//                             relabund / total_abundance_of_matching_contigs;
//                     }
//                 }
//             }
//         }
//         debug!("After E-step have contig read counts: {:?}", contig_to_read_count);

//         // M-step: Work out the relative abundance given the number of reads
//         // predicted in the E-step. Relative abundance is just the ratio of the
//         // coverages, weighted by the inverse of each contig's length.
//         //
//         // Or, for genome mode, weighted by the inverse of the sum of the
//         // genome's length.
//         //TODO: zip here instead?
//         // First determine the total scaled abundance
//         let mut total_scaling_abundance: f64 = 0.0;
//         let mut converge = true;

//         for (i, read_count) in contig_to_read_count.iter().enumerate() {
//             total_scaling_abundance += read_count / (index.seq_lengths[i] as f64)
//         }

//         // Next set the abundances so the total == 1.0
//         //
//         // Converge when all contigs with total abundance > 0.01 * kmer_coverage_total
//         // change abundance by < 1%.
//         for (i, read_count) in contig_to_read_count.iter().enumerate() {
//             let to_add = read_count / (index.seq_lengths[i] as f64);
//             let new_relabund = to_add / total_scaling_abundance;
//             if new_relabund * kmer_coverage_total > 0.01*100.0 { // Add 100 in there
//                 // as a rough mapped kmers per aligning fragment.

//                 // Enough abundance that this contig might stop convergence if it
//                 // changed by enough.
//                 let delta = new_relabund / contig_to_relative_abundance[i];
//                 debug!("For testing convergence of index {}, found delta {}", i, delta);
//                 if delta < 0.99 || delta > 1.01 {
//                     debug!("No converge for you");
//                     converge = false;
//                 }
//             }
//             contig_to_relative_abundance[i] = new_relabund;
//         }
//         debug!("At end of M-step, have relative abundances: {:?}", contig_to_relative_abundance);

//         num_covergence_steps += 1;
//         if converge {
//             info!("EM process converged after {} steps", num_covergence_steps);
//             break;
//         }
//     }

//     // Print results
//     let mut to_return = vec!();
//     for (i, total_coverage) in contig_to_read_count.iter().enumerate() {
//         if print_zero_coverage_contigs || *total_coverage > 0.0 {
//             // Print average coverage as total coverage divided by contig length.
//             to_return.push((i, total_coverage / (index.seq_lengths[i] as f64)));
//         }
//     }
//     return to_return;
// }

pub struct PseudoalignmentReadInput {
    pub forward_fastq: String,
    pub reverse_fastq: Option<String>,
    pub sample_name: String
}

// pub fn calculate_and_print_contig_kmer_coverages<K: Kmer + Send + Sync>(
//     read_inputs: &Vec<PseudoalignmentReadInput>,
//     num_threads: usize,
//     print_zero_coverage_contigs: bool,
//     index: &DebruijnIndex<K>) {

//     println!("Sample\tContig\tCoverage");
//     for read_input in read_inputs {
//         let covs = calculate_contig_kmer_coverage(
//             &read_input.forward_fastq,
//             match read_input.reverse_fastq {
//                 Some(ref s) => Some(&s),
//                 None => None
//             },
//             num_threads,
//             print_zero_coverage_contigs,
//             index);

//         for contig_res in covs {
//             println!(
//                 "{}\t{}\t{}",
//                 read_input.sample_name,
//                 index.tx_names[contig_res.0],
//                 contig_res.1);
//         }
//         info!("Finished printing contig coverages for sample {}",
//               read_input.sample_name);
//     }
//     info!("Finished printing contig coverages");
// }
