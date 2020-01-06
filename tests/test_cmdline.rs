extern crate assert_cli;
extern crate tempfile;

#[cfg(test)]
mod tests {
    use assert_cli::Assert;
    extern crate tempfile;
 
    // #[test]
    // fn test_contig_kmer_hello_world() {
    //     Assert::main_binary()
    //         .with_args(&[
    //             "contig",
    //             "--single",
    //             "tests/data/reads_for_seq1_and_seq2.1.fq",
    //             "-r",
    //             "tests/data/7seqs.fna"])
    //         .succeeds()
    //         .stdout().is(
    //             "Sample	Contig	Coverage\n\
    //              tests/data/7seqs.fna/tests/data/reads_for_seq1_and_seq2.1.fq	genome1~random_sequence_length_11000	0\n\
    //              tests/data/7seqs.fna/tests/data/reads_for_seq1_and_seq2.1.fq	genome1~random_sequence_length_11010	0\n\
    //              tests/data/7seqs.fna/tests/data/reads_for_seq1_and_seq2.1.fq	genome2~seq1	0.6\n\
    //              tests/data/7seqs.fna/tests/data/reads_for_seq1_and_seq2.1.fq	genome3~random_sequence_length_11001	0\n\
    //              tests/data/7seqs.fna/tests/data/reads_for_seq1_and_seq2.1.fq	genome4~random_sequence_length_11002	0\n\
    //              tests/data/7seqs.fna/tests/data/reads_for_seq1_and_seq2.1.fq	genome5~seq2	0.6\n\
    //              tests/data/7seqs.fna/tests/data/reads_for_seq1_and_seq2.1.fq	genome6~random_sequence_length_11003	0\n")
    //          .unwrap();
    // }

    // #[test]
    // fn test_contig_kmer_different_ref_lengths() {
    //     Assert::main_binary()
    //         .with_args(&[
    //             "contig",
    //             "-r",
    //             "tests/data/2_single_species_dummy_dataset/2genomes_different_lengths.fna",
    //             "--single",
    //             "tests/data/2_single_species_dummy_dataset/reads/2genomes_2_reads.fq",
    //         ])
    //         .succeeds()
    //         .stdout().is(
    //             "Sample	Contig	Coverage\n\
    //              tests/data/2_single_species_dummy_dataset/2genomes_different_lengths.fna/tests/data/2_single_species_dummy_dataset/reads/2genomes_2_reads.fq	genome1	0.002844711575150617\n\
    //              tests/data/2_single_species_dummy_dataset/2genomes_different_lengths.fna/tests/data/2_single_species_dummy_dataset/reads/2genomes_2_reads.fq	genome2	0.37477019228323294\n")
    //         .unwrap()
    // }

    fn make_index(
        clades_file: &str,
        genome_definition_file: &str,
        extra_arguments: &Vec<&str>,
    ) -> tempfile::NamedTempFile {
        let tf = tempfile::NamedTempFile::new().unwrap();
        let mut args = vec![
            "index",
                "--index",
                tf.path().to_str().unwrap(),
                "--clades",
                clades_file,
                "--genome-definition",
                genome_definition_file];
        args.extend(extra_arguments.iter());

        Assert::main_binary()
            .with_args(&args)
            .succeeds()
            .unwrap();
        return tf;
    }


    #[test]
    fn test_bad_clade_file() {
        let tf = tempfile::NamedTempFile::new().unwrap();
        let clades_file = "tests/data/clade_file_reverse";
        let genome_definition_file = "tests/data/2_single_species_dummy_dataset/single_genome_example_tsv";
        Assert::main_binary()
            .with_args(&[
                "index",
                "--index",
                tf.path().to_str().unwrap(),
                "--clades",
                clades_file,
                "--genome-definition",
                genome_definition_file])
            .fails()
            .stderr()
            .contains("if that helps")
            .unwrap();
    }

    #[test]
    fn test_genome_kmer_one_genome() {
        let index = make_index(
            "tests/data/2_single_species_dummy_dataset/2genomes_same_genome.clades",
            "tests/data/2_single_species_dummy_dataset/single_genome_example_tsv", &vec![]
        );
        
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--single",
                "tests/data/2_single_species_dummy_dataset/reads/2genomes_2_reads.fq",
                "--index",
                index.path().to_str().unwrap()
                ])
            .succeeds()
            // TODO: The test value here is actually wrong - this seems to be a legitimate problem. Genome length is incorrect because repeats are only counted once.
            .stdout()
            .is(
                format!(
                    "Sample	Genome	Coverage\n\
                    {}/tests/data/2_single_species_dummy_dataset/reads/2genomes_2_reads.fq	g	0.2488986784140969\n",
                    index.path().to_str().unwrap())
                    .as_str())
            .unwrap()
    }

    #[test]
    fn test_genome_kmer_two_genomes_single_input() {
        let index = make_index(
            "tests/data/2_single_species_dummy_dataset/2genomes_different_lengths.different_clades",
            "tests/data/2_single_species_dummy_dataset/two_genomes_tsv", &vec![]
        );
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--single",
                "tests/data/2_single_species_dummy_dataset/reads/2genomes_2_reads.fq",
                "-d",
                index.path().to_str().unwrap()])
            .succeeds()
            .stdout().is(format!(
                "Sample	Genome	Coverage\n\
                 {}/tests/data/2_single_species_dummy_dataset/reads/2genomes_2_reads.fq	g1	0.0033675020949463837\n\
                 {}/tests/data/2_single_species_dummy_dataset/reads/2genomes_2_reads.fq	g2	0.40579044089961913\n",
                 index.path().to_str().unwrap(),
                 index.path().to_str().unwrap())
                 .as_str())
            .unwrap()
    }

    #[test]
    fn test_genome_kmer_two_genomes_paired_input() {
        let index = make_index(
            "tests/data/2_single_species_dummy_dataset/2genomes_different_lengths.different_clades",
            "tests/data/2_single_species_dummy_dataset/two_genomes_tsv", &vec![]
        );
        Assert::main_binary()
            .with_args(&[
                "genome",
                "-c",
                "tests/data/2_single_species_dummy_dataset/genome1_read_over100bp.1.fq",
                "tests/data/2_single_species_dummy_dataset/genome1_read_over100bp_one_bp_shorter.1.fq",
                "-d",
                index.path().to_str().unwrap()])
            .succeeds()
            .stdout().is(format!(
                "Sample	Genome	Coverage\n\
                 {}/tests/data/2_single_species_dummy_dataset/genome1_read_over100bp.1.fq	g1	0.7853107344632768\n\
                 {}/tests/data/2_single_species_dummy_dataset/genome1_read_over100bp.1.fq	g2	0\n",
                 index.path().to_str().unwrap(),
                 index.path().to_str().unwrap())
                 .as_str())
            .unwrap()
    }

    #[test]
    fn test_screen_positive_single() {
        Assert::main_binary()
            .with_args(&[
                "screen",
                "--mash-screen-db",
                "tests/data/screen/1/genome1_genome2.msh",
                "--single",
                "tests/data/screen/1/genome1.50paired_reads.fa.gz",
            ])
            .succeeds()
            .stdout().is(
                "genome1.fna	0.97991\n")
            .unwrap()
    }

    #[test]
    fn test_screen_negative() {
        Assert::main_binary()
            .with_args(&[
                "screen",
                "--mash-screen-db",
                "tests/data/screen/1/genome1_genome2.msh",
                "--single",
                "tests/data/screen/1/genome1.50paired_reads.fa.gz",
                "--min-identity",
                "0.999"
            ])
            .succeeds()
            .stdout().is("")
            .unwrap()
    }

    #[test]
    fn test_screen_positive_paired() {
        Assert::main_binary()
            .with_args(&[
                "screen",
                "--mash-screen-db",
                "tests/data/screen/1/genome1_genome2.msh",
                "--coupled",
                "tests/data/screen/1/genome1.1paired_reads.fa.gz",
                "tests/data/screen/1/genome1.another1paired_reads.fa.gz",
            ])
            .succeeds()
            .stdout().is(
                "genome1.fna	0.896575\n")
            .unwrap()
    }

    #[test]
    fn test_coverage_capping_hello_world() {
        let index = make_index(
            "tests/data/joel_test_data/2.clades",
            "tests/data/joel_test_data/2.definition", &vec![]
        );
        
        // Mapping a single read should work.
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--single",
                "tests/data/mapping_bugs/bad4.fq",
                "--index",
                index.path().to_str().unwrap()
                ])
            .succeeds()
            .stdout()
            .is(
                format!(
                    "Sample	Genome	Coverage\n\
                    {}/tests/data/mapping_bugs/bad4.fq	NC_023740.1	0.0021965791940018746\n\
                    {}/tests/data/mapping_bugs/bad4.fq	NC_028907.1	0\n",
                    index.path().to_str().unwrap(),
                    index.path().to_str().unwrap())
                    .as_str())
            .unwrap();

        // Mapping the same read pair twice should give the same result.
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--single",
                "tests/data/mapping_bugs/bad4_twice.fq",
                "--index",
                index.path().to_str().unwrap()
                ])
            .succeeds()
            .stdout()
            .is(
                format!(
                    "Sample	Genome	Coverage\n\
                    {}/tests/data/mapping_bugs/bad4_twice.fq	NC_023740.1	0.0021965791940018746\n\
                    {}/tests/data/mapping_bugs/bad4_twice.fq	NC_028907.1	0\n",
                    index.path().to_str().unwrap(),
                    index.path().to_str().unwrap())
                    .as_str())
            .unwrap()
    }

    #[test]
    fn test_no_core_genome() {
        let index = make_index(
            "tests/data/2_single_species_dummy_dataset/two_genomes_clades_same",
            "tests/data/2_single_species_dummy_dataset/two_genomes_tsv",
            &vec!["--no-core-genome"]
        );
        
        // Mapping a single read should work.
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--single",
                "tests/data/2_single_species_dummy_dataset/reads/2genomes_3_reads_one_not_core.fq",
                "--index",
                index.path().to_str().unwrap()
                ])
            .succeeds()
            .stdout()
            .is(
                format!(
                    "Sample	Genome	Coverage\n\
                    {}/tests/data/2_single_species_dummy_dataset/reads/2genomes_3_reads_one_not_core.fq	g1	0.003139807446870145\n\
                    {}/tests/data/2_single_species_dummy_dataset/reads/2genomes_3_reads_one_not_core.fq	g2	0.46730777646896743\n",
                    index.path().to_str().unwrap(),
                    index.path().to_str().unwrap())
                    .as_str())
            .unwrap();
    }
}


// TODO: Add mismatching bases test
// TODO: Filter fails when reference sequences are duplicated?
// TODO: Filter should spit things out if no thresholds are specified.
