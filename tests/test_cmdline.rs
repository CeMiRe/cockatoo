extern crate assert_cli;
extern crate tempfile;

#[cfg(test)]
mod tests {
    use assert_cli::Assert;
    extern crate tempfile;
 
    #[test]
    fn test_contig_kmer_hello_world() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "--single",
                "tests/data/reads_for_seq1_and_seq2.1.fq",
                "-r",
                "tests/data/7seqs.fna"])
            .succeeds()
            .stdout().is(
                "Sample	Contig	Coverage\n\
                 tests/data/7seqs.fna/tests/data/reads_for_seq1_and_seq2.1.fq	genome1~random_sequence_length_11000	0\n\
                 tests/data/7seqs.fna/tests/data/reads_for_seq1_and_seq2.1.fq	genome1~random_sequence_length_11010	0\n\
                 tests/data/7seqs.fna/tests/data/reads_for_seq1_and_seq2.1.fq	genome2~seq1	0.6\n\
                 tests/data/7seqs.fna/tests/data/reads_for_seq1_and_seq2.1.fq	genome3~random_sequence_length_11001	0\n\
                 tests/data/7seqs.fna/tests/data/reads_for_seq1_and_seq2.1.fq	genome4~random_sequence_length_11002	0\n\
                 tests/data/7seqs.fna/tests/data/reads_for_seq1_and_seq2.1.fq	genome5~seq2	0.6\n\
                 tests/data/7seqs.fna/tests/data/reads_for_seq1_and_seq2.1.fq	genome6~random_sequence_length_11003	0\n")
             .unwrap();
    }

    #[test]
    fn test_contig_kmer_different_ref_lengths() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "-r",
                "tests/data/2_single_species_dummy_dataset/2genomes_different_lengths.fna",
                "--single",
                "tests/data/2_single_species_dummy_dataset/reads/2genomes_2_reads.fq",
            ])
            .succeeds()
            .stdout().is(
                "Sample	Contig	Coverage\n\
                 tests/data/2_single_species_dummy_dataset/2genomes_different_lengths.fna/tests/data/2_single_species_dummy_dataset/reads/2genomes_2_reads.fq	genome1	0.002844711575150617\n\
                 tests/data/2_single_species_dummy_dataset/2genomes_different_lengths.fna/tests/data/2_single_species_dummy_dataset/reads/2genomes_2_reads.fq	genome2	0.37477019228323294\n")
            .unwrap()
    }

    #[test]
    fn test_genome_kmer_one_genome() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--single",
                "tests/data/2_single_species_dummy_dataset/reads/2genomes_2_reads.fq",
                "--clades",
                "tests/data/2_single_species_dummy_dataset/2genomes_same_genome.clades",
                "-r",
                "tests/data/2_single_species_dummy_dataset/2genomes_different_lengths.fna",
                "--genome-definition",
                "tests/data/2_single_species_dummy_dataset/single_genome_example_tsv"])
            .succeeds()
            // This seems to be a legitimate problem. Genome length is incorrect, maybe because repeats are only counted once?
            .stdout().is(
                "Sample	Genome	Coverage\n\
                 tests/data/2_single_species_dummy_dataset/2genomes_different_lengths.fna/tests/data/2_single_species_dummy_dataset/reads/2genomes_2_reads.fq	g	0.226\n")
            .unwrap()
    }

    #[test]
    fn test_genome_kmer_two_genomes_single_input() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--single",
                "tests/data/2_single_species_dummy_dataset/reads/2genomes_2_reads.fq",
                "-r",
                "tests/data/2_single_species_dummy_dataset/2genomes_different_lengths.fna",
                "--genome-definition",
                "tests/data/2_single_species_dummy_dataset/two_genomes_tsv",
                "--clades",
                "tests/data/2_single_species_dummy_dataset/2genomes_different_lengths.different_clades",
                ])
            .succeeds()
            .stdout().is(
                "Sample	Genome	Coverage\n\
                 tests/data/2_single_species_dummy_dataset/2genomes_different_lengths.fna/tests/data/2_single_species_dummy_dataset/reads/2genomes_2_reads.fq	g1	0.002844711575150617\n\
                 tests/data/2_single_species_dummy_dataset/2genomes_different_lengths.fna/tests/data/2_single_species_dummy_dataset/reads/2genomes_2_reads.fq	g2	0.37477019228323294\n")
            .unwrap()
    }

    #[test]
    fn test_genome_kmer_two_genomes_paired_input() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "-c",
                "tests/data/2_single_species_dummy_dataset/genome1_read_over100bp.1.fq",
                "tests/data/2_single_species_dummy_dataset/genome1_read_over100bp_one_bp_shorter.1.fq",
                "-r",
                "tests/data/2_single_species_dummy_dataset/2genomes_different_lengths.fna",
                "--genome-definition",
                "tests/data/2_single_species_dummy_dataset/two_genomes_tsv"])
            .succeeds()
            .stdout().is(
                "Sample	Genome	Coverage\n\
                 tests/data/2_single_species_dummy_dataset/2genomes_different_lengths.fna/tests/data/2_single_species_dummy_dataset/genome1_read_over100bp.1.fq	g1	0.35\n\
                 tests/data/2_single_species_dummy_dataset/2genomes_different_lengths.fna/tests/data/2_single_species_dummy_dataset/genome1_read_over100bp.1.fq	g2	0\n")
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

}


// TODO: Add mismatching bases test
// TODO: Filter fails when reference sequences are duplicated?
// TODO: Filter should spit things out if no thresholds are specified.
