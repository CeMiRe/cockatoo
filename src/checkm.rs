use std;
use csv;
use log::*;

pub struct CheckMTabTable {
}

impl CheckMTabTable {
    pub fn good_quality_genome_names(file_path: &str, min_completeness: f32, max_contamination: f32) -> Vec<String> {
        let rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path(std::path::Path::new(file_path));

        let mut passes = vec![];
        let mut total_seen = 0usize;
        for result in rdr
            .expect(&format!("Failed to parse CheckM tab table {}", file_path))
            .records() {

            let res = result.expect("Parsing error in CheckM tab table file");
            if res.len() != 14 {
                error!("Parsing error in CheckM tab table file - didn't find 13 columns in line {:?}", res);
                std::process::exit(1);
            }
            let completeness: f32 = res[11].parse::<f32>().expect("Error parsing completeness in checkm tab table");
            let contamination: f32 = res[12].parse::<f32>().expect("Error parsing contamination in checkm tab table");
            trace!("For {}, found completeness {} and contamination {}", &res[0], completeness, contamination);
            if completeness >= min_completeness && contamination <= max_contamination {
                passes.push(res[0].to_string());
            }
            total_seen += 1;
        }
        debug!("Read in {} genomes from {}, {} passed the quality thresholds", total_seen, file_path, passes.len());
        return passes;
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }
    
    #[test]
    fn test_good_quality_genome_names() {
        init();
        assert_eq!(
            vec!["GUT_GENOME006390.gff","GUT_GENOME011264.gff","GUT_GENOME011296.gff"], 
            CheckMTabTable::good_quality_genome_names(&"tests/data/checkm.tsv",56.,1.))
    }
}