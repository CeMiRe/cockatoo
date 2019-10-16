// NOTE: Much of this code is a direct copy from coverm, but there are some differences.

use std::process;

use tempfile::NamedTempFile;

#[derive(Clone)]
pub enum ReadFormat {
    Coupled,
    Interleaved,
    Single,
}

pub struct MappingParameters<'a> {
    references: Vec<&'a str>,
    threads: u16,
    read1: Vec<&'a str>,
    read2: Vec<&'a str>,
    interleaved: Vec<&'a str>,
    unpaired: Vec<&'a str>,
    iter_reference_index: usize,
}

impl<'a> MappingParameters<'a> {
    pub fn generate_from_clap(
        m: &'a clap::ArgMatches,
        reference_tempfile: &'a Option<NamedTempFile>)
        -> MappingParameters<'a> {

        let mut read1: Vec<&str> = vec!();
        let mut read2: Vec<&str> = vec!();
        let mut interleaved: Vec<&str> = vec!();
        let mut unpaired: Vec<&str> = vec!();

        if m.is_present("read1") {
            read1 = m.values_of("read1").unwrap().collect();
            read2 = m.values_of("read2").unwrap().collect();
            if read1.len() != read2.len() {
                error!("When specifying paired reads with the -1 and -2 flags, \
                        there must be equal numbers specified. Instead found \
                        {} and {} respectively", read1.len(), read2.len());
                process::exit(1);
            }
        }

        // Parse --coupled
        if m.is_present("coupled") {
            let coupled: Vec<&str> = m.values_of("coupled").unwrap().collect();
            if coupled.len() % 2 != 0 {
                error!(
                    "The --coupled flag must be set with pairs of read \
                     sets, but an odd number ({}) was specified",
                    coupled.len()
                );
                process::exit(1);
            }
            let mut i = 0;
            while i < coupled.len() {
                read1.push(coupled[i]);
                read2.push(coupled[i+1]);
                i += 2;
            }
        }

        if m.is_present("interleaved") {
            interleaved = m.values_of("interleaved").unwrap().collect();
        }
        if m.is_present("single") {
            unpaired = m.values_of("single").unwrap().collect();
        }

        return MappingParameters {
            references: match reference_tempfile {
                Some(r) => vec!(r.path().to_str().unwrap()),
                None => match m.values_of("reference") {
                    Some(refs) => refs.collect(),
                    None => vec![]
                }
            },
            threads: m.value_of("threads").unwrap().parse::<u16>()
                .expect("Failed to convert threads argument into integer"),
            read1: read1,
            read2: read2,
            interleaved: interleaved,
            unpaired: unpaired,
            iter_reference_index: 0,
        }
    }

    // Return a Vec of str + Option<str> where each entry is a read pair or
    // single with None.
    pub fn readsets(&self) -> Vec<(&str, Option<&str>)> {
        let mut to_return: Vec<(&str, Option<&str>)> = vec!();

        for (ref r1, ref r2) in self.read1.iter().zip(self.read2.iter()) {
            to_return.push((r1, Some(r2)))
        }
        for ref s in self.unpaired.iter() {
            to_return.push((&s, None))
        }
        return to_return
    }
}

pub struct SingleReferenceMappingParameters<'a> {
    pub reference: &'a str,
    threads: u16,
    read1: Vec<&'a str>,
    read2: Vec<&'a str>,
    interleaved: Vec<&'a str>,
    unpaired: Vec<&'a str>,

    iter_read_pair_index: usize,
    iter_interleaved_index: usize,
    iter_unpaired_index: usize,
}

impl<'a> SingleReferenceMappingParameters<'a> {
    pub fn len(&self) -> usize {
        self.read1.len() + self.interleaved.len() + self.unpaired.len()
    }
}

impl<'a> Iterator for MappingParameters<'a> {
    type Item = SingleReferenceMappingParameters<'a>;

    fn next(&mut self) -> Option<SingleReferenceMappingParameters<'a>> {
        if self.iter_reference_index < self.references.len() {
            let i = self.iter_reference_index;
            self.iter_reference_index += 1;
            return Some(SingleReferenceMappingParameters {
                reference: self.references[i],
                threads: self.threads,
                read1: self.read1.clone(),
                read2: self.read2.clone(),
                interleaved: self.interleaved.clone(),
                unpaired: self.unpaired.clone(),
                iter_read_pair_index: 0,
                iter_interleaved_index: 0,
                iter_unpaired_index: 0,
            })
        } else {
            return None
        }
    }
}

impl<'a> Iterator for SingleReferenceMappingParameters<'a> {
    type Item = OneSampleMappingParameters<'a>;

    fn next(&mut self) -> Option<OneSampleMappingParameters<'a>> {
        if self.iter_read_pair_index < self.read1.len() {
            let i = self.iter_read_pair_index;
            self.iter_read_pair_index += 1;
            return Some(OneSampleMappingParameters {
                reference: self.reference,
                read_format: ReadFormat::Coupled,
                read1: self.read1[i],
                read2: Some(self.read2[i]),
                threads: self.threads,
            })
        } else if self.iter_interleaved_index < self.interleaved.len() {
            let i = self.iter_interleaved_index;
            self.iter_interleaved_index += 1;
            return Some(OneSampleMappingParameters {
                reference: self.reference,
                read_format: ReadFormat::Interleaved,
                read1: self.interleaved[i],
                read2: None,
                threads: self.threads,
            })
        } else if self.iter_unpaired_index < self.unpaired.len() {
            let i = self.iter_unpaired_index;
            self.iter_unpaired_index += 1;
            return Some(OneSampleMappingParameters {
                reference: self.reference,
                read_format: ReadFormat::Single,
                read1: self.unpaired[i],
                read2: None,
                threads: self.threads,
            })
        } else {
            return None
        }
    }
}

pub struct OneSampleMappingParameters<'a> {
    pub reference: &'a str,
    pub read_format: ReadFormat,
    pub read1: &'a str,
    pub read2: Option<&'a str>,
    pub threads: u16,
}
