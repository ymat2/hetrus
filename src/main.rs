use clap::Parser;
use rust_htslib::bcf::{Reader, Read, Record};
use std::path::PathBuf;
use std::convert::TryFrom;


#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    #[arg(short, long, value_hint = clap::ValueHint::FilePath, help = "VCF, bgzipped VCF, and BCF.")]
    input: PathBuf,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let mut bcf = Reader::from_path(args.input)?;

    for record_result in bcf.records() {
        let record = record_result?;
        let pos = record.pos();
        // number of sample in the vcf
        let sample_count = usize::try_from(record.sample_count()).unwrap();

        if !is_biallelic(&record) {
            // println!("Locus {} is not biallelic site.", pos);
            continue;
        }

        let gt_count = count_genotypes(&record, sample_count)?;
        // haploid が含まれてしまっているが、"all missing" を除外することでひとまず対処？
        let observed_het = calc_observed_heterozygosity(gt_count, true)?;
        println!("Locus: {}, GTs: {:?}, Ho: {:.4}", pos, gt_count, observed_het)
    }

    Ok(())
}


fn is_biallelic(record: &Record) -> bool {
    record.alleles().len() == 2
}


fn count_genotypes(record: &Record, n_samples: usize) -> Result<(u32, u32, u32), Box<dyn std::error::Error>> {
    let gts = record.genotypes()?;

    let mut n_homo = 0;
    let mut n_hetero = 0;
    let mut n_missing = 0;

    for sample_index in 0..n_samples {
        let mut ref_count = 0;
        let mut alt_count = 0;
        let mut missing_count = 0;

        for gta in gts.get(sample_index).iter() {
            match gta.index() {
                Some(0) => ref_count += 1,
                Some(_) => alt_count += 1,
                None => missing_count += 1,
            }
        }

        // diploid 前提
        match (ref_count, alt_count, missing_count) {
            (2, 0, 0) | (0, 2, 0) => n_homo += 1,
            (1, 1, 0) => n_hetero += 1,
            _ => n_missing += 1, // missing or partial missing
        }
    }

    Ok((n_homo, n_hetero, n_missing))
}

fn calc_observed_heterozygosity(gt_count: (u32, u32, u32), ignore_missing: bool) -> Result<f32, Box<dyn std::error::Error>> {
    let ho: f32;
    let (n_homo, n_hetero, n_missing) = gt_count;
    if ignore_missing {
        if n_homo + n_hetero == 0 {
            ho = f32::NAN
        } else {
            ho = (n_hetero) as f32 / (n_homo + n_hetero) as f32;
        }
    } else {
        ho = (n_hetero) as f32 / (n_homo + n_hetero + n_missing) as f32;
    }

    Ok(ho)
}
