use clap::Parser;
use rust_htslib::bcf::{Read, Reader, Record};
use std::collections::HashMap;
use std::convert::TryFrom;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    #[arg(short, long, value_hint = clap::ValueHint::FilePath, help = "VCF, bgzipped VCF, and BCF.")]
    input: PathBuf,
    #[arg(short, long, value_hint = clap::ValueHint::FilePath, help = "Output file path.")]
    out: PathBuf,
    #[arg(long = "window-size")]
    window_size: i64,
    #[arg(long = "window-step")]
    window_step: i64,
}

#[derive(Clone, Default, Debug)]
struct WindowBin {
    count: i64,
}
type ChromBins = Vec<WindowBin>;
type Bins = HashMap<String, ChromBins>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let mut bcf = Reader::from_path(args.input)?;
    let outfile = File::create(args.out)?;
    let mut out = BufWriter::new(outfile);
    let header = bcf.header().clone();
    let mut bins: Bins = HashMap::new();
    let window_size: i64 = args.window_size;
    let window_step: i64 = args.window_step;

    for record_result in bcf.records() {
        let record = record_result?;
        let rid = record.rid().expect("RID missing");
        let chr = std::str::from_utf8(header.rid2name(rid)?)?;
        let pos = record.pos() + 1;
        let sample_count = usize::try_from(record.sample_count()).unwrap();

        if !is_biallelic(&record) {
            continue;
        }

        let gt_count = count_genotypes(&record, sample_count)?;
        // haploid が含まれてしまっているが、"all missing" を除外することでひとまず対処？
        let observed_het = calc_observed_heterozygosity(gt_count, true)?;
        let expected_het = calc_expected_heterozygosity(gt_count, true)?;
        writeln!(
            out,
            "{}\t{}\t{:.4}\t{:.4}",
            chr, pos, observed_het, expected_het
        )?;

        add_snp_to_bins(&mut bins, chr, pos, window_size, window_step);
    }

    let mut chroms: Vec<_> = bins.keys().collect();
    chroms.sort();
    println!("CHROM\tBIN_START\tBIN_END\tN_VARIANTS");
    for chrom in chroms {
        let chrom_bins = &bins[chrom];
        for (idx, bin) in chrom_bins.iter().enumerate() {
            if bin.count == 0 {
                continue;
            }
            let start: i64 = idx as i64 * window_step + 1;
            let end: i64 = start + window_size - 1;
            println!("{}\t{}\t{}\t{}", chrom, start, end, bin.count);
        }
    }

    Ok(())
}

fn is_biallelic(record: &Record) -> bool {
    record.alleles().len() == 2
}

fn count_genotypes(
    record: &Record,
    n_samples: usize,
) -> Result<(i64, i64, i64, i64), Box<dyn std::error::Error>> {
    let gts = record.genotypes()?;

    let mut n_ref_homo = 0;
    let mut n_hetero = 0;
    let mut n_alt_homo = 0;
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

        match (ref_count, alt_count, missing_count) {
            (2, 0, 0) => n_ref_homo += 1,
            (0, 2, 0) => n_alt_homo += 1,
            (1, 1, 0) => n_hetero += 1,
            _ => n_missing += 1, // missing, partial missing, or not diploid
        }
    }

    Ok((n_ref_homo, n_hetero, n_alt_homo, n_missing))
}

fn calc_observed_heterozygosity(
    gt_count: (i64, i64, i64, i64),
    ignore_missing: bool,
) -> Result<f64, Box<dyn std::error::Error>> {
    let ho: f64;
    let (n_ref_homo, n_hetero, n_alt_homo, n_missing) = gt_count;
    if ignore_missing {
        if n_ref_homo + n_hetero + n_alt_homo == 0 {
            ho = f64::NAN
        } else {
            ho = (n_hetero) as f64 / (n_ref_homo + n_hetero + n_alt_homo) as f64;
        }
    } else {
        ho = (n_hetero) as f64 / (n_ref_homo + n_hetero + n_alt_homo + n_missing) as f64;
    }

    Ok(ho)
}

fn calc_expected_heterozygosity(
    gt_count: (i64, i64, i64, i64),
    ignore_missing: bool,
) -> Result<f64, Box<dyn std::error::Error>> {
    let af: f64;
    let he: f64;
    let (n_ref_homo, n_hetero, n_alt_homo, n_missing) = gt_count;
    if ignore_missing {
        if n_ref_homo + n_hetero + n_alt_homo == 0 {
            he = f64::NAN
        } else {
            af = (n_hetero + 2 * n_alt_homo) as f64
                / (2 * (n_ref_homo + n_hetero + n_alt_homo)) as f64;
            he = 2.0 * af * (1.0 - af)
        }
    } else {
        af = (n_hetero + 2 * n_alt_homo) as f64
            / (2 * (n_ref_homo + n_hetero + n_alt_homo + n_missing)) as f64;
        he = 2.0 * af * (1.0 - af)
    }

    Ok(he)
}

fn add_snp_to_bins(bins: &mut Bins, chrom: &str, pos: i64, size: i64, step: i64) {
    let mut first: usize = 0;
    if pos >= size {
        first = ((pos - size) as f64 / step as f64).ceil() as usize;
    };
    let last = (pos as f64 / step as f64).ceil() as usize;

    let chrom_bins = bins.entry(chrom.to_string()).or_default();
    if chrom_bins.len() < last {
        chrom_bins.resize_with(last, WindowBin::default);
    }

    for bin in chrom_bins.iter_mut().take(last).skip(first) {
        bin.count += 1;
    }
}
