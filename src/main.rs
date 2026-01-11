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
    window_size: u64,
    #[arg(long = "window-step")]
    window_step: u64,
}

#[derive(Clone, Default, Debug)]
struct WindowBin {
    count: u64,
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
    let window_size: u64 = args.window_size;
    let window_step: u64 = args.window_step;

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
    for chrom in chroms {
        let chrom_bins = &bins[chrom];
        for (i, bin) in chrom_bins.iter().enumerate() {
            if bin.count == 0 {
                continue;
            }
            let start: u64 = i as u64 * window_step + 1;
            let end: u64 = start + window_size - 1;
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
) -> Result<(u32, u32, u32, u32), Box<dyn std::error::Error>> {
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
    gt_count: (u32, u32, u32, u32),
    ignore_missing: bool,
) -> Result<f32, Box<dyn std::error::Error>> {
    let ho: f32;
    let (n_ref_homo, n_hetero, n_alt_homo, n_missing) = gt_count;
    if ignore_missing {
        if n_ref_homo + n_hetero + n_alt_homo == 0 {
            ho = f32::NAN
        } else {
            ho = (n_hetero) as f32 / (n_ref_homo + n_hetero + n_alt_homo) as f32;
        }
    } else {
        ho = (n_hetero) as f32 / (n_ref_homo + n_hetero + n_alt_homo + n_missing) as f32;
    }

    Ok(ho)
}

fn calc_expected_heterozygosity(
    gt_count: (u32, u32, u32, u32),
    ignore_missing: bool,
) -> Result<f32, Box<dyn std::error::Error>> {
    let af: f32;
    let he: f32;
    let (n_ref_homo, n_hetero, n_alt_homo, n_missing) = gt_count;
    if ignore_missing {
        if n_ref_homo + n_hetero + n_alt_homo == 0 {
            he = f32::NAN
        } else {
            af = (n_hetero + 2 * n_alt_homo) as f32
                / (2 * (n_ref_homo + n_hetero + n_alt_homo)) as f32;
            he = 2.0 * af * (1.0 - af)
        }
    } else {
        af = (n_hetero + 2 * n_alt_homo) as f32
            / (2 * (n_ref_homo + n_hetero + n_alt_homo + n_missing)) as f32;
        he = 2.0 * af * (1.0 - af)
    }

    Ok(he)
}

fn window_range(pos: i64, size: u64, step: u64) -> (usize, usize) {
    let mut first = ((pos - size as i64) as f64 / step as f64).ceil() as isize;
    if first < 0 {
        first = 0;
    }
    let last = (pos as f64 / step as f64).ceil() as usize;
    (first as usize, last)
}

fn add_snp_to_bins(bins: &mut Bins, chrom: &str, pos: i64, size: u64, step: u64) {
    let (first, last) = window_range(pos, size, step);
    let chrom_bins = bins.entry(chrom.to_string()).or_default();

    if chrom_bins.len() < last {
        chrom_bins.resize_with(last, WindowBin::default);
    }

    for idx in first..last {
        let bin = &mut chrom_bins[idx];
        bin.count += 1;
    }
}
