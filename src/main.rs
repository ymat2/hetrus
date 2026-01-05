use clap::Parser;
use rust_htslib::bcf::{Reader, Read};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    #[arg(short, long, help = "Provide your name.", default_value = "world")]
    name: String,

    #[args(short, long, value_hint = clap::ValueHint::FilePath)]
    input: PathBuf,
}

fn main() {
    let args = Args::parse();
    println!("Hello, {}!", args.name);
    let mut bcf = Reader::from_path(args.input).expect("Error opening file.");

    for record in bcf.records() {
        let record = record?;
        let pos = record.pos();
        println!("POS: {}", pos);
    }
}
