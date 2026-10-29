use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use anyhow::Result;
use clap::{Parser, Subcommand};
use wcl_koppen::koppen::{self, CLASSES, Rules};
use wcl_koppen::wcl::Export;
use wcl_koppen::{NO_CLASS, classify_grid};

/// Makes Köppen-Geiger climate zones from a World Climate Lab export.
#[derive(Parser)]
#[command(version)]
struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Classifies the export on its own grid and compares the result with
    /// its `climate_class.png`.
    Check {
        /// The export folder, with `climate_data.json`.
        export: PathBuf,
        /// The temperature of the coldest month in °C that separates C from D.
        #[arg(long, default_value_t = 0.0, allow_hyphen_values = true)]
        cd_boundary: f32,
    },
}

fn main() -> Result<()> {
    match Args::parse().command {
        Command::Check {
            export,
            cd_boundary,
        } => check(&export, Rules { cd_boundary }),
    }
}

fn check(dir: &Path, rules: Rules) -> Result<()> {
    let export = Export::open(dir)?;
    let climate = export.climate()?;
    let reference = export.classes()?;
    let ours = classify_grid(&climate, rules);

    // The class of each Köppen-Geiger class number of the export.
    let theirs: BTreeMap<u32, koppen::Class> = reference
        .codes
        .iter()
        .filter(|(_, _, system)| system == "Koppen-Geiger")
        .filter_map(|(value, code, _)| koppen::class(code).map(|c| (*value, c)))
        .collect();

    let mut total = 0usize;
    let mut same = 0usize;
    let mut other_system = 0usize;
    let mut confusion: BTreeMap<(koppen::Class, koppen::Class), usize> = BTreeMap::new();
    for (i, &c) in ours.iter().enumerate() {
        if c == NO_CLASS {
            continue;
        }
        let Some(&t) = theirs.get(&(reference.grid.values[i] as u32)) else {
            other_system += 1;
            continue;
        };
        total += 1;
        if c == t {
            same += 1;
        } else {
            *confusion.entry((t, c)).or_default() += 1;
        }
    }
    println!(
        "{same} of {total} land cells match ({:.2} %)",
        100.0 * same as f64 / total.max(1) as f64
    );
    if other_system > 0 {
        println!("{other_system} land cells have a class of another system in the export");
    }
    let mut worst: Vec<_> = confusion.into_iter().collect();
    worst.sort_by_key(|&(_, n)| std::cmp::Reverse(n));
    for ((t, c), n) in worst.into_iter().take(15) {
        println!(
            "{n:>8}  export {:<4} here {}",
            CLASSES[t as usize].0, CLASSES[c as usize].0
        );
    }
    Ok(())
}
