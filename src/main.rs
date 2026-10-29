use std::collections::BTreeMap;
use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result, bail};
use clap::{Args as ClapArgs, Parser, Subcommand};
use wcl_koppen::downscale::{self, Downscale};
use wcl_koppen::koppen::{self, CLASSES, Rules};
use wcl_koppen::trace::{self, Shape};
use wcl_koppen::wcl::Export;
use wcl_koppen::{NO_CLASS, classify_grid, filter, geojson, height};

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
    /// Writes the climate zones at the detail of a heightmap as GeoJSON.
    Zones(Zones),
}

#[derive(ClapArgs)]
struct Zones {
    /// The export folder, with `climate_data.json`.
    export: PathBuf,
    /// The heightmap: a GeoTIFF from Vale, or a 16-bit greyscale PNG. Each
    /// is an equirectangular map of the whole globe. Land is above 0 m.
    #[arg(long)]
    heightmap: PathBuf,
    /// The heights of levels 0 and 65,535 in m, for example -6000,6000. A
    /// PNG needs it, and for a GeoTIFF it replaces the values of the file.
    #[arg(long, value_name = "MIN,MAX", allow_hyphen_values = true)]
    height_range: Option<String>,
    /// The GeoJSON file to write.
    #[arg(long)]
    out: PathBuf,
    /// The cells of the heightmap along the edge of one cell of the zones.
    #[arg(long, default_value_t = 1)]
    step: usize,
    /// The fall of the temperature with height, in °C per km.
    #[arg(long, default_value_t = 6.5)]
    lapse_rate: f32,
    /// The temperature of the coldest month in °C that separates C from D.
    #[arg(long, default_value_t = 0.0, allow_hyphen_values = true)]
    cd_boundary: f32,
    /// The smallest zone in km². A smaller zone joins the zones around it.
    /// The default is 25 cells at the equator.
    #[arg(long)]
    min_area_km2: Option<f64>,
    /// The width of the smoothing of the boundaries, in km. The default is
    /// one cell at the equator.
    #[arg(long)]
    smooth_km: Option<f64>,
    /// The largest distance of a removed point from the simpler line, in
    /// km. The default is a quarter of a cell at the equator.
    #[arg(long)]
    tolerance_km: Option<f64>,
    /// The radius of the world in km.
    #[arg(long, default_value_t = 6371.0)]
    radius_km: f64,
}

fn main() -> Result<()> {
    match Args::parse().command {
        Command::Check {
            export,
            cd_boundary,
        } => check(&export, Rules { cd_boundary }),
        Command::Zones(args) => zones(&args),
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

fn parse_range(text: &str) -> Result<(f64, f64)> {
    let (a, b) = text
        .split_once(',')
        .context("the height range is two numbers with a comma, for example -6000,6000")?;
    let (min, max): (f64, f64) = (a.trim().parse()?, b.trim().parse()?);
    if min.partial_cmp(&max) != Some(std::cmp::Ordering::Less) {
        bail!("the first height of the range must be below the second");
    }
    Ok((min, max))
}

fn zones(args: &Zones) -> Result<()> {
    let start = Instant::now();
    let export = Export::open(&args.export)?;
    let climate = export.climate()?;
    let range = args.height_range.as_deref().map(parse_range).transpose()?;
    let heights = height::read(&args.heightmap, range)?;
    let step = args.step.max(1);
    let (width, rows) = (heights.width / step, heights.height / step);
    eprintln!(
        "Read the export and a heightmap of {} × {} in {:.1} s",
        heights.width,
        heights.height,
        start.elapsed().as_secs_f64()
    );

    let slope = downscale::elevation_slope(&climate, &heights);
    if slope.is_finite() {
        eprintln!("The elevation of the export is {slope:.3} times the heightmap");
        if (slope - 1.0).abs() > 0.1 {
            eprintln!(
                "Warning: the two heights do not have the same scale. The export \
                 used other mountains, so the correction is large. Check the height \
                 range, or give World Climate Lab the same peak height."
            );
        }
    }

    let cell_km = std::f64::consts::TAU * args.radius_km / width as f64;
    let min_area = args.min_area_km2.unwrap_or(25.0 * cell_km * cell_km);
    let shape = Shape {
        smooth_km: args.smooth_km.unwrap_or(cell_km),
        tolerance_km: args.tolerance_km.unwrap_or(cell_km / 4.0),
        radius_km: args.radius_km,
    };

    let t = Instant::now();
    let params = Downscale {
        lapse_rate: args.lapse_rate,
        step,
        rules: Rules {
            cd_boundary: args.cd_boundary,
        },
    };
    let mut grid = downscale::classify(&climate, &heights, params);
    drop(heights);
    eprintln!(
        "Classified {width} × {rows} cells of {cell_km:.2} km in {:.1} s",
        t.elapsed().as_secs_f64()
    );

    let t = Instant::now();
    let merged = filter::merge_small(&mut grid, min_area, args.radius_km);
    eprintln!(
        "Merged {merged} zones below {min_area:.0} km² in {:.1} s",
        t.elapsed().as_secs_f64()
    );

    let t = Instant::now();
    let polygons = trace::polygons(&grid, shape);
    let points: usize = polygons
        .iter()
        .map(|p| p.exterior.len() + p.holes.iter().map(Vec::len).sum::<usize>())
        .sum();
    eprintln!(
        "Traced {} polygons with {points} points in {:.1} s",
        polygons.len(),
        t.elapsed().as_secs_f64()
    );

    let file = File::create(&args.out).with_context(|| format!("{}", args.out.display()))?;
    let mut out = BufWriter::new(file);
    geojson::write(&mut out, &polygons, args.radius_km)?;
    eprintln!(
        "Wrote {} in {:.1} s in all",
        args.out.display(),
        start.elapsed().as_secs_f64()
    );
    Ok(())
}
