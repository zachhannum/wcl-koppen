//! The reader of a World Climate Lab export folder.
//!
//! An export folder has the manifest `climate_data.json`. Each grid is a
//! greyscale PNG with a JSON sidecar that gives its encoding.

use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use serde::Deserialize;

use crate::grid::Grid;

const FORMAT: &str = "wcl-climate-data";
const FORMAT_VERSION: u32 = 1;
const MANIFEST: &str = "climate_data.json";

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("{path}: {message}")]
    Io { path: PathBuf, message: String },
    #[error("{path}: {message}")]
    Json { path: PathBuf, message: String },
    #[error("{path}: the format is {found}, and this tool reads {FORMAT} version {FORMAT_VERSION}")]
    Format { path: PathBuf, found: String },
    #[error("the manifest has no grid {0}")]
    Missing(String),
    #[error("{path}: {message}")]
    Png { path: PathBuf, message: String },
    #[error("{path}: the grid is {found}, and the other grids are {expected}")]
    Size {
        path: PathBuf,
        found: String,
        expected: String,
    },
}

#[derive(Debug, Deserialize)]
struct Manifest {
    format: String,
    format_version: u32,
    files: Vec<Entry>,
}

#[derive(Debug, Deserialize)]
struct Entry {
    file: String,
    metadata: String,
}

#[derive(Debug, Deserialize)]
struct Sidecar {
    format: String,
    format_version: u32,
    grid: GridInfo,
    encoding: Encoding,
}

#[derive(Debug, Deserialize)]
struct GridInfo {
    width: usize,
    height: usize,
    longitude_left_edge: f64,
    longitude_right_edge: f64,
    latitude_top_edge: f64,
    latitude_bottom_edge: f64,
}

#[derive(Debug, Deserialize)]
#[serde(tag = "type", rename_all = "lowercase")]
enum Encoding {
    Linear {
        bits: u8,
        encoding_min: f64,
        encoding_max: f64,
    },
    Mask,
}

/// The grids that the classification reads.
#[derive(Clone, Debug)]
pub struct Climate {
    /// The monthly mean temperatures in °C, at the height of the ground of
    /// `elevation`, for months 1 to 12.
    pub temperature: Vec<Grid>,
    /// The monthly precipitation totals in mm, for months 1 to 12.
    pub precipitation: Vec<Grid>,
    /// The height of the ground in m. A sea cell is 0.
    pub elevation: Grid,
    /// 1 for land and 0 for sea.
    pub land: Grid,
}

/// A class map of the export, with the code of each class number.
#[derive(Clone, Debug)]
pub struct Classes {
    pub grid: Grid,
    /// The code of each class number, for example "Cfb", with its system.
    pub codes: Vec<(u32, String, String)>,
}

#[derive(Debug, Deserialize)]
struct ClassList {
    classes: Vec<ClassEntry>,
}

#[derive(Debug, Deserialize)]
struct ClassEntry {
    value: u32,
    system: String,
    code: String,
}

/// An export folder with its manifest.
pub struct Export {
    dir: PathBuf,
    entries: Vec<Entry>,
}

impl Export {
    pub fn open(dir: &Path) -> Result<Self, Error> {
        let path = dir.join(MANIFEST);
        let manifest: Manifest = read_json(&path)?;
        check_format(&path, &manifest.format, manifest.format_version)?;
        Ok(Self {
            dir: dir.to_path_buf(),
            entries: manifest.files,
        })
    }

    /// Reads the grid with the file stem `stem`, for example
    /// `temperature_month_01`.
    pub fn grid(&self, stem: &str) -> Result<Grid, Error> {
        let file = format!("{stem}.png");
        let entry = self
            .entries
            .iter()
            .find(|e| e.file == file)
            .ok_or_else(|| Error::Missing(stem.to_string()))?;
        let sidecar_path = self.dir.join(&entry.metadata);
        let sidecar: Sidecar = read_json(&sidecar_path)?;
        check_format(&sidecar_path, &sidecar.format, sidecar.format_version)?;
        let g = &sidecar.grid;
        let global = g.longitude_left_edge == -180.0
            && g.longitude_right_edge == 180.0
            && g.latitude_top_edge == 90.0
            && g.latitude_bottom_edge == -90.0;
        if !global {
            return Err(Error::Json {
                path: sidecar_path,
                message: "the grid does not cover the whole globe".to_string(),
            });
        }
        let path = self.dir.join(&entry.file);
        let (width, height, levels) = read_png(&path)?;
        if (width, height) != (g.width, g.height) {
            return Err(Error::Size {
                path,
                found: format!("{width} × {height}"),
                expected: format!("{} × {}", g.width, g.height),
            });
        }
        let values = match sidecar.encoding {
            Encoding::Linear {
                bits,
                encoding_min,
                encoding_max,
            } => {
                let top = ((1u32 << bits) - 1) as f64;
                let span = encoding_max - encoding_min;
                levels
                    .iter()
                    .map(|&l| (encoding_min + l as f64 / top * span) as f32)
                    .collect()
            }
            Encoding::Mask => levels
                .iter()
                .map(|&l| if l > 0 { 1.0 } else { 0.0 })
                .collect(),
        };
        Ok(Grid::new(width, height, values))
    }

    /// Reads the monthly temperatures and precipitation, the elevation, and
    /// the land mask. Each grid has the same size.
    pub fn climate(&self) -> Result<Climate, Error> {
        let months = |name: &str| -> Result<Vec<Grid>, Error> {
            (1..=12)
                .map(|m| self.grid(&format!("{name}_month_{m:02}")))
                .collect()
        };
        let climate = Climate {
            temperature: months("temperature")?,
            precipitation: months("precipitation")?,
            elevation: self.grid("elevation")?,
            land: self.grid("land_mask")?,
        };
        let size = (climate.land.width, climate.land.height);
        let all = climate
            .temperature
            .iter()
            .chain(&climate.precipitation)
            .chain([&climate.elevation]);
        for g in all {
            if (g.width, g.height) != size {
                return Err(Error::Size {
                    path: self.dir.clone(),
                    found: format!("{} × {}", g.width, g.height),
                    expected: format!("{} × {}", size.0, size.1),
                });
            }
        }
        Ok(climate)
    }

    /// Reads `climate_class.png` and the names of its classes.
    pub fn classes(&self) -> Result<Classes, Error> {
        let path = self.dir.join("climate_classes.json");
        let list: ClassList = read_json(&path)?;
        let (width, height, levels) = read_png(&self.dir.join("climate_class.png"))?;
        let values = levels.iter().map(|&l| l as f32).collect();
        Ok(Classes {
            grid: Grid::new(width, height, values),
            codes: list
                .classes
                .into_iter()
                .map(|c| (c.value, c.code, c.system))
                .collect(),
        })
    }
}

fn read_json<T: for<'de> Deserialize<'de>>(path: &Path) -> Result<T, Error> {
    let text = std::fs::read_to_string(path).map_err(|e| Error::Io {
        path: path.to_path_buf(),
        message: e.to_string(),
    })?;
    serde_json::from_str(&text).map_err(|e| Error::Json {
        path: path.to_path_buf(),
        message: e.to_string(),
    })
}

fn check_format(path: &Path, format: &str, version: u32) -> Result<(), Error> {
    if format == FORMAT && version == FORMAT_VERSION {
        Ok(())
    } else {
        Err(Error::Format {
            path: path.to_path_buf(),
            found: format!("{format} version {version}"),
        })
    }
}

/// Reads a greyscale PNG of 8 or 16 bits. Returns the width, the height, and
/// the levels.
fn read_png(path: &Path) -> Result<(usize, usize, Vec<u16>), Error> {
    let png_error = |message: String| Error::Png {
        path: path.to_path_buf(),
        message,
    };
    let file = File::open(path).map_err(|e| Error::Io {
        path: path.to_path_buf(),
        message: e.to_string(),
    })?;
    let mut reader = png::Decoder::new(BufReader::new(file))
        .read_info()
        .map_err(|e| png_error(e.to_string()))?;
    let size = reader
        .output_buffer_size()
        .ok_or_else(|| png_error("the image is too large".to_string()))?;
    let mut buf = vec![0; size];
    let info = reader
        .next_frame(&mut buf)
        .map_err(|e| png_error(e.to_string()))?;
    if info.color_type != png::ColorType::Grayscale {
        return Err(png_error(format!(
            "the color type is {:?}, not greyscale",
            info.color_type
        )));
    }
    let bytes = &buf[..info.buffer_size()];
    let levels = match info.bit_depth {
        png::BitDepth::Eight => bytes.iter().map(|&b| b as u16).collect(),
        png::BitDepth::Sixteen => bytes
            .as_chunks::<2>()
            .0
            .iter()
            .map(|&b| u16::from_be_bytes(b))
            .collect(),
        depth => {
            return Err(png_error(format!(
                "the bit depth is {depth:?}, not 8 or 16"
            )));
        }
    };
    Ok((info.width as usize, info.height as usize, levels))
}

#[cfg(test)]
mod tests {
    use super::*;

    const GRID: &str = r#"{"width": 2, "height": 1, "longitude_left_edge": -180.0,
        "longitude_right_edge": 180.0, "latitude_top_edge": 90.0, "latitude_bottom_edge": -90.0}"#;

    fn folder(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("wcl-koppen-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn write_png(path: &Path, depth: png::BitDepth, data: &[u8]) {
        let file = File::create(path).unwrap();
        let mut encoder = png::Encoder::new(file, 2, 1);
        encoder.set_color(png::ColorType::Grayscale);
        encoder.set_depth(depth);
        encoder
            .write_header()
            .unwrap()
            .write_image_data(data)
            .unwrap();
    }

    fn write(dir: &Path, name: &str, text: &str) {
        std::fs::write(dir.join(name), text).unwrap();
    }

    fn manifest(dir: &Path, stems: &[&str]) {
        let files: Vec<String> = stems
            .iter()
            .map(|s| format!(r#"{{"file": "{s}.png", "metadata": "{s}.json"}}"#))
            .collect();
        let text = format!(
            r#"{{"format": "wcl-climate-data", "format_version": 1, "files": [{}]}}"#,
            files.join(",")
        );
        write(dir, MANIFEST, &text);
    }

    #[test]
    fn linear_and_mask() {
        let dir = folder("grids");
        manifest(&dir, &["temperature_month_01", "land_mask"]);
        let linear = r#""encoding": {"type": "linear", "bits": 16, "encoding_min": -150.0, "encoding_max": 150.0}"#;
        write(
            &dir,
            "temperature_month_01.json",
            &format!(
                r#"{{"format": "wcl-climate-data", "format_version": 1, "grid": {GRID}, {linear}}}"#
            ),
        );
        write_png(
            &dir.join("temperature_month_01.png"),
            png::BitDepth::Sixteen,
            &[0, 0, 0xff, 0xff],
        );
        let mask = r#""encoding": {"type": "mask", "bits": 8}"#;
        write(
            &dir,
            "land_mask.json",
            &format!(
                r#"{{"format": "wcl-climate-data", "format_version": 1, "grid": {GRID}, {mask}}}"#
            ),
        );
        write_png(&dir.join("land_mask.png"), png::BitDepth::Eight, &[0, 255]);

        let export = Export::open(&dir).unwrap();
        let t = export.grid("temperature_month_01").unwrap();
        assert_eq!((t.width, t.height), (2, 1));
        assert_eq!(t.values, [-150.0, 150.0]);
        assert_eq!(export.grid("land_mask").unwrap().values, [0.0, 1.0]);
        assert!(matches!(export.grid("elevation"), Err(Error::Missing(_))));
    }

    #[test]
    fn unknown_version() {
        let dir = folder("version");
        write(
            &dir,
            MANIFEST,
            r#"{"format": "wcl-climate-data", "format_version": 2, "files": []}"#,
        );
        assert!(matches!(Export::open(&dir), Err(Error::Format { .. })));
    }
}
