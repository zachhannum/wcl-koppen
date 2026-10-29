//! The reader of the heightmap of the world.
//!
//! A heightmap is an equirectangular greyscale image of 16-bit levels over
//! the whole globe. A GeoTIFF from Vale holds the offset and the scale of its
//! levels in the GDAL metadata. A PNG holds no values, so the range comes
//! from the caller.

use std::path::{Path, PathBuf};

use tiff::decoder::{Decoder, DecodingResult, Limits};
use tiff::tags::Tag;

/// The tag of the GDAL metadata.
const GDAL_METADATA: Tag = Tag::Unknown(42112);

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("{path}: {message}")]
    Read { path: PathBuf, message: String },
    #[error("{0}: a PNG heightmap needs --height-range, because the file holds no values")]
    NoRange(PathBuf),
    #[error(
        "{path}: the image is {width} × {height}, and a world map is twice as wide as it is tall"
    )]
    Shape {
        path: PathBuf,
        width: usize,
        height: usize,
    },
}

/// A heightmap of 16-bit levels. Row 0 is at the north, and column 0 is at
/// the west edge, at -180°. Each level is the mean over its cell.
pub struct Heightmap {
    pub width: usize,
    pub height: usize,
    pub levels: Vec<u16>,
    /// The height of level 0, in m.
    pub offset: f64,
    /// The height of one level step, in m.
    pub scale: f64,
}

impl Heightmap {
    /// The height of a cell in m.
    pub fn get(&self, col: usize, row: usize) -> f32 {
        let level = self.levels[row * self.width + col];
        (self.offset + level as f64 * self.scale) as f32
    }

    /// The mean height of the block of `step` × `step` cells with its top
    /// left cell at (`col`, `row`).
    pub fn block_mean(&self, col: usize, row: usize, step: usize) -> f32 {
        let mut sum = 0u64;
        for r in row..row + step {
            let line = &self.levels[r * self.width + col..r * self.width + col + step];
            sum += line.iter().map(|&l| l as u64).sum::<u64>();
        }
        let mean = sum as f64 / (step * step) as f64;
        (self.offset + mean * self.scale) as f32
    }
}

/// Reads a GeoTIFF or a PNG heightmap. `range` gives the heights of levels
/// 0 and 65,535, and it replaces the values of a GeoTIFF.
pub fn read(path: &Path, range: Option<(f64, f64)>) -> Result<Heightmap, Error> {
    let error = |message: String| Error::Read {
        path: path.to_path_buf(),
        message,
    };
    let bytes = std::fs::read(path).map_err(|e| error(e.to_string()))?;
    let is_png = bytes.starts_with(b"\x89PNG");
    let (width, height, levels, stored) = if is_png {
        let (w, h, levels) = read_png(&bytes).map_err(error)?;
        (w, h, levels, None)
    } else {
        read_tiff(&bytes).map_err(error)?
    };
    if width != 2 * height {
        return Err(Error::Shape {
            path: path.to_path_buf(),
            width,
            height,
        });
    }
    let (offset, scale) = match (range, stored) {
        (Some((min, max)), _) => (min, (max - min) / 65535.0),
        (None, Some(values)) => values,
        (None, None) => return Err(Error::NoRange(path.to_path_buf())),
    };
    Ok(Heightmap {
        width,
        height,
        levels,
        offset,
        scale,
    })
}

fn read_png(bytes: &[u8]) -> Result<(usize, usize, Vec<u16>), String> {
    let mut decoder = png::Decoder::new(std::io::Cursor::new(bytes));
    decoder.set_ignore_text_chunk(true);
    let mut reader = decoder.read_info().map_err(|e| e.to_string())?;
    let info = reader.info();
    if info.color_type != png::ColorType::Grayscale || info.bit_depth != png::BitDepth::Sixteen {
        return Err("the image is not 16-bit greyscale".to_string());
    }
    let (width, height) = (info.width as usize, info.height as usize);
    let mut levels = Vec::with_capacity(width * height);
    let mut row = vec![0u8; width * 2];
    while let Some(r) = reader.next_row().map_err(|e| e.to_string())? {
        row.copy_from_slice(r.data());
        levels.extend(
            row.as_chunks::<2>()
                .0
                .iter()
                .map(|&b| u16::from_be_bytes(b)),
        );
    }
    Ok((width, height, levels))
}

type Tiff = (usize, usize, Vec<u16>, Option<(f64, f64)>);

fn read_tiff(bytes: &[u8]) -> Result<Tiff, String> {
    let file = std::io::Cursor::new(bytes);
    let mut decoder = Decoder::new(file)
        .map_err(|e| e.to_string())?
        .with_limits(Limits::unlimited());
    let (width, height) = decoder.dimensions().map_err(|e| e.to_string())?;
    let metadata = decoder
        .get_tag_ascii_string(GDAL_METADATA)
        .ok()
        .and_then(|xml| gdal_values(&xml));
    let levels = match decoder.read_image().map_err(|e| e.to_string())? {
        DecodingResult::U16(levels) => levels,
        _ => return Err("the image does not have 16-bit levels".to_string()),
    };
    Ok((width as usize, height as usize, levels, metadata))
}

/// The offset and the scale of the GDAL metadata.
fn gdal_values(xml: &str) -> Option<(f64, f64)> {
    let item = |role: &str| -> Option<f64> {
        let at = xml.find(&format!("role=\"{role}\""))?;
        let rest = &xml[at..];
        let start = rest.find('>')? + 1;
        let end = rest.find('<')?;
        rest[start..end].trim().parse().ok()
    };
    Some((item("offset")?, item("scale")?))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;

    #[test]
    fn gdal_metadata() {
        let xml = "<GDALMetadata>\n  <Item name=\"OFFSET\" sample=\"0\" role=\"offset\">-6000</Item>\n  <Item name=\"SCALE\" sample=\"0\" role=\"scale\">0.18310826276035706</Item>\n</GDALMetadata>\n";
        let (offset, scale) = gdal_values(xml).unwrap();
        assert_eq!(offset, -6000.0);
        assert!((scale * 65535.0 - 12000.0).abs() < 1e-6);
    }

    #[test]
    fn png_needs_a_range() {
        let dir = std::env::temp_dir().join(format!("wcl-koppen-height-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("h.png");
        let mut encoder = png::Encoder::new(File::create(&path).unwrap(), 2, 1);
        encoder.set_color(png::ColorType::Grayscale);
        encoder.set_depth(png::BitDepth::Sixteen);
        encoder
            .write_header()
            .unwrap()
            .write_image_data(&[0, 0, 0xff, 0xff])
            .unwrap();
        assert!(matches!(read(&path, None), Err(Error::NoRange(_))));
        let map = read(&path, Some((-6000.0, 6000.0))).unwrap();
        assert_eq!(map.get(0, 0), -6000.0);
        assert_eq!(map.get(1, 0), 6000.0);
    }
}
