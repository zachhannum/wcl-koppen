//! Classification on the grid of a finer heightmap.

use std::collections::VecDeque;

use rayon::prelude::*;

use crate::NO_CLASS;
use crate::grid::Grid;
use crate::height::Heightmap;
use crate::koppen::{self, Class, Rules};
use crate::wcl::Climate;

/// A class for each cell of an equirectangular grid. Row 0 is at the north.
/// A sea cell holds [`NO_CLASS`].
#[derive(Clone, Debug, PartialEq)]
pub struct ClassGrid {
    pub width: usize,
    pub height: usize,
    pub classes: Vec<Class>,
}

impl ClassGrid {
    /// The latitude of the center of a row, in degrees.
    pub fn latitude(&self, row: usize) -> f64 {
        90.0 - (row as f64 + 0.5) * 180.0 / self.height as f64
    }
}

/// The parameters of the downscaling.
#[derive(Clone, Copy, Debug)]
pub struct Downscale {
    /// The fall of the temperature with height, in °C per km.
    pub lapse_rate: f32,
    /// The cells of the heightmap along the edge of one output cell.
    pub step: usize,
    pub rules: Rules,
}

/// Classifies each land cell of the heightmap. A cell is land where its
/// height is above 0 m. The temperatures of the export change by the lapse
/// rate times the difference between the heightmap and the elevation of the
/// export. Only the land cells of the export give values. A place with no
/// land cell of the export around it takes the nearest land cell.
pub fn classify(climate: &Climate, heights: &Heightmap, params: Downscale) -> ClassGrid {
    let step = params.step.max(1);
    let width = heights.width / step;
    let height = heights.height / step;
    let land = &climate.land;
    let nearest = nearest_land(land);
    let lapse = params.lapse_rate / 1000.0;
    let mut classes = vec![NO_CLASS; width * height];
    classes
        .par_chunks_mut(width)
        .enumerate()
        .for_each(|(row, out)| {
            let lat = 90.0 - (row as f64 + 0.5) * 180.0 / height as f64;
            for (col, class) in out.iter_mut().enumerate() {
                let h = heights.block_mean(col * step, row * step, step);
                if h <= 0.0 {
                    continue;
                }
                let lon = -180.0 + (col as f64 + 0.5) * 360.0 / width as f64;
                let stencil = land
                    .stencil_where(lon, lat, |i| land.values[i] > 0.0)
                    .unwrap_or_else(|| {
                        let i = nearest[land.cell(lon, lat)];
                        [(i, 1.0), (i, 0.0), (i, 0.0), (i, 0.0)]
                    });
                let dh = h - climate.elevation.apply(&stencil);
                let t =
                    std::array::from_fn(|m| climate.temperature[m].apply(&stencil) - lapse * dh);
                let p = std::array::from_fn(|m| climate.precipitation[m].apply(&stencil));
                *class = koppen::classify(&t, &p, params.rules);
            }
        });
    ClassGrid {
        width,
        height,
        classes,
    }
}

/// The slope of the elevation of the export against the mean height of
/// the heightmap over each land cell of the export. A slope of 1 means that
/// the two have the same scale.
pub fn elevation_slope(climate: &Climate, heights: &Heightmap) -> f64 {
    let grid = &climate.elevation;
    let k = heights.width / grid.width;
    if k == 0 || heights.height / grid.height != k {
        return f64::NAN;
    }
    let (mut sx, mut sy, mut sxx, mut sxy, mut n) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for row in 0..grid.height {
        for col in 0..grid.width {
            let i = row * grid.width + col;
            let x = heights.block_mean(col * k, row * k, k) as f64;
            if climate.land.values[i] > 0.0 && x > 0.0 {
                let y = grid.values[i] as f64;
                sx += x;
                sy += y;
                sxx += x * x;
                sxy += x * y;
                n += 1.0;
            }
        }
    }
    (n * sxy - sx * sy) / (n * sxx - sx * sx)
}

/// The nearest land cell of each cell, in steps across the grid. The steps
/// wrap at 180°. A grid with no land gives each cell itself.
fn nearest_land(land: &Grid) -> Vec<usize> {
    let (w, h) = (land.width, land.height);
    let mut nearest = vec![usize::MAX; w * h];
    let mut queue = VecDeque::new();
    for (i, &v) in land.values.iter().enumerate() {
        if v > 0.0 {
            nearest[i] = i;
            queue.push_back(i);
        }
    }
    while let Some(i) = queue.pop_front() {
        let (col, row) = (i % w, i / w);
        let mut next = vec![row * w + (col + 1) % w, row * w + (col + w - 1) % w];
        if row > 0 {
            next.push(i - w);
        }
        if row + 1 < h {
            next.push(i + w);
        }
        for n in next {
            if nearest[n] == usize::MAX {
                nearest[n] = nearest[i];
                queue.push_back(n);
            }
        }
    }
    nearest
        .into_iter()
        .enumerate()
        .map(|(i, n)| if n == usize::MAX { i } else { n })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn uniform(width: usize, height: usize, v: f32) -> Grid {
        Grid::new(width, height, vec![v; width * height])
    }

    fn climate(t: f32, elevation: f32) -> Climate {
        Climate {
            temperature: (0..12).map(|_| uniform(4, 2, t)).collect(),
            precipitation: (0..12).map(|_| uniform(4, 2, 100.0)).collect(),
            elevation: uniform(4, 2, elevation),
            land: uniform(4, 2, 1.0),
        }
    }

    fn heightmap(levels: Vec<u16>, width: usize, height: usize) -> Heightmap {
        // One level is one meter.
        Heightmap {
            width,
            height,
            levels,
            offset: 0.0,
            scale: 1.0,
        }
    }

    #[test]
    fn mountains_are_colder() {
        let c = climate(27.0, 0.0);
        // Sea, low land, and a peak of 3,000 m.
        let map = heightmap(vec![0, 100, 3000, 100, 0, 100, 3000, 100], 4, 2);
        let params = Downscale {
            lapse_rate: 6.5,
            step: 1,
            rules: Rules::default(),
        };
        let out = classify(&c, &map, params);
        let code = |i: usize| koppen::CLASSES[out.classes[i] as usize].0;
        assert_eq!(out.classes[0], NO_CLASS);
        assert_eq!(code(1), "Af");
        // 27 - 6.5 × 3 = 7.5 °C in each month.
        assert_eq!(code(2), "ET");
    }

    #[test]
    fn heights_of_the_export_cancel() {
        let c = climate(7.5, 3000.0);
        let map = heightmap(vec![3000; 8], 4, 2);
        let params = Downscale {
            lapse_rate: 6.5,
            step: 1,
            rules: Rules::default(),
        };
        let out = classify(&c, &map, params);
        assert_eq!(koppen::CLASSES[out.classes[0] as usize].0, "ET");
    }

    #[test]
    fn sea_of_the_export_gives_no_values() {
        let mut c = climate(27.0, 0.0);
        // Column 1 of the export is land at 0 °C, and the rest is warm sea.
        for col in 0..4 {
            for row in 0..2 {
                let i = row * 4 + col;
                c.land.values[i] = if col == 1 { 1.0 } else { 0.0 };
                for m in 0..12 {
                    c.temperature[m].values[i] = if col == 1 { -5.0 } else { 27.0 };
                }
            }
        }
        let map = heightmap(vec![10; 8], 4, 2);
        let params = Downscale {
            lapse_rate: 6.5,
            step: 1,
            rules: Rules::default(),
        };
        let out = classify(&c, &map, params);
        assert!(
            out.classes
                .iter()
                .all(|&k| koppen::CLASSES[k as usize].0 == "EF")
        );
    }

    #[test]
    fn nearest_land_wraps() {
        let mut land = uniform(4, 1, 0.0);
        land.values[0] = 1.0;
        assert_eq!(nearest_land(&land), [0, 0, 0, 0]);
    }
}
