//! Köppen-Geiger climate zones from a World Climate Lab export.

pub mod downscale;
pub mod grid;
pub mod height;
pub mod koppen;
pub mod wcl;

use grid::Grid;
use koppen::{Class, Rules};
use wcl::Climate;

/// The class of a cell with no class.
pub const NO_CLASS: Class = Class::MAX;

/// Classifies each land cell of the export on its own grid.
pub fn classify_grid(climate: &Climate, rules: Rules) -> Vec<Class> {
    let land = &climate.land;
    let mut out = vec![NO_CLASS; land.width * land.height];
    for row in 0..land.height {
        for col in 0..land.width {
            let i = row * land.width + col;
            if land.values[i] == 0.0 {
                continue;
            }
            let t = months(&climate.temperature, i);
            let p = months(&climate.precipitation, i);
            out[i] = koppen::classify(&t, &p, rules);
        }
    }
    out
}

fn months(grids: &[Grid], i: usize) -> [f32; 12] {
    std::array::from_fn(|m| grids[m].values[i])
}
