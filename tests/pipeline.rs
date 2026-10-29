//! Runs the whole chain on a small synthetic world.

use wcl_koppen::downscale::{self, Downscale};
use wcl_koppen::grid::Grid;
use wcl_koppen::height::Heightmap;
use wcl_koppen::koppen::{CLASSES, Rules};
use wcl_koppen::trace::{self, Shape};
use wcl_koppen::wcl::Climate;
use wcl_koppen::{NO_CLASS, filter, geojson};

const RADIUS: f64 = 6371.0;

/// A climate that is warm at the equator and cold at the poles, with the
/// same rain everywhere, on land everywhere.
fn climate() -> Climate {
    let (w, h) = (64, 32);
    let temperature = |row: usize| 30.0 - 50.0 * ((row as f32 + 0.5) / h as f32 - 0.5).abs() * 2.0;
    let t: Vec<f32> = (0..w * h).map(|i| temperature(i / w)).collect();
    Climate {
        temperature: (0..12).map(|_| Grid::new(w, h, t.clone())).collect(),
        precipitation: (0..12)
            .map(|_| Grid::new(w, h, vec![80.0; w * h]))
            .collect(),
        elevation: Grid::new(w, h, vec![0.0; w * h]),
        land: Grid::new(w, h, vec![1.0; w * h]),
    }
}

/// A continent with a high range, in a sea.
fn heightmap() -> Heightmap {
    let (w, h) = (512, 256);
    let mut levels = vec![0u16; w * h];
    for row in 40..216 {
        for col in 100..400 {
            let range = (col as i64 - 250).abs() < 8;
            levels[row * w + col] = if range { 4000 } else { 200 };
        }
    }
    Heightmap {
        width: w,
        height: h,
        levels,
        offset: 0.0,
        scale: 1.0,
    }
}

#[test]
fn zones_cover_the_land() {
    let heights = heightmap();
    let params = Downscale {
        lapse_rate: 6.5,
        step: 1,
        rules: Rules::default(),
    };
    let mut grid = downscale::classify(&climate(), &heights, params);
    let cell_km = std::f64::consts::TAU * RADIUS / grid.width as f64;
    filter::merge_small(&mut grid, 4.0 * cell_km * cell_km, RADIUS);

    // The range is colder than the land around it.
    let at = |col: usize, row: usize| CLASSES[grid.classes[row * grid.width + col] as usize].0;
    assert_eq!(at(150, 128), "Af");
    assert_ne!(at(250, 128), "Af");

    let areas = filter::row_areas(grid.height, grid.width, RADIUS);
    let land: f64 = (0..grid.classes.len())
        .filter(|&i| grid.classes[i] != NO_CLASS)
        .map(|i| areas[i / grid.width])
        .sum();
    let shape = Shape {
        smooth_km: cell_km,
        tolerance_km: cell_km / 4.0,
        radius_km: RADIUS,
    };
    let polygons = trace::polygons(&grid, shape);
    let total: f64 = polygons.iter().map(|p| geojson::area_km2(p, RADIUS)).sum();
    assert!((total - land).abs() / land < 0.01, "{total} of {land} km²");
    assert!(polygons.iter().any(|p| CLASSES[p.class as usize].0 == "Af"));
}
