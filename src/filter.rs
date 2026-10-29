//! The merge of small zones into their neighbors.

use std::collections::VecDeque;

use crate::NO_CLASS;
use crate::downscale::ClassGrid;
use crate::koppen::Class;

/// A cell of a small zone while the merge runs.
const OPEN: Class = Class::MAX - 1;

/// The area of each cell of a row, in km², on a sphere of `radius_km`.
pub fn row_areas(height: usize, width: usize, radius_km: f64) -> Vec<f64> {
    let dlon = std::f64::consts::TAU / width as f64;
    (0..height)
        .map(|row| {
            let north = (90.0 - row as f64 * 180.0 / height as f64).to_radians();
            let south = (90.0 - (row + 1) as f64 * 180.0 / height as f64).to_radians();
            radius_km * radius_km * dlon * (north.sin() - south.sin())
        })
        .collect()
}

/// The four neighbors of a cell. The neighbors wrap at 180°, and a cell at a
/// pole has three.
fn neighbors(i: usize, width: usize, height: usize) -> impl Iterator<Item = usize> {
    let (col, row) = (i % width, i / width);
    let east = row * width + (col + 1) % width;
    let west = row * width + (col + width - 1) % width;
    let north = (row > 0).then(|| i - width);
    let south = (row + 1 < height).then(|| i + width);
    [Some(east), Some(west), north, south].into_iter().flatten()
}

/// Merges each zone below `min_area_km2` into the zones around it. A zone
/// is a set of cells of one class that join at their sides. Each cell of a
/// small zone takes the class of the nearest cell of a large zone on the same
/// land. A small zone with no large zone on its land keeps its class.
/// Returns the count of merged zones.
pub fn merge_small(grid: &mut ClassGrid, min_area_km2: f64, radius_km: f64) -> usize {
    let (w, h) = (grid.width, grid.height);
    let areas = row_areas(h, w, radius_km);
    let mut seen = vec![false; w * h];
    let mut small: Vec<(u32, Class)> = Vec::new();
    let mut merged = 0;
    let mut queue = VecDeque::new();
    for start in 0..w * h {
        let class = grid.classes[start];
        if seen[start] || class == NO_CLASS {
            continue;
        }
        seen[start] = true;
        queue.push_back(start);
        let mut area = 0.0;
        let mut cells = Vec::new();
        while let Some(i) = queue.pop_front() {
            area += areas[i / w];
            if area < min_area_km2 {
                cells.push(i as u32);
            } else if !cells.is_empty() {
                cells = Vec::new();
            }
            for n in neighbors(i, w, h) {
                if !seen[n] && grid.classes[n] == class {
                    seen[n] = true;
                    queue.push_back(n);
                }
            }
        }
        if area < min_area_km2 {
            merged += 1;
            small.extend(cells.into_iter().map(|i| (i, class)));
        }
    }
    drop(seen);
    for &(i, _) in &small {
        grid.classes[i as usize] = OPEN;
    }

    // The cells of the large zones next to an open cell grow into the open cells.
    let mut queue: VecDeque<usize> = VecDeque::new();
    for &(i, _) in &small {
        for n in neighbors(i as usize, w, h) {
            let c = grid.classes[n];
            if c != OPEN && c != NO_CLASS {
                queue.push_back(n);
            }
        }
    }
    while let Some(i) = queue.pop_front() {
        let class = grid.classes[i];
        for n in neighbors(i, w, h) {
            if grid.classes[n] == OPEN {
                grid.classes[n] = class;
                queue.push_back(n);
            }
        }
    }
    for (i, class) in small {
        if grid.classes[i as usize] == OPEN {
            grid.classes[i as usize] = class;
        }
    }
    merged
}

#[cfg(test)]
mod tests {
    use super::*;

    const S: Class = NO_CLASS;

    fn grid(width: usize, classes: Vec<Class>) -> ClassGrid {
        ClassGrid {
            width,
            height: classes.len() / width,
            classes,
        }
    }

    #[test]
    fn areas_add_up_to_the_sphere() {
        let total: f64 = row_areas(90, 180, 1.0).iter().sum::<f64>() * 180.0;
        assert!((total - 4.0 * std::f64::consts::PI).abs() < 1e-9);
    }

    #[test]
    fn small_zone_joins_its_neighbor() {
        let mut g = grid(4, vec![1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]);
        // The sphere has a radius of 1, so the area of the sphere is 4π.
        let n = merge_small(&mut g, 2.0, 1.0);
        assert_eq!(n, 1);
        assert!(g.classes.iter().all(|&c| c == 1));
    }

    #[test]
    fn island_keeps_its_class() {
        let mut g = grid(4, vec![S, S, S, S, S, 2, S, S, S, S, S, S, S, S, S, S]);
        merge_small(&mut g, 2.0, 1.0);
        assert_eq!(g.classes[5], 2);
    }

    #[test]
    fn zone_across_the_seam_is_one_zone() {
        // Class 2 is in the first and the last column of each row. Together
        // the cells are half of the sphere.
        let mut g = grid(4, vec![2, 1, 1, 2, 2, 1, 1, 2]);
        merge_small(&mut g, 5.0, 1.0);
        assert_eq!(g.classes, [2, 1, 1, 2, 2, 1, 1, 2]);
    }
}
