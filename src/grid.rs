//! Global equirectangular grids.

/// Four cells and their weights.
pub type Stencil = [(usize, f32); 4];

/// A grid of values over the whole globe. Row 0 is at the north, and column 0
/// is at the west edge, at -180°. Each value is the mean over its cell.
#[derive(Clone, Debug, PartialEq)]
pub struct Grid {
    pub width: usize,
    pub height: usize,
    pub values: Vec<f32>,
}

impl Grid {
    pub fn new(width: usize, height: usize, values: Vec<f32>) -> Self {
        assert_eq!(
            values.len(),
            width * height,
            "the grid has the wrong length"
        );
        Self {
            width,
            height,
            values,
        }
    }

    pub fn get(&self, col: usize, row: usize) -> f32 {
        self.values[row * self.width + col]
    }

    /// The longitude and the latitude of the center of a cell, in degrees.
    pub fn center(&self, col: usize, row: usize) -> (f64, f64) {
        let lon = -180.0 + (col as f64 + 0.5) * 360.0 / self.width as f64;
        let lat = 90.0 - (row as f64 + 0.5) * 180.0 / self.height as f64;
        (lon, lat)
    }

    /// The bilinear blend of the four cell centers around a place. The blend
    /// wraps at 180° and clamps at the rows next to the poles.
    pub fn sample(&self, lon: f64, lat: f64) -> f32 {
        let (x, y) = self.position(lon, lat);
        let (c0, c1, fx) = self.columns(x);
        let (r0, r1, fy) = self.rows(y);
        let top = lerp(self.get(c0, r0), self.get(c1, r0), fx);
        let bottom = lerp(self.get(c0, r1), self.get(c1, r1), fx);
        lerp(top, bottom, fy)
    }

    /// The bilinear blend of the four cells around a place, with weight only
    /// on the cells where `keep` is true. `None`: none of the four cells is
    /// kept.
    pub fn sample_where(&self, lon: f64, lat: f64, keep: impl Fn(usize) -> bool) -> Option<f32> {
        self.stencil_where(lon, lat, keep).map(|s| self.apply(&s))
    }

    /// The cells and the weights of [`Grid::sample_where`]. The weights add
    /// up to 1. A grid of the same size takes the same stencil.
    pub fn stencil_where(
        &self,
        lon: f64,
        lat: f64,
        keep: impl Fn(usize) -> bool,
    ) -> Option<Stencil> {
        let (x, y) = self.position(lon, lat);
        let (c0, c1, fx) = self.columns(x);
        let (r0, r1, fy) = self.rows(y);
        let mut cells = [
            (r0 * self.width + c0, (1.0 - fx) * (1.0 - fy)),
            (r0 * self.width + c1, fx * (1.0 - fy)),
            (r1 * self.width + c0, (1.0 - fx) * fy),
            (r1 * self.width + c1, fx * fy),
        ];
        let mut weight = 0.0;
        for (i, w) in &mut cells {
            if keep(*i) {
                weight += *w;
            } else {
                *w = 0.0;
            }
        }
        if weight > 0.0 {
            for (_, w) in &mut cells {
                *w /= weight;
            }
            return Some(cells);
        }
        // A kept cell with a weight of zero still gives the place a value.
        let i = cells.iter().map(|&(i, _)| i).find(|&i| keep(i))?;
        Some([(i, 1.0), (i, 0.0), (i, 0.0), (i, 0.0)])
    }

    /// The blend of the cells of a stencil.
    pub fn apply(&self, stencil: &Stencil) -> f32 {
        stencil.iter().map(|&(i, w)| self.values[i] * w).sum()
    }

    /// The cell that holds a place.
    pub fn cell(&self, lon: f64, lat: f64) -> usize {
        let col = ((lon + 180.0) / 360.0 * self.width as f64).floor() as i64;
        let row = ((90.0 - lat) / 180.0 * self.height as f64).floor() as i64;
        let col = col.rem_euclid(self.width as i64) as usize;
        let row = row.clamp(0, self.height as i64 - 1) as usize;
        row * self.width + col
    }

    /// The place in cell units, where the center of cell (0, 0) is at (0, 0).
    fn position(&self, lon: f64, lat: f64) -> (f64, f64) {
        let x = (lon + 180.0) * self.width as f64 / 360.0 - 0.5;
        let y = (90.0 - lat) * self.height as f64 / 180.0 - 0.5;
        (x, y)
    }

    fn columns(&self, x: f64) -> (usize, usize, f32) {
        let x0 = x.floor();
        let w = self.width as i64;
        let c0 = (x0 as i64).rem_euclid(w) as usize;
        let c1 = (x0 as i64 + 1).rem_euclid(w) as usize;
        (c0, c1, (x - x0) as f32)
    }

    fn rows(&self, y: f64) -> (usize, usize, f32) {
        let last = self.height as f64 - 1.0;
        let y = y.clamp(0.0, last);
        let y0 = y.floor().min(last - 1.0).max(0.0);
        let r0 = y0 as usize;
        let r1 = (r0 + 1).min(self.height - 1);
        (r0, r1, (y - y0) as f32)
    }
}

fn lerp(a: f32, b: f32, t: f32) -> f32 {
    a + (b - a) * t
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ramp() -> Grid {
        // The value is the column number.
        let values = (0..4 * 2).map(|i| (i % 4) as f32).collect();
        Grid::new(4, 2, values)
    }

    #[test]
    fn center_of_a_cell_gives_its_value() {
        let g = ramp();
        for col in 0..4 {
            let (lon, lat) = g.center(col, 1);
            assert_eq!(g.sample(lon, lat), col as f32);
        }
    }

    #[test]
    fn blend_wraps_at_the_seam() {
        let g = ramp();
        // Halfway between the last column (3) and the first column (0).
        assert_eq!(g.sample(180.0, 0.0), 1.5);
        assert_eq!(g.sample(-180.0, 0.0), 1.5);
    }

    #[test]
    fn blend_clamps_at_the_poles() {
        let g = Grid::new(1, 2, vec![10.0, 20.0]);
        assert_eq!(g.sample(0.0, 90.0), 10.0);
        assert_eq!(g.sample(0.0, -90.0), 20.0);
        assert_eq!(g.sample(0.0, 0.0), 15.0);
    }

    #[test]
    fn blend_where_skips_cells() {
        let g = ramp();
        let (lon, lat) = g.center(1, 0);
        // Between columns 1 and 2, with column 2 left out.
        let v = g.sample_where(lon + 45.0, lat, |i| i % 4 != 2);
        assert_eq!(v, Some(1.0));
        assert_eq!(g.sample_where(lon, lat, |_| false), None);
    }
}
