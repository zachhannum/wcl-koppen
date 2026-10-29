//! The polygons of the zones of a class grid.
//!
//! The tracer follows the cell edges between two classes, or between a class
//! and the sea. Each run of edges between two nodes is an arc. A node is a
//! corner where three labels meet, a corner where two labels touch across a
//! diagonal, or a corner on the edge of the grid. Both zones of a shared arc
//! take the same smoothed and simplified points, so the zones have no gap
//! and no overlap.
//!
//! The grid does not wrap at 180°. A zone across the 180° meridian becomes
//! two polygons, with edges along the meridian, as RFC 7946 asks. A zone at a
//! pole has an edge along the latitude of the pole.

use std::collections::HashMap;

use crate::NO_CLASS;
use crate::downscale::ClassGrid;
use crate::koppen::Class;

/// The label of a place outside the grid.
const OUTSIDE: u16 = 0x100;

/// The parameters of the smoothing and the simplification.
#[derive(Clone, Copy, Debug)]
pub struct Shape {
    /// The width of the Gaussian smoothing, as a standard deviation in km.
    pub smooth_km: f64,
    /// The largest distance of a removed point from the simpler line, in km.
    pub tolerance_km: f64,
    pub radius_km: f64,
}

/// A polygon of one zone, in longitude and latitude. The exterior turns
/// counterclockwise and each hole turns clockwise. Each ring is closed.
#[derive(Clone, Debug, PartialEq)]
pub struct Polygon {
    pub class: Class,
    pub exterior: Vec<[f64; 2]>,
    pub holes: Vec<Vec<[f64; 2]>>,
}

type Corner = (u32, u32);

/// One closed boundary of a zone, with the zone on the left.
struct Ring {
    class: Class,
    corners: Vec<Corner>,
    /// The label at the right of each edge, from corner i to corner i + 1.
    others: Vec<u16>,
    /// The cell at the left of the first edge, as (column, row).
    inside: (u32, u32),
}

/// Directions: 0 east, 1 north, 2 west, 3 south. Y grows to the north.
const STEP: [(i64, i64); 4] = [(1, 0), (0, 1), (-1, 0), (0, -1)];

struct Tracer<'a> {
    grid: &'a ClassGrid,
    w: i64,
    h: i64,
}

impl Tracer<'_> {
    /// The label of the cell at column `c` and row `r`, where row 0 is at
    /// the north.
    fn label(&self, c: i64, r: i64) -> u16 {
        if c < 0 || r < 0 || c >= self.w || r >= self.h {
            OUTSIDE
        } else {
            self.grid.classes[(r * self.w + c) as usize] as u16
        }
    }

    /// The cells at the left and at the right of the edge from corner (x, y)
    /// in direction `d`, as (column, row).
    fn sides(&self, x: i64, y: i64, d: usize) -> ((i64, i64), (i64, i64)) {
        let h = self.h;
        match d {
            0 => ((x, h - 1 - y), (x, h - y)),
            1 => ((x - 1, h - 1 - y), (x, h - 1 - y)),
            2 => ((x - 1, h - y), (x - 1, h - 1 - y)),
            _ => ((x, h - y), (x - 1, h - y)),
        }
    }

    /// True if the edge from (x, y) in direction `d` has `class` at its left
    /// and another label at its right.
    fn is_boundary(&self, x: i64, y: i64, d: usize, class: u16) -> bool {
        let ((lc, lr), (rc, rr)) = self.sides(x, y, d);
        self.label(lc, lr) == class && self.label(rc, rr) != class
    }

    /// True if the corner (x, y) ends arcs.
    fn is_node(&self, x: i64, y: i64) -> bool {
        if x == 0 || y == 0 || x == self.w || y == self.h {
            return true;
        }
        let r = self.h - y;
        let nw = self.label(x - 1, r - 1);
        let ne = self.label(x, r - 1);
        let sw = self.label(x - 1, r);
        let se = self.label(x, r);
        let mut labels = [nw, ne, sw, se];
        labels.sort_unstable();
        let distinct = 1 + labels.windows(2).filter(|p| p[0] != p[1]).count();
        distinct >= 3 || (nw == se && ne == sw && nw != ne)
    }
}

/// Traces each zone of the grid. A cell of [`NO_CLASS`] is sea.
pub fn polygons(grid: &ClassGrid, shape: Shape) -> Vec<Polygon> {
    let tracer = Tracer {
        grid,
        w: grid.width as i64,
        h: grid.height as i64,
    };
    let rings = trace_rings(&tracer);
    let mut arcs = ArcCache::default();
    let shaped: Vec<Vec<[f64; 2]>> = rings
        .iter()
        .map(|ring| shape_ring(&tracer, ring, shape, &mut arcs))
        .collect();
    assemble(rings, shaped, grid.height as u32)
}

fn trace_rings(t: &Tracer) -> Vec<Ring> {
    let (w, h) = (t.w, t.h);
    // Bit d of a cell marks the edge in direction d with the cell at its left.
    let mut seen = vec![0u8; (w * h) as usize];
    let mut rings = Vec::new();
    for r in 0..h {
        for c in 0..w {
            let class = t.label(c, r);
            if class == NO_CLASS as u16 {
                continue;
            }
            let yb = h - 1 - r;
            // The four sides of the cell, as the edges with the cell at the left.
            let starts = [(c, yb, 0), (c + 1, yb, 1), (c + 1, yb + 1, 2), (c, yb + 1, 3)];
            for (x, y, d) in starts {
                if seen[(r * w + c) as usize] & (1 << d) != 0 || !t.is_boundary(x, y, d, class) {
                    continue;
                }
                rings.push(trace_ring(t, &mut seen, x, y, d, class));
            }
        }
    }
    rings
}

fn trace_ring(t: &Tracer, seen: &mut [u8], x0: i64, y0: i64, d0: usize, class: u16) -> Ring {
    let mut corners = Vec::new();
    let mut others = Vec::new();
    let ((ic, ir), _) = t.sides(x0, y0, d0);
    let (mut x, mut y, mut d) = (x0, y0, d0);
    loop {
        let ((lc, lr), (rc, rr)) = t.sides(x, y, d);
        seen[(lr * t.w + lc) as usize] |= 1 << d;
        corners.push((x as u32, y as u32));
        others.push(t.label(rc, rr));
        x += STEP[d].0;
        y += STEP[d].1;
        // A left turn first keeps two zones that touch at a corner apart.
        d = [(d + 1) % 4, d, (d + 3) % 4]
            .into_iter()
            .find(|&n| t.is_boundary(x, y, n, class))
            .expect("a boundary goes on at each corner");
        if (x, y, d) == (x0, y0, d0) {
            break;
        }
    }
    Ring {
        class: class as Class,
        corners,
        others,
        inside: (ic as u32, ir as u32),
    }
}

/// The shaped arcs between two zones, until the second zone takes them.
#[derive(Default)]
struct ArcCache {
    arcs: HashMap<(Corner, Corner, Corner), Vec<[f64; 2]>>,
}

fn shape_ring(t: &Tracer, ring: &Ring, shape: Shape, cache: &mut ArcCache) -> Vec<[f64; 2]> {
    let n = ring.corners.len();
    let nodes: Vec<usize> = (0..n)
        .filter(|&i| {
            let (x, y) = ring.corners[i];
            t.is_node(x as i64, y as i64)
        })
        .collect();
    let mut out: Vec<[f64; 2]> = Vec::new();
    if nodes.is_empty() {
        let points = shape_loop(t, &ring.corners, ring.others[0], shape, cache);
        out.extend(points);
        out.push(out[0]);
        return out;
    }
    for (k, &start) in nodes.iter().enumerate() {
        let end = nodes[(k + 1) % nodes.len()];
        let len = (end + n - start) % n;
        let len = if len == 0 { n } else { len };
        let corners: Vec<Corner> = (0..=len).map(|j| ring.corners[(start + j) % n]).collect();
        let points = shape_arc(t, &corners, ring.others[start], shape, cache);
        let skip = usize::from(!out.is_empty());
        out.extend(points.into_iter().skip(skip));
    }
    out
}

/// The points of an arc between two nodes, from its first corner to its last.
fn shape_arc(
    t: &Tracer,
    corners: &[Corner],
    other: u16,
    shape: Shape,
    cache: &mut ArcCache,
) -> Vec<[f64; 2]> {
    let n = corners.len();
    if other == OUTSIDE {
        // An arc along the 180° meridian or a pole is straight.
        return vec![lonlat(t, corners[0]), lonlat(t, corners[n - 1])];
    }
    if other == NO_CLASS as u16 {
        return shape_open(t, corners, shape);
    }
    let forward = (corners[0], corners[1], corners[n - 1]);
    let backward = (corners[n - 1], corners[n - 2], corners[0]);
    let reversed = backward < forward;
    let key = if reversed { backward } else { forward };
    let points = match cache.arcs.remove(&key) {
        Some(points) => points,
        None => {
            let mut canonical = corners.to_vec();
            if reversed {
                canonical.reverse();
            }
            let points = shape_open(t, &canonical, shape);
            cache.arcs.insert(key, points.clone());
            points
        }
    };
    if reversed {
        points.into_iter().rev().collect()
    } else {
        points
    }
}

/// The points of a closed boundary with no node, in the order of `corners`,
/// without the closing point.
fn shape_loop(
    t: &Tracer,
    corners: &[Corner],
    other: u16,
    shape: Shape,
    cache: &mut ArcCache,
) -> Vec<[f64; 2]> {
    let n = corners.len();
    // The canonical loop starts at the least corner and goes to the lesser
    // of its two neighbors, so both zones of the loop get the same points.
    let first = (0..n).min_by_key(|&i| corners[i]).unwrap();
    let next = corners[(first + 1) % n];
    let prev = corners[(first + n - 1) % n];
    let reversed = prev < next;
    let canonical: Vec<Corner> = if reversed {
        (0..n).map(|j| corners[(first + n - j) % n]).collect()
    } else {
        (0..n).map(|j| corners[(first + j) % n]).collect()
    };
    let key = (canonical[0], canonical[1], canonical[0]);
    let shared = other != NO_CLASS as u16;
    let points = match shared.then(|| cache.arcs.remove(&key)).flatten() {
        Some(points) => points,
        None => {
            let points = shape_closed(t, &canonical, shape);
            if shared {
                cache.arcs.insert(key, points.clone());
            }
            points
        }
    };
    // The ring may start at any point, so only the direction matters.
    if reversed {
        points.into_iter().rev().collect()
    } else {
        points
    }
}

fn lonlat(t: &Tracer, (x, y): Corner) -> [f64; 2] {
    [
        -180.0 + x as f64 * 360.0 / t.w as f64,
        -90.0 + y as f64 * 180.0 / t.h as f64,
    ]
}

fn to_vector([lon, lat]: [f64; 2]) -> [f64; 3] {
    let (lon, lat) = (lon.to_radians(), lat.to_radians());
    [lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin()]
}

/// The longitude and latitude of a vector, with the longitude near `near`.
fn from_vector(v: [f64; 3], near: [f64; 2]) -> [f64; 2] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let lat = (v[2] / len).clamp(-1.0, 1.0).asin().to_degrees();
    let mut lon = v[1].atan2(v[0]).to_degrees();
    if lon - near[0] > 180.0 {
        lon -= 360.0;
    } else if near[0] - lon > 180.0 {
        lon += 360.0;
    }
    [lon.clamp(-180.0, 180.0), lat.clamp(-90.0, 90.0)]
}

fn distance_km(a: [f64; 3], b: [f64; 3], radius_km: f64) -> f64 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() * radius_km
}

/// Smooths and simplifies an open arc. The ends stay where they are.
fn shape_open(t: &Tracer, corners: &[Corner], shape: Shape) -> Vec<[f64; 2]> {
    let raw: Vec<[f64; 2]> = corners.iter().map(|&c| lonlat(t, c)).collect();
    let n = raw.len();
    let v: Vec<[f64; 3]> = raw.iter().map(|&p| to_vector(p)).collect();
    let mut s = vec![0.0; n];
    for i in 1..n {
        s[i] = s[i - 1] + distance_km(v[i - 1], v[i], shape.radius_km);
    }
    let sigma = shape.smooth_km;
    let reach = 3.0 * sigma;
    let mut smooth = raw.clone();
    if sigma > 0.0 {
        for i in 1..n - 1 {
            let mut sum = [0.0; 3];
            let mut add = |p: [f64; 3], ds: f64| {
                let w = (-(ds * ds) / (2.0 * sigma * sigma)).exp();
                for k in 0..3 {
                    sum[k] += w * p[k];
                }
            };
            add(v[i], 0.0);
            let mut j = i;
            while j > 0 && s[i] - s[j - 1] <= reach {
                j -= 1;
                add(v[j], s[i] - s[j]);
            }
            // Beyond an end, the arc goes on as its reflection through the end.
            if j == 0 {
                let mut k = 1;
                while k < n && s[i] + s[k] <= reach {
                    add(reflect(v[0], v[k]), s[i] + s[k]);
                    k += 1;
                }
            }
            let mut j = i;
            while j + 1 < n && s[j + 1] - s[i] <= reach {
                j += 1;
                add(v[j], s[j] - s[i]);
            }
            if j == n - 1 {
                let end = s[n - 1];
                let mut k = n - 2;
                while (end - s[i]) + (end - s[k]) <= reach {
                    add(reflect(v[n - 1], v[k]), (end - s[i]) + (end - s[k]));
                    if k == 0 {
                        break;
                    }
                    k -= 1;
                }
            }
            smooth[i] = from_vector(sum, raw[i]);
        }
    }
    simplify(&smooth, shape)
}

/// The reflection of `p` through `center`.
fn reflect(center: [f64; 3], p: [f64; 3]) -> [f64; 3] {
    [
        2.0 * center[0] - p[0],
        2.0 * center[1] - p[1],
        2.0 * center[2] - p[2],
    ]
}

/// Smooths and simplifies a closed loop. Returns the points without the
/// closing point.
fn shape_closed(t: &Tracer, corners: &[Corner], shape: Shape) -> Vec<[f64; 2]> {
    let raw: Vec<[f64; 2]> = corners.iter().map(|&c| lonlat(t, c)).collect();
    let n = raw.len();
    let v: Vec<[f64; 3]> = raw.iter().map(|&p| to_vector(p)).collect();
    let step: Vec<f64> = (0..n)
        .map(|i| distance_km(v[i], v[(i + 1) % n], shape.radius_km))
        .collect();
    let total: f64 = step.iter().sum();
    let sigma = shape.smooth_km;
    let reach = (3.0 * sigma).min(total / 2.0);
    let mut smooth = raw.clone();
    if sigma > 0.0 {
        for i in 0..n {
            let mut sum = [0.0; 3];
            for k in 0..3 {
                sum[k] += v[i][k];
            }
            let mut add = |j: usize, ds: f64| {
                let w = (-(ds * ds) / (2.0 * sigma * sigma)).exp();
                for k in 0..3 {
                    sum[k] += w * v[j][k];
                }
            };
            let (mut ds, mut j) = (0.0, i);
            for _ in 0..n / 2 {
                ds += step[j];
                j = (j + 1) % n;
                if ds > reach {
                    break;
                }
                add(j, ds);
            }
            let (mut ds, mut j) = (0.0, i);
            for _ in 0..(n - 1) / 2 {
                j = (j + n - 1) % n;
                ds += step[j];
                if ds > reach {
                    break;
                }
                add(j, ds);
            }
            smooth[i] = from_vector(sum, raw[i]);
        }
    }
    // The loop splits at its first point and at the point farthest from it.
    let far = (0..n)
        .max_by(|&a, &b| {
            let d = |i: usize| planar(smooth[0], smooth[i], shape.radius_km);
            let (da, db) = (d(a), d(b));
            (da[0] * da[0] + da[1] * da[1]).total_cmp(&(db[0] * db[0] + db[1] * db[1]))
        })
        .unwrap();
    let mut first: Vec<[f64; 2]> = smooth[..=far].to_vec();
    let mut second: Vec<[f64; 2]> = smooth[far..].to_vec();
    second.push(smooth[0]);
    first = simplify(&first, shape);
    second = simplify(&second, shape);
    let mut out = first;
    out.extend(second.into_iter().skip(1));
    out.pop();
    if out.len() < 4 && n >= 4 {
        // A tiny loop keeps a shape with an area.
        return (0..4).map(|k| smooth[k * n / 4]).collect();
    }
    out
}

/// The offset of `p` from `origin` in km, on a plane at the latitude of
/// `origin`.
fn planar(origin: [f64; 2], p: [f64; 2], radius_km: f64) -> [f64; 2] {
    let k = radius_km * std::f64::consts::PI / 180.0;
    let cos = origin[1].to_radians().cos();
    [(p[0] - origin[0]) * k * cos, (p[1] - origin[1]) * k]
}

/// The Douglas-Peucker simplification of an open line. The ends stay.
fn simplify(points: &[[f64; 2]], shape: Shape) -> Vec<[f64; 2]> {
    let n = points.len();
    if n <= 2 {
        return points.to_vec();
    }
    let mut keep = vec![false; n];
    keep[0] = true;
    keep[n - 1] = true;
    let mut stack = vec![(0, n - 1)];
    while let Some((a, b)) = stack.pop() {
        if b <= a + 1 {
            continue;
        }
        let end = planar(points[a], points[b], shape.radius_km);
        let len2 = end[0] * end[0] + end[1] * end[1];
        let (mut worst, mut at) = (0.0, a);
        for (i, &p) in points.iter().enumerate().take(b).skip(a + 1) {
            let q = planar(points[a], p, shape.radius_km);
            let d2 = if len2 == 0.0 {
                q[0] * q[0] + q[1] * q[1]
            } else {
                let cross = q[0] * end[1] - q[1] * end[0];
                cross * cross / len2
            };
            if d2 > worst {
                worst = d2;
                at = i;
            }
        }
        if worst > shape.tolerance_km * shape.tolerance_km {
            keep[at] = true;
            stack.push((a, at));
            stack.push((at, b));
        }
    }
    points
        .iter()
        .zip(keep)
        .filter_map(|(&p, k)| k.then_some(p))
        .collect()
}

/// The signed area of a ring in longitude and latitude units.
fn signed_area(points: &[[f64; 2]]) -> f64 {
    points
        .windows(2)
        .map(|p| p[0][0] * p[1][1] - p[1][0] * p[0][1])
        .sum::<f64>()
        / 2.0
}

/// True if a point is inside the corners of a ring, by the even-odd rule.
fn contains(corners: &[Corner], x: f64, y: f64) -> bool {
    let n = corners.len();
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = (corners[i].0 as f64, corners[i].1 as f64);
        let (xj, yj) = (corners[j].0 as f64, corners[j].1 as f64);
        if (yi > y) != (yj > y) && x < (xj - xi) * (y - yi) / (yj - yi) + xi {
            inside = !inside;
        }
        j = i;
    }
    inside
}

/// The corners of a ring where it turns.
fn turns(corners: &[Corner]) -> Vec<Corner> {
    let n = corners.len();
    (0..n)
        .filter(|&i| {
            let (a, b, c) = (corners[(i + n - 1) % n], corners[i], corners[(i + 1) % n]);
            (a.0 == b.0) != (b.0 == c.0)
        })
        .collect::<Vec<_>>()
        .into_iter()
        .map(|i| corners[i])
        .collect()
}

/// Joins each hole to the smallest exterior of its class around it.
fn assemble(rings: Vec<Ring>, shaped: Vec<Vec<[f64; 2]>>, height: u32) -> Vec<Polygon> {
    struct Outer {
        polygon: Polygon,
        turns: Vec<Corner>,
        bbox: (u32, u32, u32, u32),
        area: f64,
    }
    let mut outers: Vec<Outer> = Vec::new();
    let mut holes: Vec<(Class, (u32, u32), Vec<[f64; 2]>)> = Vec::new();
    for (ring, points) in rings.into_iter().zip(shaped) {
        let area = signed_area(&points);
        if area > 0.0 {
            let turns = turns(&ring.corners);
            let bbox = turns.iter().fold((u32::MAX, u32::MAX, 0, 0), |b, &(x, y)| {
                (b.0.min(x), b.1.min(y), b.2.max(x), b.3.max(y))
            });
            outers.push(Outer {
                polygon: Polygon {
                    class: ring.class,
                    exterior: points,
                    holes: Vec::new(),
                },
                turns,
                bbox,
                area,
            });
        } else {
            holes.push((ring.class, ring.inside, points));
        }
    }
    let mut order: Vec<usize> = (0..outers.len()).collect();
    order.sort_by(|&a, &b| outers[a].area.total_cmp(&outers[b].area));
    for (class, (col, row), points) in holes {
        // The center of a cell of the zone, in corner units with Y to the north.
        let x = col as f64 + 0.5;
        let y = height as f64 - row as f64 - 0.5;
        let owner = order.iter().copied().find(|&k| {
            let o = &outers[k];
            if o.polygon.class != class {
                return false;
            }
            let (x0, y0, x1, y1) = o.bbox;
            x >= x0 as f64 && x <= x1 as f64 && y >= y0 as f64 && y <= y1 as f64 && contains(&o.turns, x, y)
        });
        if let Some(k) = owner {
            outers[k].polygon.holes.push(points);
        }
    }
    outers.into_iter().map(|o| o.polygon).collect()
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

    const SHARP: Shape = Shape {
        smooth_km: 0.0,
        tolerance_km: 0.0,
        radius_km: 6371.0,
    };

    const SMOOTH: Shape = Shape {
        smooth_km: 300.0,
        tolerance_km: 1.0,
        radius_km: 6371.0,
    };

    fn sorted(points: &[[f64; 2]]) -> Vec<(i64, i64)> {
        let mut v: Vec<(i64, i64)> = points
            .iter()
            .map(|p| ((p[0] * 1e9).round() as i64, (p[1] * 1e9).round() as i64))
            .collect();
        v.sort_unstable();
        v.dedup();
        v
    }

    #[test]
    fn square_in_the_sea() {
        #[rustfmt::skip]
        let g = grid(8, vec![
            S, S, S, S, S, S, S, S,
            S, S, 1, 1, S, S, S, S,
            S, S, 1, 1, S, S, S, S,
            S, S, S, S, S, S, S, S,
        ]);
        let p = polygons(&g, SHARP);
        assert_eq!(p.len(), 1);
        assert_eq!(p[0].class, 1);
        assert!(p[0].holes.is_empty());
        let ring = &p[0].exterior;
        assert_eq!(ring.first(), ring.last());
        assert!(signed_area(ring) > 0.0, "the exterior turns counterclockwise");
        assert_eq!(
            sorted(ring),
            sorted(&[[-90.0, 45.0], [0.0, 45.0], [0.0, -45.0], [-90.0, -45.0]])
        );
    }

    #[test]
    fn shared_boundary_has_the_same_points() {
        #[rustfmt::skip]
        let g = grid(8, vec![
            S, S, S, S, S, S, S, S,
            S, 1, 1, 1, 2, 2, S, S,
            S, 1, 1, 2, 2, 2, S, S,
            S, 1, 2, 2, 2, 2, S, S,
            S, S, S, S, S, S, S, S,
            S, S, S, S, S, S, S, S,
        ]);
        let p = polygons(&g, SMOOTH);
        assert_eq!(p.len(), 2);
        let a: std::collections::BTreeSet<_> = sorted(&p[0].exterior).into_iter().collect();
        let b: std::collections::BTreeSet<_> = sorted(&p[1].exterior).into_iter().collect();
        // The shared arc has more than its two nodes.
        assert!(a.intersection(&b).count() > 2);
    }

    #[test]
    fn island_is_a_hole() {
        #[rustfmt::skip]
        let g = grid(8, vec![
            1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 2, 2, 1, 1, 1,
            1, 1, 1, 2, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1,
        ]);
        let p = polygons(&g, SMOOTH);
        let outer = p.iter().find(|p| p.class == 1).unwrap();
        let island = p.iter().find(|p| p.class == 2).unwrap();
        assert_eq!(outer.holes.len(), 1);
        assert!(signed_area(&outer.holes[0]) < 0.0, "a hole turns clockwise");
        assert_eq!(sorted(&outer.holes[0]), sorted(&island.exterior));
    }

    #[test]
    fn hole_goes_to_the_zone_around_it() {
        // Zone 1 holds a lake, and an island of zone 1 is in the lake.
        #[rustfmt::skip]
        let g = grid(8, vec![
            1, 1, 1, 1, 1, 1, 1, 1,
            1, S, S, S, S, S, 1, 1,
            1, S, 1, 1, S, S, 1, 1,
            1, S, S, S, S, S, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1,
        ]);
        let p = polygons(&g, SHARP);
        assert_eq!(p.len(), 2);
        let big = p.iter().find(|p| p.holes.len() == 1).unwrap();
        let small = p.iter().find(|p| p.holes.is_empty()).unwrap();
        assert!(signed_area(&big.exterior) > signed_area(&small.exterior));
    }

    #[test]
    fn zone_across_the_seam_splits() {
        #[rustfmt::skip]
        let g = grid(8, vec![
            S, S, S, S, S, S, S, S,
            1, S, S, S, S, S, S, 1,
            1, S, S, S, S, S, S, 1,
            S, S, S, S, S, S, S, S,
        ]);
        let p = polygons(&g, SMOOTH);
        assert_eq!(p.len(), 2);
        for poly in &p {
            let lons: Vec<f64> = poly.exterior.iter().map(|q| q[0]).collect();
            let west = lons.iter().all(|&l| l <= -90.0);
            let east = lons.iter().all(|&l| l >= 90.0);
            assert!(west || east, "{lons:?}");
            assert!(lons.iter().any(|&l| l.abs() == 180.0));
        }
    }

    #[test]
    fn corner_touch_gives_two_zones() {
        #[rustfmt::skip]
        let g = grid(8, vec![
            S, S, S, S, S, S, S, S,
            S, 1, 2, S, S, S, S, S,
            S, 2, 1, S, S, S, S, S,
            S, S, S, S, S, S, S, S,
        ]);
        let p = polygons(&g, SHARP);
        assert_eq!(p.iter().filter(|p| p.class == 1).count(), 2);
        assert_eq!(p.iter().filter(|p| p.class == 2).count(), 2);
    }

    #[test]
    fn zone_at_the_pole() {
        #[rustfmt::skip]
        let g = grid(8, vec![
            1, 1, 1, 1, 1, 1, 1, 1,
            S, S, S, S, S, S, S, S,
            S, S, S, S, S, S, S, S,
            S, S, S, S, S, S, S, S,
        ]);
        let p = polygons(&g, SMOOTH);
        assert_eq!(p.len(), 1);
        let ring = &p[0].exterior;
        assert!(signed_area(ring) > 0.0);
        assert!(ring.contains(&[-180.0, 90.0]) && ring.contains(&[180.0, 90.0]));
        // No point of the ring is south of the zone.
        assert!(ring.iter().all(|q| q[1] >= 45.0 - 1e-9));
    }
}
