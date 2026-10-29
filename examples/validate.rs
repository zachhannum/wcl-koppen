//! Counts the crossings of the edges of the zones of a GeoJSON file: within
//! a polygon, and between the polygons. A valid file has none.
//!
//! ```sh
//! cargo run --release --example validate -- zones.geojson
//! ```

use std::collections::HashMap;

type P = [f64; 2];

fn cross(o: P, a: P, b: P) -> f64 {
    (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
}

/// True if the two segments cross at a point inside both.
fn proper(a: P, b: P, c: P, d: P) -> bool {
    let (d1, d2) = (cross(c, d, a), cross(c, d, b));
    let (d3, d4) = (cross(a, b, c), cross(a, b, d));
    d1 * d2 < 0.0 && d3 * d4 < 0.0
}

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: validate ZONES.geojson");
    let json: serde_json::Value =
        serde_json::from_reader(std::fs::File::open(path).unwrap()).unwrap();
    let mut segments: Vec<(P, P, usize)> = Vec::new();
    for (f, feature) in json["features"].as_array().unwrap().iter().enumerate() {
        for ring in feature["geometry"]["coordinates"].as_array().unwrap() {
            let pts: Vec<P> = ring
                .as_array()
                .unwrap()
                .iter()
                .map(|p| [p[0].as_f64().unwrap(), p[1].as_f64().unwrap()])
                .collect();
            for w in pts.windows(2) {
                segments.push((w[0], w[1], f));
            }
        }
    }
    let cell = 0.1;
    let mut bins: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
    for (i, &(a, b, _)) in segments.iter().enumerate() {
        let (x0, x1) = (
            (a[0].min(b[0]) / cell).floor() as i64,
            (a[0].max(b[0]) / cell).floor() as i64,
        );
        let (y0, y1) = (
            (a[1].min(b[1]) / cell).floor() as i64,
            (a[1].max(b[1]) / cell).floor() as i64,
        );
        for x in x0..=x1 {
            for y in y0..=y1 {
                bins.entry((x, y)).or_default().push(i);
            }
        }
    }
    let mut found = std::collections::HashSet::new();
    for list in bins.values() {
        for (k, &i) in list.iter().enumerate() {
            for &j in &list[k + 1..] {
                let (a, b, _) = segments[i];
                let (c, d, _) = segments[j];
                if proper(a, b, c, d) {
                    found.insert((i.min(j), i.max(j)));
                }
            }
        }
    }
    let within = found
        .iter()
        .filter(|&&(i, j)| segments[i].2 == segments[j].2)
        .count();
    println!(
        "{} segments; {} crossings within a polygon, {} between polygons",
        segments.len(),
        within,
        found.len() - within
    );
    for &(i, _) in found.iter().take(5) {
        println!("  near {:?}", segments[i].0);
    }
}
