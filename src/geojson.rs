//! The GeoJSON output.

use std::io::{self, Write};

use crate::koppen::{CLASSES, COLORS};
use crate::trace::Polygon;

/// The area of a ring in km², positive for a counterclockwise ring.
pub fn ring_area_km2(ring: &[[f64; 2]], radius_km: f64) -> f64 {
    // A ring that is straight on a cylindrical equal-area map.
    let sum: f64 = ring
        .windows(2)
        .map(|p| {
            let (x0, y0) = (p[0][0].to_radians(), p[0][1].to_radians().sin());
            let (x1, y1) = (p[1][0].to_radians(), p[1][1].to_radians().sin());
            x0 * y1 - x1 * y0
        })
        .sum();
    sum / 2.0 * radius_km * radius_km
}

/// The area of a polygon in km².
pub fn area_km2(polygon: &Polygon, radius_km: f64) -> f64 {
    ring_area_km2(&polygon.exterior, radius_km)
        + polygon
            .holes
            .iter()
            .map(|h| ring_area_km2(h, radius_km))
            .sum::<f64>()
}

/// Writes a FeatureCollection with one Polygon feature for each zone.
pub fn write(out: &mut impl Write, polygons: &[Polygon], radius_km: f64) -> io::Result<()> {
    writeln!(out, "{{\"type\":\"FeatureCollection\",\"features\":[")?;
    for (i, polygon) in polygons.iter().enumerate() {
        let (code, name) = CLASSES[polygon.class as usize];
        let group = &code[..1];
        let [r, g, b] = COLORS[polygon.class as usize];
        let area = area_km2(polygon, radius_km);
        write!(
            out,
            "{{\"type\":\"Feature\",\"properties\":{{\"code\":\"{code}\",\"name\":{},\
             \"system\":\"Koppen-Geiger\",\"group\":\"{group}\",\"color\":\"#{r:02x}{g:02x}{b:02x}\",\
             \"area_km2\":{area:.1}}},\
             \"geometry\":{{\"type\":\"Polygon\",\"coordinates\":[",
            serde_json::to_string(name).expect("a name is a string"),
        )?;
        let rings = std::iter::once(&polygon.exterior).chain(&polygon.holes);
        for (k, ring) in rings.enumerate() {
            if k > 0 {
                write!(out, ",")?;
            }
            write!(out, "[")?;
            for (j, p) in ring.iter().enumerate() {
                if j > 0 {
                    write!(out, ",")?;
                }
                write!(out, "[{},{}]", round(p[0]), round(p[1]))?;
            }
            write!(out, "]")?;
        }
        let end = if i + 1 < polygons.len() { "," } else { "" };
        writeln!(out, "]}}}}{end}")?;
    }
    writeln!(out, "]}}")
}

/// A coordinate with six decimals, about 0.1 m.
fn round(v: f64) -> f64 {
    (v * 1e6).round() / 1e6
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hemisphere_area() {
        let ring = [
            [-180.0, 0.0],
            [180.0, 0.0],
            [180.0, 90.0],
            [-180.0, 90.0],
            [-180.0, 0.0],
        ];
        let area = ring_area_km2(&ring, 1.0);
        assert!((area - 2.0 * std::f64::consts::PI).abs() < 1e-9);
    }

    #[test]
    fn valid_json() {
        let polygon = Polygon {
            class: crate::koppen::class("Cfb").unwrap(),
            exterior: vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 0.0]],
            holes: vec![],
        };
        let mut out = Vec::new();
        write(&mut out, &[polygon.clone(), polygon], 6371.0).unwrap();
        let json: serde_json::Value = serde_json::from_slice(&out).unwrap();
        let features = json["features"].as_array().unwrap();
        assert_eq!(features.len(), 2);
        assert_eq!(features[0]["properties"]["code"], "Cfb");
        assert_eq!(features[0]["properties"]["name"], "Oceanic");
        assert_eq!(features[0]["properties"]["color"], "#64ff50");
        assert_eq!(features[0]["geometry"]["coordinates"][0][2][1], 1.0);
    }
}
