//! Draws a GeoJSON file of zones to an equirectangular PNG, for a look at
//! the output.
//!
//! ```sh
//! cargo run --release --example render -- zones.geojson zones.png 4096 [WEST,SOUTH,EAST,NORTH]
//! ```

use std::fs::File;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 4 {
        eprintln!("usage: render ZONES.geojson OUT.png WIDTH [WEST,SOUTH,EAST,NORTH]");
        std::process::exit(2);
    }
    let json: serde_json::Value =
        serde_json::from_reader(File::open(&args[1]).expect("the zones file opens")).unwrap();
    let bbox: [f64; 4] = match args.get(4) {
        Some(b) => {
            let v: Vec<f64> = b.split(',').map(|x| x.parse().unwrap()).collect();
            [v[0], v[1], v[2], v[3]]
        }
        None => [-180.0, -90.0, 180.0, 90.0],
    };
    let width: usize = args[3].parse().unwrap();
    let height = ((width as f64) * (bbox[3] - bbox[1]) / (bbox[2] - bbox[0])).round() as usize;
    let to_px = |p: &serde_json::Value| -> (f64, f64) {
        let (lon, lat) = (p[0].as_f64().unwrap(), p[1].as_f64().unwrap());
        (
            (lon - bbox[0]) / (bbox[2] - bbox[0]) * width as f64,
            (bbox[3] - lat) / (bbox[3] - bbox[1]) * height as f64,
        )
    };
    let mut img = vec![[20u8, 40, 70]; width * height];
    let features = json["features"].as_array().unwrap();
    let mut outlines = Vec::new();
    for f in features {
        let color = f["properties"]["color"].as_str().unwrap();
        let rgb = [1, 3, 5].map(|i| u8::from_str_radix(&color[i..i + 2], 16).unwrap());
        let rings: Vec<Vec<(f64, f64)>> = f["geometry"]["coordinates"]
            .as_array()
            .unwrap()
            .iter()
            .map(|r| r.as_array().unwrap().iter().map(to_px).collect())
            .collect();
        // Even-odd scanline fill of all the rings.
        for y in 0..height {
            let yc = y as f64 + 0.5;
            let mut xs = Vec::new();
            for ring in &rings {
                for e in ring.windows(2) {
                    let ((x0, y0), (x1, y1)) = (e[0], e[1]);
                    if (y0 > yc) != (y1 > yc) {
                        xs.push(x0 + (yc - y0) / (y1 - y0) * (x1 - x0));
                    }
                }
            }
            xs.sort_by(f64::total_cmp);
            for pair in xs.as_chunks::<2>().0 {
                let a = (pair[0] - 0.5).ceil().max(0.0) as usize;
                let b = ((pair[1] - 0.5).floor() as i64).min(width as i64 - 1);
                for x in a as i64..=b {
                    img[y * width + x as usize] = rgb;
                }
            }
        }
        outlines.extend(rings);
    }
    for ring in &outlines {
        for e in ring.windows(2) {
            let ((x0, y0), (x1, y1)) = (e[0], e[1]);
            let n = ((x1 - x0).abs().max((y1 - y0).abs()).ceil() as usize).max(1);
            for k in 0..=n {
                let t = k as f64 / n as f64;
                let (x, y) = (x0 + t * (x1 - x0), y0 + t * (y1 - y0));
                if x >= 0.0 && y >= 0.0 && (x as usize) < width && (y as usize) < height {
                    img[y as usize * width + x as usize] = [30, 30, 30];
                }
            }
        }
    }
    let out = File::create(&args[2]).unwrap();
    let mut encoder = png::Encoder::new(out, width as u32, height as u32);
    encoder.set_color(png::ColorType::Rgb);
    encoder.set_depth(png::BitDepth::Eight);
    let data: Vec<u8> = img.into_iter().flatten().collect();
    encoder
        .write_header()
        .unwrap()
        .write_image_data(&data)
        .unwrap();
}
