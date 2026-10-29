use wcl_koppen::height;
use wcl_koppen::wcl::Export;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let e = Export::open(a[1].as_ref()).unwrap();
    let elev = e.grid("elevation").unwrap();
    let land = e.grid("land_mask").unwrap();
    let t0 = std::time::Instant::now();
    let h = height::read(a[2].as_ref(), Some((-6000.0, 6000.0))).unwrap();
    eprintln!("read {} x {} in {:?}", h.width, h.height, t0.elapsed());
    let k = h.width / elev.width;
    let (mut sx, mut sy, mut sxx, mut sxy, mut n) = (0f64, 0f64, 0f64, 0f64, 0f64);
    let (mut both, mut vale_only, mut wcl_only) = (0, 0, 0);
    let mut hmax = f32::MIN;
    for r in 0..elev.height { for c in 0..elev.width {
        let i = r * elev.width + c;
        let v = h.block_mean(c * k, r * k, k);
        hmax = hmax.max(v);
        let wl = land.values[i] > 0.0;
        match (wl, v > 0.0) { (true, true) => both += 1, (false, true) => vale_only += 1, (true, false) => wcl_only += 1, _ => {} }
        if wl && v > 0.0 {
            let x = v as f64; let y = elev.values[i] as f64;
            sx += x; sy += y; sxx += x * x; sxy += x * y; n += 1.0;
        }
    }}
    let slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    let icpt = (sy - slope * sx) / n;
    println!("land both {both}, vale only {vale_only}, export only {wcl_only}");
    println!("export elevation = {slope:.3} x vale + {icpt:.1}; vale block max {hmax}");
    let lv = h.levels.iter().copied().max().unwrap();
    let lmin = h.levels.iter().copied().min().unwrap();
    println!("levels {lmin}..{lv}");
}
