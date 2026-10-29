//! Compares the classes with the classes of a real export. The variable
//! `WCL_EXPORT` names the export folder. Without it, the test does nothing.

use wcl_koppen::koppen::{self, Rules};
use wcl_koppen::wcl::Export;
use wcl_koppen::{NO_CLASS, classify_grid};

#[test]
fn matches_the_classes_of_the_export() {
    let Some(dir) = std::env::var_os("WCL_EXPORT") else {
        eprintln!("WCL_EXPORT is not set, so the test does nothing");
        return;
    };
    let export = Export::open(dir.as_ref()).unwrap();
    let climate = export.climate().unwrap();
    let reference = export.classes().unwrap();
    let ours = classify_grid(&climate, Rules::default());
    let (mut total, mut same) = (0, 0);
    for (i, &c) in ours.iter().enumerate() {
        let value = reference.grid.values[i] as u32;
        let theirs = reference
            .codes
            .iter()
            .find(|(v, _, system)| *v == value && system == "Koppen-Geiger")
            .and_then(|(_, code, _)| koppen::class(code));
        if c == NO_CLASS || theirs.is_none() {
            continue;
        }
        total += 1;
        same += usize::from(Some(c) == theirs);
    }
    let share = same as f64 / total as f64;
    assert!(share >= 0.995, "{same} of {total} land cells match");
}
