//! The Köppen-Geiger classification.
//!
//! The thresholds are those of Peel et al. (2007), with the As class. The
//! order of the groups and the rule for a dry summer are those of the
//! classes of World Climate Lab: B comes before E, and a dry summer needs
//! 40 % of the precipitation in winter.

/// The 31 Köppen-Geiger classes, with their codes and names.
pub const CLASSES: [(&str, &str); 31] = [
    ("Af", "Tropical Rainforest"),
    ("Am", "Tropical Monsoon"),
    ("Aw", "Savanna, Dry Winter"),
    ("As", "Savanna, Dry Summer"),
    ("BWh", "Hot Desert"),
    ("BWk", "Cold Desert"),
    ("BSh", "Hot Steppe"),
    ("BSk", "Cold Steppe"),
    ("Cfa", "Humid Subtropical"),
    ("Cfb", "Oceanic"),
    ("Cfc", "Subpolar Oceanic"),
    ("Csa", "Hot Mediterranean"),
    ("Csb", "Warm Mediterranean"),
    ("Csc", "Cold Mediterranean"),
    ("Cwa", "Dry-Winter Subtropical"),
    ("Cwb", "Subtropical Highland"),
    ("Cwc", "Cold Highland"),
    ("Dfa", "Hot Humid Continental"),
    ("Dfb", "Warm Humid Continental"),
    ("Dfc", "Subarctic"),
    ("Dfd", "Extremely Cold Subarctic"),
    ("Dsa", "Hot Dry-Summer Continental"),
    ("Dsb", "Warm Dry-Summer Continental"),
    ("Dsc", "Dry-Summer Subarctic"),
    ("Dsd", "Extremely Cold Dry-Summer Subarctic"),
    ("Dwa", "Hot Dry-Winter Continental"),
    ("Dwb", "Warm Dry-Winter Continental"),
    ("Dwc", "Dry-Winter Subarctic"),
    ("Dwd", "Extremely Cold Dry-Winter Subarctic"),
    ("ET", "Tundra"),
    ("EF", "Ice Cap"),
];

/// The usual map color of each class of [`CLASSES`], from Beck et al. (2018).
pub const COLORS: [[u8; 3]; 31] = [
    [0, 0, 255],
    [0, 120, 255],
    [70, 170, 250],
    [120, 200, 250],
    [255, 0, 0],
    [255, 150, 150],
    [245, 165, 0],
    [255, 220, 100],
    [200, 255, 80],
    [100, 255, 80],
    [50, 200, 0],
    [255, 255, 0],
    [200, 200, 0],
    [150, 150, 0],
    [150, 255, 150],
    [100, 200, 100],
    [50, 150, 50],
    [0, 255, 255],
    [55, 200, 255],
    [0, 125, 125],
    [0, 70, 95],
    [255, 0, 255],
    [200, 0, 200],
    [150, 50, 150],
    [150, 100, 150],
    [170, 175, 255],
    [90, 120, 220],
    [75, 80, 180],
    [50, 0, 135],
    [178, 178, 178],
    [102, 102, 102],
];

/// An index into [`CLASSES`].
pub type Class = u8;

/// The parameters of the classification.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rules {
    /// The temperature of the coldest month in °C that separates C from D.
    pub cd_boundary: f32,
}

impl Default for Rules {
    fn default() -> Self {
        Self { cd_boundary: 0.0 }
    }
}

/// Returns the index of the class with `code`.
pub fn class(code: &str) -> Option<Class> {
    CLASSES
        .iter()
        .position(|(c, _)| *c == code)
        .map(|i| i as Class)
}

/// Classifies one place from its monthly mean temperatures in °C and its
/// monthly precipitation in mm, for months 1 to 12. Summer is the warmer of
/// the two halves of the year, months 4 to 9 or months 10 to 3.
pub fn classify(t: &[f32; 12], p: &[f32; 12], rules: Rules) -> Class {
    let t_hot = t.iter().copied().fold(f32::MIN, f32::max);
    let t_cold = t.iter().copied().fold(f32::MAX, f32::min);
    let t_mean = t.iter().sum::<f32>() / 12.0;
    let p_total: f32 = p.iter().sum();
    let p_dry = p.iter().copied().fold(f32::MAX, f32::min);

    let april_warmer = t[3..9].iter().sum::<f32>() >= t.iter().sum::<f32>() / 2.0;
    let summer = |m: usize| (3..9).contains(&m) == april_warmer;
    let (mut p_summer, mut p_winter) = (0.0, 0.0);
    let (mut s_dry, mut s_wet) = (f32::MAX, f32::MIN);
    let (mut w_dry, mut w_wet) = (f32::MAX, f32::MIN);
    for (m, &v) in p.iter().enumerate() {
        if summer(m) {
            p_summer += v;
            s_dry = s_dry.min(v);
            s_wet = s_wet.max(v);
        } else {
            p_winter += v;
            w_dry = w_dry.min(v);
            w_wet = w_wet.max(v);
        }
    }

    let threshold = if p_winter >= 0.7 * p_total {
        2.0 * t_mean
    } else if p_summer >= 0.7 * p_total {
        2.0 * t_mean + 28.0
    } else {
        2.0 * t_mean + 14.0
    };
    if p_total < 10.0 * threshold {
        let code = match (p_total < 5.0 * threshold, t_mean >= 18.0) {
            (true, true) => "BWh",
            (true, false) => "BWk",
            (false, true) => "BSh",
            (false, false) => "BSk",
        };
        return code_class(code);
    }
    if t_hot < 10.0 {
        return code_class(if t_hot > 0.0 { "ET" } else { "EF" });
    }
    if t_cold >= 18.0 {
        let code = if p_dry >= 60.0 {
            "Af"
        } else if p_dry >= 100.0 - p_total / 25.0 {
            "Am"
        } else if s_dry < w_dry {
            "As"
        } else {
            "Aw"
        };
        return code_class(code);
    }

    let group = if t_cold > rules.cd_boundary { 'C' } else { 'D' };
    // A dry summer also needs 40 % of the precipitation in winter, and it wins
    // over a dry winter.
    let dry_summer = s_dry < 40.0 && s_dry < w_wet / 3.0 && p_winter >= 0.4 * p_total;
    let dry_winter = w_dry < s_wet / 10.0;
    let season = if dry_summer {
        's'
    } else if dry_winter {
        'w'
    } else {
        'f'
    };
    let warm_months = t.iter().filter(|&&v| v >= 10.0).count();
    let heat = if t_hot >= 22.0 {
        'a'
    } else if warm_months >= 4 {
        'b'
    } else if group == 'D' && t_cold < -38.0 {
        'd'
    } else {
        'c'
    };
    code_class(&format!("{group}{season}{heat}"))
}

fn code_class(code: &str) -> Class {
    class(code).expect("each code of the rules is a class")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn code(t: [f32; 12], p: [f32; 12]) -> &'static str {
        CLASSES[classify(&t, &p, Rules::default()) as usize].0
    }

    #[test]
    fn rainforest() {
        assert_eq!(code([27.0; 12], [200.0; 12]), "Af");
    }

    #[test]
    fn hot_desert() {
        assert_eq!(code([30.0; 12], [2.0; 12]), "BWh");
    }

    #[test]
    fn ice_cap_and_tundra() {
        assert_eq!(code([-20.0; 12], [10.0; 12]), "EF");
        let mut t = [-20.0; 12];
        t[6] = 5.0;
        assert_eq!(code(t, [30.0; 12]), "ET");
    }

    // London, Heathrow, 1991 to 2020.
    #[test]
    fn oceanic() {
        let t = [
            5.2, 5.5, 7.6, 10.0, 13.3, 16.5, 18.7, 18.5, 15.7, 12.0, 8.0, 5.5,
        ];
        let p = [
            56.0, 41.0, 37.0, 42.0, 46.0, 46.0, 45.0, 50.0, 47.0, 66.0, 62.0, 56.0,
        ];
        assert_eq!(code(t, p), "Cfb");
    }

    // A Mediterranean climate in the south swaps the seasons.
    #[test]
    fn summer_follows_the_temperature() {
        let t = [
            22.0, 22.0, 20.0, 17.0, 14.0, 12.0, 11.0, 12.0, 14.0, 16.0, 19.0, 21.0,
        ];
        let p = [
            10.0, 10.0, 20.0, 50.0, 90.0, 100.0, 100.0, 80.0, 50.0, 30.0, 15.0, 10.0,
        ];
        assert_eq!(code(t, p), "Csa");
        // With the rain in summer, the same total makes the place dry.
        let mut wet_summer = p;
        wet_summer.rotate_left(6);
        assert_eq!(code(t, wet_summer), "BSk");
    }

    #[test]
    fn dry_summer_needs_rain_in_winter() {
        let t = [
            -10.9, -8.4, -4.5, -2.3, 7.6, 17.5, 22.5, 17.3, 8.2, 0.3, -6.1, -7.4,
        ];
        let p = [
            11.3, 14.8, 11.2, 22.7, 50.9, 28.2, 5.4, 3.5, 13.0, 15.0, 14.4, 13.6,
        ];
        assert_eq!(code(t, p), "Dfa");
    }

    #[test]
    fn cd_boundary() {
        let mut t = [
            -2.0, 0.0, 5.0, 10.0, 15.0, 19.0, 21.0, 20.0, 15.0, 10.0, 4.0, 0.0,
        ];
        let p = [50.0; 12];
        assert_eq!(code(t, p), "Dfb");
        let rules = Rules { cd_boundary: -3.0 };
        assert_eq!(CLASSES[classify(&t, &p, rules) as usize].0, "Cfb");
        t[0] = -40.0;
        assert_eq!(code(t, p), "Dfb");
    }
}
