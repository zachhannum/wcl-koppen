# wcl-koppen

`wcl-koppen` reads a climate export of World Climate Lab and classifies it into Köppen-Geiger climate zones. The goal is a vector file of the zones in longitude and latitude, at the detail of a finer heightmap, that a GIS such as Vale imports as a polygon layer.

## Current state

The tool reads an export and classifies each land cell on the grid of the export. The command `check` compares the result with the `climate_class.png` of the export:

```sh
cargo run --release -- check path/to/climate-data
```

On a test export, 99.6 % of the land cells match.

The correction of the temperature for a finer heightmap and the vector output are not done yet.

## Input

An export folder has the manifest `climate_data.json`. Each grid is an equirectangular greyscale PNG for the whole globe, with a JSON sidecar that gives its encoding. The tool reads these grids:

- `temperature_month_01.png` to `temperature_month_12.png`: the monthly mean temperature in °C, at the height of the ground of the export.
- `precipitation_month_01.png` to `precipitation_month_12.png`: the monthly precipitation in mm.
- `elevation.png`: the height of the ground in m.
- `land_mask.png`: land and sea.

## Classification

The thresholds are those of Peel et al. (2007), with the As class. These rules match the classes of World Climate Lab:

- Summer is the warmer half of the year, months 4 to 9 or months 10 to 3.
- The B group comes before the E group.
- A dry summer (s) needs 40 % of the precipitation in winter, and it comes before a dry winter (w).

The C and D groups split at a coldest month of 0 °C. The option `--cd-boundary -3` uses −3 °C.

## Tests

The variable `WCL_EXPORT` names an export folder for the test against a real export:

```sh
WCL_EXPORT=path/to/climate-data cargo test --release
```

## License

MIT OR Apache-2.0.
