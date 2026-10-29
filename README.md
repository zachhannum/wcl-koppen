# wcl-koppen

`wcl-koppen` makes Köppen-Geiger climate zones from a climate export of World Climate Lab. It classifies the climate at the detail of a finer heightmap and writes the zones as polygons in longitude and latitude. A GIS such as Vale or QGIS imports the GeoJSON file as a polygon layer.

## Usage

This command writes the zones of an export for a heightmap from Vale:

```sh
cargo run --release -- zones path/to/climate-data --heightmap world.tif --out zones.geojson
```

The heightmap is an equirectangular map of the whole globe, twice as wide as it is tall, with 16-bit levels. Land is above 0 m. A GeoTIFF from the Vale command `--export-heightmap` holds the heights of its levels. A PNG holds no heights, so the option `--height-range` gives the heights of levels 0 and 65,535:

```sh
cargo run --release -- zones path/to/climate-data --heightmap world.png --height-range=-6000,6000 --out zones.geojson
```

The command `check` classifies the export on its own grid and compares the result with the `climate_class.png` of the export:

```sh
cargo run --release -- check path/to/climate-data
```

## Steps

1. The tool reads the monthly temperature and precipitation, the elevation, and the land mask of the export.
2. For each cell of the heightmap, the tool blends the four nearest land cells of the export. A cell with no land cell of the export around it takes the nearest land cell.
3. The tool changes each monthly temperature by the lapse rate times the difference between the heightmap and the elevation of the export. The default lapse rate is 6.5 °C per km. The precipitation does not change.
4. The tool classifies each land cell.
5. A zone below the minimum area joins the zones around it. A zone is a set of cells of one class that join at their sides.
6. The tool traces the boundaries of the zones. It smooths each boundary between two zones one time, so the two zones share the same points.
7. The tool simplifies the boundaries and writes one GeoJSON feature for each polygon.

A 32,768 × 16,384 heightmap takes about 20 seconds and 2 GB of memory. The option `--step` makes one cell of the zones from a block of cells of the heightmap.

## Output

The file is a FeatureCollection of Polygon features, as RFC 7946 describes. Each feature has these properties:

| Property | Example |
| --- | --- |
| `code` | `Cfb` |
| `name` | `Oceanic` |
| `system` | `Koppen-Geiger` |
| `group` | `C` |
| `color` | `#64ff50`, the usual map color of the class |
| `area_km2` | `12345.6` |

A zone across the 180° meridian is two polygons, with edges along the meridian. A zone at a pole has an edge along the latitude of the pole.

## Height scale

The temperatures of the export include the effect of the heights that World Climate Lab used. The command `zones` compares the elevation of the export with the heightmap. If the scales are not the same, it shows a warning. For example, World Climate Lab can map the brightest level of a heightmap to 8,849 m, but Vale maps it to the top of the range of the height layer. The tool then corrects the temperatures to the heights of the heightmap, but the winds and the precipitation of the export still come from the other heights.

## Classification

The thresholds are those of Peel et al. (2007), with the As class. These rules match the classes of World Climate Lab:

- Summer is the warmer half of the year, months 4 to 9 or months 10 to 3.
- The B group comes before the E group.
- A dry summer (s) needs 40 % of the precipitation in winter, and it comes before a dry winter (w).

On a test export, 99.6 % of the land cells match the classes of the export. The C and D groups split at a coldest month of 0 °C. The option `--cd-boundary -3` uses −3 °C.

## Options

| Option | Default |
| --- | --- |
| `--lapse-rate` | 6.5 °C per km |
| `--cd-boundary` | 0 °C |
| `--min-area-km2` | 25 cells at the equator |
| `--smooth-km` | 1 cell at the equator |
| `--tolerance-km` | 0.25 cells at the equator |
| `--step` | 1 |
| `--radius-km` | 6,371 km |

## Development

The example `render` draws a GeoJSON file of zones to a PNG. The example `validate` counts the crossings of the edges of the zones.

```sh
cargo run --release --example render -- zones.geojson zones.png 4096
cargo run --release --example validate -- zones.geojson
```

The variable `WCL_EXPORT` names an export folder for the test against a real export:

```sh
WCL_EXPORT=path/to/climate-data cargo test --release
```

## License

MIT OR Apache-2.0.
