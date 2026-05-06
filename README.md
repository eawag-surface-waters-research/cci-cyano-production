# CCI Cyano Production

![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)
![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)

A processing pipeline for extracting per-lake time series and computing phenology metrics from ESA CCI Lakes satellite data. The pipeline downscales global NetCDF imagery to individual lake polygons and fits cubic splines to derive green-up, green-down, peak, and trough dates at the pixel level.

## Overview

The pipeline runs in two sequential stages:

1. **Extract** — clips global satellite NetCDF files to each lake's bounding box and raster mask, producing a per-lake time series NetCDF containing the target variable and quality flags.
2. **Phenology** — reads the extracted time series, fits a smoothing cubic spline per pixel, and extracts phenology metrics (peaks, troughs, green-up/green-down onset, mid, and advanced dates, and data gaps).
3. **Visualization** — loads extracted and phenology outputs to produce spatial metric maps (R², MAD, RMSE, correlation), per-pixel time series plots overlaying raw observations, fitted splines, and cross-dataset comparison; it can also be used directly in Jupyter notebooks for more detailed or custom interactive analysis.

Supported datasets:
- **v2.1** — ESA CCI Lakes v2.1, variable `chla_mean`, QA flag `lwlr_quality_flag`
- **v3.1** — ESA CCI Lakes v3.1, variables `chla` and `phycocyanin`, QA flag `lwlr_quality_flags`

## Setup

### Conda (recommended)

```bash
conda env create -f environment.yml
conda activate cci
```

The environment uses Python 3.11 and installs all required packages via conda-forge, including `csaps`, `scipy`, `geopandas`, `rasterio`, `netCDF4`, `xarray`, and `tqdm`.

## Usage

Run from the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/main.py -f <arg_file> [options]
```

Setting `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1` avoids thread oversubscription when using the `-p pixels` parallelisation mode.

### Arguments

| Flag | Long | Description | Default |
|------|------|-------------|---------|
| `-f` | `--file` | Argument file name (without extension) from `args/` | required |
| `-l` | `--logs` | Write logs to a timestamped file in `logs/` | off |
| `-t` | `--threads` | Number of parallel worker threads | `1` |
| `-p` | `--parallel` | Parallelise over `lakes` or `pixels` | `lake` |
| `-b` | `--batch-size` | Pixels per batch for phenology I/O | `100` |

### Parallelisation modes

- **`-p pixels`** (recommended for large lakes): processes each lake serially but parallelises pixel batches within a lake using `ProcessPoolExecutor`. Workers inherit the in-memory data array via fork copy-on-write, avoiding redundant I/O.
- **`-p lakes`**: processes lakes in parallel, each lake on a single thread. Better for datasets with many small lakes.

### Example commands

Run phenology on the v2.1 chlorophyll-a dataset using 50 pixel threads:
```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/main.py -f v2_chla -p pixels -t 50
```

Run full extraction and phenology on v3.1 phycocyanin with logging:
```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/main.py -f v3_phycocyanin_new -p pixels -t 50 -l
```

## Configuration Files

Each JSON file in `args/` controls one run. All keys except `variable`, `qa`, `shapefile`, `images`, and `out_folder` have defaults and can be omitted.

```json
{
  "variable": "chla_mean",
  "qa": "lwlr_quality_flag",
  "shapefile": "/path/to/lakes.shp",
  "images": "/path/to/NetCDF/imagery/",
  "out_folder": "/path/to/output/",
  "lakes": [5, 15, 6],
  "extract": true,
  "phenology": true,
  "analysis": false,
  "maps": false,
  "pixel_plots": false,
  "comparison": false,
  "comparison_classes": ["chla21", "chla31", "phycocyanin31"],
  "comparison_plot_types": ["chla21 vs chla31", "triple"],
  "background_pts": true,
  "purple_chla21": false,
  "time_splits": [[0, 9999]],
  "aggregation": true,
  "qa_filter": true,
  "spline_min_phase_length": 14,
  "spline_min_relative_amplitude": 0,
  "spline_min_phase_data": 0,
  "spline_data_gap_size": 31,
  "spline_data_gap_size_buffer": 0,
  "spline_subs_peak_win_size": 365,
  "spline_subs_peak_ampl_frac": 0.05
}
```

### Pipeline keys

| Key | Description |
|-----|-------------|
| `variable` | NetCDF variable to extract (`chla_mean`, `chla`, `phycocyanin`) |
| `qa` | Quality flag variable name |
| `shapefile` | Path to the lake boundary shapefile |
| `images` | Root directory of input NetCDF files (searched recursively for `*.nc`) |
| `out_folder` | Root directory for output files |
| `lakes` | List of lake IDs to process; omit or leave empty to process all lakes in the shapefile |
| `start_index` / `end_index` | Process a slice of the shapefile by row index (alternative to `lakes`) |
| `extract` | Run the extraction stage |
| `phenology` | Run the phenology stage |
| `qa_filter` | Exclude pixels where QA flag ≠ 0 (Good) |
| `spline_min_phase_length` | Minimum phase length in days for a valid spline peak/trough |
| `spline_min_relative_amplitude` | Minimum relative amplitude (0–1) for a phase to be retained |
| `spline_min_phase_data` | Minimum number of observations within a phase |
| `spline_data_gap_size` | Minimum gap length (days) to flag as a data gap |
| `spline_data_gap_size_buffer` | Buffer (days) added around flagged data gaps |
| `spline_subs_peak_win_size` | Window size (days) for the substantial-peak amplitude check |
| `spline_subs_peak_ampl_frac` | Amplitude fraction threshold for the substantial-peak check (0.05 retains smaller peaks; 0.35 filters them) |

### Analysis and visualization keys

| Key | Default | Description |
|-----|---------|-------------|
| `analysis` | `false` | Run the post-processing analysis stage |
| `maps` | `false` | Save metric maps (R², MAD, RMSE, correlation, values per pixel) as PNG files |
| `pixel_plots` | `false` | Path to a JSON pixel dictionary file (see below), or `false` to skip pixel plots |
| `comparison` | `false` | Run cross-version comparison plots |
| `comparison_classes` | `["chla21","chla31","phycocyanin31"]` | Dataset labels used for comparison; each must correspond to a configured `PhenologyVisualization` instance |
| `comparison_plot_types` | (all combinations) | Which pairwise or triple comparisons to generate; valid values include `"chla21 vs chla31"`, `"chla21 vs phyco"`, `"chla31 vs phyco"`, `"triple"`, and `*_split` variants |
| `background_pts` | `true` | Show raw or aggregated observations as a scatter background behind peak/trough stem plots |
| `purple_chla21` | `false` | Colour v2.1 chl-a data in purple instead of light green when overlaying with other datasets |
| `time_splits` | `[[0, 9999]]` | List of `[start, end]` year ranges to compute metrics and generate plots for. `0` means the earliest available year; `9999` means the latest. A single entry produces one panel per output; two entries produce a side-by-side split layout; more entries create a grid. |
| `aggregation` | `true` | Use 3×3 neighbourhood median values as the scatter background instead of the raw pixel time series |

#### Pixel dictionary format

`pixel_plots` must point to a JSON file in `pixel_dicts/`. The file maps lake IDs (as strings) to lists of `[i, j]` pixel index pairs:

```json
{
  "327":        [[23, 30]],
  "300013962":  [[8, 3]],
  "300013966":  [[7, 3]]
}
```

If a lake ID is not present in the dictionary, pixel plots for that lake are skipped with a warning. The `comparison` stage also requires `pixel_plots` to be set.

## Visualization and Analysis (`scripts/visualization.py`)

`visualization.py` provides the `PhenologyVisualization` class for interactive and static exploration of extracted and phenology NetCDF outputs in Jupyter notebooks. It mirrors the `PhenologyEDA` class in `scripts/analysis.py` and exposes the same interface.

### Setup

```python
from scripts.visualization import PhenologyVisualization

PhenologyVisualization.set_shapefile_path("/path/to/lakes.shp")

vis = PhenologyVisualization(
    extract_path="/path/to/extract/chla/12345.nc",
    phenology_path="/path/to/phenology/chla/12345.nc",
)
```

`set_shapefile_path` must be called once at the class level before instantiating any object.

### Key methods

| Method | Description |
|--------|-------------|
| `load_pixel_data(i, j)` | Return a datetime-indexed Series of valid observations for pixel `(i, j)` |
| `pixel_map(lat, lon, ax)` | Grayscale coverage map with a selected pixel marked |
| `interactive_pixel_map(ax)` | Clickable coverage map that prints `(lat_idx, lon_idx)` on click (requires `%matplotlib widget`) |
| `single_plot(lat, lon, ax, start, end)` | Observations, smoothed spline, and all phenological events for one pixel |
| `split_plot(lat, lon, ax0, ax1, start0, end0, start1, end1)` | Two side-by-side `single_plot` panels for different year windows |
| `full_plot(lat, lon, ax)` | `single_plot` over the complete available time range |
| `extrema_plot(lat, lon, ax, peak, start, end)` | Stem plot of peaks or troughs; pass `peak=False` for troughs |
| `extrema_comparison(other1, lat, lon, ax, …, other2)` | Overlay extrema plots from two or three `PhenologyVisualization` objects |
| `single_plot_background(lat, lon, ax, fig, …)` | Like `single_plot` with scatter coloured by QA flag and a QA colorbar |
| `metric_map(metric_scores, metric_str, fig, ax, …)` | Spatial heatmap of any `{(i,j): value}` metric dict; pixels outside the 1 km-inset boundary are masked |
| `interactive_metric_map(metric_scores, metric_str, fig, ax)` | Clickable version of `metric_map` |
| `time_map(fig, ax, year, peaks, max)` | Map of peak or green-up day-of-year for a given year; pixels outside the 1 km-inset boundary are masked |
| `single_day_map(date)` | Spatial map of the target variable for a single observation date; pixels outside the 1 km-inset boundary are masked. `date` must be a UTC-aware `datetime` matching an entry in the extract time axis — use `series.idxmax()` or another index value from `load_pixel_data` |
| `r2_scores(time_split)` | `{(i,j): R²}` for all valid pixels; cached to CSV |
| `MAD_scores(time_split)` | `{(i,j): MAD}` for all valid pixels; cached to CSV |
| `RMSE_scores(time_split)` | `{(i,j): RMSE}` for all valid pixels; cached to CSV |
| `correlation_scores(time_split)` | `{(i,j): Pearson r}` for all valid pixels; cached to CSV |
| `values_per_pixel(time_split)` | `{(i,j): count}` of valid observations; cached to CSV |
| `spatial_aggregation()` | Compute 3×3 neighbourhood medians for all pixels and timesteps; cached to CSV |

The metric methods (`r2_scores`, `MAD_scores`, `RMSE_scores`, `correlation_scores`, `values_per_pixel`) each accept a `time_split` argument: a single-element list containing one `[start, end]` year pair, e.g. `[[2003, 2012]]` or `[[0, 9999]]` for the full series. Results are cached to CSV; the cache filename encodes the time window so different windows are stored independently.

All metric and aggregation computations are cached to CSV on first call and loaded from cache on subsequent calls. Interactive plot methods require an interactive Matplotlib backend (`%matplotlib widget`).

## Output Format

### Extraction output (`extract/{variable}/{lake_id}.nc`)

| Variable | Dimensions | Description |
|----------|-----------|-------------|
| `{variable}` | `(time, lat, lon)` | Extracted pixel values; fill value `-9999` |
| `{qa}` | `(time, lat, lon)` | QA flag: 0=Good, 1=Fair, 2=Poor, 3=No data |
| `summary` | `(lat, lon)` | Count of valid (non-fill) observations per pixel |
| `time` | `(time,)` | Unix timestamps |
| `lat`, `lon` | 1-D | Coordinates of the lake sub-grid |

### Phenology output (`phenology/{variable}/{lake_id}.nc`)

All variables are shaped `(lat, lon, record)` where the `record` dimension is unlimited and grows to the maximum number of events across all pixels.

| Variable | Description |
|----------|-------------|
| `smoothing_parameter` | Optimal csaps smoothing parameter per pixel |
| `pks_x` / `pks_y` / `pks_qa` | Peak time (Unix), value, and QA |
| `trgs_x` / `trgs_y` / `trgs_qa` | Trough time (Unix), value, and QA |
| `green_up_onset_x/y`, `green_up_mid_x/y`, `green_up_advanced_x/y` | Green-up onset, mid, and advanced dates and values |
| `green_down_onset_x/y`, `green_down_mid_x/y`, `green_down_advanced_x/y` | Green-down onset, mid, and advanced dates and values |
| `data_gap_start` / `data_gap_end` | Data gap start and end (Unix) |

Global attributes store the run parameters used to produce the file.

### Visualization and analysis output

Metric CSVs and spatial aggregation values are written relative to the phenology file's root (`out_folder`). Analysis plots are written to a separate lake analysis folder passed at runtime.

```
out_folder/
├── extract/
│   └── {variable}/
│       └── {lake_id}.nc
├── phenology/
│   └── {variable}/
│       ├── {lake_id}.nc
│       └── checkpoints/
│           └── {lake_id}/
│               └── bs{batch_size}/
│                   └── *.npy          # temporary; deleted after successful write
└── calculated_values/
    ├── metrics/
    │   └── {metric_name}/             # r2 | MAD | RMSE | correlation | values_per_pixel
    │       ├── full_ts.csv                          # time_split [0, 9999]
    │       ├── ts_end_{year}.csv                    # time_split [0, year]
    │       ├── ts_start_{year}.csv                  # time_split [year, 9999]
    │       └── ts_start_{year}_end_{year}.csv       # time_split [start, end]
    └── spatial_aggregation_values/
        └── aggregation_background_values.csv

lake_analysis_folder/
└── {lake_str}/
    └── plots/
        ├── metric_maps/
        │   ├── {variable}_v{version}_{metric}_full_ts.png      # single time_split [0, 9999]
        │   ├── {variable}_v{version}_{metric}_split_ts.png     # two time_splits
        │   └── {variable}_v{version}_{metric}_{n}_split_ts.png # n > 2 time_splits
        └── pixel_plots/
            ├── aggregated/
            │   └── {i}_{j}/
            │       ├── location.png
            │       ├── {variable}_v{version}_full_ts.png
            │       ├── {variable}_v{version}_split_ts.png
            │       └── {variable}_v{version}_peaks_*.png
            ├── not_aggregated/
            │   └── {i}_{j}/           # same structure as aggregated/
            └── summaries/
                └── {i}_{j}/
                    └── summary.txt
```

## Fault Tolerance and Restart

The phenology stage saves intermediate results as `.npy` checkpoint files in `data/{version}/phenology/{variable}/checkpoints/{lake_id}/bs{batch_size}/`. If a run is interrupted, re-running the same command will skip completed batches and resume from where it left off. Checkpoints are deleted after the final NetCDF is successfully written.

Both stages skip lakes whose output file already exists, so it is safe to rerun the command after adding new lakes to the `lakes` list.
