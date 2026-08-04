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
  "timing_plots": false,
  "comparison": false,
  "comparison_classes": ["chla21", "chla3", "phycocyanin3"],
  "comparison_plot_types": ["chla21 vs chla3", "chla21 vs phyco", "chla3 vs phyco", "triple"],
  "background_pts": true,
  "purple_chla21": true,
  "ratio_qa_source": "self",
  "kde_qa": null,
  "time_splits": [[0, 9999]],
  "provenance": false,
  "aggregation": true,
  "aggregation_format": "csv",
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
| `timing_plots` | `false` | Save lake-wide bloom-timing outputs to `plots/timing_plots/`: DOY panels (peaks and green-up mid), a lake-wide year×quarter heatmap, one KDE bloom-timing plot per `kde_qa` entry, and lake-wide QA boxplots (see [Bloom-timing KDE](#bloom-timing-kde)) |
| `comparison` | `false` | Run cross-version comparison plots |
| `comparison_classes` | `["chla21","chla3","phycocyanin3"]` | Dataset labels used for comparison; each must correspond to a configured `PhenologyVisualization` instance. Currently `chla3`/`phycocyanin3` resolve to the `v3.0` folder (see `_class_paths` in `main.py`) |
| `comparison_plot_types` | `["chla21 vs chla3","chla21 vs phyco","chla3 vs phyco","triple"]` | Which pairwise or triple comparisons to generate |
| `background_pts` | `true` | Show raw or aggregated observations as a scatter background behind peak/trough stem plots |
| `purple_chla21` | `true` | Colour v2.1 chl-a data in purple instead of light green when overlaying with other datasets |
| `ratio_qa_source` | `"self"` | For ratio boxplots (`qa_boxplot`/`qa_boxplot_lake` with `other` set, and comparison plots): which side's QA flag to group by — `"self"`, `"other"`, or `"matched"` (keeps only pairs where both sides agree on QA) |
| `kde_qa` | `null` | List of QA filters, one KDE plot per entry, e.g. `[[0], [0,1], null]` (`null` = all QA levels). `null`/unset produces a single all-levels plot. Only used by the `timing_plots` stage |
| `time_splits` | `[[0, 9999]]` | List of `[start, end]` year ranges to compute metrics and generate plots for. `0` means the earliest available year; `9999` means the latest. A single entry produces one panel per output; two entries produce a side-by-side split layout; more entries create a grid. |
| `aggregation` | `true` | Use 3×3 neighbourhood median values as the scatter background instead of the raw pixel time series |
| `aggregation_format` | `"csv"` | On-disk cache format for `spatial_aggregation()`: `"csv"` (long-format table, one row per timestep × pixel) or `"netcdf"` (compressed `(time, pixel)` array, chunked for fast single-pixel reads — smaller and faster for large lakes). Set globally via `PhenologyVisualization.set_aggregation_format()` when not going through `main.py` |
| `provenance` | `false` | Append a record of this run to `provenance_logs/` (see [Provenance](#provenance)) |

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
| `single_plot_insitu(lat, lon, ax, insitu_df, …)` | `single_plot` overlaid with in-situ observations (e.g. Diemark data) matched by date and station |
| `single_years_plot(lat, lon, years, ncol, nrow, annotation, ylim)` | Grid of `single_plot` panels, one per calendar year, for a single pixel |
| `qa_boxplot(lat, lon, ax, metric, start, end, other, tolerance_days, qa_source)` | Boxplot of peak/trough values grouped by QA level for one pixel; pass `other` to plot the value-ratio between two instances instead. Annotated with pairwise Mann-Whitney U significance brackets (`*`/`**`/`***`/`ns`) |
| `qa_boxplot_lake(ax, metric, start, end, other, tolerance_days, qa_source)` | Same as `qa_boxplot`, pooled across every valid pixel in the lake |
| `pair_phenology_events(other, lat, lon, metric, tolerance_days)` | Match this instance's peak/trough events against another instance's at the same pixel, within `tolerance_days` |
| `plot_background_ratio_timeseries(other, ax, lat, lon, color)` | Scatter of the per-observation ratio (self / other) between two instances at one pixel, over time |
| `plot_background_ratio_v_self(other, ax, lat, lon, color)` | Same ratio, plotted against self's own value instead of time |
| `yearly_heatmap_pixel(lat, lon, color_scheme, show_gaps, qa, mask_cells, mask_color)` | Bivariate year×quarter heatmap of peak/trough counts for one pixel |
| `yearly_heatmap_lake(color_scheme, qa, mask_cells, mask_color)` | Bivariate year×quarter heatmap of peak/trough fractions pooled across the whole lake |
| `metric_map(metric_scores, metric_str, fig, ax, …)` | Spatial heatmap of any `{(i,j): value}` metric dict; pixels outside the 1 km-inset boundary are masked |
| `interactive_metric_map(metric_scores, metric_str, fig, ax)` | Clickable version of `metric_map` |
| `time_map(fig, ax, year, peaks, max)` | Map of peak or green-up day-of-year for a given year; pixels outside the 1 km-inset boundary are masked |
| `time_map_panel(years, nrow, ncol, peaks, max)` | Grid of `time_map` panels, one per year, with a shared colorbar and lake-outline legend |
| `single_day_map(date)` | Spatial map of the target variable for a single observation date; pixels outside the 1 km-inset boundary are masked. `date` must be a UTC-aware `datetime` matching an entry in the extract time axis — use `series.idxmax()` or another index value from `load_pixel_data` |
| `r2_scores(time_split)` | `{(i,j): R²}` for all valid pixels; cached to CSV |
| `MAD_scores(time_split)` | `{(i,j): MAD}` for all valid pixels; cached to CSV |
| `RMSE_scores(time_split)` | `{(i,j): RMSE}` for all valid pixels; cached to CSV |
| `correlation_scores(time_split)` | `{(i,j): Pearson r}` for all valid pixels; cached to CSV |
| `values_per_pixel(time_split)` | `{(i,j): count}` of valid observations; cached to CSV |
| `spatial_aggregation()` | Compute 3×3 neighbourhood medians for all pixels and timesteps; cached to CSV or NetCDF depending on `aggregation_format` (see below) |
| `set_aggregation_format(fmt)` *(classmethod)* | Set the on-disk cache format used by `spatial_aggregation()` for all instances: `"csv"` or `"netcdf"` |
| `lake_bloom_kde(ax, qa_value, start_year, end_year, plt_kwargs, probability, interval, resolution, x_max, y_max)` | Lake-wide 2D KDE of bloom timing (green-up advance DOY vs. green-down onset DOY), or — with `probability=True` — a bloom-window probability contour. See [Bloom-timing KDE](#bloom-timing-kde) |
| `calculate_bloom_probabilities_from_kde(qa_value, start_year, end_year, interval, x_max, y_max, resolution, save_path)` | Probability of an `interval`-day bloom window under the fitted KDE, evaluated on a `resolution`-day grid; returns a DataFrame and optionally caches the grid to NetCDF via `save_path` |

The metric methods (`r2_scores`, `MAD_scores`, `RMSE_scores`, `correlation_scores`, `values_per_pixel`) each accept a `time_split` argument: a single-element list containing one `[start, end]` year pair, e.g. `[[2003, 2012]]` or `[[0, 9999]]` for the full series. Results are cached to CSV; the cache filename encodes the time window so different windows are stored independently.

All metric and aggregation computations are cached to CSV on first call and loaded from cache on subsequent calls. Interactive plot methods require an interactive Matplotlib backend (`%matplotlib widget`).

### Bloom-timing KDE

`lake_bloom_kde` fits a 2D Gaussian KDE (via `scipy.stats.gaussian_kde`) to bloom events pooled across every valid pixel in the lake, pairing each detected peak with its bracketing green-up advanced (start) and green-down onset (end) events. Building this event set is expensive (`assemble_kde_data` + `prep_kde_data` scan every pixel's phenology arrays in parallel), so the bracketed events are cached once per lake/variable/version to `kde_events.csv` (see [Output Format](#output-format) below) and reused across calls with different `qa_value` or year-range filters.

Two plotting modes, selected by `probability`:
- **Density (default, `probability=False`)** — a filled contour of the raw KDE density over green-up advance DOY (x) vs. green-down onset DOY (y).
- **Probability (`probability=True`)** — calls `calculate_bloom_probabilities_from_kde` to compute, via a 2D summed-area table over the KDE density, the probability that a bloom falls within every `interval`-day window (default 21 days, `resolution`-day grid steps). Windows overlap, so probabilities do not sum to 1 — each is `P(bloom in *this* window)`. A red contour line marks the 0.05 probability threshold.

`qa_value` filters which peak QA levels are included (e.g. `{0}` for Good only); `start_year`/`end_year` filter by the bloom's overlap with that range. When run through the pipeline, the `timing_plots` stage produces one KDE plot per entry in `kde_qa` (see [Configuration Files](#configuration-files)) for every configured `time_splits` window.

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

Extraction and phenology NetCDFs live under `out_folder` (`v{version}/`). Everything `PhenologyVisualization` computes or plots — cached metric/aggregation/KDE data *and* all analysis plots — is written under a separate **`lake_analysis/{lake_str}/`** folder, one per lake, keyed off `out_folder`'s parent directory (`lake_str` = `"ID{lake_id}_{lake_name}"`).

This applies whenever `out_folder` is reassigned to the lake analysis path, which is exactly what `main.py`'s `analysis` stage does for every lake before calling any scoring/plotting method — so this is the layout you'll see after a normal pipeline run. (If you instantiate `PhenologyVisualization` directly, e.g. in a notebook, without reassigning `out_folder`, cached metrics and spatial-aggregation values fall back to `out_folder/calculated_values/...` instead — but the KDE event cache always resolves to `lake_analysis/{lake_str}/calculated_values/kde_data/`, since `build_kde_path` derives that location independently of `out_folder`.)

> **Naming requirement:** `out_folder` itself must be a directory literally named `v{version}` (e.g. `v2.1`, `v3.1`), and it must sit alongside any other version folders under a shared parent directory. `PhenologyVisualization` derives `self.version` by reading the basename of `out_folder` from the `phenology_path` it is given (stripping the leading `v`), and `main.py`'s `comparison` stage locates other versions as siblings of `out_folder` (`os.path.join(os.path.dirname(out_folder), other_version)`). Renaming or nesting `out_folder` differently breaks both the version label shown in plots and the comparison stage.

```
{parent}/                                # shared parent of all version folders
├── {data_folder}/                       # typically 'merged_product'
│   ├── v3.0/                            # = out_folder for a run
│   │   ├── extract/
│   │   │   └── {variable}/
│   │   │       └── {lake_id}.nc
│   │   └── phenology/
│   │       └── {variable}/
│   │           ├── {lake_id}.nc
│   │           └── checkpoints/
│   │               └── {lake_id}/
│   │                   └── bs{batch_size}/
│   │                       └── *.npy          # temporary; deleted after successful write
│   └── v2.1/                            # = out_folder for a v2.1 run; same layout as v3.0/
│
└── lake_analysis/                # = dirname(dirname(out_folder))/lake_analysis
    └── {lake_str}/                      #   lake_str = "ID{lake_id}_{lake_name}"
        ├── calculated_values/
        │   ├── metrics/
        │   │   └── {metric_name}/             # r2 | MAD | RMSE | correlation | values_per_pixel
        │   │       └── v{version}/
        │   │           └── {variable}/
        │   │               ├── full_ts.csv                       # time_split [0, 9999]
        │   │               ├── ts_2002_to_{end}.csv               # time_split [0, end]
        │   │               ├── ts_{start}_to_2024.csv             # time_split [start, 9999]
        │   │               └── ts_{start}_to_{end}.csv            # time_split [start, end]
        │   ├── spatial_aggregation_values/
        │   │   └── v{version}/
        │   │       └── {variable}/
        │   │           └── aggregation_background_values.{csv|nc} # format set by aggregation_format
        │   └── kde_data/
        │       └── v{version}/
        │           └── {variable}/
        │               └── kde_events.csv        # bracketed peak events cached by build_kde_path
        └── plots/
            ├── metric_maps/
            │   ├── {variable}_v{version}_{metric}_full_ts.png      # single time_split [0, 9999]
            │   ├── {variable}_v{version}_{metric}_split_ts.png     # two time_splits
            │   └── {variable}_v{version}_{metric}_{n}_split_ts.png # n > 2 time_splits
            ├── pixel_plots/
            │   ├── aggregated/
            │   │   └── {i}_{j}/
            │   │       ├── location.png
            │   │       ├── {variable}_v{version}_full_ts.png
            │   │       ├── {variable}_v{version}_split_ts.png
            │   │       ├── {variable}_v{version}_peaks_*.png
            │   │       ├── {variable}_v{version}_{peaks|troughs}_qa_boxplot_*.png
            │   │       └── {variable}_v{version}_heatmap_{i}_{j}.png
            │   ├── not_aggregated/
            │   │   └── {i}_{j}/         # same structure as aggregated/
            │   ├── comparisons/                 # only when comparison: true
            │   │   └── {i}_{j}/
            │   │       └── comparison_pks_{pair_label}_{ts_suffix}[_agg].png
            │   └── summaries/
            │       └── {i}_{j}/
            │           └── summary_{variable}.txt
            └── timing_plots/                    # only when timing_plots: true
                ├── {variable}_v{version}_doy_peaks_*.png
                ├── {variable}_v{version}_doy_green_up_mid_*.png
                ├── {variable}_v{version}_lake_heatmap.png
                ├── {variable}_v{version}_kde_qa{qa_suffix}_*.png
                └── {variable}_v{version}_lake_{peaks|troughs}_qa_boxplot_*.png
```

## Fault Tolerance and Restart

The phenology stage saves intermediate results as `.npy` checkpoint files in `data/{version}/phenology/{variable}/checkpoints/{lake_id}/bs{batch_size}/`. If a run is interrupted, re-running the same command will skip completed batches and resume from where it left off. Checkpoints are deleted after the final NetCDF is successfully written.

Both stages skip lakes whose output file already exists, so it is safe to rerun the command after adding new lakes to the `lakes` list.

## Provenance

Each invocation of `main.py` writes one file, `{out_folder}/provenance_logs/provenance_{run_id}.json` (`run_id` is the UTC timestamp of that invocation), containing a record for every pipeline stage that actually executes during that run (`extract`, `phenology`, `analysis`). Stages that are skipped (`"extract": false`, etc.) do not add an entry, so a run like an analysis-only rerun against existing extract/phenology data still leaves a record of when it ran and with what settings. Requires `"provenance": true`.

Each entry contains:

| Field | Description |
|-------|-------------|
| `stage` | `"extract"`, `"phenology"`, or `"analysis"` |
| `timestamp` | UTC timestamp when the stage was launched |
| `git_commit` | Short hash of the checked-out commit, or `null` if not in a git repo |
| `args_file` | Name of the `args/*.json` file the run was launched from |
| `args` | The full resolved parameters dict (after `parse_args` defaults are applied) |
| `lakes` | IDs of the lakes processed in this stage |
| `threads`, `parallel`, `batch_size` | Parallelisation settings (extract/phenology only) |

Within one invocation's file, stage entries are appended, never overwritten — so if a single `main.py` call runs `extract` → `phenology` → `analysis` (or resumes an interrupted run), all of that call's stages land in the same file rather than erasing each other. Each separate `main.py` invocation gets its own timestamped file, so a full history across many runs means listing `provenance_logs/`, not reading one ever-growing file:

```json
{
  "out_folder": "v3.1",
  "runs": [
    { "stage": "extract",   "timestamp": "...", "git_commit": "d979727", "args_file": "v3_chla", "args": {...}, "lakes": [5, 15, 6], "threads": 50, "parallel": "pixels" },
    { "stage": "phenology", "timestamp": "...", "git_commit": "d979727", "args_file": "v3_chla", "args": {...}, "lakes": [5, 15, 6], "threads": 50, "parallel": "pixels", "batch_size": 100 },
    { "stage": "analysis",  "timestamp": "...", "git_commit": "facdadb", "args_file": "v3_chla_analysis", "args": {...}, "lakes": [12] }
  ]
}
```
