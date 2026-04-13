# CCI Cyano Production

![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)
![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)

A processing pipeline for extracting per-lake time series and computing phenology metrics from ESA CCI Lakes satellite data. The pipeline downscales global NetCDF imagery to individual lake polygons and fits cubic splines to derive green-up, green-down, peak, and trough dates at the pixel level.

## Overview

The pipeline runs in two sequential stages:

1. **Extract** — clips global satellite NetCDF files to each lake's bounding box and raster mask, producing a per-lake time series NetCDF containing the target variable and quality flags.
2. **Phenology** — reads the extracted time series, fits a smoothing cubic spline per pixel, and extracts phenology metrics (peaks, troughs, green-up/green-down onset, mid, and advanced dates, and data gaps).

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

## Fault Tolerance and Restart

The phenology stage saves intermediate results as `.npy` checkpoint files in `data/{version}/phenology/{variable}/checkpoints/{lake_id}/bs{batch_size}/`. If a run is interrupted, re-running the same command will skip completed batches and resume from where it left off. Checkpoints are deleted after the final NetCDF is successfully written.

Both stages skip lakes whose output file already exists, so it is safe to rerun the command after adding new lakes to the `lakes` list.
