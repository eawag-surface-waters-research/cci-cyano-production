# Synthetic lake scenarios

Synthetic chlorophyll-a scenarios for testing the phenology and postprocessing
pipeline. Each scenario has a known, noise-free reference signal. Gaussian
noise is added to the generated data, and the pipeline results are compared
with the reference.

## Status

- Implemented: `9999001` to `9999008` and `9999013`.
- Planned: `9999009` to `9999012` and `9999014` to `9999016`.

## Files

- `common.py`: shared functions (grid, time axis, signals, noise, output).
- `default.json`: shared settings (grid, time range, mean, amplitude, noise).
- `generate_<id>.py`: one generator per scenario.
- `check_reference_events.py`: compares fitted peaks and troughs with the
  reference.
- `data/`, `lake_analysis/`, `validation/`, `backups/`: generated or local
  files. They are not committed to Git.

## Scenarios

| id | name | description |
|---|---|---|
| 9999001 | SYNTH_full_constant_1 | 1 seasonal bloom, constant in time and space |
| 9999002 | SYNTH_inc_variance_1 | as 001, noise increases linearly over time |
| 9999003 | SYNTH_snsr_bias_1 | mean shifts at 1/3 and 2/3 of the time series |
| 9999004 | SYNTH_snsr_amp_1 | amplitude changes at 1/3 and 2/3 of the time series |
| 9999005 | SYNTH_lon_shift_grad_1 | bloom shifts later across longitude |
| 9999006 | SYNTH_lat_shift_grad_1 | bloom shifts later across latitude |
| 9999007 | SYNTH_lon_amp_grad_1 | amplitude increases across longitude |
| 9999008 | SYNTH_lat_amp_grad_1 | amplitude increases across latitude |
| 9999013 | SYNTH_full_constant_2 | 2 seasonal blooms, constant in time and space |

The full catalogue of 16 scenarios, including the planned ones, is defined in
the metadata sheet `synth_descriptions` of `synth_lakes.xlsx`.

## Noise

- The reference stays noise-free. Independent Gaussian noise is added to the
  generated data.
- Scenarios `001`, `003` to `008` and `013`: constant standard deviation
  (`noise.std_ug_l` in `default.json`).
- Scenario `002`: standard deviation increases from 0.3 to 0.7 ug/L.
- The random seed comes from `default.json` and is shared by all scenarios.
- Negative values stop the generator and are never clipped. Choose the mean
  high enough that mean minus maximum amplitude stays above the noise range.

## Workflow

1. Generate the input data (`<id>` is e.g. `9999001`):

   ```powershell
   python .\synthetic\generate_<id>.py
   ```

2. Run phenology and postprocessing with `args/synthetic_fit.json`. First run
   with `phenology` true and `postprocess` false, then the reverse:

   ```powershell
   python .\scripts\main.py -f synthetic_fit
   ```

3. Create plots (set `plot` and `timingplots` to true; the file name has no
   `.json` ending):

   ```powershell
   python .\scripts\plot.py --file synthetic_fit
   ```

4. Check the fitted events against the reference:

   ```powershell
   python .\synthetic\check_reference_events.py 9999001 9999002 9999003 9999004 9999005 9999006 9999007 9999008 9999013
   ```

After changing generator settings, delete the old outputs
(`data/v3.0/phenology`, `data/v3.0/calculated_values`, `lake_analysis`,
`validation`) before re-running.

## Outputs

- Input data: `data/v3.0/extract/chla/<id>.nc`
- Reference: `data/v3.0/reference/<id>/` with `reference.csv` (time-series
  scenarios `001` to `004`, `013`) or `reference.nc` (spatial scenarios `005`
  to `008`), plus `resolved_params.json` and `preview.png`.
- Phenology results: `data/v3.0/phenology/chla/`
- Postprocessing results: `data/v3.0/calculated_values/`
- Plots and summaries: `lake_analysis/`
- Validation tables: `validation/<id>/`

## Validation results (with noise)

The checker excludes 60 days at each boundary and matches events one-to-one
within 60 days. This matching distance is not an acceptance threshold.
All scenarios had 96 of 96 valid pixels.

| id | matched peaks / troughs | missing (pks / trgs) | extra (pks / trgs) | time MAE in days (pks / trgs) | value RMSE in ug/L (pks / trgs) |
|---|---|---|---|---|---|
| 9999001 | 2016 / 1920 | 0 / 0 | 0 / 0 | 6.98 / 6.60 | 0.283 / 0.280 |
| 9999002 | 2016 / 1920 | 0 / 0 | 0 / 0 | 3.84 / 3.70 | 0.098 / 0.096 |
| 9999003 | 2016 / 1920 | 0 / 0 | 0 / 0 | 7.06 / 8.11 | 0.285 / 0.313 |
| 9999004 | 2016 / 1920 | 96 / 96 | 4 / 4 | 6.91 / 8.13 | 0.269 / 0.370 |
| 9999005 | 2016 / 1920 | 0 / 0 | 0 / 0 | 6.49 / 6.41 | 0.284 / 0.276 |
| 9999006 | 2016 / 1920 | 0 / 0 | 0 / 0 | 6.72 / 6.58 | 0.283 / 0.278 |
| 9999007 | 2016 / 1920 | 0 / 0 | 0 / 0 | 6.35 / 6.03 | 0.294 / 0.290 |
| 9999008 | 2016 / 1920 | 0 / 0 | 0 / 0 | 6.60 / 6.28 | 0.291 / 0.288 |
| 9999013 | 3936 / 3936 | 0 / 0 | 0 / 0 | 3.40 / 3.44 | 0.355 / 0.350 |

Scenario `9999004` has missing and additional events. The cause has not been
analysed yet.

Gradient checks (event means per grid column or row):

- `9999005`: mean peak day of year rises from 183.2 to 200.2 along longitude.
- `9999006`: mean peak day of year rises from 183.2 to 194.2 along latitude.
- `9999007`: peak values rise from 21.0 to 23.2 ug/L along longitude.
- `9999008`: peak values rise from 21.0 to 22.4 ug/L along latitude.

The fit follows the reference in all four cases (within 1.6 days or 0.1 ug/L).

## Notes

- Peak day-of-year maps only show peaks between day 160 and 250 (default of
  the plotting function). Green-up maps are therefore empty, and single pixels
  can be blank in the peak maps.
- Pending: full-curve validation, analytical generator timing,
  green-up/green-down validation and boundary-event analysis.