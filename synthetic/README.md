## Scenarios

The catalogue contains 16 synthetic lake scenarios.
Scenarios `9999001` through `9999008` have been implemented and
processed with the phenology and postprocessing pipeline.
Scenarios `9999009` through `9999016` are planned.

The descriptions below define the intended scenarios. Parameters are
defined in the generator scripts and `default.json`.

All paths below are relative to the `synthetic/` directory:

- Generated extracts: `data/v3.0/extract/chla/`
- Synthetic references: `data/v3.0/reference/<lake_id>/`
- Spatial references for scenarios `9999005` through `9999008`:
  `data/v3.0/reference/<lake_id>/reference.nc`
- Generated plots and summaries: `lake_analysis/`
- Reference validation tables: `validation/<lake_id>/`

Generated data, plots and validation tables are excluded from Git.

## Reference event validation

Run from the repository root:

```powershell
python synthetic/check_reference_events.py 9999005 9999006 9999007 9999008
```

The checker compares stored peak and trough events with extrema of the
daily sampled, noise-free reference. It excludes 60 days at each
boundary and matches events one-to-one within a maximum distance of
60 days. This matching distance is not an acceptance threshold.

For the tested outputs of scenarios `9999005` through `9999008`:

- All 96 pixels per scenario had valid smoothing parameters.
- Each scenario had 2016 matched peaks and 1920 matched troughs
  within the evaluation window.
- No missing or additional extrema were found within that window.
- Timing errors were zero at the reported precision.
- Maximum absolute event-value errors were below 0.000003 ug/L.

Full-curve validation, analytical generator timing, green-up/green-down
validation and boundary-event analysis remain pending.

## Scenario catalogue

| id | short_name | name | tag | Description | Base chl-a distribution |
|---|---|---|---|---|---|
| 9999001 | SYNTH9999001 | SYNTH_full_constant_1 | full_constant_1 | 1 seasonal bloom throughout time and space | annual sinusoid |
| 9999002 | SYNTH9999002 | SYNTH_inc_variance_1 | inc_variance_1 | 1 seasonal bloom throughout time and space, noise increases linearly over time | annual sinusoid |
| 9999003 | SYNTH9999003 | SYNTH_snsr_bias_1 | snsr_bias_1 | 1 seasonal bloom throughout time and space, mean shifts at 1/3 and 2/3 of the time series | annual sinusoid |
| 9999004 | SYNTH9999004 | SYNTH_snsr_amp_1 | snsr_amp_1 | 1 seasonal bloom throughout time and space, amplitude stretched at 1/3 and 2/3 of the time series | annual sinusoid |
| 9999005 | SYNTH9999005 | SYNTH_lon_shift_grad_1 | lon_shift_grad_1 | 1 seasonal bloom throughout time and space, bloom shifted later in time linearly across longitude | annual sinusoid |
| 9999006 | SYNTH9999006 | SYNTH_lat_shift_grad_1 | lat_shift_grad_1 | 1 seasonal bloom throughout time and space, bloom shifted later in time linearly across latitude | annual sinusoid |
| 9999007 | SYNTH9999007 | SYNTH_lon_amp_grad_1 | lon_amp_grad_1 | 1 seasonal bloom throughout time and space, bloom amplitude increases linearly across longitude | annual sinusoid |
| 9999008 | SYNTH9999008 | SYNTH_lat_amp_grad_1 | lat_amp_grad_1 | 1 seasonal bloom throughout time and space, bloom amplitude increases linearly across latitude | annual sinusoid |
| 9999009 | SYNTH9999009 | SYNTH_winter_ice_1 | winter_ice_1 | 1 seasonal bloom, no data for December–May, data approaches and reaches 0 at data gaps | annual Gaussian |
| 9999010 | SYNTH9999010 | SYNTH_dec_ice_1 | dec_ice_1 | 1 seasonal bloom, no-data period decreases from December–May to January–April, data approaches and reaches 0 at data gaps | annual Gaussian |
| 9999011 | SYNTH9999011 | SYNTH_var_ice_1 | var_ice_1 | 1 seasonal bloom, no data for 4 winter months with timing varying year-to-year, data approaches and reaches 0 at data gaps | annual Gaussian |
| 9999012 | SYNTH9999012 | SYNTH_ice_hide_1 | ice_hide_1 | 1 seasonal bloom, no data for December–May, data approaches but does not reach 0 at data gaps | annual Gaussian |
| 9999013 | SYNTH9999013 | SYNTH_full_constant_2 | full_constant_2 | 2 seasonal blooms throughout time and space | semi-annual sinusoid |
| 9999014 | SYNTH9999014 | SYNTH_1_to_2_abrupt | 1_to_2_abrupt | Shifts from 1 to 2 seasonal blooms halfway through the time series, abrupt frequency change | mixed frequency sinusoid |
| 9999015 | SYNTH9999015 | SYNTH_1_to_2_grad | 1_to_2_grad | Shifts from 1 to 2 seasonal blooms linearly through the time series, constant frequency for the first and last 5 years | mixed frequency sinusoid |
| 9999016 | SYNTH9999016 | SYNTH_constant_poly | constant_poly | Multiple non-seasonal blooms throughout time and space | non-seasonal sinusoid |
