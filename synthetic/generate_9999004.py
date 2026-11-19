"""Generate scenario 9999004: annual sinusoid with amplitude changes."""

import numpy as np

from common import (
    load_defaults,
    load_scenario_metadata,
    create_grid_and_time,
    create_sinusoid,
    create_daily_sections,
    broadcast_reference_to_grid,
    save_results,
)

def main():
    p = load_defaults()

    # No additional noise: isolate the amplitude changes.
p["noise"]["std_ug_l"] = 0.0

# Load the scenario metadata from the authoritative Excel source.
lake = load_scenario_metadata(9999004, p)

# Create the daily time axis, pixel coordinates, and lake boundary.
dates, lat, lon, geometry = create_grid_and_time(p)

# Generate one noise-free baseline concentration value for each date.
baseline = create_sinusoid(dates, p)

mean = float(p["baseline"]["mean_ug_l"])
base_amplitude = float(p["baseline"]["amplitude_ug_l"])

# Split the daily time series into three nearly equal sections.
n = len(dates)
first_change, second_change = create_daily_sections(n)

# Define the amplitude factors for the three sections.
factors = np.ones(n, dtype=float)
factors[first_change:second_change] = 0.5
factors[second_change:] = 1.5

# Scale the oscillation around the mean, not the mean itself.
seasonal_signal = baseline - mean
reference = mean + factors * seasonal_signal

# Add the scenario-specific provenance fields to the Excel metadata.
lake["amplitude_factors"] = [1.0, 0.5, 1.5]
lake["scaling_reference"] = "deviation from baseline mean"
lake["split_method"] = "daily observation index"
lake["change_indices"] = [first_change, second_change]
lake["change_dates"] = [
    dates[first_change].isoformat(),
    dates[second_change].isoformat(),
]

# Broadcast the amplitude-modified reference across the full pixel grid.
values = broadcast_reference_to_grid(reference, lat, lon)

if np.any(values < 0):
        raise ValueError(
            "Negative concentrations generated. Review the baseline "
            "and amplitude factors; values have not been clipped."
        )

save_results(
    lake=lake,
    dates=dates,
    lat=lat,
    lon=lon,
    geometry=geometry,
    values=values,
    reference=reference,
    p=p,
)

print(f"\nUnchanged baseline mean: {mean:g} ug/L")
print("Amplitude changes:")

for label, factor in [
    (f"Before {dates[first_change].date()}", 1.0),
    (f"From   {dates[first_change].date()}", 0.5),
    (f"From   {dates[second_change].date()}", 1.5),
]:
    amplitude = base_amplitude * factor
    print(
        f"  {label}: factor {factor:g}, "
        f"amplitude {amplitude:g} ug/L, "
        f"theoretical range {mean - amplitude:g}"
        f" to {mean + amplitude:g} ug/L"
    )

if __name__ == "__main__":
    main()