"""Generate scenario 9999003: annual sinusoid with two mean shifts."""

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

    # Isolate the mean shifts: no additional noise in this scenario.
p["noise"]["std_ug_l"] = 0.0

# Load the scenario metadata from the authoritative Excel source.
lake = load_scenario_metadata(9999003, p)

# Create the daily time axis, pixel coordinates, and lake boundary.
dates, lat, lon, geometry = create_grid_and_time(p)

# Generate one noise-free baseline concentration value for each date.
baseline = create_sinusoid(dates, p)

# Split the daily time series into three nearly equal sections.
n = len(dates)
first_change, second_change = create_daily_sections(n)

# Define the mean shifts for the three sections.
offsets = np.zeros(n, dtype=float)
offsets[first_change:second_change] = -2.0
offsets[second_change:] = 2.0

# Apply the offsets to the original baseline.
reference = baseline + offsets

# Add the scenario-specific provenance fields to the Excel metadata.
lake["offsets_ug_l"] = [0.0, -2.0, 2.0]
lake["offset_reference"] = "original baseline"
lake["split_method"] = "daily observation index"
lake["change_indices"] = [first_change, second_change]
lake["change_dates"] = [
    dates[first_change].isoformat(),
    dates[second_change].isoformat(),
]

# Broadcast the shifted reference across the full pixel grid.
values = broadcast_reference_to_grid(reference, lat, lon)

if np.any(values < 0):
        raise ValueError(
            "Negative concentrations generated. "
            "Review the baseline and offsets; values have not been clipped."
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

print("\nMean shifts:")
print(f"  Before {dates[first_change].date()}: offset 0 ug/L")
print(f"  From   {dates[first_change].date()}: offset -2 ug/L")
print(f"  From   {dates[second_change].date()}: offset +2 ug/L")


if __name__ == "__main__":
    main()