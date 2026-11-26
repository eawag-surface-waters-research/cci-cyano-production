"""Generate scenario 9999013: two seasonal blooms per year."""

import numpy as np

from common import (
    load_defaults,
    load_scenario_metadata,
    create_grid_and_time,
    add_gaussian_noise,
    save_results,
)


def main():
    p = load_defaults()
    lake = load_scenario_metadata(9999013, p)

    dates, lat, lon, geometry = create_grid_and_time(p)

    baseline = p["baseline"]
    mean = float(baseline["mean_ug_l"])
    amplitude = float(baseline["amplitude_ug_l"])
    phase_fraction = float(baseline["phase_months"]) / 12.0

    year_days = np.where(dates.is_leap_year, 366.0, 365.0)
    year_fraction = (dates.dayofyear.to_numpy() - 1) / year_days

    reference = mean + amplitude * np.sin(
        2 * 2 * np.pi * (year_fraction - phase_fraction)
    )

    values = np.broadcast_to(
        reference[:, None, None],
        (len(dates), len(lat), len(lon)),
    ).copy()

    values, noise_type = add_gaussian_noise(
        values,
        dates,
        p,
        std_ug_l=p["noise"]["std_ug_l"],
    )

    if np.any(values < 0):
        raise ValueError(
            "Negative concentrations generated. Review mean, amplitude "
            "and noise settings; values have not been clipped."
        )

    print(f"Scenario: {lake['id']} - {lake['name']}")
    print(f"Noise type: {noise_type}")
    print(f"Blooms per year: 2")
    print(f"Range: {reference.min():.3f} to {reference.max():.3f} ug/L")

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


if __name__ == "__main__":
    main()