"""Generate scenario 9999004: annual sinusoid with amplitude changes."""

import numpy as np

from common import (
    load_defaults,
    load_scenario_metadata,
    create_grid_and_time,
    create_sinusoid,
    create_daily_sections,
    broadcast_reference_to_grid,
    add_gaussian_noise,
    save_results,
)


def main():
    p = load_defaults()

    lake = load_scenario_metadata(9999004, p)
    dates, lat, lon, geometry = create_grid_and_time(p)
    baseline = create_sinusoid(dates, p)

    mean = float(p["baseline"]["mean_ug_l"])
    base_amplitude = float(p["baseline"]["amplitude_ug_l"])

    n = len(dates)
    first_change, second_change = create_daily_sections(n)

    factors = np.ones(n, dtype=float)
    factors[first_change:second_change] = 0.5
    factors[second_change:] = 1.5

    seasonal_signal = baseline - mean
    reference = mean + factors * seasonal_signal

    lake["amplitude_factors"] = [1.0, 0.5, 1.5]
    lake["scaling_reference"] = "deviation from baseline mean"
    lake["split_method"] = "daily observation index"
    lake["change_indices"] = [first_change, second_change]
    lake["change_dates"] = [
        dates[first_change].isoformat(),
        dates[second_change].isoformat(),
    ]

    values = broadcast_reference_to_grid(reference, lat, lon)
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

    print(f"\nNoise type: {noise_type}")
    print(f"Unchanged baseline mean: {mean:g} ug/L")
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