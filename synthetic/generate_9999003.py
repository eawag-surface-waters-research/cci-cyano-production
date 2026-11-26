"""Generate scenario 9999003: annual sinusoid with two mean shifts."""

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

    lake = load_scenario_metadata(9999003, p)
    dates, lat, lon, geometry = create_grid_and_time(p)
    baseline = create_sinusoid(dates, p)

    n = len(dates)
    first_change, second_change = create_daily_sections(n)

    offsets = np.zeros(n, dtype=float)
    offsets[first_change:second_change] = -2.0
    offsets[second_change:] = 2.0

    reference = baseline + offsets

    lake["offsets_ug_l"] = [0.0, -2.0, 2.0]
    lake["offset_reference"] = "original baseline"
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
    print("Mean shifts:")
    print(f"  Before {dates[first_change].date()}: offset 0 ug/L")
    print(f"  From   {dates[first_change].date()}: offset -2 ug/L")
    print(f"  From   {dates[second_change].date()}: offset +2 ug/L")


if __name__ == "__main__":
    main()