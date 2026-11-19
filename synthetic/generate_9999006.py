"""Generate scenario 9999006: bloom timing shifts south to north."""

import numpy as np

from common import (
    load_defaults,
    load_scenario_metadata,
    create_grid_and_time,
    create_lat_phase_gradient,
    add_gaussian_noise,
    save_results,
)


def main():
    p = load_defaults()
    lake = load_scenario_metadata(9999006, p)

    dates, lat, lon, geometry = create_grid_and_time(p)

    reference = create_lat_phase_gradient(dates, lat, lon, p)

    expected_shape = (len(dates), len(lat), len(lon))
    if reference.shape != expected_shape:
        raise ValueError(
            f"Expected reference shape {expected_shape}, "
            f"got {reference.shape}."
        )

    values, noise_type = add_gaussian_noise(
        reference.copy(),
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

    phase_gradient = float(
        p["spatial"]["phase_shift_months_per_pixel_lat"]
    )
    total_shift = (len(lat) - 1) * phase_gradient
    print(f"South-to-north phase shift: {total_shift:.3f} months")

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