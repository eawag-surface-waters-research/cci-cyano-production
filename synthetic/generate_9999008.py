"""Generate scenario 9999008: amplitude increases south to north."""

import numpy as np

from common import (
    load_defaults,
    load_scenario_metadata,
    create_grid_and_time,
    create_lat_amplitude_gradient,
    add_gaussian_noise,
    save_results,
)


def main():
    p = load_defaults()
    lake = load_scenario_metadata(9999008, p)

    dates, lat, lon, geometry = create_grid_and_time(p)

    reference = create_lat_amplitude_gradient(dates, lat, lon, p)

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

    base_amplitude = float(p["baseline"]["amplitude_ug_l"])
    amplitude_gradient = float(
        p["spatial"]["amplitude_increase_per_pixel_lat"]
    )
    north_amplitude = (
        base_amplitude + (len(lat) - 1) * amplitude_gradient
    )

    print(f"Scenario: {lake['id']} - {lake['name']}")
    print(f"Noise type: {noise_type}")
    print(f"South amplitude: {base_amplitude:.3f} ug/L")
    print(f"North amplitude: {north_amplitude:.3f} ug/L")

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