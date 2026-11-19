"""Generate scenario 9999001: a spatially uniform annual sinusoid."""

import numpy as np

from common import (
    load_defaults,
    load_scenario_metadata,
    create_grid_and_time,
    create_sinusoid,
    broadcast_reference_to_grid,
    add_gaussian_noise,
    save_results,
)


def main():
    # Load the shared settings from default.json.
    p = load_defaults()

    # Load the scenario metadata from the authoritative Excel source.
    lake = load_scenario_metadata(9999001, p)

    # Create the daily time axis, pixel coordinates, and lake boundary.
    dates, lat, lon, geometry = create_grid_and_time(p)

    # Generate one noise-free concentration value for each date.
    reference = create_sinusoid(dates, p)

    # Broadcast the reference across the full pixel grid.
    values = broadcast_reference_to_grid(reference, lat, lon)

    # Add independent Gaussian noise using the shared noise function.
    values, noise_type = add_gaussian_noise(
        values=values,
        dates=dates,
        p=p,
        std_ug_l=p["noise"]["std_ug_l"],
    )

    # Keep the resolved noise type in the saved parameters.
    p["noise"]["type"] = noise_type

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


if __name__ == "__main__":
    main()