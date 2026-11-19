"""Generate scenario 9999005: bloom timing shifts west to east."""

import numpy as np

from common import (
    load_defaults,
    create_grid_and_time,
    create_lon_phase_gradient,
    add_gaussian_noise,
    save_results,
)


def main():
    p = load_defaults()

    lake = {
        "id": 9999005,
        "short_name": "SYNTH9999005",
        "name": "SYNTH_lon_shift_grad_1",
        "tag": "lon_shift_grad_1",
        "description": (
            "1 seasonal bloom throughout time and space, bloom shifted "
            "later in time linearly across longitude"
        ),
    }

    dates, lat, lon, geometry = create_grid_and_time(p)

    reference = create_lon_phase_gradient(dates, lat, lon, p)

    values = reference.copy()

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


if __name__ == "__main__":
    main()