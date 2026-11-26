from pathlib import Path
import argparse

import netCDF4
import numpy as np
import pandas as pd
from scipy.signal import find_peaks


ROOT = Path(__file__).resolve().parent
DAY = 86400.0
EDGE_DAYS = 60
MATCH_DAYS = 60


def read_array(variable):
    return np.ma.asarray(variable[:], dtype=float).filled(np.nan)


def match_events(reference_times, fitted_times):
    candidates = []
    for r, rt in enumerate(reference_times):
        for f, ft in enumerate(fitted_times):
            distance = abs(ft - rt) / DAY
            if distance <= MATCH_DAYS:
                candidates.append((distance, r, f))

    used_r, used_f, pairs = set(), set(), []
    for _, r, f in sorted(candidates):
        if r not in used_r and f not in used_f:
            pairs.append((r, f))
            used_r.add(r)
            used_f.add(f)

    return pairs, used_r, used_f


def check(lake_id):
    ref_dir = ROOT / f"data/v3.0/reference/{lake_id}"
    ref_nc_path = ref_dir / "reference.nc"
    ref_csv_path = ref_dir / "reference.csv"
    fit_path = ROOT / f"data/v3.0/phenology/chla/{lake_id}.nc"
    output = ROOT / f"validation/{lake_id}"
    output.mkdir(parents=True, exist_ok=True)

    rows, counts = [], []

    with netCDF4.Dataset(fit_path) as fit:
        if ref_nc_path.exists():
            ref = netCDF4.Dataset(ref_nc_path)
            ref_is_nc = True
        elif ref_csv_path.exists():
            ref = pd.read_csv(ref_csv_path)
            ref_is_nc = False
        else:
            raise FileNotFoundError(
                f"Neither reference.nc nor reference.csv found in {ref_dir}"
            )

        if ref_is_nc:
            for coord in ("lat", "lon"):
                a = read_array(ref.variables[coord])
                b = read_array(fit.variables[coord])
                if a.shape != b.shape or not np.allclose(
                    a, b, rtol=0, atol=1e-5
                ):
                    raise ValueError(f"Koordinaten stimmen nicht: {coord}")

            tv = ref.variables["time"]
            dates = netCDF4.num2date(
                tv[:], tv.units,
                calendar=getattr(tv, "calendar", "standard")
            )
            time = np.asarray(netCDF4.date2num(
                dates,
                "seconds since 1970-01-01 00:00:00",
                calendar="standard"
            ), dtype=float)

            n_time = len(time)
            n_lat = fit.dimensions["lat"].size
            n_lon = fit.dimensions["lon"].size
            ref_values = read_array(ref.variables["reference_chla_ug_l"])

            if ref_values.ndim == 1:
                if ref_values.shape != (n_time,):
                    raise ValueError(
                        "1D reference does not match time dimension."
                    )
                ref_values = np.broadcast_to(
                    ref_values[:, None, None],
                    (n_time, n_lat, n_lon),
                ).copy()
        else:
            dates = pd.to_datetime(ref["date"], utc=True)
            time = np.asarray(
                [d.timestamp() for d in dates], dtype=float
            )
            n_time = len(time)
            n_lat = fit.dimensions["lat"].size
            n_lon = fit.dimensions["lon"].size
            series = ref["reference_chla_ug_l"].to_numpy(dtype=float)

            if series.shape != (n_time,):
                raise ValueError(
                    "CSV reference does not match time dimension."
                )
            ref_values = np.broadcast_to(
                series[:, None, None],
                (n_time, n_lat, n_lon),
            ).copy()

        if not np.all(np.diff(time) > 0):
            raise ValueError("Referenzzeit muss streng ansteigend sein.")

        step_days = np.median(np.diff(time)) / DAY
        lower = time[0] + EDGE_DAYS * DAY
        upper = time[-1] - EDGE_DAYS * DAY
        valid = np.isfinite(
            read_array(fit.variables["smoothing_parameter"])
        )

        print(f"\nSzenario: {lake_id}")
        print(f"Referenzformat: {'NetCDF' if ref_is_nc else 'CSV'}")
        print(f"Gültige Fit-Pixel: {valid.sum()}/{valid.size}")
        print(f"Medianes Referenzintervall: {step_days:.6f} Tage")
        print(f"Randausschluss: {EDGE_DAYS} Tage je Seite")
        print(f"Maximaler Zuordnungsabstand: {MATCH_DAYS} Tage")

        for i, j in np.argwhere(valid):
            signal = ref_values[:, i, j]

            if not np.all(np.isfinite(signal)):
                raise ValueError(f"Referenzlücken bei Pixel {(i, j)}")

            for kind, sign in (("pks", 1), ("trgs", -1)):
                indices, _ = find_peaks(sign * signal)
                reference_t = time[indices]
                reference_y = signal[indices]
                keep = (reference_t >= lower) & (reference_t <= upper)
                reference_t = reference_t[keep]
                reference_y = reference_y[keep]

                fitted_t = np.ma.asarray(
                    fit.variables[f"{kind}_x"][i, j, :],
                    dtype=float
                ).filled(np.nan)
                fitted_y = np.ma.asarray(
                    fit.variables[f"{kind}_y"][i, j, :],
                    dtype=float
                ).filled(np.nan)

                finite_t = np.isfinite(fitted_t)
                inside = finite_t & (fitted_t >= lower) & (fitted_t <= upper)
                boundary = int(np.sum(finite_t & ~inside))
                invalid_y = int(np.sum(inside & ~np.isfinite(fitted_y)))
                keep = inside & np.isfinite(fitted_y)
                fitted_t, fitted_y = fitted_t[keep], fitted_y[keep]

                pairs, used_r, used_f = match_events(
                    reference_t, fitted_t
                )

                counts.append({
                    "lake_id": lake_id, "i": i, "j": j,
                    "event": kind,
                    "reference_count": len(reference_t),
                    "fit_count": len(fitted_t),
                    "matched": len(pairs),
                    "missing": len(reference_t) - len(used_r),
                    "extra": len(fitted_t) - len(used_f),
                    "fit_boundary_excluded": boundary,
                    "fit_invalid_value": invalid_y,
                })

                for r, f in pairs:
                    rows.append({
                        "lake_id": lake_id, "i": i, "j": j,
                        "event": kind,
                        "reference_unix": reference_t[r],
                        "fit_unix": fitted_t[f],
                        "time_error_days":
                            (fitted_t[f] - reference_t[r]) / DAY,
                        "reference_ug_l": reference_y[r],
                        "fit_ug_l": fitted_y[f],
                        "value_error_ug_l":
                            fitted_y[f] - reference_y[r],
                    })

        if ref_is_nc:
            ref.close()

    count_df = pd.DataFrame(counts)
    event_df = pd.DataFrame(rows)
    count_df.to_csv(output / "event_counts.csv", index=False)
    event_df.to_csv(output / "event_errors.csv", index=False)

    print("\nEreigniszahlen im Prüfzeitraum:")
    columns = [
        "reference_count", "fit_count", "matched", "missing", "extra",
        "fit_boundary_excluded", "fit_invalid_value"
    ]
    if not count_df.empty:
        print(count_df.groupby("event")[columns].sum().to_string())

    if event_df.empty:
        print("\nKeine zugeordneten Ereignisse.")
        return

    summary = []
    for kind, group in event_df.groupby("event"):
        dt = group["time_error_days"].to_numpy()
        dy = group["value_error_ug_l"].to_numpy()
        summary.append({
            "event": kind,
            "matched": len(group),
            "time_bias_days": dt.mean(),
            "time_mae_days": np.abs(dt).mean(),
            "time_max_abs_days": np.abs(dt).max(),
            "value_bias_ug_l": dy.mean(),
            "value_rmse_ug_l": np.sqrt(np.mean(dy ** 2)),
            "value_max_abs_ug_l": np.abs(dy).max(),
        })

    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(output / "event_summary.csv", index=False)
    print("\nFehler der zugeordneten Ereignisse:")
    print(summary_df.to_string(index=False, float_format=lambda x: f"{x:.6f}"))
    print(f"\nPrüftabellen: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("lake_ids", nargs="+", type=int)
    args = parser.parse_args()
    for lake_id in args.lake_ids:
        check(lake_id)