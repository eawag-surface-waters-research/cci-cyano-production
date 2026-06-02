import netCDF4
import argparse
import numpy as np
from csaps import csaps
from matplotlib import pyplot as plt
from functions import unix_to_datetime, unix_to_datenum, datenum_to_datetime, remove_nan

def plot(e_file, p_file):
    with netCDF4.Dataset(e_file) as nc:
        summary = np.array(nc.variables["summary"][:, :])   # (lat=37, lon=95)
        lat = np.array(nc.variables["lat"])
        lon = np.array(nc.variables["lon"])
        t_all = unix_to_datenum(nc.variables["time"])

    valid_coords = np.argwhere(summary > 2)
    selected = list(valid_coords[0])

    fig, (ax_grid, ax_ts) = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle("Click a cell on the left to view its time series", fontsize=11)

    masked_summary = np.ma.masked_where(summary <= 2, summary)
    ax_grid.set_facecolor("#cccccc")
    im = ax_grid.imshow(masked_summary, cmap="viridis", aspect="auto", origin="lower")
    plt.colorbar(im, ax=ax_grid, label="Observation count")
    ax_grid.set_title("Cell grid")
    ax_grid.set_xlabel("Lon index")
    ax_grid.set_ylabel("Lat index")

    sel_marker, = ax_grid.plot([], [], "r*", markersize=14, zorder=5, label="Selected")
    ax_grid.legend(loc="upper right", fontsize=8)

    def update_time_series(x, y):
        """Reload and redraw the right-hand time-series panel for cell (x=lat_idx, y=lon_idx)."""
        ax_ts.clear()
        try:
            with netCDF4.Dataset(p_file) as nc:
                smoothing = float(nc.variables["smoothing_parameter"][x, y])

                pks_x_raw = np.array(nc.variables["pks_x"][x, y, :])
                pk_mask = ~np.isnan(pks_x_raw)
                pks_x = unix_to_datetime(pks_x_raw[pk_mask])
                pks_y = np.array(nc.variables["pks_y"][x, y, :])[pk_mask]
                pks_qa = np.array(nc.variables["pks_qa"][x, y, :])[pk_mask]

                trgs_x_raw = np.array(nc.variables["trgs_x"][x, y, :])
                trg_mask = ~np.isnan(trgs_x_raw)
                trgs_x = unix_to_datetime(trgs_x_raw[trg_mask])
                trgs_y = np.array(nc.variables["trgs_y"][x, y, :])[trg_mask]
                trgs_qa = np.array(nc.variables["trgs_qa"][x, y, :])[trg_mask]

                gap_starts = unix_to_datetime(remove_nan(nc.variables["data_gap_start"][x, y, :]))
                gap_ends = unix_to_datetime(remove_nan(nc.variables["data_gap_end"][x, y, :]))


            with netCDF4.Dataset(e_file) as nc:
                variable = getattr(nc, "variable")
                values = np.array(nc.variables[variable][:, x, y])
                mask = (values != -9999) & (np.array(nc.variables[getattr(nc, 'qa')][:, x, y]) == 0)
                values_m = values[mask]
                time_m = t_all[mask]

            if len(values_m) > 1:
                smooth_x = np.arange(t_all.min(), t_all.max() + 1, 1)
                smooth_y = csaps(time_m, values_m, smooth_x, smooth=smoothing)
                # Shade data gaps
                for gs, ge in zip(gap_starts, gap_ends):
                    ax_ts.axvspan(gs, ge, color="orange", alpha=0.15, zorder=0)
                if len(gap_starts) > 0:
                    ax_ts.axvspan(gap_starts[0], gap_ends[0], color="orange",
                                  alpha=0.15, zorder=0, label="Data gap")

                ax_ts.scatter(datenum_to_datetime(time_m), values_m,
                              color="grey", alpha=0.3, s=10, label="Data")
                ax_ts.plot(datenum_to_datetime(smooth_x), smooth_y,
                           color="blue", linewidth=1, label="Spline")

                # Peaks/troughs colored by QA flag (0=Good, 1=Fair, 2=Poor)
                qa_colors = {0: "green", 1: "orange", 2: "red"}
                qa_labels = {0: "Good", 1: "Fair", 2: "Poor"}
                for qa in (0, 1, 2):
                    pm = pks_qa == qa
                    if pm.any():
                        ax_ts.scatter(pks_x[pm], pks_y[pm], color=qa_colors[qa], s=60,
                                      marker="^", edgecolors="black", linewidths=0.5,
                                      zorder=4, label=f"Peak ({qa_labels[qa]})")
                    tm = trgs_qa == qa
                    if tm.any():
                        ax_ts.scatter(trgs_x[tm], trgs_y[tm], color=qa_colors[qa], s=60,
                                      marker="v", edgecolors="black", linewidths=0.5,
                                      zorder=4, label=f"Trough ({qa_labels[qa]})")

                # Scale y-axis to the peak/trough range (with a small margin)
                feat_y = np.concatenate([pks_y, trgs_y])
                if feat_y.size > 0:
                    lo, hi = feat_y.min(), feat_y.max()
                    pad = (hi - lo) * 0.1 or abs(hi) * 0.1 or 1.0
                    ax_ts.set_ylim(lo - pad, hi + pad)
                ax_ts.legend(fontsize=8)
            else:
                ax_ts.text(0.5, 0.5, "No valid data for this cell",
                           transform=ax_ts.transAxes, ha="center", va="center")

        except Exception as exc:
            ax_ts.text(0.5, 0.5, f"Error: {exc}",
                       transform=ax_ts.transAxes, ha="center", va="center", wrap=True)

        ax_ts.set_title(
            f"lat_idx={x}  lon_idx={y}  |  lat={lat[x]:.3f}°  lon={lon[y]:.3f}°"
        )
        ax_ts.set_xlabel("Date")
        ax_ts.set_ylabel(variable)
        fig.canvas.draw_idle()


    def on_click(event):
        if event.inaxes is not ax_grid or event.xdata is None:
            return
        lon_idx = int(round(event.xdata))
        lat_idx = int(round(event.ydata))
        lat_idx = max(0, min(lat_idx, summary.shape[0] - 1))
        lon_idx = max(0, min(lon_idx, summary.shape[1] - 1))
        if summary[lat_idx, lon_idx] <= 2:
            return  # invalid cell — ignore click
        selected[0], selected[1] = lat_idx, lon_idx
        sel_marker.set_data([lon_idx], [lat_idx])
        update_time_series(lat_idx, lon_idx)


    fig.canvas.mpl_connect("button_press_event", on_click)

    sel_marker.set_data([selected[1]], [selected[0]])
    update_time_series(selected[0], selected[1])

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot results of phenology computation')
    parser.add_argument('--extract_file', '-e', help='Absolute path of extract file')
    parser.add_argument('--phenology_file', '-p', help='Absolute path of phenology file')
    args = parser.parse_args()
    plot(args.extract_file, args.phenology_file)
