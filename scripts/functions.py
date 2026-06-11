import os
import sys
import logging
import netCDF4
import numpy as np
from csaps import csaps
from scipy.signal import find_peaks
from datetime import datetime, timezone
import warnings
from scipy.sparse import SparseEfficiencyWarning
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from shapely.geometry import Point
import pandas as pd
import plotly.graph_objects as go
from matplotlib.patches import Rectangle



warnings.filterwarnings("ignore", category=SparseEfficiencyWarning)

def parse_args(args):
    default_args = {
        "variable": "chla_mean", #
        "qa": "lwlr_quality_flag",
        "shapefile": "metadata.shp",
        "start_index": False,
        "end_index": False,
        "lakes": [],
        "out_folder": "",
        "images": "",  #
        "extract": True,
        "phenology": True,
        "analysis": False,
        "maps": False,
        "pixel_plots": False,
        "comparison": False,
        "comparison_classes": ["chla21", "chla31", "phycocyanin31"],
        "comparison_plot_types": ["chla21 vs chla31", "chla21 vs phyco", "chla31 vs phyco", "triple"],
        "background_pts": True,
        "purple_chla21": True,
        "time_splits" : [(0,9999)], 
        # "start":0,
        # "end":9999,
        # "split_start": 2016,
        # "split_end": 2012,
        "aggregation": True,
        "qa_filter": True, # Only accept qa_flag = 0
        "spline_min_phase_length": 14,
        "spline_min_relative_amplitude": 0,
        "spline_min_phase_data": 0,
        "spline_data_gap_size": 31,
        "spline_data_gap_size_buffer": 0,
        "spline_subs_peak_win_size": 365,
        "spline_subs_peak_ampl_frac": 0.05, # 0.35 removes lower peaks & 0.05 adds lower peaks
    }
    return default_args | args

def chunked(iterable, n):
    for i in range(0, len(iterable), n):
        yield iterable[i:i + n]

def set_logging(save):
    if save:
        filename = datetime.now().strftime("%Y%m%d_%H%M%S") + ".log"
        repo = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
        logging.basicConfig(
            filename=os.path.join(repo, "logs", filename),
            filemode='a',
            format='%(asctime)s - %(levelname)s - %(message)s',
            level=logging.INFO
        )
    else:
        logging.basicConfig(
            format='%(asctime)s - %(levelname)s - %(message)s',
            level=logging.INFO
        )

def verify_arg_file(value):
    arg_folder = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "../args"))
    for file in os.listdir(arg_folder):
        if os.path.splitext(file)[0] == value or file == value:
            return os.path.join(arg_folder, file)
    raise ValueError("Argument file {} not found in the args folder.".format(value))


def _datenum_to_unix(datenum):
    """Convert MATLAB-style datenum values to Unix timestamps."""
    # datenum is days since Jan 0, year 0. Ordinal is days since Jan 1, year 1.
    # datenum = ordinal + 366
    datenum = np.asarray(datenum, dtype=np.float64)
    ordinals = datenum - 366
    base = np.floor(ordinals).astype(int)
    frac = ordinals - base
    return np.array([datetime.fromordinal(int(b)).replace(tzinfo=timezone.utc).timestamp() + f * 86400
                     for b, f in zip(base, frac)])


def unix_to_datenum(arr):
    return np.array([datetime.utcfromtimestamp(int(ts)).toordinal() + 366 for ts in np.array(arr)])

def datenum_to_datetime(arr):
    return np.array([datetime.fromordinal(int(dn) - 366).replace(tzinfo=timezone.utc) for dn in np.array(arr)])

def unix_to_datetime(arr):
    return np.array([datetime.fromtimestamp(int(ts), tz=timezone.utc) for ts in np.array(arr)])

def remove_nan(arr):
    arr = np.array(arr)
    return arr[~np.isnan(arr)]

# define year range for plotting 
def define_year_range(start, end, years):
                return(years.min() if start == 0 else start, 
                       years.max() if end == 9999 else end)

# find closest factor for nice plotting
def close_factors(number):
    ''' 
    find the closest pair of factors for a given number
    '''
    factor1 = 0
    factor2 = number
    while factor1 +1 <= factor2:
        factor1 += 1
        if number % factor1 == 0:
            factor2 = number // factor1
        
    return factor1, factor2

def plot_lake_outline(geometry, ax):
    if geometry.geom_type == "Polygon":
        x, y = geometry.exterior.xy
        ax.plot(x, y, color="black", linewidth=1, label="Lake Outline")
    elif geometry.geom_type == "MultiPolygon":
        for i, poly in enumerate(geometry.geoms):
            x, y = poly.exterior.xy
            ax.plot(x, y, color="black", linewidth=1, label="Lake Outline" if i == 0 else None)

def grab_metrics(e_path, metric_scores, buffered_geom_prep):
    with netCDF4.Dataset(e_path) as nc:
        summary = np.array(nc.variables["summary"][:, :])
        lats = nc.variables["lat"][:]
        lons = nc.variables["lon"][:]
        map_data = np.full(summary.shape, np.nan)
        extent = [lons.min(), lons.max(), lats.min(), lats.max()]
        for (i, j), r2 in metric_scores.items():
            lon = lons[j]
            lat = lats[i]
            if buffered_geom_prep.contains(Point(lon, lat)):
                map_data[i, j] = r2
        return map_data, extent

def plot_map_data(colormap, map_data, extent, ax, cmap_extent=None):
    cmap = colormap if colormap else "RdYlBu"
    if cmap_extent is not None:
        im = ax.imshow(map_data, cmap=cmap, aspect="auto", origin="lower",
                       vmin=cmap_extent[0], vmax=cmap_extent[1], extent=extent)
    else:
        im = ax.imshow(map_data, cmap=cmap, aspect="auto", origin="lower", extent=extent)
    return im

def set_labels(ax, fig, im, title, colorbar_label, colorbar_ticks=None):
    ax.set_xlabel("Lon index")
    ax.set_ylabel("Lat index")
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax, label=colorbar_label)
    if colorbar_ticks is not None:
        cbar.set_ticks(colorbar_ticks)
    ax.legend()

def grab_time_data(e_path, p_path, valid_coords, buffered_geom_prep,
                   var_x, var_y, year, use_max):
    with netCDF4.Dataset(e_path) as nc:
        summary = np.array(nc.variables["summary"][:, :])
        lats = nc.variables["lat"][:]
        lons = nc.variables["lon"][:]

    map_data = np.full(summary.shape, np.nan)
    neg_warned = False
    extent = [lons.min(), lons.max(), lats.min(), lats.max()]

    with netCDF4.Dataset(p_path) as nc_p:
        x_full = nc_p.variables[var_x][:, :, :]
        y_full = nc_p.variables[var_y][:, :, :]

    for (i, j) in valid_coords:
        lon, lat = lons[j], lats[i]
        if not buffered_geom_prep.contains(Point(lon, lat)):
            continue
        x_arr = np.array(unix_to_datetime(remove_nan(x_full[i, j, :])))
        y_arr = np.array(remove_nan(y_full[i, j, :]))
        if len(x_arr) == 0:
            continue
        year_arr = np.array([d.year for d in x_arr])
        doys = np.array([d.timetuple().tm_yday for d in x_arr])
        mask = (year_arr == year) & (doys >= 160) & (doys <= 250)
        x_sub, y_sub = x_arr[mask], y_arr[mask]
        if len(x_sub) == 0:
            continue
        if (y_sub < 0).any() and not neg_warned:
            warnings.warn(f"Negative values in {year}", Warning)
            neg_warned = True
        if use_max:
            map_data[i, j] = float(x_sub[int(np.argmax(y_sub))].timetuple().tm_yday)
        else:
            map_data[i, j] = float(x_sub[0].timetuple().tm_yday)

    return map_data, extent



def save_maps(eda_instance, lake_analysis_folder, lake_str,  time_splits, metric= ["R2", "MAD", "RMSE", "correlation", "values_per_pixel"]):
    for m in metric:
        text_strings = []
        num_splits = len(time_splits)
        out_path = os.path.join(lake_analysis_folder, lake_str, "plots", "metric_maps")
        os.makedirs(out_path, exist_ok= True)
    
        if m == "R2":
            score_fn = eda_instance.r2_scores
            metric_str = "R$^2$"
            cmap = "RdYlBu"
            colorbar_extent = [0, 1]
        elif m == "MAD":
            score_fn = eda_instance.MAD_scores
            metric_str = "MAD"
            cmap = "RdYlBu_r"
            colorbar_extent = [0, 1]
        elif m == "RMSE":
            score_fn = eda_instance.RMSE_scores
            metric_str = "RMSE"
            cmap = "RdYlBu_r"
            colorbar_extent = [0, 1]
        elif m == "correlation":
            score_fn = eda_instance.correlation_scores
            metric_str = "correlation"
            cmap = "RdYlBu"
            colorbar_extent = [0, 1]
        elif m == "values_per_pixel":
            score_fn = eda_instance.values_per_pixel
            metric_str = "values_per_pixel"
            cmap = "winter"
            colorbar_extent = None
        else:
            raise ValueError("Please provide a valid metric string")

        for start, end in time_splits:
            if start > end:
                raise ValueError("Beginning of time split cannot be larger than the end")

        rows, cols = close_factors(num_splits)
        fig, axs = plt.subplots(rows, cols, constrained_layout=True)

        for i, (start, end) in enumerate(time_splits):
            if start == 0 and end == 9999:
                text_strings.append(f"{metric_str}-scores {eda_instance.variable} {eda_instance.version} full ts")
            elif start == 0:
                text_strings.append(f"{metric_str}-scores {eda_instance.variable} {eda_instance.version} {2002}-{end}")
            elif end == 9999:
                text_strings.append(f"{metric_str}-scores {eda_instance.variable} {eda_instance.version} {start}-{2024}")
            else:
                text_strings.append(f"{metric_str}-scores {eda_instance.variable} {eda_instance.version} {start}-{end}")

            scores = score_fn([(start, end)])

            if num_splits == 1:
                eda_instance.metric_map(scores, metric_str, fig, axs, cmap, colorbar_extent=colorbar_extent)
                axs.set_title(text_strings[0])
            else:
                eda_instance.metric_map(scores, metric_str, fig, axs[i], cmap, colorbar_extent=colorbar_extent)
                axs[i].set_title(text_strings[i])

        if num_splits == 1 and time_splits[0] == (0, 9999):
            file_name = f"{eda_instance.variable}_v{eda_instance.version.replace('.', '')}_{m}_full_ts.png"
        elif num_splits == 2:
            file_name = f"{eda_instance.variable}_v{eda_instance.version.replace('.', '')}_{m}_split_ts.png"
        else:
            file_name = f"{eda_instance.variable}_v{eda_instance.version.replace('.', '')}_{m}_{num_splits}_split_ts.png"
        fig.savefig(os.path.join(out_path, file_name), dpi=600)
        plt.close(fig)
       

def save_pixel_plots(eda_instance, pixels, lake_analysis_folder, lake_str, time_splits, aggregation = True):
    for i,j in pixels:
        if aggregation:
            out_path = os.path.join(lake_analysis_folder, lake_str, "plots", "pixel_plots", "aggregated", f"{i}_{j}")
        else:
            out_path = os.path.join(lake_analysis_folder, lake_str, "plots", "pixel_plots", "not_aggregated", f"{i}_{j}")
        os.makedirs(out_path, exist_ok=True)
        fig, ax = plt.subplots(1,1)
        eda_instance.pixel_map(i,j, ax)
        
        fig.savefig(os.path.join(out_path, "location.png"),dpi=600)
        plt.close(fig)

        num_splits = len(time_splits)
        rows, cols = close_factors(num_splits)

        # Single plots
        fig, axs = plt.subplots(rows, cols, constrained_layout=True)
        for num, (start, end) in enumerate(time_splits):
            if start > end:
                raise ValueError("Beginning of time split cannot be larger than the end")
            ax = axs if num_splits == 1 else np.atleast_1d(axs).flatten()[num]
            eda_instance.single_plot(ax=ax, latitude=i, longitude=j, aggregation=aggregation, start=start, end=end)
            ax.set_ylim(bottom=-0.5)
            

        if num_splits == 1:
            start, end = time_splits[0]
            if start == 0 and end == 9999:
                file_name = f"{eda_instance.variable}_v{eda_instance.version.replace('.', '')}_full_ts.png"
            elif start == 0:
                file_name = f"{eda_instance.variable}_v{eda_instance.version.replace('.', '')}_ts_{2002}_to_{end}.png"
            elif end == 9999:
                file_name = f"{eda_instance.variable}_v{eda_instance.version.replace('.', '')}_ts_{start}_to_{2024}.png"
            else:
                file_name = f"{eda_instance.variable}_v{eda_instance.version.replace('.', '')}_ts_{start}_to_{end}.png"
        elif num_splits == 2:
            file_name = f"{eda_instance.variable}_v{eda_instance.version.replace('.', '')}_split_ts.png"
        else:
            file_name = f"{eda_instance.variable}_v{eda_instance.version.replace('.', '')}_{num_splits}_split_ts.png"
        fig.savefig(os.path.join(out_path, file_name), dpi=600)
        plt.close(fig)

        # Peaks plots
        fig, axs = plt.subplots(rows, cols, constrained_layout=True)
        for num, (start, end) in enumerate(time_splits):
            if start > end:
                raise ValueError("Beginning of time split cannot be larger than the end")
            ax = axs if num_splits == 1 else np.atleast_1d(axs).flatten()[num]
            eda_instance.extrema_plot(ax=ax, latitude=i, longitude=j, aggregation=aggregation, start=start, end=end, peak=True)
            ax.set_ylim(bottom=-0.5)
         

        if num_splits == 1:
            start, end = time_splits[0]
            if start == 0 and end == 9999:
                file_name = f"{eda_instance.variable}_v{eda_instance.version.replace('.', '')}_peaks_full_ts.png"
            elif start == 0:
                file_name = f"{eda_instance.variable}_v{eda_instance.version.replace('.', '')}_peaks_ts_{2002}_to_{end}.png"
            elif end == 9999:
                file_name = f"{eda_instance.variable}_v{eda_instance.version.replace('.', '')}_peaks_ts_{start}_to_{2024}.png"
            else:
                file_name = f"{eda_instance.variable}_v{eda_instance.version.replace('.', '')}_peaks_ts_{start}_to_{end}.png"
        elif num_splits == 2:
            file_name = f"{eda_instance.variable}_v{eda_instance.version.replace('.', '')}_peaks_split_ts.png"
        else:
            file_name = f"{eda_instance.variable}_v{eda_instance.version.replace('.', '')}_{num_splits}_peaks_split_ts.png"
        fig.savefig(os.path.join(out_path, file_name), dpi=600)
        plt.close(fig)



def create_summary(eda_instance, pixels, lake_analysis_folder, lake_str, time_splits, summary_types=["n_peaks", "r2", "rmse", "mad", "correlation", "values", "qa"]):

    with netCDF4.Dataset(eda_instance.p_path) as nc:
        smoothing_arr = nc.variables["smoothing_parameter"][:]

    for i, j in pixels:
        out_path = os.path.join(lake_analysis_folder, lake_str, "plots", "pixel_plots", "summaries", f"{i}_{j}")
        os.makedirs(out_path, exist_ok=True)
        txt_path = os.path.join(out_path, f"summary_{eda_instance.variable}.txt")

        with open(txt_path, "w") as file:
            file.write(f"{eda_instance.variable} {eda_instance.version}:\n")
            file.write(f"    Smoothing Parameter: {float(smoothing_arr[i, j])}\n")
            for start, end in time_splits:
                if start == 0 and end == 9999:
                    label = "full ts"
                elif start == 0:
                    label = f"2002-{end}"
                elif end == 9999:
                    label = f"{start}-2024"
                else:
                    label = f"{start}-{end}"

                file.write(f"    {label}:\n")
                if "n_peaks" in summary_types:
                    file.write(f"        Number of Peaks: {eda_instance.count_extrema(i, j, start=start, end=end, peaks = True)}\n")
                if "r2" in summary_types:
                    file.write(f"        R2: {round(eda_instance.pixel_r2(i, j, start=start, end=end), 4)}\n")
                if "rmse" in summary_types:
                    file.write(f"        RMSE: {round(eda_instance.pixel_rmse(i, j, start=start, end=end), 4)}\n")
                if "mad" in summary_types:
                    file.write(f"        MAD: {round(eda_instance.pixel_mad(i, j, start=start, end=end), 4)}\n")
                if "correlation" in summary_types:
                    file.write(f"        Correlation: {round(eda_instance.pixel_correlation(i, j, start=start, end=end), 4)}\n")
                if "values" in summary_types:
                    file.write(f"        Values: {round(eda_instance.pixel_values(i, j, start=start, end=end), 4)}\n")

            if "qa" in summary_types:
                px = eda_instance._load_pixel_data(i, j)
                pks_x  = px["pks_x"]
                pks_qa = px["pks_qa"]

                trgs_x  = px["trgs_x"]
                trgs_qa = px["trgs_qa"]

                qa_periods = [
                    ("full ts",   0,    9999),
                    ("2002-2012", 2002, 2012),
                    ("2016-2024", 2016, 2024),
                ]
                file.write(f"    Peak QA:\n")
                file.write(f"        {'Period':<14}  {'Good':>6}  {'Fair':>6}  {'Poor':>6}\n")
                for period_label, s, e in qa_periods:
                    mask   = np.array([(d.year >= s) & (d.year <= e) for d in pks_x])
                    qa_sub_pks = pks_qa[mask]
                    file.write(f"        {period_label:<14}  {int(np.sum(qa_sub_pks == 0)):>6}  {int(np.sum(qa_sub_pks == 1)):>6}  {int(np.sum(qa_sub_pks == 2)):>6}\n")

                file.write(f"    Trough QA:\n")
                file.write(f"        {'Period':<14}  {'Good':>6}  {'Fair':>6}  {'Poor':>6}\n")
                for period_label, s, e in qa_periods:
                    mask   = np.array([(d.year >= s) & (d.year <= e) for d in trgs_x])
                    qa_sub_trgs = trgs_qa[mask]
                    file.write(f"        {period_label:<14}  {int(np.sum(qa_sub_trgs == 0)):>6}  {int(np.sum(qa_sub_trgs == 1)):>6}  {int(np.sum(qa_sub_trgs == 2)):>6}\n")

            file.write("\n")




def save_comparison_plots(instances, pixels, lake_analysis_folder, lake_str,
                          time_splits, comparison_plot_types,
                          aggregation=True, background_pts=False, purple_chla21=False):
    chla21 = instances.get("chla21")
    chla31 = instances.get("chla31")
    phyco  = instances.get("phycocyanin31")
    agg = "_agg" if aggregation else ""

    for start, end in time_splits:
        if start > end:
            raise ValueError("Beginning of time split cannot be larger than the end")

    pair_plots = []
    if "chla21 vs chla31" in comparison_plot_types and chla21 and chla31:
        pair_plots.append((chla21, chla31, None, "chla21_chla31"))
    if "chla21 vs phyco" in comparison_plot_types and chla21 and phyco:
        pair_plots.append((chla21, phyco, None, "chla21_phyco"))
    if "chla31 vs phyco" in comparison_plot_types and chla31 and phyco:
        pair_plots.append((chla31, phyco, None, "chla31_phyco"))
    if "triple" in comparison_plot_types and chla21 and chla31 and phyco:
        pair_plots.append((chla21, chla31, phyco, "chla21_chla31_phyco"))

    n = len(time_splits)

    def _year(v, zero_val=2002, max_val=2024):
        return zero_val if v == 0 else (max_val if v == 9999 else v)

    if n == 1 and time_splits[0] ==[0, 9999]:
        ts_suffix = "full_ts"
    elif n == 2:
        s1, e1 = time_splits[0]
        s2, e2 = time_splits[1]
        ts_suffix = f"{_year(s1)}_{_year(e1)}_vs_{_year(s2)}_{_year(e2)}"
    else:
        ts_suffix = f"{n}_split_ts"

    for i, j in pixels:
        save_path = os.path.join(lake_analysis_folder, lake_str, "plots", "pixel_plots", "comparisons", f"{i}_{j}")
        os.makedirs(save_path, exist_ok=True)

        for inst1, inst2, inst3, label in pair_plots:
            file_name = f"comparison_pks_{label}_{ts_suffix}{agg}.png"

            rows, cols = close_factors(n)
            fig, axs = plt.subplots(rows, cols, figsize=(15 * cols, 5 * rows), constrained_layout=True)
            axs_flat = np.atleast_1d(axs).flatten()

            for k, (start, end) in enumerate(time_splits):
                inst1.extrema_comparison(
                    inst2, i, j, axs_flat[k],
                    aggregation=aggregation, start=start, end=end,
                    background_pts=background_pts,
                    other2=inst3,
                    purple_chla21=purple_chla21
                )

            fig.savefig(os.path.join(save_path, file_name), dpi=600)
            plt.close(fig)


def prep_dimark_data(insitu_df,start=0, end=9999,  insitu_date_col="datetime", insitu_value_col="chlorophyll_a", insitu_station_col=None, station_id=None, max_depth = 5):
        insitu = insitu_df.copy()
        insitu[insitu_date_col] = pd.to_datetime(
                insitu[insitu_date_col], format = "mixed", dayfirst = True)


        # optional station filtering
        if (insitu_station_col is not None) and (station_id is not None):
                insitu = insitu[
                insitu[insitu_station_col] == station_id
                ]

        # time filtering
        insitu = insitu[
                (insitu[insitu_date_col].dt.year >= start) &
                (insitu[insitu_date_col].dt.year <= end)
        ]

        # remove invalid values
        insitu = insitu[
                np.isfinite(insitu[insitu_value_col])]
        insitu[insitu_date_col] = pd.to_datetime(insitu["datetime"]).dt.date
        insitu = insitu[insitu["depth"]< max_depth]
        insitu_mean = (insitu.groupby("datetime", as_index=False)[insitu_value_col].mean())
        return insitu_mean





def bivariate_legend( ax, color_set):

    labels = ["0", "1", "2", "2+"]

    cell_size = 1

    for row in range(4):
        for col in range(4):

            idx = row * 4 + col

            ax.add_patch(
                Rectangle(
                    (col, row),          # x,y
                    cell_size,           # width
                    cell_size,           # height
                    facecolor=color_set[idx],
                    edgecolor="black"
                )
            )

    ax.set_xlim(0, 4)
    ax.set_ylim(0, 4)

    ax.set_xticks([0.5, 1.5, 2.5, 3.5])
    ax.set_xticklabels(labels)

    ax.set_yticks([0.5, 1.5, 2.5, 3.5])
    ax.set_yticklabels(labels)

    ax.set_xlabel("# Troughs")
    ax.set_ylabel("# Peaks")

    ax.set_aspect("equal")



def create_empty_heatmap(nrows=4, ncols=4):

    data = np.full((nrows, ncols), np.nan)

    fig, ax = plt.subplots(figsize=(5, 5))

    ax.set_xlim(0, 3)
    ax.set_ylim(int(2001), int(2024))

    # Draw grid lines
    ax.set_xticks(np.arange(0, 5, 1), minor=True)
    ax.grid(which='minor')
    ax.set_yticks(np.arange(2001, 2024), minor = True)


    ax.set_xticks(np.arange(0.5, 4, 1))
    ax.set_yticks(np.arange(2001.5, 2024.5, 1))
    months = ["Jan-Mar", "Apr-Jun", "Jul-Sep", "Oct-Dez"]
    years = ["2002", "2003", "2004", "2005", 
             "2006", "2007", "2008", "2009", "2010", 
             "2011", "2012", "2013", "2014", "2015", 
             "2016", "2017", "2018", "2019","2020", 
             "2021", "2022", "2023", "2024"]
    ax.set_xticklabels(months, horizontalalignment = "center")
    ax.set_yticklabels(years, horizontalalignment = "center")

    # Hide tick marks
    ax.tick_params(axis='x', length=0)
    ax.tick_params(axis='y', length=0, pad = 20)



    return fig, ax, data






def init_phenology_output(out, lat, lon, p=None):
    """Initialise a streaming phenology NetCDF file.

    Variables are shaped (lat, lon, record) where 'record' is an unlimited
    dimension.  Pixel identity is encoded by array position so no index
    variables are needed.  The record dimension grows automatically as pixels
    are written; unused slots for a given pixel are filled with the variable's
    fill value.

    Parameters
    ----------
    out : netCDF4.Dataset
        Open, writable NetCDF4 dataset.
    lat : array-like
        1-D array of latitude values for the grid.
    lon : array-like
        1-D array of longitude values for the grid.
    p : dict, optional
        Parameters dict.  If provided, spline parameters are written as global
        attributes so the settings used to produce the file are self-documented.
    """
    if p is not None:
        out.setncatts({
            'variable':                    p['variable'],
            'qa':                          p['qa'],
            'qa_filter':                   int(p['qa_filter']),
            'spline_min_phase_length':     p['spline_min_phase_length'],
            'spline_min_relative_amplitude': p['spline_min_relative_amplitude'],
            'spline_min_phase_data':       p['spline_min_phase_data'],
            'spline_data_gap_size':        p['spline_data_gap_size'],
            'spline_data_gap_size_buffer': p['spline_data_gap_size_buffer'],
            'spline_subs_peak_win_size':   p['spline_subs_peak_win_size'],
            'spline_subs_peak_ampl_frac':  p['spline_subs_peak_ampl_frac'],
        })
    out.createDimension('lat', len(lat))
    out.createDimension('lon', len(lon))
    out.createDimension('record', None)  # unlimited

    lat_var = out.createVariable('lat', 'f4', ('lat',))
    lat_var[:] = lat
    lat_var.units = 'degrees_north'
    lon_var = out.createVariable('lon', 'f4', ('lon',))
    lon_var[:] = lon
    lon_var.units = 'degrees_east'

    sm = out.createVariable('smoothing_parameter', 'f8', ('lat', 'lon'), fill_value=np.nan,
                            zlib=True, complevel=4)
    sm[:] = np.nan

    v = out.createVariable('pks_x', 'f8', ('lat', 'lon', 'record'), fill_value=np.nan,
                           zlib=True, complevel=4)
    v.long_name = 'Peak time (Unix)'
    v = out.createVariable('pks_y', 'f8', ('lat', 'lon', 'record'), fill_value=np.nan,
                           zlib=True, complevel=4)
    v.long_name = 'Peak value'
    v = out.createVariable('pks_qa', 'i2', ('lat', 'lon', 'record'), fill_value=-1,
                           zlib=True, complevel=4)
    v.long_name = 'Peak QA (0=Good, 1=Fair, 2=Poor)'

    v = out.createVariable('trgs_x', 'f8', ('lat', 'lon', 'record'), fill_value=np.nan,
                           zlib=True, complevel=4)
    v.long_name = 'Trough time (Unix)'
    v = out.createVariable('trgs_y', 'f8', ('lat', 'lon', 'record'), fill_value=np.nan,
                           zlib=True, complevel=4)
    v.long_name = 'Trough value'
    v = out.createVariable('trgs_qa', 'i2', ('lat', 'lon', 'record'), fill_value=-1,
                           zlib=True, complevel=4)
    v.long_name = 'Trough QA (0=Good, 1=Fair, 2=Poor)'

    for name, desc in [('green_up_onset_x',    'Green-up onset time (Unix)'),
                       ('green_up_onset_y',    'Green-up onset value'),
                       ('green_up_mid_x',      'Green-up mid time (Unix)'),
                       ('green_up_mid_y',      'Green-up mid value'),
                       ('green_up_advanced_x', 'Green-up advanced time (Unix)'),
                       ('green_up_advanced_y', 'Green-up advanced value'),
                       ('green_down_onset_x',    'Green-down onset time (Unix)'),
                       ('green_down_onset_y',    'Green-down onset value'),
                       ('green_down_mid_x',      'Green-down mid time (Unix)'),
                       ('green_down_mid_y',      'Green-down mid value'),
                       ('green_down_advanced_x', 'Green-down advanced time (Unix)'),
                       ('green_down_advanced_y', 'Green-down advanced value')]:
        v = out.createVariable(name, 'f8', ('lat', 'lon', 'record'), fill_value=np.nan,
                               zlib=True, complevel=4)
        v.long_name = desc

    v = out.createVariable('data_gap_start', 'f8', ('lat', 'lon', 'record'), fill_value=np.nan,
                           zlib=True, complevel=4)
    v.long_name = 'Data gap start (Unix)'
    v = out.createVariable('data_gap_end', 'f8', ('lat', 'lon', 'record'), fill_value=np.nan,
                           zlib=True, complevel=4)
    v.long_name = 'Data gap end (Unix)'


def append_pixel_phenology(out, x, y, result, pheno):
    """Write one pixel's phenology results into an open streaming NetCDF file.

    Each feature type is written into its (lat, lon, record) slice starting at
    index 0.  The unlimited record dimension expands automatically whenever a
    pixel has more records than the current dimension size.

    Parameters
    ----------
    out : netCDF4.Dataset
        Open, writable dataset previously initialised with init_phenology_output.
    x : int
        Lat index of this pixel.
    y : int
        Lon index of this pixel.
    result : dict
        Output from smooth_cubic_spline (must have 'smoothing' key).
    pheno : dict
        Output from extract_phenology_metrics.
    """
    out.variables['smoothing_parameter'][x, y] = result['smoothing']

    n_pk = len(pheno['pks_x'])
    if n_pk > 0:
        out.variables['pks_x'][x, y, :n_pk] = _datenum_to_unix(pheno['pks_x'])
        out.variables['pks_y'][x, y, :n_pk] = pheno['pks_y']
        out.variables['pks_qa'][x, y, :n_pk] = pheno['pks_qa']

    n_trg = len(pheno['trgs_x'])
    if n_trg > 0:
        out.variables['trgs_x'][x, y, :n_trg] = _datenum_to_unix(pheno['trgs_x'])
        out.variables['trgs_y'][x, y, :n_trg] = pheno['trgs_y']
        out.variables['trgs_qa'][x, y, :n_trg] = pheno['trgs_qa']

    n_gu = len(pheno['green_up'])
    if n_gu > 0:
        gu_all_x = np.array([[g['onset_x'], g['mid_x'], g['advanced_x']] for g in pheno['green_up']])
        gu_all_x_unix = _datenum_to_unix(gu_all_x.ravel()).reshape(n_gu, 3)
        out.variables['green_up_onset_x'][x, y, :n_gu] = gu_all_x_unix[:, 0]
        out.variables['green_up_onset_y'][x, y, :n_gu] = [g['onset_y'] for g in pheno['green_up']]
        out.variables['green_up_mid_x'][x, y, :n_gu] = gu_all_x_unix[:, 1]
        out.variables['green_up_mid_y'][x, y, :n_gu] = [g['mid_y'] for g in pheno['green_up']]
        out.variables['green_up_advanced_x'][x, y, :n_gu] = gu_all_x_unix[:, 2]
        out.variables['green_up_advanced_y'][x, y, :n_gu] = [g['advanced_y'] for g in pheno['green_up']]

    n_gd = len(pheno['green_down'])
    if n_gd > 0:
        gd_all_x = np.array([[g['onset_x'], g['mid_x'], g['advanced_x']] for g in pheno['green_down']])
        gd_all_x_unix = _datenum_to_unix(gd_all_x.ravel()).reshape(n_gd, 3)
        out.variables['green_down_onset_x'][x, y, :n_gd] = gd_all_x_unix[:, 0]
        out.variables['green_down_onset_y'][x, y, :n_gd] = [g['onset_y'] for g in pheno['green_down']]
        out.variables['green_down_mid_x'][x, y, :n_gd] = gd_all_x_unix[:, 1]
        out.variables['green_down_mid_y'][x, y, :n_gd] = [g['mid_y'] for g in pheno['green_down']]
        out.variables['green_down_advanced_x'][x, y, :n_gd] = gd_all_x_unix[:, 2]
        out.variables['green_down_advanced_y'][x, y, :n_gd] = [g['advanced_y'] for g in pheno['green_down']]

    n_gap = len(pheno['data_gap_start_days'])
    if n_gap > 0:
        out.variables['data_gap_start'][x, y, :n_gap] = _datenum_to_unix(pheno['data_gap_start_days'])
        out.variables['data_gap_end'][x, y, :n_gap] = _datenum_to_unix(pheno['data_gap_end_days'])


def smooth_cubic_spline(time, values, spline_min_phase_length, spline_min_relative_amplitude, spline_min_phase_data,
                        smoothing_min=0.0, smoothing_max=1.0, smoothing_change=1e-12, max_iterations=1e5, smooth_x_axis=None):
    """Find the optimal smoothing parameter for a cubic spline using binary search.

    Searches for the highest smoothing parameter that still satisfies the minimum
    phase length, relative amplitude, and data points per phase conditions.

    Parameters
    ----------
    time : array-like
        Numeric time values of observations.
    values : array-like
        Observed data values.
    spline_min_phase_length : int
        Minimum allowed phase length in days.
    spline_min_relative_amplitude : float
        Minimum relative amplitude between consecutive peaks/troughs (0-1).
    spline_min_phase_data : int
        Minimum number of data points per phase.
    smoothing_min : float, optional
        Lower bound for smoothing parameter (default 0.0).
    smoothing_max : float, optional
        Upper bound for smoothing parameter (default 1.0).
    smoothing_change : float, optional
        Convergence threshold for the smoothing parameter (default 1e-12).
    max_iterations : float, optional
        Maximum number of binary search iterations (default 1e5).

    Returns
    -------
    dict
        Output from cubic_spline with the optimal smoothing parameter, or a dict
        with None values and conditions_satisfied=False if no valid fit is found.

    Raises
    ------
    ValueError
        If max_iterations is reached without convergence.
    """
    smoothing = smoothing_min + 0.5 * (smoothing_max - smoothing_min)
    smoothing_old = smoothing
    iteration = 0

    while True:
        # Binary search for smoothing parameter
        if iteration >= max_iterations:
            raise ValueError('No fit found: maximum number of iterations reached')

        result = cubic_spline(
            time,
            values,
            spline_min_phase_length,
            spline_min_relative_amplitude,
            spline_min_phase_data,
            smoothing,
            smooth_x_axis
        )

        if result['conditions_satisfied']:
            if abs(smoothing_old - smoothing) < smoothing_change:
                return result
            smoothing_min = smoothing
        else:
            smoothing_max = smoothing

        smoothing_old = smoothing
        smoothing = smoothing_min + 0.5 * (smoothing_max - smoothing_min)

        if smoothing_max <= smoothing_change:
            # Flat line
            return {
                'smooth_y_data': None,
                'min_phase_length': None,
                'min_rel_amp': None,
                'min_phase_nr_data': None,
                'x_pks': None,
                'y_pks': None,
                'x_trgs': None,
                'y_trgs': None,
                'smoothing': None,
                'conditions_satisfied': False
            }
        iteration += 1


def cubic_spline(time, values, cond_min_phase_length, cond_min_rel_amp, cond_min_phase_nr_data, smoothing, smooth_x_axis=None):
    """Fit a smoothing cubic spline and evaluate phase conditions.

    Fits a cubic smoothing spline to the data, finds peaks and troughs, and checks
    whether the resulting phases satisfy the minimum phase length, relative amplitude,
    and data points per phase conditions.

    Parameters
    ----------
    time : array-like
        Numeric time values of observations.
    values : array-like
        Observed data values.
    cond_min_phase_length : int
        Minimum allowed phase length in days.
    cond_min_rel_amp : float
        Minimum relative amplitude between consecutive peaks/troughs (0-1).
    cond_min_phase_nr_data : int
        Minimum number of data points per phase.
    smoothing : float
        Smoothing parameter for csaps (0 = interpolation, 1 = least-squares line).

    Returns
    -------
    dict
        smooth_x_axis : ndarray - regular daily grid used for spline evaluation
        smooth_y_data : ndarray - spline values on the grid
        min_phase_length : float or None - shortest phase length found
        min_rel_amp : float or None - smallest relative amplitude found
        min_phase_nr_data : int - fewest data points in any phase
        x_pks, y_pks : ndarray - peak positions and values
        x_trgs, y_trgs : ndarray - trough positions and values
        smoothing : float - the smoothing parameter used
        conditions_satisfied : bool - whether all conditions are met
    """
    time = np.asarray(time)
    values = np.asarray(values)
    if smooth_x_axis is None:
        smooth_x_axis = np.arange(time.min(), time.max() + 1, 1)

    # Fit smoothing spline using numeric values
    smooth_y_data = csaps(time, values, smooth_x_axis, smooth=smoothing)

    # Find peaks and troughs
    pks_indices, _ = find_peaks(smooth_y_data)
    y_pks = smooth_y_data[pks_indices]
    x_pks = smooth_x_axis[pks_indices]

    trgs_indices, _ = find_peaks(-smooth_y_data)
    y_trgs = smooth_y_data[trgs_indices]
    x_trgs = smooth_x_axis[trgs_indices]

    # Initialize outputs
    min_phase_length = None
    min_rel_amp = None
    min_phase_nr_data = len(time)
    conditions_satisfied = False

    # Combine and sort peaks and troughs
    x_pks_and_trgs = np.concatenate([x_pks, x_trgs])
    y_pks_and_trgs = np.concatenate([y_pks, y_trgs])

    if len(x_pks_and_trgs) >= 2:
        sort_indices = np.argsort(x_pks_and_trgs)
        x_pks_and_trgs_sorted = x_pks_and_trgs[sort_indices]
        y_pks_and_trgs_sorted = y_pks_and_trgs[sort_indices]

        # Minimum phase length (in days)
        min_phase_length = np.min(np.diff(x_pks_and_trgs_sorted))

        # Minimum relative amplitude
        y_range = np.max(smooth_y_data) - np.min(smooth_y_data)
        if y_range > 0:
            min_rel_amp = np.min(np.abs(np.diff(y_pks_and_trgs_sorted))) / y_range
        else:
            min_rel_amp = 0

        # Minimum data points per phase
        for i in range(len(x_pks_and_trgs_sorted) - 1):
            count = np.sum(
                (time >= x_pks_and_trgs_sorted[i]) &
                (time <= x_pks_and_trgs_sorted[i + 1])
            )
            min_phase_nr_data = min(min_phase_nr_data, count)

        # Check conditions
        if (min_phase_length >= cond_min_phase_length and
                min_rel_amp >= cond_min_rel_amp and
                min_phase_nr_data >= cond_min_phase_nr_data):
            conditions_satisfied = True

    return {
        'smooth_x_axis': smooth_x_axis,
        'smooth_y_data': smooth_y_data,
        'min_phase_length': min_phase_length,
        'min_rel_amp': min_rel_amp,
        'min_phase_nr_data': min_phase_nr_data,
        'x_pks': x_pks,
        'y_pks': y_pks,
        'x_trgs': x_trgs,
        'y_trgs': y_trgs,
        'smoothing': smoothing,
        'conditions_satisfied': conditions_satisfied
    }


def remove_non_substantial_phases(x_pks, y_pks, x_trgs, y_trgs, win_size, ampl_frac):
    """Iteratively remove non-substantial peaks and their neighboring troughs.

    Uses the MCD-style algorithm: repeatedly finds the lowest-amplitude peak,
    compares its amplitude to the overall amplitude within a surrounding window,
    and removes it (along with its highest neighboring trough) if the peak is
    not substantial (amplitude <= ampl_frac * window_amplitude).

    Parameters
    ----------
    x_pks : array-like
        X positions of peaks.
    y_pks : array-like
        Y values of peaks.
    x_trgs : array-like
        X positions of troughs.
    y_trgs : array-like
        Y values of troughs.
    win_size : int
        Window size (in days) added on each side of the peak-trough pair to
        compute the reference amplitude.
    ampl_frac : float
        Fraction of window amplitude below which a peak is considered
        non-substantial and removed.

    Returns
    -------
    tuple of ndarray
        (x_pks, y_pks, x_trgs, y_trgs) with non-substantial phases removed.
    """
    x_pks = list(x_pks)
    y_pks = list(y_pks)
    x_trgs = list(x_trgs)
    y_trgs = list(y_trgs)

    while True:
        nr_peaks_initial = len(x_pks)
        nr_peaks_remain = 0

        while True:
            if nr_peaks_remain >= len(x_pks):
                break

            if len(x_trgs) == 0 or len(x_pks) == 0:
                x_pks = []
                y_pks = []
                x_trgs = []
                y_trgs = []
                break

            # Find lowest peak that was not skipped
            y_sort_indices = np.argsort(y_pks)
            pk_idx = y_sort_indices[nr_peaks_remain]

            # Find highest corresponding neighboring trough
            starts_with_pk = x_pks[0] < x_trgs[0]
            ends_with_pk = x_pks[-1] > x_trgs[-1]

            if starts_with_pk:
                if ends_with_pk:
                    if pk_idx == 0:
                        trg_idx = 0
                    elif pk_idx == len(x_pks) - 1:
                        trg_idx = pk_idx - 1
                    else:
                        trg_idx = pk_idx - 1 if y_trgs[pk_idx - 1] >= y_trgs[pk_idx] else pk_idx
                else:
                    if pk_idx == 0:
                        trg_idx = 0
                    else:
                        trg_idx = pk_idx - 1 if y_trgs[pk_idx - 1] >= y_trgs[pk_idx] else pk_idx
            else:
                if ends_with_pk:
                    if pk_idx == len(x_pks) - 1:
                        trg_idx = pk_idx
                    else:
                        trg_idx = pk_idx if y_trgs[pk_idx] >= y_trgs[pk_idx + 1] else pk_idx + 1
                else:
                    trg_idx = pk_idx if y_trgs[pk_idx] >= y_trgs[pk_idx + 1] else pk_idx + 1

            # Phase amplitude
            phase_ampl = y_pks[pk_idx] - y_trgs[trg_idx]

            # Amplitude in window
            window_x_start = min(x_pks[pk_idx], x_trgs[trg_idx]) - win_size
            window_x_end = max(x_pks[pk_idx], x_trgs[trg_idx]) + win_size

            trgs_in_window = [y_trgs[i] for i in range(len(x_trgs)) if window_x_start <= x_trgs[i] <= window_x_end]
            pks_in_window = [y_pks[i] for i in range(len(x_pks)) if window_x_start <= x_pks[i] <= window_x_end]
            window_amplitude = max(pks_in_window) - min(trgs_in_window)

            # Remove when not substantial
            if phase_ampl <= ampl_frac * window_amplitude:
                del x_pks[pk_idx]
                del y_pks[pk_idx]
                del x_trgs[trg_idx]
                del y_trgs[trg_idx]
            else:
                nr_peaks_remain += 1

        if nr_peaks_remain >= nr_peaks_initial:
            break

    return np.array(x_pks), np.array(y_pks), np.array(x_trgs), np.array(y_trgs)


def extract_phenology_metrics(x_pks, y_pks, x_trgs, y_trgs,
                              time, smooth_x_axis, smooth_y_data,
                              cond_subs_peak_win_size, cond_subs_peak_ampl_frac,
                              cond_data_gap_size, cond_data_gap_size_buffer):
    """Extract phenology metrics from cubic spline peaks and troughs.

    Removes non-substantial phases, then extracts green-up and green-down
    phenology metrics (onset at 25%, mid at 50%, advanced at 75% thresholds).
    Also detects data gaps and assigns QA flags to peaks and troughs.

    Parameters
    ----------
    x_pks : ndarray
        X positions of peaks from the spline fit.
    y_pks : ndarray
        Y values of peaks from the spline fit.
    x_trgs : ndarray
        X positions of troughs from the spline fit.
    y_trgs : ndarray
        Y values of troughs from the spline fit.
    time : array-like
        Numeric time values of the original observations.
    smooth_x_axis : ndarray
        Regular daily grid the spline was evaluated on.
    smooth_y_data : ndarray
        Spline values on the grid.
    cond_subs_peak_win_size : int
        Window size (days) for substantial peak check.
    cond_subs_peak_ampl_frac : float
        Amplitude fraction threshold for substantial peak check.
    cond_data_gap_size : int
        Minimum gap size (days) between observations to flag as a data gap.
    cond_data_gap_size_buffer : int
        Buffer (days) added on each side of detected data gaps.

    Returns
    -------
    dict
        pks_x, pks_y : ndarray - substantial peak positions and values
        trgs_x, trgs_y : ndarray - substantial trough positions and values
        pks_qa, trgs_qa : ndarray - QA flags (0=Good, 1=Fair, 2=Poor)
        green_up : list of dict - green-up metrics per phase, each with keys:
            pk_nr, trg_nr, ampl, onset_x, onset_y, mid_x, mid_y,
            advanced_x, advanced_y
        green_down : list of dict - green-down metrics per phase (same keys)
        data_gap_start_days, data_gap_end_days : ndarray - gap boundaries
        smooth_qa_data : ndarray - per-day QA (0=Good, 2=in data gap)
    """
    # Remove non-substantial phases
    pks_x, pks_y, trgs_x, trgs_y = remove_non_substantial_phases(
        x_pks, y_pks, x_trgs, y_trgs,
        cond_subs_peak_win_size, cond_subs_peak_ampl_frac
    )

    nr_pks = len(pks_x)
    nr_trgs = len(trgs_x)

    if nr_pks == 0 or nr_trgs == 0:
        return {
            'pks_x': pks_x, 'pks_y': pks_y,
            'trgs_x': trgs_x, 'trgs_y': trgs_y,
            'pks_qa': np.zeros(nr_pks, dtype=int), 'trgs_qa': np.zeros(nr_trgs, dtype=int),
            'green_up': [], 'green_down': [],
            'data_gap_start_days': np.array([]),
            'data_gap_end_days': np.array([]),
            'smooth_qa_data': np.zeros(len(smooth_x_axis), dtype=int)
        }

    # Combine and sort peaks and troughs
    all_x = np.concatenate([pks_x, trgs_x])
    all_y = np.concatenate([pks_y, trgs_y])
    sort_indices = np.argsort(all_x)
    x_sorted = all_x[sort_indices]
    y_sorted = all_y[sort_indices]

    # Extract phenology metrics per phase
    nr_phases = len(x_sorted) - 1
    green_up = []
    green_down = []

    # Track which peak/trough index each sorted position maps to
    # In the concatenated array: 0..nr_pks-1 are peaks, nr_pks..end are troughs
    green_up_pk_nrs = []
    green_up_trg_nrs = []
    green_down_pk_nrs = []
    green_down_trg_nrs = []

    for phase_nr in range(nr_phases):
        phase_start_x = x_sorted[phase_nr]
        phase_end_x = x_sorted[phase_nr + 1]
        phase_start_y = y_sorted[phase_nr]
        phase_end_y = y_sorted[phase_nr + 1]

        # Determine peak/trough indices from the sort mapping
        start_orig_idx = sort_indices[phase_nr]
        end_orig_idx = sort_indices[phase_nr + 1]

        if phase_start_y < phase_end_y:
            # Green-up: trough -> peak
            pk_nr = end_orig_idx  # index in peaks array
            trg_nr = start_orig_idx - nr_pks  # index in troughs array

            ampl = phase_end_y - phase_start_y

            # Get spline segment for this phase
            phase_mask = (smooth_x_axis >= trgs_x[trg_nr]) & (smooth_x_axis <= pks_x[pk_nr])
            phase_sx = smooth_x_axis[phase_mask]
            phase_sy = smooth_y_data[phase_mask]

            # Onset: first day above 25% threshold
            threshold_25 = phase_start_y + 0.25 * ampl
            onset_idx = np.argmax(phase_sy >= threshold_25)
            onset_x = phase_sx[onset_idx]
            onset_y = phase_sy[onset_idx]

            # Mid: mean of first day above 50% and last day below 50%
            threshold_50 = phase_start_y + 0.5 * ampl
            mid_one_idx = np.argmax(phase_sy >= threshold_50)
            mid_two_idx = len(phase_sy) - 1 - np.argmax(phase_sy[::-1] <= threshold_50)
            mid_x = (phase_sx[mid_one_idx] + phase_sx[mid_two_idx]) / 2.0
            mid_y = (phase_sy[mid_one_idx] + phase_sy[mid_two_idx]) / 2.0

            # Advanced: last day below 75% threshold
            threshold_75 = phase_start_y + 0.75 * ampl
            advanced_idx = len(phase_sy) - 1 - np.argmax(phase_sy[::-1] <= threshold_75)
            advanced_x = phase_sx[advanced_idx]
            advanced_y = phase_sy[advanced_idx]

            green_up.append({
                'pk_nr': pk_nr, 'trg_nr': trg_nr, 'ampl': ampl,
                'onset_x': onset_x, 'onset_y': onset_y,
                'mid_x': mid_x, 'mid_y': mid_y,
                'advanced_x': advanced_x, 'advanced_y': advanced_y
            })
            green_up_pk_nrs.append(pk_nr)
            green_up_trg_nrs.append(trg_nr)
        else:
            # Green-down: peak -> trough
            pk_nr = start_orig_idx  # index in peaks array
            trg_nr = end_orig_idx - nr_pks  # index in troughs array

            ampl = phase_start_y - phase_end_y

            # Get spline segment for this phase
            phase_mask = (smooth_x_axis >= pks_x[pk_nr]) & (smooth_x_axis <= trgs_x[trg_nr])
            phase_sx = smooth_x_axis[phase_mask]
            phase_sy = smooth_y_data[phase_mask]

            # Onset: first day dropping to 75% threshold
            threshold_75 = phase_end_y + 0.75 * ampl
            onset_idx = np.argmax(phase_sy <= threshold_75)
            onset_x = phase_sx[onset_idx]
            onset_y = phase_sy[onset_idx]

            # Mid: mean of first day dropping to 50% and last day above 50%
            threshold_50 = phase_end_y + 0.5 * ampl
            mid_one_idx = np.argmax(phase_sy <= threshold_50)
            mid_two_idx = len(phase_sy) - 1 - np.argmax(phase_sy[::-1] >= threshold_50)
            mid_x = (phase_sx[mid_one_idx] + phase_sx[mid_two_idx]) / 2.0
            mid_y = (phase_sy[mid_one_idx] + phase_sy[mid_two_idx]) / 2.0

            # Advanced: last day above 25% threshold
            threshold_25 = phase_end_y + 0.25 * ampl
            advanced_idx = len(phase_sy) - 1 - np.argmax(phase_sy[::-1] >= threshold_25)
            advanced_x = phase_sx[advanced_idx]
            advanced_y = phase_sy[advanced_idx]

            green_down.append({
                'pk_nr': pk_nr, 'trg_nr': trg_nr, 'ampl': ampl,
                'onset_x': onset_x, 'onset_y': onset_y,
                'mid_x': mid_x, 'mid_y': mid_y,
                'advanced_x': advanced_x, 'advanced_y': advanced_y
            })
            green_down_pk_nrs.append(pk_nr)
            green_down_trg_nrs.append(trg_nr)

    # Detect data gaps
    time_sorted = np.sort(time)
    gap_sizes = np.diff(time_sorted)
    large_gap_mask = gap_sizes >= cond_data_gap_size
    gap_start_days = time_sorted[:-1][large_gap_mask] - cond_data_gap_size_buffer
    gap_end_days = time_sorted[1:][large_gap_mask] + cond_data_gap_size_buffer

    # Merge overlapping gaps
    if len(gap_start_days) > 1:
        overlap = gap_start_days[1:] <= gap_end_days[:-1] + 1
        gap_start_days = np.delete(gap_start_days, np.where(overlap)[0] + 1)
        gap_end_days = np.delete(gap_end_days, np.where(overlap)[0])

    # Find peaks and troughs in gaps, build smooth QA data
    pk_nrs_in_gap = []
    trg_nrs_in_gap = []
    smooth_qa_data = np.zeros(len(smooth_x_axis), dtype=int)

    for gap_nr in range(len(gap_start_days)):
        gs = gap_start_days[gap_nr]
        ge = gap_end_days[gap_nr]

        pk_nrs_in_gap.extend(np.where((pks_x >= gs) & (pks_x <= ge))[0].tolist())
        trg_nrs_in_gap.extend(np.where((trgs_x >= gs) & (trgs_x <= ge))[0].tolist())

        smooth_qa_data[(smooth_x_axis >= gs) & (smooth_x_axis <= ge)] = 2

    pk_nrs_in_gap = set(pk_nrs_in_gap)
    trg_nrs_in_gap = set(trg_nrs_in_gap)

    # Assign QA flags to peaks
    # QA: 0 = Good, 1 = Fair (neighboring trough in gap), 2 = Poor (peak itself in gap)
    pks_qa = np.zeros(nr_pks, dtype=int)
    for pk_nr in range(nr_pks):
        # Check if a neighboring trough is in a data gap -> Fair
        for gu in green_up:
            if gu['pk_nr'] == pk_nr and gu['trg_nr'] in trg_nrs_in_gap:
                pks_qa[pk_nr] = 1
        for gd in green_down:
            if gd['pk_nr'] == pk_nr and gd['trg_nr'] in trg_nrs_in_gap:
                pks_qa[pk_nr] = 1
        # Peak itself in gap -> Poor (overrides Fair)
        if pk_nr in pk_nrs_in_gap:
            pks_qa[pk_nr] = 2

    # Assign QA flags to troughs
    trgs_qa = np.zeros(nr_trgs, dtype=int)
    for trg_nr in range(nr_trgs):
        # Check if a neighboring peak is in a data gap -> Fair
        for gu in green_up:
            if gu['trg_nr'] == trg_nr and gu['pk_nr'] in pk_nrs_in_gap:
                trgs_qa[trg_nr] = 1
        for gd in green_down:
            if gd['trg_nr'] == trg_nr and gd['pk_nr'] in pk_nrs_in_gap:
                trgs_qa[trg_nr] = 1
        # Trough itself in gap -> Poor (overrides Fair)
        if trg_nr in trg_nrs_in_gap:
            trgs_qa[trg_nr] = 2 # There is a bug in matlab code

    return {
        'pks_x': pks_x, 'pks_y': pks_y,
        'trgs_x': trgs_x, 'trgs_y': trgs_y,
        'pks_qa': pks_qa, 'trgs_qa': trgs_qa,
        'green_up': green_up, 'green_down': green_down,
        'data_gap_start_days': gap_start_days,
        'data_gap_end_days': gap_end_days,
        'smooth_qa_data': smooth_qa_data
    }