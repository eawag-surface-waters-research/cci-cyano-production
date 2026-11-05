import pandas as pd
import netCDF4
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.lines as mlines
import matplotlib.colors as mcolors
import numpy as np
from pathlib import Path
import warnings
from matplotlib.colors import ListedColormap, BoundaryNorm
import functions as f
from shapely.geometry import Point
import colorcet as cc
from visualization.base import PhenologyVisualizationBase


class PixelVisualization(PhenologyVisualizationBase):
    """Pixel-level phenology graphs.

    This class groups the public methods that render or support graphs for a
    single pixel.  It intentionally inherits the shared implementation so
    existing callers can migrate incrementally.
    """

    graph_type = "pixel"
    graph_methods = frozenset({
        "pixel_map",
        "single_plot",
        "split_plot",
        "full_plot",
        "single_years_plot",
        "single_plot_insitu",
        "extrema_plot",
        "qa_boxplot",
        "yearly_heatmap_pixel",
    })
    

    def plot_background_pts(self, ax, latitude_idx, longitude_idx, masked_values, masked_time, aggregation = False, **style_kwargs):
        base_style = {"alpha":0.3, "color":"grey","s":10}
        style = {**base_style,**style_kwargs}

        if aggregation:
            fn = f.construct_background_agg_filename(self.lakeID)
            fp = f.construct_background_agg_dirpath(self.data_folder,self.variable)
            time, values = f.read_background_agg_pixel(fp/fn, latitude_idx, longitude_idx)

            x = f.datenum_to_datetime(time)
            y = values
        else:
            x = f.datenum_to_datetime(masked_time)
            y = masked_values

        sc = ax.scatter(x, y, label="Data", **style)
        return sc


    def plot_data_gaps(self, ax, pixel_data):
        gap_starts = pixel_data["gap_starts"]
        gap_ends   = pixel_data["gap_ends"]
        for gs, ge in zip(gap_starts, gap_ends):
            ax.axvspan(gs, ge, color="orange", alpha=0.15, zorder=0)
        if len(gap_starts) > 0:
            ax.axvspan(gap_starts[0], gap_ends[0], color="orange", alpha=0.15, zorder=0, label="Data gap")

            
    def plot_pheno_metrics(self, ax, plotting_data, spline_x, spline_y, time_frame, variables = None):

        if variables is None:
            variables = ["pks", "trgs", "midUP", "midDOWN"]

        ax.plot(f.datenum_to_datetime(spline_x), spline_y, color="black", linewidth=1, label="Spline")

        used_labels = set()
        qa_colors = {0: "blue", 1: "orange", 2: "red"}
        qa_labels = {0: "Good", 1: "Fair", 2: "Poor"}
        qa_vars = {'pks','trgs'}
        for var in variables:
            if var in qa_vars:
                qa_vals = plotting_data[var][2]

                for qa in self.QA_LEVELS:
                    pm =  qa_vals == qa
                    if not pm.any():
                        continue
                    label = qa_labels[qa]
                    ax.scatter(plotting_data[var][0][pm], plotting_data[var][1][pm], color=qa_colors[qa], s=50,
                                    marker="o", edgecolors="black", linewidths=0.5,
                                    zorder=4, label=label if label not in used_labels else None)
                    used_labels.add(label)

            else:
                ax.scatter(plotting_data[var][0], plotting_data[var][1], s=30, zorder=4, **self.METRIC_CONFIG[var]['style'])
                used_labels.add(label)
                # labeling not working right now

        neg_values_sub = f.mark_negative_values_timeseries_plot(ax,plotting_data,time_frame)
        return neg_values_sub


    def yearly_cubic_spline(self, ax, latitude_idx, longitude_idx, years=None):
        """Overlay csaps splines for multiple years on a common fractional-month x-axis.

        Fits a single spline over the full valid time series for the pixel at
        (latitude_idx, longitude_idx), then slices it year by year and plots each slice
        against fractional month (1–12) using a distinct colour from the
        cc.glasbey_light palette. Detected peaks and troughs are overlaid as
        scatter markers coloured by QA level.

        A vertical colorbar maps year indices to their assigned colours, and two
        separate legends show QA marker styles (Good / Fair / Poor) and event
        types (Peak / Trough).

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes on which to draw the overlaid splines and markers.
        latitude_idx : int
            Row (lat) index of the pixel.
        longitude_idx : int
            Column (lon) index of the pixel.
        years : list of str, optional
            Calendar years to include. Each element must be a string (e.g. '2005').
            Defaults to 2002–2024. Years with fewer than 3 spline points are
            skipped with a warning.

        Returns
        -------
        None
        """
        if years is None:
            years = range(2002, 2025)

        g = self._load_extracted_globals()
        px = self._load_pixel_data(latitude_idx, longitude_idx)
        lat, lon, t_all = g["lat"], g["lon"], g["t_all"]
        smoothing = px["smoothing"]
        lat_val = round(float(lat[latitude_idx]), 4)
        lon_val = round(float(lon[longitude_idx]), 4)

        # QA legend (markers only, no year lines)
        qa_markers = {0: "o", 1: ".", 2: "x"}
        qa_labels = {0: "Good", 1: "Fair", 2: "Poor"}
        qa_handles = [
            mlines.Line2D(
                [], [],
                color="black",
                marker=qa_markers[qa],
                linestyle="None",
                markersize=6,
                label=qa_labels[qa]
            )
            for qa in self.QA_LEVELS
        ]
        
        cmap = cc.cm.rainbow
        year_colors = {year: cmap(i / max(len(years) - 1, 1)) for i, year in enumerate(years)}
        cmap_discrete = ListedColormap([year_colors[y] for y in years])

        # Fit one spline over all valid data
        mask_all = (px["values"] != -9999) & (px["qa"] == 0)
        values_m_all = px["values"][mask_all]
        time_m_all = t_all[mask_all]
        smooth_x_all, smooth_y_all = f.calculate_spline(
            whole_timeframe=t_all, masked_values=values_m_all,
            masked_time=time_m_all, smoothing_parameter=smoothing
        )
        if smooth_x_all is None:
            warnings.warn("No data to plot")
            return
        smooth_dates_all = np.array(f.datenum_to_datetime(smooth_x_all))

        for year in years:
            plotting_data = f.grab_plotting_variables(
                start=year,
                end=year,
                pixel_data=px,
                variables=["pks", "trgs"]
            )
            pks_x_sub, pks_y_sub, pks_qa_sub   = plotting_data["pks"]
            trgs_x_sub, trgs_y_sub, trgs_qa_sub = plotting_data["trgs"]

            # Subset the pre-fitted spline to this year and convert to fractional month (1–12)
            mask_year = np.array([d.year == year for d in smooth_dates_all])
            smooth_dates_year = smooth_dates_all[mask_year]
            smooth_y  = smooth_y_all[mask_year]

            if len(smooth_dates_year) > 1:
                smooth_x_month = f.to_frac_month(smooth_dates_year)
                ax.plot(smooth_x_month, smooth_y, color=year_colors[year], linewidth=1, label=str(year))

                for qa in self.QA_LEVELS:
                    pm = pks_qa_sub == qa
                    tm = trgs_qa_sub == qa
                    if pm.any():
                        ax.scatter(f.to_frac_month(pks_x_sub[pm]), pks_y_sub[pm], color=year_colors[year], s=50,
                            marker=qa_markers[qa], edgecolors="black", linewidths=0.5,
                            zorder=4, label= qa_labels[qa] if year == years[0] else None)
                    if tm.any():
                        ax.scatter(f.to_frac_month(trgs_x_sub[tm]), trgs_y_sub[tm], color=year_colors[year], s=50,
                            marker=qa_markers[qa], edgecolors="black", linewidths=0.5,
                            zorder=4, label=qa_labels[qa] if (year == years[0] and not pm.any()) else None)
            else:
                warnings.warn(f"Not enough data to plot for year {year}")

        # One-time axis setup after all years are plotted
        int_years = [int(y) for y in years]
        # derive months from data
        all_months = sorted(set(d.month for d in smooth_dates_all if d.year in int_years))
        month_names = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
        ax.set_xlim(min(all_months), max(all_months))
        ax.set_xticks(all_months)
        ax.set_xticklabels([month_names[m - 1] for m in all_months])
        all_pks_y = px["pks_y"][np.array([d.year in int_years for d in px["pks_x"]])]
        all_trgs_y = px["trgs_y"][np.array([d.year in int_years for d in px["trgs_x"]])]
        if len(all_pks_y) > 0 and len(all_trgs_y) > 0:
            pks_sorted = sorted(all_pks_y)
            ymax = pks_sorted[-2] if pks_sorted[-1] > 10 and len(pks_sorted) > 1 else pks_sorted[-1]
            ax.set_ylim(sorted(all_trgs_y)[0] - 0.5, ymax + 0.5)
        textstr = f"{self.variable}\n Lake ID:{self.lakeID}\n lat, lon: {lat_val}, {lon_val}"
        ax.set_title(textstr)
        ax.grid(axis="x", linewidth=0.5)
        ax.grid(axis="y")
        ax.set_ylabel("[ug/L]")

        # Color legend for pks and trgs
        type_handles = [
            mlines.Line2D([], [], color="black", marker="o",
                        linestyle="None", markersize=8, label="Peak"),
            mlines.Line2D([], [], color="darkgray", marker="o",
                        linestyle="None", markersize=8, label="Trough"),
        ]

        leg1 = ax.legend(handles=qa_handles, loc="upper left")
        ax.add_artist(leg1)
        ax.legend(handles=type_handles, loc="upper right")

        # Colorbar for year colors
        norm = mcolors.BoundaryNorm(boundaries=range(len(years) + 1), ncolors=len(years))
        sm   = plt.cm.ScalarMappable(cmap=cmap_discrete, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, orientation="vertical", pad=0.01, aspect=30)
        cbar.set_ticks([i + 0.5 for i in range(len(years))])
        cbar.set_ticklabels(years)


    def annotations_and_limits(self, ax, plotting_data, metrics_dict, time_frame, lat_val, lon_val, neg_values_sub, annotation = None):
        start, end = time_frame[0], time_frame[1]

        ax.legend(loc="upper left", ncol= 2)
        textstr = f"{self.variable}\n Lake ID:{self.lakeID}\n lat, lon: {lat_val:.4f}, {lon_val:.4f}\n Total RMSE, R$^2$, MAD: {metrics_dict['rmse'][1]:.4f}, {metrics_dict['r2'][1]:.4f}, {metrics_dict['mad'][1]:.4f}"
        ax.set_title(textstr)
        ax.xaxis.set_minor_locator(mdates.YearLocator())
        ax.grid(axis="x", which="minor", linewidth=0.5)
        ax.grid(axis="x", which="major", linewidth=0.5)
        ax.grid(axis="y")
        ax.set_ylabel("[ug/L]")
        ax.set_xlim(pd.to_datetime('01-01-' + str(start), format='%d-%m-%Y') , pd.to_datetime('31-12-' + str(end), format='%d-%m-%Y'))
        trgs_y_sub = plotting_data["trgs"][1]

        pks_lim_sub = sorted(plotting_data["pks"][1])
        # if max(pks_lim_sub)> 10:
        #         ax.set_ylim(sorted(trgs_y_sub)[0]-0.5, pks_lim_sub[-2]+0.5)
        # else:
        #         ax.set_ylim(sorted(trgs_y_sub)[0]-0.5, pks_lim_sub[-1]+0.5)

        trgs_lim_sub = sorted(plotting_data["trgs"][1])

        if len(pks_lim_sub) > 0 and len(trgs_lim_sub) > 0:
            ymax = pks_lim_sub[-1]
            if ymax > 10 and len(pks_lim_sub) > 1:
                ymax = pks_lim_sub[-2]
            ymin = trgs_lim_sub[0]

            ax.set_ylim(ymin - 0.5, ymax + 0.5)

        else:
            warnings.warn(
                f"No peaks/troughs available for {start}-{end}; using automatic y-limits."
            )
        if not annotation:
            if neg_values_sub:
                ax.text(0.99,0.99,f"# Neg.values: {sum(neg_values_sub)} \n RMSE: {metrics_dict['rmse'][0]:.3f}\n R$^2$: {metrics_dict['r2'][0]:.3f}\n MAD: {metrics_dict['mad'][0]:.0f}", transform = ax.transAxes,   ha= "right", va= "top", zorder = 10)
            else:
                ax.text(0.99,0.99,f"RMSE:{metrics_dict['rmse'][0]:.3f} \n R$^2$: {metrics_dict['r2'][0]:.3f}\n MAD: {metrics_dict['mad'][0]:.3f}", transform = ax.transAxes,   ha= "right", va= "top", zorder = 10)
        else:
            lines = []
            if "R2" in annotation:
                lines.append(f"R$^2$: {round(metrics_dict['r2'][0], 3)}")
            if "RMSE" in annotation:
                lines.append(f"RMSE: {round(metrics_dict['rmse'][0], 3)}")
            if "MAD" in annotation:
                lines.append(f"MAD: {round(metrics_dict['mad'][0], 3)}")
            if "neg" in annotation and neg_values_sub:
                lines.append(f"# Neg.values: {sum(neg_values_sub)}")
            if lines:
                ax.text(0.99, 0.99, "\n".join(lines),
                    transform=ax.transAxes, ha="right", va="top", zorder = 10)


    def single_plot(self, latitude_idx, longitude_idx, ax, aggregation = False, start= 0, end= 9999, annotation = None,variables = None, check_buffer = False):
        """Plot raw observations, the smoothed spline, and all phenological events for a pixel.

        Displays a scatter of valid (QA==0) observations or 3×3 aggregated values,
        overlaid with the csaps spline and scatter markers for peaks, troughs,
        green-up midpoints, and green-down midpoints. Negative values are flagged
        with red crosses and trigger a warning. Fit metrics (RMSE, R², MAD) are
        annotated on the plot.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel.
        longitude_idx : int
            Column (lon) index of the pixel.
        ax : matplotlib.axes.Axes
            Axes on which to draw the plot.
        aggregation : bool, optional
            If True, replace raw scatter with the 3×3 neighbourhood spatial median.
            Requires spatial_aggregation() to have been called or available cache.
        start : int, optional
            First year to display (inclusive). 0 = earliest in the series.
        end : int, optional
            Last year to display (inclusive). 9999 = latest in the series.
        check_buffer : bool, optional
            If True, warn when the requested pixel falls outside the lake's 1 km
            inward-shrunk boundary (self.prepped_geom). Default False - the pixel
            is not checked against the buffer.

        Returns
        -------
        None
        """
        g  = self._load_extracted_globals()
        pixel_data = self._load_pixel_data(latitude_idx, longitude_idx)
        lat, lon, t_all = g["lat"], g["lon"], g["t_all"]
        lat_val = float(lat[latitude_idx])
        lon_val = float(lon[longitude_idx])
        smoothing = pixel_data["smoothing"]

        if check_buffer and not self.prepped_geom.contains(Point(lon_val, lat_val)):
            warnings.warn(
                f"Pixel ({latitude_idx}, {longitude_idx}) is outside the 1 km lake buffer "
                f"for lake ID {self.lakeID}."
            )

        plotting_data = f.grab_plotting_variables(start = start, end = end, pixel_data=pixel_data, variables=variables)

        mask     = (pixel_data["values"] != -9999) & (pixel_data["qa"] == 0)
        values_m = pixel_data["values"][mask]
        time_m   = t_all[mask]

        if len(values_m) == 0:
            warnings.warn("No data to plot")
            return

        smooth_x, smooth_y = f.calculate_spline(whole_timeframe= t_all, masked_values=values_m, masked_time= time_m, smoothing_parameter=smoothing)

        metrics_dict, plot_time_frame = f.calculate_metrics_to_plot(start = start, end = end, masked_values= values_m, masked_time=time_m, smoothing_parameter=smoothing)

        if metrics_dict is None:
            return

        self.plot_background_pts(ax = ax, latitude_idx= latitude_idx, longitude_idx = longitude_idx, masked_values=values_m, masked_time=time_m, aggregation=aggregation)
        self.plot_data_gaps(ax = ax, pixel_data = pixel_data)
        neg_values_sub = f.plot_variables(ax = ax, plotting_data= plotting_data, spline_x= smooth_x, spline_y= smooth_y, time_frame= plot_time_frame, variables= variables)
        self.annotations_and_limits(ax = ax, plotting_data= plotting_data, metrics_dict= metrics_dict, time_frame=plot_time_frame, lat_val = lat_val, lon_val = lon_val, neg_values_sub=neg_values_sub, annotation = annotation)


    def split_plot(self, latitude_idx, longitude_idx, ax0, ax1, aggregation = False, start0= 0, end0= 9999, start1= 0, end1=9999):
        """Plot two year-windowed single_plots side by side for the same pixel.

        Calls single_plot twice — once on ax0 with [start0, end0] and once on ax1
        with [start1, end1]. Intended for comparing two non-overlapping time periods
        at the same pixel.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel.
        longitude_idx : int
            Column (lon) index of the pixel.
        ax0 : matplotlib.axes.Axes
            Axes for the first time window.
        ax1 : matplotlib.axes.Axes
            Axes for the second time window.
        aggregation : bool, optional
            If True, use 3×3 spatial median scatter on both panels.
        start0 : int, optional
            First year of the first window. 0 = earliest in the series.
        end0 : int, optional
            Last year of the first window. 9999 = latest in the series.
        start1 : int, optional
            First year of the second window. 0 = earliest in the series.
        end1 : int, optional
            Last year of the second window. 9999 = latest in the series.

        Returns
        -------
        None
        """
        if (start0 == 0) and (start1 == 0) and (end0 == 9999) and (end1 == 9999):
            warnings.warn("split_plot needs a least end0 and start1 parameter, otherwise use full_plot")

        self.single_plot(latitude_idx = latitude_idx, longitude_idx= longitude_idx, ax=ax0, aggregation = aggregation, start= start0, end = end0)
        self.single_plot(latitude_idx = latitude_idx, longitude_idx= longitude_idx, ax=ax1, aggregation = aggregation, start=start1, end=end1)


    def full_plot(self, latitude_idx, longitude_idx, ax, aggregation = False):
        """Plot the complete valid time series for a pixel.

        Auto-detects the first and last years with valid observations and
        delegates to single_plot with those bounds.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel.
        longitude_idx : int
            Column (lon) index of the pixel.
        ax : matplotlib.axes.Axes
            Axes on which to draw the plot.
        aggregation : bool, optional
            If True, use 3×3 neighbourhood spatial median scatter.

        Returns
        -------
        None
        """
        with netCDF4.Dataset(self.e_path) as nc:
            t_all = f.unix_to_datenum(nc.variables["time"])
            variable = getattr(nc, "variable")
            values = np.array(nc.variables[variable][:, latitude_idx, longitude_idx])
            mask = (values != -9999) & (np.array(nc.variables[getattr(nc, 'qa')][:, latitude_idx, longitude_idx]) == 0)
            values_m = values[mask]
            time_m = t_all[mask]

        if len(values_m) > 1:
            limits = sorted(f.datenum_to_datetime(time_m))
            full_plot_start = min(limits).year
            full_plot_end = max(limits).year
            self.single_plot(latitude_idx= latitude_idx, longitude_idx= longitude_idx, ax = ax, aggregation = aggregation, start=full_plot_start, end= full_plot_end)

        else:
            warnings.warn("No data to plot")


    def single_years_plot(self, latitude_idx, longitude_idx, years, ncol, nrow, annotation, ylim=None):
        """Plot one panel per year in a grid, each showing phenology for a single pixel.

        Creates a figure with ``nrow × ncol`` subplots. Each subplot calls
        :meth:`single_plot` for one year in *years*, with enlarged markers and
        month-number x-axis labels. Panels beyond ``len(years)`` are hidden.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel.
        longitude_idx : int
            Column (lon) index of the pixel.
        years : list of int
            Calendar years to display, one per panel.
        ncol : int
            Number of subplot columns.
        nrow : int
            Number of subplot rows. Must satisfy ``ncol * nrow >= len(years)``.
        annotation : list of str or None
            Passed to :meth:`single_plot`. Controls which fit metrics are shown
            (e.g. ``["R2", "RMSE", "MAD"]``). ``None`` shows all metrics.
        ylim : tuple of (float, float) or None, optional
            If provided, sets the y-axis limits as ``(bottom, top)`` for every
            panel, overriding the automatic limits set by :meth:`single_plot`.

        Returns
        -------
        matplotlib.figure.Figure
        """
        _MARKER_SIZES = {"Data": 50, "Peaks": 150, "Troughs": 150, "Mid Up": 150, "Mid Down": 150}

        fig, axs = plt.subplots(nrow, ncol, constrained_layout=True, squeeze=False, figsize=(ncol * 5, nrow * 4))
        for year, ax in zip(years, axs.flatten()):
            self.single_plot(latitude_idx, longitude_idx, ax, start=year, end=year, annotation=annotation)

            for col in ax.collections:
                if col.get_label() in _MARKER_SIZES:
                    col.set_sizes([_MARKER_SIZES[col.get_label()]])

            if ax.texts:
                ax.texts[-1].set_fontsize(15)

            legend = ax.get_legend()
            if legend is not None:
                legend.remove()
            if ylim is not None:
                ax.set_ylim(ylim)
            else:
                ax.set_ylim(bottom=-0.5,top = ax.get_ylim()[1]*1.5)
            ax.set_title(str(year), fontsize=20)
            ax.set_ylabel("[ug/L]", fontsize=15)
            ax.xaxis.set_major_locator(mdates.MonthLocator())
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%#m'))
            ax.tick_params(labelsize=15)

        for ax in axs.flatten()[len(years):]:
            ax.set_visible(False)


    def single_plot_insitu(self, latitude_idx,longitude_idx,ax,insitu_df,aggregation=False, start = 0, end = 9999, insitu_date_col="datetime", insitu_value_col="chlorophyll_a", insitu_station_col=None, station_id=None, max_depth = 5):
        """
        Plot satellite observations + spline + phenology + in situ overlay.
        """

        # -------------------------------------------------
        # FIRST: draw the original plot
        # -------------------------------------------------

        self.single_plot(
            latitude_idx=latitude_idx,
            longitude_idx=longitude_idx,
            ax=ax,
            aggregation=aggregation,
            start=start,
            end=end)

        # -------------------------------------------------
        # PREPARE IN SITU DATA
        # -------------------------------------------------

        insitu_mean = f.prep_dimark_data(insitu_df = insitu_df,start=start, end=end,  insitu_date_col=insitu_date_col, insitu_value_col=insitu_value_col, insitu_station_col=insitu_station_col, station_id=station_id, max_depth = max_depth)
        
        # -------------------------------------------------
        # OVERLAY IN SITU DATA
        # -------------------------------------------------

        ax.scatter(
            insitu_mean[insitu_date_col],
            insitu_mean[insitu_value_col],
            color="red",
            marker="D",
            s=20,
            edgecolor="black",
            linewidth=0.5,
            zorder=5,
            alpha= 0.5,
            label=f"In Situ (<{max_depth}m)"
        )

        # optional connecting line
        # ax.plot(
        #         insitu[insitu_date_col],
        #         insitu[insitu_value_col],
        #         color="red",
        #         alpha=0.5,
        #         linewidth=1
        # )

        # -------------------------------------------------
        # UPDATE LEGEND
        # -------------------------------------------------

        handles, labels = ax.get_legend_handles_labels()

        by_label = dict(zip(labels, handles))

        ax.legend(
            by_label.values(),
            by_label.keys(),
            loc="upper left",
            ncol=2
        )


    def extrema_plot(self, latitude_idx, longitude_idx, ax,  peak = True, aggregation= False,
                     start = 0, end = 9999, background_pts = True, purple_chla21= False, show_legend = True):
        """Plot detected peaks or troughs as a stem plot with optional background scatter.

        Displays summer peaks or winter troughs for the pixel at (latitude_idx, longitude_idx)
        as vertical stems, with each extremum marker shaped by its QA flag
        while keeping the product's colour. Background observations may be shown as
        a raw scatter (QA==0 only) or 3×3 spatial median (when aggregation=True). Negative
        values are flagged with red crosses and trigger a warning.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel.
        longitude_idx : int
            Column (lon) index of the pixel.
        ax : matplotlib.axes.Axes
            Axes on which to draw the plot.
        peak : bool, optional
            If True (default), plot summer peaks. If False, plot troughs.
        aggregation : bool, optional
            If True, show the 3×3 neighbourhood median instead of raw scatter.
            Requires spatial_aggregation() to have been called or available cache.
        start : int, optional
            First year to display (inclusive). 0 = earliest in the series.
        end : int, optional
            Last year to display (inclusive). 9999 = latest in the series.
        background_pts : bool, optional
            If True (default), show scatter background observations.
            Cannot be False when aggregation is True.
        purple_chla21 : bool, optional
            Unused colour-override flag kept for API compatibility.

        Returns
        -------
        float or None
            The upper y-axis limit set for the plot, or None if there is no data.

        Raises
        ------
        ValueError
            If background_pts is False and aggregation is True simultaneously.
        """
        g = self._load_extracted_globals()
        px = self._load_pixel_data(latitude_idx, longitude_idx)
        lat, lon, t_all = g["lat"], g["lon"], g["t_all"]
        smoothing = px["smoothing"]
        lat_val = float(lat[latitude_idx])
        lon_val = float(lon[longitude_idx])

        var = "pks" if peak else "trgs"
        plotting_data = f.grab_plotting_variables(start=start, end=end, pixel_data=px, variables=[var])
        x_sub, y_sub, qa_sub = plotting_data[var]

        mask     = (px["values"] != -9999) & (px["qa"] == 0)
        values_m = px["values"][mask]
        time_m   = t_all[mask]

        
        if len(values_m) <= 1:
            warnings.warn("No data to plot (check valid indices)")
            return None

        limits = sorted(f.datenum_to_datetime(time_m))
        function_start = min(limits).year if start <= min(limits).year else start
        function_end = max(limits).year if end >= max(limits).year else end

        phenology_name = self.variable
        
        var_cfg = self.get_plot_config(
            "var",
            phenology_name,
            use_alt=purple_chla21
        )
        
        var_label = var_cfg["label"]
        var_style = var_cfg["style"]
        background_style = {"alpha": 0.1, "s": 10}
    
        if background_pts:
            self.plot_background_pts(
                ax=ax,
                latitude_idx = latitude_idx,
                longitude_idx = longitude_idx,
                masked_values = values_m,
                masked_time = time_m,
                aggregation= aggregation,
                **{**var_style, **background_style}
            )
        elif not background_pts and aggregation:
            warnings.warn(f"Aggregation ignored for lake ID {self.lakeID} since backrgound_pts turned off.")

        self.plot_data_gaps(ax=ax, pixel_data=px)

        ax.stem(x_sub, y_sub, markerfmt=" ", basefmt = " ",linefmt = var_style['color'])
        
        seen_labels = set()
        for q in self.QA_LEVELS:
            qm = qa_sub == q
            if not qm.any():
                continue

            qa_cfg = self.get_plot_config("qa", q, use_alt = True)
            qa_style = qa_cfg["style"]

            qa_label = qa_cfg["label"]
            if qa_label in seen_labels:
                qa_label = None
            else:
                seen_labels.add(qa_label)

            combined_style = {
                **var_style,            # color
                "marker": qa_style["marker"],
                "s": 50,
                "edgecolors": var_style['color'],
                "linewidths": 2,
                "zorder": 4,
            }

            ax.scatter(
                x_sub[qm],
                y_sub[qm],
                label=qa_label,
                **combined_style,
            )

        if (y_sub < 0).any():
            mask =  y_sub<0
            ax.scatter(x_sub[mask], y_sub[mask], color="red", s=50, marker="x", zorder=6, label="Negative value")
            warnings.warn(f"Negative Peak(s) in time period {start}-{end}", Warning)

        if show_legend:
            ax.legend(loc="upper left", ncol= 2)
        
        plot_var = "Peak" if peak else "Trough"
        textstr = f"{plot_var} Comparison\n Lake ID:{self.lakeID}\n lat, lon: {lat_val:.4f}, {lon_val:.4f}"

        ax.set_title(textstr)
        ax.xaxis.set_minor_locator(mdates.YearLocator())
        ax.grid(axis="x", which="minor", linewidth=0.5)
        ax.grid(axis="x", which="major", linewidth=0.5)
        ax.grid(axis="y", linewidth=0.5)
        ax.set_ylabel("[ug/L]")

        ax.set_xlim(
            pd.to_datetime('01-01-' + str(function_start), format='%d-%m-%Y'),
            pd.to_datetime('31-12-' + str(function_end), format='%d-%m-%Y')
        )
        pks_lim_sub = sorted(y_sub)
        if max(pks_lim_sub)> 10:
            ymax = pks_lim_sub[-2]+0.5
            ax.set_ylim(-0.5, ymax)
        else:
            ymax = pks_lim_sub[-1]+0.5
            ax.set_ylim(-0.5, ymax)
        return ymax


    def single_plot_background(self, latitude_idx, longitude_idx, ax, fig, aggregation = False, start= 0, end= 9999):
        """Plot a pixel time series with QA-coloured scatter, spline, and phenological events.

        Like single_plot, but colours each background scatter point by its QA flag
        using a discrete colormap and adds a QA colorbar to the figure. Peaks,
        troughs, green-up and green-down midpoints are overlaid as scatter markers.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel.
        longitude_idx : int
            Column (lon) index of the pixel.
        ax : matplotlib.axes.Axes
            Axes on which to draw the plot.
        fig : matplotlib.figure.Figure
            Figure used to attach the QA colorbar.
        aggregation : bool, optional
            If True, replace raw scatter with the 3×3 neighbourhood median.
            Requires spatial_aggregation() to have been called or available cache.
        start : int, optional
            First year to display (inclusive). 0 = earliest in the series.
        end : int, optional
            Last year to display (inclusive). 9999 = latest in the series.

        Returns
        -------
        None
        """
        g  = self._load_extracted_globals()
        pixel_data = self._load_pixel_data(latitude_idx, longitude_idx)
        lat, lon, t_all = g["lat"], g["lon"], g["t_all"]
        smoothing = pixel_data["smoothing"]

        plotting_data = f.grab_plotting_variables(start=start, end=end, pixel_data=pixel_data)

        # No QA==0 filter here — all non-fill values kept for QA-coloured scatter
        mask     = (pixel_data["values"] != -9999)
        values_m = pixel_data["values"][mask]
        time_m   = t_all[mask]
        qa_mask  = pixel_data["qa"][mask]

        if len(values_m) == 0:
            warnings.warn("No data to plot")
            return

        smooth_x, smooth_y = f.calculate_spline(
            whole_timeframe=t_all, masked_values=values_m,
            masked_time=time_m, smoothing_parameter=smoothing
        )
        if smooth_x is None:
            warnings.warn("No data to plot")
            return

        metrics_dict, plot_time_frame = f.calculate_metrics_to_plot(
            start=start, end=end, masked_values=values_m,
            masked_time=time_m, smoothing_parameter=smoothing
        )
        if metrics_dict is None:
            return

        # QA-coloured background scatter (unique to this method)
        qa_unique = sorted(np.unique(qa_mask))
        cmap = plt.cm.get_cmap("tab10", len(qa_unique))
        cmap_new = ListedColormap(cmap(np.arange(len(qa_unique))))
        qa_to_idx = {qa: i for i, qa in enumerate(qa_unique)}
        qa_idx = np.array([qa_to_idx[q] for q in qa_mask])
        bounds = np.arange(-0.5, len(qa_unique) + 0.5, 1)

        norm = BoundaryNorm(bounds, cmap_new.N)

        sc = self.plot_background_pts(ax, latitude_idx, longitude_idx, values_m, time_m, aggregation = aggregation)

        neg_values_sub = f.plot_variables(
            ax=ax, plotting_data=plotting_data, spline_x=smooth_x, spline_y=smooth_y,
            time_frame=plot_time_frame
        )

        self.annotations_and_limits(
            ax=ax, plotting_data=plotting_data, metrics_dict=metrics_dict,
            time_frame=plot_time_frame, lat=lat, lon=lon,
            latitude_idx=latitude_idx, longitude_idx=longitude_idx,
            neg_values_sub=neg_values_sub
        )

        cbar = fig.colorbar(sc, ax=ax, boundaries=bounds)
        cbar.set_label("QA indicators")
        cbar.set_ticks(np.arange(len(qa_unique)))
        cbar.set_ticklabels([str(q) for q in qa_unique])

