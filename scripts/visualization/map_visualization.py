import pandas as pd
import netCDF4
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from matplotlib.colors import ListedColormap
import functions as f
from shapely.prepared import prep
from shapely.geometry import Point
from visualization.base import PhenologyVisualizationBase

class MapVisualization(PhenologyVisualizationBase):
    """Spatial and temporal map graphs."""

    graph_type = "map"
    graph_methods = frozenset({
        "interactive_pixel_map",
        "metric_map",
        "interactive_metric_map",
        "time_map",
        "time_map_panel",
        "single_day_map",
        "yearly_heatmap_lake",
    })

    
    def pixel_map(self, latitude_idx, longitude_idx, ax):
        """Plot a grayscale coverage map with the selected pixel marked.

        Reads the summary grid from the extract NetCDF and displays it in
        grayscale, masking cells with fewer than 2 valid observations. The
        selected pixel is overlaid as a red star.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel to highlight.
        longitude_idx : int
            Column (lon) index of the pixel to highlight.
        ax : matplotlib.axes.Axes
            Axes on which to draw the map.

        Returns
        -------
        matplotlib.image.AxesImage
            The imshow image object.
        """
        with netCDF4.Dataset(self.e_path) as nc:
            summary = np.array(nc.variables["summary"][:, :])

        # mask invalid cells
        masked_summary = np.ma.masked_where(summary <= 2, summary)

        cmap = ListedColormap("gray")

        # plot as geographic map
        im = ax.imshow(masked_summary, cmap=cmap, aspect="auto", origin="lower")

        # plot selected pixel
        ax.plot(longitude_idx, latitude_idx, "r*", markersize=14, zorder=5, label="Pixel")

        ax.set_xlabel("Lon index")
        ax.set_ylabel("Lat index")
        text_str = f"Pixel Location\n Lake ID:{self.lakeID}"
        ax.set_title(text_str)
        ax.legend()

        return im


    def interactive_pixel_map(self, ax):
        """Plot a clickable coverage map for pixel selection.

        Displays the summary coverage grid in grayscale. Clicking a valid cell
        prints its (lat_idx, lon_idx) and marks it with a dot and text label.
        Requires an interactive matplotlib backend (e.g. %matplotlib widget).

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes on which to draw the map.

        Returns
        -------
        int
            Matplotlib canvas connection ID, which can be passed to
            ax.figure.canvas.mpl_disconnect() to remove the click handler.
        """
        with netCDF4.Dataset(self.e_path) as nc:
            summary = np.array(nc.variables["summary"][:, :])

        masked_summary = np.ma.masked_where(summary <= 2, summary)

        cmap = ListedColormap("gray")

        ax.imshow(masked_summary, cmap=cmap, aspect="auto", origin="lower")

        ax.set_xlabel("Lon index")
        ax.set_ylabel("Lat index")
        ax.set_title(f"Pixel Location\n Lake ID: {self.lakeID}")

        def on_click(event):
            """Label the clicked valid pixel with its lat/lon indices."""
            print("click detected")

            if event.inaxes is not ax or event.xdata is None or event.ydata is None:
                    return

            lon_idx = int(round(event.xdata))
            lat_idx = int(round(event.ydata))

            lat_idx = max(0, min(lat_idx, summary.shape[0] - 1))
            lon_idx = max(0, min(lon_idx, summary.shape[1] - 1))

            if summary[lat_idx, lon_idx] <= 2:
                    print("invalid cell")
                    return

            ax.text(
            lon_idx, lat_idx,
            f"{lat_idx},{lon_idx}",
            color="black",
            fontsize=10,
            ha="center", va="center"
            )

            ax.plot(lon_idx, lat_idx, "ro", markersize=6)
            ax.figure.canvas.draw_idle()

        cid = ax.figure.canvas.mpl_connect("button_press_event", on_click)

        return cid


    def metric_map(self, metric_scores:dict, metric_str:str, fig, ax, colormap = None, colorbar_extent= None):
        """Plot a spatial heatmap of a precomputed per-pixel metric.

        Pixels outside the 1 km-inset lake boundary are masked. The lake outline
        is overlaid from the shapefile geometry.

        Parameters
        ----------
        metric_scores : dict
            Mapping of (i, j) grid index tuples to metric float values, as
            returned by r2_scores, MAD_scores, RMSE_scores, etc.
        metric_str : str
            Label used for the colorbar and plot title (e.g. 'R$^2$', 'RMSE').
        fig : matplotlib.figure.Figure
            Figure used to attach the colorbar.
        ax : matplotlib.axes.Axes
            Axes on which to draw the heatmap.
        colormap : str or matplotlib.colors.Colormap, optional
            Colormap passed to imshow. Defaults to 'RdYlBu'.
        colorbar_extent : list of float, optional
            [vmin, vmax] for the colorbar. Defaults to [0, 1].

        Returns
        -------
        matplotlib.image.AxesImage
            The imshow image object.

        Raises
        ------
        ValueError
            If the lake ID derived from p_path is not found in the shapefile.
        """
        lake_id = int(self.lakeID)
        lake_row = self.gdf[self.gdf["id"] == lake_id]
        if lake_row.empty:
            raise ValueError(f"Lake ID {lake_id} not found in shapefile.")
        
        geom = lake_row.geometry.iloc[0]
        buffered_geom = self.shrink_geometry(geom)
        buffered_geom_prepared = prep(buffered_geom)

        map_data, extent = f.grab_metrics(self.e_path, metric_scores, buffered_geom_prepared)

        im = f.plot_map_data(colormap, map_data, extent, ax, cmap_extent=colorbar_extent)
        f.plot_lake_outline(geometry=geom, ax=ax)
        f.set_labels(ax, fig, im,
            title=f"{metric_str}-Scores for Lake: ID {lake_id}",
            colorbar_label=metric_str)
        return im


    def interactive_metric_map(self, metric_scores, metric_str:str, fig, ax):
        """Plot a clickable spatial heatmap of a precomputed per-pixel metric.

        Like metric_map but with click interaction: clicking a pixel prints and
        labels its (lat_idx, lon_idx). Colorbar range is fixed to [0, 1].
        Requires an interactive matplotlib backend (e.g. %matplotlib widget).

        Parameters
        ----------
        metric_scores : dict
            Mapping of (i, j) grid index tuples to metric float values.
        metric_str : str
            Label used for the colorbar and plot title.
        fig : matplotlib.figure.Figure
            Figure used to attach the colorbar.
        ax : matplotlib.axes.Axes
            Axes on which to draw the heatmap.

        Returns
        -------
        int
            Matplotlib canvas connection ID for removing the click handler.

        Raises
        ------
        ValueError
            If the lake ID derived from p_path is not found in the shapefile.
        """
        lake_id = int(self.lakeID)

        lake_row = self.gdf[self.gdf["id"] == lake_id]
        if lake_row.empty:
                raise ValueError(f"Lake ID {lake_id} not found in shapefile.")
        geom = lake_row.geometry.iloc[0]
        buffered_geom = self.shrink_geometry(geom)
        buffered_geom_prepared = prep(buffered_geom)

        map_data, extent = f.grab_metrics(self.e_path, metric_scores, buffered_geom_prepared)
        g = self._load_extracted_globals()
        lats = g["lat"]
        lons = g["lon"]

        im = f.plot_map_data(None, map_data, extent, ax, cmap_extent=[0, 1])
        f.plot_lake_outline(geometry=geom, ax=ax)
        f.set_labels(ax, fig, im,
            title=f"{metric_str}-Scores for Lake: ID {lake_id}",
            colorbar_label=metric_str)

        def on_click(event):
            """Label the clicked metric-map pixel with its lat/lon indices."""
            if event.inaxes is not ax or event.xdata is None or event.ydata is None:
                return

            clicked_lon = event.xdata
            clicked_lat = event.ydata

            lon_idx = int(np.abs(lons - clicked_lon).argmin())
            lat_idx = int(np.abs(lats - clicked_lat).argmin())

            ax.text(
            lons[lon_idx], lats[lat_idx],
            f"{lat_idx},{lon_idx}",
            color="black",
            fontsize=10,
            ha="center", va="center"
            )
            ax.plot(lons[lon_idx], lats[lat_idx], "ro", markersize=6)
            fig.canvas.draw_idle()

        cid = fig.canvas.mpl_connect("button_press_event", on_click)
        return cid
    

    def time_map(self, fig, ax, year, peaks=True, max = True, colorbar=True, DOY_range=None):
        """Map the day-of-year of a phenological event across all pixels for one year.

        For each valid pixel within the 1 km-inset lake boundary, extracts the
        summer peak or green-up midpoint DOY (restricted to DOY 160–250) for the
        given year and displays it as a spatial heatmap with a fixed colorbar.

        Parameters
        ----------
        fig : matplotlib.figure.Figure
            Figure used to attach the colorbar.
        ax : matplotlib.axes.Axes
            Axes on which to draw the heatmap.
        year : int
            Calendar year for which to map the phenological event.
        peaks : bool, optional
            If True (default), map summer peaks (pks_x / pks_y).
            If False, map green-up midpoints (green_up_mid_x / green_up_mid_y).
        max : bool, optional
            If True (default), use the highest-amplitude peak in the DOY window.
            If False, use the first peak in the window.

        Returns
        -------
        matplotlib.image.AxesImage
            The imshow image object.
        """

        if DOY_range is None:
            DOY_range = [160,250]

        lake_id = int(self.lakeID)

        # lake_row = self.gdf[self.gdf["id"] == lake_id]
        geom = self.geometry #lake_row.geometry.iloc[0]
        # buffered_geom = self.shrink_geometry(geom)
        buffered_geom_prepared = self.prepped_geom

        var_x = "pks_x" if peaks else "green_up_mid_x"
        var_y = "pks_y" if peaks else "green_up_mid_y"
        extrema_label = "Peak" if peaks else "Green Mid Up"

        map_data, extent = f.grab_time_data(self.e_path, self.p_path, self.valid_coords,
                                            buffered_geom_prepared, var_x, var_y, year, max, 
                                            DOY_start = DOY_range[0], DOY_end = DOY_range[1])
        im = f.plot_map_data("rainbow", map_data, extent, ax, cmap_extent=DOY_range)
        f.plot_lake_outline(geometry=geom, ax=ax)
        if colorbar:
            f.set_labels(ax, fig, im,
                title=f"{extrema_label} Day of Year\n Lake ID: {lake_id}\n Year: {year}",
                colorbar_label="Day of Year",
                colorbar_ticks=np.arange(*DOY_range,30))
        else:
            ax.set_title(f"{extrema_label} Day of Year\n Lake ID: {lake_id}\n Year: {year}")
            ax.set_xlabel("Lon index")
            ax.set_ylabel("Lat index")
        return im
    

    def time_map_panel(self, years, nrow, ncol, peaks = True, max = True,DOY_range=None):
        extrema_label = "Peak" if peaks else "Green Mid Up"
        fig, axs = plt.subplots(nrow, ncol, constrained_layout=True, squeeze=False, figsize=(ncol * 5, nrow * 4))
        im = None
        outline_handle, outline_label = None, None
        for year, ax in zip(years, axs.flatten()):
            im = self.time_map(fig=fig, ax=ax, year=year, peaks=peaks, max=max, colorbar=False,DOY_range=DOY_range)
            if outline_handle is None:
                handles, labels = ax.get_legend_handles_labels()
                if handles:
                    outline_handle, outline_label = handles[0], labels[0]
            ax.set_title(str(year), fontsize=20)
            ax.set_ylabel("Lat index", fontsize=15)
            ax.set_xlabel("Lon index", fontsize=15)
            ax.tick_params(labelsize=15)

        for ax in axs.flatten()[len(years):]:
            ax.set_visible(False)

        if im is not None:
            cbar = fig.colorbar(im, ax=axs.ravel().tolist(), location="right", shrink=0.8)
            cbar.set_label("Day of Year", fontsize=20)
            cbar.set_ticks([160, 185, 215, 250])
            cbar.ax.tick_params(labelsize=15)

        # single shared legend for the lake outline, placed next to the colorbar
        # instead of repeating it inside every subplot
        if outline_handle is not None:
            fig.legend([outline_handle], [outline_label], loc="upper left",
                       bbox_to_anchor=(1.0, 1.0), fontsize=15, frameon=False)

        fig.suptitle(f"{extrema_label} Day of Year", fontsize=25)

        return fig , axs


    def single_day_map(self, date):
        """Plot a spatial map of chlorophyll values for a single observation date.

        Creates a new figure internally. Values failing the fill-value or QA filter
        are masked. The lake outline is overlaid from the shapefile geometry.

        Parameters
        ----------
        date : datetime.datetime
            A UTC-aware timestamp matching one entry in the extract time axis.
            Typically obtained from a pixel time series, e.g.:
                series = vis.load_pixel_data(i, j)
                vis.single_day_map(series.idxmax())

        Returns
        -------
        matplotlib.image.AxesImage
            The imshow image object.

        Raises
        ------
        ValueError
            If date is not found in the extract time axis.
        """
        date = pd.Timestamp(date)
        if date.tzinfo is None:
            date = date.tz_localize("UTC")
        else:
            date = date.tz_convert("UTC")
        g = self._load_extracted_globals()
        t_all = g["t_all"]
        idx = np.argwhere(f.datenum_to_datetime(t_all) == date)
        if len(idx) == 0:
            raise ValueError(f"Date {date} not found in the extract time axis.")
        time_index = idx[0, 0]

        with netCDF4.Dataset(self.e_path) as nc:
            values = np.array(nc.variables[g["variable"]][time_index, :, :])
            qa     = np.array(nc.variables[g["qa"]][time_index, :, :])
            lats   = nc.variables["lat"][:]
            lons   = nc.variables["lon"][:]

        mask = (values != -9999) & (qa == 0)
        lake_id = int(self.lakeID)
        lake_row = self.gdf[self.gdf["id"] == lake_id]
        geom = lake_row.geometry.iloc[0]
        buffered_geom = self.shrink_geometry(geom)
        buffered_geom_prepared = prep(buffered_geom)

        map_data = np.full(values.shape, np.nan)
        for i in range(values.shape[0]):
            for j in range(values.shape[1]):
                if mask[i, j] and buffered_geom_prepared.contains(Point(lons[j], lats[i])):
                    map_data[i, j] = values[i, j]
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        im = ax.imshow(map_data, cmap="winter", aspect="auto", origin="lower",
                        extent=[lons.min(), lons.max(), lats.min(), lats.max()])

        label = False
        if geom.geom_type == "Polygon":
            x, y = geom.exterior.xy
            ax.plot(x, y, color="black", linewidth=1, label="Lake Outline")
        elif geom.geom_type == "MultiPolygon":
            for poly in geom.geoms:
                x, y = poly.exterior.xy
                ax.plot(x, y, color="black", linewidth=1,
                        label="Lake Outline" if not label else None)
                label = True

        fig.colorbar(im, orientation='vertical', label=f'{self.variable} (ug/L)')
        plt.xticks([])
        plt.yticks([])
        plt.title(date.strftime("%Y-%m-%d"))
        return im
