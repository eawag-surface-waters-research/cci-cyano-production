import pandas as pd
import netCDF4
import matplotlib
matplotlib.use("Agg")
import numpy as np
from pathlib import Path
from matplotlib.patches import Rectangle, Patch
import functions as f
from visualization.base import PhenologyVisualizationBase


class HeatmapVisualization(PhenologyVisualizationBase):

    graph_type = "heatmap"
    graph_methods = frozenset({
        "create_heatmap_output",
        "yearly_heatmap_pixel",
        "yearly_heatmap_lake",
    })
    
    def create_heatmap_output(self, latitude_idx, longitude_idx, start_year=2002, end_year=2024, fraction=False, qa=None):
        """Return peak/trough counts or lake-wide fractions per year and quarter.

        When fraction=False, counts are read from the single pixel at (latitude_idx,
        longitude_idx). When fraction=True, all pixels in the lake are aggregated and
        each quarter value is expressed as the fraction of that year's total events
        occurring in that quarter (0.0 – 1.0); years with no events return 0.0.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel. Only used when fraction=False.
        longitude_idx : int
            Column (lon) index of the pixel. Only used when fraction=False.
        start_year : int, optional
            First calendar year to include (inclusive). Default 2002.
        end_year : int, optional
            Last calendar year to include (inclusive). Default 2024.
        fraction : bool, optional
            If False (default), return per-pixel integer counts.
            If True, return lake-wide fractions aggregated across all pixels.
        qa : set of int or None, optional
            QA levels to include when counting peaks and troughs (fraction=False only).
            E.g. {0} for good only, {0, 1} for good and fair. None includes all
            levels (default). Valid values: 0 (Good), 1 (Fair), 2 (Poor).

        Returns
        -------
        dict
            Keys are years (int). Values are lists of 4 tuples, one per quarter
            (Jan-Mar, Apr-Jun, Jul-Sep, Oct-Dec). Each tuple is
            (n_peaks, n_troughs) when fraction=False, or
            (peaks_fraction, troughs_fraction) when fraction=True.
        """
        quarters = [(1, 3), (4, 6), (7, 9), (10, 12)]
        result = {}
        if not fraction:
            with netCDF4.Dataset(self.p_path) as nc:
                pks_x_raw = np.array(nc.variables["pks_x"][latitude_idx, longitude_idx, :])
                pk_mask = ~np.isnan(pks_x_raw)
                pks_x = f.unix_to_datetime(pks_x_raw[pk_mask])
                trgs_x_raw = np.array(nc.variables["trgs_x"][latitude_idx, longitude_idx, :])
                trg_mask = ~np.isnan(trgs_x_raw)
                trgs_x = f.unix_to_datetime(trgs_x_raw[trg_mask])

                if qa is not None:
                    qa_set = set(qa)
                    pks_qa_arr = np.array(nc.variables["pks_qa"][latitude_idx, longitude_idx, :])[pk_mask]
                    trgs_qa_arr = np.array(nc.variables["trgs_qa"][latitude_idx, longitude_idx, :])[trg_mask]
                    pks_x = pks_x[np.isin(pks_qa_arr, list(qa_set))]
                    trgs_x = trgs_x[np.isin(trgs_qa_arr, list(qa_set))]
        
                for year in range(start_year, end_year + 1):
                    year_counts = []
                    for (q_start, q_end) in quarters:
                        n_pks = sum(
                            1 for d in pks_x
                            if d.year == year and q_start <= d.month <= q_end
                        )
                        n_trgs = sum(
                            1 for d in trgs_x
                            if d.year == year and q_start <= d.month <= q_end
                        )
                        year_counts.append((n_pks, n_trgs))
                    result[year] = year_counts
        else:
            with netCDF4.Dataset(self.p_path) as nc:
                vpx, vpqa = nc.variables["pks_x"], nc.variables["pks_qa"]
                vtx, vtqa = nc.variables["trgs_x"], nc.variables["trgs_qa"]

                # restrict to pixels inside the 1 km-inset lake boundary (self.valid_idx_prep),
                # matching metric_map/time_map/lake_bloom_kde/qa_boxplot_lake, instead of the
                # full raw grid which includes border/mixed pixels. Read per-pixel
                # (nc.variables[var][i, j, :]) rather than a bulk [:, :, :] read - netCDF4 1.7.4
                # misattributes data between pixels when read that way for files with an
                # unlimited 'record' dimension, which would corrupt this restriction
                pks_x_parts, pks_qa_parts = [], []
                trgs_x_parts, trgs_qa_parts = [], []
                for i, j in self.valid_idx_prep:
                    px_raw = np.asarray(vpx[i, j, :])
                    pmask = ~np.isnan(px_raw)
                    if pmask.any():
                        pks_x_parts.append(px_raw[pmask])
                        pks_qa_parts.append(np.asarray(vpqa[i, j, :])[pmask])

                    tx_raw = np.asarray(vtx[i, j, :])
                    tmask = ~np.isnan(tx_raw)
                    if tmask.any():
                        trgs_x_parts.append(tx_raw[tmask])
                        trgs_qa_parts.append(np.asarray(vtqa[i, j, :])[tmask])

                pks_x_raw = np.concatenate(pks_x_parts) if pks_x_parts else np.array([])
                pks_qa_arr = np.concatenate(pks_qa_parts) if pks_qa_parts else np.array([])
                trgs_x_raw = np.concatenate(trgs_x_parts) if trgs_x_parts else np.array([])
                trgs_qa_arr = np.concatenate(trgs_qa_parts) if trgs_qa_parts else np.array([])

                pks_x = f.unix_to_datetime(pks_x_raw)
                trgs_x = f.unix_to_datetime(trgs_x_raw)

                if qa is not None:
                    qa_set = set(qa)
                    pks_x = pks_x[np.isin(pks_qa_arr, list(qa_set))]
                    trgs_x = trgs_x[np.isin(trgs_qa_arr, list(qa_set))]

                for year in range(start_year, end_year + 1):
                    year_fractions = []
                    for (q_start, q_end) in quarters:
                        n_pks = sum(
                            1 for d in pks_x
                            if d.year == year and q_start <= d.month <= q_end
                        )
                        yearly_pks = sum(1 for d in pks_x if d.year == year)
                        pks_fraction = n_pks / yearly_pks if yearly_pks > 0 else 0.0

                        n_trgs = sum(
                            1 for d in trgs_x
                            if d.year == year and q_start <= d.month <= q_end
                        )
                        yearly_trgs = sum(1 for d in trgs_x if d.year == year)
                        trgs_fraction = n_trgs / yearly_trgs if yearly_trgs > 0 else 0.0

                        year_fractions.append((pks_fraction, trgs_fraction))
                    result[year] = year_fractions
        return result


    def yearly_heatmap_pixel(self, latitude_idx, longitude_idx, color_scheme='pink-blue', show_gaps=False, qa=None, mask_cells=None, mask_color='grey'):
        """Plot a bivariate heatmap of peak and trough counts or fractions by year and quarter.

        Each cell in the heatmap represents one calendar quarter of one year. The
        cell colour encodes two variables simultaneously using a 4×4 bivariate
        colour palette from color_sets_4x4.

        When whole_lake=False, counts for the single pixel at (latitude_idx, longitude_idx)
        are binned into four levels (0, 1, 2, 3+) and the cell is coloured from the
        discrete 4×4 grid using f.bivariate_legend.

        When whole_lake=True, peak and trough events are aggregated across all
        lake pixels and each cell shows the fraction of that year's total events
        falling in that quarter. Colours are interpolated continuously across the
        4×4 grid using bivariate_continuous_legend(), and the legend is a smooth 2-D gradient
        rendered by bivariate_continuous_legend.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel. Only used when whole_lake=False.
        longitude_idx : int
            Column (lon) index of the pixel. Only used when whole_lake=False.
        color_scheme : str, optional
            Key into color_sets_4x4 selecting the bivariate palette.
            One of 'pink-blue', 'teal-red', 'teal-red1', 'blue-orange'.
            Default 'pink-blue'.
        whole_lake : bool, optional
            If False (default), plot per-pixel counts with a discrete legend.
            If True, plot lake-wide fractions with a continuous gradient legend.
        show_gaps : bool, optional
            If True, overlay a grey strip on the left of each cell whose width
            is proportional to the fraction of that quarter covered by data gaps.
            Default False.
        mask_cells : list of str or None, optional
            Cells to manually override with a solid colour after all other
            drawing, regardless of the underlying data. Each entry is a string
            '<year> <quarter>', e.g. ``['2015 Apr-Jun', '2016 Jul-Sep']``.
            Valid quarter labels: 'Jan-Mar', 'Apr-Jun', 'Jul-Sep', 'Oct-Dec'.
            Default None (no masking).
        mask_color : str, optional
            Colour used to fill masked cells. Default 'grey'.

        Returns
        -------
        fig : matplotlib.figure.Figure
            The figure containing the heatmap.
        ax : matplotlib.axes.Axes
            The axes on which the heatmap is drawn.
        """
        heatmap_data = self.create_heatmap_output(latitude_idx=latitude_idx, longitude_idx=longitude_idx, fraction=False, qa=qa)
        fig, ax, _ = f.create_empty_heatmap()
        g  = self._load_extracted_globals()
        lat, lon = g["lat"], g["lon"]
        lat_val = float(lat[latitude_idx])
        lon_val = float(lon[longitude_idx])
        textstr = f"Yearly Heatmap for Pixel\n lat, lon: {lat_val:.4f}, {lon_val:.4f}\n {self.variable}, Lake ID:{self.lakeID}"
        ax.set_title(textstr)

        color_set = PhenologyVisualizationBase.color_sets_4x4[color_scheme]

        for year, quarters in heatmap_data.items():
            for q_idx, (n_pks, n_trgs) in enumerate(quarters):
                pk_bin  = min(n_pks,  3)
                trg_bin = min(n_trgs, 3)
                color = color_set[trg_bin * 4 + pk_bin]
                ax.add_patch(Rectangle((q_idx, year - 1), 1, 1, facecolor=color, edgecolor='none'))

        if show_gaps:
            px = self._load_pixel_data(latitude_idx, longitude_idx)
            gap_rects = f.gap_rectangles(px["gap_starts"], px["gap_ends"])
            gap_patch_added = False
            for year, q_idx, x_offset, x_width in gap_rects:
                if 2002 <= year <= 2024:
                    label = "Data gap" if not gap_patch_added else None
                    ax.add_patch(Rectangle((q_idx + x_offset, year - 1), x_width, 1,
                                           facecolor='grey', alpha=0.6, edgecolor='none', label=label))
                    gap_patch_added = True
            if gap_patch_added:
                ax.legend(handles=[Patch(facecolor='grey', alpha=0.6, label='Data gap')],
                          bbox_to_anchor=(1.5, 0.5),
                          bbox_transform=ax.transAxes, fontsize=8)

        if mask_cells:
            for year, q_idx in f.parse_mask_cells(mask_cells):
                ax.add_patch(Rectangle((q_idx, year - 1), 1, 1,
                                       facecolor=mask_color, edgecolor='none', zorder=4))

        ax_legend = ax.inset_axes([1.2, 0.7, 0.3, 0.3], transform = ax.transAxes)
        f.bivariate_legend(ax_legend, color_set)

        return fig, ax


    def yearly_heatmap_lake(self, color_scheme="pink-blue", qa=None, mask_cells=None, mask_color='grey'):
        """Plot a bivariate heatmap of lake-wide peak/trough fractions by year and quarter.

        Parameters
        ----------
        color_scheme : str, optional
            Key into color_sets_4x4 selecting the bivariate palette. Default 'pink-blue'.
        qa : set or list of int or None, optional
            QA levels to include when counting events. E.g. {0} for good only,
            {0, 1} for good and fair. None includes all levels (default).
        mask_cells : list of str or None, optional
            Cells to manually override with a solid colour after all other
            drawing, regardless of the underlying data. Each entry is a string
            '<year> <quarter>', e.g. ``['2015 Apr-Jun', '2016 Jul-Sep']``.
            Valid quarter labels: 'Jan-Mar', 'Apr-Jun', 'Jul-Sep', 'Oct-Dec'.
            Default None (no masking).
        mask_color : str, optional
            Colour used to fill masked cells. Default 'grey'.

        Returns
        -------
        fig : matplotlib.figure.Figure
        ax : matplotlib.axes.Axes
        """
        # lat and lon are not needed as the heatmap uses all pixels from the lake, thus they can be arbitrary
        heatmap_data = self.create_heatmap_output(latitude_idx=-1, longitude_idx=-1, fraction=True, qa=qa)
        fig, ax, _ = f.create_empty_heatmap()
        textstr = f"Yearly Heatmap for Lake ID: {self.lakeID}\n  {self.variable}"
        ax.set_title(textstr)

        color_set = PhenologyVisualizationBase.color_sets_4x4[color_scheme]

        for year, quarters in heatmap_data.items():
            for q_idx, (pks_frac, trgs_frac) in enumerate(quarters):
                color = f.interpolate_from_color_set(pks_frac, trgs_frac, color_set)
                ax.add_patch(Rectangle((q_idx, year - 1), 1, 1, facecolor=color, edgecolor='none'))

        if mask_cells:
            for year, q_idx in f.parse_mask_cells(mask_cells):
                ax.add_patch(Rectangle((q_idx, year - 1), 1, 1,
                                       facecolor=mask_color, edgecolor='none', zorder=4))

        ax_legend = ax.inset_axes([1.2, 0.7, 0.3, 0.3], transform=ax.transAxes)
        f.bivariate_continuous_legend(ax_legend, color_set)

        return fig, ax
