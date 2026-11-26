import pandas as pd
import netCDF4
import matplotlib
matplotlib.use("Agg")
import numpy as np
from pathlib import Path
import warnings
import functions as f
from visualization.base import PhenologyVisualizationBase

class BoxplotVisualization(PhenologyVisualizationBase):

    graph_type = "boxplot"
    graph_methods = frozenset({
        "pair_phenology_events",
        "qa_boxplot",
        "qa_boxplot_lake",
        "_boxplot_by_qa",
    })
    
    def pair_phenology_events(self, other, latitude_idx, longitude_idx, metric="pks", tolerance_days=4):
        """Pair phenology metric events between self and another phenology product/version.

        Each event may be used in at most one pair. Candidate pairs (any self/other
        event within `tolerance_days`) are ranked by how close they are in time and
        claimed greedily closest-first, so once an event is paired it is removed from
        contention - no event is ever paired twice. Everything is done with an all-pairs
        distance matrix and pandas/numpy set ops (drop_duplicates, boolean masks), so no
        Python loop is used over the (typically ragged, variable-length) per-pixel event
        arrays.

        Parameters
        ----------
        other : PhenologyVisualization
            Second instance to pair against (e.g. phycocyanin vs chla, or v2.1 vs v3.0),
            for the same lake/pixel.
        latitude_idx : int
            Row (lat) index of the pixel.
        longitude_idx : int
            Column (lon) index of the pixel.
        metric : str, optional
            One of 'pks', 'trgs', 'midUP', 'midDOWN', 'onsetUP', 'onsetDOWN', 'advUP',
            'advDOWN'. Defaults to 'pks' (summer peaks).
        tolerance_days : int, optional
            Maximum gap in days between two events for them to count as a pair.
            Defaults to 4.

        Returns
        -------
        pandas.DataFrame
            One row per paired event, with columns time_<label>/value_<label> for
            each side plus 'day_diff' (absolute gap in days). Unpaired events are
            not included, and no event appears in more than one row.

        Raises
        ------
        ValueError
            If `metric` is not a recognised phenology metric key.
        """
        x_key, y_key, qa_key = f"{metric}_x", f"{metric}_y", f"{metric}_qa"

        px_self = self._load_pixel_data(latitude_idx, longitude_idx)
        px_other = other._load_pixel_data(latitude_idx, longitude_idx)

        if x_key not in px_self or x_key not in px_other:
            raise ValueError(f"'{metric}' is not a valid phenology metric.")

        self_label, other_label = self.variable, other.variable
        # if self_label == other_label:
        #     self_label, other_label = f"{self_label}_v{self.version}", f"{other_label}_v{other.version}"

        has_qa = qa_key in px_self and qa_key in px_other

        columns = [f"time_{self_label}", f"value_{self_label}",
                   f"time_{other_label}", f"value_{other_label}", "day_diff"]
        if has_qa:
            columns += [f"qa_{self_label}", f"qa_{other_label}"]

        time_self = pd.to_datetime(px_self[x_key]).values
        time_other = pd.to_datetime(px_other[x_key]).values

        self_idx_matched, other_idx_matched, day_diff_matched = f.greedy_match_within_tolerance(
            time_self, time_other, tolerance_days
        )
        if len(self_idx_matched) == 0:
            return pd.DataFrame(columns=columns)

        paired_dict = {
            f"time_{self_label}": time_self[self_idx_matched],
            f"value_{self_label}": px_self[y_key][self_idx_matched],
            f"time_{other_label}": time_other[other_idx_matched],
            f"value_{other_label}": px_other[y_key][other_idx_matched],
            "day_diff": day_diff_matched,
        }
        if has_qa:
            paired_dict[f"qa_{self_label}"] = px_self[qa_key][self_idx_matched]
            paired_dict[f"qa_{other_label}"] = px_other[qa_key][other_idx_matched]

        paired = pd.DataFrame(paired_dict)

        return paired.sort_values(f"time_{self_label}").reset_index(drop=True)


    def qa_boxplot(self, latitude_idx, longitude_idx, ax, metric="pks", start=0, end=9999,
                other=None, tolerance_days=4, qa_source="self"):
        """Boxplot of a phenology metric's values grouped by QA level, for one pixel.

        X-axis groups are QA levels (Good/Fair/Poor). If `other` is None, the
        y-axis is the metric value itself (e.g. summer peak chla or phycocyanin).
        If `other` is given, self's and other's events are paired one-to-one
        (within `tolerance_days`, see pair_phenology_events) and the y-axis becomes
        the ratio value_self / value_other instead - e.g. pass other=chla_v3 on a
        phycocyanin instance to get a Phycocyanin/chla ratio boxplot. Only 'pks'
        and 'trgs' carry a QA flag per event, so those are the only supported
        metrics either way.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel.
        longitude_idx : int
            Column (lon) index of the pixel.
        ax : matplotlib.axes.Axes
            Axes on which to draw the boxplot.
        metric : str, optional
            'pks' (summer peaks, default) or 'trgs' (winter troughs).
        start : int, optional
            First year to include (inclusive; keyed on self's event time when
            `other` is given). 0 = earliest in the series.
        end : int, optional
            Last year to include (inclusive). 9999 = latest in the series.
        other : PhenologyVisualization, optional
            If given, plot the value_self / value_other ratio of paired events
            instead of self's raw metric values (e.g. self=phycocyanin,
            other=chla gives Phyco/chla).
        tolerance_days : int, optional
            Only used when `other` is given: maximum gap in days for two events
            to count as a pair. Default 4.
        qa_source : {'self', 'other', 'matched'}, optional
            Only used when `other` is given: which side's QA flag to group by.
            'matched' keeps only pairs where both sides have the same QA level
            (dropping e.g. a Good self event paired with a Poor other event)
            and groups by that shared QA value. Default 'self'.

        Returns
        -------
        dict or None
            The dict returned by ax.boxplot (boxes/medians/whiskers/...), or None
            if there is no QA-labelled data to plot for this pixel/year range.

        Raises
        ------
        ValueError
            If `metric` is not 'pks' or 'trgs', or `qa_source` is not 'self'/'other'.
        """
        if metric not in ("pks", "trgs"):
            raise ValueError("qa_boxplot only supports metrics with a QA flag: 'pks' or 'trgs'.")

        metric_label = "Peak" if metric == "pks" else "Trough"

        if other is None:
            px = self._load_pixel_data(latitude_idx, longitude_idx)
            plotting_data = f.grab_plotting_variables(start=start, end=end, pixel_data=px, variables=[metric])

            if metric not in plotting_data:
                warnings.warn(f"No {metric} data to plot for lake ID {self.lakeID}.")
                return None

            _, values, qa_values = plotting_data[metric]
            var_label = self.get_plot_config("var", self.variable)["label"]
            title = f"{var_label} {metric_label} Values by QA\n Lake ID: {self.lakeID}"
            ylabel = f"{var_label} [ug/L]"
        else:
            if qa_source not in ("self", "other", "matched"):
                raise ValueError("qa_source must be 'self', 'other', or 'matched'.")

            self_label, other_label = self.variable, other.variable
            self_pretty = self.get_plot_config("var", self.variable)["label"]
            other_pretty = other.get_plot_config("var", other.variable)["label"]

            paired = self.pair_phenology_events(other, latitude_idx, longitude_idx,
                                                metric=metric, tolerance_days=tolerance_days)
            if len(paired) == 0:
                warnings.warn(f"No paired {metric} events to plot for lake ID {self.lakeID}.")
                return None

            time_col = f"time_{self_label}"
            if start != 0:
                paired = paired[paired[time_col].dt.year >= start]
            if end != 9999:
                paired = paired[paired[time_col].dt.year <= end]
            if len(paired) == 0:
                warnings.warn(f"No paired {metric} events in {start}-{end} for lake ID {self.lakeID}.")
                return None

            if qa_source == "matched":
                # only keep pairs where both sides agree on QA - so the ratio is
                # never computed from a Good event divided by a Poor one, say
                paired = paired[paired[f"qa_{self_label}"] == paired[f"qa_{other_label}"]]
                if len(paired) == 0:
                    warnings.warn(f"No {metric} pairs with matching QA on both sides for lake ID {self.lakeID}.")
                    return None
                qa_values = paired[f"qa_{self_label}"].to_numpy()
            else:
                qa_label = self_label if qa_source == "self" else other_label
                qa_values = paired[f"qa_{qa_label}"].to_numpy()

            values = (paired[f"value_{self_label}"] / paired[f"value_{other_label}"]).to_numpy()
            title = f"{metric_label} {self_pretty}/{other_pretty} Ratio by QA\n Lake ID: {self.lakeID}"
            ylabel = f"{self_label}/{other_label} ratio"

        return self._boxplot_by_qa(ax, values, qa_values, title=title, ylabel=ylabel)


    def qa_boxplot_lake(self, ax, metric="pks", start=0, end=9999, other=None, tolerance_days=4, qa_source="self"):
        """Boxplot of a phenology metric's values grouped by QA level, pooled across the whole lake.

        Same idea as qa_boxplot, but instead of a single pixel it pools events from
        every pixel in the lake. If `other` is None, this reads the raw pks/trgs
        arrays directly (like create_heatmap_output/yearly_heatmap_lake) rather than
        looping pixel by pixel, since no per-pixel correspondence is needed. If
        `other` is given, a peak must be paired against the *same pixel's* peak in
        `other` (not any pixel), so pair_phenology_events is run per valid pixel
        (self.valid_idx_prep) and the results are pooled - this is slower for large
        lakes.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes on which to draw the boxplot.
        metric : str, optional
            'pks' (summer peaks, default) or 'trgs' (winter troughs).
        start : int, optional
            First year to include (inclusive). 0 = earliest in the series.
        end : int, optional
            Last year to include (inclusive). 9999 = latest in the series.
        other : PhenologyVisualization, optional
            If given, plot the value_self / value_other ratio of paired events
            instead of self's raw metric values (e.g. self=phycocyanin,
            other=chla gives Phyco/chla), pooled across all lake pixels.
        tolerance_days : int, optional
            Only used when `other` is given: maximum gap in days for two events
            at the same pixel to count as a pair. Default 4.
        qa_source : {'self', 'other', 'matched'}, optional
            Only used when `other` is given: which side's QA flag to group by.
            'matched' keeps only pairs where both sides have the same QA level.
            Default 'self'.

        Returns
        -------
        dict or None
            The dict returned by ax.boxplot, or None if there is no QA-labelled
            data to plot.

        Raises
        ------
        ValueError
            If `metric` is not 'pks' or 'trgs', or `qa_source` is not recognised.
        """
        if metric not in ("pks", "trgs"):
            raise ValueError("qa_boxplot_lake only supports metrics with a QA flag: 'pks' or 'trgs'.")

        metric_label = "Peak" if metric == "pks" else "Trough"

        x_key, y_key, qa_key = f"{metric}_x", f"{metric}_y", f"{metric}_qa"

        if other is None:
            values_parts, qa_parts, times_parts = [], [], []
            with netCDF4.Dataset(self.p_path) as nc:
                vx, vy, vqa = nc.variables[x_key], nc.variables[y_key], nc.variables[qa_key]
                for i, j in self.valid_idx_prep:
                    # read one pixel at a time (nc.variables[var][i, j, :]) rather than
                    # a bulk nc.variables[var][:, :, :] read - netCDF4 1.7.4 silently
                    # misaligns data read this way for files with an unlimited 'record'
                    # dimension, verified by comparing against the trusted per-pixel reads
                    # used everywhere else in this class (e.g. _load_pixel_data)
                    x_raw = np.asarray(vx[i, j, :])
                    mask = ~np.isnan(x_raw)
                    if not mask.any():
                        continue
                    values_parts.append(np.asarray(vy[i, j, :])[mask])
                    qa_parts.append(np.asarray(vqa[i, j, :])[mask])
                    times_parts.append(x_raw[mask])

            if not values_parts:
                warnings.warn(f"No {metric} data to plot for lake ID {self.lakeID}.")
                return None

            values, qa_values = np.concatenate(values_parts), np.concatenate(qa_parts)
            years = pd.to_datetime(np.concatenate(times_parts), unit="s", utc=True).year.to_numpy()

            year_mask = np.ones(len(values), dtype=bool)
            if start != 0:
                year_mask &= years >= start
            if end != 9999:
                year_mask &= years <= end
            values, qa_values = values[year_mask], qa_values[year_mask]

            if len(values) == 0:
                warnings.warn(f"No {metric} data to plot for lake ID {self.lakeID}.")
                return None

            var_label = self.get_plot_config("var", self.variable)["label"]
            title = f"{var_label} {metric_label} Values by QA\n Lake ID: {self.lakeID} (lake-wide)"
            ylabel = f"{var_label} [ug/L]"
        else:
            if qa_source not in ("self", "other", "matched"):
                raise ValueError("qa_source must be 'self', 'other', or 'matched'.")

            self_label, other_label = self.variable, other.variable
            self_pretty = self.get_plot_config("var", self.variable)["label"]
            other_pretty = other.get_plot_config("var", other.variable)["label"]

            self_vals, other_vals, self_qa_vals, other_qa_vals, self_times = [], [], [], [], []
            with netCDF4.Dataset(self.p_path) as nc_self, netCDF4.Dataset(other.p_path) as nc_other:
                svx, svy, svqa = nc_self.variables[x_key], nc_self.variables[y_key], nc_self.variables[qa_key]
                ovx, ovy, ovqa = nc_other.variables[x_key], nc_other.variables[y_key], nc_other.variables[qa_key]

                for i, j in self.valid_idx_prep:
                    # per-pixel reads only (see note above) - each file is opened
                    # once, then indexed pixel by pixel, avoiding both the bulk-read
                    # correctness bug and _load_pixel_data's per-pixel full-file cost
                    s_x_raw = np.asarray(svx[i, j, :])
                    s_mask = ~np.isnan(s_x_raw)
                    o_x_raw = np.asarray(ovx[i, j, :])
                    o_mask = ~np.isnan(o_x_raw)
                    if not s_mask.any() or not o_mask.any():
                        continue

                    s_time = pd.to_datetime(s_x_raw[s_mask], unit="s", utc=True).values
                    o_time = pd.to_datetime(o_x_raw[o_mask], unit="s", utc=True).values

                    s_idx, o_idx, _ = f.greedy_match_within_tolerance(s_time, o_time, tolerance_days)
                    if len(s_idx) == 0:
                        continue

                    self_vals.append(np.asarray(svy[i, j, :])[s_mask][s_idx])
                    other_vals.append(np.asarray(ovy[i, j, :])[o_mask][o_idx])
                    self_qa_vals.append(np.asarray(svqa[i, j, :])[s_mask][s_idx])
                    other_qa_vals.append(np.asarray(ovqa[i, j, :])[o_mask][o_idx])
                    self_times.append(s_time[s_idx])

            if not self_vals:
                warnings.warn(f"No paired {metric} events to plot for lake ID {self.lakeID}.")
                return None

            self_vals, other_vals = np.concatenate(self_vals), np.concatenate(other_vals)
            self_qa_vals, other_qa_vals = np.concatenate(self_qa_vals), np.concatenate(other_qa_vals)
            years = pd.DatetimeIndex(np.concatenate(self_times)).year.to_numpy()

            year_mask = np.ones(len(self_vals), dtype=bool)
            if start != 0:
                year_mask &= years >= start
            if end != 9999:
                year_mask &= years <= end
            self_vals, other_vals = self_vals[year_mask], other_vals[year_mask]
            self_qa_vals, other_qa_vals = self_qa_vals[year_mask], other_qa_vals[year_mask]

            if len(self_vals) == 0:
                warnings.warn(f"No paired {metric} events in {start}-{end} for lake ID {self.lakeID}.")
                return None

            if qa_source == "matched":
                match_mask = self_qa_vals == other_qa_vals
                self_vals, other_vals = self_vals[match_mask], other_vals[match_mask]
                self_qa_vals = self_qa_vals[match_mask]
                if len(self_vals) == 0:
                    warnings.warn(f"No {metric} pairs with matching QA on both sides for lake ID {self.lakeID}.")
                    return None
                qa_values = self_qa_vals
            else:
                qa_values = self_qa_vals if qa_source == "self" else other_qa_vals

            values = self_vals / other_vals
            title = f"{metric_label} {self_pretty}/{other_pretty} Ratio by QA\n Lake ID: {self.lakeID} (lake-wide)"
            ylabel = f"{self_label}/{other_label} ratio"

        return self._boxplot_by_qa(ax, values, qa_values, title=title, ylabel=ylabel)


    def _boxplot_by_qa(self, ax, values, qa_values, title, ylabel):
        """Draw a QA-grouped boxplot of `values` on `ax`. Shared by qa_boxplot and qa_boxplot_ratio.

        Also annotates every pair of QA groups with a two-sided Mann-Whitney U
        test (non-parametric, no normality assumption), drawn as a significance
        bracket above the boxes: '***' p<0.001, '**' p<0.01, '*' p<0.05, 'ns' otherwise.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes on which to draw the boxplot.
        values : array-like
            Values to group and plot (one entry per event).
        qa_values : array-like
            QA level (0/1/2) for each entry in `values`, same length.
        title : str
            Plot title.
        ylabel : str
            Y-axis label.

        Returns
        -------
        dict or None
            The dict returned by ax.boxplot, or None if no QA level has data.
        """
        values = np.asarray(values)
        qa_values = np.asarray(qa_values)

        present_levels = [q for q in self.QA_LEVELS if (qa_values == q).any()]
        if not present_levels:
            warnings.warn(f"No QA-labelled data to plot for lake ID {self.lakeID}.")
            return None

        groups = [values[qa_values == q] for q in present_levels]
        labels = [self.get_plot_config("qa", q)["label"] for q in present_levels]
        colors = [self.get_plot_config("qa", q)["style"]["color"] for q in present_levels]

        box = ax.boxplot(groups, tick_labels=labels, patch_artist=True)
        for patch, color in zip(box["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.5)
        for median in box["medians"]:
            median.set_color("black")

        ax.set_title(title)
        ax.set_xlabel("QA")
        ax.set_ylabel(ylabel)

        f.annotate_pairwise_significance(ax, groups)

        return box

