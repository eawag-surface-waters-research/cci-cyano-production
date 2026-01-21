import os
import netCDF4
import logging
import numpy as np
from csaps import csaps
from scipy.signal import find_peaks
from datetime import datetime, timedelta
from matplotlib import pyplot as plt


def phenology(lake, p):
    input_file = os.path.join(p["out_folder"], "extract", p["variable"], "{}.nc".format(lake["id"]))
    output_file = os.path.join(p["out_folder"], "phenology", p["variable"], "{}.nc".format(lake["id"]))

    if os.path.isfile(output_file):
        logging.info(f"File {lake['id']}.nc already exists")
        return

    if os.path.isfile(input_file):
        logging.info(f"Starting phenology computation for {lake['id']}")
    else:
        logging.info(f"Input file not available for {lake['id']}")
        return

    with netCDF4.Dataset(input_file) as nc:
        valid_coords = np.argwhere(nc.variables["summary"][:, :] > 2)
        valid_coords = [[27, 50]]
        for x, y in valid_coords:
            values = np.array(nc.variables[p["variable"]][:, x, y])
            time = np.array([datetime.utcfromtimestamp(int(ts)) for ts in np.array(nc.variables["time"])])
            mask = values != -9999
            if p["qa_filter"]:
                qa = nc.variables[p["variable"]][:, x, y]
                mask = qa == 0 & mask
            values = values[mask]
            time = time[mask]

            result = smooth_cubic_spline(time, values, p["spline_min_phase_length"], p["spline_min_relative_amplitude"],
                                         p["spline_min_phase_data"], smoothing_change=1e-12, max_iterations=1e5)
            fig, ax = plot_cubic_spline(time, values, result)
            plt.show()


def smooth_cubic_spline(time, values, spline_min_phase_length, spline_min_relative_amplitude, spline_min_phase_data,
                        smoothing_min=0.0, smoothing_max=1.0, smoothing_change=1e-12, max_iterations=1e5):

    smoothing = smoothing_min + 0.5 * (smoothing_max - smoothing_min)
    smoothing_old = smoothing_min
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
            smoothing
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


def cubic_spline(time, values, cond_min_phase_length, cond_min_rel_amp, cond_min_phase_nr_data, smoothing):
    time = np.asarray(time)
    values = np.asarray(values)

    # Convert datetime to numeric (days since first date)
    time_numeric = np.array([(t - time.min()).days for t in time])

    # Create regular grid for evaluation
    smooth_x_axis_numeric = np.arange(time_numeric.min(), time_numeric.max() + 1, 1)

    # Convert numeric grid back to datetime for output
    smooth_x_axis = np.array([time.min() + timedelta(days=int(d)) for d in smooth_x_axis_numeric])

    # Fit smoothing spline using numeric values
    smooth_y_data = csaps(time_numeric, values, smooth_x_axis_numeric, smooth=smoothing)

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

    # Combine and sort peaks and troughs (use numeric for calculations)
    x_pks_numeric = smooth_x_axis_numeric[pks_indices]
    x_trgs_numeric = smooth_x_axis_numeric[trgs_indices]
    x_pks_and_trgs_numeric = np.concatenate([x_pks_numeric, x_trgs_numeric])
    y_pks_and_trgs = np.concatenate([y_pks, y_trgs])

    if len(x_pks_and_trgs_numeric) >= 2:
        sort_indices = np.argsort(x_pks_and_trgs_numeric)
        x_pks_and_trgs_sorted = x_pks_and_trgs_numeric[sort_indices]
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
                (time_numeric >= x_pks_and_trgs_sorted[i]) &
                (time_numeric <= x_pks_and_trgs_sorted[i + 1])
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


def plot_cubic_spline(time, values, result, title=None):
    """
    Plot the original data and cubic spline fit with peaks and troughs.

    Parameters:
    -----------
    time : array-like
        Original time data
    values : array-like
        Original values
    result : dict
        Output from cubic_spline function
    title : str, optional
        Plot title
    """
    fig, ax = plt.subplots(figsize=(12, 6))

    # Original data
    ax.scatter(time, values, color='blue', alpha=0.5, label='Original data', zorder=3)

    # Smoothed spline
    ax.plot(result['smooth_x_axis'], result['smooth_y_data'], color='red', linewidth=2,
            label=f"Spline (smoothing={result['smoothing']:.2e})")

    # Peaks
    if len(result['x_pks']) > 0:
        ax.scatter(result['x_pks'], result['y_pks'], color='green', s=100, marker='^', label='Peaks', zorder=4)

    # Troughs
    if len(result['x_trgs']) > 0:
        ax.scatter(result['x_trgs'], result['y_trgs'], color='orange', s=100, marker='v', label='Troughs', zorder=4)

    ax.set_xlabel('Time')
    ax.set_ylabel('Value')
    ax.legend()
    ax.grid(True, alpha=0.3)

    if title:
        ax.set_title(title)
    else:
        ax.set_title(f"Cubic Spline Fit (conditions satisfied: {result['conditions_satisfied']})")

    # Rotate x-axis labels if datetime
    fig.autofmt_xdate()

    plt.tight_layout()
    return fig, ax