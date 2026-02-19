"""Compare Python NetCDF phenology output against MATLAB .mat output."""

import numpy as np
from scipy.io import loadmat
from datetime import datetime
import netCDF4


def unix_to_datenum(unix_ts):
    """Convert Unix timestamps to MATLAB datenum values."""
    return np.array([datetime.utcfromtimestamp(ts).toordinal() + 366
                     + (ts - datetime.utcfromtimestamp(ts).replace(
                         hour=0, minute=0, second=0, microsecond=0).timestamp()) / 86400
                     for ts in unix_ts])


def compare(python_file, matlab_file, input_file, matlab_data_file=None, matlab_info_file=None,
            pixel=None, rtol=1e-6, plot=False, variable=None):
    """Compare Python and MATLAB phenology outputs.

    Parameters
    ----------
    python_file : str
        Path to Python NetCDF output (new format: variables shaped (lat, lon, record)).
    matlab_file : str
        Path to MATLAB .mat output.
    input_file : str
        Path to Python input NetCDF file (for pixel coordinate mapping).
    matlab_data_file : str, optional
        Path to MATLAB input data file (for plotting raw data).
    matlab_info_file : str, optional
        Path to MATLAB input info file (for plotting raw data).
    pixel : int, optional
        Compare only this MATLAB pixel index. If None, compare all.
    rtol : float
        Relative tolerance for float comparisons.
    plot : bool
        Whether to plot mismatched pixels.
    """
    py = netCDF4.Dataset(python_file)
    try:
        mat = loadmat(matlab_file, squeeze_me=False)

        def get_cell(key, idx):
            cell = mat[key]
            val = cell[idx, 0]
            if val.size == 0:
                return np.array([])
            return val.flatten()

        # Build pixel coordinate mapping.
        # MATLAB iterates over all non-zero summary pixels in argwhere order.
        with netCDF4.Dataset(input_file) as nc:
            matlab_coords = np.argwhere(np.array(nc.variables["summary"][:, :]) != 0)

        # Determine which pixels to compare
        n_matlab_pixels = mat['cubicSpline_smPar'].shape[0]
        pixels = [pixel] if pixel is not None else range(n_matlab_pixels)

        # Summary counters
        total = 0
        matched = 0
        mismatched = 0
        python_only = 0
        matlab_only = 0

        for pix in pixels:
            # Map MATLAB pixel index to lat/lon indices via argwhere(summary != 0)
            lat_idx, lon_idx = matlab_coords[pix]

            m_sm = get_cell("cubicSpline_smPar", pix)
            has_matlab = len(m_sm) > 0 and not np.isnan(m_sm[0])

            # Check if Python has results for this pixel.
            # smoothing_parameter is NaN for unprocessed pixels (initialised to
            # all-NaN by init_phenology_output; set to the real value by
            # append_pixel_phenology for processed pixels).
            py_sm = float(np.ma.filled(py.variables['smoothing_parameter'][lat_idx, lon_idx], np.nan))
            has_python = not np.isnan(py_sm)

            if not has_matlab and not has_python:
                continue

            total += 1

            if has_python and not has_matlab:
                python_only += 1
                print(f"[PIX {pix}] ({lat_idx},{lon_idx}): Python has result, MATLAB does not")
                continue

            if has_matlab and not has_python:
                matlab_only += 1
                print(f"[PIX {pix}] ({lat_idx},{lon_idx}): MATLAB has result, Python does not")
                continue

            # Both have results — compare
            pixel_ok = True
            issues = []

            # Smoothing parameter
            if not np.isclose(m_sm[0], py_sm, rtol=rtol):
                pixel_ok = False
                issues.append(f"  Smoothing: MATLAB={m_sm[0]:.12e} Python={py_sm:.12e}")

            # Peaks
            m_pks_x = get_cell("phenoMETs_pks_X", pix)
            m_pks_y = get_cell("phenoMETs_pks_Y", pix)
            py_pks_x_unix = _get_padded(py, 'pks_x', lat_idx, lon_idx)
            py_pks_y = _get_padded(py, 'pks_y', lat_idx, lon_idx)
            py_pks_x = unix_to_datenum(py_pks_x_unix) if len(py_pks_x_unix) > 0 else np.array([])

            if len(m_pks_x) != len(py_pks_x):
                pixel_ok = False
                issues.append(f"  Peaks count: MATLAB={len(m_pks_x)} Python={len(py_pks_x)}")
            elif len(m_pks_x) > 0:
                if not np.allclose(m_pks_x, py_pks_x, rtol=rtol):
                    pixel_ok = False
                    issues.append(f"  Peaks X max diff: {np.max(np.abs(m_pks_x - py_pks_x)):.6f}")
                if not np.allclose(m_pks_y, py_pks_y, rtol=rtol):
                    pixel_ok = False
                    issues.append(f"  Peaks Y max diff: {np.max(np.abs(m_pks_y - py_pks_y)):.6e}")

            # Troughs
            m_trgs_x = get_cell("phenoMETs_trgs_X", pix)
            m_trgs_y = get_cell("phenoMETs_trgs_Y", pix)
            py_trgs_x_unix = _get_padded(py, 'trgs_x', lat_idx, lon_idx)
            py_trgs_y = _get_padded(py, 'trgs_y', lat_idx, lon_idx)
            py_trgs_x = unix_to_datenum(py_trgs_x_unix) if len(py_trgs_x_unix) > 0 else np.array([])

            if len(m_trgs_x) != len(py_trgs_x):
                pixel_ok = False
                issues.append(f"  Troughs count: MATLAB={len(m_trgs_x)} Python={len(py_trgs_x)}")
            elif len(m_trgs_x) > 0:
                if not np.allclose(m_trgs_x, py_trgs_x, rtol=rtol):
                    pixel_ok = False
                    issues.append(f"  Troughs X max diff: {np.max(np.abs(m_trgs_x - py_trgs_x)):.6f}")
                if not np.allclose(m_trgs_y, py_trgs_y, rtol=rtol):
                    pixel_ok = False
                    issues.append(f"  Troughs Y max diff: {np.max(np.abs(m_trgs_y - py_trgs_y)):.6e}")

            # Peaks QA
            m_pks_qa = get_cell("phenoMETs_pks_QA", pix)
            py_pks_qa = _get_padded(py, 'pks_qa', lat_idx, lon_idx, fill=-1)
            if len(m_pks_qa) == len(py_pks_qa) and len(m_pks_qa) > 0:
                if not np.array_equal(m_pks_qa.astype(int), py_pks_qa.astype(int)):
                    pixel_ok = False
                    mismatches = np.where(m_pks_qa.astype(int) != py_pks_qa.astype(int))[0]
                    for i in mismatches:
                        issues.append(f"  Peaks QA[{i}]: MATLAB={int(m_pks_qa[i])} Python={int(py_pks_qa[i])}")

            # Troughs QA
            m_trgs_qa = get_cell("phenoMETs_trgs_QA", pix)
            py_trgs_qa = _get_padded(py, 'trgs_qa', lat_idx, lon_idx, fill=-1)
            if len(m_trgs_qa) == len(py_trgs_qa) and len(m_trgs_qa) > 0:
                if not np.array_equal(m_trgs_qa.astype(int), py_trgs_qa.astype(int)):
                    pixel_ok = False
                    mismatches = np.where(m_trgs_qa.astype(int) != py_trgs_qa.astype(int))[0]
                    for i in mismatches:
                        issues.append(f"  Troughs QA[{i}]: MATLAB={int(m_trgs_qa[i])} Python={int(py_trgs_qa[i])} (MATLAB has known bug)")

            # Green-up metrics
            for name, suffix in [("Onset", "onset"), ("Mid", "mid"), ("Advanced", "advanced")]:
                m_x = get_cell(f"phenoMETs_GreenUp_{name}_X", pix)
                m_y = get_cell(f"phenoMETs_GreenUp_{name}_Y", pix)
                py_x_unix = _get_padded(py, f'green_up_{suffix}_x', lat_idx, lon_idx)
                py_y = _get_padded(py, f'green_up_{suffix}_y', lat_idx, lon_idx)
                py_x = unix_to_datenum(py_x_unix) if len(py_x_unix) > 0 else np.array([])

                if len(m_x) != len(py_x):
                    pixel_ok = False
                    issues.append(f"  GreenUp {name} count: MATLAB={len(m_x)} Python={len(py_x)}")
                elif len(m_x) > 0:
                    if not np.allclose(m_x, py_x, atol=1.0):
                        pixel_ok = False
                        issues.append(f"  GreenUp {name} X max diff: {np.max(np.abs(m_x - py_x)):.6f}")
                    if not np.allclose(m_y, py_y, rtol=rtol):
                        pixel_ok = False
                        issues.append(f"  GreenUp {name} Y max diff: {np.max(np.abs(m_y - py_y)):.6e}")

            # Green-down metrics
            for name, suffix in [("Onset", "onset"), ("Mid", "mid"), ("Advanced", "advanced")]:
                m_x = get_cell(f"phenoMETs_GreenDown_{name}_X", pix)
                m_y = get_cell(f"phenoMETs_GreenDown_{name}_Y", pix)
                py_x_unix = _get_padded(py, f'green_down_{suffix}_x', lat_idx, lon_idx)
                py_y = _get_padded(py, f'green_down_{suffix}_y', lat_idx, lon_idx)
                py_x = unix_to_datenum(py_x_unix) if len(py_x_unix) > 0 else np.array([])

                if len(m_x) != len(py_x):
                    pixel_ok = False
                    issues.append(f"  GreenDown {name} count: MATLAB={len(m_x)} Python={len(py_x)}")
                elif len(m_x) > 0:
                    if not np.allclose(m_x, py_x, atol=1.0):
                        pixel_ok = False
                        issues.append(f"  GreenDown {name} X max diff: {np.max(np.abs(m_x - py_x)):.6f}")
                    if not np.allclose(m_y, py_y, rtol=rtol):
                        pixel_ok = False
                        issues.append(f"  GreenDown {name} Y max diff: {np.max(np.abs(m_y - py_y)):.6e}")

            # Data gaps
            m_gap_start = get_cell("cubicSpline_dataGap_startDAYs", pix)
            m_gap_end = get_cell("cubicSpline_dataGap_endDAYs", pix)
            py_gap_start_unix = _get_padded(py, 'data_gap_start', lat_idx, lon_idx)
            py_gap_end_unix = _get_padded(py, 'data_gap_end', lat_idx, lon_idx)
            py_gap_start = unix_to_datenum(py_gap_start_unix) if len(py_gap_start_unix) > 0 else np.array([])
            py_gap_end = unix_to_datenum(py_gap_end_unix) if len(py_gap_end_unix) > 0 else np.array([])

            if len(m_gap_start) != len(py_gap_start):
                pixel_ok = False
                issues.append(f"  Data gaps count: MATLAB={len(m_gap_start)} Python={len(py_gap_start)}")
            elif len(m_gap_start) > 0:
                if not np.allclose(m_gap_start, py_gap_start):
                    pixel_ok = False
                    issues.append(f"  Gap start max diff: {np.max(np.abs(m_gap_start - py_gap_start)):.6f}")
                if not np.allclose(m_gap_end, py_gap_end):
                    pixel_ok = False
                    issues.append(f"  Gap end max diff: {np.max(np.abs(m_gap_end - py_gap_end)):.6f}")

            if pixel_ok:
                matched += 1
            else:
                mismatched += 1
                print(f"\n[PIX {pix}] ({lat_idx},{lon_idx}): MISMATCH")
                for issue in issues:
                    print(issue)

                if plot and matlab_data_file and matlab_info_file:
                    _plot_comparison(py, mat, pix, lat_idx, lon_idx,
                                     matlab_data_file, matlab_info_file,
                                     input_file=input_file, variable=variable)

        # Summary
        print(f"\n{'='*60}")
        print(f"Summary: {total} pixels compared")
        print(f"  Matched:     {matched}")
        print(f"  Mismatched:  {mismatched}")
        print(f"  Python only: {python_only}")
        print(f"  MATLAB only: {matlab_only}")
        print(f"{'='*60}")

    finally:
        py.close()


def _get_padded(nc, var_name, lat_idx, lon_idx, fill=np.nan):
    """Extract non-fill values from a (lat, lon, record) NetCDF variable.

    The new phenology format stores each feature type in a (lat, lon, record)
    variable where the record dimension is unlimited.  Unused record slots for
    a given pixel are left at the fill value (NaN for floats, -1 for QA ints).
    This helper strips those fill slots and returns only the real values.
    """
    if var_name not in nc.variables:
        return np.array([])
    data = np.ma.filled(nc.variables[var_name][lat_idx, lon_idx, :], fill)
    if np.isnan(fill):
        return data[~np.isnan(data)]
    return data[data != fill]


def _load_matlab_raw(mat, pix, matlab_data_file, matlab_info_file):
    """Load and filter raw MATLAB time-series data for a single pixel."""
    m_info = loadmat(matlab_info_file)
    m_data = loadmat(matlab_data_file)
    m_time_all = np.array(m_info["timeSeries_dateNR"][0])
    m_values_all = np.array(m_data["timeSeries_data"][pix])
    m_qa_all = np.array(m_data["timeSeries_QA"][pix])
    m_values_all[m_qa_all == 1] = np.nan
    m_time_filt = m_time_all[~np.isnan(m_values_all)]
    m_values_filt = m_values_all[~np.isnan(m_values_all)]
    m_smooth_x = mat["cubicSpline_smoothXaxis"].flatten()
    return m_smooth_x, m_time_filt, m_values_filt


def _plot_netcdf(py, lat_idx, lon_idx, input_file, variable, ax=None):
    """Plot Python NetCDF phenology results for a single pixel.

    Parameters
    ----------
    py : netCDF4.Dataset
        Open Python phenology output dataset.
    lat_idx, lon_idx : int
        Pixel coordinates.
    input_file : str
        Path to Python input NetCDF (for raw observations and spline evaluation).
    variable : str
        Variable name in input_file to plot as raw observations.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. If None, a new figure is created and shown.
    """
    from matplotlib import pyplot as plt
    from csaps import csaps

    standalone = ax is None
    if standalone:
        _, ax = plt.subplots(figsize=(16, 5))

    qa_var = getattr(py, 'qa', None)
    qa_filter = bool(getattr(py, 'qa_filter', 0))

    with netCDF4.Dataset(input_file) as nc:
        py_time_unix = np.array(nc.variables["time"][:], dtype=np.float64)
        py_values_all = np.array(nc.variables[variable][:, lat_idx, lon_idx])
        py_mask = py_values_all != -9999
        if qa_filter and qa_var and qa_var in nc.variables:
            qa = np.array(nc.variables[qa_var][:, lat_idx, lon_idx])
            py_mask = py_mask & (qa == 0)
        # Reproduce the exact arrays used by the processing code:
        # - t_all: integer-only datenum for all (unfiltered) timestamps
        # - smooth_x: daily integer grid over the FULL dataset time range
        # phenology.py builds smooth_x_axis once from all timestamps (line 57)
        # and passes it into every pixel's spline fit, so we must do the same.
        t_all = np.array([datetime.utcfromtimestamp(int(ts)).toordinal() + 366
                          for ts in py_time_unix])
    t = t_all[py_mask]
    py_values_filt = py_values_all[py_mask]
    smooth_x = np.arange(t_all.min(), t_all.max() + 1, 1)

    py_pks_x_unix = _get_padded(py, 'pks_x', lat_idx, lon_idx)
    py_pks_y = _get_padded(py, 'pks_y', lat_idx, lon_idx)
    py_trgs_x_unix = _get_padded(py, 'trgs_x', lat_idx, lon_idx)
    py_trgs_y = _get_padded(py, 'trgs_y', lat_idx, lon_idx)
    py_sm = float(np.ma.filled(py.variables['smoothing_parameter'][lat_idx, lon_idx], np.nan))
    # Peak/trough x positions are stored as unix; convert using the same integer
    # datenum formula so they map to exact grid points on smooth_x.
    _to_int_datenum = lambda ts_arr: np.array(
        [datetime.utcfromtimestamp(int(ts)).toordinal() + 366 for ts in ts_arr])
    py_pks_x = _to_int_datenum(py_pks_x_unix) if len(py_pks_x_unix) > 0 else np.array([])
    py_trgs_x = _to_int_datenum(py_trgs_x_unix) if len(py_trgs_x_unix) > 0 else np.array([])

    py_gu = {}
    py_gd = {}
    for suffix in ['onset', 'mid', 'advanced']:
        x_unix = _get_padded(py, f'green_up_{suffix}_x', lat_idx, lon_idx)
        py_gu[f'{suffix}_x'] = unix_to_datenum(x_unix) if len(x_unix) > 0 else np.array([])
        py_gu[f'{suffix}_y'] = _get_padded(py, f'green_up_{suffix}_y', lat_idx, lon_idx)
        x_unix = _get_padded(py, f'green_down_{suffix}_x', lat_idx, lon_idx)
        py_gd[f'{suffix}_x'] = unix_to_datenum(x_unix) if len(x_unix) > 0 else np.array([])
        py_gd[f'{suffix}_y'] = _get_padded(py, f'green_down_{suffix}_y', lat_idx, lon_idx)

    py_smooth_y = csaps(t, py_values_filt, smooth_x, smooth=py_sm)

    ax.scatter(t, py_values_filt, color='grey', alpha=0.3, s=10, label='Data')
    ax.plot(smooth_x, py_smooth_y, color='blue', linewidth=1, label='Python spline')
    if len(py_pks_x) > 0:
        ax.scatter(py_pks_x, py_pks_y, color='red', s=60, marker='^', zorder=4, label='Peaks')
    if len(py_trgs_x) > 0:
        ax.scatter(py_trgs_x, py_trgs_y, color='green', s=60, marker='v', zorder=4, label='Troughs')
    for suffix, color, alpha in [('onset', 'lime', 0.3), ('mid', 'lime', 0.5), ('advanced', 'lime', 0.7)]:
        for x_val in py_gu[f'{suffix}_x']:
            ax.axvline(x_val, color=color, alpha=alpha, linewidth=0.5)
    for suffix, color, alpha in [('onset', 'orange', 0.3), ('mid', 'orange', 0.5), ('advanced', 'orange', 0.7)]:
        for x_val in py_gd[f'{suffix}_x']:
            ax.axvline(x_val, color=color, alpha=alpha, linewidth=0.5)
    ax.set_title(f'Python (smoothing={py_sm:.12e})')
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(True, alpha=0.3)

    if standalone:
        plt.tight_layout()
        plt.show()


def _plot_matlab(mat, pix, matlab_data_file, matlab_info_file, ax=None):
    """Plot MATLAB phenology results for a single pixel.

    Parameters
    ----------
    mat : dict
        Loaded MATLAB .mat output.
    pix : int
        MATLAB pixel index.
    matlab_data_file : str
        Path to MATLAB input data file.
    matlab_info_file : str
        Path to MATLAB input info file.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. If None, a new figure is created and shown.
    """
    from matplotlib import pyplot as plt
    from csaps import csaps

    standalone = ax is None
    if standalone:
        _, ax = plt.subplots(figsize=(16, 5))

    m_smooth_x, m_time_filt, m_values_filt = _load_matlab_raw(mat, pix, matlab_data_file, matlab_info_file)

    def get_cell(key, idx):
        cell = mat[key]
        val = cell[idx, 0]
        if val.size == 0:
            return np.array([])
        return val.flatten()

    m_sm = get_cell("cubicSpline_smPar", pix)[0]
    m_pks_x = get_cell("phenoMETs_pks_X", pix)
    m_pks_y = get_cell("phenoMETs_pks_Y", pix)
    m_trgs_x = get_cell("phenoMETs_trgs_X", pix)
    m_trgs_y = get_cell("phenoMETs_trgs_Y", pix)

    m_gu = {}
    m_gd = {}
    for name, suffix in [("Onset", "onset"), ("Mid", "mid"), ("Advanced", "advanced")]:
        m_gu[f'{suffix}_x'] = get_cell(f"phenoMETs_GreenUp_{name}_X", pix)
        m_gu[f'{suffix}_y'] = get_cell(f"phenoMETs_GreenUp_{name}_Y", pix)
        m_gd[f'{suffix}_x'] = get_cell(f"phenoMETs_GreenDown_{name}_X", pix)
        m_gd[f'{suffix}_y'] = get_cell(f"phenoMETs_GreenDown_{name}_Y", pix)

    m_smooth_y = csaps(m_time_filt, m_values_filt, m_smooth_x, smooth=m_sm)

    ax.scatter(m_time_filt, m_values_filt, color='grey', alpha=0.3, s=10, label='Data')
    ax.plot(m_smooth_x, m_smooth_y, color='blue', linewidth=1, label='MATLAB spline')
    if len(m_pks_x) > 0:
        ax.scatter(m_pks_x, m_pks_y, color='red', s=60, marker='^', zorder=4, label='Peaks')
    if len(m_trgs_x) > 0:
        ax.scatter(m_trgs_x, m_trgs_y, color='green', s=60, marker='v', zorder=4, label='Troughs')
    for suffix, color, alpha in [('onset', 'lime', 0.3), ('mid', 'lime', 0.5), ('advanced', 'lime', 0.7)]:
        for x_val in m_gu[f'{suffix}_x']:
            ax.axvline(x_val, color=color, alpha=alpha, linewidth=0.5)
    for suffix, color, alpha in [('onset', 'orange', 0.3), ('mid', 'orange', 0.5), ('advanced', 'orange', 0.7)]:
        for x_val in m_gd[f'{suffix}_x']:
            ax.axvline(x_val, color=color, alpha=alpha, linewidth=0.5)
    ax.set_title(f'MATLAB (smoothing={m_sm:.12e})')
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(True, alpha=0.3)

    if standalone:
        plt.tight_layout()
        plt.show()


def _plot_comparison(py, mat, pix, lat_idx, lon_idx,
                     matlab_data_file, matlab_info_file, input_file=None, variable=None):
    """Plot Python vs MATLAB results for a single pixel."""
    from matplotlib import pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(16, 10), sharex=True, sharey=True)

    _plot_netcdf(py, lat_idx, lon_idx, input_file, variable, ax=axes[0])
    _plot_matlab(mat, pix, matlab_data_file, matlab_info_file, ax=axes[1])

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    python_file = "/home/runnalja/git/cci-cyano-production/data/v3.1/phenology/chla/327.nc"
    input_file = "/home/runnalja/git/cci-cyano-production/data/v3.1/extract/chla/327.nc"
    matlab_file = "/home/runnalja/git/bgb-jelle/data/CCISTEP3_phenoMETs_ss/Lake_297_327_GLWD00000327/pixFileDataList_chla_mean_731335_738886/cubicSpline_1_14_0_0_31_0_365_5_pixNRs_1_982.mat"
    matlab_data = "/home/runnalja/git/bgb-jelle/data/CCISTEP1_timeSeriesAndNullInfo/Lake_297_327_GLWD00000327/CCI_Data_chla_mean_731335_738886_pixNRs_1_982.mat"
    matlab_info = "/home/runnalja/git/bgb-jelle/data/CCISTEP1_timeSeriesAndNullInfo/Lake_297_327_GLWD00000327/CCI_NullInfo_chla_mean_731335_738886.mat"
    pixel = 600

    compare(python_file, matlab_file, input_file, matlab_data_file=matlab_data, matlab_info_file=matlab_info, pixel=pixel, plot=True, variable="chla")
