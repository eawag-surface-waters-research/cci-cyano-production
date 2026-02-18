import os
import netCDF4
import logging
import numpy as np
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import functions

def phenology(lake, p, threads=1):
    """Extract phenology metrics for a single lake from time series data.

    Reads extracted time series, applies a smoothing cubic spline, and extracts
    phenology metrics (green-up/green-down onset, mid, advanced) for each pixel.
    Results are written to a NetCDF file incrementally as each pixel is processed.

    Parameters
    ----------
    lake : dict or GeoDataFrame row
        Lake metadata containing at minimum 'id'.
    p : dict
        Parameters dict with keys:
        - out_folder : str - base output directory
        - variable : str - variable name (e.g. 'chla_mean')
        - qa : str - QA variable name
        - qa_filter : bool - whether to filter on QA flag
        - spline_min_phase_length : int - minimum phase length in days
        - spline_min_relative_amplitude : float - minimum relative amplitude (0-1)
        - spline_min_phase_data : int - minimum data points per phase
        - spline_subs_peak_win_size : int - window size for substantial peak check (days)
        - spline_subs_peak_ampl_frac : float - amplitude fraction for substantial peak check
        - spline_data_gap_size : int - minimum gap size to flag (days)
        - spline_data_gap_size_buffer : int - buffer around data gaps (days)
        - parallel : str - 'lakes' (parallelise across lakes) or 'pixels' (parallelise
          pixels within each lake); defaults to 'lakes'
        - threads : int - number of parallel workers
    """
    input_file = os.path.join(p["out_folder"], "extract", p["variable"], "{}.nc".format(lake["id"]))
    output_file = os.path.join(p["out_folder"], "phenology", p["variable"], "{}.nc".format(lake["id"]))

    if os.path.isfile(output_file):
        logging.info(f"File {lake['id']}.nc already exists")
        return

    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    if os.path.isfile(input_file):
        logging.info(f"Starting phenology computation for {lake['id']}")
    else:
        logging.info(f"Input file not available for {lake['id']}")
        return

    with netCDF4.Dataset(input_file) as nc:
        valid_coords = np.argwhere(nc.variables["summary"][:, :] > 2)
        lat = np.array(nc.variables["lat"])
        lon = np.array(nc.variables["lon"])
        t = np.array([datetime.utcfromtimestamp(int(ts)).toordinal() + 366 for ts in np.array(nc.variables["time"])])

    smooth_x_axis = np.arange(t.min(), t.max() + 1, 1)
    out = None

    if threads > 1:
        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = {executor.submit(compute_pixel, input_file, x, y, t, smooth_x_axis, p): None
                       for x, y in valid_coords}
            for future in tqdm(as_completed(futures), total=len(futures), desc=str(lake['id']), unit='px'):
                x, y, result, pheno = future.result()
                if result is not None:
                    if out is None:
                        out = netCDF4.Dataset(output_file, 'w', format='NETCDF4')
                        functions.init_phenology_output(out, lat, lon)
                    functions.append_pixel_phenology(out, x, y, result, pheno)
    else:
        for x, y in tqdm(valid_coords, desc=str(lake['id']), unit='px'):
            x, y, result, pheno = compute_pixel(input_file, x, y, t, smooth_x_axis, p)
            if result is not None:
                if out is None:
                    out = netCDF4.Dataset(output_file, 'w', format='NETCDF4')
                    functions.init_phenology_output(out, lat, lon)
                functions.append_pixel_phenology(out, x, y, result, pheno)

    if out is not None:
        out.close()
    else:
        logging.info(f"No valid phenology results for {lake['id']}")

def compute_pixel(input_file, x, y, t, smooth_x_axis, p):
    """Compute spline and phenology metrics for a single pixel.

    Opens the input file, extracts the time series for pixel (x, y), fits the
    smoothing cubic spline, and extracts phenology metrics.  Defined at module
    level so it can be pickled by ProcessPoolExecutor.

    Parameters
    ----------
    input_file : str
        Path to the extract NetCDF file for the lake.
    x : int
        Lat index of the pixel.
    y : int
        Lon index of the pixel.
    t : ndarray
        Ordinal time values (days since Jan 0, year 0) for all observations.
    smooth_x_axis : ndarray
        Regular daily grid used for spline evaluation.
    p : dict
        Parameters dict (see phenology docstring).

    Returns
    -------
    tuple
        (x, y, result, pheno) where result and pheno are None if the spline
        conditions are not satisfied for this pixel.
    """
    with netCDF4.Dataset(input_file) as nc:
        values = np.array(nc.variables[p["variable"]][:, x, y])
        mask = values != -9999
        if p["qa_filter"]:
            qa = np.array(nc.variables[p["qa"]][:, x, y])
            mask = (qa == 0) & mask
        values = values[mask]
        time = t[mask]

    result = functions.smooth_cubic_spline(
        time, values, p["spline_min_phase_length"], p["spline_min_relative_amplitude"],
        p["spline_min_phase_data"], smoothing_change=1e-12, max_iterations=1e5,
        smooth_x_axis=smooth_x_axis
    )
    if not result['conditions_satisfied']:
        return x, y, None, None

    pheno = functions.extract_phenology_metrics(
        result['x_pks'], result['y_pks'], result['x_trgs'], result['y_trgs'],
        time, result['smooth_x_axis'], result['smooth_y_data'],
        p["spline_subs_peak_win_size"], p["spline_subs_peak_ampl_frac"],
        p["spline_data_gap_size"], p["spline_data_gap_size_buffer"]
    )
    return x, y, result, pheno
