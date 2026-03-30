import os
import netCDF4
import logging
import numpy as np
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import functions

def phenology(lake, p, threads=1, batch_size=100):
    """Extract phenology metrics for a single lake from time series data.

    Reads extracted time series, applies a smoothing cubic spline, and extracts
    phenology metrics (green-up/green-down onset, mid, advanced) for each pixel.
    Pixels are processed in batches to reduce I/O: the input file is opened once
    per batch and the output file is opened at most once per lake.

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
    batch_size : int
        Number of pixels to process per batch.  Each batch opens the input file
        once.  With threads > 1 each worker receives one batch.
    """
    input_file = os.path.join(p["out_folder"], "extract", p["variable"], "{}.nc".format(lake["id"]))
    output_file = os.path.join(p["out_folder"], "phenology", p["variable"], "{}.nc".format(lake["id"]))

    if os.path.isfile(output_file):
        logging.info(f"File {lake['id']}.nc already exists")
        return

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    output_file_temp = output_file.replace(".nc", "_tmp.nc")

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
    batches = [valid_coords[i:i + batch_size] for i in range(0, len(valid_coords), batch_size)]
    out = None

    if threads > 1:
        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = {executor.submit(compute_pixel_batch, input_file, batch, t, smooth_x_axis, p): None
                       for batch in batches}
            for future in tqdm(as_completed(futures), total=len(futures), desc=str(lake['id']), unit='batch'):
                for x, y, result, pheno in future.result():
                    if result is not None:
                        if out is None:
                            out = netCDF4.Dataset(output_file_temp, 'w', format='NETCDF4')
                            functions.init_phenology_output(out, lat, lon, p=p)
                        functions.append_pixel_phenology(out, x, y, result, pheno)
    else:
        for batch in tqdm(batches, desc=str(lake['id']), unit='batch'):
            for x, y, result, pheno in compute_pixel_batch(input_file, batch, t, smooth_x_axis, p):
                if result is not None:
                    if out is None:
                        out = netCDF4.Dataset(output_file_temp, 'w', format='NETCDF4')
                        functions.init_phenology_output(out, lat, lon, p=p)
                    functions.append_pixel_phenology(out, x, y, result, pheno)

    if out is not None:
        out.close()
        os.rename(output_file_temp, output_file)
        logging.info(f"Completed phenology computation for {lake['id']}")
    else:
        logging.info(f"No valid phenology results for {lake['id']}")

def compute_pixel_batch(input_file, coords, t, smooth_x_axis, p):
    """Compute spline and phenology metrics for a batch of pixels.

    Opens the input file once for the entire batch, reads all pixel time series,
    then fits splines and extracts metrics.  Defined at module level so it can
    be pickled by ProcessPoolExecutor.

    Parameters
    ----------
    input_file : str
        Path to the extract NetCDF file for the lake.
    coords : array-like of (int, int)
        Sequence of (lat_index, lon_index) pairs to process.
    t : ndarray
        Ordinal time values (days since Jan 0, year 0) for all observations.
    smooth_x_axis : ndarray
        Regular daily grid used for spline evaluation.
    p : dict
        Parameters dict (see phenology docstring).

    Returns
    -------
    list of tuple
        Each element is (x, y, result, pheno) where result and pheno are None
        if the spline conditions are not satisfied for that pixel.
    """
    pixel_data = []
    with netCDF4.Dataset(input_file) as nc:
        for x, y in coords:
            values = np.array(nc.variables[p["variable"]][:, x, y])
            mask = values != -9999
            if p["qa_filter"]:
                qa = np.array(nc.variables[p["qa"]][:, x, y])
                mask = (qa == 0) & mask
            pixel_data.append((x, y, values[mask], t[mask]))

    batch_results = []
    for x, y, values, time in pixel_data:
        if len(time) < 2:
            batch_results.append((x, y, None, None))
            continue
        try:
            result = functions.smooth_cubic_spline(
                time, values, p["spline_min_phase_length"], p["spline_min_relative_amplitude"],
                p["spline_min_phase_data"], smoothing_change=1e-12, max_iterations=1e5,
                smooth_x_axis=smooth_x_axis
            )
            if not result['conditions_satisfied']:
                batch_results.append((x, y, None, None))
                continue
            pheno = functions.extract_phenology_metrics(
                result['x_pks'], result['y_pks'], result['x_trgs'], result['y_trgs'],
                time, result['smooth_x_axis'], result['smooth_y_data'],
                p["spline_subs_peak_win_size"], p["spline_subs_peak_ampl_frac"],
                p["spline_data_gap_size"], p["spline_data_gap_size_buffer"]
            )
            batch_results.append((x, y, result, pheno))
        except Exception as e:
            logging.warning("Pixel (%d, %d) failed: %s", x, y, e)
            batch_results.append((x, y, None, None))

    return batch_results
