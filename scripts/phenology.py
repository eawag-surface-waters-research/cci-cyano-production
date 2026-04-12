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
    Only metadata is read upfront; workers open the input file themselves and
    read their own pixel data to avoid large memory allocations before forking.

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
        Number of pixels per batch dispatched to each worker.
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

    coord_list = [(int(x), int(y)) for x, y in valid_coords]
    batches = [coord_list[i:i + batch_size] for i in range(0, len(coord_list), batch_size)]
    out = None

    if threads > 1:
        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = {executor.submit(compute_pixel_batch, batch, smooth_x_axis, p, input_file, t): None
                       for batch in batches}
            for future in tqdm(as_completed(futures), total=len(futures), desc=str(lake['id']), unit='batch'):
                try:
                    batch_result = future.result()
                except Exception as e:
                    print(f"Worker failed: {e}", flush=True)
                    continue
                for x, y, result, pheno in batch_result:
                    if result is not None:
                        if out is None:
                            out = netCDF4.Dataset(output_file_temp, 'w', format='NETCDF4')
                            functions.init_phenology_output(out, lat, lon, p=p)
                        functions.append_pixel_phenology(out, x, y, result, pheno)
    else:
        for batch in tqdm(batches, desc=str(lake['id']), unit='batch'):
            for x, y, result, pheno in compute_pixel_batch(batch, smooth_x_axis, p, input_file, t):
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

def compute_pixel_batch(coord_batch, smooth_x_axis, p, input_file, t):
    """Compute spline and phenology metrics for a batch of pixels.

    Opens the input NetCDF4 file read-only, reads only the pixel columns
    needed for this batch, then performs spline fitting and phenology
    extraction.  Defined at module level so it can be pickled by
    ProcessPoolExecutor.

    Parameters
    ----------
    coord_batch : list of (x, y)
        Pixel coordinates to process.
    smooth_x_axis : ndarray
        Regular daily grid used for spline evaluation.
    p : dict
        Parameters dict (see phenology docstring).
    input_file : str
        Path to the input NetCDF4 file (opened read-only).
    t : ndarray
        Pre-computed time ordinal array for all time steps.

    Returns
    -------
    list of tuple
        Each element is (x, y, result, pheno) where result and pheno are None
        if the spline conditions are not satisfied for that pixel.
    """
    batch_results = []
    xs = [xy[0] for xy in coord_batch]
    ys = [xy[1] for xy in coord_batch]
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    with netCDF4.Dataset(input_file, 'r') as nc:
        var_slice = np.array(nc.variables[p["variable"]][:, x_min:x_max + 1, y_min:y_max + 1])
        qa_slice = np.array(nc.variables[p["qa"]][:, x_min:x_max + 1, y_min:y_max + 1]) if p["qa_filter"] else None
    for x, y in coord_batch:
        values = var_slice[:, x - x_min, y - y_min]
        mask = values != -9999
        if p["qa_filter"]:
            mask = (qa_slice[:, x - x_min, y - y_min] == 0) & mask
        values = values[mask]
        time = t[mask]

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
