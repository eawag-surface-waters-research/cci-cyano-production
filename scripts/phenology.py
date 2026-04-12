import os
import sys
import glob
import shutil
import signal
import faulthandler
import netCDF4
import logging
import numpy as np
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import functions

# Module-level globals set in the main process before forking workers.
# Forked workers inherit these via copy-on-write without pickling.
_VAR_DATA = None
_QA_DATA = None

def _worker_init():
    """Initialiser run once in each worker process at startup."""
    faulthandler.enable()

    def _sigterm_handler(signum, frame):
        print(f"Worker {os.getpid()} received SIGTERM", flush=True)
        sys.exit(1)

    signal.signal(signal.SIGTERM, _sigterm_handler)


def _write_from_checkpoints(checkpoint_dir, output_file_temp, output_file, lat, lon, p, lake_id):
    """Write the final netCDF4 output from saved batch checkpoints.

    Reads each .npy checkpoint in order and writes all valid pixel results
    to the output file in a single sequential pass.

    Parameters
    ----------
    checkpoint_dir : str
        Directory containing batch_NNNN.npy checkpoint files.
    output_file_temp : str
        Temporary output path (renamed to output_file on success).
    output_file : str
        Final output netCDF4 path.
    lat, lon : ndarray
        Coordinate arrays for initialising the output file.
    p : dict
        Parameters dict.
    lake_id : str or int
        Lake identifier used for tqdm label.
    """
    checkpoint_files = sorted(glob.glob(os.path.join(checkpoint_dir, 'batch_*.npy')))
    out = None
    for cp_file in tqdm(checkpoint_files, desc=f"{lake_id} writing", unit='batch'):
        for x, y, smoothing, pheno in np.load(cp_file, allow_pickle=True):
            if out is None:
                out = netCDF4.Dataset(output_file_temp, 'w', format='NETCDF4')
                functions.init_phenology_output(out, lat, lon, p=p)
            functions.append_pixel_phenology(out, x, y, {'smoothing': smoothing}, pheno)
    if out is not None:
        out.close()
        os.rename(output_file_temp, output_file)
        logging.info(f"Completed phenology computation for {lake_id}")
    else:
        logging.info(f"No valid phenology results for {lake_id}")


def phenology(lake, p, threads=1, batch_size=100):
    """Extract phenology metrics for a single lake from time series data.

    Reads extracted time series, applies a smoothing cubic spline, and extracts
    phenology metrics (green-up/green-down onset, mid, advanced) for each pixel.
    Batch results are saved as .npy checkpoints so the run can be resumed if
    interrupted. The final netCDF4 is written in one sequential pass after all
    batches complete.

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
    global _VAR_DATA, _QA_DATA

    input_file = os.path.join(p["out_folder"], "extract", p["variable"], "{}.nc".format(lake["id"]))
    output_file = os.path.join(p["out_folder"], "phenology", p["variable"], "{}.nc".format(lake["id"]))
    checkpoint_dir = os.path.join(p["out_folder"], "phenology", p["variable"], "checkpoints", str(lake["id"]), f"bs{batch_size}")

    if os.path.isfile(output_file):
        logging.info(f"File {lake['id']}.nc already exists")
        return

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    os.makedirs(checkpoint_dir, exist_ok=True)
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
        _VAR_DATA = np.array(nc.variables[p["variable"]])
        _QA_DATA = np.array(nc.variables[p["qa"]]) if p["qa_filter"] else None

    smooth_x_axis = np.arange(t.min(), t.max() + 1, 1)

    coord_list = [(int(x), int(y)) for x, y in valid_coords]
    batches = [coord_list[i:i + batch_size] for i in range(0, len(coord_list), batch_size)]

    # Find already completed batches for restart
    completed = {int(os.path.splitext(os.path.basename(f))[0].split('_')[1])
                 for f in glob.glob(os.path.join(checkpoint_dir, 'batch_*.npy'))}
    logging.info(f"Lake {lake['id']}: {len(completed)}/{len(batches)} batches already complete")

    if threads > 1:
        with ProcessPoolExecutor(max_workers=threads, initializer=_worker_init) as executor:
            futures = {executor.submit(compute_pixel_batch, batch, smooth_x_axis, p, t): i
                       for i, batch in enumerate(batches) if i not in completed}
            for future in tqdm(as_completed(futures), total=len(batches),
                               initial=len(completed), desc=str(lake['id']), unit='batch'):
                batch_idx = futures[future]
                try:
                    batch_result = future.result()
                except Exception as e:
                    print(f"Worker failed on batch {batch_idx}: {e}", flush=True)
                    continue
                checkpoint = [(x, y, result['smoothing'], pheno)
                              for x, y, result, pheno in batch_result if result is not None]
                np.save(os.path.join(checkpoint_dir, f'batch_{batch_idx:04d}.npy'),
                        checkpoint, allow_pickle=True)
    else:
        for i, batch in enumerate(tqdm(batches, desc=str(lake['id']), unit='batch')):
            if i in completed:
                continue
            batch_result = compute_pixel_batch(batch, smooth_x_axis, p, t)
            checkpoint = [(x, y, result['smoothing'], pheno)
                          for x, y, result, pheno in batch_result if result is not None]
            np.save(os.path.join(checkpoint_dir, f'batch_{i:04d}.npy'),
                    checkpoint, allow_pickle=True)

    _VAR_DATA = None
    _QA_DATA = None

    _write_from_checkpoints(checkpoint_dir, output_file_temp, output_file, lat, lon, p, lake["id"])
    shutil.rmtree(checkpoint_dir)


def compute_pixel_batch(coord_batch, smooth_x_axis, p, t):
    """Compute spline and phenology metrics for a batch of pixels.

    Reads pixel data from the module-level _VAR_DATA and _QA_DATA arrays
    inherited from the parent process via fork copy-on-write, then performs
    spline fitting and phenology extraction. Defined at module level so it
    can be pickled by ProcessPoolExecutor.

    Parameters
    ----------
    coord_batch : list of (x, y)
        Pixel coordinates to process.
    smooth_x_axis : ndarray
        Regular daily grid used for spline evaluation.
    p : dict
        Parameters dict (see phenology docstring).
    t : ndarray
        Pre-computed time ordinal array for all time steps.

    Returns
    -------
    list of tuple
        Each element is (x, y, result, pheno) where result and pheno are None
        if the spline conditions are not satisfied for that pixel.
    """
    batch_results = []
    for x, y in coord_batch:
        values = _VAR_DATA[:, x, y]
        mask = values != -9999
        if p["qa_filter"]:
            mask = (_QA_DATA[:, x, y] == 0) & mask
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
