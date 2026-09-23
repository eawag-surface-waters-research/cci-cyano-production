import os
from pathlib import Path
import logging
import netCDF4
import numpy as np
import pandas as pd
import datetime
import multiprocessing
from functools import partial
from shapely.geometry import Point
import warnings
from contextlib import contextmanager
import functions as f
from pixel_metric_calcs.base import PixelCalcBase

_GLOBALS = {}

REQ_COLS = ["primary", "qa_column", "secondary", "i", "j", "year.DOY"]


DTYPES = {
    "primary": "i1",
    "secondary": "i1",
    "i": "i4",
    "j": "i4",
    "qa_column": "i4",
    "year.DOY": "f8",
}

def _init_bloom_assemble_worker(p_path, var_names):
    """Initialise per-process globals for bloom pixel extraction.

    Called once per worker process by multiprocessing.Pool. Preloads all
    required phenology arrays as numpy arrays so per-pixel work is pure
    in-memory indexing with no NetCDF I/O in the hot path.

    Parameters
    ----------
    p_path : str
        Path to the phenology NetCDF file.
    var_names : list of str
        NetCDF variable names to preload (determined by assemble_bloom_data).
    """
    with netCDF4.Dataset(p_path) as nc:
        for var in var_names:
            v = nc.variables[var]
            nlat, nlon, nrec = v.shape
            arr = np.empty((nlat, nlon, nrec), dtype=v.dtype)
            # read pixel by pixel (v[i, j, :]) rather than a bulk v[:] read - netCDF4
            # 1.7.4 silently misattributes data between pixels when read this way for
            # files with an unlimited 'record' dimension, verified against the trusted
            # per-pixel access pattern used elsewhere in this class (e.g. _load_pixel_data)
            for i in range(nlat):
                for j in range(nlon):
                    arr[i, j, :] = v[i, j, :]
            _GLOBALS[f"bloom_{var}"] = arr


class BloomAssemble(PixelCalcBase):
    """Everything here used to be called 'KDE calculations'.
    If you run into a 'KDE' reference somewhere, 
    it is (likely) related to bloom assembly."""
    def __init__(self, *args, **kwargs):

        kwargs = kwargs | {'metric_name': "bloom_assembly"}
        super().__init__(*args, **kwargs)
        
    def build_path(self):
        """Return the bloom cache directory and file path."""
        out_dir = os.path.join(
            self.out_folder,
            "calculated_values",
            "bloom_assembly_data",
            self.variable,
        )
        os.makedirs(out_dir, exist_ok=True)

        filename = f"ID{self.lakeID}_{self.start_year}_{self.end_year}.nc"
        self.save_fp = os.path.join(out_dir, filename)

    def read_input(self):
        """Resolve the source files used to calculate bloom events."""
        self.p_path = Path(self.p_path).resolve()
        self.e_path = Path(self.e_path).resolve()

        if not self.p_path.is_file():
            raise FileNotFoundError(f"Phenology file not found: {self.p_path}")

        if not self.e_path.is_file():
            raise FileNotFoundError(f"Extract file not found: {self.e_path}")

        return self.p_path, self.e_path

    def calculate(
        self,
        primary_vars=None,
        secondary_vars=None,
        qa_var="pks",
    ):
        """Collect bloom events for all valid pixels."""
        self.bloom_assembly_df = self.assemble_bloom_data(
            primary_vars=primary_vars,
            secondary_vars=secondary_vars,
            qa_var=qa_var,
        )
        return self.bloom_assembly_df

    def _create_output(self, ds, n_events=None, block_size=1000):
        """
        Create bloom event output schema.

        Parameters
        ----------
        ds : netCDF4.Dataset
            Open output dataset.
        n_events : int or None
            Number of events. Use None for unlimited dimension.
        block_size : int
            Chunk size along event dimension.
        """

        ds.createDimension("event", n_events)

        output_vars = {}

        for column in REQ_COLS:

            chunksizes = None
            if n_events != 0:
                chunksizes = (block_size,)

            output_vars[column] = ds.createVariable(
                column,
                DTYPES[column],
                ("event",),
                zlib=True,
                complevel=4,
                chunksizes=chunksizes,
            )

        ds.lake_id = str(self.lakeID)
        ds.version = str(self.version)
        ds.variable = str(self.variable)
        ds.metric_name = self.metric_name
        ds.start_year = self.start_year
        ds.end_year = self.end_year

        return output_vars

    
    def write_output(self):
        """Write calculated bloom events to the cache."""
        if not hasattr(self, "bloom_assembly_df") or self.bloom_assembly_df is None:
            raise ValueError("No bloom data available to write.")

        with netCDF4.Dataset(self.save_fp, "w") as ds:
            output_vars = self._create_output(ds,
                                              n_events=len(self.bloom_assembly_df),
                                              block_size=len(self.bloom_assembly_df))

            for column in REQ_COLS:
                output_vars[column][:] = (
                    self.bloom_assembly_df[column].to_numpy()
                )

    def calculate_and_write_chunked(
        self,
        primary_vars=None,
        secondary_vars=None,
        qa_var="pks",
        pixel_block_size=1000,
        block_size = 1000
    ):
        """Calculate and write bloom events incrementally for large lakes."""

        with netCDF4.Dataset(self.save_fp, "w") as ds:
            output_vars = self._create_output(
                ds,
                n_events = None,
                block_size = block_size
            )

            event_start = 0

            for start in range(0, len(self.valid_coords), pixel_block_size):
                stop = min(
                    start + pixel_block_size,
                    len(self.valid_coords),
                )

                coords = self.valid_coords[start:stop]

                bloom_assembly_df = self.assemble_bloom_data(
                    primary_vars=primary_vars,
                    secondary_vars=secondary_vars,
                    qa_var=qa_var,
                    coords=coords,
                )

                if bloom_assembly_df.empty:
                    continue

                event_stop = event_start + len(bloom_assembly_df)

                for column in REQ_COLS:
                    output_vars[column][event_start:event_stop] = (
                        bloom_assembly_df[column].to_numpy()
                    )

                event_start = event_stop


    def read_cached_metric(self):
        """Read cached bloom event data from NetCDF."""
        with netCDF4.Dataset(self.save_fp, "r") as ds:
            required = set(REQ_COLS)
            missing = required.difference(ds.variables)

            if missing:
                raise ValueError(
                    f"Bloom cache {self.save_fp} is missing variables: "
                    f"{sorted(missing)}"
                )

            event_df = pd.DataFrame({col : np.asarray(ds.variables[col][:]) for col in required})
            return event_df[REQ_COLS]

    @staticmethod
    def _extract_pixel_bloom_events(nc, i, j,
                                   primary_vars=None, secondary_vars=None,
                                   qa_var="pks", arrays=None):
        """Extract bracketing and peak events for a single pixel.

        Parameters
        ----------
        nc : netCDF4.Dataset or None
            Open phenology dataset. Ignored when arrays is provided.
        i, j : int
            Pixel row and column indices.
        primary_vars : list of str, optional
            NetCDF variable names whose events mark the start of a bloom bracket
            (e.g. green-up). Defaults to ['green_up_advanced'].
        secondary_vars : list of str, optional
            NetCDF variable names whose events mark the end of a bloom bracket
            (e.g. green-down onset). Defaults to ['green_down_onset'].
        qa_var : str, optional
            Short name for the peak QA variable (passed through parse_qa_var_from_str).
            Defaults to 'pks'.
        arrays : dict of str -> np.ndarray, optional
            Preloaded full arrays keyed by NetCDF variable name. When provided,
            pixel slices are taken from these in-memory arrays instead of reading
            from nc, which is the fast path used by _bloom_pixel_worker.

        Returns
        -------
        pandas.DataFrame
            Columns: year.DOY, i, j, green_up_advanced, peaks, green_down_onset.
            Empty DataFrame if the pixel has no events.
        """
        if primary_vars is None:
            primary_vars = ["green_up_advanced"]
        if secondary_vars is None:
            secondary_vars = ["green_down_onset"]

        def _get(var_name):
            if arrays is not None:
                return arrays[var_name][i, j, :]
            return nc.variables[var_name][i, j, :]

        frames = []

        for var in primary_vars:
            var_x = f.coerce_varname_to_var_x(var)
            raw = f.remove_nan(_get(var_x))
            if len(raw) == 0:
                continue
            dt = pd.to_datetime(raw, unit="s", utc=True)
            frames.append(pd.DataFrame({
                "year.DOY": dt.year + dt.day_of_year / 1000,
                "i": i, "j": j,
                "primary": True,
                "qa_column": np.nan,
                "secondary": False,
            }))

        qa_var_parsed = f.parse_qa_var_from_str(qa_var)
        if qa_var_parsed is not None:
            qa_x = f.coerce_varname_to_var_x(qa_var_parsed)
            qa_x_raw = np.array(_get(qa_x))
            pk_mask = ~np.isnan(qa_x_raw)
            if pk_mask.any():
                pks_dt = pd.to_datetime(qa_x_raw[pk_mask], unit="s", utc=True)
                frames.append(pd.DataFrame({
                    "year.DOY": pks_dt.year + pks_dt.day_of_year / 1000,
                    "i": i, "j": j,
                    "primary": False,
                    "qa_column": np.array(_get(qa_var_parsed))[pk_mask].astype(int),
                    "secondary": False,
                }))

        for var in secondary_vars:
            var_x = f.coerce_varname_to_var_x(var)
            raw = f.remove_nan(_get(var_x))
            if len(raw) == 0:
                continue
            dt = pd.to_datetime(raw, unit="s", utc=True)
            frames.append(pd.DataFrame({
                "year.DOY": dt.year + dt.day_of_year / 1000,
                "i": i, "j": j,
                "primary": False,
                "qa_column": np.nan,
                "secondary": True,
            }))

        if not frames:
            return pd.DataFrame(columns=REQ_COLS)
        return pd.concat(frames, ignore_index=True)

    @staticmethod
    def _bloom_pixel_worker(coord, primary_vars, secondary_vars, qa_var):
        """Multiprocessing worker for a single pixel's bloom event extraction.

        Uses preloaded numpy arrays from _GLOBALS (populated by _init_bloom_worker)
        for pure in-memory indexing — no NetCDF I/O in the hot path.

        Parameters
        ----------
        coord : tuple of int
            (i, j) grid index pair.
        primary_vars, secondary_vars : list of str
            Passed through to the extraction logic.
        qa_var : str
            Short name for the peak QA variable.

        Returns
        -------
        pandas.DataFrame
            Same schema as _extract_pixel_bloom_events; empty if no events found.
        """
        i, j = coord
        arrays = {k[4:]: v for k, v in _GLOBALS.items() if k.startswith("bloom_")}
        return BloomAssemble._extract_pixel_bloom_events(
            nc=None, i=i, j=j, primary_vars = primary_vars, secondary_vars = secondary_vars, qa_var = qa_var, arrays=arrays
        )


    def assemble_bloom_data(self, primary_vars=None, secondary_vars=None, qa_var="pks",coords=None):
        """Collect bracketing and peak events lake-wide into a DataFrame.

        Iterates all valid pixels inside the 1 km-inset lake boundary and gathers
        events for each variable type. Each row represents one event occurrence.

        Parameters
        ----------
        primary_vars : list of str, optional
            Variables marking the start of a bloom bracket. Defaults to ['green_up_advanced'].
        secondary_vars : list of str, optional
            Variables marking the end of a bloom bracket. Defaults to ['green_down_onset'].
        qa_var : str, optional
            Short name for the peak QA variable. Defaults to 'pks'.

        Returns
        -------
        pandas.DataFrame
            Index named 'year.DOY' — a float of the form year + DOY/1000
            (e.g. 2005.150 = year 2005, day-of-year 150).
            Columns: i, j, primary (bool), qa_column (float), secondary (bool).
            Row count equals the total number of events across all variables and pixels.
        """

        print(f"assemble bloom events started at: {datetime.datetime.now()}")
        g = self._load_extracted_globals()
        lats = g["lat"]
        lons = g["lon"]

        if coords is None:
            coords = self.valid_coords

        inset_coords = [
            (i, j) for (i, j) in coords
            if self.prepped_geom.contains(Point(lons[j], lats[i]))
        ]

        # Compute the full set of NetCDF variable names needed so the initializer
        # can preload them as numpy arrays — eliminates per-pixel disk reads.
        var_names_set = set()
        _pv = primary_vars if primary_vars is not None else ["green_up_advanced"]
        _sv = secondary_vars if secondary_vars is not None else ["green_down_onset"]
        for var in _pv:
            var_names_set.add(f.coerce_varname_to_var_x(var))
        qa_var_parsed = f.parse_qa_var_from_str(qa_var)
        if qa_var_parsed is not None:
            var_names_set.add(f.coerce_varname_to_var_x(qa_var_parsed))
            var_names_set.add(qa_var_parsed)
        for var in _sv:
            var_names_set.add(f.coerce_varname_to_var_x(var))
        var_names = list(var_names_set)

        worker = partial(
            BloomAssemble._bloom_pixel_worker,
            primary_vars=primary_vars, secondary_vars=secondary_vars, qa_var=qa_var,
        )
        n_workers = min(10, os.cpu_count() or 4)
        chunksize = max(1, len(inset_coords) // (n_workers * 4))
        with multiprocessing.Pool(
            initializer=_init_bloom_assemble_worker, initargs=(self.p_path, var_names), processes=n_workers
        ) as pool:
            results = list(pool.imap_unordered(worker, inset_coords, chunksize=chunksize))

        frames = [df for df in results if not df.empty]
        if not frames:
            return pd.DataFrame(columns=REQ_COLS)
        print(f"assemble bloom events finished at: {datetime.datetime.now()}")
        return pd.concat(frames, ignore_index=True)