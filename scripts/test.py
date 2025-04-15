import os
import time
import glob
import netCDF4
import numpy as np
import xarray as xr
import geopandas as gpd
from rasterio import features, transform
from concurrent.futures import ProcessPoolExecutor
import matplotlib.pyplot as plt

def process_lake(args):
    start = time.time()

    print(f"Starting {args['lake']['id']}.")
    output_path = os.path.join(args["out_folder"], "{}.nc".format(args["lake"]["id"]))
    output_temp = os.path.join(args["out_folder"], "{}_tmp.nc".format(args["lake"]["id"]))
    if os.path.isfile(output_path):
        print(f"File {args['lake']['id']}.nc already exists")
        return
    if os.path.isfile(output_temp):
        os.remove(output_temp)

    with xr.open_dataset(args["files"][0]) as ds:
        lat = ds["lat"].values
        lon = ds["lon"].values

    geom = args["lake"].geometry
    minx, miny, maxx, maxy = geom.bounds
    lat_mask = (lat >= miny) & (lat <= maxy)
    lon_mask = (lon >= minx) & (lon <= maxx)

    lat_indices = np.where(lat_mask)[0]
    lon_indices = np.where(lon_mask)[0]

    lat_sub = lat[lat_mask]
    lon_sub = lon[lon_mask]
    lon_sub_grid, lat_sub_grid = np.meshgrid(lon_sub, lat_sub)

    if lat_sub.shape[0] == 1 and lon_sub.shape[0] == 1:
        mask = np.array([[True]])
    else:
        t = transform.from_bounds(
            west=lon_sub_grid.min(), east=lon_sub_grid.max(),
            south=lat_sub_grid.min(), north=lat_sub_grid.max(),
            width=lon_sub_grid.shape[1], height=lat_sub_grid.shape[0]
        )
        mask = np.flipud(features.rasterize(
            [(geom, 1)],
            out_shape=lon_sub_grid.shape,
            transform=t,
            fill=0,
            all_touched=True,
            dtype="uint8"
        ).astype(bool))

    lat_dim = lat_sub.shape[0]
    lon_dim = lon_sub.shape[0]

    with netCDF4.Dataset(output_temp, "w", format="NETCDF4") as out_nc:
        out_nc.createDimension("time", None)
        out_nc.createDimension("lat", lat_dim)
        out_nc.createDimension("lon", lon_dim)

        time_var = out_nc.createVariable("time", "f4", ("time",))
        lat_var = out_nc.createVariable("lat", "f4", ("lat",))
        lon_var = out_nc.createVariable("lon", "f4", ("lon",))
        data_var = out_nc.createVariable(args["variable"], "f4", ("time", "lat", "lon"), fill_value=-9999, zlib=True, complevel=4, shuffle=True)

        lat_var[:] = lat_sub
        lon_var[:] = lon_sub

        index = 0
        start_time = time.time()
        for i, file in enumerate(args["files"]):
            with netCDF4.Dataset(file) as nc:
                t = int(nc.variables["time"][0])
                data_slice = np.array(nc.variables[args["variable"]][0, lat_indices[0]:lat_indices[-1]+1, lon_indices[0]:lon_indices[-1]+1])
                data_slice[data_slice > 10e6] = np.nan
                data_slice[~mask] = np.nan
                if not np.all(np.isnan(data_slice)):
                    data_slice[np.isnan(data_slice)] = -9999
                    data_var[index, :, :] = data_slice
                    time_var[index] = t
                    index = index + 1
            if (i + 1) % 1000 == 0:
                elapsed = time.time() - start_time
                avg_time = (elapsed / (i + 1) * 1000) / args["threads"]
                print(f"Average time per file: {avg_time:.2f} ms")

    os.rename(output_temp, output_path)
    elapsed = time.time() - start
    print(f"Completed {args['lake']['id']} in {elapsed:.2f} seconds.")


folder = "/media/runnalja/My Passport/Documents/Phenology/data/CCI_lakes_nc/v2.1.0/**/*.nc"
shapefile = "/media/runnalja/My Passport/Documents/Phenology/data/CCI_lakes_ShpFls/lakes_cci_v2.0.2_data_availability_shp/lakescci_v2.0.2_data-availability.shp"
out_folder = "data"
variable = "chla_mean"
n = 10

files = glob.glob(folder, recursive=True)
files.sort()
gdf = gpd.read_file(shapefile)
lakes = [{"lake": row, "variable": variable, "files": files, "out_folder": out_folder, "threads": n} for i, (_, row) in enumerate(gdf.iterrows())]
with ProcessPoolExecutor(max_workers=n) as executor:
    executor.map(process_lake, lakes)

