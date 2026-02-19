import numpy as np
from scipy.io import loadmat
from datetime import datetime, timezone
import netCDF4
from matplotlib import pyplot as plt
from csaps import csaps

p_file = "/home/runnalja/git/cci-cyano-production/data/v3.1/phenology/chla/327.nc"
e_file = "/home/runnalja/git/cci-cyano-production/data/v3.1/extract/chla/327.nc"
variable = "chla"
qa_name = "lwlr_quality_flags"
x, y = 20, 20

def unix_to_datenum(arr):
    return np.array([datetime.utcfromtimestamp(int(ts)).toordinal() + 366 for ts in np.array(arr)])

def datenum_to_datetime(arr):
    return np.array([datetime.fromordinal(int(dn) - 366).replace(tzinfo=timezone.utc) for dn in np.array(arr)])

def unix_to_datetime(arr):
    return np.array([datetime.fromtimestamp(int(ts), tz=timezone.utc) for ts in np.array(arr)])

def remove_nan(arr):
    arr = np.array(arr)
    return arr[~np.isnan(arr)]

with netCDF4.Dataset(e_file) as nc:
    valid_coords = np.argwhere(nc.variables["summary"][:, :] > 2)
    lat = np.array(nc.variables["lat"])
    lon = np.array(nc.variables["lon"])
    t = unix_to_datenum(nc.variables["time"])
    values = np.array(nc.variables[variable][:, x, y])
    mask = values != -9999
    qa = np.array(nc.variables[qa_name][:, x, y])
    mask = (qa == 0) & mask
    values = values[mask]
    time = t[mask]

with netCDF4.Dataset(p_file) as nc:
    smoothing = float(nc.variables["smoothing_parameter"][x, y])
    pks_x = unix_to_datetime(remove_nan(nc.variables["pks_x"][x, y, :]))
    pks_y = remove_nan(nc.variables["pks_y"][x, y, :])
    trgs_x = unix_to_datetime(remove_nan(nc.variables["trgs_x"][x, y, :]))
    trgs_y = remove_nan(nc.variables["trgs_y"][x, y, :])

smooth_x_axis = np.arange(t.min(), t.max() + 1, 1)
smooth_y_data = csaps(time, values, smooth_x_axis, smooth=smoothing)

plt.scatter(datenum_to_datetime(time), values, color='grey', alpha=0.3, s=10, label='Data')
plt.plot(datenum_to_datetime(smooth_x_axis), smooth_y_data, color='blue', linewidth=1, label='MATLAB spline')
plt.scatter(pks_x, pks_y, color='red', s=60, marker='^', zorder=4, label='Peaks')
plt.scatter(trgs_x, trgs_y, color='green', s=60, marker='v', zorder=4, label='Troughs')

plt.show()
