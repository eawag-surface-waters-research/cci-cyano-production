import os
import netCDF4
import numpy as np
import pandas as pd
from scipy.io import loadmat
from datetime import datetime
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

file = "../data/v3.1/extract/chla/327.nc"
base = "/home/runnalja/git/bgb-jelle/data/CCISTEP1_timeSeriesAndNullInfo/Lake_297_327_GLWD00000327/"

index = 586

with netCDF4.Dataset(file) as nc:
    coords = np.argwhere(nc.variables["summary"][:, :] != 0)
    x = nc.variables["chla"][:, coords[index][0], coords[index][1]]
    t = [datetime.utcfromtimestamp(int(ts)) for ts in np.array(nc.variables["time"])]

info = loadmat(base + "CCI_NullInfo_chla_mean_731335_738886.mat")
t_mat = pd.to_datetime(np.array(info["timeSeries_dateNR"][0]) - 719529, unit='D', origin='unix')

data = loadmat(base + "CCI_Data_chla_mean_731335_738886_pixNRs_1_982.mat")
x_mat = np.array(data["timeSeries_data"][index])


result = loadmat("/home/runnalja/git/bgb-jelle/data/CCISTEP3_phenoMETs/Lake_297_327_GLWD00000327/pixFileDataList_chla_mean_731335_738886/cubicSpline_1_14_0_0_31_0_365_5_pixNRs_1_982.mat")

print(result.keys())
exit()

#plt.scatter(t, x, color="black")
plt.plot(t_mat, x_mat, color="red", marker="x")

print(len(t_mat), len(t))

plt.show()

x1 = x[~np.isnan(x)]
x2 = x_mat[~np.isnan(x_mat)]

print(len(x1), len(x2))
