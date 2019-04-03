import numpy as np
import multiprocessing
import os
from math import ceil
import subprocess

stations = ['03020', '23050', '23047', '23023']
latitudes = [31.8559, 35.0844, 35.2220, 31.9973]
latitudes = list(map(np.deg2rad, latitudes))

station_latitude = {}
for i in range(len(stations)):
	station_latitude[stations[i]] = latitudes[i]

path_to_radiosonde_data = '/Users/willfried/Desktop/radiosonde_data/core_4/'
path_to_csv_files = '/Users/willfried/Desktop/test/Joint-Probability-Distribution-of-Atmospheric-Gravity-Wave-Parameters/'




