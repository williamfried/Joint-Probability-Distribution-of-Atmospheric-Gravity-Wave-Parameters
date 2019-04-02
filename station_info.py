import numpy as np
import multiprocessing
import os
from math import ceil
import subprocess

stations = ['13985', '13996']
latitudes = [37.7528, 39.0473]
latitudes = list(map(np.deg2rad, latitudes))

station_latitude = {}
for i in range(len(stations)):
	station_latitude[stations[i]] = latitudes[i]

path_to_radiosonde_data = '/Users/willfried/Desktop/trial/'
path_to_csv_files = '/Users/willfried/Desktop/data/'



