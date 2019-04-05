import numpy as np
import multiprocessing
import os
from math import ceil
import subprocess

stations = []
latitudes = []
latitudes = list(map(np.deg2rad, latitudes))

station_latitude = {}
for i in range(len(stations)):
	station_latitude[stations[i]] = latitudes[i]

path_to_radiosonde_data = ''
path_to_csv_files = ''




