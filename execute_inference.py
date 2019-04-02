import numpy as np
import sys
import os
from inference import *
from station_info import *
import json
from multiprocessing import Pool

# apply hodograph method to all soundings from upper air station
def execute_inference(station):
	path = path_to_radiosonde_data + station
	os.chdir(path)
	files = os.listdir(path)
	previous_file = None

	dat_files = sorted([x for x in files if x.endswith('.dat') and x[0].isnumeric()])

	for filename in dat_files:
		print(filename)
		previous_file = hodograph(filename, 18, 25, 50, previous_file)
	with open('zodograph.txt','w') as f:
		f.write(json.dumps(hodograph_dict))
	with open('years_dict.txt','w') as f:
		f.write(json.dumps(years_dict))

	return(0)

p = Pool()
p.map(execute_inference, stations)



