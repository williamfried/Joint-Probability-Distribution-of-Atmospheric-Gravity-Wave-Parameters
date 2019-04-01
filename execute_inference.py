import numpy as np
import sys
from inference import *
from station_info import *
import json

station = sys.argv[1]

stations_4 = set(['03020','23050','23047','23023'])
stations_7 = set(['03160', '04105', '53103', '23160', '24127', '23062', '23066'])

def station_directory(station):
	if station in stations_4:
		return('core_4')
	elif station in stations_7:
		return('extra_7')
	else:
		return(-1)

# apply hodograph method to all soundings from upper air station
def execute():
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

execute()



