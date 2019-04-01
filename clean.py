import numpy as np 
import os
from station_info import *
from preprocessing import *
import sys

station_code = sys.argv[1]

# filter out soundings that either don't ascend to 25km or have too many missing values
def clean(filename):

	data = extract(filename, True)

	if data == 'short':
		os.rename(filename, 'short' + filename)
		return(0)

	data = data[data[:,13] >= 18*1000]
	data = data[data[:,13] <= 25.1*1000]

	# extract altitude as well as velocity and wind profiles
	alt = data[:,13]
	alt /= 1000

	U_velocity = data[:,4]
	V_velocity = data[:,5]
	temp = data[:,1]

	alt_filtered = alt[(U_velocity < 900) & (V_velocity < 900) & (temp < 900)]

	if len(alt_filtered) == 0 or max(alt_filtered) < 25:
		os.rename(filename, 'below25' + filename)
		return(0)

	x = np.count_nonzero(U_velocity > 900) + np.count_nonzero(V_velocity > 900) + np.count_nonzero(temp > 900)
	if x > 500: 
		os.rename(filename, 'five_hundred' + filename)
	elif x > 200: 
		os.rename(filename, 'two_hundred' + filename)
	elif x > 100: 
		os.rename(filename, 'one_hundred' + filename)
	elif x > 50: 
		os.rename(filename, 'fifty' + filename)
	elif x > 10: 
		os.rename(filename, 'ten' + filename)
	
	return(0)

# undo effect of 'clean' function above
def undo(filename):
	if filename.startswith('two_hundred'):
		os.rename(filename, filename[len('two_hundred')::])
	elif filename.startswith('ten'):
		os.rename(filename, filename[len('ten')::])
	elif filename.startswith('one_hundred'):
		os.rename(filename, filename[len('one_hundred')::])
	elif filename.startswith('five_hundred'):
		os.rename(filename, filename[len('five_hundred')::])
	elif filename.startswith('fifty'):
		os.rename(filename, filename[len('fifty')::])
	elif filename.startswith('below25'):
		os.rename(filename, filename[len('below25')::])	

	return(0)

def clean_files():
	path = path_to_radiosonde_data + station_code
	os.chdir(path)
	files = os.listdir(path)
	dat_files = sorted([x for x in files if x.endswith('.dat') and x[0].isnumeric()])

	for filename in dat_files:
		print(filename)
		clean(filename)

	return(0)

def undo_files():
	path = path_to_radiosonde_data + station_code
	os.chdir(path)
	files = os.listdir(path)
	dat_files = sorted([x for x in files if x.endswith('.dat') and x[0].isalpha()])
	for filename in dat_files:
		print(filename)
		undo(filename)
	return(0)

def execute(action):
	if action == 'clean':
		clean_files()
		return(0)
	elif action == 'undo':
		undo_files()
		return(0)
	else:
		return('invalid selection')

execute(sys.argv[2])


