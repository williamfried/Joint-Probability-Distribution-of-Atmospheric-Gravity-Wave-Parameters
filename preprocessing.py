import numpy as np 
import os
from station_info import *

# extract data from raw DAT file
def extract(filename, clean):

	# rows before first combination of two rows into one
	f = open(filename,'r')
	try:
		data0 = np.genfromtxt((" ".join(ln.split()[3:17]) for ln in f if len(ln.split()) == 23))
	except:
		return('short')

	data0 = data0[data0[:,4] < 500] 

	# rows after first combination of two rows into one
	f = open(filename,'r')
	data1 = np.genfromtxt((" ".join(ln.split()[2:16]) for ln in f if len(ln.split()) == 22))

	# rows after second combination of two rows into one
	f = open(filename,'r')
	data2 = np.genfromtxt((" ".join(ln.split()[1:15]) for ln in f if len(ln.split()) == 21)) 

	# combine into one array
	if len(data1.shape) == 1 or len(data1) == 0:
		return('short')

	if len(data2) != 0 and len(data2.shape) != 1:
		data = np.concatenate((data0, data1, data2), axis=0)	
	else:
		data = np.concatenate((data0, data1), axis=0)

	data = data[~np.isnan(data).any(axis=1)]

	if not clean:
		data = data[(data[:,4] < 900)]  
		data = data[(data[:,5] < 900)]
		data = data[(data[:,1] < 900)]

	return(data)



	
