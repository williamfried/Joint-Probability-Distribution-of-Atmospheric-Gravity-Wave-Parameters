import numpy as np 
import matplotlib.pyplot as plt
import matplotlib
import os
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.stats import norm
from scipy.stats import gaussian_kde
from copy import deepcopy
import csv
from station_info import *

years = [str(x + 1998) for x in range(11)]
months = ['0' + str(x + 1) for x in range(9)] + ['10', '11', '12']
days = [str((x+1)) for x in range(31)]
times = ['0' + str(x) for x in range(10)] + [str(x) for x in range(10,24)]
quantities = ['vertical_wavelength','freq','U_wavelength','V_wavelength','U_amp','V_amp','T_amp']

winter_months = ['01', '02', '03', '04', '11', '12']
summer_months = ['06', '07', '08', '09']
intermediate_months = ['05', '10']

# store data from all stations in one dictionary
data_dict = {}
for station in stations:

	os.chdir(path_to_radiosonde_data + station)

	s = open('zodograph.txt', 'r').read()
	hodograph_dict = eval(s)

	s = open('years_dict.txt', 'r').read()
	years_dict = eval(s)

	data_dict[station] = hodograph_dict
	data_dict[station + '1'] = years_dict

os.chdir(path_to_csv_files)

# store data by year
def yearly_readings():
	for quantity in quantities:
		cnt = 0
		for year in years:
			l = []
			for station in stations:
				for month in months:
					for time in times:
						readings = data_dict[station + '1'][year][month][time]
						for reading in readings:
							if quantity in set(['U_wavelength','V_wavelength']):
								l.append(readings[reading][quantity] / 1000)
							else:
								l.append(readings[reading][quantity])

			mode = 'w' if cnt == 0 else 'a'
			with open(quantity + '_year.csv', mode) as f:
				writer = csv.writer(f)
				writer.writerow(l)
			cnt += 1

# store data by month 
def monthly_readings():
	for quantity in quantities:
		cnt = 0
		for month in months:
			l = []
			for station in stations:
				for year in years:
					for time in times:
						readings = data_dict[station + '1'][year][month][time]
						for reading in readings:
							if quantity in set(['U_wavelength','V_wavelength']):
								l.append(readings[reading][quantity] / 1000)
							else:
								l.append(readings[reading][quantity])

			mode = 'w' if cnt == 0 else 'a'
			with open(quantity + '_month.csv', mode) as f:
				writer = csv.writer(f)
				writer.writerow(l)
			cnt += 1

# store data by station 
def station_readings():
	for quantity in quantities:
		cnt = 0
		for station in stations:
			l = []
			for year in years:
				for month in months:
					for time in times:
						readings = data_dict[station + '1'][year][month][time]
						for reading in readings:
							if quantity in set(['U_wavelength','V_wavelength']):
								l.append(readings[reading][quantity] / 1000)
							else:
								l.append(readings[reading][quantity])

			mode = 'w' if cnt == 0 else 'a'
			with open(quantity + '_station.csv', mode) as f:
				writer = csv.writer(f)
				writer.writerow(l)
			cnt += 1

station_readings()

# group soundings by season (winter months (November-April), summer months (June-September) or all months of year)
def season_readings(season):
	season_dict = {'winter' : winter_months,
				   'summer' : summer_months,
				   'all' : months}
	cnt = 0
	for quantity in quantities:
		l = []
		for month in season_dict[season]:
			for station in stations:
				for year in years:
					for time in times:
						readings = data_dict[station + '1'][year][month][time]
						for reading in readings:
							if quantity in set(['U_wavelength','V_wavelength']):
								l.append(readings[reading][quantity] / 1000)
							else:
								l.append(readings[reading][quantity])

		mode = 'w' if cnt == 0 else 'a'
		with open(season + '.csv', mode) as f:
			writer = csv.writer(f)
			writer.writerow(l)
		cnt += 1

yearly_readings()
monthly_readings()
season_readings('winter')
season_readings('summer')
season_readings('all')
