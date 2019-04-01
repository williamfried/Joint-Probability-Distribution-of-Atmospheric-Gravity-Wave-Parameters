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
quantities_full_names = ['vertical wavelength','frequency','zonal wavelength','meridional wavelength','zonal wind amplitude','meridional wind amplitude','temperature amplitude']
quantities_units = ['(km)', '($10^{-4}$$s^{-1}$)', '(km)', '(km)', '(m/s)', '(m/s)', '(K)']

quantities_info = {}
for i in range(len(quantities)):
	quantities_info[quantities[i]] = {'name' : quantities_full_names[i], 'units' : quantities_units[i]}

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

parameter_dict = {}
for quantity in quantities:
	parameter_dict[quantity] = []

# graph of sounding to sounding difference in parameter values based on number of missing soundings
def reading_to_reading(num):

	missing_color = {0 : '#99000d',
				 	 1 : '#cb181d',
				 	 2 : '#ef3b2c',
				 	 3 : '#fb6a4a',
				 	 4 : '#fc9272',
				 	 5 : '#fcbba1'}

	for station in stations:
		for quantity in quantities:
			prev = None
			cnt = 0
			for reading in data_dict[station]:
				if cnt < num:
					if prev == None:
						prev = data_dict[station][reading][quantity]
					else:
						reading_data = data_dict[station][reading]
						if reading_data['missing'] in missing_color:
							plt.plot([cnt, cnt + 1], [prev, reading_data[quantity]], c=missing_color[reading_data['missing']])
						else:
							plt.plot([cnt, cnt + 1], [prev, reading_data[quantity]], c='#fee5d9')
	
						prev = reading_data[quantity]
						cnt += 1

			plt.xlabel('sounding number')
			plt.ylabel(full_name_dict[quantity] + ' ' + units_dict[quantity])
			plt.title(full_name_dict[quantity] + ' by number of missing soundings')	
			plt.show()

# create bar graphs that illusrate average difference between consecutive soundings as a function of the number of hours between consecutive soundings
def missing_diff():

	missing_dict = {}
	for i in range(6):
		missing_dict[i] = deepcopy(parameter_dict)
	missing_dict['greater'] = deepcopy(parameter_dict)
	categories = ['12', '24', '36', '48', '60', '72', '> 72']
	locations = [1.2*i for i in range(7)]

	prev_dict = {}
	for quantity in quantities:
		prev_dict[quantity] = None

	for station in stations:
		for reading in data_dict[station]:
			for quantity in quantities:
				if prev_dict[quantity] == None:
					prev_dict[quantity] = data_dict[station][reading][quantity]
				else:
					reading_data = data_dict[station][reading]
					if reading_data['missing'] in missing_dict:
						missing_dict[reading_data['missing']][quantity].append(abs(reading_data[quantity] - prev_dict[quantity]))
					else:
						if reading_data['missing'] > 5: 
							missing_dict['greater'][quantity].append(abs(reading_data[quantity] - prev_dict[quantity]))
						elif reading_data['missing'] < 0:
							missing_dict[0][quantity].append(abs(reading_data[quantity] - prev_dict[quantity]))

					prev_dict[quantity] = reading_data[quantity]

		for quantity in quantities:
			prev_dict[quantity] = None

	ll = []
	for quantity in quantities:
		avg_diff_list = []
		for missing_val in missing_dict:
			avg_diff_list.append(np.mean(np.array(missing_dict[missing_val][quantity])))

		overall_mean = np.mean(np.array(avg_diff_list))
		dev = ((avg_diff_list[0] - overall_mean) / overall_mean) * 100

		plt.bar(locations, avg_diff_list, align='center')
		matplotlib.rcParams.update({'font.size': 16})
		plt.xticks(locations, categories)
		plt.xlabel('hours')
		plt.ylabel('average difference between consecutive soundings ' + quantities_info[quantity]['units'])
		plt.title(quantities_info[quantity]['name'] + ': ' + str(round(dev,1)) + '%')
		plt.show()

# compare marginal distributions of parameters inferred at noon and midnight 
def day_vs_night(type):
	matplotlib.rcParams.update({'font.size': 14})
	times = ['00', '12']
	time_dict = {}
	for time in times:
		time_dict[time] = deepcopy(parameter_dict)

	for quantity in quantities:
		for time in times:
			for station in stations:
				for year in years:
					for month in months:
						soundings = data_dict[station + '1'][year][month][time]
						for sounding in soundings:
							time_dict[time][quantity].append(soundings[sounding][quantity])

		if type == 'marginal':
			midnight_density = gaussian_kde(time_dict['00'][quantity])
			noon_density = gaussian_kde(time_dict['12'][quantity])
			xs = np.linspace(0,max(np.quantile(time_dict['00'][quantity], 0.999), np.quantile(time_dict['12'][quantity], 0.999)), num=200)
			plt.plot(xs,midnight_density(xs), label='midnight')
			plt.plot(xs,noon_density(xs), label='noon')
			plt.xlabel(quantities_info[quantity]['name'] + ' ' + quantities_info[quantity]['units'])
			plt.ylabel('density')
			plt.title('Differences in marginal distribution by time of day')
			plt.legend(loc='best')
			plt.show()

	if type == 'correlation':
		for i in range(len(quantities)-1):
			for j in range((i+1),len(quantities)):
				plt.scatter(time_dict['00'][quantities[i]], time_dict['00'][quantities[j]], color='b', s=5)
				plt.scatter(time_dict['12'][quantities[i]], time_dict['12'][quantities[j]], color='r', s=5)
				plt.xlabel(quantities_info[quantities[i]]['name'] + ' ' + quantities_info[quantities[i]]['units'])
				plt.ylabel(quantities_info[quantities[j]]['name'] + ' ' + quantities_info[quantities[j]]['units'])
				plt.title('Spearman\'s rank correlation coefficient: ' + str(round(spearmanr(time_dict['00'][quantities[i]], time_dict['00'][quantities[j]])[0],2)) + ' ' + str(round(spearmanr(time_dict['12'][quantities[i]], time_dict['12'][quantities[j]])[0],2)))
				plt.show()

missing_diff()
day_vs_night('marginal')

