import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import warnings
import os
import sys
from inference_functions import *
from station_info import *
from copy import deepcopy
warnings.filterwarnings('ignore', category=UserWarning) 

# quantities to interpolate 
quantities = ['P','T','U','V','alt','dens']

def energy(filename, lower, upper, interval):

	# specify altitudes that are equally spaced at 50m intervals
	alt_interval = np.linspace(lower, upper, ((1000 / interval) * (upper - lower)) + 1)

	# extract data from file
	data = extract(filename, False)

	# limit analysis to lower stratosphere 
	data = data[data[:,13] >= lower * 1000 - interval * 2]
	data = data[data[:,13] <= upper * 1000 + interval * 2]

	# extract relevant data
	alt = data[:,13]
	alt /= 1000
	temp = data[:,1]
	U_velocity = data[:,4]
	V_velocity = data[:,5]
	pressure = data[:,0]
	RH = data[:,3]

	# higher order polynomial fit
	raw_data = np.concatenate((np.array([U_velocity]).T, np.array([V_velocity]).T, np.array([temp]).T), axis=1)
	profile_fit = polyfit(alt, raw_data, 2)

	# residuals of 3 quantities
	residual = raw_data - profile_fit 
	U_perturb = residual[:,0]
	V_perturb = residual[:,1]
	temp_perturb = residual[:,2]
	
	# interpolate quantities
	U_interpolated = interpolate_array(alt, U_perturb, alt_interval)
	V_interpolated = interpolate_array(alt, V_perturb, alt_interval)
	temp_interpolated = interpolate_array(alt, temp_perturb, alt_interval)
	background_temp = interpolate_array(alt, profile_fit[:,2], alt_interval)
	pressure_interpolated = interpolate_array(alt, pressure, alt_interval)
	RH_interpolated = interpolate_array(alt, RH, alt_interval)
	temp_interpolated = interpolate_array(alt, temp, alt_interval)

	# calculate vertical velocity fluctuation using hydrostatic altitude
	hydro_alt = hydrostatic_alt_diff_array(temp, pressure, RH)
	hydro_velo = hydro_alt / 12
	hydro_velo_interpolated = interpolate_array(alt[1:-1], hydro_velo, alt_interval)

	# remove 5 km running average from vertical velocity flunctation profile
	w_perturb_temp = np.zeros(len(hydro_velo_interpolated))
	for i in range(len(w_perturb_temp)):
		lower_end = max(0, i-50)
		upper_end = min(len(w_perturb_temp),i+50)
		w_perturb_temp[i] = hydro_velo_interpolated [i] - np.mean(hydro_velo_interpolated[lower_end:upper_end])
	w_perturb = w_perturb_temp

	# plot vertical velocity flunctation as a function of altitude
	def plot_vertical_velocity():
		plt.rcParams.update({'font.size': 12})
		plt.plot(w_perturb, alt_interval, 'g')
		plt.xlabel('vertical wind perturbation (m/s)')
		plt.ylabel('altitude (km)')
		plt.title('vertical wind perturbation at 50m intervals in lower stratosphere')
		plt.show()
	
	# estimate Brunt-Vaisala frequency in lower stratosphere 
	N = brunt_vaisala(pressure_interpolated, temp_interpolated, alt_interval)

	# calculate energy density associated with atmospheric gravity wave
	E = energy_density(U_interpolated, V_interpolated, W_interpolated, background_temp, temp_interpolated, N)

	def plot_energy_density():
		plt.rcParams.update({'font.size': 18})
		plt.plot(E, alt_interval)
		plt.xlabel('specific energy (J/kg)')
		plt.ylabel('altitude (km)')
		plt.title('Specific energy at 50m intervals in lower stratosphere')
		plt.show()
		print(np.mean(E))
		print(np.std(E))

	return(0)

# create dictionary to store gravity wave parameters
time_dict, month_dict, years_dict = {}, {}, {} 
for time in times:
	time_dict[time] = {}
for month in months:
	month_dict[month] = deepcopy(time_dict)
for year in years:
	years_dict[year] = deepcopy(month_dict)

hodograph_dict = {}

# infer 7 gravity waves parameters using hodograph method
def hodograph(filename, lower, upper, interval, prev_file):

	# extract data from file
	data = extract(filename, False)

	# specify altitudes that are equally spaced at 50m intervals
	alt_interval = np.linspace(lower, upper, ((1000 / interval) * (upper - lower)) + 1)

	# limit analysis to lower stratosphere 
	data = data[data[:,13] >= lower * 1000 - interval * 2]
	data = data[data[:,13] <= upper * 1000 + interval * 2]

	# extract relevant data
	alt = data[:,13]
	alt /= 1000
	temp = data[:,1]
	U_velocity = data[:,4]
	V_velocity = data[:,5]

	# higher order polynomial fit
	raw_data = np.concatenate((np.array([U_velocity]).T, np.array([V_velocity]).T, np.array([temp]).T), axis=1)
	profile_fit = polyfit(alt, raw_data, 2)

	# residuals of 3 quantities
	residual = raw_data - profile_fit 
	U_perturb = residual[:,0]
	V_perturb = residual[:,1]
	temp_perturb = residual[:,2]

	# linearly interpolate quantities at 50m intervals between 18-25 km
	U_interpolated = interpolate_array(alt, U_perturb, alt_interval)
	V_interpolated = interpolate_array(alt, V_perturb, alt_interval)
	temp_interpolated = interpolate_array(alt, temp_perturb, alt_interval)

	# plot different polynomial fit
	def graph_poly():
		plt.rcParams.update({'font.size': 12})
		colors = ['#a50f15', '#cb181d', '#ef3b2c', '#fb6a4a', '#fc9272', '#fcbba1']
		plt.plot(V_velocity, alt, 'k-', label='raw data', zorder=10)
		for i in range(2, 8):
			profile_fit_iter = polyfit(alt, raw_data, i)
			plt.plot(profile_fit_iter[:,1], alt, colors[i-2], label=str(i), zorder=10-i)
		plt.xlabel('meridional wind velocity (m/s)')
		plt.ylabel('altitude (km)')
		plt.title('Comparison of different degree polynomial fits')
		plt.legend(loc='top right')
		plt.show()

	# fit sine curve (amplitude, frequency, phase) to each quantity 	
	parameters1, wave1 = fit_sin(alt_interval, U_interpolated)
	parameters2, wave2 = fit_sin(alt_interval, V_interpolated)
	parameters3, wave3 = fit_sin(alt_interval, temp_interpolated)

	#plot best fit sine curve
	def graph_sine():
		A,w,p = parameters1
		aa = A*np.sin(alt_interval*w + p)
		plt.plot(aa, alt_interval, 'g')
		plt.scatter(U_interpolated, alt_interval, s=1)
		plt.show()

	# calculate relative standard error of wavelengths of 3 quantities 
	wavelengths = np.array([wave1, wave2, wave3])
	vert_wavelength = np.mean(wavelengths)
	diff = wavelengths - vert_wavelength
	rel_standard_error = (np.sqrt(sum(diff*diff) / 2) / np.sqrt(3)) / vert_wavelength

	# only infer soundings where relative standard error is below 0.2
	if rel_standard_error < 0.2:

		# refit sine curve to 3 quantities using average wavelength calculated above
		def sine_fit (x, phase, amp):
			return(amp * (np.sin(((2*np.pi) / vert_wavelength) * x + phase)))

		param_U, _ = curve_fit(sine_fit, alt_interval, U_interpolated)
		param_V, _ = curve_fit(sine_fit, alt_interval, V_interpolated)
		param_temp, _ = curve_fit(sine_fit, alt_interval, temp_interpolated)

		# inferred values of amplitudes of 3 quantities
		(_, U_amp), (_, V_amp), (_, T_amp) = param_U, param_V, param_temp

		# fit hodograph
		x_data = np.linspace(0, vert_wavelength, 50)

		u = sine_fit(x_data, *param_U)
		v = sine_fit(x_data, *param_V)

		l1, l2, angle = calc_ellipse(u, v, 1, float('-Inf'))

		if np.isreal(l2):
			
			# calculate ratio between major and minor axis
			major = l1 if l1 > l2 else l2
			minor = l2 if l1 > l2 else l1
			axis_ratio = np.real(major / minor)
			
			# only infer soundings where ratio between major and minor axis is below 10 
			if 0 < axis_ratio < 10:

				# compute Coriolis frequency value based on geographical location of upper air station
				cor = 2 * coriolis * np.sin(station_latitude[filename[0:5]])

				# infer frequency of atmospheric gravity wave
				freq = axis_ratio * cor

				# estimate Brunt-Vaisala frequency in lower stratosphere 
				abs_temp_interpolated = interpolate_array(alt, temp, alt_interval)
				pressure_interpolated = interpolate_array(alt, data[:,0], alt_interval)
				N = brunt_vaisala(pressure_interpolated, abs_temp_interpolated, alt_interval)

				# infer zonal and meridional wavelengths using dispersion equation for inertial gravity waves
				vert_wavenumber = 1 / vert_wavelength
				hor_wavelength = 1 / (np.sqrt((((vert_wavenumber) ** 2) * ((freq ** 2) - (cor ** 2))) / (N ** 2)))
				U_wavelength = abs(hor_wavelength * np.cos(angle))
				V_wavelength = abs(hor_wavelength * np.sin(angle))

				#print(between_soundings(prev_file, filename))

				# combine all 7 inferred atmospheric gravity wave parameters 
				# determine number of missing soundings between given sounding and previous sounding
				inferred_data = {'freq' : freq * (10**4), 'vertical_wavelength' : vert_wavelength, 'horizontal_wavelength' : hor_wavelength, 'U_wavelength' : U_wavelength, 'V_wavelength' : V_wavelength, 'missing' : between_soundings(prev_file, filename), 'U_amp' : abs(U_amp), 'V_amp' : abs(V_amp), 'T_amp' : abs(T_amp)} 

				# organize inferred data by date of sounding 
				years_dict[filename[6:10]][filename[10:12]][filename[14:16]][filename] = inferred_data

				# organize inferred data by name of sounding 
				hodograph_dict[filename] = inferred_data

				return(filename)
	return(prev_file)
