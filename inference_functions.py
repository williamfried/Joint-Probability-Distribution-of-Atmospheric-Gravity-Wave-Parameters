from numpy.linalg import eig, inv
from scipy.optimize import curve_fit
import scipy 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import time
import warnings
import os
import sys
from preprocessing import *
warnings.filterwarnings('ignore', category=UserWarning) 
warnings.filterwarnings('ignore', category=FutureWarning) 

# calculate best fit ellipse from hodograph of zonal and meridional winds (http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html)
# returns length of major and minor axis and angle between them 
def calc_ellipse (u, v, num, prev_offset):

	# compute angle that major axis makes to horizontal  
	def slope(xx, yy):
		index = None
		maximum = float('-inf')
		for i in range(len(xx)):
			dist = (xx[i])**2 + (yy[i])**2
			if dist > maximum:
				maximum = dist
				index = i
		
		return(np.arctan(np.real(yy[index] / xx[index])))

	# deviation of hodograph from fitted ellipse 
	def offset(u, v, xx, yy):

		def closest_point(x1, y1, xx, yy):
			minimum = float('inf')
			for i in range(len(xx)):
				dist = (xx[i] - x1)**2 + (yy[i] - y1)**2
				if dist < minimum:
					minimum = dist
			return(np.sqrt(minimum))

		cnt = 0
		for i in range(len(u)):
			cnt += closest_point(u[i],v[i], xx, yy)
		return(cnt)

	def fitEllipse(x,y):
	    x = x[:,np.newaxis]
	    y = y[:,np.newaxis]
	    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
	    S = np.dot(D.T,D)
	    C = np.zeros([6,6])
	    C[0,2] = C[2,0] = 2; C[1,1] = -1
	    try:
	    	E, V =  eig(np.dot(inv(np.round(S,10)), C))
	    except:
	    	return('singular matrix')
	    n = np.argmax(np.abs(E))
	    a = V[:,n]
	    return a

	def ellipse_center(a):
	    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
	    num = b*b-a*c
	    x0=(c*d-b*f)/num
	    y0=(a*f-b*d)/num
	    return np.array([x0,y0])

	def ellipse_angle_of_rotation2(a):
	    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
	    if b == 0:
	        if a > c:
	            return 0
	        else:
	            return np.pi/2
	    else:
	        if a > c:
	            return np.arctan(2*b/(a-c))/2
	        else:
	            return np.pi/2 + np.arctan(2*b/(a-c))/2

	def ellipse_axis_length(a):
	    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
	    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
	    down1=(b*b-a*c)*( num*(c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
	    down2=(b*b-a*c)*( num*(a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
	    res1=np.sqrt(up/down1)
	    res2=np.sqrt(up/down2)
	    return np.array([res1, res2])

	arc = 2
	R = np.arange(0,arc*np.pi, 0.01)

	a = fitEllipse(u,v)
	if (a == 'singular matrix'):
		return(100, 1, 0)
	axes = ellipse_axis_length(a)
	center = ellipse_center(a)
	phi = ellipse_angle_of_rotation2(a)

	a, b = axes

	xx = center[0] + a*np.cos(R)*np.cos(phi) - b*np.sin(R)*np.sin(phi)
	yy = center[1] + a*np.cos(R)*np.sin(phi) + b*np.sin(R)*np.cos(phi)

	offset_val = np.real(offset(u, v, xx, yy))
	if prev_offset == float('-Inf'):
		return(calc_ellipse(u, v, -1 * num, offset_val))
	elif offset_val > prev_offset:
		return(calc_ellipse(u, v, -1 * num, offset_val))

	angle = slope(xx, yy)

	def plot_ellipse():
		matplotlib.rcParams.update({'font.size': 12})
		plt.plot(u,v,'.')
		plt.plot(xx,yy)
		plt.axis('equal')
		plt.xlabel('zonal wind (m/s)')
		plt.ylabel('meridional wind (m/s)')
		plt.title('Typical hodograph diagram')
		val = 2.42
		plt.plot([-val*np.cos(angle), val*np.cos(angle),],[-val*np.sin(angle),val*np.sin(angle)],'g--')
		plt.show()

	return(a,b,angle)

# column indices corresponding to quantities in DAT files 
index = {'P' : 0,
		 'T' : 1,
		 'U' : 4,
		 'V' : 5,
		 'alt' : 13}

# rotation rate of earth (rad/s)
coriolis = 7.2921e-5

# calendar info 
years = [str(x + 1998) for x in range(11)]
months = ['0' + str(x + 1) for x in range(9)] + ['10', '11', '12']
times = ['0' + str(x) for x in range(10)] + [str(x) for x in range(10,24)]

# linearly interpolate value between two points
def interpolate(x1, x2, y1, y2, x):
	if x2-x1 != 0:
		return(y1 + ((x - x1) / (x2 - x1)) * (y2 - y1))
	return((y1 + y2) / 2)

# find closest value in array to given value and return its index
def find_nearest(array, value):
    return((np.abs(array - value)).argmin())

# compute density given temperature and pressure using ideal gas equation of state
def density(T, P):
	return((0.1 * P) / (0.287 * (T + 273.15)))

# calculate energy density associated with gravity wave at particular altitude
def energy_density(u, v, w, T_avg, T_per, N):
	KE = (u*u + v*v + w*w) / 2
	PE = (((9.81**2) * (T_per ** 2)) / ((N**2) * (T_avg**2))) / 2
	return(KE + PE)

# calculate saturation pressure of water vapor as a function of temperature
def virtual_temp(T, P, RH):
	T += 273.15
	p_sat = np.exp(54.842763 - (6763.22 / T) - 4.21*np.log(T) + 0.000367*T + np.tanh(0.0415*(T - 218.8)) * (53.878 - (1331.22 / T) - 9.44523 * np.log(T) + 0.014025*T))
	pp_wv = p_sat * RH
	mr = pp_wv / (100*P)
	return(T*(1+(0.61*mr)))

# calculate difference in hydrostatic altitude between two points
def hydrostatic_alt_diff(Tl, Pl, RHl, Th, Ph, RHh):
	VTl = virtual_temp(Tl, Pl, RHl)
	VTh = virtual_temp(Th, Ph, RHh)
	return((287/9.81)*np.mean([VTl, VTh])*np.log(Pl/Ph))

# calculate difference in hydrostatic altitude between many pairs of points
def hydrostatic_alt_diff_array(T, P, RH):
	values = np.zeros(len(T)-2)
	for i in range(1,len(T)-1):
		values[i-1] = hydrostatic_alt_diff(T[i-1], P[i-1], RH[i-1], T[i+1], P[i+1], RH[i+1])
	return(values)

# linearly interpolate data at equally spaced altitudes
# alt_data: raw altitude
# x_data: data to be linearly interpolated
# altitudes: altitudes for interpolation
def interpolate_array(alt_data, x_data, altitudes):
	
	store_array = np.zeros(len(altitudes))
	
	for i in range(len(altitudes)):
		equal = None
		idx = find_nearest(alt_data, altitudes[i])
		if alt_data[idx] > altitudes[i]:
			above = idx
			below = idx - 1
		elif alt_data[idx] < altitudes[i]:
			above = idx + 1
			below = idx 
		else:
			equal = idx

		if equal == None:
			store_array[i] = interpolate(alt_data[below], alt_data[above], x_data[below], x_data[above], altitudes[i])
		else:
			store_array[i] = x_data[equal]

	return(store_array)

# higher order least squares polynomial fit
def polyfit(alt, raw_data, degree):

	fit = np.polyfit(alt, raw_data, degree)

	alt_matrix = np.ones((len(alt), 1))
	for i in range(1, degree+1):
		alt_matrix = np.concatenate((np.array([list(map(lambda x: x**i, alt))]).T, alt_matrix), axis=1)

	profile_fit = np.dot(alt_matrix, fit)

	return(profile_fit)

# least squares sine curve fit (https://stackoverflow.com/questions/16716302/how-do-i-fit-a-sine-curve-to-my-data-with-pylab-and-numpy)
def fit_sin(tt, yy):

	def sinfunc(t, A, w, p):  
		return (A * np.sin(w*t + p))

	ff = np.fft.fftfreq(len(tt), d=(tt[1]-tt[0]))
	Fyy = abs(np.fft.fft(yy))
	guess_freq = abs(ff[np.argmax(Fyy[1:])+1])
	guess_amp = np.std(yy) * 2.**0.5
	guess = np.array([guess_amp, 2.*np.pi*guess_freq, 0.])
	popt, _ = curve_fit(sinfunc, tt, yy, p0=guess, maxfev=20000)
	A, w, p = popt
	wave = (2.*np.pi)/w
	return(np.array([A, w, p]), wave)


# calendar infomation (needed for between_soundings function below)
time_dict = {'month_no_leap' : {1 : 31, 2 : 28, 3 : 31, 4 : 30, 5 : 31, 6 : 30, 7 : 31, 8 : 31, 9 : 30, 10 : 31, 11 : 30, 12 : 31},
			 'month_leap' : {1 : 31, 2 : 29, 3 : 31, 4 : 30, 5 : 31, 6 : 30, 7 : 31, 8 : 31, 9 : 30, 10 : 31, 11 : 30, 12 : 31}}

# determine number of missing soundings between two consecutive soundings 
def between_soundings(sounding1, sounding2):

	if sounding1 == None:
		return(0)
	
	r1, r2 = sounding1[6:-4], sounding2[6:-4]
	r1_year, r2_year, r1_month, r2_month, r1_day, r2_day, r1_hour, r2_hour = r1[0:4], r2[0:4], r1[4:6], r2[4:6], r1[6:8], r2[6:8], r1[8:10], r2[8:10]
	
	if (int(r2_year) - int(r1_year)) == 1:
		return(1)
	cnt = 0
	cnt += (int(r2_hour) - int(r1_hour)) / 12
	cnt += 2 * (int(r2_day) - int(r1_day))
	if (int(r2_month) - int(r1_month)) != 0:
		if int(r2_year) % 4 == 0:
			cnt += 2 * time_dict['month_leap'][int(r1_month)]
		else:
			cnt += 2 * time_dict['month_no_leap'][int(r1_month)]
	return(round((cnt - 1), 2))

# compute Brunt Vaisala frequency for a given sounding by averaging over values across lower stratosphere
def brunt_vaisala(pressure, temp, alt):

	def pot_temp(t, p):
		return(t * ((1000 / p) ** 0.286))

	temp_K = temp + 273.15

	BV_array = np.zeros(len(alt) - 5)

	for i in range(len(BV_array)):
		pt_lower = pot_temp(temp_K[i], pressure[i])
		pt_upper = pot_temp(temp_K[i+5], pressure[i+5])
		pt_avg = (pt_lower + pt_upper) / 2
		grad = (pt_upper - pt_lower) / (alt[i+5] - alt[i])

		if grad > 0:
			BV_array[i] = np.sqrt(9.81e-3 * grad / pt_avg) 
		else:
			BV_array[i] = 0

	return(np.mean(BV_array))

# scatterplot of temperature as a function of altitude
def temp_profile(filename):

	data = extract(filename, False)

	t = data[:,1]
	a = data[:,13]
	a /= 1000

	plt.scatter(t,a, s=1)
	plt.xlabel('temperature (C)')
	plt.ylabel('altitude (km)')
	plt.title('temperature as a function of altitude')
	plt.axhline(y=a[np.argmin(t)], color = 'r', )
	plt.show()

	return(0)

