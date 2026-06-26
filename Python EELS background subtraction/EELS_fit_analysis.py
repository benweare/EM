'''
Python version of EELS_fit_analysis.m script.

Author: E Weare
Location: nmRC
Contact: benjamin.weare1@nottingham.ac.uk

Please cite if you found this script useful.
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import trapezoid

# Define functions.

# Import from .msa file as [x, y].
def _import_data_from_file( datafile, delim ):
    data = np.genfromtxt(datafile, delimiter=delim, comments='#')
    array = np.array( data )
    return array


# Extract spectrum region of interest.
def _extract_roi( inarr, start, end ):
	if inarr.ndim == 2:
		rows = np.where( inarr[:,0] > start )
		temp = inarr[rows]
		rows = np.where( temp[:,0] < end )
		outarr = temp[rows]
	if inarr.ndim == 1:
		rows = np.where( inarr > start )
		temp = inarr[rows]
		rows = np.where( temp < end )
		outarr = temp[rows]
	return outarr


# Normalsie range of array.
def _normalise_data_range( data, dmin=0, dmax=1 ):
    return ((data-np.min(data))/(np.max(data)-np.min(data)))*( dmax - dmin )


# Print the results of the fitting.
def _print_best_fit( popt, p_sigma, sigma ):
	print('\nmu +/- ' + str(sigma) + 'sigma:')
	if len(popt) == 4:
		print( 'a = ' + str(popt[0]) + ' +/- ' + str(p_sigma[0]*sigma) )
		print( 'b = ' + str(popt[1]) + ' +/- ' + str(p_sigma[1]*sigma) )
		print( 'c = ' + str(popt[2]) + ' +/- ' + str(p_sigma[2]*sigma) )
		print( 'd = ' + str(popt[0]) + ' +/- ' + str(p_sigma[3]*sigma) )
	if len(popt) == 3:
		print( 'a = ' + str(popt[0]) + ' +/- ' + str(p_sigma[0]*sigma) )
		print( 'b = ' + str(popt[1]) + ' +/- ' + str(p_sigma[1]*sigma) )
		print( 'c = ' + str(popt[2]) + ' +/- ' + str(p_sigma[2]*sigma) )
	if len(popt) == 2:
		print( 'a = ' + str(popt[0]) + ' +/- ' + str(p_sigma[0]*sigma) )
		print( 'b = ' + str(popt[1]) + ' +/- ' + str(p_sigma[1]*sigma) )
	else:
		print( 'Could not print optimised parameters.' )
	return


# Import .msa data exported from Digital Micrograph 3.
filename = 'file_name.msa'; # This is your file name.

data1 = _import_data_from_file( filename, ',' )

# There may be situations where the data of interest is a smaller section
# of the total spectrum. The code below allows the user to define the start
# and end of that data. A larger background before the start of the edge
# tends to lead to better fitting of the background.

# Extracting the edge. Define the start of the edge.
startedge = 175
endedge = 380
data2 = _extract_roi( data1, startedge, endedge )

# Analysing exponential fitting
# Exclude data from table where xdata is above i eV and change i to the
# appropriate value.
i = 260
fitmask = np.where( data2[:,0] < i )

# Define fit - use one of the fits in single quotations below:
# 'exp1' is a one-term exponential model f(x) = a*exp(b*x)
# 'exp2' is a two-term exponential model f(x) = a*exp(b*x) + c*exp(d*x)
# 'power1' is a one-term power model f(x) = a*x^b
# 'power2' is a two-term power model f(x) = a*x^b + c

# If any other fits are needed, you need to define the fittype as they're
# not included in the MATLAB Curve Fitting Toolbox e.g. exp3 is a
# three-term exponential model f(x) = a*exp(b*x) + c*exp(d*x) + e*exp(f*x).
# exp3 = fittype('a*exp(b*x) + c*exp(d*x)+
# e*exp(f*x)','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b','c','d','e','f'});
# You then need to define reasonable start points if exp3 is to be used and
# input exp3 and not 'exp3'. This type of fitting wasn't done as part of
# this script because defining the starting points of six different
# coefficients gets rather complicated and, for the most part, the fits
# available in the Curve Fitting Toolbox are sufficient.

def _exp1( x, a, b ):
	func = a*np.exp(b*x)
	return func 

def _exp2( x, a, b, c, d ):
	func = a*np.exp(b*x) + c*np.exp(d*x)
	return func

def _power1( x, a, b ):
	func = a*np.power(x, b)
	return func

def _power2( x, a, b, c ):
	func = a*np.power(x, b) + c
	return func

func = _power1

# Fit the new excluded data xdata2 (eV) - the fit in single quotations can
# be changed as appropriate.
popt, pcov = curve_fit( func, data2[:,0][fitmask], data2[:,1][fitmask]  )

# Get residuals from fit
residuals = data2[:,1] - func(data2[:,0], *popt)

# Fitting statistics.
ss_res = np.sum( np.square(residuals) )
ss_tot = np.sum( data2[:,1] - np.square( np.mean(data2[:,1]) ) )
r_square = 1 - (ss_res/ss_tot)
p_sigma = np.sqrt(np.diag(pcov))

_print_best_fit( popt, p_sigma, 1 )
print( 'R2 = ' + str(r_square) )

# Approximate the signal integral using the trapezoidal rule. This can be
# used to quantify the signal-to-noise ratio (SNR) as well as in further
# quantification of chemical composition. The partial inelastic scattering
# cross-section needed for absolute quantification depends on the system
# being studied. The window of integration can be defined by the user.
# Define the start of integration for signal integral.

startint = 200
endint = 400

dataint = data2.copy()
dataint[:,1] = residuals
dataint = _extract_roi( data2, startint, endint )

# Integrate signal
ik = trapezoid( dataint[0], dataint[1] )

# Integrate background
bkgint = dataint.copy()
bkgint[:,1] = func( bkgint[:,0], *popt )

ib = trapezoid( bkgint[:,0], bkgint[:,0] )

# Variance in background signal.
varib = np.var( bkgint[:,1] )
# H parameter.
h = (ib+varib)/ib
# Signal-to-noise ratio.
snr = np.sqrt( ik/((ik+(h*ib)) ) )


# Store the results as a dictionary.
results_dict = { 'h' : h,\
				'snr' : snr,\
				'Background variance' : ib,\
				'R square' : r_square,\
				'Sigma' : p_sigma,\
				'Fitting parameters' : popt}

## Figure plotting
# Plot fit with confidence bounds, original EELS data, and the subtracted
# spectrum. The confidence bounds for the fitted coefficients determine the
# accuracy of the fit. Bounds that are far apart indicate uncertainty in
# the fit.


fig, axs = plt.subplots( 1,2 )
# Plot fit with confidence bounds and original data

#axs.plot(data1[:,0], data1[:,1], lw=2, label='data1');
axs[0].plot(data2[:,0],\
			data2[:,1],\
			lw=1,\
			label='Orignal EELS data',\
			color='k');
axs[0].plot(data2[:,0],\
			residuals,\
			label='Subtracted spectrum',\
			lw=1,\
			color='r')
axs[0].plot( data2[:,0],\
			func( data2[:,0], *popt ),\
			lw=1,\
			label='Fitted curve',\
			color='b'  )

axs[1].plot(data2[:,0][fitmask],\
			_normalise_data_range(residuals[fitmask],-1,1),\
			'.',\
			color='b',\
			lw = 1,\
			ms=2,\
			label='Residuals of fit')
axs[1].axhline(0, color='red',linestyle='--', label='Zero line')

#axs.plot(residuals[:,0], residuals[:,1], label='residuals')


#axs.plot( dataint[:,0], dataint[:,1], label='dataint' )
#axs.plot( bkgint[:,0], bkgint[:,1], label='bkgint' )

#axs.plot(ib, data2[:,0], label='ib')

#axs.set_ylim([2e7, 5e7])

for ax in axs:
	ax.legend(loc='upper right')
	ax.set_xlabel('Energy loss / eV')
	ax.set_ylabel('Counts')
	ax.grid(True)
	ax.set_box_aspect(1)

axs[0].set_xlim( [data2[0,0], data2[-1,0]] )
axs[1].set_xlim( [startedge, i] )

axs[0].set_ylim( [None, np.max(data2[:,1])] )
axs[1].set_ylim( [None, None] )

axs[0].set_title( ('Background fit of data below ' + str(i) + ' eV') )
axs[0].set_title( ('Residuals of fit') )

plt.show() 