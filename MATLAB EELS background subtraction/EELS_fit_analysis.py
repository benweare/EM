'''
Python version of EELS fit analysis script.

Author: E Weare
Location: nmRC
Contact: benjamin.weare1@nottingham.ac.uk
'''

import numpy as np

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


def _calc_goodness_of_fit( pcov ):
	return

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
i = 350
data3 = _extract_roi( data2, data2[0,0], i )

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
from scipy.optimize import curve_fit

popt, pcov = curve_fit( func, data3[:,0], data3[:,1]  )

# Get residuals from fit
residuals = data2.copy()
residuals[:,1] = data2[:,1] - func(data2[:,0], *popt)

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

from scipy.integrate import trapezoid

startint = 200
endint = 400

dataint = _extract_roi( data2, startint, endint )
residualsint = _extract_roi( residuals, startint, endint )

# Integrate signal
ik = trapezoid( dataint, residualsint )
# Integrate background
fityvalues = data2.copy()
fityvalues[:,1] = func( data2[:,0], *popt )
bkgint = _extract_roi( fityvalues, startint, endint )
ib = trapezoid( dataint, bkgint )

# Variance in background signal.
varib = np.var( bkgint )
# H parameter.
h = (ib+varib)/ib
# Signal-to-noise ratio.
snr = np.sqrt( ik/((ik+(h*ib)) ) )


## Figure plotting
# Plot fit with confidence bounds, original EELS data, and the subtracted
# spectrum. The confidence bounds for the fitted coefficients determine the
# accuracy of the fit. Bounds that are far apart indicate uncertainty in
# the fit.

import matplotlib.pyplot as plt

fig, axs = plt.subplots( 1,2 )
# Plot fit with confidence bounds and original data

#axs.plot(data1[:,0], data1[:,1], lw=2, label='data1');
axs[0].plot(data2[:,0], data2[:,1], lw=2, label='orignal EELS data');
axs[0].plot(residuals[:,0], residuals[:,1], label='Subtracted spectrum')
axs[0].plot( data2[:,0], func( data2[:,0], *popt ), lw=2, label='fitted curve'  )

axs[1].plot(residuals[:,0], residualsint[:,1], '.', label='Residuals of fit')
axs[1].axhline(0, color='red', label='Zero line')

#axs.plot(residuals[:,0], residuals[:,1], label='residuals')


#axs.plot( dataint[:,0], dataint[:,1], label='dataint' )
#axs.plot( bkgint[:,0], bkgint[:,1], label='bkgint' )

#axs.plot(ib, data2[:,0], label='ib')

#axs.set_ylim([2e7, 5e7])

axs[0].legend()
axs[1].legend()

plt.show()

'''
# Plot subtracted spectrum
axs[0].plot(data2, residuals,'r', lw=2);

# Define characteristics of fit and original data axes
axs[0].set_xlim(None, None) # Limits of x-axis
axs[0].set_ylim(None, None) # Limits of y-axis
#axs[0].FontName = 'Calibri';
#axs[0].setTickDir = 'out';
#axs[0].TickLength = [0.005 0.005];
axs[0].grid(true)
#axs[0].Layer = 'bottom';
axs[0].set_title( ('Background of fit data below ' + str(i) + ' eV'), fontsize=30 )
axs[0].set_xlabel('eV')
axs[0].set_ylabel('Counts')
axs[0].legend('Original data','Fitted curve','Upper prediction bounds','Lower prediction bounds','Subtracted spectrum')


# Residuals from a fitted model are the differences between the original
# data and the fit to the original data. Assuming the model is correct, the
# residuals approximate the random errors. The model fits the data well if
# the residuals appear to behave randomly. If the residuals have a
# systematic pattern, the model does not fit the data very well. Plot
# residuals

axs[1].plot(f, data3, data3,'Residuals', lw=2);

# Define characteristics for residuals axes
axs[1].set_xlim(None, None) # Limits of x-axis
axs[1].set_ylim(None, None) # Limits of y-axis
#axs[1].TickDir = 'out';
#axs[1].TickLength = [0.005 0.005];
axs[1].grid(true)
axs[1].Layer = 'bottom'
axs[1].set_title( 'Residuals of fit',fontsize=30 )
axs[1].set_xlabel('eV')
axs[1].set_ylabel('Counts')
axs[1].legend('Residuals', 'Zero line')

plt.show()
# End of script.'''