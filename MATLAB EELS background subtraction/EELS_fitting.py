'''
Python version of EELS fiting script.

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
	rows = np.where( inarr[:,0] > start )
	temp = inarr[rows]
	rows = np.where( temp[:,0] < end )
	outarr = temp[rows]
	return outarr


# Normalsie range of array.
def _normalise_data_range( data, dmin=0, dmax=1 ):
    return ((data-np.min(data))/(np.max(data)-np.min(data)))*( dmax - dmin )

# Import .msa data exported from Digital Micrograph 3.
filename = 'file_name.msa'; # This is your file name.

data = _import_data_from_file( filename, ',' )

# There may be situations where the data of interest is a smaller section
# of the total spectrum. The code below allows the user to define the start
# and end of that data. A larger background before the start of the edge
# tends to lead to better fitting of the background.

# Extracting the edge.
# Define the start of the edge.
startedge = 176
endedge = 381
data2 = _extract_roi( data, 176, 381 )

# Fitting 'for' loop
# Define fit - use one of the fits in single quotations below:

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

# (i = start:increment:end for excluded data points)
i = range(200, 280, 10)

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

fig, ax = plt.subplots()

for j in i:
	# Define fit used and 'Exclude' data points in xdata (eV) for fitting
	fitdata = _extract_roi( data2, j, endedge )
	popt, pcov = curve_fit( func, fitdata[:,0], fitdata[:,1]  )
	# Get residuals from fit to plot subtracted spectra
	residuals = fitdata[:,1] - func(fitdata[:,0], *popt)


	ax.plot( fitdata[:,0],\
			_normalise_data_range(residuals,dmin=0, dmax=np.max(data2[:,1])),\
			label=str(j) )


ax.plot( data2[:,0],\
		_normalise_data_range(data2[:,1],dmin=0,dmax=np.max(data2[:1])),\
		color='k',\
		label='Orignal EELS data' )

ax.legend()

ax.set_xlim(startedge, endedge)

plt.show()

# End of script.