'''
Python version of EELS_fitting.m script.

Author: E Weare
Location: nmRC
Contact: benjamin.weare1@nottingham.ac.uk

Please cite if you found this script usesful.
'''

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


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
i = range(280, 230, -10)

fig, ax = plt.subplots()

ax.plot( data2[:,0],\
		data2[:,1],\
		color='k',\
		label='Orignal EELS data' )

co = [[1.00 , 0.40, 0.40],
	[1.00, 0.70, 0.40],
	[1.00, 1.00, 0.40],
	[0.70, 1.00, 0.40],
	[0.40, 1.00, 0.40],
	[0.40, 1.00, 0.70],
	[0.40, 1.00, 1.00],
	[0.40, 0.70, 1.00],
	[0.40, 0.40, 1.00],
	[0.70, 0.40, 1.00],
	[1.00, 0.40, 1.00],
	[1.00, 0.40, 0.70]]

n = 0
for j in i:
	fitmask = np.where( data2[:,0] < j )

	# Do the curve fit.
	try:
		popt, pcov = curve_fit( func, data2[:,0][fitmask], data2[:,1][fitmask]  )
	except:
		print('Could not fit.')
		continue

	# Get residuals from fit to plot subtracted spectra
	residuals = data2[:,1] - func(data2[:,0], *popt)
	pltlabel = str(j) + ' eV'
	try:
		ax.plot( data2[:,0],\
				residuals,\
				color=co[n],\
				label=pltlabel )
		n = n+1
	except:
		print('Out of custom colours.')
		ax.plot( data2[:,0],\
				residuals,\
				label=pltlabel )


ax.legend(bbox_to_anchor=(1.05, 1),\
			title='Fits excluding data above i eV',\
            loc='upper left', borderaxespad=0.)
ax.set_xlabel('Energy loss / eV')
ax.set_ylabel('Normalised Counts')
ax.grid(True)
ax.set_title('EELS spectra after subtracting fitted curves')

ax.set_box_aspect(1)
ax.set_xlim([startedge, endedge])

plt.show()

# End of script.