'''
Python version of EELS fit analysis script.

Need to convert to eV for correct array indexing in functions below.
'''

import numpy as np

# Define functions.

# Import from .msa file.
def _import_data_from_file( datafile, delim ):
    data = genfromtxt(datafile, delimiter=delim, comments='#')
    array = np.array( data[:,0:1] )
    return array


# Extract spectrum region of interest.
def _extract_roi( x1, start, end ):
	x2 = x1[start:end]
	return x2


# Get the energy scale from the msa
def _get_scale( xdata ):
	scale = xdata[1] - xdata[0]
	return scale

# Import .msa data exported from Digital Micrograph 3.
filename = 'file_name.msa'; # This is your file name.

xdata, ydata = _import_data_from_file( filename, ',' )

# There may be situations where the data of interest is a smaller section
# of the total spectrum. The code below allows the user to define the start
# and end of that data. A larger background before the start of the edge
# tends to lead to better fitting of the background.

# Extracting the edge. Define the start of the edge.
startedge = 176
endedge = 381
xdata2 = _extract_roi( xdata1, startedge, endedge )
ydata2 = _extract_roi( ydata1, startedge, endedge )


# Analysing exponential fitting
# Exclude data from table where xdata is above i eV and change i to the
# appropriate value.
i = 280
xdata3, ydata3 = _extract_roi( xdata2, ydata2, 0, i )

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
	return = a*np.exp(b*x)

def _exp2( x, a, b, c, d ):
	return = a*np.exp(b*x) + c*np.exp(d*x)

def _power1( x, a, b ):
	return = a*np.power(x, b)

def _power2( x, a, b, c ):
	return = a*np.power(x, b) + c


# Fit the new excluded data xdata2 (eV) - the fit in single quotations can
# be changed as appropriate.
from scipy.optimise import curve_fit
#[f,gof,output] = 
popt, pcov = curve_fit( _exp1, xdata3, ydata3  )
# Get residuals from fit
residuals = #ydata2 - f(xdata2); f(xdata2, popt) - ydata
# Get fit options
fitoptions = fitoptions(f)

# Approximate the signal integral using the trapezoidal rule. This can be
# used to quantify the signal-to-noise ratio (SNR) as well as in further
# quantification of chemical composition. The partial inelastic scattering
# cross-section needed for absolute quantification depends on the system
# being studied. The window of integration can be defined by the user.
# Define the start of integration for signal integral.
startint = xdata2 > 284;
xdataint1 = xdata2(startint);
residualsint1 = residuals(startint);
# Define end of integration for signal integral.
endint = xdataint1 < 300;
xdataint2 = xdataint1(endint);
residualsint2 = residualsint1(endint);
# Integrate signal
ik = trapz(xdataint2,residualsint2)
# Integrate background
fityvalues = f(xdata2);
bkgint1 = fityvalues(startint);
bkgint2 = bkgint1(endint);
ib = trapz(xdataint2,bkgint2)
# Calculate variance in the background integral
varib = var(bkgint2)
# h parameter
h = (ib+varib)/ib
# Signal-to-noise ratio (SNR)
snr = ik/((ik+(h*ib))^0.5)

from scipy.integrate import trapezoid

xdataint = _extract_roi( xdata, 284, 300 )
residualsint = _extract_roi( residuals, 284, 300 )

ik = trapezoid( xdataint, residualsint )# check this one
fityvalues = 
bkgint = _extract_roi( fityvalues, 284, 300 )
ib = trapezoid( xdataint, bkgint )

varib = np.var( bkgint )
h = (ib+varib)/ib
snr = np.sqrt( ik/((ik+(h*ib)) ) )


## Figure plotting
# Plot fit with confidence bounds, original EELS data, and the subtracted
# spectrum. The confidence bounds for the fitted coefficients determine the
# accuracy of the fit. Bounds that are far apart indicate uncertainty in
# the fit.

import matplotlib.pyplot as plt

fig, axs = plt.subplots( (1,2) )
# Plot fit with confidence bounds and original data
axs[0].plot(f,'b',xdata2,ydata2,'k','predobs', lw=2);

# Plot subtracted spectrum
axs[0].plot(xdata2,residuals,'r', lw=2);

# Define characteristics of fit and original data axes
axs[0].XLim = [-inf inf]; % Limits of x-axis
axs[0].YLim = [-inf 500000]; % Limits of y-axis
axs[0].FontName = 'Calibri';
#axs[0].setTickDir = 'out';
#axs[0].TickLength = [0.005 0.005];
axs[0].grid(true)
#axs[0].Layer = 'bottom';
axs[0].set_title( ('Background of fit data below ' + str(i) + ' eV'), fontsize=30 )
axs[0].set_xlabel('eV')
axs[0].set_ylabel('Counts')
axs[1].legend('Original data','Fitted curve','Upper prediction bounds','Lower prediction bounds','Subtracted spectrum')


# Residuals from a fitted model are the differences between the original
# data and the fit to the original data. Assuming the model is correct, the
# residuals approximate the random errors. The model fits the data well if
# the residuals appear to behave randomly. If the residuals have a
# systematic pattern, the model does not fit the data very well. Plot
# residuals

axs[1].plot(f,xdata3,ydata3,'Residuals', lw=2);

# Define characteristics for residuals axes
axs[1].set_xlim = [-inf inf]; % Limits of x-axis
axs[1].set_ylim = [-inf inf]; % Limits of y-axis
#axs[1].TickDir = 'out';
axs[1].TickLength = [0.005 0.005];
axs[1].grid(true)
axs[1].Layer = 'bottom'
axs[1].set_title( 'Residuals of fit',fontsize=30 )
axs[1].set_xlabel('eV')
axs[1].set_ylabel('Counts')
axs[1].legend('Residuals', 'Zero line')

# End of script.