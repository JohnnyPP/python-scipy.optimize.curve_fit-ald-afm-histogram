import numpy as np 
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit

filename = '/home/kolan/mycode/python/ald/1cycle2_bcr_hist.asc'

def FindHeaderLength():
    """
    Finds the positionon the last header item line 
    This is then used by numpy.loadtxt thats wrides items to np.array
    """

    lookup = '#------------------------------'
    
    with open(filename) as myFile:
        for FoundPosition, line in enumerate(myFile, 1):
            if lookup in line:
                #print 'Scan Data found at line:', FoundPosition
                break
    
    return FoundPosition


x=np.loadtxt(filename,dtype=float,delimiter='\t',skiprows=FindHeaderLength(),usecols=(0,))
y=np.loadtxt(filename,dtype=float,delimiter='\t',skiprows=FindHeaderLength(),usecols=(1,))


######################################
##Fitting
######################################

# Define model functions to be used to fit to the data above:
def gauss(x, A, mu, sigma):
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def doubleGauss(x, A, A2, mu, mu2, sigma, sigma2):
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))+A2*np.exp(-(x-mu2)**2/(2.*sigma2**2))


#p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
#p0 = [0.01, -0.006,  0.4]
#p0 is the initial guess for the fitting coefficients (A, A2, mu, mu2, sigma, sigma2 above)
p0 = [0.005, 0.006,  -0.3, 0.3, 0.2, 0.35]

# Fit an the gauss curve to the data
#popt, pcov = curve_fit(gauss, x, y, p0=p0)
#popt, pcov = curve_fit(gauss, x, y)

#print popt
#print pcov

# Fit an the doubleGauss curve to the data
popt2, pcov2 = curve_fit(doubleGauss, x, y, p0=p0)
print popt2
print pcov2

# Get the fitted curve
#yFit = gauss(x, *popt)
yFit = doubleGauss(x, *popt2)


def rsquared(x, y):
    """ Return R^2 where x and y are array-like."""

    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    return r_value**2
    
    
def gauss1(x):
    return popt2[0]*np.exp(-(x-popt2[2])**2/(2.*popt2[4]**2))

def gauss2(x):
    return popt2[1]*np.exp(-(x-popt2[3])**2/(2.*popt2[5]**2))


plt.figure('1st ALD cycle')
plt.plot(x, y, 'bo', label='Data')
plt.plot(x, yFit, linewidth=2, label='Fit', color='red')
plt.plot(x, gauss1(x), label='SiO$2$', color='green')
plt.plot(x, gauss2(x), label='HfO$2$', color='blue')
plt.fill_between(x, gauss2(x), color='blue', alpha=0.3)
plt.fill_between(x, gauss1(x), color='green', alpha=0.3)
plt.title('Raw data plot')
plt.xlabel('Height [nm]')
plt.ylabel('Frequency [counts]')
plt.grid(True)
plt.legend()
print rsquared(y, yFit)
plt.show()