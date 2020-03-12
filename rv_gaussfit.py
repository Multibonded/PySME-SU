from galah_sp_part2 import total_ccd_wavelength, total_ccd_sob_flux
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Broad lines are the halpha and hbeta lines, wide and always visible.
# Issue is with the optional values, how to do em? Currently we require all arguments.@@@@@@@
def gauss_function(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def fit_func(x, a0, a1, a2, a3, a4, a5):
    z = (x-a1)/a2
    y =  a0*np.exp(-z**2/2) + a3
    print(y)
    return y

def rv_gaussfit(total_ccd_wavelength, total_ccd_sob_flux, broad_lines, mrv = 0 , srv = 0, erv = 0, rv_range = [2.5, 2.5], auto = False):
    rvtab = np.zeros(len(broad_lines))

    print(rvtab)

    for i in range(0, len(broad_lines)):
        while True:
            # If rv_range is not equal to the default set value do this, otherwise if it is(one time only I guess) do the other.
            if rv_range != [2.5, 2.5]:
                broad_line_index = np.where(np.logical_and(total_ccd_wavelength > broad_lines[i]-rv_range[0], total_ccd_wavelength < broad_lines[i]+rv_range[1]))[0]
            else:
                broad_line_index = np.where(np.logical_and(total_ccd_wavelength > broad_lines[i]-2.5, total_ccd_wavelength < broad_lines[i]+2.5))[0]

            # When can it ever be a value of -1 hmm
            if broad_line_index[0] == -1:
                if auto == False:
                    print("Line not covered in spectrum")
                    rv = -1

            else:
                # Check if this is good enough replacement for gaussfit from idl
                xdata= total_ccd_wavelength[broad_line_index]
                ydata = total_ccd_sob_flux[broad_line_index]
                n = len(xdata)  # the number of data
                mean = sum(xdata * ydata) / n  # note this correction
                sigma = sum(ydata * (xdata - mean) ** 2) / n  # note this correction
                parameters, covariance = curve_fit(fit_func, xdata, ydata)
                print("parameters:", parameters, "\n|||||\n", "covariance:", covariance)

                fitdata = fit_func(xdata, *parameters)
                print(fitdata)
                plt.plot(xdata, ydata)
                plt.plot(xdata, fitdata)

                plt.show()
range_def = [15,15]
lines = [4861.3230,6562.7970]
auto = True
rv_gaussfit(total_ccd_wavelength, total_ccd_sob_flux, lines ,auto)