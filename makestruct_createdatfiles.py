"""We create dat files from the fits files of the object ccd1-4 files containing columns of wl, sob, uob,
# (observed signal flux, flux uncertainty) then concatenate them into one to loop through to compare to segm file to
# see which lines are ones we're interested in. The bigger cannon files has the variable data like grav, NOT this stuff.
We also create and use the autonorm function, which basically uses a polyfit to fit a polynomial, then remove points
that are too far away from this; the outliers"""




import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

# Import the wavelength, flux, and error from the ccd readings
from galah_sp_part2 import total_ccd_wavelength, total_ccd_sob_flux, total_ccd_flux_error_uob

# Using the resolution file to grab the linemask we've made
from makestruct_resolution import linemask_sp, segmask_segm, c
# Taking radial velocity from the file which takes it from cannon fits file. Used for doppler shift.
from galah_sp_part3 import radial_velocity_global

# Shifts the wavelengths to account for the doppler effect.
total_ccd_wavelength_dopplered = (total_ccd_wavelength)/((radial_velocity_global[0]*(1E3/c)) + 1E0)

# Skipping the "if not normalising, make sure there are lines in the segment" part as I think we can do that here
# Unsure how many there will be, so arrays can't be used. Instead we use a temp list.
ccd_wavelengths_inside_segmask = []
ccd_flux_inside_segmask = []
ccd_flux_error_inside_segmask = []
# array for the first and final wavelength indexes of each segment respectively. Remember SME only wants the one col.
wavelength_start_end_index = np.zeros((len(segmask_segm["Wavelength_Start"]), 2))

# For each segment in segmask, find the values of dopplered wavelength (and associated flux from indexing) that are inside.
# This is the ;select wavelength regions inside segments      part of makestruct
for segment in range(0,len(segmask_segm["Wavelength_Start"])):
    # Beginning wavelength and end of that segment. Put as variables here for readability.
    seg_start = (pd.to_numeric(segmask_segm["Wavelength_Start"][segment]))
    seg_stop = (pd.to_numeric(segmask_segm["Wavelength_End"][segment]))
    # Finding the index of values inside the segment, using logical and is a neccesity.
    wavelength_inside_segmask_index = np.where(np.logical_and(seg_stop >= total_ccd_wavelength_dopplered, total_ccd_wavelength_dopplered >= seg_start))
    ccd_wavelengths_inside_segmask.extend(total_ccd_wavelength_dopplered[wavelength_inside_segmask_index])
    ccd_flux_inside_segmask.extend(total_ccd_sob_flux[wavelength_inside_segmask_index])
    ccd_flux_error_inside_segmask.extend(total_ccd_flux_error_uob[wavelength_inside_segmask_index])
    # Numpy array of indexes of the first and final wavelengths per segment with column 0 being the first.
    # the wind SME wants seems to just be the final values so keep a note of that @@@@@@@@@@@@@@@
    wavelength_start_end_index[segment, 0]= (wavelength_inside_segmask_index[0][0])
    wavelength_start_end_index[segment, 1]= (wavelength_inside_segmask_index[-1][-1])

# Turning lists into arrays for numpy indexing with np.where. Will this slow things down too much by having both list
# and array? Delete one? Or does the creation take too long?
ccd_wavelengths_inside_segmask_array = np.array(ccd_wavelengths_inside_segmask)
ccd_flux_inside_segmask_array = np.array(ccd_flux_inside_segmask)
ccd_flux_error_inside_segmask_array = np.array(ccd_flux_error_inside_segmask)


# Number of segments that contain visible spectra. Wind (wavelength_start_end_index)
# is the array of indexes of first and final wavelengths per segment.
number_of_segments = len(wavelength_start_end_index)
if number_of_segments == 0:
    print("No observations in segment mask")
    exit()

# Starts at 1 and only changes in sp after normalization (to -3) for some damn reason. Very convoluted. @@@@@@@@@@@@@
continuum_scale_flag = 1
# selects contiuum points or something. Something to do with running_snr. Placeholder before the cscale_flag changes
if continuum_scale_flag >= 0:
    pass

""" I can't find this stuff anywhere in make struct or galahsp?! @@@@@@@@@@@
# We want to set the initial scale to 1 (2 values for each segment, for flag = 1) what does this mean @@@@@@@@@@@@@@@@@@
if continuum_scale_flag == 1:
    continuum_scale = np.ones((2, number_of_segments))
# When would it EVER be 0?? @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
elif continuum_scale_flag == 0:
    continuum_scale = np.ones((number_of_segments))
# Why would it not need an array? What changes when it's been normalized? @@@@@@@@@@@@@@@@@@@@@@@@@@@
elif continuum_scale_flag < 0:
    continuum_scale_flag = 1
"""
# wran is segment_begin_end
# An array with two columns, the first and last recorded wavelength in each segment
# Why do we have both the index AN D the values in two separate arrays? @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
segment_begin_end = np.copy(wavelength_start_end_index) # copy to avoid overwriting wind. Just using it for its size.
for windex_row in range(0, number_of_segments):
    # At indexes 0,0 and 0,1 (and then 1,0etc) of the index array 'wavelength_start_end_index' we take the value and
    # apply it to the wavelngtharray as the values we have taken are indexes of the first and last wavelength of each
    # segment. windrow, 0 is the segment beginning. , 1 is the end.
    segment_begin_end[windex_row, 0] = total_ccd_wavelength_dopplered[int(wavelength_start_end_index[windex_row, 0])]
    segment_begin_end[windex_row, 1] = total_ccd_wavelength_dopplered[int(wavelength_start_end_index[windex_row, 1])]

import matplotlib.pyplot as plt
# Need to check with Karin and Sven that this works appropriately @@@@@@@@@@@@@@@@@@@@@
                        # wavelegnth x   flux y  order of fn( 0=linr)
                        # We use fit_param to modify the index selections, but it is a tuple here. So we use its 1st value
def autonormalisation(wavelength_array, flux_array, polynomial_order, fit_parameters):

    if polynomial_order == 0: polynomial_order = 2
    if fit_parameters == 0: fit_parameters = 1.5
    # Using np.polyfit to replcae the idl Robust_poly_fit, is this a fair replacement? @@@@@@@@@@@@@@@@
    for polyfit_loop in range(0, 99):

        # Gets the coefficients for a fit of the order polynomial_order
        polyfit_coefficients = np.polynomial.polynomial.polyfit(wavelength_array,flux_array,polynomial_order)

        # Uses these to get an array of the y values of this line
        fitted_flux = np.polynomial.polynomial.polyval(wavelength_array, polyfit_coefficients)

              # Creates an array of the error (sigma) of the line compared to original data, to find outliers by getting the
        # standard deviation of the difference
        fitted_sigma = np.std(flux_array - fitted_flux)

        # Find the fluxes that exist below the linear fit + error but above the lower error boundary (* p or something)
        # So not outliers, but inliers :) we take the first value from the output which is the index array itself
        inlier_index = (np.where(np.logical_and(flux_array < (fitted_flux+(2*fitted_sigma)),
                                               (flux_array > (fitted_flux-(fit_parameters*fitted_sigma))))))[0]
        # cont flux is c in idl.
        if polynomial_order == 2:
            continuous_function = polyfit_coefficients[0] + (polyfit_coefficients[1]*ccd_wavelengths_inside_segmask_array[segment_indexes])\
                                                    +(polyfit_coefficients[2]*ccd_wavelengths_inside_segmask_array[segment_indexes]**2)
        if polynomial_order == 1:
            continuous_function = polyfit_coefficients[0] + (polyfit_coefficients[1]*ccd_wavelengths_inside_segmask_array[segment_indexes])
        # Stops when no more convergence occurs I suppose, breaks the loop. Again, np where gives a tuple with dtype
        # the second condition uses the original non edited wavelength array
        if len(inlier_index) == len(wavelength_array) or len(inlier_index)/len(ccd_wavelengths_inside_segmask_array[segment_indexes]) <= 0.1:
            print("Converged")
            break

        # Replace the array with only values that lie inside the error boundaries we have.
        wavelength_array = wavelength_array[inlier_index]
        flux_array = flux_array[inlier_index]

    # This is 'stopp' in idl.
    if polynomial_order == 2:
        # co in idl, compared to the continuous flux of c. These variables, man.. Does not get returned in make struct
        # Additional just means there's something different that I don't know.
        continuous_function_additional = polyfit_coefficients[0] + polyfit_coefficients[1]*wavelength_array[inlier_index] +\
                                     polyfit_coefficients[2]*wavelength_array[inlier_index]**2
    elif polynomial_order == 1:
        continuous_function_additional = polyfit_coefficients[0] + polyfit_coefficients[1]*wavelength_array[inlier_index]

    # Something to do with further sigma calculations again  maybe, idl has it as sn
    sigma_flux_comparison = 1/(np.std(flux_array[inlier_index]/[continuous_function_additional]))
    #plt.plot(wavelength_array, flux_array)
    #plt.plot(wavelength_array, fitted_flux)
    #plt.show()
    return polyfit_coefficients, continuous_function, sigma_flux_comparison

# Pre-normalization steps. ";Performs first-guess normalisation by robustly converging straight line fit to high pixels"
# I removed the requirement for sob > 0. @@@@@@@@@@
for segment_band in range(0,number_of_segments):
    # Finds the index where the wavelengths are between the start and end of each segment, to be able to loop each seg as i
    # We do this in wavelength_inside_segmask_index already, but we know that the IDL code is so damn convoluted
    # then maybe a case occurs where we didn't make that. Doesn't make any sense to me, but is possible.
    # If this is slow, we can optimise somehow probably. @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    segment_indexes = (np.where(np.logical_and(ccd_wavelengths_inside_segmask >= segment_begin_end[segment_band, 0],
                                               ccd_wavelengths_inside_segmask <= segment_begin_end[segment_band, 1])))[0]

    # If count is greater than 20, do this. len segindex is the number of values that fit our criteria.
    if len(segment_indexes) > 20:

        # Take the coefficients of the polyfit, and the flux of it too using the equation from IDL
        polyfit_coefficients, continuous_function, sigma_flux_comparison = autonormalisation(ccd_wavelengths_inside_segmask_array[segment_indexes],
                                                                  ccd_flux_inside_segmask_array[segment_indexes], 1, 0)
        # Puts it to a relative flux out of 1 that we see in the abundance charts.
        ccd_flux_inside_segmask_array[segment_indexes] = ccd_flux_inside_segmask_array[segment_indexes]/continuous_function
        ccd_flux_error_inside_segmask_array[segment_indexes] = ccd_flux_error_inside_segmask_array[segment_indexes]/continuous_function
        #print(ccd_flux_inside_segmask_array[segment_indexes])

    # If we don't have enough points, we just use the mean value instead. Numpy mean wasn#t working. Still need to
    # Double check this is correct, as it was wrong for seemingly no reason before @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    else:
        # np.mean had issues with our results. Making a variable for readability.
        flux_mean = (sum(ccd_wavelengths_inside_segmask_array[segment_indexes]))/len(ccd_wavelengths_inside_segmask_array[segment_indexes])

        # Gotta be ordered correctly or we modify the sob before we use it to modify uob!
        # Is this the right way around? @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        ccd_flux_error_inside_segmask_array[segment_indexes] = \
            ccd_flux_error_inside_segmask_array[segment_indexes]/(flux_mean)

        ccd_flux_inside_segmask_array[segment_indexes] = \
                (ccd_flux_inside_segmask_array[segment_indexes])/(flux_mean)

# end of prenormalisation! Phew!



