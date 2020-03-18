"""
We read the segment mask segm, removing the comment lines and columns.
We proceed to open the ccd1-4 files and take their wavelength and resolution out to plot a resolution mask and
interpolating to be able to find resolution at peak wavelengths that we calculate at that point of the segment mask.
We open the 4 ccd fits files to produce a telephonenumber_Sp file that contains Wavelength array, observed flux sign (sob)
uncertainity in flux (uob)(essentially error), smod (not important now) and mod (Idk) that we call in make_struct
to modify the wavelength array by the radial velocity doppler shift and make some sme variables.
Running create_structure should be the only thing needed.
"""

from astropy.io import fits
import galah_sp_part1
import galah_sp_part3
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
# Import the wavelength, flux, and error from the ccd readings
import galah_sp_part2
import pickle
from datetime import date
import sme_rdlin
import makestruct_abund

# ccd_flux_inside_segmask_array is sob
# ccd_wavelengths_inside_segmask_array is wave

# Need to check with Karin and Sven that this works appropriately @@@@@@@@@@@@@@@@@@@@@
# wavelegnth x   flux y  order of fn( 0=linr)
# We use fit_param to modify the index selections, but it is a tuple here. So we use its 1st value
# We're inputting the wavelengths and fluxes of the readings in the first segment
def autonormalisation(wavelength_array, flux_array, polynomial_order, fit_parameters):
    if polynomial_order == 0: polynomial_order = 2
    if fit_parameters == 0: fit_parameters = 1.5
    # Using numpy polyfit to replcae the idl Robust_poly_fit, is this a fair replacement? @@@@@@@@@@@@@@@@
    for polyfit_loop in range(0, 99):
        # Gets the coefficients for a fit of the order polynomial_order
        polyfit_coefficients = np.polynomial.polynomial.polyfit(wavelength_array, flux_array, polynomial_order)

        # Uses these to get an array of the y values of this line
        fitted_flux = np.polynomial.polynomial.polyval(wavelength_array, polyfit_coefficients)

        # Creates an array of the error (sigma) of the line compared to original data, to find outliers by getting the
        # standard deviation of the difference
        fitted_sigma = np.std(flux_array - fitted_flux)

        # Find the fluxes that exist below the linear fit + error but above the lower error boundary (* p or something)
        # So not outliers, but inliers :) we take the first value from the output which is the index array itself
        inlier_index = (np.where(np.logical_and(flux_array < (fitted_flux + (2 * fitted_sigma)),
                                                (flux_array > (fitted_flux - (fit_parameters * fitted_sigma))))))[0]
        # cont flux is c in idl.
        if polynomial_order == 2:
            continuous_function = polyfit_coefficients[0] + (polyfit_coefficients[1] * wavelength_array) \
                                  + (polyfit_coefficients[2] * wavelength_array ** 2)
        if polynomial_order == 1:
            continuous_function = polyfit_coefficients[0] + (polyfit_coefficients[1] * wavelength_array)
        else:
            continuous_function = 1
        # Stops when no more convergence occurs I suppose, breaks the loop. Again, np where gives a tuple with dtype
        # the second condition uses the original non edited wavelength array
        if len(inlier_index) == len(wavelength_array) or len(inlier_index) / len(wavelength_array) <= 0.1:
            print("Converged after", polyfit_loop, "loops")
            break
        if polyfit_loop >= 98:
            print("Did not converge")
            break
        # Replace the array with only values that lie inside the error boundaries we have.
        wavelength_array = wavelength_array[inlier_index]
        flux_array = flux_array[inlier_index]
    # This is 'stopp' in idl.
    if polynomial_order == 2:
        # co in idl, compared to the continuous flux of c. These variables, man.. Does not get returned in make struct
        # Additional just means there's something different that I don't know.
        continuous_function_additional = polyfit_coefficients[0] + polyfit_coefficients[1] * wavelength_array[inlier_index] \
                                         + polyfit_coefficients[2] * wavelength_array[inlier_index] ** 2
    elif polynomial_order == 1:
        continuous_function_additional = polyfit_coefficients[0] + polyfit_coefficients[1] * wavelength_array[
            inlier_index]
    else:
        continuous_function_additional = 1
    # Something to do with further sigma calculations again  maybe, idl has it as sn
    sigma_flux_comparison = 1 / (np.std(flux_array[inlier_index] / [continuous_function_additional]))
    # import matplotlib.pyplot as plt
    # plt.plot(wavelength_array, flux_array)
    # plt.plot(wavelength_array, fitted_flux)
    # plt.show()
    return polyfit_coefficients, continuous_function, sigma_flux_comparison, wavelength_array


# Produces the arrays with the stellar information from the master file for every wavelength. such as line_atomic etc.
def indexed_stellar_information():
    # The reduced amount of spectral lines we need to look at, applied to masterfile, brings it down from 300k to 20k. much faster.
    # need int as it's an index array which requires that type. Un-needed if the data file is already created with rdlin
    data_index = np.loadtxt('alllinesindex.csv', delimiter=" ").astype(int)
    try:
        print("Line_merge data file found, this will be faster!")
        line_atomic, lande_mean, depth, data_reference_array, species, j_e_array, lower_level, \
        upper_level, lu_lande = pickle.load(open("galah_master_file_arrays", "rb"))
        line_atomic, lande_mean, depth, data_reference_array, species, j_e_array, lower_level, \
        upper_level, lu_lande = line_atomic[data_index], lande_mean[data_index], depth[data_index], \
                                data_reference_array[data_index], species[data_index], j_e_array[data_index], \
                                lower_level[data_index], upper_level[data_index], lu_lande[data_index]

    except FileNotFoundError:
        print("No line_merge data file created previously, running a line merger. This could take up to 2 minutes, and will"
              "create a data file for later use for any star with the ", galah_sp_part1.setup_for_obs, "release")
        # Will create the data file for use if it doesn't exist.
        sme_rdlin.run_merger()

        # Now we know it exists, we can do what we were trying to before.
        line_atomic, lande_mean, depth, data_reference_array, species, j_e_array, lower_level, \
        upper_level, lu_lande = pickle.load(open("galah_master_file_arrays", "rb"))
        line_atomic, lande_mean, depth, data_reference_array, species, j_e_array, lower_level, \
        upper_level, lu_lande = line_atomic[data_index], lande_mean[data_index], depth[data_index], \
                                data_reference_array[data_index], species[data_index], j_e_array[data_index], \
                                lower_level[data_index], upper_level[data_index], lu_lande[data_index]

    return line_atomic, lande_mean, depth, data_reference_array, species, j_e_array, lower_level, \
        upper_level, lu_lande


# Creates the segment mask and the resolution of it from the csv file that we name in the variable input.
def segment_mask_creation(segment_mask):


    # What is status. At one point during make_obs, status it says if status = 0 it's an error.
    # It gets set to 1 only at 1 point after selecting wavelength regions inside segments, before
    # starting the pre normalisation procedure.
    status = 0

    # Segm_mask is _Segm.data, unsurprisingly. It takes the start, end wavelength, and the resolution base guess of 35000
    # which is totally wrong. That's something to fix/modify when we're streamlining this.
    # The seperator separates via 2 spaces or a space and ; or . and space. We use the overflow column to account for lines
    # that begin with ; (commented out) which we immediately delete afterwards. Hard to separate comments when ; was forgotten.
    # Can't do lower case potentials because it randomly selects lower case letters in the middle of sentences. Something is wrong.
    # Same as in galah part 4, this should be a dummy file with only one segment, but we're trying to avoid that so we
    # are using the full regular segment mask NOT obsname.
    segment_mask_data_with_res = pd.read_csv(
        "DATA/" + segment_mask, sep='[ ]{2,}| ;| [A-Z]', header=None,
        names=["Wavelength_Start", "Wavelength_End", "Resolution", "comment", "overflow"], engine='python',
        skipinitialspace=True,
        usecols=["Wavelength_Start", "Wavelength_End",
                 "Resolution"])  # Here we auto delete the comment and overflow columns
    # ~ asks for negation, removes any row that starts with ; in the first column. Unsure what the str does but required.
    segment_mask_data_with_res = segment_mask_data_with_res[
        ~segment_mask_data_with_res['Wavelength_Start'].str.startswith(";")]
    # Reset the index to account for the now missing values
    segment_mask_data_with_res = segment_mask_data_with_res.reset_index(drop=True)
    # Next is: ;This sorts out segments 35 \n i=sort(seg_st) & seg_st=seg_st[i] & seg_en=seg_en[i] & seg_ipres=seg_ipres[i]
    # Sorts in ascending order of wavelength of starting wavelength, and re-orders the others.
    # But it's already ordered? weird. Will leave it in case data comes in unordered
    segment_mask_data_with_res.sort_values(by=['Wavelength_Start', 'Wavelength_End'])
    # It is a dummy file as of 27-02-2020 so it should only be one segment. (Why do we make atomic data for the whole thing
    # with each loop though? in sme rdlin )@ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #segment_mask_data_with_res = segment_mask_data_with_res.head(1)
    # segment_mask_data_with_res = segment_mask
    return segment_mask_data_with_res


# Creates an interpolation equation based on the resolutions we already have and their wavelgnths.
def resolution_interpolation(object_pivot):
    # mrdfits reads fits file of ccd1_piv etc. and saves the data as res1. The piv files are essentially the resolution files
    # I believe
    ccd1_res_file = fits.open(r'ccd1_piv.fits', ext=0)
    ccd1_res_data = (ccd1_res_file[0].data)
    ccd2_res_file = fits.open(r'ccd2_piv.fits', ext=0)
    ccd2_res_data = (ccd2_res_file[0].data)
    ccd3_res_file = fits.open(r'ccd3_piv.fits', ext=0)
    ccd3_res_data = (ccd3_res_file[0].data)
    ccd4_res_file = fits.open(r'ccd4_piv.fits', ext=0)
    ccd4_res_data = (ccd4_res_file[0].data)

    # We're making an array of the resolution of the (actually x) axis
    # Extracting the row of data of ccd1 that matches the piv number (-1 as piv starts at 1)
    # Piv is the 'strand'(?) of the ccd that we're looking at. @@
    ccd1_obj_piv_y = (ccd1_res_data[object_pivot - 1])
    # Creates a wavelength list from starting CRVAL1 (4700) in steps of CRDELT1 (0.1) until it matches the array len of NAXIS1
    ccd1_wavelength_x = ccd1_res_file[0].header['CRVAL1'] + (
            ccd1_res_file[0].header['CDELT1'] * np.arange(ccd1_res_file[0].header['NAXIS1']))
    #                               start                        steps                            number of data points
    # Remember these are all arrays.
    ccd2_obj_piv_y = (ccd2_res_data[object_pivot - 1])
    ccd2_wavelength_x = ccd2_res_file[0].header['CRVAL1'] + (
            ccd2_res_file[0].header['CDELT1'] * np.arange(ccd2_res_file[0].header['NAXIS1']))
    ccd3_obj_piv_y = (ccd3_res_data[object_pivot - 1])
    ccd3_wavelength_x = ccd3_res_file[0].header['CRVAL1'] + (
            ccd3_res_file[0].header['CDELT1'] * np.arange(ccd3_res_file[0].header['NAXIS1']))
    ccd4_obj_piv_y = (ccd4_res_data[object_pivot - 1])
    ccd4_wavelength_x = ccd4_res_file[0].header['CRVAL1'] + (
            ccd4_res_file[0].header['CDELT1'] * np.arange(ccd4_res_file[0].header['NAXIS1']))

    # combines the arrays together to interpolate them
    wavelength_res_x_collection = np.concatenate(
        [ccd1_wavelength_x, ccd2_wavelength_x, ccd3_wavelength_x, ccd4_wavelength_x])
    wavelength_piv_y_collection = np.concatenate([ccd1_obj_piv_y, ccd2_obj_piv_y, ccd3_obj_piv_y, ccd4_obj_piv_y])

    # iterpolate the res collections to find wlpeak for them later
    interpolation = interp1d(wavelength_res_x_collection, wavelength_piv_y_collection)
    return interpolation


# Uses our previously created resolution interpolation to interpolate at the peak wavelength (which we create)
def interpolate_peak_resolution(segment_mask_data_with_res, interpolation):
    # Loop through each wavelength range (consisting of a start and end) in the data and take the peak value and its res factor.
    # The idl code seems to ignore dum1 column, and instead calls seg_st and en
    # We make a list to avoid calling pandas data frame repeatedly, I believe this is faster and avoids any copy errors.
    temporarywavelengthpeaklist = []
    print("Here", segment_mask_data_with_res['Resolution'])
    for wavelengthband in range(0, len(segment_mask_data_with_res['Resolution'])):
        # Calculating the peak wavelength
        wlpeak = 0.5 * (float(segment_mask_data_with_res.loc[wavelengthband, 'Wavelength_Start'])
                        + float(segment_mask_data_with_res.loc[wavelengthband, 'Wavelength_End']))
        # Appending it to a list to add as a column to segment_mask_data_with_res
        temporarywavelengthpeaklist.append(wlpeak)
        # Interpolating the resolution at the wavelength peak and replacing the resolution of that index with it
        segment_mask_data_with_res.loc[wavelengthband, 'Resolution'] = interpolation(
            wlpeak) * galah_sp_part2.resolution_factor
    # We insert the peak list into the dataframe as a new column. I believe faster and easier than inserting a new row
    # each loop.
    segment_mask_data_with_res.insert(0, "Wavelength_Peak", temporarywavelengthpeaklist)
    return segment_mask_data_with_res


# Adjusts the wavelengths according to the doppler shift and returns it.
def doppler_wavelengths():
    # Setting speed of light
    c = 2.99792458E8

    # Shifts the wavelengths to account for the doppler effect.
    total_ccd_wavelength_dopplered = (galah_sp_part2.total_ccd_wavelength) / (
            (galah_sp_part3.radial_velocity_global[0] * (1E3 / c)) + 1E0)

    return total_ccd_wavelength_dopplered


# Finds the wavelengths that we have that are also inside the segments.
def wavelengths_flux_inside_segments(segment_mask_data_with_res, total_ccd_wavelength_dopplered):
    # Skipping the "if not normalising, make sure there are lines in the segment" part as I think we can do that here
    # Unsure how many there will be, so arrays can't be used. Instead we use a temp list.
    ccd_wavelengths_inside_segmask = []
    ccd_flux_inside_segmask = []
    ccd_flux_error_inside_segmask = []
    ccd_resolution_of_segment = []
    # array for the first and final wavelength indexes of each segment respectively. Remember SME only wants the one col.
    wavelength_start_end_index = np.zeros((len(segment_mask_data_with_res["Wavelength_Start"]), 2))

    # For each segment in segmask, find the values of dopplered wavelength (and associated flux from indexing) that are inside.
    # This is the ;select wavelength regions inside segments      part of makestruct
    # Despite having this array we still use np where each time to find the wavelengths in the segments
    for segment in range(0, len(segment_mask_data_with_res["Wavelength_Start"])):
        # Beginning wavelength and end of that segment. Put as variables here for readability.
        seg_start = (pd.to_numeric(segment_mask_data_with_res["Wavelength_Start"][segment]))
        seg_stop = (pd.to_numeric(segment_mask_data_with_res["Wavelength_End"][segment]))
        # Finding the index of values inside the segment, using logical and is a neccesity.
        wavelength_inside_segmask_index = np.where(
            np.logical_and(seg_stop >= total_ccd_wavelength_dopplered, total_ccd_wavelength_dopplered >= seg_start))
        ccd_wavelengths_inside_segmask.extend(total_ccd_wavelength_dopplered[wavelength_inside_segmask_index])
        ccd_flux_inside_segmask.extend(galah_sp_part2.total_ccd_sob_flux[wavelength_inside_segmask_index])
        ccd_flux_error_inside_segmask.extend(galah_sp_part2.total_ccd_flux_error_uob[wavelength_inside_segmask_index])
        # Numpy array of indexes of the first and final wavelengths per segment with column 0 being the first.
        # the wind SME wants seems to just be the final values so keep a note of that @@@@@@@@@@@@@@@
        wavelength_start_end_index[segment, 0] = (wavelength_inside_segmask_index[0][0])
        wavelength_start_end_index[segment, 1] = (wavelength_inside_segmask_index[-1][-1])
        ccd_resolution_of_segment.append(segment_mask_data_with_res['Resolution'][segment])
    # Turning lists into arrays for numpy indexing with np.where. Will this slow things down too much by having both list
    # and array? Delete one? Or does the creation take too long? no. probs not.
    ccd_wavelengths_inside_segmask_array = np.array(ccd_wavelengths_inside_segmask)
    ccd_flux_inside_segmask_array = np.array(ccd_flux_inside_segmask)
    ccd_flux_error_inside_segmask_array = np.array(ccd_flux_error_inside_segmask)
    ccd_resolution_of_segment_array = np.array(ccd_resolution_of_segment)
    return ccd_wavelengths_inside_segmask_array, ccd_flux_inside_segmask_array, ccd_flux_error_inside_segmask_array,\
           ccd_resolution_of_segment_array, wavelength_start_end_index


# Finds the first and final wavelengths of the segments that we have.
def find_segment_limits(wavelength_start_end_index, total_ccd_wavelength_dopplered):
    number_of_segments = len(wavelength_start_end_index)

    # wran is segment_begin_end
    # An array with two columns, the first and last recorded wavelength in each segment
    # Why do we have both the index AN D the values in two separate arrays? @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    segment_begin_end = np.copy(wavelength_start_end_index)  # copy to avoid overwriting wind. Just using it for its size.
    for windex_row in range(0, number_of_segments):
        # At indexes 0,0 and 0,1 (and then 1,0etc) of the index array 'wavelength_start_end_index' we take the value and
        # apply it to the wavelngtharray as the values we have taken are indexes of the first and last wavelength of each
        # segment. windrow, 0 is the segment beginning. , 1 is the end.
        segment_begin_end[windex_row, 0] = total_ccd_wavelength_dopplered[
            int(wavelength_start_end_index[windex_row, 0])]
        segment_begin_end[windex_row, 1] = total_ccd_wavelength_dopplered[
            int(wavelength_start_end_index[windex_row, 1])]
    return segment_begin_end

# make_obs in idl, sets up the important wavelength, flux, and error of the important segments of wavelengths we have.
# This is the first one to run I think.
def object_array_setup():

    # We use the full setup file instead of the dummy/temp obsname file as we're trying to avoid constand looping.
    segment_mask = galah_sp_part1.setup_for_obs + '_Segm.dat'
    segment_mask_data_with_res = segment_mask_creation(segment_mask)

    # Takes the last three digits and turns into an int. The last 3 represent the fibre of the ccd used that we want to
    # look at. Then we produce an interpolation equation for the resolutions.
    interpolation = resolution_interpolation(int(str(galah_sp_part1.object_for_obs)[-3:]))

    # Uses our resolution interpolatino to find the resolution at the peak wavelengths that we also find.
    segment_mask_data_with_res = interpolate_peak_resolution(segment_mask_data_with_res, interpolation)

    # Checks for a negative range which I assume would cause issues. Check it works @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    if min(segment_mask_data_with_res['Wavelength_End'] - pd.to_numeric(
            segment_mask_data_with_res['Wavelength_Start'])) <= 0:
        print("Segment has a negative range!")
        return

    # Checks for overlapping segments if there's more than one. Double check this works when you have data @@@@@@@@@@@@@@
    if len(segment_mask_data_with_res['Wavelength_End']) > 1:
        if max(segment_mask_data_with_res['Wavelength_End'][0:len(segment_mask_data_with_res['Wavelength_End'])]
               - pd.to_numeric(segment_mask_data_with_res['Wavelength_Start'][1:len(segment_mask_data_with_res['Wavelength_Start'])])) > 0:
            print("Overlapping segments")
            return

    # Adjust the wavelengths for the doppler shift.
    total_ccd_wavelength_dopplered = doppler_wavelengths()

    ccd_wavelengths_inside_segmask_array, ccd_flux_inside_segmask_array, ccd_flux_error_inside_segmask_array, \
    ccd_resolution_of_segment_array, wavelength_start_end_index =\
        wavelengths_flux_inside_segments(segment_mask_data_with_res, total_ccd_wavelength_dopplered)

    # Number of segments that contain visible spectra. Wind (wavelength_start_end_index)
    # is the array of indexes of first and final wavelengths per segment.
    number_of_segments = len(wavelength_start_end_index)
    if number_of_segments == 0:
        print("No observations in segment mask")
        return

    # Creates an array with the beginning and end wavelenghts of the segments. different to the start end array as that's indexes.
    segment_begin_end = find_segment_limits(wavelength_start_end_index, total_ccd_wavelength_dopplered)

    # No clue. I think just confirms it reached the end of the function? Seems ridiculous. Removed. @@
    #status = 1

     # Returning the variables used for prenormalisation and the rest of makestruct. The first three are the important
    # data. Segment begin end is just made from the data (beginning and end wavelengths), but I'd rather not remake it each time.
    return ccd_wavelengths_inside_segmask_array, ccd_flux_inside_segmask_array, ccd_flux_error_inside_segmask_array,\
           segment_begin_end, wavelength_start_end_index, ccd_resolution_of_segment_array


# This ist he function to prenormalise the line. Calls on autonormalise, removes data that is too far away from the
# point and uses the close data to normalised itself. Doesn't affect the non-inlier data which is confusing and weird
# as it is not normalised for them.?
def pre_normalise(ccd_wavelengths_inside_segmask_array, ccd_flux_inside_segmask_array, ccd_flux_error_inside_segmask_array,
                  segment_begin_end):
    # Number of segments we have. Made each function to avoid potential length errors.
    number_of_segments = len(segment_begin_end)

    # Pre-normalization steps. ";Performs first-guess normalisation by robustly converging straight line fit to high pixels"
    # I removed the requirement for sob > 0. @@@@@@@@@@
    for segment_band in range(0, number_of_segments):
        # Finds the index where the wavelengths are between the start and end of each segment, to be able to loop each seg as i
        # We do this in wavelength_inside_segmask_index already, but we know that the IDL code is so damn convoluted
        # then maybe a case occurs where we didn't make that. Doesn't make any sense to me, but is possible.
        # If this is slow, we can optimise somehow probably. @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        segment_indexes = (np.where(np.logical_and(ccd_wavelengths_inside_segmask_array >= segment_begin_end[segment_band, 0],
                                                   ccd_wavelengths_inside_segmask_array <= segment_begin_end[
                                                       segment_band, 1])))[0]

        # If count is greater than 20, do this. len segindex is the number of values that fit our criteria.
        if np.size(segment_indexes) > 20:

            # Take the coefficients of the polyfit, and the flux of it too using the equation from IDL and then applies
            # a ploynomial fit, removing outlying values until it no longer changes. The outlying value limit is hardcoded
            # Wew then use the equation fit to it to normalise it
            polyfit_coefficients, continuous_function, sigma_flux_comparison, wavelength_array = autonormalisation(
                ccd_wavelengths_inside_segmask_array[segment_indexes],
                ccd_flux_inside_segmask_array[segment_indexes], 1, 0)
            # The index of the wavelengths we want to keep as they are not outliers.
            inlier_index = np.isin(ccd_wavelengths_inside_segmask_array[segment_indexes], wavelength_array)
            # Puts it to a relative flux out of 1 that we see in the abundance charts.
            # Have to take relative indexes for cont_func to have the correct shape.
            ccd_flux_inside_segmask_array[segment_indexes[inlier_index]] = ccd_flux_inside_segmask_array[
                                                                 segment_indexes[inlier_index]] / continuous_function
            ccd_flux_error_inside_segmask_array[segment_indexes[inlier_index]] = ccd_flux_error_inside_segmask_array[
                                                                       segment_indexes[inlier_index]] / continuous_function
            #print(wavelength_array, "\n",ccd_flux_inside_segmask_array[segment_indexes])
            # print(ccd_flux_inside_segmask_array[segment_indexes])
            # import matplotlib.pyplot as plt
            # plt.plot(ccd_wavelengths_inside_segmask_array[segment_indexes], ccd_flux_inside_segmask_array[segment_indexes])
            # plt.show()
        # If we don't have enough points, we just use the mean value instead. Numpy mean wasn#t working. Still need to
        # Double check this is correct, as it was wrong for seemingly no reason before @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        else:

            # np.mean had issues with our results. Making a variable for readability.
            flux_mean = (sum(ccd_wavelengths_inside_segmask_array[segment_indexes])) / len(
                ccd_wavelengths_inside_segmask_array[segment_indexes])

            # Gotta be ordered correctly or we modify the sob before we use it to modify uob!
            # Is this the right way around? @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            ccd_flux_error_inside_segmask_array[segment_indexes] = \
                ccd_flux_error_inside_segmask_array[segment_indexes] / (flux_mean)

            ccd_flux_inside_segmask_array[segment_indexes] = \
                (ccd_flux_inside_segmask_array[segment_indexes]) / (flux_mean)
            # import matplotlib.pyplot as plt
            # plt.plot(ccd_wavelengths_inside_segmask_array[segment_indexes], ccd_flux_inside_segmask_array[segment_indexes])
            # plt.show()

    # end of prenormalisation! Phew! Again we return the modified observational data that has been normalised (and some
    # of it has not?) @@
    return ccd_wavelengths_inside_segmask_array, ccd_flux_inside_segmask_array, ccd_flux_error_inside_segmask_array


# Open up the line and conituum masks for making mob/flagged fluxes
def open_masks():
    # of opening it we're just gonna import it from sp4 as an array
    line_mask = galah_sp_part1.setup_for_obs + "_Sp.dat"  # Line masks for log(g) sensitive lines, i.e Fe Ti Sc

    # Reads out the columns which are centre/peak wavelength, start and end of wavelength peak (all simulated), and atomic number
    # Careful with the python enginer, it's slower. If we are looking at BIG data files this might be bad.
    # line0, line_st, and line_en in idl. they ignore atomic number for some reason @@@@@@@@@@@@@@@@@@@
    linemask_data = pd.read_csv("DATA/" + line_mask, sep="[ ]{2,}", header=None, engine='python',
        names=["Sim_Wavelength_Peak", "Sim_Wavelength_Start", "Sim_Wavelength_End", "Atomic_Number"])
    # Removes all rows that begin with ; as it's commented out. Required to be able to modify the data by column rather than
    # by row
    linemask_data = linemask_data[~linemask_data['Sim_Wavelength_Peak'].str.startswith(";")]
    # Reset the index to account for the now missing values
    linemask_data = linemask_data.reset_index(drop=True)

    continuum_mask = str(galah_sp_part1.setup_for_obs) + "_Cont.dat"
    continuum_mask_data = pd.read_csv("DATA/" + continuum_mask, sep="[ ]{1,}", header=None, engine='python',
                                      names=["Continuum_Start", "Continuum_End"])
    return continuum_mask_data, linemask_data


# Flag the fluxes that correspond to the peaks
def flag_flux_peaks(ccd_wavelengths_inside_segmask_array, ccd_flux_inside_segmask_array, segment_begin_end):

    # Number of segments we have. Made each function to avoid potential length errors.
    number_of_segments = len(segment_begin_end)
    # Opens the data files containing the line and continuum data. we use these limits to flag the continuum
    # and peak fluxes. We take only the 2nd variable as this is the linemask peak function, not the continuum.
    linemask_data = open_masks()[1]


    """this stuff might work or it might not?! CHECK IT SOMEHOW! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ vvvvv"""
    # mob in idl. Mask of observed pixels
    flagged_fluxes = np.zeros(len(ccd_wavelengths_inside_segmask_array))

    # We flag the fluxes that are probably peaks inside our segments that we care about that are probably atomic absorb.
    #for line_loop in range(0,len(linemask_data['Sim_Wavelength_Start'])):
    for line_loop in range(0, 5):  # to avoid long waits when testing

        # Usual case of making sure the wavelengths we want to use are in the lines. (They're already made sure to be in
        # the segement mask. What's the diff? It's to separate by segment for each loop.
        wavelengths_inside_linemask_index = np.where(np.logical_and(
            ccd_wavelengths_inside_segmask_array >= linemask_data['Sim_Wavelength_Start'][line_loop],
            ccd_wavelengths_inside_segmask_array <= linemask_data['Sim_Wavelength_End'][line_loop]))
        # running_snr in idl, sets values to 1 if they're below the max noise spike. This means we flag all good points
        # at 1 that are below spikes in noise to be used.
        signal_to_noise = []
        # We're trying to find signal to noise ratios I guess, where 1.5 is the limit for our noise? Averages over +/-4 values. Unsure about the 11 though
        for flux_row in range(0, len(ccd_flux_inside_segmask_array)):  # Should this not be uob? for 2nd len()
            #              Indexes the obs flux from ii-4 to ii + 4 (or limits of length)
            signal_to_noise.append(max([1 + 10 / np.mean(ccd_flux_inside_segmask_array[
                                                         max(0, flux_row - 4):min(flux_row + 4,
                                                                                  len(ccd_flux_inside_segmask_array))] /
                                                         galah_sp_part2.total_ccd_flux_error_uob[max(0, flux_row - 4):min(flux_row + 4,
                                                                            len(ccd_flux_inside_segmask_array))]), 1.5]))

        signal_to_noise = (np.array(signal_to_noise))
        if len(wavelengths_inside_linemask_index[0]) > 0:
            # print("flux: ", max(ccd_flux_inside_segmask_array[wavelengths_inside_linemask_index[0]]), "\n Snr:", max(signal_to_noise[wavelengths_inside_linemask_index[0]]))
            # If the flux exists and is less than the noise, set a marker to 1 to indicate this.
            # Is this working? @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            if min(ccd_flux_inside_segmask_array[wavelengths_inside_linemask_index[0]]) > 0 and \
                    max(ccd_flux_inside_segmask_array[wavelengths_inside_linemask_index[0]]) < max(
                signal_to_noise[wavelengths_inside_linemask_index[0]]):
                # 1 is a good thing! it means it fits nicely in the peak and the noise is nothing to worry about
                flagged_fluxes[wavelengths_inside_linemask_index[0]] = 1
    # Return our array of which wavelengths are flagged as possible peaks to be modified further in contiuum and more.
    return flagged_fluxes


# And the same as flagging but for continuum
def flag_continuum(ccd_wavelengths_inside_segmask_array, ccd_flux_inside_segmask_array, flagged_fluxes):
    """; cont mob  - continuum points selected between 0 and 1.2 or 1+3*sigma, where there are no line masks
    ;             avoid buffer zone at edges
    """
    # Opens the data files containing the line and continuum data. we use these limits to flag the continuum
    # and peak fluxes. We take only the 1st variable as this is the continuum.
    continuum_mask_data = open_masks()[0]



    # If the flag is above 0 we're resetting signal to noise. Seems silly to reset it.. Maybe skip for fluxrow in the
    # function above? @@@@@@@@@@@@@@@@@@@@@@@@@@@
    # Basically same as above but for finding contiuum fluxes instead of peak fluxes.
    if galah_sp_part1.continuum_scale_flag >= 0:

        # For each segment in the continuum file. Currently only one massive range so meh.
        for continuum_loop in range(0, len(continuum_mask_data['Continuum_Start'])):
            # A list to append to with the signal to noise ratios, either 1.5 or 1 + 10/mean(flux/error) from i to ii -/+4.
            # So the higher the
            signal_to_noise = []
            for flux_row in range(0, len(ccd_flux_inside_segmask_array)):
                signal_to_noise.append(max([1 + 10 /
                                            np.mean(ccd_flux_inside_segmask_array[max(0, flux_row - 4):min(flux_row + 4, len(ccd_flux_inside_segmask_array))]
                                                    /  galah_sp_part2.total_ccd_flux_error_uob[max(0, flux_row - 4):min(flux_row + 4, len(ccd_flux_inside_segmask_array))]), 1.5]))
            # Indexes where the wavelengths are inside the continuum.
            wavelengths_inside_continuum_index = np.where((np.logical_and(
                ccd_wavelengths_inside_segmask_array >= continuum_mask_data['Continuum_Start'][continuum_loop],
                ccd_wavelengths_inside_segmask_array <= continuum_mask_data['Continuum_End'][continuum_loop])))
            # This is cleaner code (albeit slower I imagine) than having 5 np logical ands in the above statement.
            # If we haven't already flagged the flux as peak flux, and it's less than the noise (?) then we flag it with '2'
            # to mark it as continuum
            if len(wavelengths_inside_continuum_index[0]) > 0:
                for continuum_index in wavelengths_inside_continuum_index[0]:
                    if flagged_fluxes[continuum_index] != 1 and ccd_flux_inside_segmask_array[continuum_index] > 0 and \
                            ccd_flux_inside_segmask_array[continuum_index] < signal_to_noise[continuum_index]:
                        flagged_fluxes[continuum_index] = 2
        print("Before cont2")
        return flagged_fluxes


# This function removes the lowest 70% of fluxes while retaining enough points on both sides of peak to have a continuum
def cutting_low_flux(ccd_wavelengths_inside_segmask_array, ccd_flux_inside_segmask_array, segment_begin_end, flagged_fluxes):
    # Number of segments we have. Made each function to avoid potential length errors.
    number_of_segments = len(segment_begin_end)
    for segment_band in range(0, number_of_segments):

        ";Deselect 70% lowest continuum points in each segment using synthesis. Ensure both ends have continuum points"
        # The fraction of points we want to remove. Potentially modified in the while loop hence kept outside of it.
        fraction = 0.7
        # print("new segment")
        # While the fraction to use is not 0 we continue looping. Python is setting it to E-17 in loop so this fixes it
        # by using 0.01 instead of 0.
        while fraction > 0.01:
            # We take the wavelength indexes of the wavelengths that exist in the segments as always. We just repeat it
            # each time to be able to work on each segment individually I suppose? @@@@@@@@@
            wavelength_inside_segment_index = np.where(
                np.logical_and(ccd_wavelengths_inside_segmask_array >= segment_begin_end[segment_band, 0],
                               segment_begin_end[segment_band, 1] >= ccd_wavelengths_inside_segmask_array))
            # The value of the fraction of the fluxes we chose, so how many IS 70% for example.
            value_of_fraction = int(len(ccd_flux_inside_segmask_array[wavelength_inside_segment_index]) * fraction)
            # Takes the index of the 70%th lowest value (our cut off point) from a sorted list. Sorting takes a long
            # time but this doesn't appear to be the bottleneck. @@@
            cutting_flux_value = sorted(ccd_flux_inside_segmask_array[wavelength_inside_segment_index])[value_of_fraction]
            # We take a list of indexes where the flux in the segment is below our cut off point.
            cutting_flux_index = np.where(
                ccd_flux_inside_segmask_array[wavelength_inside_segment_index] < cutting_flux_value)
            # We need to count how many values there are at the extreme ends of the segment, as we need a continuum on both
            # Here we see how many are in the first 1/3, and how many in the last 2/3.
            # Also makes sure that the flux they have is flagged as continuum
            # Also have to be flagged as 2 - so not peaks.
            # Is this ok? @@@@@@@@@@@@@@@@@@@@@@@
            low_continuum_index = np.where(
                np.logical_and(ccd_wavelengths_inside_segmask_array[cutting_flux_index] <=
                               (ccd_wavelengths_inside_segmask_array)[
                                   int(len(wavelength_inside_segment_index[0]) / 3)],
                               (flagged_fluxes[cutting_flux_index] == 2)))
            high_continuum_index = np.where(
                np.logical_and(ccd_wavelengths_inside_segmask_array[cutting_flux_index] >=
                               (ccd_wavelengths_inside_segmask_array)[
                                   int(len(wavelength_inside_segment_index[0]) * 2 / 3)],
                               (flagged_fluxes[cutting_flux_index] == 2)))

            # print("number of low high seg", len(low_continuum_index[0]), len(high_continuum_index[0]))
            # print("wavelengths",sorted(ccd_wavelengths_inside_segmask_array[wavelength_inside_segment_index][cutting_flux_index[:value_of_fraction]]))

            # If we don't have enough points, decrease the fraction we remove.
            if len(low_continuum_index[0]) < 5 or len(high_continuum_index[0]) < 5:
                fraction -= 0.1
                # print("Adjusting fraction down 10% to", fraction)
                # if fraction < 0.1:
                # print("number of low high seg", len(low_continuum_index[0]), len(high_continuum_index[0]))

                # print(len(wavelength_inside_segment_index[0]))
            # If we have enough points on both sides, we can continue and remove them by looping through the indexes of
            # the low fluxes. These are the indexes OF the wavelength indexes so we then have to apply them to the wave
            # inside segment as seen below.
            else:
                for index_loop in wavelength_inside_segment_index[0][cutting_flux_index]:
                    # Checks if it's 2, as we don't want to remove spectra values that are at 1.
                    if flagged_fluxes[index_loop] == 2:
                        # print("2")
                        flagged_fluxes[index_loop] = 0

                """print("wl indexes",wavelength_inside_segment_index[0])
                print("cutting fluxindexes", sorted(cutting_flux_index[0]))
                print("applying those indexes to the wl indexeslol", wavelength_inside_segment_index[0][cutting_flux_index])
                print("the flagged flux of those", flagged_fluxes[wavelength_inside_segment_index[0][cutting_flux_index]])
                print("What is printed", flagged_fluxes[wavelength_inside_segment_index])
                print("high flux", flagged_fluxes[cutting_flux_index[value_of_fraction:]])
                print("low flux", flagged_fluxes[cutting_flux_index[:value_of_fraction]])"""

                """import matplotlib.pyplot as plt
    
                plt.plot(ccd_wavelengths_inside_segmask_array[wavelength_inside_segment_index],
                         flagged_fluxes[wavelength_inside_segment_index])
                plt.plot(ccd_wavelengths_inside_segmask_array[wavelength_inside_segment_index],
                         ccd_flux_inside_segmask_array[wavelength_inside_segment_index])
    
                plt.show()"""
                # Signal to break the while loop and continue iteration
                fraction = 0

                """this stuff might work or it might not?! CHECK IT SOMEHOW! @@@@@@@@@@@@@@@@@@@@@@@@@@@@^^^^^^^^^^^^"""

    return flagged_fluxes



# Honestly unsure what this does. soemthing to do with cores and "avoiding strong NLTE-dominated cores in mask" @@@@
def removing_nlte(ccd_wavelengths_inside_segmask_array, ccd_flux_inside_segmask_array, flagged_fluxes):
    # ;Avoid strong NLTE-dominated cores in mask, what. @@@@@@@@@@@@@
    core_minimum = 0.6
    if galah_sp_part3.metalicity < -2:
        core_minimum = 0.72
    elif galah_sp_part3.metalicity < -1:
        core_minimum = 0.65

    # Checks for line_cores existing from galahsp1, if no variable exists then we can't do this can we.
    try:
        print("Using the line cores:", galah_sp_part1.line_cores)
        # Not sure what line cores are, but basically if it's near the value of the wavelength, and the flux is low we're setting
        # it to 0.
        for line_core_loop in range(0, len(galah_sp_part1.line_cores)):
            line_core_wavelength_index = np.where(
                np.logical_and(
                    abs(ccd_wavelengths_inside_segmask_array - galah_sp_part1.line_cores[line_core_loop]) < 4,
                    ccd_flux_inside_segmask_array < core_minimum))
            flagged_fluxes[line_core_wavelength_index] = 0
    except AttributeError:
        print("No line_cores")
    return flagged_fluxes


# Sums the flagged flux creation functions in one to produce a final product. make_mob in idl
def create_observational_definitions(ccd_wavelengths_inside_segmask_array, ccd_flux_inside_segmask_array, segment_begin_end):
    # Open masks is called in flag functions so no need to do it here for input.
    flagged_fluxes = flag_flux_peaks(ccd_wavelengths_inside_segmask_array, ccd_flux_inside_segmask_array, segment_begin_end)
    flagged_fluxes = flag_continuum(ccd_wavelengths_inside_segmask_array, ccd_flux_inside_segmask_array, flagged_fluxes)
    flagged_fluxes = cutting_low_flux(ccd_wavelengths_inside_segmask_array, ccd_flux_inside_segmask_array, segment_begin_end, flagged_fluxes)
    flagged_fluxes = removing_nlte(ccd_wavelengths_inside_segmask_array, ccd_flux_inside_segmask_array, flagged_fluxes)
    return flagged_fluxes


# Compare atomic_line array to segments we have to see which fit in. But we already indexed atomic array against segments
# in galahpart4 so this might be useless@@@@
def atomic_lines_in_segments(desired_atomic_lines_index, segment_begin_end, line_atomic, depth):
    number_of_segments = len(segment_begin_end)
    # The list to add the indexes that are in the segments and linemask. Overall it seems to inclujde almost everything.
    # currently we take every single value out of the 24k - I guess because we already detailed the important values
    # in rdlin when making the index to create line_atomic. @@@@@@@@@@@@@@@@@@@@
    # Selecting lines within wavelength segments (We've done that SO MUCH!!) and deep than a given depth
    # this is de in idl, I think it's just a buffer?
    buffer = 0.7 #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    # Is this even needed? We index lineatomic in smerdlin using the index of important wavelengths in galahpart4 @@
    # So we already know they're in the segments. @@@@
    for segment in range(0, number_of_segments):

        # Here we reference the data created in sme_rdlin, which takes about 70 seconds for the full 300k lines.
        # we test to see which parts are in the segments we have in segment_begin_end which is segment mask.
        desired_atomic_lines_index.extend(np.where(np.logical_and(np.logical_and(line_atomic[:, 2] > segment_begin_end[segment, 0] - buffer,
                                                   line_atomic[:, 2] < segment_begin_end[segment, 1] + buffer),
                                                   depth > galah_sp_part1.depthmin))[0])
        # If broad lines are near (but not inside), select those too.
        for broad_line_single in galah_sp_part1.broad_lines:
            # If any line is within 100a of a broadline we'll include it
            if np.any(abs(broad_line_single - segment_begin_end[segment]) < 100):
                # Where does this broad line exist in our atomic line array? rounds automatically in np where for the array.
                # We have duplicates in line_atiomic and (therefore?) the d_a_l_index, do we want to remove those?
                # We add all(?) the lineatomic to the index, but that makes sense as it is preselected to include
                # the indexes of the wavelengths inside our segments. guess it's more useful for if we run through
                # a single wavelength segment like in idl.
                desired_atomic_lines_index.extend(np.where(line_atomic[:, 2] == broad_line_single)[0])

    return desired_atomic_lines_index


# Now we see which are in the line list, which is the list of important wavelengths to look out for. This is if we aren't
# doing atomic+lnes)in)segments
def atomic_lines_in_linelist(desired_atomic_lines_index, line_atomic, depth, linemask_data):
    buffer = 0.7
    # Used to take the certain number of linemask indexes (roughly 20ish ) what does it stand for? @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    nrline = 20 + ((8000 - galah_sp_part3.effective_temperature) / 1.3E3) ** 4

    # We see which parts of our atomic line data lies within the linelist (the important wavelengths) as well as identifying
    # the peak wavelength if we can, and the broad lines. We get a lot (12kish) duplicates here
    # I don't see why we do this, we've already indexed to only include wavelengths in SegM in galah part 4, so when
    # woudl we have data that is out of that and in the linemask? @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    for line in range(0, len(linemask_data['Sim_Wavelength_Peak'])):
        # Ok.. again we find the wavelengths within the ranges, but this time it'sthe linemask ranges, not segments.
        inside_linemask_index = np.where(
            np.logical_and(np.logical_and(line_atomic[:, 2] > linemask_data['Sim_Wavelength_Start'][line] - buffer,
                                          line_atomic[:, 2] < linemask_data['Sim_Wavelength_End'][line] + buffer),
                           depth > galah_sp_part1.depthmin))
        # We reverse it to take the the largest indexes? But that doesn't necessarily mean the strongest lines
        # does it? @@@@@@@@@@@@@@@@@@@@
        inside_linemask_index = np.flip(inside_linemask_index)
        # Take either the last nrline number or all the index? How odd, why is nrline a thing@@@@@@@@@@@@@@@@@@@@@@@@@@
        desired_atomic_lines_index.extend((inside_linemask_index[0:min([nrline, len(inside_linemask_index)])])[0])
        # ;always select the main line if it's present
        peak_index = (np.where(line_atomic[:, 2] == float(linemask_data['Sim_Wavelength_Peak'][line])))[0]
        if peak_index.size != 0:
            desired_atomic_lines_index.extend(peak_index)
        else:
            print("No peak line (", linemask_data['Sim_Wavelength_Peak'][line],
                  ") available in the atomic line list for line",
                  line, "(", linemask_data['Sim_Wavelength_Start'][line], "to",
                  linemask_data['Sim_Wavelength_End'][line], ")")
        # And of course, always select broad lines when close just like before (but that was broad lines in segments?) @@@@@@
        # What's the diff betwee nthe broad lines? And why ARE we using the line list now and adding it to the same
        # list as the segments information? @@@@@@@@@@@@@@@@@@@@@@@@@@
        for broad_line_single in galah_sp_part1.broad_lines:
            # If any line is within 100a of a broadline we'll include it
            if np.any(abs(broad_line_single - float(linemask_data['Sim_Wavelength_Peak'][line])) < 100):
                # Where does this broad line exist in our atomic line array? rounds automatically in np where for the array.
                # We have duplicates in line_atiomic and (therefore?) the d_a_l_index, do we want to remove those?
                # We add all(?) the lineatomic to the index, but that makes sense as it is preselected to include
                # the indexes of the wavelengths inside our segments. guess it's more useful for if we run through
                # a single wavelength segment like in idl.
                desired_atomic_lines_index.extend(np.where(line_atomic[:, 2] == broad_line_single)[0])
    # but why are we adding linemask AND segment mask lines into the same array and do they ever combine
    return desired_atomic_lines_index


# I think this is make_line in idl. We produce an indexed version of the data from indexed_stellar_information that
# takes only the atomic lines that are either in the segments (but they are already?? @@) or in the linemask, depending
# on what we send the variable 'run' to.
def produce_indexed_atomic_data(segment_begin_end, linemask_data, run = 2):
    """ Here's where they run sme rdlin, but we do that at the beginning and only if we don't have the premade data file"""
    line_atomic, lande_mean, depth, data_reference_array, species, j_e_array, lower_level, \
    upper_level, lu_lande = indexed_stellar_information()

    # If run is 2 we want the indexes of wavelengths in segments, else we want them from the linelist ranges.
    desired_atomic_lines_index = []
    if run == 2:
        desired_atomic_lines_index = atomic_lines_in_segments(desired_atomic_lines_index, segment_begin_end, line_atomic, depth)
    else:
        desired_atomic_lines_index = atomic_lines_in_linelist(desired_atomic_lines_index, line_atomic, depth, linemask_data)

    # Is now an array as we're doing adding to it, and we take only the unique values using Unique as there's a lot of
    # overlap between the segment mask and linemask
    desired_atomic_lines_index = np.unique(desired_atomic_lines_index)
    nselect = len(desired_atomic_lines_index)
    # Why do we have 12 THOUSAND duplicates?
    line_atomic = line_atomic[desired_atomic_lines_index]

    # Now we also sort them according to wavelength. Not sure how needed this is. @@@@@@@@@@@
    sort_line_index = np.argsort(line_atomic[:,2])

    # So now we apply these indexes to the information taken from smerdlin. As these indexes were taken useing line atomic
    # then the position of the index shoudl be fine.
    species = species[desired_atomic_lines_index][sort_line_index]
    line_atomic = line_atomic[sort_line_index]
    lande_mean = lande_mean[desired_atomic_lines_index][sort_line_index]
    depth = depth[desired_atomic_lines_index][sort_line_index]
    data_reference_array = data_reference_array[desired_atomic_lines_index][sort_line_index]
    lower_level = lower_level[desired_atomic_lines_index][sort_line_index]
    upper_level = upper_level[desired_atomic_lines_index][sort_line_index]
    j_e_array = j_e_array[desired_atomic_lines_index][sort_line_index]
    lu_lande = lu_lande[desired_atomic_lines_index][sort_line_index]
    print(nselect, "spectral lines are selected within wavelength segments")

    return line_atomic, lande_mean, depth, data_reference_array, species, j_e_array, lower_level, \
    upper_level, lu_lande


# runs all the rest of 'em and creates the sme structure!
"""Run this one!"""
def create_structure():

    # ; set remaining variables, usually fixed?
    # Have I set these elsewhere with better names? @@@@@@@@@@@@@
    # Global correction factor to all van der Waals damping constants. Values of
    # 1.5 to 2.5 are sometimes used for iron.
    global_waals_correction = 1 #gam6
    # Minimum accuracy for linear spectrum interpolation vs. wavelength.
    wavelength_interpolation_accuracy = 0.005 # accwi

    # Minimum accuracy for sme.sint (Specific intensities on an irregular wavelength grid given in sme.wint.)
    # at wavelength grid points in sme.wint. (Irregularly spaced wavelengths for specific intensities in sme.sint.)
    # Values below 10-4 are not meaningful.
    specific_intensity_accuracy_min = 0.005 # accrt

    # If continuum points are not specied in the mask sme.mob (flagged fluxes) (where sme.mob
    # eq 2), points within sme.clim of the maximum in the observation are con-
    # sidered continuum points and sme.mob is updated.
    backup_contiuum = 0.01 # clim

    # Fractional change in sme.chisq (Chi-square weighted by the observed spectrum
    # flux for each iteration below which convergence is assumed.) below which convergence is assumed.
    chi_square_convergence = 0.001 # chirat

    # Number of "equal-area"  angles at which to calculate specic intensity.
    specific_intensity_angles = 7 #nmu

    # Dunno
    obs_type = 3

    #Type of prole used for instrumental broadening. Possible values are gauss,
    # sinc, or table. See Section 3.4.
    broadening_profile = "gauss" #iptype

    # The equal-area midpoints of each equal-area annulus for which specic inten-
    # sities were calculated.  values for Gaussian quadrature are not conducive
    # to subsequent disk integration (???).
    intensity_midpoints = np.flip(np.sqrt(0.5*(2*np.arange(1,specific_intensity_angles+1))/specific_intensity_angles)) #mu
    id = date.today()
    # This is meant to be taken from gala part 4 as it's -1 for normalisation, but not for iteration.
    # Array of strings, specifying which of the global parameters are free pa-
    # rameters. Possible values are 'TEFF', 'GRAV', 'FEH', 'VMIC', 'VMAC',
    # 'VSINI', 'GAM6', and 'VRAD' (if sme.vrad flag eq -1), and 'CSCALE'
    # (if sme.cscale flag eq -1).
    free_global_parameters = '-1' #glob_free
    # Is this still right?
    atmo_pro = 'interp_atmo_grid'
    try:
        atmosphere_grid_file = galah_sp_part1.atmosphere_grid_file
    except AttributeError or NameError:
        # Why are we using 20212? @@@@@@@@@@@@
        atmosphere_grid_file = 'marcs2012.sav'
    # atmo_grid_vers
    atmosphere_grid_version = 4.0

    # This is set in galahsp4 or whatever, but it IS modified in ab, so we're testing for it here.
    # Array of flags (0 or 1), specifying for which of the 99 atomic elements the
    # abundances are free parameters.
    try:
        galah_sp_part1.atomic_abundances_free_parameters #ab_free
    except AttributeError or NameError:
        atomic_abundances_free_parameters = np.zeros(99)

    # Produce the basic blocks of which wavelengths etc are inside our segment masks.
    ccd_wavelengths_inside_segmask_array, ccd_flux_inside_segmask_array, ccd_flux_error_inside_segmask_array, \
    segment_begin_end, wavelength_start_end_index, ccd_resolution_of_segment_array = object_array_setup()

    number_of_segments = len(wavelength_start_end_index)

    continuum_mask_data, linemask_data = open_masks()

    if galah_sp_part1.radial_velocity_flag == 0:
        # vrad. If flag is 0, we have individual radial velocities for each segment, else, it's global and thus not
        # an array.
        radial_velocity = np.zeros(number_of_segments)
    else:
        radial_velocity = 0

    # We want the conituum scale as 1 (2 values for each segment). 0 means one value per segment, negative means global.
    if galah_sp_part1.continuum_scale_flag == 1:
        continuum_scale = np.ones((2, number_of_segments))
    elif galah_sp_part1.continuum_scale_flag == 0:
        continuum_scale = np.ones(number_of_segments)
    else:
        continuum_scale = 1

    # Only normalise the data if that's how this file is called.
    if 'normalise_flag' in locals() and galah_sp_part1.continuum_scale_flag >= 0:
        ccd_wavelengths_inside_segmask_array, ccd_flux_inside_segmask_array, ccd_flux_error_inside_segmask_array  =\
            pre_normalise(ccd_wavelengths_inside_segmask_array, ccd_flux_inside_segmask_array,
                          ccd_flux_error_inside_segmask_array, segment_begin_end)

    # Produce the array that tells us what fluxes to ignore, which are contiuum or peak. make_mob
    flagged_fluxes = create_observational_definitions(ccd_wavelengths_inside_segmask_array,
                                                      ccd_flux_inside_segmask_array, segment_begin_end)

    # make_line, we get the indexed atomic data and all.
    line_atomic, lande_mean, depth, data_reference_array, species, j_e_array, lower_level, \
    upper_level, lu_lande = produce_indexed_atomic_data(segment_begin_end, linemask_data)
    # Array of ags (0 or 1), specifying for which of the spectral lines the gf values
    # are free parameters.
    spectral_lines_free_parameters = np.zeros(len(species)) # gf_free
    # Tau is the optical depth at some wavelength reference for continuum to determine the deepest point of interest.
    # You can plot T against log tau and when x is 0, it is opaque so anything below that we are viewing.
    # Rhox is the other option of column mass (accumulated mass) of the atmosphere. These are both pretty bad as
    # they cause some kind of spike in the abundance of stars at temperatures but people are working on replacing them
    # with a neural network.
    # Something to do with setdefaultvalue but I can't figure out if it does anything else but set to this string.
    # They're meant to only be set if newsme but I can't find them being made elsewhere
    atmosphere_depth = "TAU"
    atmosphere_interpolation = "RHOX"
    # Plane parallel
    atmoshpere_geometry = "PP"
    """; Temporary fix for stagger grid, which has only EITHER tau or rhox depth scale available:
  if strmatch(atmogrid_file, '*stagger-t*') then atmo_depth = 'TAU'
  if strmatch(atmogrid_file, '*stagger-r*') then atmo_interp = 'RHOX'"""
    # Always two for fits files I think.
    short_format = 2

    # finally we set the dictionary! Contains the info for sme!
    sme = {'version':atmosphere_grid_version, 'id': id, 'teff': galah_sp_part3.effective_temperature, 'grav': galah_sp_part3.log_surface_gravity,
           'feh': galah_sp_part3.metalicity, 'vmic': galah_sp_part3.microturbulence_velocity, 'vmac': galah_sp_part3.macrotublence_velocity,
           'vsini': galah_sp_part3.rotational_velocity, 'vrad': radial_velocity, 'vrad_flag': galah_sp_part1.radial_velocity_flag,
           'cscale': continuum_scale, 'cscale_flag': galah_sp_part1.continuum_scale_flag, 'gam6': global_waals_correction,
           'accwi': wavelength_interpolation_accuracy, 'accrt': specific_intensity_accuracy_min, 'clim': backup_contiuum,
           'maxiter': max(galah_sp_part3.iterations), 'chirat': chi_square_convergence, 'nmu': specific_intensity_angles,
           'abund': makestruct_abund.abundances, 'mu': intensity_midpoints, 'atmo': {'source': atmosphere_grid_file,
                                                                                     'method': 'grid', 'depth': atmosphere_depth,
                                                                                     'interp': atmosphere_interpolation,
                                                                                     'geom': atmoshpere_geometry},
           'sob': ccd_flux_inside_segmask_array, 'uob': ccd_flux_error_inside_segmask_array, 'obs_name': galah_sp_part1.obs_name,
           'obs_type': obs_type, 'iptype': broadening_profile, 'glob_free': free_global_parameters, 'ab_free': atomic_abundances_free_parameters,
           'gf_free': spectral_lines_free_parameters, 'species': species, 'atomic': line_atomic, 'lande': lande_mean,
           'depth': depth, 'lineref': data_reference_array, 'line_term_low': lower_level, 'line_term_upper': upper_level,
           'short_format': short_format, 'line_extra': j_e_array, 'line_lulande': lu_lande, 'nseg': number_of_segments,
           'wran': segment_begin_end, 'wave': ccd_wavelengths_inside_segmask_array, 'wind': wavelength_start_end_index,
           'mob': flagged_fluxes, 'ipres': ccd_resolution_of_segment_array, 'auto_alpha': galah_sp_part1.auto_alpha
           }

    # Only do so if nlte is set on (1)
    if galah_sp_part1.nonlte_flag:
        nltestruct = {'nlte_pro': 'sme_nlte', 'nlte_elem_flags': galah_sp_part1.nonlte_element_flags, 'nlte_subgrud_size':
                      [3, 3, 3, 3], 'nlte_grids': galah_sp_part1.nonlte_atmosphere_grids, 'nlte_debug': 1}
        # In idl this is made super weirdly, there's a 'NLTE' in there if it's newsme? Is it a dictionary with create_struct
        # or not? It says they are assigned values but it's so weird. what I've done is the literal translation of
        # the code but it's IDL so who fricking knows.
        sme['NLTE'] = nltestruct

    store_sme_input(sme)


# IDL saves it as an inp file, but..
# Called from make_struct function
def store_sme_input(sme):
    input_file = open(r'OUTPUT/'+galah_sp_part1.obs_name+'_SME.pkl', 'wb')
    pickle.dump(sme,input_file)
    input_file.close()

# Checks there are enough points left in line makss before running
def check_line_points(sme):
    i = np.where(sme['sme'])


def run_sme(run = 0):
    if run:
        r"D:\Work\SME_574D:\Work\SME_574\sme_main.pro"
        pass



# this would be the reduced data array but we have it saved already. To do: write the code that reduces it and calls
# on the reduction.
segment_mask = galah_sp_part1.obs_name + '_Segm.dat' # is this the wrong one..
create_structure()