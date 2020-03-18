"""
We read the segment mask segm, removing the comment lines and columns.
We proceed to open the ccd1-4 files and take their wavelength and resolution out to plot a resolution mask and
interpolating to be able to find resolution at peak wavelengths that we calculate at that point of the segment mask.
We open the 4 ccd fits files to produce a telephonenumber_Sp file that contains Wavelength array, observed flux sign (sob)
uncertainity in flux (uob)(essentially error), smod (not important now) and mod (Idk) that we call in make_struct
to modify the wavelength array by the radial velocity doppler shift and make some sme variables.
"""

from astropy.io import fits
import galah_sp_part1
import galah_sp_part2
import galah_sp_part3
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
# Import the wavelength, flux, and error from the ccd readings
from galah_sp_part2 import total_ccd_wavelength, total_ccd_sob_flux, total_ccd_flux_error_uob, resolution_factor
import pickle

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


def structure_creation(segment_mask):
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
              " create a data file for later use for any star with the ", galah_sp_part1.setup_for_obs, "release")
        # Will create the data file for use if it doesn't exist.
        import sme_rdlin
        # Now we know it exists, we can do what we were trying to before.
        line_atomic, lande_mean, depth, data_reference_array, species, j_e_array, lower_level, \
        upper_level, lu_lande = pickle.load(open("galah_master_file_arrays", "rb"))
        line_atomic, lande_mean, depth, data_reference_array, species, j_e_array, lower_level, \
        upper_level, lu_lande = line_atomic[data_index], lande_mean[data_index], depth[data_index], \
                                data_reference_array[data_index], species[data_index], j_e_array[data_index], \
                                lower_level[data_index], upper_level[data_index], lu_lande[data_index]

    #line_list = fits.open(r'C:\Users\jama2357\Documents\Galafiles\GALAH\LINELIST/%s' % galah_sp_part1.line_list)[1]

    # Takes the last three digits and turns into an int. The last 3 represent the fibre of the ccd used that we want to
    # look at.
    object_pivot = int(str(galah_sp_part1.object_for_obs)[-3:])
    # This is the temporary segment mask file with only one segment inside it. Better way to do this.....
    # Initially the segmask is only the setup name (DR2_segm.dat) which is then changed when we create our own using
    # obsname.
    """segment_mask = galah_sp_part1.obs_name + '_Segm.dat'"""
    segment_mask = galah_sp_part1.setup_for_obs + '_Segm.dat'

    # Setting speed of light
    c = 2.99792458E8

    # What is status. At one point during make_obs, status it says if status = 0 it's an error.
    # It gets set to 1 only at 1 point after selecting wavelength regions inside segments, before
    # starting the pre normalisation procedure.
    status = 0

    # Segm_mask is _Segm.data, unsurprisingly. It takes the start, end wavelength, and the resolution base guess of 35000
    # which is totally wrong. That's something to fix/modify when we're streamlining this.
    # The seperator separates via 2 spaces or a space and ; or . and space. We use the overflow column to account for lines
    # that begin with ; (commented out) which we immediately delete afterwards. Hard to separate comments when ; was forgotten.
    # Can't do lower case potentials because it randomly selects lower case letters in the middle of sentences. Something is wrong.
    # Same as in galah part 4, this should be a dummy file with only one segment, so why the bloody hell do
    # we loop through it so much rater than just using its valjues? @@@@@@@@@@@@@@@@@
    segment_mask_data_with_res = pd.read_csv(
         segment_mask, sep='[ ]{2,}| ;| [A-Z]', header=None,
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
    print(segment_mask_data_with_res)
    ccd1_res_file = fits.open(r'ccd1_piv.fits', ext=0)
    ccd1_res_data = (ccd1_res_file[0].data)
    ccd2_res_file = fits.open(r'ccd2_piv.fits', ext=0)
    ccd2_res_data = (ccd2_res_file[0].data)
    ccd3_res_file = fits.open(r'ccd3_piv.fits', ext=0)
    ccd3_res_data = (ccd3_res_file[0].data)
    ccd4_res_file = fits.open(r'ccd4_piv.fits', ext=0)
    ccd4_res_data = (ccd4_res_file[0].data)

    # Takes the last three digits and turns into an int. The last 3 represent the fibre of the ccd used that we want to
    # look at.
    object_pivot = int(str(galah_sp_part1.object_pivot)[-3:])
    # We're making an array of the resolution of the (actually x) axis
    # Extracting the row of data of ccd1 that matches the piv number (-1 as piv starts at 1)
    # What IS piv why are we interested in the piv valued row of the data? @@@@@@@@@@@@@@@@@
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
    wavelength_res_y_collection = np.concatenate([ccd1_obj_piv_y, ccd2_obj_piv_y, ccd3_obj_piv_y, ccd4_obj_piv_y])
    # iterpolate the res collections to find wlpeak for them later
    interpolation = interp1d(wavelength_res_x_collection, wavelength_res_y_collection)

    # Loop through each wavelength range (consisting of a start and end) in the data and take the peak value and its res factor.
    # The idl code seems to ignore dum1 column, and instead calls seg_st and en
    # We make a list to avoid calling pandas data frame repeatedly, I believe this is faster and avoids any copy errors.
    temporarywavelengthpeaklist = []
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

    # Checks for a negative range which I assume would cause issues. Check it works @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    if min(segment_mask_data_with_res['Wavelength_End'] - pd.to_numeric(
            segment_mask_data_with_res['Wavelength_Start'])) <= 0:
        print("Segment has a negative range!")
        exit()

    # Checks for overlapping segments if there's more than one. Double check this works when you have data @@@@@@@@@@@@@@
    if len(segment_mask_data_with_res['Wavelength_End']) > 1:
        if max(segment_mask_data_with_res['Wavelength_End'][0:len(segment_mask_data_with_res['Wavelength_End'])]
               - pd.to_numeric(
            segment_mask_data_with_res['Wavelength_Start'][1:len(segment_mask_data_with_res['Wavelength_Start'])])) > 0:
            print("Overlapping segments")
            exit()

    # Ok we are opening the 4 ccds and assigning their values to new columns as lists (wl, sob, uob)
    # object_id_desired = input("Gimme the id you want: ")
    object_id_desired = galah_sp_part1.object_for_obs

    # Shifts the wavelengths to account for the doppler effect.
    total_ccd_wavelength_dopplered = (total_ccd_wavelength) / (
                (galah_sp_part3.radial_velocity_global[0] * (1E3 / c)) + 1E0)

    # Skipping the "if not normalising, make sure there are lines in the segment" part as I think we can do that here
    # Unsure how many there will be, so arrays can't be used. Instead we use a temp list.
    ccd_wavelengths_inside_segmask = []
    ccd_flux_inside_segmask = []
    ccd_flux_error_inside_segmask = []
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
        ccd_flux_inside_segmask.extend(total_ccd_sob_flux[wavelength_inside_segmask_index])
        ccd_flux_error_inside_segmask.extend(total_ccd_flux_error_uob[wavelength_inside_segmask_index])
        # Numpy array of indexes of the first and final wavelengths per segment with column 0 being the first.
        # the wind SME wants seems to just be the final values so keep a note of that @@@@@@@@@@@@@@@
        wavelength_start_end_index[segment, 0] = (wavelength_inside_segmask_index[0][0])
        wavelength_start_end_index[segment, 1] = (wavelength_inside_segmask_index[-1][-1])

    # Turning lists into arrays for numpy indexing with np.where. Will this slow things down too much by having both list
    # and array? Delete one? Or does the creation take too long? no. probs not.
    ccd_wavelengths_inside_segmask_array = np.array(ccd_wavelengths_inside_segmask)
    ccd_flux_inside_segmask_array = np.array(ccd_flux_inside_segmask)
    ccd_flux_error_inside_segmask_array = np.array(ccd_flux_error_inside_segmask)

    # Number of segments that contain visible spectra. Wind (wavelength_start_end_index)
    # is the array of indexes of first and final wavelengths per segment.
    number_of_segments = len(wavelength_start_end_index)
    if number_of_segments == 0:
        print("No observations in segment mask")
        exit()

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
    segment_begin_end = np.copy(
        wavelength_start_end_index)  # copy to avoid overwriting wind. Just using it for its size.
    for windex_row in range(0, number_of_segments):
        # At indexes 0,0 and 0,1 (and then 1,0etc) of the index array 'wavelength_start_end_index' we take the value and
        # apply it to the wavelngtharray as the values we have taken are indexes of the first and last wavelength of each
        # segment. windrow, 0 is the segment beginning. , 1 is the end.
        segment_begin_end[windex_row, 0] = total_ccd_wavelength_dopplered[
            int(wavelength_start_end_index[windex_row, 0])]
        segment_begin_end[windex_row, 1] = total_ccd_wavelength_dopplered[
            int(wavelength_start_end_index[windex_row, 1])]

    # Pre-normalization steps. ";Performs first-guess normalisation by robustly converging straight line fit to high pixels"
    # I removed the requirement for sob > 0. @@@@@@@@@@
    for segment_band in range(0, number_of_segments):
        # Finds the index where the wavelengths are between the start and end of each segment, to be able to loop each seg as i
        # We do this in wavelength_inside_segmask_index already, but we know that the IDL code is so damn convoluted
        # then maybe a case occurs where we didn't make that. Doesn't make any sense to me, but is possible.
        # If this is slow, we can optimise somehow probably. @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        segment_indexes = (np.where(np.logical_and(ccd_wavelengths_inside_segmask >= segment_begin_end[segment_band, 0],
                                                   ccd_wavelengths_inside_segmask <= segment_begin_end[
                                                       segment_band, 1])))[0]

        # If count is greater than 20, do this. len segindex is the number of values that fit our criteria.
        if len(segment_indexes) > 20:

            # Take the coefficients of the polyfit, and the flux of it too using the equation from IDL and then applies
            # a ploynomial fit, removing outlying values until it no longer changes. The outlying value limit is hardcoded
            # Wew then use the equation fit to it to normalise it
            polyfit_coefficients, continuous_function, sigma_flux_comparison, wavelength_array = autonormalisation(
                ccd_wavelengths_inside_segmask_array[segment_indexes],
                ccd_flux_inside_segmask_array[segment_indexes], 1, 0)
            # The index of the wavelengths we want to keep as they are not outliers.
            inlier_index = np.isin(ccd_wavelengths_inside_segmask_array[segment_indexes], wavelength_array)
            # Puts it to a relative flux out of 1 that we see in the abundance charts.
            # Have to take relative indexes for cont_func to have the correct shape. Is the wavelengths outside of it
            # not being modified an issue? Should we remove the other wavelengtghs? @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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

    # end of prenormalisation! Phew!

    # of opening it we're just gonna import it from sp4 as an array
    line_mask = galah_sp_part1.setup_for_obs + "_Sp.dat"  # Line masks for log(g) sensitive lines, i.e Fe Ti Sc

    # Reads out the columns which are centre/peak wavelength, start and end of wavelength peak (all simulated), and atomic number
    # Careful with the python enginer, it's slower. If we are looking at BIG data files this might be bad.
    # line0, line_st, and line_en in idl. they ignore atomic number for some reason @@@@@@@@@@@@@@@@@@@
    linemask_data = pd.read_csv(
        line_mask, sep="[ ]{2,}", header=None, engine='python',
        names=["Sim_Wavelength_Peak", "Sim_Wavelength_Start", "Sim_Wavelength_End", "Atomic_Number"])
    # Removes all rows that begin with ; as it's commented out. Required to be able to modify the data by column rather than
    # by row
    linemask_data = linemask_data[~linemask_data['Sim_Wavelength_Peak'].str.startswith(";")]
    # Reset the index to account for the now missing values
    linemask_data = linemask_data.reset_index(drop=True)

    continuum_mask = str(galah_sp_part1.setup_for_obs) + "_Cont.dat"
    continuum_mask_data = pd.read_csv(continuum_mask, sep="[ ]{1,}", header=None, engine='python',
                                      names=["Continuum_Start", "Continuum_End"])

    """There doesn't seem to be any need to repeat Signal to noise each loop. Especially if we want to recreate it again
    in a second if statement."""

    """this stuff might work or it might not?! CHECK IT SOMEHOW! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ vvvvv"""
    # mob in idl. Mask of observed pixels
    flagged_fluxes = np.zeros(len(ccd_wavelengths_inside_segmask_array))
    x2 = 0
    # for line_loop in range(0,len(linemask_data['Sim_Wavelength_Start'])):
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
                                                         total_ccd_flux_error_uob[max(0, flux_row - 4):min(flux_row + 4,
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

    """; cont mob  - continuum points selected between 0 and 1.2 or 1+3*sigma, where there are no line masks
    ;             avoid buffer zone at edges
    """
    print("Before cont")
    # If the flag is above 0 we're resetting signal to noise. Seems silly to reset it.. Maybe skip for fluxrow in the
    # function above? @@@@@@@@@@@@@@@@@@@@@@@@@@@
    # Basically same as above but for finding contiuum fluxes instead of peak fluxes.
    if galah_sp_part1.continuum_scale_flag >= 0:

        for continuum_loop in range(0, len(continuum_mask_data['Continuum_Start'])):
            signal_to_noise = []
            for flux_row in range(0, len(ccd_flux_inside_segmask_array)):
                signal_to_noise.append(max([1 + 10 / np.mean(ccd_flux_inside_segmask_array[
                                                             max(0, flux_row - 4):min(flux_row + 4, len(
                                                                 ccd_flux_inside_segmask_array))] /
                                                             total_ccd_flux_error_uob[
                                                             max(0, flux_row - 4):min(flux_row + 4, len(
                                                                 ccd_flux_inside_segmask_array))]), 1.5]))
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
        for segment_band in range(0, number_of_segments):

            ";Deselect 70% lowest continuum points in each segment using synthesis. Ensure both ends have continuum points"
            # The fraction of points we want to remove. Potentially modified in the while loop hence kept outside of it.
            fraction = 0.7
            # print("new segment")
            # While the fraction to use is not 0 we continue looping. Python is setting it to E-17 in loop so this fixes it.
            while fraction > 0.01:
                # We take the wavelength indexes of the wavelengths that exist in the segments as always. We just repeat it
                # each time to be able to work on each segment individually I suppose? @@@@@@@@@
                wavelength_inside_segment_index = np.where(
                    np.logical_and(ccd_wavelengths_inside_segmask_array >= segment_begin_end[segment_band, 0],
                                   segment_begin_end[segment_band, 1] >= ccd_wavelengths_inside_segmask_array))
                # The value of the fraction of the fluxes we chose, so how many IS 70% for example.
                value_of_fraction = int(len(ccd_flux_inside_segmask_array[wavelength_inside_segment_index]) * fraction)
                # Takes the index of the 70%th lowest value (our cut off point) from a sorted list. Sorting takes a long
                # time but this if statement doesn't appear to be the bottleneck.
                cutting_flux_value = sorted(ccd_flux_inside_segmask_array[wavelength_inside_segment_index])[
                    value_of_fraction]
                # print((ccd_flux_inside_segmask_array[wavelength_inside_segment_index]))
                # print("cfv ", cutting_flux_value)
                # We take a list of indexes where the flux in the segment is below our cut off point.
                cutting_flux_index = np.where(
                    ccd_flux_inside_segmask_array[wavelength_inside_segment_index] < cutting_flux_value)
                # print("cfi", cutting_flux_index)
                # We need to count how many values there are at the extreme ends of the segment, as we need a continuum on both
                # Here we see how many are in the first 1/3, and how many in the last 2/3.
                # Also makes sure that the flux they have is flagged as continuum
                # currently having issues where some can never find enough points, as there simply AREN'T 10 points.
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

    """ checking the number of 0 to non 0 points.
    print("len", len(flagged_fluxes))
    print("max", max(flagged_fluxes))
    x0 = 0
    x1 = 0
    for x in flagged_fluxes:
        if x > 0:
            x1 +=1

        else:
            x0 +=1
    print("x0, x1", x0, x1)"""

    """this stuff might work or it might not?! CHECK IT SOMEHOW! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@^^^^^^^^^^^^"""

    # ;Avoid strong NLTE-dominated cores in mask, what. @@@@@@@@@@@@@
    core_minimum = 0.6
    if galah_sp_part3.metalicity < -2:
        core_minimum = 0.72
    elif galah_sp_part3.metalicity < -1:
        core_minimum = 0.65

    # Checks for line_cores existing from galahsp1, if no variable exists then we can't do this can we.
    try:
        galah_sp_part1.line_cores
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

    # temporary, should be in part 4
    run = 2
    newsme = True
    """ Here's where they run sme rdlin, but we do that at the beginning and only if we don't have the premade data file"""

    # Selecting lines within wavelength segments (We've done that SO MUCH!!) and deep than a given depth
    # de = 0.7
    nrline = 20 + ((8000 - galah_sp_part3.effective_temperature) / 1.3E3) ** 4
    desired_atomic_lines_index = []
    if run == 2:
        # Wh yare we looping therough the segments AGAIN when we're already looping through them in galahpart 4 when calling
        # make struct @@@@@@@@@@@@@@     Should just be 1???
        for i in range(0, number_of_segments):
            print("wwww")


            # Here we reference the data created in sme_rdlin, which takes about 70 seconds for the full 300k lines.
            # However, it creates a file with all the data so we test to see if it exists, and if it does we just use
            # that and take the index from the beginning that is important lines that we care about.
            # Haven't we done this shit so many times at this point?
            # Has a buffer of 0.7 on either side, for some reason.
            desired_atomic_lines_index.extend(np.where(np.logical_and(np.logical_and(line_atomic[:, 2] > segment_begin_end[i, 0] - 0.7,
                                                       line_atomic[:, 2] < segment_begin_end[i, 1] + 0.7),
                                                       depth > galah_sp_part1.depthmin))[0])
            # If broad lines are near (but not inside), select those too.
            for broad_line_single in galah_sp_part1.broad_lines:
                # If any line is within 100a of a broadline we'll include it
                if np.all(abs((broad_line_single - segment_begin_end[i])) < 100):
                    print("big tru", abs(broad_line_single - segment_begin_end[i]))

                    # Where does this broad line exist in our atomic line array? rounds automatically in np where for the array.
                    desired_atomic_lines_index.extend(np.where(line_atomic[:, 2] == broad_line_single)[0])


    print(len(set(desired_atomic_lines_index)))
    print(len(line_atomic))

# this would be the reduced data array but we have it saved already. To do: write the code that reduces it and calls
# on the reduction.
segment_mask = galah_sp_part1.obs_name + '_Segm.dat'
structure_creation(segment_mask)