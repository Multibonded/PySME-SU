

"""read in spectra nad possibly perform skyline and tellurics correction"""

"Read in 4 band spectra (ext. 0)"

# This assumes spectrum is not a multivisit/stacked spectrum ???@@@@@@@@@@@@@@@@@@
# and what does it stand for @@@@@@@@@@@@@
com = 'com'
import galah_sp_part1
import numpy as np
from astropy.io import fits
from scipy.interpolate import interp1d

# The digit represents ... something. ask. @@@@@@@@@@@@@@@@@@@@@@@@@
if str(galah_sp_part1.object_for_obs)[11] == 1: com = 'com2'

# defines a function that we can use to set cluster if its keyword pre-exists. This kind of thing is a little wrong
# until we find out how we'll be calling these things, because idl works difo.
# What is the cluster? Just the group of stars? @@@@@@@@
def set_cluster(cluster):

    if cluster == 'avatar':
        # Obviously need to change these locations to something realistic. @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        spectra = r'/avatar/buder/trunk/GALAH/SPECTRA/'

    elif cluster == 'gemini12':
        spectra = '/shared-storage/buder/svn-repos/trunk/GALAH/SPECTRA'
    else:
        #spectra = r'SPECTRA/'
        #temporary use for my computer
        spectra = r'C:\Users\jama2357\Documents\Galafiles\GALAH\SPECTRA'


'the next part takes six numbers from the object to find the right folder. Another case of changing after we finish ' \
'testing. God it is so damn frustrating to not have a working example directory but this crap'

# We open the four ccd files as we do in the makestruct_createdatfiles.


# The location of these files will CHANGE! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# We use primary here for flux, and [1] (1st extention) for serr/uob. other ones are for more unimportant values.
# calls for flux from f1, 0, h1 and serr from f1, 1, h1. (printing to a list called h1 basically.) We change serr to uob
# later. The wavefits just makes the wavelength list.
#ccd1_object_line_mask = fits.open(r"C:\Users\jama2357\Documents\Galafiles\GALAH\SPECTRA\test_lbol\1504270008010491.fits")
ccd1_object_line_mask = fits.open(r"1504270008010491.fits")

# Creating the wavelength array, the flux, and the relative error. The + * arrange just plots the wavelength array.
ccd1_wavelength_x = ccd1_object_line_mask[0].header['CRVAL1'] +\
                    (ccd1_object_line_mask[0].header['CDELT1'] * np.arange(ccd1_object_line_mask[0].header['NAXIS1']))
ccd1_sob_flux_y = ccd1_object_line_mask[0].data
ccd1_relative_flux_error_y = ccd1_object_line_mask[1].data

#ccd2_object_line_mask = fits.open(r"C:\Users\jama2357\Documents\Galafiles\GALAH\SPECTRA\test_lbol\1504270008010492.fits")
ccd2_object_line_mask = fits.open(r"1504270008010492.fits")

ccd2_wavelength_x = ccd2_object_line_mask[0].header['CRVAL1'] + \
                    (ccd2_object_line_mask[0].header['CDELT1'] * np.arange(ccd2_object_line_mask[0].header['NAXIS1']))
ccd2_sob_flux_y = ccd2_object_line_mask[0].data
ccd2_relative_flux_error_y = ccd2_object_line_mask[1].data

#ccd3_object_line_mask = fits.open(r"C:\Users\jama2357\Documents\Galafiles\GALAH\SPECTRA\test_lbol\1504270008010493.fits")
ccd3_object_line_mask = fits.open(r"1504270008010493.fits")

ccd3_wavelength_x = ccd3_object_line_mask[0].header['CRVAL1'] + \
                    (ccd3_object_line_mask[0].header['CDELT1'] * np.arange(ccd3_object_line_mask[0].header['NAXIS1']))
ccd3_sob_flux_y = ccd3_object_line_mask[0].data
ccd3_relative_flux_error_y = ccd3_object_line_mask[1].data

#ccd4_object_line_mask = fits.open(r"C:\Users\jama2357\Documents\Galafiles\GALAH\SPECTRA\test_lbol\1504270008010494.fits")
ccd4_object_line_mask = fits.open(r"1504270008010494.fits")

ccd4_wavelength_x = ccd4_object_line_mask[0].header['CRVAL1'] + \
                    (ccd4_object_line_mask[0].header['CDELT1'] * np.arange(ccd4_object_line_mask[0].header['NAXIS1']))
ccd4_sob_flux_y = ccd4_object_line_mask[0].data
ccd4_relative_flux_error_y = ccd4_object_line_mask[1].data

# We combine the wavelength, flux, and relative flux error into one array each. Represented by dum1, etc in makestruct
total_ccd_wavelength = np.concatenate(
    (ccd1_wavelength_x, ccd2_wavelength_x, ccd3_wavelength_x, ccd4_wavelength_x), axis=0)
total_ccd_sob_flux = np.concatenate(
    (ccd1_sob_flux_y, ccd2_sob_flux_y, ccd3_sob_flux_y, ccd4_sob_flux_y), axis=0)
total_ccd_relative_flux_error = np.concatenate(
    (ccd1_relative_flux_error_y, ccd2_relative_flux_error_y, ccd3_relative_flux_error_y, ccd4_relative_flux_error_y), axis=0)
# Do we need to check for 0 flux? @@@@@@@@@@@@@@@@@


# It seems we're now determining the centre of the peaks and then applying the resolution factor (which is a stated variable)
# Idl just uses fxpar so I am not sure which file it's using to determine slitmask resolution? Are they all the same @@@@@@@
if ccd1_object_line_mask[0].header['SLITMASK'] == 'IN      ': # High res
    resolution_factor = 1.789
else: # low
    resolution_factor = 1.0

"resolution time here. Same as from makestruct at the beginning, getting the concatented resolutions @@@@@@@@@@@"
#mrdfits reads fits file of ccd1_piv etc. and saves the data as res1. The piv files are essentially the resolution files
# I believe
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
object_pivot=int(str(galah_sp_part1.object_pivot)[-3:])
# We're making an array of the resolution of the (actually x) axis
# Extracting the row of data of ccd1 that matches the piv number (-1 as piv starts at 1)
# What IS piv why are we interested in the piv valued row of the data? @@@@@@@@@@@@@@@@@
ccd1_obj_piv_y = (ccd1_res_data[object_pivot-1])
# Creates a wavelength list from starting CRVAL1 (4700) in steps of CRDELT1 (0.1) until it matches the array len of NAXIS1
ccd1_wavelength_x = ccd1_res_file[0].header['CRVAL1'] + (ccd1_res_file[0].header['CDELT1'] * np.arange(ccd1_res_file[0].header['NAXIS1']))
#                               start                        steps                            number of data points

# Remember these are all arrays.
ccd2_obj_piv_y = (ccd2_res_data[object_pivot-1])
ccd2_wavelength_x = ccd2_res_file[0].header['CRVAL1'] + (ccd2_res_file[0].header['CDELT1'] * np.arange(ccd2_res_file[0].header['NAXIS1']))
ccd3_obj_piv_y = (ccd3_res_data[object_pivot-1])
ccd3_wavelength_x = ccd3_res_file[0].header['CRVAL1'] + (ccd3_res_file[0].header['CDELT1'] * np.arange(ccd3_res_file[0].header['NAXIS1']))
ccd4_obj_piv_y = (ccd4_res_data[object_pivot-1])
ccd4_wavelength_x = ccd4_res_file[0].header['CRVAL1'] + (ccd4_res_file[0].header['CDELT1'] * np.arange(ccd4_res_file[0].header['NAXIS1']))

# combines the arrays together to interpolate them
wavelength_res_x_collection = np.concatenate([ccd1_wavelength_x, ccd2_wavelength_x, ccd3_wavelength_x, ccd4_wavelength_x])
wavelength_res_y_collection = np.concatenate([ccd1_obj_piv_y, ccd2_obj_piv_y, ccd3_obj_piv_y, ccd4_obj_piv_y])
# iterpolate the res collections to find wlpeak for them later
interpolation = interp1d(wavelength_res_x_collection, wavelength_res_y_collection)

" Read in reduction pipeline output"
# What the heck is iraf_dr53
# Why is galah_sp calling on BailerJones k2seis fits file? What does it have that we don't in cannon? Just extra stuff I guess
# Must changed iraf and stuff variable names once I find out what they stand for.
reduction_and_analysis=fits.open(r'sobject_iraf_53'
                                 r'_2MASS_GaiaDR2_WISE_PanSTARRSDR1_BailerJones_K2seis_small.fits')
# Finds the index of where we can find our object id.
reduction_and_analysis_index = np.where(reduction_and_analysis[1].data['sobject_id']==galah_sp_part1.object_for_obs)
reduction_and_analysis_data = reduction_and_analysis[1].data[reduction_and_analysis_index]
# Checks the index exists therefore the object id exists in the place we're looking for.
if not reduction_and_analysis_index[0].size:
    print("Does not exist")
    exit()
# Check this is the right name @@@@@@@@@@@@@@
vbary = reduction_and_analysis_data['v_bary']
velocity_barycenter = vbary
#print(repr(reduction_and_analysis[1].header))

# We're interpolating the telluric data on to the grid dimensions of the spectra we read.
"Telluric correction (Incr errors), removest he wavelengths that the earths atmosphere produces"

telluric = fits.open(r'telluric_noao_21k.fits')
# Taking the wavelengths from the telluric fits file
telluric_wavelengths = telluric[1].data['wave']/(1-(velocity_barycenter/299792.458))
# Made so we can append a new flux to it. Change it if we don't need the two if loops directly below. (max < max)
telluric_flux = telluric[1].data['flux']
# Is interp1d the correct choice?
# The flux is already normalised at 1 or below I think
# Also need to extend the edges of telluric if the current range isn't enough to fit the wavelength!
# scipy interp extends automatically, but numpy is faster. Gotta pick one eventually though.
# ah shooot is doing this taking longer than just using scipy though?
if max(telluric_wavelengths) < max(total_ccd_wavelength):
    # Adds a new max wavelength to telluric at the end, and sets the flux of it to 1.
    telluric_wavelengths = np.append(telluric_wavelengths, max(total_ccd_wavelength))
    telluric_flux = np.append(telluric_flux, 1)

if min(telluric_wavelengths) > min(total_ccd_wavelength):
    # Adds a minimum wavelength at the beginning of the array, and sets flux to 1.
    telluric_wavelengths = np.insert(telluric_wavelengths, 0, min(total_ccd_wavelength))
    telluric_flux = np.insert(telluric_flux,0, 1)
# Finally make the interpolation function coefficient thing
telluric_interpolation = interp1d(telluric_wavelengths, telluric_flux)

# Issue with telluric correction in that the telluric wavelengths start at 4716A, but our wavelengths start at 4713..
# So can't interpolate it. But why aren't we just like.. removing the same wavelegnths? Probably best to double check
# how it works again @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
telluric_interpolated_ccd_wavelength_array = telluric_interpolation(total_ccd_wavelength)
telluric_below_zero = np.where(telluric_interpolated_ccd_wavelength_array < 0.81)
telluric_above_one = np.where(telluric_interpolated_ccd_wavelength_array > 0.998)
# Prevents the equation below dividing by 0 and stuff. Should make this better @@@@@@@@@@@@@@@@@@@@@@@@@
telluric_interpolated_ccd_wavelength_array[telluric_below_zero] = 0.81
telluric_interpolated_ccd_wavelength_array[telluric_above_one] = 1.0
# Increasing the error due to this, not very scientific.
total_ccd_relative_flux_error = total_ccd_relative_flux_error/(telluric_interpolated_ccd_wavelength_array*5 -4 )



"Skyline correction (incr. errors)"
# Again need to have a case for if the wavelengths we input are out of the range of the sky mask @@@@@@@@@@@@@@@@@@@@@@@
sky_mask = fits.open(r'Skyspectrum_161105.fits')
# Adjusting for... sometrhing @@@@@@@@@@@@
sky_mask_wavelengths = ((sky_mask[1].data['wave'])/(1-(velocity_barycenter/299792.458)))
sky_mask_interpolation = interp1d(sky_mask_wavelengths, sky_mask[1].data['sky'])
# 'sky' is... a bunch of zeroes for object 150427000801049?? @@@@@@@@@@@@@@@@@@@
sky_mask_interpolated_wavelengths = sky_mask_interpolation(total_ccd_wavelength)
total_ccd_relative_flux_error = total_ccd_relative_flux_error+sky_mask_interpolated_wavelengths


"Final uncertainties UOB"
# Relative error to actual error.
total_ccd_flux_error_uob = total_ccd_relative_flux_error * total_ccd_sob_flux
"There's a part here about making sure sob is finite, is this important? Does our code work when any values are infinite?" \
"Should check for int overflow maybe? @@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
# Saves the arrays as columns seperated by two spaces (delimeter)
np.savetxt(galah_sp_part1.obs_file, np.c_[total_ccd_wavelength,total_ccd_sob_flux,total_ccd_flux_error_uob], delimiter='  ')
print(total_ccd_wavelength,total_ccd_sob_flux,total_ccd_flux_error_uob)
# The smooth(sob, 10, /NaN) in the idl code boxcar smoothing to fabricate a new point for it
# we do the same here with convolve. We don't need it. ignore it.







