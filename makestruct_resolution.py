"""
We read the segment mask segm, removing the comment lines and columns.
We proceed to open the ccd1-4 files and take their wavelength and resolution out to plot a resolution mask and
interpolating to be able to find resolution at peak wavelengths that we calculate at that point of the segment mask.
We open the 4 ccd fits files to produce a telephonenumber_Sp file that contains Wavelength array, observed flux sign (sob)
uncertainity in flux (uob)(essentially error), smod (not important now) and mod (Idk) that we call in make_struct
to modify the wavelength array by the radial velocity doppler shift and make some sme variables.
"""

import pandas as pd

from astropy.io import fits

import numpy as np

from scipy.interpolate import interp1d
resolution = 0 #Low
resolution = 1 #High resolution

# Takes a 3 digit number instead of looping back to the beginning. Str to subscript it later.
#objectinput=str(input("Object: "))
object_input = "150427000801049"
# Takes the last three digits and turns into an int. The last 3 represent the fibre of the ccd used that we want to
# look at.
object_pivot=int(object_input[-3:])

# Setting speed of light
c = 2.99792458E8

# What is status. At one point during make_obs, status it says if status = 0 it's an error.
# It gets set to 1 only at 1 point after selecting wavelength regions inside segments, before
# starting the pre normalisation procedure.
status = 0

# Next it reads the 'SPECTRA/'+obs_file,dum1,dum2,dum3,dum5,/silent,format='d,d,d,d'
# So (readcol) it reads the obs_file and outputs the column as dum1 etc?
# So dum 1 = sme.wave, dum2 = sme.sob, dum 3 = sme.uob, dum 4 = sme.smod, dum 5 = sme.mob
# But DR2_SP only has 4 columns....... perhaps number 4 doesn't exist, then

# Reads out the columns which are centre/peak wavelength, start and end of wavelength peak (all simulated), and atomic number
# Careful with the python enginer, it's slower. If we are looking at BIG data files this might be bad.
linemask_sp=pd.read_csv(
    r"C:\Users\jama2357\Documents\Galafiles\GALAH\DATA\DR2_Sp.dat",sep="[ ]{2,}",header=None,engine='python',
    names=["Sim_Wavelength_Peak","Sim_Wavelength_Start","Sim_Wavelength_End", "Atomic_Number"]
)
# Removes all rows that begin with ; as it's commented out. Required to be able to modify the data by column rather than
# by row
linemask_sp=linemask_sp[~linemask_sp['Sim_Wavelength_Peak'].str.startswith(";")]
# Reset the index to account for the now missing values
linemask_sp = linemask_sp.reset_index(drop=True)


# We're reading the segment mask created in galah sp
# Next line is ;Read segments and sort in increasing wavelength. [But it Doesn't seem to sort it at all? It's already ordered]
# readcol,'DATA/'+segm_mask,seg_st,seg_en,seg_ipres,comment=';',/silent,format='d,d,f'
# Segm_mask is _Segm.data, unsurprisingly. It takes the start, end wavelength, and the resolution base guess of 35000
# which is totally wrong. That's something to fix/modify when we're streamlining this.
# The seperator separates via 2 spaces or a space and ; or . and space. We use the overflow column to account for lines
# that begin with ; (commented out) which we immediately delete afterwards. Hard to separate comments when ; was forgotten.
# Can't do lower case potentials because it randomly selects lower case letters in the middle of sentences. Something is wrong.
segmask_segm=pd.read_csv(
    r"C:\Users\jama2357\Documents\Galafiles\GALAH\DATA\DR2_Segm.dat",sep='[ ]{2,}| ;| [A-Z]',header=None,
    names=["Wavelength_Start","Wavelength_End","Resolution", "comment", "overflow"],engine='python', skipinitialspace=True,
    usecols=["Wavelength_Start","Wavelength_End","Resolution"]) # Here we auto delete the comment and overflow columns

# ~ asks for negation, removes any row that starts with ; in the first column. Unsure what the str does but required.
segmask_segm=segmask_segm[~segmask_segm['Wavelength_Start'].str.startswith(";")]
# Reset the index to account for the now missing values
segmask_segm = segmask_segm.reset_index(drop=True)

# Next is: ;This sorts out segments 35 \n i=sort(seg_st) & seg_st=seg_st[i] & seg_en=seg_en[i] & seg_ipres=seg_ipres[i]
# Sorts in ascending order of wavelength of starting wavelength, and re-orders the others.
# But it's already ordered? weird. Will leave it in case data comes in unordered
segmask_segm.sort_values(by=['Wavelength_Start', 'Wavelength_End'])


# Next is: ;Adjust resolution to resolution maps
# fix converts floats to int wtf why
# strmid(x,y,z) extracts a string from string x at position y for z characters long. It comes from galah sp, and is
# a command line argument takeaway. At one point strmid(object,11,2)=2 is important. So is it 2?
# It seems to start 3 characters before the end, but then takes the full string length?
"""piv=fix(strmid(object_pivot,strlen(object_pivot)-3,strlen(object_pivot)))
res1=mrdfits('DATA/ccd1_piv.fits',0,res_h1,/silent)
res2=mrdfits('DATA/ccd2_piv.fits',0,res_h2,/silent)
res3=mrdfits('DATA/ccd3_piv.fits',0,res_h3,/silent)
res4=mrdfits('DATA/ccd4_piv.fits',0,res_h4,/silent)
resy1=res1[0:fxpar(res_h1,'NAXIS1')-1,piv-1]
resx1=fxpar(res_h1,'CRVAL1')+fxpar(res_h1,'CDELT1')*indgen(fxpar(res_h1,'NAXIS1'))
resy2=res2[0:fxpar(res_h2,'NAXIS1')-1,piv-1]
resx2=fxpar(res_h2,'CRVAL1')+fxpar(res_h2,'CDELT1')*indgen(fxpar(res_h2,'NAXIS1'))
resy3=res3[0:fxpar(res_h3,'NAXIS1')-1,piv-1]
resx3=fxpar(res_h3,'CRVAL1')+fxpar(res_h3,'CDELT1')*indgen(fxpar(res_h3,'NAXIS1'))
resy4=res4[0:fxpar(res_h4,'NAXIS1')-1,piv-1]
resx4=fxpar(res_h4,'CRVAL1')+fxpar(res_h4,'CDELT1')*indgen(fxpar(res_h4,'NAXIS1'))
resx=[resx1,resx2,resx3,resx4]
resy=[resy1,resy2,resy3,resy4]  """



#mrdfits reads fits file of ccd1_piv etc. and saves the data as res1. The piv files are essentially the resolution files
# I believe
ccd1_res_file = fits.open(r'C:\Users\jama2357\Documents\Galafiles\GALAH\DATA\ccd1_piv.fits', ext=0)
ccd1_res_data = (ccd1_res_file[0].data)
ccd2_res_file = fits.open(r'C:\Users\jama2357\Documents\Galafiles\GALAH\DATA\ccd2_piv.fits', ext=0)
ccd2_res_data = (ccd2_res_file[0].data)
ccd3_res_file = fits.open(r'C:\Users\jama2357\Documents\Galafiles\GALAH\DATA\ccd3_piv.fits', ext=0)
ccd3_res_data = (ccd3_res_file[0].data)
ccd4_res_file = fits.open(r'C:\Users\jama2357\Documents\Galafiles\GALAH\DATA\ccd4_piv.fits', ext=0)
ccd4_res_data = (ccd4_res_file[0].data)

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
wavelength_piv_y_collection = np.concatenate([ccd1_obj_piv_y, ccd2_obj_piv_y, ccd3_obj_piv_y, ccd4_obj_piv_y])
# iterpolate the res collections to find wlpeak for them later
interpolation = interp1d(wavelength_res_x_collection, wavelength_piv_y_collection)

""";This part is to get the resolution from the resolution map
for i=0,n_elements(seg_ipres)-1 do begin
    seg_mid=0.5*(seg_st[i]+seg_en[i])
    seg_ipres[i]=interpol(resy,resx,seg_mid)*resolution_factor # Does resx seem to have no effect? x seems to be in interpol
endfor
print,'Resolution: ',seg_ipres

if min(seg_en-seg_st) le 0 then begin
    print,'Segment has negative range'
    return
endif

if n_elements(seg_en) gt 1 then begin
    if max(seg_en[0:*]-seg_st[1:*]) gt 0 then begin
        print,'Overlapping segments'
        return
    endif
endif"""
# It seems we're now determining the centre of the peaks and then applying the resolution factor (which is a stated variable)
if resolution == 1: # High res
    resolution_factor = 1.789
elif resolution == 0: # low
    resolution_factor = 1.0
# Loop through each wavelength range (consisting of a start and end) in the data and take the peak value and its res factor.
# The idl code seems to ignore dum1 column, and instead calls seg_st and en


# We make a list to avoid calling pandas data frame repeatedly, I believe this is faster and avoids any copy errors.
# this loop is in part 4 n ow? :@@@@@@@@@@@@@@@@@@@@@@
temporarywavelengthpeaklist = []
for wavelengthband in range(0, len(segmask_segm['Resolution'])):
    # Calculating the peak wavelength
    wlpeak = 0.5*(float(segmask_segm.loc[wavelengthband, 'Wavelength_Start'])
                  +float(segmask_segm.loc[wavelengthband, 'Wavelength_End']))
    # Appending it to a list to add as a column to segmask_segm
    temporarywavelengthpeaklist.append(wlpeak)
    # Interpolating the resolution at the wavelength peak and replacing the resolution of that index with it
    segmask_segm.loc[wavelengthband,'Resolution']= interpolation(wlpeak)*resolution_factor
# We insert the peak list into the dataframe
segmask_segm.insert(0, "Wavelength_Peak", temporarywavelengthpeaklist)

# Checks for a negative range which I assume would cause issues. Check it works @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if min(segmask_segm['Wavelength_End']-pd.to_numeric(segmask_segm['Wavelength_Start'])) <= 0:
    print("Segment has a negative range!")
    exit()

# Checks for overlapping segments if there's more than one. Double check this works when you have data @@@@@@@@@@@@@@
if len(segmask_segm['Wavelength_End']) > 1:
    if max(segmask_segm['Wavelength_End'][0:len(segmask_segm['Wavelength_End'])]
           -pd.to_numeric(segmask_segm['Wavelength_Start'][1:len(segmask_segm['Wavelength_Start'])])) > 0:
        print("Overlapping segments")
        exit()



# Ok we are opening the 4 ccds and assigning their values to new columns as lists (wl, sob, uob)
#object_id_desired = input("Gimme the id you want: ")
object_id_desired = 150427000801049




# Reads the DR2 Sp data file. Changes to a  ''SPECTRA/'+obs_file' kinda thing later in galahsp.-
# Reads out the columns which are centre/peak wavelength, start and end of wavelength peak (all simulated), and atomic number
# Careful with the python enginer, it's slower. If we are looking at BIG data files this might be bad.
linemask_sp=pd.read_csv(
    r"C:\Users\jama2357\Documents\Galafiles\GALAH\DATA\DR2_Sp.dat",sep="[ ]{2,}",header=None,engine='python',
    names=["Sim_Wavelength_Peak","Sim_Wavelength_Start","Sim_Wavelength_End", "Atomic_Number"]
)
# Removes all rows that begin with ; as it's commented out. Required to be able to modify the data by column rather than
# by row
linemask_sp=linemask_sp[~linemask_sp['Sim_Wavelength_Peak'].str.startswith(";")]
# Reset the index to account for the now missing values
linemask_sp = linemask_sp.reset_index(drop=True)


