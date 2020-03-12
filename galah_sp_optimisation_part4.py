import numpy as np
from astropy.io import fits
# we could get feh from when it is created (makestruct cannon parameters), but idl calls it from the abundances
# despite it making no changes to it there
from makestruct_abund import abundances
import makestruct_cannon_parameters
from galah_sp_part1 import setup_for_obs, field_for_obs, object_for_obs, ab_free, continuum_scale_flag, line_list, obs_name, broad_lines, depthmin, version
from galah_sp_part2 import wavelength_res_x_collection, wavelength_res_y_collection, interpolation, resolution_factor
from galah_sp_part3 import iterations, metalicity, log_surface_gravity, effective_temperature, ebv
import pandas as pd

class read_iso():
    def __init__(self):
        self.num_cols=4
        self.columns = ['M_Mo', 'logTeff', 'logG', 'logL_Lo']
        self.num_ages = len(age)
        self.ages = age

    def fill_chemistry(self, m_h, fe_h, alpha_fe):
        self.FeH = fe_h
        self.Z = 10**m_h*0.0152
        self.aFe = alpha_fe

    def fill_iso(self, iso_input):
        self.data = iso_input
# from update_logg import optimize_logg_with_extra_input




"Main part 4, main loop: Optimisation. Normalisation into iterative optimisation."
for optimisation_loop in range(0, len(iterations)):

    # Still keeping the 3rd try coded for future if we want it
    if optimisation_loop == 0: print("STARTING LOOP 1 -> Pre-Normalisation")
    if optimisation_loop == 1: print("STARTING LOOP 2 -> First iterations (max. 2 iterations)")
    if optimisation_loop == 2: print("STARTING LOOP 3 -> Normalization")
    if optimisation_loop == 3 and len(iterations) == 6: print("STARTING LOOP 4 -> Second iterations (max. 20 iterations)")
    if optimisation_loop == 3 and len(iterations) == 4: print("STARTING LOOP 4 -> Final iterations (max. 20 iterations)")
    if optimisation_loop == 4: print("STARTING LOOP 5 -> Normalization")
    if optimisation_loop == 5: print("STARTING LOOP 6 -> Final iterations (max. 20 iterations)")

    maximum_iteration = max(iterations)
    # A continuum mask covering the full segment. Why are we setting it each loop? Its never reset. And it's only
    # ever used ONCE and that's in makestruct
    continuum_mask = str(setup_for_obs)+"_Cont.dat"
    # The segment definitions containing all lines
    segment_mask = str(setup_for_obs)+"_Segm.dat"

    # Define initial abundances, initialising the MARCS model abundances
    if optimisation_loop == 0:
        # Load the abundances here in idl, but we import earlier

        # Element 6 or 7? Array starts at 0, so [6] is actually the 7th element @@@@@@@@@@@@@@@
        if log_surface_gravity < 2: abundances[6] += 0.4
    else:
    # Restoring the variables and updating sme parameters @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    # Restore first loop and normalised spectra
        pass



    "Loop step 4.1: Segment selection and normalisation"
    if iterations[optimisation_loop] <= 1:
        # No free parameters except normalisation using straight lines
        glob_free = '-1'
        ab_free.fill(0)
        continuum_scale_flag = 1 #cscale

        "Linelist, Depth-Grid, etc."
        # Read line-list and segments
        if optimisation_loop == 0:
            # We start with the galah master file, and change it to our own fits file after (I assume) creating it? @@@
            master_line_list = fits.open(r'C:\Users\jama2357\Documents\Galafiles\GALAH\LINELIST/%s' % line_list)[1]

        print("Using trilinearly interpolated depth")
        # Here we want to call bc_ks (depth in idl) from update_logg, but it's a nested function so we have to do the
        #whole thing, and the actual value itself is not output.
        #depth = optimize_logg_with_extra_input.bc_interpolation(effective_temperature, log_surface_gravity, metalicity, ebv)

        # Why are we replacing the depth info?
        #master_line_list.header['depth'] = [depth[100471:162356],depth[318432:435589],depth[652638:748075],depth[851743:884262]]
        line_list = str(obs_name)+".fits" # Why are we replacing galah master.fits now?
        line_mask = setup_for_obs+"_Sp.dat" # Line masks for log(g) sensitive lines, i.e Fe Ti Sc

        segment_mask_data_with_res = pd.read_csv(
            "DATA/"+segment_mask, sep='[ ]{2,}| ;| [A-Z]', header=None, names=["Wavelength_Start", "Wavelength_End",
                                                                               "Resolution", "comment", "overflow"],
            engine='python', skipinitialspace=True, usecols=["Wavelength_Start", "Wavelength_End", "Resolution"])  # Here we auto delete the comment and overflow columns

        # These lines I guess are not needed when not using IDL created files? @@@@@@@@@@@@@@@@@@@@@@@@@
        # ~ asks for negation, removes any row that starts with ; in the first column. Unsure what the str does but required.
        segment_mask_data_with_res = segment_mask_data_with_res[~segment_mask_data_with_res['Wavelength_Start'].str.startswith(";")]
        # Reset the index to account for the now missing values
        segment_mask_data_with_res = segment_mask_data_with_res.reset_index(drop=True)

        # Next is: ;This sorts out segments 35 \n i=sort(seg_st) & seg_st=seg_st[i] & seg_en=seg_en[i] & seg_ipres=seg_ipres[i]
        # Sorts in ascending order of wavelength of starting wavelength, and re-orders the others.
        # But it's already ordered? weird. Will leave it in case data comes in unordered
        segment_mask_data_with_res.sort_values(by=['Wavelength_Start', 'Wavelength_End'])

        # Now we get the resolution from the resolution map.
        # We make a list to avoid calling pandas data frame repeatedly, I believe this is faster and avoids any copy errors.
        temporarywavelengthpeaklist = []
        for wavelength_band_resolution in range(0, len(segment_mask_data_with_res['Resolution'])):
            # Calculating the peak wavelength
            temporarywavelengthpeaklist.append(0.5 * (float(segment_mask_data_with_res.loc[wavelength_band_resolution, 'Wavelength_Start'])
                            + float(segment_mask_data_with_res.loc[wavelength_band_resolution, 'Wavelength_End'])))
            # Interpolating the resolution at the wavelength peak and replacing the resolution of that index with it
            segment_mask_data_with_res.loc[wavelength_band_resolution, 'Resolution'] = interpolation(temporarywavelengthpeaklist[wavelength_band_resolution])*resolution_factor

        broad_line_index = []
        if 'broad_lines' in locals():
            # list over numpy array to store the indexes of the broad lines where they equal the linelist. So the
            # broad lines that were find? @@@@@@@@@@@@@@@@@@@@@@
            # Not sure why we get the index yet, or do we just want the value of the broad lines?
            for broadline in broad_lines:
                print(np.where(broadline == master_line_list.data['lambda'])[0])

                broad_line_index.extend((np.where(broadline == master_line_list.data['lambda']))[0])
            print(broad_line_index)

        "Loop over segments, normalise them, and update the line list"
        # the list that we turn into an np array containing the indexes of the parts of the master file we care about
        # Must apply it to sme rdlin to save time.
        all_lines_index = []
        print(all_lines_index,"\nyes")
        for wavelength_band in range(0, len(segment_mask_data_with_res['Wavelength_Start'])):


            #print("Segment Number:", wavelength_band,"/", len(segment_mask_data_with_res['Wavelength_Start']))
            # "Print a dummy fits-file, covering only one segment" IDL doesn't specify the depth we use
            # are we supposed to use the [wavelength band] one? IDL just seems to do the whole array @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            # Double 'and' as there are 3 required variables to be true
            #print(master_line_list.data['depth'][wavelength_band])
            single_line = np.where(np.logical_and(np.logical_and(
                master_line_list.data['lambda'] >= float(segment_mask_data_with_res['Wavelength_Start'][wavelength_band]),
                master_line_list.data['lambda'] <= float(segment_mask_data_with_res['Wavelength_End'][wavelength_band])),
                np.logical_or(master_line_list.data['depth'] > depthmin,
                 str(master_line_list.data['name'][0]).strip() == 'H'))) # Above min depth OR is hydrogen I guess
            # If there are no non broad lines, all_lines_index are just broad, else combine the two. These are the INDEXES
            # but we turn it into the real thing when creating the smaller linelist of obsname.fits for makestruct
            # all_lines_index is plines in idl
            #print("single", single_line)
            all_lines_index.extend(single_line[0])
        # If we have broad lines in the local variable definitions we want to add them. Out of loop to prevent repeated
        # adding of the same ones.
        if 'broad_lines' in locals():
            # all lines is plines in idl. Contains the regular lines in the wavelength bands, and the broad ones
            # that impact it but with peaks that are out of the range.
            # So theoretically, it could try to concatenate b l i if it doesn't exist if the previous if statement
            # is skipped, but it can't happen as they have the same if requirements, so what's the issue?
            # np.sort to keep in numerical order.
            all_lines_index.extend(broad_line_index)
        # Avoid pesky duplicates of broad lines! Which we otherwise get. This is plines.
        all_lines_index = np.unique(np.asarray(all_lines_index))
        # So we have a csv file containing the indexes of the master line list that are within our wavelengths.
        # Idl writes this to linelist/+line_list
        np.savetxt('alllinesindex.csv',(all_lines_index), delimiter=" ")
        #print("ayeee", (master_line_list.data[all_lines_index]))
        # Writing the wavelengths to a fits file. Skip this and just make an array as we're all using the same files!
        """columndata = fits.Column(name='Obs_data', array=master_line_list.data[all_lines_index], format='E')
        tableoflines = fits.BinTableHDU.from_columns([columndata])
        tableoflines.writeto('LINELIST/'+line_list)
        exit()
        reduced_fits_file = master_line_list.data[all_lines_index]

        # New dummy file with wavelength stars and end and res for the next loop of makestruct. This name becomes
        # obs_name
        np.savetxt('DATA/'+obs_name+'_Segm.dat', np.c_(segment_mask_data_with_res['Wavelength_Start'][wavelength_band],
                                                       segment_mask_data_with_res['Wavelength_End'][wavelength_band],
                                                       segment_mask_data_with_res['Resolution'][wavelength_band]),
                   delimiter='  ')
        print("Wavelength start, end, and resolution:\n", segment_mask_data_with_res['Wavelength_Start'][wavelength_band],
              segment_mask_data_with_res['Wavelength_End'][wavelength_band], segment_mask_data_with_res['Resolution'][wavelength_band])
        # Segment mask reset to obsname in makestruct to use the file with only one segment. Seems ridiculous to do
        # it that way, why not just do [wavelengthband]

        # Run SME for the segment and restore the output
        if optimisation_loop == 0:
            # Prenormalise for initial loop
            normalisation = 1
        else:
            normalisation = 0

        if version >= 5.00:
            run = 2
            newsme = True
            print("Starting version:", version)
        else:
            run = 2
            newsme = False
            print("Starting version:", version)

        # restore OUTPUT/obsname_Sme.out
        # Guess we need to run sme now to get this?

        # Create dummy vectors to collect line-lists and normalised spectra."""
