import numpy as np
from astropy.io import fits
# we could get feh from when it is created (makestruct cannon parameters), but idl calls it from the abundances
# despite it making no changes to it there
from makestruct_abund import abundances
import makestruct_cannon_parameters
from galah_sp_part1 import setup_for_obs, field_for_obs, object_for_obs, ab_free, cscale_flag, line_list, obs_name, broad_lines, depthmin
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
from update_logg import optimize_logg_with_extra_input




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
    # A continuum mask covering the full segment
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


        cscale_flag = 1

        "Linelist, Depth-Grid, etc."
        # Read line-list and segments
        if optimisation_loop == 0:
            linefull = fits.open(r'C:\Users\jama2357\Documents\Galafiles\GALAH\LINELIST/%s' % line_list)[1]

        print("Using trilinearly interpolated depth")
        # Here we want to call bc_ks (depth in idl) from update_logg, but it's a nested function so we have to do the
        #whole thing, and the actual value itself is not output.
        #depth = optimize_logg_with_extra_input.bc_interpolation(effective_temperature, log_surface_gravity, metalicity, ebv)

        # Why are we replacing the depth info?
        #linefull.header['depth'] = [depth[100471:162356],depth[318432:435589],depth[652638:748075],depth[851743:884262]]
        line_list = str(obs_name)+".fits" # Why are we replacing galah master.fits now?
        line_mask = setup_for_obs+"_Sp.dat" # Line masks for log(g) sensitive lines, i.e Fe Ti Sc

        segmask_segm = pd.read_csv(
            r"C:\Users\jama2357\Documents\Galafiles\GALAH\DATA\DR2_Segm.dat", sep='[ ]{2,}| ;| [A-Z]', header=None,
            names=["Wavelength_Start", "Wavelength_End", "Resolution", "comment", "overflow"], engine='python',
            skipinitialspace=True,
            usecols=["Wavelength_Start", "Wavelength_End",
                     "Resolution"])  # Here we auto delete the comment and overflow columns

        # These lines I guess are not needed when not using IDL created files? @@@@@@@@@@@@@@@@@@@@@@@@@
        # ~ asks for negation, removes any row that starts with ; in the first column. Unsure what the str does but required.
        segmask_segm = segmask_segm[~segmask_segm['Wavelength_Start'].str.startswith(";")]
        # Reset the index to account for the now missing values
        segmask_segm = segmask_segm.reset_index(drop=True)

        # Next is: ;This sorts out segments 35 \n i=sort(seg_st) & seg_st=seg_st[i] & seg_en=seg_en[i] & seg_ipres=seg_ipres[i]
        # Sorts in ascending order of wavelength of starting wavelength, and re-orders the others.
        # But it's already ordered? weird. Will leave it in case data comes in unordered
        segmask_segm.sort_values(by=['Wavelength_Start', 'Wavelength_End'])

        # Now we get the resolution from the resolution map.
        # We make a list to avoid calling pandas data frame repeatedly, I believe this is faster and avoids any copy errors.
        temporarywavelengthpeaklist = []
        for wavelength_band in range(0, len(segmask_segm['Resolution'])):
            # Calculating the peak wavelength
            temporarywavelengthpeaklist.append(0.5 * (float(segmask_segm.loc[wavelength_band, 'Wavelength_Start'])
                            + float(segmask_segm.loc[wavelength_band, 'Wavelength_End'])))
            # Interpolating the resolution at the wavelength peak and replacing the resolution of that index with it
            segmask_segm.loc[wavelength_band, 'Resolution'] = interpolation(temporarywavelengthpeaklist[wavelength_band])*resolution_factor
        # list over numpy array to store the indexes of the broad lines where they equal the linelist. So the
        # broad lines that were find? @@@@@@@@@@@@@@@@@@@@@@
        # Not sure why we get the index yet, or do we just want the value of the broad lines?
        broad_line_index = []
        for broadline in broad_lines:
            broad_line_index.append((np.where(broadline == linefull.data['lambda'])))

        print("ayeee", broad_line_index)
        exit()



        "Loop over segments, normalise them, and update the line list"
        for wavelength_band in range(0, len(segmask_segm['Wavelength_Start'])):
            print("Segment Number:", wavelength_band,"/", len(segmask_segm['Wavelength_Start']))

            # "Print a dummy fits-file, covering only one segment"
            single_line = np.where(np.logical_and(linefull.data['lambda'] >= segmask_segm['Wavelength_Start'][wavelength_band],
                                                  linefull.data['lambda']) >= segmask_segm['Wavelength_Start'][wavelength_band],
                                   (linefull.data['depth'] > depthmin or str(linefull.data['name'][0]).strip() == 'H')) # Above min depth OR is hydrogen I guess
            # If the array is empty, just set to int of 0
            if single_line[0].size == 0: single_line = 0
            if 'broad_lines' in locals(): pl