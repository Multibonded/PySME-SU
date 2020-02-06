import numpy as np
from astropy.io import fits

# we could get feh from when it is created (makestruct cannon parameters), but idl calls it from the abundances
# despite it making no changes to it there
from makestruct_abund import abundances, metalicity_modified
import makestruct_cannon_parameters
from galah_sp_part1 import setup_for_obs, field_for_obs, object_for_obs, ab_free, cscale_flag, line_list
from galah_sp_part3 import iterations, metalicity, log_surface_gravity

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
    # A continuum mask covering the full segment @@@@@@@@@@
    continuum_mask = str(setup_for_obs)+"_Cont.dat"
    # The segment definitions containing all lines @@@@@@@@@@
    segment_mask = str(setup_for_obs)+"_Segm.dat"

    # Define initial abundances, initialising the MARCS model abundances
    if optimisation_loop == 0:
        # Element 6 or 7? Array starts at 0, so [6] is actually the 7th element @@@@@@@@@@@@@@@
        if log_surface_gravity < 2: abundances[6] += 0.4
    else:
    # Restoring the variables and updating sme parameters

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