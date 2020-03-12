
"""This file just collects variables from other scripts and calls the most basic constants to plug into a dictionary
to input to SME at some point."""


"""So high alpha elements:Fe ratio occurs from type 2 SN as it produces alpha elements. This occured a lot
during the early universe. After a delay, 1a SN occured (as they live longer) and so increased the metalicity
and DECREASING alpha:Fe. This is a partial link back in time? But the radial distance from the galactic centre
also increases the metalicity because (?) they're less dense out there. It follows a knee pattern of alpha:Fe against
Fe:H. We see an abundance of stars outside of this pattern indicating a merger. There is a thick disk (high alpha) disk
around the galaxy of older stars that were the originals. The merger would have dispersed some of this into the thin
(low alpha) disk that has a wider radius. We also see a banana shape from where  the lower metalic galaxy influenced our
star formations."""

# Version of SME? Make a text box.
#'version = input("What version of SME are we chilling with today boys: ")
version = 5.36

# This is our version of the "id" which IDL calls with Systime()
from datetime import date
import numpy as np

currentdate = date.today()

# Effective temperature
import makestruct_cannon_parameters
# What does this even mean @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#radial_velocity_flag = input("Input value for vrad_flag. -1 for global, 0 for each wavelength segment")
radial_velocity_flag = -1
# What does it stand for? It may be changed during iterative optimisation to -3, so.. what the hell do I put it as?
# Keeping it as 1, but will have to be changed at some point. Also the ORDER of all things is important. Gotta set the flag
# before optimisation for example. @@@@@@
cscale_flag = 1
# apparently usually fixed. Find out names? @@@@@@@@@@@@@@@
gam6 = 1
accwi = 0.005
accrt = 0.005
clim = 0.01
chirat = 0.001
nmu = 7
obs_type = 3
iptype = "gauss"
# Check this with sven
mu = np.flip(np.sqrt(0.5*(2*np.arange(1,nmu+1))/nmu))
# = 0 in galah_ab, I wonder why. Also how does this relate to the alpha keyword
auto_alpha = 1

#Obs name and file setup
field_for_obs = 1 #input("Field, please.")
object_for_obs = 2 #input("Object, please.")
setup_for_obs = 3 #input("Setup, please.")

obs_name = str(field_for_obs)+"_"+str(object_for_obs)+"_"+str(setup_for_obs)+"_Sp"
obs_file = str(obs_name+".dat")

"""Here we do iterations which change variables such as vradglob and cscale flag. It's important to 
remember that they WILL change and NOT be the variables that are above or imported.
Vradglob apparently is not used in sme. I am not sure where some variables should be (before or after) so do it
on a case by case basis when reviewing this part."""

if radial_velocity_flag == 0:
       pass

import makestruct_createdatfiles

# The iteration values, something like normalise 1, iterate 2, normalise 1, iterate 20? Change later @@@@@@@@@@@@@
# These are used for maxiter, which takesthe iteration value depending on what part of the loop it's on
# Not sure how to code this currently, as we are NOT coding the loop just yet. @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
iteration_list = [1, 2, 1, 20]
current_loop_main = 1
max_iteration = iteration_list[current_loop_main]

import makestruct_abund
import makestruct_atmosphere
# Dictionary uploaded to sme
sme = {'version': version, "id":currentdate, "teff": makestruct_cannon_parameters.effective_temperature_cannon,
       "grav": makestruct_cannon_parameters.log_surface_gravity_cannon, "feh": makestruct_cannon_parameters.metalicity_cannon,
       "vmic": makestruct_cannon_parameters.microturbulence_velocity_cannon,
       "vmac": makestruct_cannon_parameters.macrotublence_velocity_cannon,
       "vsini": makestruct_cannon_parameters.rotational_velocity_cannon,
       "vrad": makestruct_cannon_parameters.radial_velocity_cannon_global,
       "vradflag": radial_velocity_flag, "cscale_flag":makestruct_createdatfiles.continuum_scale_flag,
       "cscale": makestruct_createdatfiles.continuum_scale, "gam6": gam6, "accwi": accwi, "accrt": accrt, "clim": clim,
       "maxiter": max_iteration, "chirat": chirat, "nmu": nmu, "abund": makestruct_abund.abundances, "mu":mu,
       'atmo': {'source': makestruct_atmosphere.atmosphere_grid_file, 'geom': makestruct_atmosphere.atmoshpere_geometry,
                'depth': makestruct_atmosphere.atmosphere_depth,'interp': makestruct_atmosphere.atmosphere_interpolation,
                'method': makestruct_atmosphere.atmosphere_abundance_grid},
       'sob': makestruct_createdatfiles.ccd_flux_inside_segmask_array,
       "uob": makestruct_createdatfiles.ccd_flux_error_inside_segmask_array, "obs_name": obs_name, "obs_type": obs_type,
       "iptype": iptype, "glob_free": "-1", "ab_free": np.zeros(99), "gf_free": 0
}