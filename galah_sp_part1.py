import numpy as np








#Obs name and file setup
field_for_obs = "M67lbol" #input("Field, please.")
object_for_obs = 150427000801049 #input("Object, please.")
setup_for_obs = "dr2" #input("Setup, please.")

obs_name = str(field_for_obs)+"_"+str(object_for_obs)+"_"+str(setup_for_obs)+"_Sp"
obs_file = str(obs_name+".dat")
print("obs name", obs_name)
# Not sure what this is, just leavin' it here.
specformat     = '(d9.4,2x,3d25.8)'
infoformat     = '(a35,7f8.2,i4,a4)'
segmformat     = '(2d10.2,f10.0)'
maskformat     = '(d10.4,2d10.2)'


version        = 5.36
#  flag for the radial velocity to determine how to process it
# -1 is a global vrad, 0 is a separate vrad for each wavelength segment. Not sure where it's change.d
radial_velocity_flag      = -1
# Flag for the scale of the continuum line, cscale_flag
continuum_scale_flag    = 1
# visini is, I believe a replacement for rotational velocity essentially
rotational_velocity          = 2.0
radial_velocity_global       = 0.0
# Array of flags (0 or 1), specifying for which of the 99 atomic elements the
    # abundances are free parameters.
atomic_abundances_free_parameters        = np.zeros(99) # ab_free
auto_alpha     = 1
depthmin       = 0.1
print('Using minimum depth ',depthmin)
broad_lines    = [4861.3230,6562.7970, 5167.3216,5172.6843,5183.6042] # Hbeta, Halpha, Mg Ib triplet
line_cores     = [8498.0200]
#Literally just the same as object for obs...
object_pivot   = object_for_obs
# what the f? It wants an imaginary number...? (SQRT -1) @@@@@@@@@@@@@@@@@@@@@@@
NaN            = abs(np.sqrt(1))
LBOL           = NaN
# Gotta figure these ones out lul. @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
b_mag = NaN
v_mag = NaN
g_mag = NaN
j_mag = NaN
h_mag = NaN
k_mag = NaN
w2_mag= NaN
plx   = NaN
dist  = NaN
e_dist= NaN
lbol  = NaN
e_lbol= NaN
ebv   = NaN

"""Choosing atmoshpere, grids, and NLTE."""
# These files give us the atmosphere grids to use. What do atmogrids actually DO? @@@@@@@@@@@@@@@
atmosphere_grid_file = 'marcs2014.sav'
# grid
atmosphere_abundance_grid = ['marcs2012_H2018.grd', 'marcs2012_Fe2016.grd']
line_list = 'galah_master_v5.1.fits'

"So we use these element numbers to run through the element flag array and change its values? Why not just make an array" \
"this size to begin with? @@@@@@@@@@@@@@@@@@@@@@@@@@@"
# elem. I believe the number of elements we're looking at? I am assuming it is actually 0 to 26, not just 0 and 26
# ask sven and karin @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
element_numbers = [0, 26]
nonlte_element_flags = np.zeros(99)
# In IDL it's [elem-1] whats with the -1? Why do they want a pointer (1b).. @@@@@@@@@@@@@@@@@@@@@@@
for element_flag in range(0,max(element_numbers)): nonlte_element_flags[element_flag] = 1

# nlte_grids. Object to allow us to parse in lists as he object for the array. @@@@@@@@@@@@@@@@@@@
nonlte_atmosphere_grids = np.zeros(99, dtype=object)
for element_flag in range(0,max(element_numbers)): nonlte_atmosphere_grids[element_flag] = atmosphere_abundance_grid

# Just nlte in idl. 1 is on, 0 is off.
nonlte_flag = 1
print("Using atmoshpere grid ", atmosphere_grid_file, "\nUsing linelist ", line_list, "\nUsing NLTE grid ", atmosphere_abundance_grid)
