from astropy.io import fits
import galah_sp_part1
import galah_sp_part2
import numpy as np
reduction_and_analysis_data = galah_sp_part2.reduction_and_analysis_data
"""Taking the cannon file and finding the object id desired inside it, then taking those parameters and 
modifying them based on cannon guess - if it was wrong we adjust it to a baseline - and output the parameters
This file currently just gives the parameters such as Effective Temp. Surface Grav. Metalicity. Micro Turblnce. 
Rot Vel. Radial Velocity. Macro Turb"""


# This finds Temp, Grav, and metalicity for cool dwarves. What are ai1 and stuff and where do they come from? Will
# rename after I find out details
# Should these functions be in different files? Something to think about
# Micro and macro turbulence effect the flux we see, changing our guesses on gravity and T etc
# But macro turbulence effects are so similar to radial velocity effects that our resolution cannot distinguish them
# hence we scrap the macro and just use vsini (radial velocity)
def microturbulence_macroturbulence(temperature, gravity, metalicity):
    # Prevents warning on reference being assigning, despite it not being able to happen unless there's something WRONG.
    ai1 = 0; ai2 = 0; ai3 = 0; aa1=0; ba1=0;ca1=0
    if gravity <= 4.2 or temperature >= 5500:
        ai1 = 1.1; ai2 = 1E-4; ai3 = 4E-7
        aa1 = 1.5; ba1 = 0; ca1 = 1E-6
        # ai1 = 1; ai2 = -2.5E-4; ai3 = 6.5E-7
    elif gravity >= 4.2 and temperature < 5500:
        ai1 = 1.1; ai2 = 1.6E-4; ai3 = 0
        aa1 = 1.5; ba1 = 0.2E-3; ca1 = 0
    temperature_base = 5500
    gravity_base = 4.0
    #;vmic =      ai1 + ai2*(t -t0) + ai3 * (t-t0)^2 + bi1 + bi2 * (g-g0) + bi3 * (g-g0)^2 + ci1 + ci2 * z + ci3 * z^2
    #;vmac = 3 * (aa1 + ba1*(t -t0) + ca1 * (t-t0)^2 + ba2*(g-g0) + ca2*(g-g0)^2 + ba3*z + ca3*z^2) + 2.0
    microturbulence_velocity = ai1 + ai2*(temperature-temperature_base) + \
                                      ai3 * (temperature-temperature_base)**2
    macrotublence_velocity = 3 * (aa1 + ba1*(temperature-temperature_base) + ca1
                                         * (temperature-temperature_base)**2)

    return(microturbulence_velocity, macrotublence_velocity)


'1st try Cannon, the machine learnign software that gives reasonable first guesses. Except at the edges? It gets a bit ' \
'funky. How, exactly? @@@@@@@@ what are the limitations @@@@@@@@@@'
def cannon_attempt_1(object_cannon_data):
    effective_temperature = object_cannon_data['Teff_cannon']
    # galahsp has a grav_input = this grav value but does nothing with it..? Guess it's the original gravity
    log_surface_gravity = object_cannon_data['Logg_cannon']
    origianl_log_surface_gravity = log_surface_gravity
    metalicity = object_cannon_data['Feh_cannon']
    # Same as gravity, this is unedited later.
    original_metalicity = metalicity
    microturbulence_velocity = object_cannon_data['Vmic_cannon']
    rotational_velocity = object_cannon_data['Vsini_cannon']
    macrotublence_velocity = 0
    # Ooh maybe I should keep the original names but also reassign them immediately? Idk
    # This one comes from the iraf/bailer jones, not cannon. is that right? @@@@@@@@@@@@@@@@@@@@@@@
    rv_guess = reduction_and_analysis_data['rv_guess']
    # Ooh maybe I should keep the original names but also reassign them immediately? Idk
    radial_velocity_global = rv_guess

    # These checks are for if the Cannon made a poor guess of the gravitational effect and it was off the main sequence HR
    # evolution. If that occurs we reset them to reasonable numbers.
    if effective_temperature < 5250 and log_surface_gravity > (4.25-(4.25-3.5)*(5250-effective_temperature)/1250):

        print("Cool dwarf (Upturn) [what's upturn] found, BOYS. Adjusting that sweet sweet log(g) and metalicity [Fe/H]")
        log_surface_gravity = (4.5+(0.2*(5250 - effective_temperature)/1250))
        metalicity = 0
        # Here is the vmic_vmac pro file called, translated to python earlier. Returns using the microturbulence_macro...
        # function to give the micro and macro turbulence for cool dwarves. Metalicity is unneeded with the current method.

        turbulences_list = (microturbulence_macroturbulence(effective_temperature, log_surface_gravity, metalicity))
        microturbulence_velocity = turbulences_list[0]
        macrotublence_velocity = turbulences_list[1]
        rotational_velocity = macrotublence_velocity
        # We seem to turn vmac to 0 and instead assign vsini to its value, what the what?
        macrotublence_velocity = 0.0

    # Possible parameters to distinguish between giants and dwarfs
    # Estimated with Dartmouth isochrones for metalicity = -0.5dex, age = 15Gyr 4500K, 2 dex and 3650 K, 0.5 dex
    if effective_temperature < 4500 and metalicity < -0.75 and \
            log_surface_gravity > (2-((2-0.5)*(4500-effective_temperature)/850)):
        print('Possible cool giant (TiO) identified, adjusting gravity and metalicity')
        #log_surface_gravity = 4.5 + (0.2*(5250-effective_temperature)/1250)
        metalicity = 0
        turbulences_list = (microturbulence_macroturbulence(effective_temperature, log_surface_gravity, metalicity))
        microturbulence_velocity=turbulences_list[0]
        macrotublence_velocity=(turbulences_list[1])
        # Can't differ between the two with low res, so we just stick with rot.
        rotational_velocity = macrotublence_velocity
        macrotublence_velocity = 0.0

    if effective_temperature < 4250 and log_surface_gravity < 2 and metalicity < -0.75:
        print ("Giant at end of Cannon training set identified, adjusting metalicity")
        metalicity = 0.0
    # radial_velocity_cannon_global is in sp part 2
    print("Effective Temp. Surface Grav. Metalicity. Micro Turblnce. Rot Vel. Macro Turb, Rad Vel")
    print(effective_temperature, log_surface_gravity, metalicity, microturbulence_velocity,
          rotational_velocity, macrotublence_velocity, radial_velocity_global)

    iterations = [1, 2, 1, 20]
    return effective_temperature, log_surface_gravity, origianl_log_surface_gravity, metalicity, original_metalicity, \
               microturbulence_velocity, rotational_velocity, macrotublence_velocity, radial_velocity_global, iterations

' 2nd try GUESS'
def guess_attempt_2(reduction_and_analysis_data):
    if (reduction_and_analysis_data['flag_guess']) != 0:
        print("Guess quality too low.")
        return

    effective_temperature = reduction_and_analysis_data['teff_guess']
    # galahsp has a grav_input = this grav value but does nothing with it..? Guess it's the original gravity
    log_surface_gravity = reduction_and_analysis_data['logg_guess']
    origianl_log_surface_gravity = log_surface_gravity
    metalicity = reduction_and_analysis_data['feh_guess']
    # Same as gravity, this is unedited later.
    original_metalicity = metalicity
    # vrad here
    radial_velocity_global = reduction_and_analysis_data['rv_guess']
    # spaces lines them up
    print("Guess completed\n     Object,        Teff,    Grav, Metalicity,   Vrad\n", object_id_desired,
          effective_temperature, log_surface_gravity, metalicity, radial_velocity_global)
    # Modify gravity why? @@@@@@
    if log_surface_gravity > 3.5:
        log_surface_gravity = log_surface_gravity + 0.5
        print("NB: Offset to initial GUESS log(g) by +0.5 for log(g) > 3.5")
    # Who can know why eh? @@@@@@@@
    metalicity = metalicity + 0.15
    print("NB: Offset to initial GUESS [Fe/H] by +0.15")

    # Doesn't work so well on these stars so we modify them slightly? @@@@@@@@@@@@
    if effective_temperature < 5250 and log_surface_gravity > (
            4.25 - (4.25 - 3.5) * (5250 - effective_temperature) / 1250):
        print(
            "Cool dwarf (Upturn) found. Adjusting log(g) and metalicity [Fe/H]") #what's upturn @@@@@@@
        log_surface_gravity = (4.5 + (0.2 * (5250 - effective_temperature) / 1250))
        metalicity = 0

    if effective_temperature < 4250 and log_surface_gravity < 2.0 and metalicity < -0.75:
        print("Giant below 4250K identified, adjusting [Fe/H] to 0") # Why? @@@@@@@@@@@@
        metalicity = 0

    turbulences_list = (microturbulence_macroturbulence(effective_temperature, log_surface_gravity, metalicity))
    microturbulence_velocity = turbulences_list[0]
    macrotublence_velocity = turbulences_list[1]
    rotational_velocity = macrotublence_velocity
    # We seem to turn vmac to 0 and instead assign vsini to its value, what the what?
    macrotublence_velocity = 0.0

    # Four SME calls : normalise, iterate, normalise, iterate
    iterations = [1, 2, 1, 20]

    if reduction_and_analysis_data['red_flag'] != 0:
        #bit mask for 1 to 8, summing um CCD problem, e.g. (0011) = 1+2=3 for BG @@@@@@@@@@ what @@@@@@@@@@@@
        print("Reduction issue! Red flag: ", reduction_and_analysis_data['red_flag'])
        if reduction_and_analysis_data['red_flag'] == 14:
            print("Skyflat!")

    return effective_temperature, log_surface_gravity, origianl_log_surface_gravity, metalicity, original_metalicity, \
           microturbulence_velocity, rotational_velocity, macrotublence_velocity, radial_velocity_global, iterations


# Skipping for now
def gaussian_attempt_3():
    print("Neither Cannon nor GUESS appropriate, trying Gaussian approach")
    # We get the wavelength and flux arrays from gala sp part 2




# Input because I hate the idea of having this be permanent and wrong
# cannon_version = input("Cannon version: ")
cannon_version = 'iDR2'
#object_id_desired = input("Gimme the id you want: ")
object_id_desired = 150427000801049

cannon_data = fits.open(r"sobject_iraf_iDR2_cannon_small.fits")
#print(repr(cannon_data[1].header))

# We want to isolate the sobject_id which is under TTYPE1
# We will use cannon_data[1].data to take the column we want. straight data[0] just takes the first ROW.
# print(cannon_data[1].columns) shows all column headers and their format. We specifcy not the row number index but the
# string title to print and take its column data.
# Find the index of when the id is what we want
object_index = np.where(cannon_data[1].data['sobject_id']==object_id_desired)
# This gives the row data when the object id is equal to the desired int. Only took an hour to figure out.
object_cannon_data = (cannon_data[1].data[object_index])
# If .size is bigger than 0, it exists. If it == 0, can't find it and so we stop the script.
if not object_index[0].size :
    print("Does not exist")
    exit()

# Now we wan to find a few variables from our row. We could take the index of the columns we want and apply it to our
# data OR we could take our row index and apply it to the dataframe. I have chosen the 2nd option currently.
# Should put these as options for user input?

# this is the quality of the cannon approx. Above 1 and we want to do a different method. 0 is ideal.
cannon_quality = (object_cannon_data['flag_cannon'])



print("Cannon flag quality is at ",cannon_quality)
print("Guess flag quality is at ", reduction_and_analysis_data['flag_guess'])

#temp
#radial_velocity_global = -10

if cannon_quality <= 1:
    # Takes these variables from the function
    effective_temperature, log_surface_gravity, origianl_log_surface_gravity, metalicity, original_metalicity, \
    microturbulence_velocity, rotational_velocity, macrotublence_velocity, radial_velocity_global, iterations = \
        cannon_attempt_1(object_cannon_data)

elif cannon_quality > 1:
    print("Cannon approx. quality is TOO LOW.")
    # Takes these variables from the function

    effective_temperature, log_surface_gravity, origianl_log_surface_gravity, metalicity, original_metalicity, \
    microturbulence_velocity, rotational_velocity, macrotublence_velocity, radial_velocity_global, iterations =\
        guess_attempt_2(object_cannon_data)




' Main part 3.2/5'
"""Account for special modes: LBOL, SEIS, FIXED
If Seis, prepare asteroseismic info
If LBOL, prepare photomettric/astrometric info
if Lbol/Seis, update logG based on provided info."""
print(str(reduction_and_analysis_data['ph_qual_tmass']))
# If LBOL mode: Galah + GaiaDR2
if (galah_sp_part1.field_for_obs[-4:]).lower() == 'lbol':
    if reduction_and_analysis_data['r_est'] > 0:
        print("Star in Gaia DR2 (And possible 2MASS and WISE")
        # What are the real names of these variables
        j_mag = reduction_and_analysis_data['j_m']
        e_j_mag = reduction_and_analysis_data['j_msigcom']
        h_mag = reduction_and_analysis_data['h_m']
        e_h_mag = reduction_and_analysis_data['h_msigcom']
        k_mag = reduction_and_analysis_data['ks_m']
        e_k_mag = reduction_and_analysis_data['ks_msigcom']
        w2_mag = reduction_and_analysis_data['w2mpro']
        e_w2_mag = reduction_and_analysis_data['w2mpro_error']
        ebv = reduction_and_analysis_data['ebv']
        a_k = 0.918*(h_mag - w2_mag - 0.08)
        # why

        if a_k < 0:
            a_k = 0
        print("A_Ks from RJCE:", str(a_k), '+/-', str(e_h_mag**2 + e_w2_mag**2))
        print("Quality flags (2MASS+WISE)", str(reduction_and_analysis_data['ph_qual_tmass']), '+&',
                                                 str(reduction_and_analysis_data['ph_qual_wise']))
        print("A_Ks from 0.36*E(B-V)", str(0.36*ebv))
        # What this  Need [0] first as is a list and we take the first value (only value), is [1] correct? @@@@@@~~@@@@
        if np.isinf(float(h_mag)) or np.isinf(float(w2_mag)) or \
                    str(reduction_and_analysis_data['ph_qual_tmass'][0][1]) != 'A' or \
                    str(reduction_and_analysis_data['ph_qual_wise'][0][1]) != 'A':
            a_k = 0.36*ebv

        if 0.36*ebv > 3*a_k:
            ebv = 2.78*a_k
            print("E(B-V) not trustworthy, setting to 2.78*A_Ks =", str(ebv), "instead of",
                  str(reduction_and_analysis_data['ebv']))

        # Why isn't this with the others?
        plx = reduction_and_analysis_data['parallax']
        e_plx = reduction_and_analysis_data['parallax_error']
        dist = reduction_and_analysis_data['r_est']
        print("Distance", str(dist), "with limits", str(reduction_and_analysis_data['r_lo']), "to",
                                                          str(reduction_and_analysis_data['r_hi']))
        # uhh this is 1000 over parallax, not 1 over it... @@@@@@@@@@@@@
        print("1/parallax:", str(1000/plx), "and parallax:", str(plx), "with", str(100*e_plx/plx), "% uncertainty")
        e_dist = 0.5*(reduction_and_analysis_data['r_hi'] - reduction_and_analysis_data['r_lo'])
        if dist < 100:
            print("Star within 100pc, setting A_K and EBV to 0")
            a_k = 0
            ebv = 0
    else:
        print("Star not in GAaia DR2, but lbol-keyword activated. Cancel.")
        exit()

"If LBOL mode: Sun at 10pc"
# Could just do the last 8 letters?
if galah_sp_part1.field_for_obs == 'sun_lbol' or galah_sp_part1.field_for_obs == 'RV_sun_lbol':
    # "Here we assume the sun is at 10pc, duhh!" <- Really guys, THAT'S the one comment you want to put in the last
    # 100 lines of code?
    k_mag = 3.28 - 5 * np.log(0.1) - 5
    h_mag = 3.32
    plx = 100 # "mas" ??? @@@@
    e_plx = 0
    dist = 10 # assuming 10 pc
    e_dist = 0.01
    ebv = 0

"If SEIS mode: get numax. (What's numax?)"
if galah_sp_part1.field_for_obs[-4:] == 'seis':
    try:
        numax = (reduction_and_analysis_data['numax'])
        e_numax = 0
    except KeyError as missingvalue:
        print(missingvalue)
        exit()

"Update initial LogG if LBOL/SEIS mode"
if galah_sp_part1.field_for_obs[-4:] == 'lbol':
    log_surface_gravity = log_surface_gravity #update log g here
    log_surface_gravity = log_surface_gravity[0] # Why are we doing these [0] NOW and not before? no issue with them@@@@@@@@@
    turbulences_list = microturbulence_macroturbulence(effective_temperature, log_surface_gravity, metalicity)
    microturbulence_velocity = turbulences_list[0]
    radial_velocity_global = turbulences_list[1]
    macrotublence_velocity = 0
    print("Running in LBOL mode, glob_free are TEF, FEH, VRAD, (VSINI), GRAV adjusted from TEFF, KMAG, PLX")
    print("LBOL =", str(galah_sp_part1.lbol), "Lbol_sun, KMAG =", str(k_mag[0]), "mag, DISTANCE =", str(dist[0]),
          ",PLX =", str(plx[0]), "mas, E_PLX/PLX =", str(e_plx[0]/plx[0]))
if galah_sp_part1.field_for_obs[-4:] == 'seis':
    # We don't even do the [0] here!
    log_surface_gravity = log_surface_gravity  # update log g here
    turbulences_list = microturbulence_macroturbulence(effective_temperature, log_surface_gravity, metalicity)
    microturbulence_velocity = turbulences_list[0]
    radial_velocity_global = turbulences_list[1]
    macrotublence_velocity = 0
    # Nu is the frequency of the oscillation/
    print("Running in SEISMIC mode, glob_free are TEF, FEH, VRAD, (VSINI), GRAV adjusted from TEFF, NU_MAX")
    print("NU_MAX =", str(numax), ", E_NU_MAX/NU_MAX =", str(e_numax/numax))

# Don#t need this part I believe @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#"Update TEFF/GRAV/FEH for FIXED mode"
#if 'fixed' in galah_sp_part1.field_for_obs:

"Main part 3.3/5" \
"Ensure reasonable parameters"
# setting the boundaries of reasonable parameters and resetting the variables
# print(effective_temperature)
effective_temperature = np.clip(effective_temperature,6000,7500)
log_surface_gravity = np.clip(log_surface_gravity, 0.5, 5)
metalicity = np.clip(metalicity, -3, 0.5)
# Adjust teff and logg to marcs lower grid limits. It seems odd that we're so specific and hard coded?@@@@@@@@
if galah_sp_part1.atmosphere_grid_file == 'marcs2014.sav':
    # Sets grav to the maximum of either 2 or gravity (basically clipping with only a minimum value, no max)
    if 7500 >effective_temperature > 6500: log_surface_gravity = np.maximum(2, log_surface_gravity)
    elif 6500 >effective_temperature > 5750: log_surface_gravity = np.maximum(1.5, log_surface_gravity)
    elif 5750 >effective_temperature > 5000: log_surface_gravity = np.maximum(1, log_surface_gravity)
    elif 5000 >effective_temperature > 4000: log_surface_gravity = np.maximum(0.5, log_surface_gravity)

# Do the same but with the stagger file limits
elif galah_sp_part1.atmosphere_grid_file == 'stagger-t5havni_marcs64.sav' or \
        galah_sp_part1.atmosphere_grid_file == 'stagger-t5havni_marcs.sav' or \
        galah_sp_part1.atmosphere_grid_file == 'stagger-tAA.sav':

    # Sets grav to the maximum of either 2 or gravity (basically clipping with only a minimum value, no max)
    if 4750 > effective_temperature: log_surface_gravity = np.maximum(1.75, log_surface_gravity)
    elif 5250 > effective_temperature > 4750: log_surface_gravity = np.maximum(2.75, log_surface_gravity)
    elif 5750 > effective_temperature > 5250: log_surface_gravity = np.maximum(3.25, log_surface_gravity)
    elif 6250 > effective_temperature > 5750: log_surface_gravity = np.maximum(3.75, log_surface_gravity)
    elif effective_temperature >= 6250: log_surface_gravity = np.maximum(4.15, log_surface_gravity)

# What is this and why are we printing it exactly? @@@@@@@@@@@@
print(galah_sp_part1.field_for_obs, "_", galah_sp_part1.object_for_obs,effective_temperature ,log_surface_gravity,
      metalicity, radial_velocity_global, "Initial", "Sp")
