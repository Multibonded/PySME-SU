from astropy.io import fits

import numpy as np






# Input because I hate the idea of having this be permanent and wrong
# cannon_version = input("Cannon version: ")
cannon_version = 'iDR2'
#object_id_desired = input("Gimme the id you want: ")
object_id_desired = 150427000801049

cannon_data = fits.open(r"C:\Users\jama2357\Documents\Galafiles\GALAH\DATA\sobject_iraf_iDR2_cannon_small.fits")
print(repr(cannon_data[1].header))

# We want to isolate the sobject_id which is under TTYPE1
# We will use cannon_data[1].data to take the column we want. straight data[0] just takes the first ROW.
# print(cannon_data[1].columns) shows all column headers and their format. We specifcy not the row number index but the
# string title to print and take its column data.
# Find the index of when the id is what we want
index = np.where(cannon_data[1].data['sobject_id']==object_id_desired)
# This gives the row data when the object id is equal to the desired int. Only took an hour to figure out.
object_cannon_data = (cannon_data[1].data[index])
print(index[0].size)
# If .size is bigger than 0, it exists. If it == 0, can't find it and so we stop the script.
if not index[0].size :
    print("Does not exist")
    exit()

# Now we wan to find a few variables from our row. We could take the index of the columns we want and apply it to our
# data OR we could take our row index and apply it to the dataframe. I have chosen the 2nd option currently.
# Should put these as options for user input?

# this is the quality of the cannon approx. Above 1 and we want to do a different method, but I won't implement that yet.
cannon_quality = (cannon_data[1].data['flag_cannon'][index])
if cannon_quality > 1:
    print("Cannon approx. quality is TOO LOW. TRY SUMMIN ELSE.")
    exit()
print("Cannon flag quality is at ",cannon_quality)


effective_temperature_cannon = cannon_data[1].data['Teff_cannon'][index]
log_surface_gravity_cannon = cannon_data[1].data['Logg_cannon'][index]
metalicity_cannon = cannon_data[1].data['Feh_cannon'][index]
microturbulence_velocity_cannon = cannon_data[1].data['Vmic_cannon'][index]
rotational_velocity_cannon = cannon_data[1].data['Vsini_cannon'][index]
macrotublence_velocity_cannon = 0


# Why is galah_sp calling on BailerJones k2seis fits file? What does it have that we don't in cannon?
# Must changed iraf and stuff variable names once I find out what they stand for.
iraf_dr53=fits.open(r'C:\Users\jama2357\Documents\Galafiles\GALAH\DATA\sobject_iraf_53_2MASS_GaiaDR2_WISE_PanSTARRSDR1_'
                    r'BailerJones_K2seis_small.fits')
# Finds the index of where we can find our object id.
iraf_index = np.where(cannon_data[1].data['sobject_id']==object_id_desired)
# Checks the index exists therefore the object id exists in the place we're looking for.
if not index[0].size :
    print("Does not exist")
    exit()

# Change vbary name! Oh wait we don't need it...
vbary = iraf_dr53[1].data['v_bary'][iraf_index]
# There was a choice from 1-4, an e and a shift as well. I went with rv_guess because that's what the IDL calls for.
# Change the name! Radial velocity? But then what's glob? (vradglob)
rv_guess = iraf_dr53[1].data['rv_guess'][iraf_index]
vradglob = rv_guess

# Just viewing them

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
    microturbulence_velocity_cannon = ai1 + ai2*(temperature-temperature_base) + \
                                      ai3 * (temperature-temperature_base)**2
    macrotublence_velocity_cannon = 3 * (aa1 + ba1*(temperature-temperature_base) + ca1
                                         * (temperature-temperature_base)**2)

    return(microturbulence_velocity_cannon, macrotublence_velocity_cannon)

# These checks are for if the Cannon made a poor guess of the gravitational effect and it was off the main sequence HR
# evolution. If that occurs we reset them to reasonable numbers.
if effective_temperature_cannon <= 5250 and \
        log_surface_gravity_cannon >= (4.25-(4.25-3.5)*(5250-effective_temperature_cannon)/1250):

    print("Cool dwarf (Upturn) [what's upturn] found, BOYS. Adjusting that sweet sweet log(g) and metalicity [Fe/H]")
    log_surface_gravity_cannon = (4.5+(0.2*(5250 - effective_temperature_cannon)/1250))
    metalicity_cannon = 0
    # Here is the vmic_vmac pro file called, translated to python earlier. Returns using the microturbulence_macro...
    # function to give the micro and macro turbulence for cool dwarves. Metalicity is unneeded with the current method.

    turbulences_list = (microturbulence_macroturbulence(
        effective_temperature_cannon, log_surface_gravity_cannon, metalicity_cannon))
    microturbulence_velocity_cannon=turbulences_list[0]
    macrotublence_velocity_cannon=(turbulences_list[1])
    rotational_velocity_cannon = macrotublence_velocity_cannon

    # We seem to turn vmac to 0 and instead assign vsini to its value, what the what?
    vsini = macrotublence_velocity_cannon
    macrotublence_velocity_cannon = 0.0

# Possible parameters to distinguish between giants and dwarfs
# Estimated with Dartmouth isochrones for metalicity = -0.5dex, age = 15Gyr 4500K, 2 dex and 3650 K, 0.5 dex
if effective_temperature_cannon < 4500 and metalicity_cannon < -0.75 and \
        log_surface_gravity_cannon > (2-((2-0.5)*(4500-effective_temperature_cannon)/850)):
    print('Possible cool giant (TiO) identified, adjusting gravity and metalicity')
    #log_surface_gravity_cannon = 4.5 + (0.2*(5250-effective_temperature_cannon)/1250)
    metalicity_cannon = 0
    turbulences_list = (microturbulence_macroturbulence(
        effective_temperature_cannon, log_surface_gravity_cannon, metalicity_cannon))
    microturbulence_velocity_cannon=turbulences_list[0]
    macrotublence_velocity_cannon=(turbulences_list[1])
    # Can't differ between the two with low res, so we just stick with rot.
    rotational_velocity_cannon = macrotublence_velocity_cannon
    macrotublence_velocity_cannon = 0.0

if effective_temperature_cannon < 4250 and log_surface_gravity_cannon < 2 and metalicity_cannon < -0.75:
    print ("Giant at end of Cannon training set identified, adjusting metalicity")
    metalicity_cannon = 0.0

print("Effective Temp. Surface Grav. Metalicity. Micro Turblnce. Rot Vel. Radial Velocity. Macro Turb")
print(effective_temperature_cannon, log_surface_gravity_cannon, metalicity_cannon, microturbulence_velocity_cannon,
      rotational_velocity_cannon, vradglob, macrotublence_velocity_cannon)
