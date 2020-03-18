

"""it seems that we only need to do the bottom 10% where it says #we have a fits file"""

import numpy as np
from scipy.stats import norm
from numpy import random as rand
import matplotlib.pyplot as plt
import time
from astropy.io import fits
from galah_sp_part1 import line_list # can change to from 4 import master line list
import pickle

class LineMerging(object):


    def __init__(self, master_file):
        # names in idl, from isotope mix
        self.element_list = np.array(
            ('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl',
             'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As',
             'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
             'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb',
             'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl',
             'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk',
             'Cf', 'Es'))
        self.master_line_list = master_file
        # The important arrays that this class was created to produce
        # Line atomic has duplicates, why?
        self.line_atomic, self.lande_mean, self.depth, self.data_reference_array, self.species, self.j_e_array, self.lower_level, \
        self.upper_level, self.lu_lande = self.readline(self.master_line_list)

        # Saving the data produced here, as it is NOT reliant on anything but the master file, so no point recreating it every damn time.
        # Ah but it does depend on the line index we're looking at, so maybe we should do this for everything and then
        # apply the index we want afterwards.
        pickle.dump([self.line_atomic, self.lande_mean, self.depth, self.data_reference_array, self.species,
                     self.j_e_array, self.lower_level, self.upper_level, self.lu_lande], open("galah_master_file_arrays", "wb"))
    # this is isotopic_mix.pro
    def isotopic_mix(self, species, master_line_list, rows_in_line_list):
        # isotope number and isotope fraction
        isotope_number = np.zeros((len(self.element_list), 10))
        isotope_fraction = np.zeros((len(self.element_list), 10))

        # This could be a dictionary with the names of the elements associated with the values, buuuut I guess that's not
        # how it's done.
        #Li
        isotope_number[3 - 1, [0, 1]] = [6, 7]
        isotope_fraction[3 - 1, [0, 1]] = [7.59, 92.41]
        # Light element in stars is actgually not very abundant so we replace it with 0%
        isotope_fraction[3 - 1, [0, 1]] = [0, 100.]
        #C
        isotope_number[6 - 1, [0, 1]] = [12, 13]
        isotope_fraction[6 - 1, [0, 1]] = [98.8938, 1.1062]
        #N
        isotope_number[7 - 1, [0, 1]] = [14, 15]
        isotope_fraction[7 - 1, [0, 1]] = [99.771, 0.229]
        #O
        isotope_number[8 - 1, [0, 1, 2]] = [16, 17, 18]
        isotope_fraction[8 - 1, [0, 1, 2]] = [99.7621, 0.0379, 0.2000]
        #Mg
        isotope_number[12 - 1, [0, 1, 2]] = [24, 25, 26]
        isotope_fraction[12 - 1, [0, 1, 2]] = [78.99, 10.00, 11.01]
        #Si
        isotope_number[14 - 1, [0, 1, 2]] = [28, 29, 30]
        isotope_fraction[14 - 1, [0, 1, 2]] = [92.2297, 4.6832, 3.0872]
        #Ti
        isotope_number[22 - 1, [0, 1, 2, 3, 4]] = [46, 47, 48, 49, 50]
        isotope_fraction[22 - 1, [0, 1, 2, 3, 4]] = [8.25, 7.44, 73.72, 5.41, 5.18]
        #Cu
        isotope_number[29 - 1, [0, 1]] = [63, 65]
        isotope_fraction[29 - 1, [0, 1]] = [69.17, 30.83]
        #Zr
        isotope_number[40 - 1, [0, 1, 2, 3, 4]] = [90, 91, 92, 94, 96]
        isotope_fraction[40 - 1, [0, 1, 2, 3, 4]] = [51.45, 11.22, 17.15, 17.38, 2.80]
        #Ba
        isotope_number[56 - 1, [0, 1, 2, 3, 4, 5, 6]] = [130, 132, 134, 135, 136, 137, 138]
        isotope_fraction[56 - 1, [0, 1, 2, 3, 4, 5, 6]] = [0.106, 0.101, 2.417, 6.592, 7.854, 11.232, 71.698]
        #La
        isotope_number[57 - 1, [0, 1]] = [138, 139]
        isotope_fraction[57 - 1, [0, 1]] = [0.091, 99.909]
        #Pr 2

        isotope_number[59 - 1, [0]] = [141]
        isotope_fraction[59 - 1, [0]] = [100.]
        #Nd 2

        isotope_number[60 - 1, [0, 1, 2, 3, 4, 5, 6]] = [142, 143, 144, 145, 146, 148, 150]
        isotope_fraction[60 - 1, [0, 1, 2, 3, 4, 5, 6]] = [27.044, 12.023, 23.729, 8.763, 17.130, 5.716, 5.596]
        #Sm 2

        isotope_number[62 - 1, [0, 1, 2, 3, 4, 5, 6]] = [144, 147, 148, 149, 150, 152, 154]
        isotope_fraction[62 - 1, [0, 1, 2, 3, 4, 5, 6]] = [3.07, 14.99, 11.24, 13.82, 7.38, 26.75, 22.75]
        #Eu 2

        isotope_number[63 - 1, [0, 1]] = [151, 153]
        isotope_fraction[63 - 1, [0, 1]] = [47.81, 52.19]
        import time
        t1 = time.time()
        tee = 0
        tff = 0
        tisin=0
        tapp = 0
        tarr = 0
        # A list to recreate the compound species in one word (TiO, SiH etc)
        # we repeat for each row in the line list galah file
        #for master_list_row in range(10000, 20000):
        for master_list_row in range(0, rows_in_line_list):
            t00 = time.time()
            # A list to recreate the compound species in one word (TiO, SiH etc). Faster than for i in.
            species.append("".join(master_line_list.data['name'][master_list_row]))
            t01 = time.time()
            tapp += t01 - t00

            t3 = time.time()
            # we strip the whitespace to be able to match it to the element list.
            # This takes about 1/3 of the time? @@@@@@@@@@@@@@
            t4 = time.time()
            tee += t4-t3

            # if there are any elements in the line.
            # the index of the elements in the element list, which we use to index through isotope_number/fraction
            # to get their values. Still think a dictionary would be better. @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            t7 = time.time()
            element_index = np.isin(self.element_list, master_line_list.data['name'][master_list_row])
            # The fits file (master line list) contains the isotopes that are used for this row. So it will contain 16
            # if we have O16 rather than O17, we use this index to find its fraction percentage that exists.
            isotope_index = np.isin(isotope_number[element_index], master_line_list.data[master_list_row]['isotope'])
            t8 = time.time()
            tisin += t8-t7
            """
            debugging 
            print("element list with element index", self.element_list[element_index])
            print("isotope numbers",isotope_number[element_index])
            print("isotope frac",isotope_fraction[element_index])
    
            print("master line list isotope values", master_line_list.data[i]['isotope'])
            print("Isotope index?", isotope_index)
            print("Isotope fractions of those indices?", isotope_fraction[element_index][isotope_index])
    
            print(master_line_list.data[i]['log_gf'])
            """
            # The isotope fractions appear to have identical effects on the log gf so we loop through them and apply
            # them to the log gf of that line.
            # This... just multiplies by the fraction. the rest does nothing?! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            # searches for the index of the element, and which isotope is desired, then takes the fraction it takes up
            # (/100 for the percantage)
            t5 = time.time()
            # for each fraction in the isotope fraction index that we've compared with the fits file
            for fraction in isotope_fraction[element_index][isotope_index]:
                if fraction:
                    # Updates log gf according to fraction. log10(10**) makes NO difference though so that's weird. @@@@
                    master_line_list.data[master_list_row]['log_gf'] = np.log10(10**(master_line_list.data[master_list_row]['log_gf'])*fraction/1E2)
            t6 = time.time()
            tff += t6-t5
        t2 = time.time()
        print("total time",t2-t1,"\ntee:", tee, "\ntff:", tff, "\ntisin", tisin)
        print("tapp", tapp, tarr)
        return species

    # eleparse in idl. It takes in the names of the molecules/elements from the galah file from 'species' as names and
    # if it's a molecule it's assigned a value of 100 and ionisation of 0, if there's more than one species it's separated
    # and the number afterwards is taken as a value of something or other.
    def element_parsing(self, species_input, atom_index, master_line_list):
        # Similar as created in iso mix but with diff names related to what it says in the idl code. This one is elemparse
        # Is it possibly calling a different name? Because they test for molecules with more than one element...
        atom_number = np.zeros((len(species_input)))
        ion_stage = np.zeros((len(species_input)))


        """ We set the atomic number and ion value of each atom in master file, or to 100 and 0 if it's a molecule.
        Definitely feel like we could combine with iso mix, and it would benefit from a dictionary. Although I guess the 
        indexing itself is basically like a list of 0-100"""
        for i in range(0, len(species_input)):
            # It's a molecule, so we set to default setting of 100 and 0. The [0] is to take the actual indexes
            if i not in atom_index[0]:
                atom_number[i] = 100
                ion_stage[i] = 0
            else:
                # Set the element name to the first value (as it's an atom ,so only one exists)
                element_name = master_line_list.data[i]['name'][0]
                atomic_value_index = np.where(element_name == self.element_list)
                # If the element doesn't exist in our list we flag it as a bad one with -1
                if atomic_value_index[0].size == 0:
                    atom_number[i] = -1
                    ion_stage[i] = -1
                else:
                    # We take the index value and add one, as indexing starts at 0 but Hydrogen of course has an atomic number of 1
                    atom_number[i] = atomic_value_index[0][0] + 1
                    ion_stage[i] = master_line_list.data[i]['ion']

        return atom_number, ion_stage

    # produces a sepcies list then runs through the other two functiosn to produce an array with all information on
    # each species we find.
    def readline(self, master_line_list):
        # nLines
        # shortening for debugging
        #master_line_list.data = master_line_list.data[10000:100000]
        rows_in_line_list = len(master_line_list.data['lambda'])

        species = []
        # If any isotopes are not.. uh.. 0? Is that hydrogen, or is 0 the lack of any element?
        if (np.any(master_line_list.data['isotope']) > 0):
            # Return the species list
            species = self.isotopic_mix(species, master_line_list, rows_in_line_list)
            # Yeah we could just makei t an array to begin with, but this was kinda cleaner.
            species = np.asarray(species)
        # If we aren't doing that function, we make the species list here, it only takes a few seconds.
        else:
            # For each list of elements involved in the line species, we join them together ([C, N] becomes CN)
            for elements in master_line_list.data['name']:
                species.append("".join(elements))
                species = np.asarray(species)
        print(len(species))
        # Replace CC with C2 denomination.
        species[np.where(species == 'CC')] = 'C2'
        # Finding indexes that contain only one species. Np.asarray required as the fits files are given in chararray which
        # is outdated and only there for backwards compatibility. 'not' means is empty so there are no other species beyond
        # the first. Could put this in the other masterlinelist loop? @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        # We use it to input values in atomnumb array in element parsing.
        atom_index = np.where(np.asarray([not x[1] for x in master_line_list.data['name']]))

        # It seems to think some of the element list might have isotopes, so I'm using an input for if we want to use a list
        # that is variable like that
        atom_number, ion_stage = self.element_parsing(species, atom_index, master_line_list)
        # We just plug the data from the galah file into columns, basically. Necesary input for sme

        data_array = np.transpose(np.array((atom_number, ion_stage, master_line_list.data['lambda'], master_line_list.data['e_low'],
                                             master_line_list.data['log_gf'], master_line_list.data['rad_damp'],
                                             master_line_list.data['stark_damp'], master_line_list.data['vdw_damp'])))
        lande_mean = np.transpose(master_line_list.data['lande_mean'])
        depth = np.transpose(master_line_list.data['depth'])
        lu_lande = np.transpose(np.array((master_line_list.data['lande_low'], master_line_list.data['lande_up'])))
        j_e_array = np.transpose(
            np.array((master_line_list.data['j_low'], master_line_list.data['e_up'], master_line_list.data['j_up'])))
        lower_level = np.transpose(np.array((master_line_list.data['label_low']).replace('\'', '').strip()))
        upper_level = np.transpose(np.array((master_line_list.data['label_up']).replace('\'', '').strip()))
        data_reference_array = np.transpose(
            np.array([master_line_list.data['lambda_ref'], master_line_list.data['e_low_ref'],
                      master_line_list.data['log_gf_ref'], master_line_list.data['rad_damp_ref'],
                      master_line_list.data['stark_damp_ref'],
                      master_line_list.data['vdw_damp_ref'],
                      master_line_list.data['lande_ref']]))
        # Not really totally sure wheret he 0.3 comes from, or why we use data_array at all when we could just call on lambda.
        #wavelength_range = np.array((min(data_array[:, 2]) - 0.3, max(data_array[:, 2]) + 0.3))

        bad_element_index = np.where(atom_number < 0)
        bad_ion_index = np.where(ion_stage < 0)
        if bad_element_index[0].size > 0:
            bad_element_list = species[bad_element_index]
            print("Unkown species:", bad_element_list)
        if bad_ion_index[0].size > 0:
            bad_ion_list = species[bad_ion_index]
            print("Unkown ion charges:", bad_ion_list)
        # And then we just trim off the lande and depth columns from data! what was the point in any of that?! @@@@@@@@@@@@@@
        # Oh but dataref still has the lande ref! @@@@@@@@@@@@@@@
    
        return data_array, lande_mean, depth, data_reference_array, species, j_e_array, lower_level, upper_level, lu_lande




"""at this point there's a lot of looping through each file to add on their own data array etc to the array. It doesn't
make alot of sense if we only have one file so I'm ignoring that (oh and it just calls abuff without defining it. I think
it must be data array).
Also more info on if not short etc. but we always have short set to 2 so that's useless.
Also a lot of transposing and reforming which appears to not do anything for us. So linemerge is a useless file that just
calls readline."""

# Idl sets range = reform(fltarr(2, nfile),2,nfile). Short is always 2 so not sure why we'd call it but.
# line_atomic, line_species, line_lande, line_depth, line_ref,term_low=line_term_low,term_upp=line_term_upp,short_
# format=short_format,extra=line_extra,lulande=line_lulande, line j e is extra, wavelength range is not used


def run_merger():
    # Why is lulande only set for newsme?
    # temporary hard code
    # The reduced amount of spectral lines we need to look at, applied to masterfile
    data_index = np.loadtxt('alllinesindex.csv', delimiter=" ").astype(int)
    #data_index= data_index.astype(int)
    master_line_list = fits.open(line_list)[1]
    #master_line_list.data = master_line_list.data[data_index]
    rows_in_line_list = len(master_line_list.data['lambda'])
    #line_atomic, line_lande, line_depth, line_data_references, line_species, line_j_e_array, lower_level, upper_level, lu_lande = line_merging(master_line_list)
    #data_array = (LineMerging(master_line_list).data_array)
    #print((data_array[:,0]))
    LineMerging(master_line_list)