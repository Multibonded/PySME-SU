


# So we're using marcs2014.sav, but only to specify that we use specific .grd files for the grid to begin with.
# And the fits file for the line_list

from scipy import io
import numpy as np

# So do we need to read these out properly, or are they fine in binary form for sme? @@@@@@@@@@
# makestruct says if atmogrid_file not set, do 2012 not 2014
atmosphere_grid_file = io.readsav(r"C:\Users\jama2357\Documents\Galafiles\sme_536\atmospheres\marcs2014.sav", idict=None, python_dict=False,
 uncompressed_file_name=None, verbose=False)
# The grid has ~3000 data points on a 3d grid of T, log g, and metalicity which we use to interpolate stars at their
# specific values that are between the points.
hydrogen_grid = open(r'C:\Users\jama2357\Documents\Galafiles\sme_536\NLTE\marcs2012_H2018.grd ', 'rb')
iron_grid = open(r'C:\Users\jama2357\Documents\Galafiles\sme_536\NLTE\marcs2012_Fe2016.grd ', 'rb')
atmosphere_abundance_grid = np.array([hydrogen_grid, iron_grid])
# Tau is the optical depth at some wavelength reference for continuum to determine the deepest point of interest.
# You can plot T against log tau and when x is 0, it is opaque so anything below that we are viewing.
# Rhox is the other option of column mass (accumulated mass) of the atmosphere. These are both pretty bad as
# they cause some kind of spike in the abundance of stars at temperatures but people are working on replacing them
# with a neural network.
atmosphere_depth = "TAU"
atmosphere_interpolation = "RHOX"
# Plane parallel
atmoshpere_geometry = "PP"
