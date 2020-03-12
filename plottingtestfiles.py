import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)

#from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
import numpy as np

#image_file= get_pkg_data_filename('tutorials/FITS-images/HorseHead.fits')


#print(fits.info(image_file))

image_data = fits.open(r"C:\Users\jama2357\Documents\Galafiles\Fits-files\benchmarks\1504050009013781.fits", ext = 0)
# Taking the first
print(image_data[0])
print (image_data.info())
print("2")
# Repr adds line breaks between headers making it readable
print(repr(image_data[1].header))
listed = (repr(image_data[1].header))

# the .data values of the file which is the flux
flux=(image_data[0].data)

# Takes the value of CRVAL1 which is the initial angstrom wavelength
starter =  image_data[0].header['CRVAL1']
# The step in wavelength per data point
steps = image_data[0].header['CDELT1']
# Making a list with the interval of steps
wlist = starter + (steps*np.arange(image_data[0].header['NAXIS1']))

print (listed.find("COMB3"))
for line in listed.split("\n"):
    if "COMB3" in line:
        print(line)
plt.figure()
plt.plot(wlist,flux,linewidth=0.1)
plt.xlabel( 'Wavelength (A)')
plt.ylabel( "Flux (Relative)")
plt.show()