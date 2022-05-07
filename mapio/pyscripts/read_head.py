# Load the WCS information from a fits header, and use it
# to convert pixel coordinates to world coordinates.

from __future__ import division, print_function
import math
import numpy
from astropy import wcs
from astropy.io import fits
import sys

def load_wcs_from_file(filename):
    # Load the FITS hdulist using astropy.io.fits
    hdulist = fits.open(filename, memmap=True)
    hdulist.info()
#    hdulist[0].header['CTYPE1']
#    hdulist.writeto('backup.fits')
    # Parse the WCS keywords in the primary HDU
    w = wcs.WCS(hdulist[0].header)
    # Print out the "name" of the WCS, as defined in the FITS header
    #print(w.wcs.name)
    # Print out all of the settings that were parsed from the header
    #w.wcs.print_contents()

    # Some pixel coordinates of interest.
    pixcrd = numpy.array([[0., 0.], [1., 1.]])
    
    # Convert pixel coordinates to world coordinates
    # The second argument is "origin" -- in this case we're declaring we
    # have 1-based (Fortran-like) coordinates.
    world = w.wcs_pix2world(pixcrd, 1)
    print(convert(pixcrd))
    print(world )
    # Convert the same coordinates back to pixel coordinates.
#    pixcrd2 = w.wcs_world2pix(world, 1)
#    print(pixcrd2)

    # These should be the same as the original pixel coordinates, modulo
    # some floating-point error.
  #  assert numpy.max(numpy.abs(pixcrd - pixcrd2)) < 1e-6

def convert(list):
    for x in list:
        x[0] = (x[0]-1030.)*(-0.008333333333333333333)
        if(x[0] < 0):
            x[0] = x[0] + 360
        x[1] = (x[1]-808.)*(0.00833333333333333333333)
        x[1] = math.asin(x[1]*(math.pi/180.))*180./math.pi
    return list

hdu=fits.PrimaryHDU()
hdu.header


if __name__ == '__main__':
    load_wcs_from_file(sys.argv[-1])
