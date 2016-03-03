program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  use coop_fitsio_mod
  implicit none
#include "constants.h"
  COOP_INT,parameter::n1 = 20 , n2 = 10
  COOP_REAL::image(n1, n2)
  call random_number(image)
  call coop_fits_file_write_2d_image(image, "actpost/test.fits")
end program test
