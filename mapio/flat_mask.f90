program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  type(coop_fits_image_cea)::map, mask
  COOP_STRING::fmap, fmask, fout
  COOP_REAL threshold
  if(iargc() .lt. 4)then
     write(*,*) "Syntax: ./FMask MAP MASK THRESHOLD OUTPUT"
     stop
  endif
  call coop_get_Input(1, fmap)
  call coop_get_Input(2, fmask)
  call coop_get_Input(3, threshold)
  call coop_get_Input(4, out)
  call map%open(fmap)
  call mask%open(fmask)
  if(map%npix .ne. mask%npix)then
     write(*,*) "mask and map are not of the same size"
     stop
  endif
  where(mask%image .lt. threshold)
     map%image = 0.
  end where
  call map%write(fout)

end program test
