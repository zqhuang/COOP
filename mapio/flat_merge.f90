program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  character(LEN=*),parameter::mapdir = "act/"
  type(coop_fits_image_cea)::cf(0:3), mask(0:3)
  integer i
  do i=0,3
     call cf(i)%open(mapdir//"U5"//trim(coop_num2str(i))//".fits")
     call mask(i)%open(mapdir//"W_U5"//trim(coop_num2str(i))//".fits")
     cf(i)%image = cf(i)%image*mask(i)%image
  enddo
  do i=1,3
     mask(0)%image = mask(0)%image + mask(i)%image
     cf(0)%image = cf(0)%image + cf(i)%image
  enddo

  where(mask(0)%image .gt. 0.)
     cf(0)%image = cf(0)%image/mask(0)%image
  elsewhere
     cf(0)%image = 0.d0
  end where
  call cf(0)%write(mapdir//"U5.fits")
  call mask(0)%write(mapdir//"WU5.fits")
end program test
