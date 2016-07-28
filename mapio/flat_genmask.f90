program GMask
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  type(coop_fits_image_cea)::map, map2
  COOP_STRING::fmap, fout
  COOP_REAL::cut, rms
  COOP_INT::i, nn
  if(iargc() .lt. 2)then
     write(*,*) "Syntax: ./FGenMask -map1 ... -cut1 ... -map2 ... -cut2 ... -map3 ... -cut3 ... -out ..."
     write(*,*) "cuts are in the reange 0-1 (*rms), default 0.2"
     stop
  endif
  call coop_get_command_line_argument(key = "map1", arg = fmap)
  call coop_get_command_line_argument(key = "cut1", arg = cut, default = 0.2d0)
  
  call map%open(fmap)
  map2 = map
  nn = count( map2%image .gt. 0.d0)
  if(nn .eq. 0) stop "Error: map = 0"
  rms = sqrt(sum(map2%image**2, mask = map2%image .gt. 0.d0)/nn)
  where(map2%image .gt. rms*cut)
     map%image = 1.d0
  elsewhere
     map%image = 0.d0
  end where
  i = 2
  do 
     call coop_get_command_line_argument(key = "map"//COOP_STR_OF(i), arg = fmap, default = "")
     if(trim(fmap).eq."")exit
     call coop_get_command_line_argument(key = "cut"//COOP_STR_OF(i), arg = cut, default = 0.2d0)     
     call map2%open(fmap)
     if(map2%npix .ne. map%npix) stop "Error: the maps are of different sizes"
     nn = count( map2%image .gt. 0.d0)
     if(nn .eq. 0) stop "Error: map = 0"
     rms = sqrt(sum(map2%image**2, mask = map2%image .gt. 0.d0)/nn)
     where(map2%image .lt. rms*cut)
        map%image = 0.d0
     end where
     i = i + 1
  enddo
  call coop_get_command_line_argument(key = "out", arg = fout)
  call map%write(fout)

end program GMask
