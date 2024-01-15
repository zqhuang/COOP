program fsm
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  COOP_STRING::fmap, fout
  COOP_REAL::ramin, ramax, decmin, decmax
  COOP_INT::ileft, iright, itop, ibottom
  type(coop_fits_image_cea)::map
  logical virtual
  COOP_REAL, dimension(:,:),allocatable::map_2d
  if(iargc().lt.10)then
     write(*,*) "Syntax:"
     write(*,*) "./FTrim -map ... -out ... -ramin ... -ramax ... -decmin ... -decmax ... [-virtual F/T]"
     stop
  endif
  call coop_get_command_line_argument(key = "map", arg = fmap)
  call coop_get_command_line_argument(key = "out", arg = fout)
  call map%read(fmap)
  call coop_get_command_line_argument(key = "ramin", arg = ramin)
  call coop_get_command_line_argument(key = "ramax", arg = ramax)
  call coop_get_command_line_argument(key = "decmin", arg = decmin)
  call coop_get_command_line_argument(key = "decmax", arg = decmax)
  call coop_get_command_line_argument(key = "virtual", arg=virtual, default=.false.)
  if(map%dx .gt. 0.d0)then
     ileft = floor((ramin - map%radec_center(1))*coop_SI_degree/map%dx+ map%center(1)) 
     iright = ceiling((ramax - map%radec_center(1))*coop_SI_degree/map%dx+ map%center(1)) 
  else
     iright = ceiling((ramin - map%radec_center(1))*coop_SI_degree/map%dx+ map%center(1)) 
     ileft = floor((ramax - map%radec_center(1))*coop_SI_degree/map%dx+ map%center(1))
  endif

  ileft = max(ileft, 1)
  iright = min(iright, map%nside(1))

  if(map%dy .gt. 0.d0)then
     ibottom = floor((decmin - map%radec_center(2))*coop_SI_degree/map%dy+ map%center(2)) 
     itop = ceiling((decmax - map%radec_center(2))*coop_SI_degree/map%dy+ map%center(2)) 
  else
     itop = ceiling((decmin - map%radec_center(2))*coop_SI_degree/map%dy+ map%center(2)) 
     ibottom = floor((decmax - map%radec_center(2))*coop_SI_degree/map%dy+ map%center(2))
  endif
 
  ibottom = max(ibottom, 1)
  itop = min(itop, map%nside(2))
  map%center(1) = map%center(1) + 1- ileft
  map%center(2) = map%center(2) + 1 - ibottom

  allocate(map_2d(map%nside(1), map%nside(2)))
  call coop_array_copy_real(map%image, map_2d, map%npix)
  if(virtual)then
     map_2d(1:ileft-1,:)=0.
     map_2d(iright+1:map%nside(1),:) = 0.
     map_2d(:, 1:ibottom-1) = 0.
     map_2d(:, itop+1:map%nside(2)) = 0.
     call coop_array_copy_real(map_2d, map%image, map%npix)
  else
     map%nside(1) = iright - ileft + 1
     map%nside(2) = itop - ibottom + 1
     print*, "Input npix = ", map%npix
     map%npix = map%nside(1)*map%nside(2)
     print*, "Output npix = ", map%npix
     deallocate(map%image)
     allocate(map%image(0:map%npix-1))
     call coop_array_copy_real(map_2d(ileft:iright, ibottom:itop), map%image, map%npix)
     call map%header%insert("CRPIX1", COOP_STR_OF(nint(map%center(1))), overwrite=.true.)
     call map%header%insert("CRPIX2", COOP_STR_OF(nint(map%center(2))), overwrite=.true.)
     call map%header%insert("NAXIS1", COOP_STR_OF(map%nside(1)), overwrite=.true.)
     call map%header%insert("NAXIS2", COOP_STR_OF(map%nside(2)), overwrite=.true.)

  endif
  call map%write(fout)
end program fsm
