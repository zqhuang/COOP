program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools

  implicit none
#include "constants.h"

  COOP_STRING::map_file, output_file
  COOP_UNKNOWN_STRING, parameter::imask_file = "planck14/dx11_v2_common_int_mask_010a_1024.fits"
  COOP_UNKNOWN_STRING, parameter::polmask_file = "planck14/dx11_v2_common_pol_mask_010a_1024.fits"
  type(coop_healpix_maps)::map, imask, polmask
  type(coop_file)::fp
  integer l
  map_file = trim(coop_InputArgs(1))
  if(trim(map_file).eq."")then
     write(*,*) "Syntax:"
     write(*,*) "./GetCl map [output]"
     stop 
  endif
  if(iargc() .ge. 2)then
     output_file = trim(coop_InputArgs(2))
  else
     output_file = coop_str_replace(map_file, ".fits", "_cls.txt")
  endif
  call map%read(map_file)
  select case(map%nmaps)
  case(1)
     if(imask_file .ne. "")then
        call imask%read(imask_file, spin = (/ 0 /) )
        map%map(:, 1) = map%map(:, 1) * imask%map(:, 1)
     endif
  case(2)
     if(polmask_file .ne. "")then
        call polmask%read(polmask_file, spin = (/ 0 /) )
        map%map(:, 1) = map%map(:, 1) * polmask%map(:, 1)
        map%map(:, 2) = map%map(:, 2) * polmask%map(:, 1)
     endif
  case(3)
     if(imask_file .ne. "")then
        call imask%read(imask_file, spin = (/ 0 /))
        map%map(:, 1) = map%map(:, 1) * imask%map(:, 1)
     endif
     if(polmask_file .ne. "")then
        call polmask%read(polmask_file, spin = (/ 0 /))
        map%map(:, 2) = map%map(:, 2) * polmask%map(:, 1)
        map%map(:, 3) = map%map(:, 3) * polmask%map(:, 1)
     endif     
  end select
  call fp%open(output_file)
  call map%map2alm()
  do l = 0, map%lmax
     write(fp%unit, "(I5, 6E16.7)") l, map%cl(l, :)
  enddo
  call fp%close

end program test
