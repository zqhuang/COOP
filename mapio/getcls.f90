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

  COOP_STRING::imap_file, polmap_file, output_file
  COOP_INT, parameter::lmax = 2000
  COOP_INT, parameter::smooth_delta_ell = 20
  
  COOP_UNKNOWN_STRING, parameter::imask_file = "planck14/dx11_v2_common_int_mask_010a_1024.fits"
  COOP_UNKNOWN_STRING, parameter::polmask_file = "planck14/dx11_v2_common_pol_mask_010a_1024.fits"
  type(coop_healpix_maps)::map, imask, polmask
  type(coop_file)::fp
  integer l, i
  imap_file = trim(coop_InputArgs(1))
  polmap_file = trim(coop_InputArgs(2))
  if(iargc() .ge. 3)then
     output_file = trim(coop_InputArgs(3))
  else
     output_file = coop_str_replace(imap_file, ".fits", "_cls.txt")
  endif
 ! if(.not. coop_file_exists(output_file))then
     call map%read(imap_file, nmaps_wanted = 3)
     call map%import(polmap_file, index_start = 2, index_end = 3, spin = (/ 2, 2 /))
     call map%map2alm(lmax = lmax)
     call map%get_cls()
     do l=0, map%lmax
        map%cl(l, :) = map%cl(l, :)*(l*(l+1)/coop_2pi*1.e12)
     enddo
     do i=1, 6
        call coop_smooth_data(map%lmax+1, map%cl(0:map%lmax, i), smooth_delta_ell)
     enddo
     call fp%open(output_file, "w")
     do l = 0, map%lmax
        write(fp%unit, "(I5, 6E16.7)") l, map%cl(l,:)
     enddo
     call fp%close()
 ! endif
end program test
