program test
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_fitswrap_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map, mask
  COOP_REAL::fwhm
  COOP_INT::l
  call map%read("lowl/commander_I_n0128_60a.fits")
  call mask%read("lowl/commander_mask_n0128_60a.fits")
  write(*, "(A)") "enter l = "
  read(*,*) l
  fwhm =  asin(0.5d0/(l+0.5d0)/coop_sigma_by_fwhm)*2.d0
  write(*,"(A)") "fwhm = "//COOP_STR_OF(fwhm/coop_SI_degree)//" deg"
  call map%smooth(fwhm = fwhm)
  call map%rotate_coor(l_deg = 212.d0, b_deg = -13.d0)
  call mask%rotate_coor(l_deg = 212.d0, b_deg = -13.d0)
  call map%apply_mask(mask, bad_data = .true.)
  call map%draw_latitude_line(30.d0, 0.5d0)
  call map%draw_latitude_line(60.d0, 0.5d0)
  call map%draw_latitude_line(-30.d0, 0.5d0)
  call map%draw_latitude_line(-60.d0, 0.5d0)
  call map%draw_latitude_line(0.d0, 0.5d0)
  call map%write("map_asym.fits")
  call system("map2gif -inp map_asym.fits -out map_asym_masked_l"//COOP_STR_OF(l)//"smooth.gif -bar T -ttl 'Gaussian smoothing l = "//COOP_STR_OF(l)//" (FWHM "//trim(coop_num2str(fwhm/coop_SI_degree, "(F10.1)"))//"deg)'") 
  
end program test
