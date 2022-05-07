program test
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_fitswrap_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map, mask, map2, map1
  COOP_REAL::fwhm, Tmax
  COOP_INT::l1, l2
  call map%read("lowl/commander_I_n0128_60a.fits")
  call mask%read("lowl/commander_mask_n0128_60a.fits")
  map1 = map
  map2 = map
  write(*, "(A)") "enter l1, l2 = "
  read(*,*) l1, l2
  fwhm =  1.d0/(l1+0.5d0)/coop_sigma_by_fwhm
  call map1%smooth(fwhm = fwhm)
  fwhm =  1.d0/(l2+0.5d0)/coop_sigma_by_fwhm
  call map2%smooth(fwhm = fwhm)
!  write(*,"(A)") "fwhm = "//COOP_STR_OF(fwhm/coop_SI_degree)//" deg"
!  call map%smooth(fwhm = fwhm)
  call map1%rotate_coor(l_deg = 212.d0, b_deg = -13.d0)
  call map2%rotate_coor(l_deg = 212.d0, b_deg = -13.d0)  
  call mask%rotate_coor(l_deg = 212.d0, b_deg = -13.d0)
  call map%convert2ring()
  call map1%convert2ring()
  call map2%convert2ring()  
  map%map = map2%map - map1%map
  call map%apply_mask(mask, bad_data = .true.)
  call map%draw_latitude_line(30.d0, 0.5d0)
  call map%draw_latitude_line(60.d0, 0.5d0)
  call map%draw_latitude_line(-30.d0, 0.5d0)
  call map%draw_latitude_line(-60.d0, 0.5d0)
  call map%draw_latitude_line(0.d0, 0.5d0)
  call map%write("map_asym.fits")
  Tmax = maxval(map%map, mask = mask%map .lt. 0.5)
  call system( "map2gif -inp map_asym.fits -out map_asym_masked_diff_l"//COOP_STR_OF(l1)//"_l"//COOP_STR_OF(l2)//".gif -bar T -ttl 'Difference map between l_smooth = "//COOP_STR_OF(l2)//" and l_smooth ="//COOP_STR_OF(l1)//"'" ) 
  
end program test
