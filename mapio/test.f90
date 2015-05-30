program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
#ifdef HAS_HEALPIX
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
  use udgrade_nr
  use coord_v_convert,only:coordsys2euler_zyz
#endif  
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::hmap, lmap, hcpy
  COOP_INT,parameter::lmax = 500
  COOP_REAL::Cls(0:lmax), sqrtCls(0:lmax)
  type(coop_file)::fp
  COOP_INT::l, basenside, i
  Cls(0:1) = 0.d0
  if(iargc().ge.1)then
     basenside = coop_str2int(coop_InputArgs(1))
  else
     print*, "enter nside"
     read(*,*) basenside
  endif
  call fp%open("planck14best_lensedCls.dat", "r") 
  do l=2, lmax
     read(fp%unit, *) i, Cls(l)
     if(i.ne.l) stop "error in cl file"
     Cls(l) = Cls(l)*coop_2pi/l/(l+1.d0)*coop_gaussian_filter(10.d0, l)**2
  enddo
  sqrtCls = sqrt(Cls)
  call fp%close()
  call lmap%init(nside = basenside,  nmaps = 1, genre = "TEMPERATURE",nested = .true.)
  call hmap%init(nside = basenside*8, nmaps = 1, genre = "TEMPERATURE")
  call hcpy%init(nside = basenside*8, nmaps = 1, genre = "TEMPERATURE", nested = .true.)
  call hmap%simulate_Tmaps(hmap%nside, lmax, sqrtCls)
  call hmap%convert2nested()
  call coop_healpix_maps_ave_udgrade(hmap, lmap)
  call coop_healpix_maps_ave_udgrade(lmap, hcpy)
  hmap%map = hmap%map - hcpy%map
  print*, sum(hmap%map*hcpy%map)/sqrt(sum(hmap%map**2)*sum(hcpy%map**2))
end program test
