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
#endif  
  implicit none
#include "constants.h"
  COOP_INT,parameter::lmax = 2000
  COOP_INT::i, j, l, nside, il, ih
  COOP_REAL::Cls(0:lmax), sqrtCls(0:lmax), theta(8), ells(0:lmax)
  type(coop_healpix_maps)::lmap, hmap
  type(coop_file)::fp
  Cls(0:1) = 0.d0
  call fp%open("planck14best_lensedCls.dat", "r")
  ells(1) = 1.
  ells(0) = 0.
  do l=2, lmax
     read(fp%unit, *) i, Cls(l)
     ells(l) = l
     if(i.ne.l) stop "error in Cl file"
     Cls(l) = Cls(l)*coop_2pi/l/(l+1.d0)*coop_gaussian_filter(10.d0, l)**2
  enddo  
  call fp%close()
  sqrtCls = sqrt(Cls)  
  call hmap%simulate_Tmaps(nside = 2048, lmax = lmax, sqrtCls = sqrtCls)
  call hmap%map2alm()
  do l=2, lmax
     hmap%alm(l, :, :) = hmap%alm(l, :, :)*sqrt(Cls(l)/hmap%cl(l,1))
  enddo
  call hmap%alm2map()
  do i=8, 1, -1
     call lmap%init(nside = 2**i, nmaps = 1, genre="T")
     call coop_healpix_maps_ave_udgrade(hmap, lmap)
     call coop_healpix_maps_ave_udgrade(lmap, hmap)     
     call hmap%map2alm()
     il = 2
     ih = lmap%nside+2
     theta(i) = sqrt( sum(-log(hmap%cl(il:ih,1)/cls(il:ih))/ells(il:ih)/(ells(il:ih)+1.d0))/(ih-il+1.d0) )
     write(*,*) 2**i, theta(i)*lmap%nside
  enddo
  
end program test
