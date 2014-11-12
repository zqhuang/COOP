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
  type(coop_healpix_maps)::map, imask, polmask
  integer l, m, il
  type(coop_file)::fp
  COOP_REAL::beam_fwhm = 10.*coop_SI_arcmin
  COOP_REAL sigma, w
  call coop_random_init()
  call map%init(nside = 1024, nmaps=3, spin = (/ 0, 2, 2 /))
  call map%allocate_alms(lmax=2000)
  map%Cl = 0.
  sigma = coop_sigma_by_fwhm * beam_fwhm
  call fp%open("camb_cls.dat", "r")
  do l=2, map%lmax
     read(fp%unit, *) il, map%cl(l,coop_healpix_index_TT), map%cl(l, coop_healpix_index_EE), map%cl(l, coop_healpix_index_BB), map%cl(l, coop_healpix_index_TE)
     map%cl(l, :) = map%cl(l, :)/(1.+l*(l+1.d0)*sigma**2)
     if(il.ne. l) stop "wrong index"
     map%cl(l,:) = map%cl(l,:)*(coop_2pi/l/(l+1.d0)/1.e12)
  enddo
  do l = 2, 20
     map%cl(l, coop_healpix_index_EE) = 0.
     map%cl(l, coop_healpix_index_BB)= 0.
     map%cl(l, coop_healpix_index_TE) = 0.
  enddo
  do l = 21, 39
     w = sin(coop_pio2*(l-20.d0)/20.d0)**2
     map%cl(l, coop_healpix_index_EE) =  map%cl(l, coop_healpix_index_EE)*w**2
     map%cl(l, coop_healpix_index_BB)= map%cl(l, coop_healpix_index_BB)*w**2
     map%cl(l, coop_healpix_index_TE) =  map%cl(l, coop_healpix_index_TE)*w    
  enddo
  call map%simulate()
  call map%write("simu/simu_int_010a_n1024.fits", index_list = (/ 1 /) )
  call map%write("simu/simu_pol_hp_20_40_010a_n1024.fits", index_list = (/2, 3/) )

end program test
