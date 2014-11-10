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
  COOP_REAL sigma
  call coop_random_init()
  call map%init(nside = 1024, nmaps=3, spin = (/ 0, 2, 2 /))
  call map%allocate_alms(lmax=2000)
  sigma = coop_sigma_by_fwhm * beam_fwhm
  call fp%open("planck14/planck14_smica_cls.txt", "r")
  do l=0, map%lmax
     read(fp%unit, *) il, map%cl(l,:)
     !!coop_healpix_index_TT), map%cl(l, coop_healpix_index_EE), map%cl(l, coop_healpix_index_BB), map%cl(l, coop_healpix_index_TE)
     !!map%cl(l, :) = map%cl(l, :)*exp(-l*(l+1.d0)*sigma**2)
     if(il.ne. l) stop "wrong index"
     if(l.ge. 1)then
        map%cl(l,:) = map%cl(l,:)*(coop_2pi/1.e12/l/(l+1.d0))
     endif
  enddo
  map%cl(0:1,:) = 0.
  
  do il = 1, 10
     call map%simulate()
     !  call map%map2alm()
     write(*,*) il
     call map%write("simu/cmb/int/dx11_v2_smica_int_cmb_mc_"//trim(coop_Ndigits(il-1, 5))//"_010a_1024.fits", index_list = (/ 1 /) )
     call map%write("simu/cmb/pol/dx11_v2_smica_pol_case1_cmb_mc_"//trim(coop_Ndigits(il-1, 5))//"_hp_20_40_010a_1024.fits", index_list = (/2, 3/) )
  enddo

end program test
