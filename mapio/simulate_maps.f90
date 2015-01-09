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
  COOP_INT,parameter::lmax = 1600
  type(coop_healpix_maps)::map, imask, polmask
  integer l, m, il
  type(coop_file)::fp
  COOP_REAL::beam_fwhm = 30.
  COOP_REAL sigma, w
  call coop_random_init()

  call map%init(nside = 1024, nmaps=3, spin = (/ 0, 2, 2 /))

  call map%allocate_alms(lmax=lmax)
  map%Cl = 0.
  sigma = coop_sigma_by_fwhm * beam_fwhm * coop_SI_arcmin
  
  call fp%open("planckbest_lensedtotCls.dat", "r")
  
  do l=2, map%lmax
     read(fp%unit, *) il, map%cl(l,coop_healpix_index_TT), map%cl(l, coop_healpix_index_EE), map%cl(l, coop_healpix_index_BB), map%cl(l, coop_healpix_index_TE)
     map%cl(l, :) = map%cl(l, :)*exp(-l*(l+1.d0)*sigma**2)
     if(il.ne. l) stop "wrong index"
     map%cl(l,:) = map%cl(l,:)*(coop_2pi/l/(l+1.d0)/1.e12)
  enddo

  do il = 0, 99
     if(coop_file_exists("massive/simu_TQTUT_"//trim(coop_ndigits(il, 5))//"_0"//COOP_STR_OF(nint(beam_fwhm))//"a_n1024.fits"))cycle
     call map%simulate()
     call map%write("massive/simu_int_"//trim(coop_ndigits(il, 5))//"_0"//COOP_STR_OF(nint(beam_fwhm))//"a_n1024.fits", index_list = (/ 1 /) )
     call map%write("massive/simu_pol_"//trim(coop_ndigits(il, 5))//"_0"//COOP_STR_OF(nint(beam_fwhm))//"a_n1024.fits", index_list = (/2, 3/) )

     call map%iqu2TQTUT()
     call map%write("massive/simu_TQTUT_"//trim(coop_ndigits(il, 5))//"_0"//COOP_STR_OF(nint(beam_fwhm))//"a_n1024.fits")
  enddo

end program test
