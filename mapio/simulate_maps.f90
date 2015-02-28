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
  COOP_INT,parameter::lmax = 1500
  type(coop_healpix_maps)::map, imask, polmask, mapcopy
  integer l, m, il
  type(coop_file)::fp
  COOP_UNKNOWN_STRING, parameter::prefix = "simu/simu"
  logical,parameter::do_highpass = .true.
  COOP_REAL::beam_fwhm = 5.
  COOP_INT,parameter::llow = 230, lupper=270
  COOP_REAL sigma, w
  call coop_random_init()
  if(do_highpass)then
     print*, "Warning: high-pass filter is on"
  endif
  call map%init(nside = 2048, nmaps=3, spin = (/ 0, 2, 2 /))

  call map%allocate_alms(lmax=lmax)
  map%Cl = 0.
  sigma = coop_sigma_by_fwhm * beam_fwhm * coop_SI_arcmin
  
  call fp%open("planck14best_lensedCls.dat", "r")
  
  do l=2, map%lmax
     read(fp%unit, *) il, map%cl(l,coop_healpix_index_TT), map%cl(l, coop_healpix_index_EE), map%cl(l, coop_healpix_index_BB), map%cl(l, coop_healpix_index_TE)
     map%cl(l, :) = map%cl(l, :)*exp(-l*(l+1.d0)*sigma**2)
     if(il.ne. l) stop "wrong index"
     if(do_highpass)then
        map%cl(l,:) = map%cl(l,:)*(coop_2pi/l/(l+1.d0)/1.e12*coop_highpass_filter(llow, lupper, l)**2 )        
     else
        map%cl(l,:) = map%cl(l,:)*(coop_2pi/l/(l+1.d0)/1.e12)
     endif
  enddo
  mapcopy = map
!  do il = 0, 0
!     print*, "simulating map #", il
!     if(coop_file_exists(trim(prefix)//"_TQTUT_"//trim(coop_ndigits(il, 5))//"_0"//COOP_STR_OF(nint(beam_fwhm))//"a_0512.fits"))cycle
     call map%simulate()
     if(do_highpass)then
        call map%write(trim(prefix)//"_i_hp_"//COOP_STR_OF(llow)//"_"//COOP_STR_OF(lupper)//"_smoothed_fwhm"//COOP_STR_OF(nint(beam_fwhm))//"arcmin.fits", index_list = (/ 1 /) )                
        call map%write(trim(prefix)//"_pol_hp_"//COOP_STR_OF(llow)//"_"//COOP_STR_OF(lupper)//"_smoothed_fwhm"//COOP_STR_OF(nint(beam_fwhm))//"arcmin.fits", index_list = (/2, 3/) )        
     else
        call map%write(trim(prefix)//"_i_hp_"//COOP_STR_OF(llow)//"_"//COOP_STR_OF(lupper)//"_smoothed_fwhm"//COOP_STR_OF(nint(beam_fwhm))//"arcmin.fits", index_list = (/ 1 /) )                        
        call map%write(trim(prefix)//"_i_smoothed_fwhm"//COOP_STR_OF(nint(beam_fwhm))//"arcmin.fits", index_list = (/ 1 /) )                
     endif

!     map%cl = mapcopy%cl
!  enddo

end program test
