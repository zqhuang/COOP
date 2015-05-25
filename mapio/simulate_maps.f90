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
  COOP_INT,parameter::nside = 1024
  COOP_REAL,parameter::beam_fwhm = 5.  
  type(coop_healpix_maps)::map
  integer l, m, il
  type(coop_file)::fp
  COOP_UNKNOWN_STRING, parameter::prefix = "simu/simu"
  logical,parameter::do_highpass = .false.
  COOP_REAL:: sigma
  COOP_INT,parameter::llow = 230, lupper=270
  call coop_MPI_Init()
  call coop_random_init()
  if(do_highpass)then
     print*, "Warning: high-pass filter is on"
  endif
  call map%init(nside = 1024, nmaps=1, genre="TEMPERATURE", lmax = lmax)
  map%Cl = 0.
  sigma = coop_sigma_by_fwhm * beam_fwhm * coop_SI_arcmin
  
  call fp%open("planck14best_lensedCls.dat", "r")
  
  do l=2, map%lmax
     read(fp%unit, *) il, map%cl(l,coop_healpix_index_TT), map%cl(l, coop_healpix_index_EE), map%cl(l, coop_healpix_index_BB), map%cl(l, coop_healpix_index_TE)
     map%cl(l, :) = map%cl(l, :)*exp(-l*(l+1.d0)*sigma**2)
     if(il.ne. l) stop "wrong index"
     if(do_highpass)then
        map%cl(l,:) = map%cl(l,:)*(coop_2pi/l/(l+1.d0)*coop_highpass_filter(llow, lupper, l)**2 )        
     else
        map%cl(l,:) = map%cl(l,:)*(coop_2pi/l/(l+1.d0))
     endif
  enddo
  write(*,*) "simulating maps:"
  call map%simulate()
  write(*,*) "maps are simulated"
  if(do_highpass)then
     call map%write(trim(prefix)//"_i_hp_"//COOP_STR_OF(llow)//"_"//COOP_STR_OF(lupper)//"_"//COOP_STR_OF(nside)//"_"//COOP_STR_OF(nint(beam_fwhm))//"a.fits", index_list = (/ 1 /) )                        
   !  call map%write(trim(prefix)//"_pol_hp_"//COOP_STR_OF(llow)//"_"//COOP_STR_OF(lupper)//"_"//COOP_STR_OF(nside)//"_"//COOP_STR_OF(nint(beam_fwhm))//"a.fits", index_list = (/2, 3/) )        
  else
     call map%write(trim(prefix)//"_i_"//COOP_STR_OF(nside)//"_"//COOP_STR_OF(nint(beam_fwhm))//"a.fits", index_list = (/ 1 /) )                        
     !call map%write(trim(prefix)//"_pol_"//COOP_STR_OF(nside)//"_"//COOP_STR_OF(nint(beam_fwhm))//"a.fits", index_list = (/2, 3/) )        
  endif
  call coop_MPI_Finalize()
end program test
