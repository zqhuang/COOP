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
  COOP_INT,parameter::lmax = 500
  COOP_INT,parameter::nside = 16
  COOP_REAL,parameter::beam_fwhm = 440.  
  type(coop_healpix_maps)::map, mcopy
  integer l, m, il, i
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
  call map%init(nside = nside, nmaps=1, genre="TEMPERATURE", lmax = lmax)
  map%Cl = 0.
  sigma = coop_sigma_by_fwhm * beam_fwhm * coop_SI_arcmin
  
  call fp%open("planck14best_lensedCls.dat", "r")
  
  do l=2, map%lmax
     read(fp%unit, *) il, map%cl(l,coop_healpix_index_TT)
     map%cl(l, :) = map%cl(l, :)*exp(-l*(l+1.d0)*sigma**2)
     if(il.ne. l) stop "wrong index"
     if(do_highpass)then
        map%cl(l,:) = map%cl(l,:)*(coop_2pi/l/(l+1.d0)*coop_highpass_filter(llow, lupper, l)**2 )        
     else
        map%cl(l,:) = map%cl(l,:)*(coop_2pi/l/(l+1.d0))
     endif
  enddo
  mcopy = map
  do i=1, 10
     write(*,*) "simulating maps", i     
     call map%simulate()
     call map%write(trim(prefix)//"_i_"//COOP_STR_OF(nside)//"_"//COOP_STR_OF(nint(beam_fwhm))//"a_"//COOP_STR_OF(i)//".fits", index_list = (/ 1 /) )
     map%cl = mcopy%cl
  enddo

  call coop_MPI_Finalize()
end program test
