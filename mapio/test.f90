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
  COOP_INT,parameter::lmax = 1000
  COOP_INT,parameter::nside = 512
  COOP_REAL,parameter::beam_fwhm = 20.
  type(coop_healpix_maps)::map
  integer l, m, il, i
  type(coop_file)::fp
  COOP_UNKNOWN_STRING, parameter::prefix = "simmaps/sim_zeta"
  COOP_REAL:: sigma, Cls(0:lmax), junk
  COOP_STRING::line
  call coop_MPI_Init()  
  call coop_random_init()
  !!load Cls
  Cls = 0.d0
  sigma = coop_sigma_by_fwhm * beam_fwhm * coop_SI_arcmin  
  call fp%open("../firstorder/output/lcdm_zetaCls.dat", "r")
  read(fp%unit, '(a)') line
  do l=2, lmax
     read(fp%unit, *) il, junk, junk, Cls(l)
     if(il.ne. l) stop "Cl file error"     
     Cls(l) = Cls(l)*(coop_2pi/l/(l+1.d0))*coop_gaussian_filter(beam_fwhm, l)**2/2.726**2/100.
  enddo
  call fp%close()
  call map%init(nside = nside, nmaps=1, genre="ZETA", lmax = lmax)
  map%Cl = 0.
  map%Cl(0:lmax,1) = Cls
  write(*,*) size(map%Cl), map%lmax  
  call map%simulate()
  call map%write(trim(prefix)//"_"//COOP_STR_OF(map%nside)//"_"//COOP_STR_OF(nint(beam_fwhm))//"a.fits")
end program test
