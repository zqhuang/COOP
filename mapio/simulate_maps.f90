program simmaps
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
  COOP_REAL::fwhm_arcmin
  COOP_INT::nside, lmax, nmaps
  type(coop_healpix_maps)::map
  type(coop_file)::fp
  COOP_STRING::fcl, fmap
  COOP_REAL:: Cls(lmin:lmax, 6)
  COOP_INT::l, il
  call coop_random_init()
  call coop_get_command_line_argument(key = "cl", arg = fcl)
  call coop_get_command_line_argument(key = "map", arg = fmap)
  call coop_get_command_line_argument(key = "FWHM", arg = fwhm_arcmin )
  call fp%open_skip_comments(fcl)
  Cls = 0.d0
  do l=lmin, lmax
     read(fp%unit, *) il, Cls(l,coop_TEB_index_TT), Cls(l,coop_TEB_index_EE), Cls(l,coop_TEB_index_BB), Cls(l,coop_TEB_index_TE)
     if(il.ne. l) stop "Cl file error"     
     Cls(l,:) = Cls(l,:)*(coop_2pi/l/(l+1.d0))*coop_Gaussian_filter(fwhm_arcmin = fwhm_arcmin, l = l)**2
     if(Cls(l,coop_TEB_index_TE)**2 .gt. Cls(l,coop_TEB_index_TT)*Cls(l,coop_TEB_index_EE))then
        write(*,*) "l = ", l, " TE^2 > TT * EE"
        stop
     endif
  enddo
  call fp%close()
  call map%init(nside = nside, nmaps=3, genre="TEB", lmax = lmax)
  map%Cl(lmin:lmax, 1:6) = Cls(lmin:lmax,1:6)   
  call map%simulate()
  call map%write(trim(prefix)//"_TEB.fits")
end program simmaps
