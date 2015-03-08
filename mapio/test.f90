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
  type(coop_file)::fp
  type(coop_healpix_maps)::hp, hpc, im
  COOP_INT i, list(8), nneigh, nmax, nmin
  call im%read("simu/simu_i_smoothed_fwhm20arcmin.fits")
  call im%convert2nested()  
  call hp%read("planck14/dx11_v2_common_int_mask_010a_1024.fits")  
  call hp%convert2nested()
  hpc = hp
  nmax = 0
  nmin = 0
  do i=0, hp%npix-1
     if(hpc%map(i,1).gt.0.5)then
        call neighbours_nest(hpc%nside, i, list, nneigh)
        if(any(hpc%map(list(1:nneigh), 1) .lt. 0.5))then
           hp%map(i,1) = 0.
        else
           if( im%map(i,1).gt.0. .and. all(im%map(list(1:nneigh), 1) .le. im%map(i, 1)) )then
              nmax = nmax + 1
           elseif( im%map(i,1).lt.0. .and. all(im%map(list(1:nneigh), 1) .ge. im%map(i, 1)))then
              nmin = nmin + 1
           endif
           
        endif
     endif
  enddo

  
  print*, sum(dble(hp%map))/hp%npix*2417.6*coop_4pi
  print*, nmax, nmin
end program test
