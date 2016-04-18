program shells
  use coop_wrapper_firstorder
  use coop_zeta3d_mod
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::base_imap, imap, base_iqu_hm1, base_iqu_hm2, iqu_hm1, iqu_hm2
  COOP_INT,dimension(:),allocatable::listpix, ind_start
  COOP_INT::npix, i, j, ngroups
  COOP_SINGLE::cut = 50.
  call base_imap%read("dust/dust_i_10a_hp"
  npix = count(map%map(:,1).ge. cut)
  allocate(listpix(npix), ind_start(npix))
  call map%convert2nested()
  j = 0
  do i = 0, map%npix-1
     if(map%map(i,1).ge. cut)then
        j = j + 1
        listpix(j) = i
     endif
  enddo
  call coop_healpix_group_connected_pixels(nside = map%nside, listpix = listpix, ngroups = ngroups, ind_start = ind_start)
  call mask%init(nside = map%nside, nmaps = 1, genre = "MASK", nested = .true.)
  mask%map = 0.
  print*, "ngroups = ", ngroups
  do i =  1, min(ngroups, 200)
     print*, "group ", i
     mask%map( listpix(ind_start(i):ind_start(i+1)-1),1) = 1.
     call mask%write("mask.fits")
     call system("map2gif -inp mask.fits -out mask"//COOP_STR_OF(i)//".gif -bar T")
  enddo


end program shells
