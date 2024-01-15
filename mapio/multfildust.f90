program fd
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
  type(coop_healpix_maps)::gmap, mask, map, supermask
  COOP_SINGLE::nu = 1.
  COOP_SINGLE::threshold, summask, mean, rms, meanicut
  type(coop_function)::mapping, invmapping
  COOP_INT::lcut, i
  call map%read("dust/dust_i_7a_n1024.fits")
  gmap = map
  call gmap%gaussianize(imap = 1, mapping = mapping, invmapping = invmapping)
  call gmap%map2alm()


  call mask%read("planck15/mask_lat30_n1024.fits")
  summask = sum(mask%map(:,1))
  supermask = mask


!!$  call gmap%convert2nested()
!!$  call mask%convert2nested()

!!$  mean = sum(gmap%map(:,1)*mask%map(:,1))/summask
!!$  rms = sum((gmap%map(:,1)-mean)*mask%map(:,1))/summask
!!$  meanicut = invmapping%eval(dble(mean + rms*nu))
!!$  gmap%map = map%map  !!gmap will still be used to check the mean luminosity

  do lcut = 500, 50, -50
     print*, lcut
     call band_filter(lcut-40, lcut+40)
     call supermask%write("dust/mask_l"//COOP_STR_OF(lcut)//".fits")
     call system("map2gif -inp dust/mask_l"//COOP_STR_OF(lcut)//".fits -out mask"//COOP_STR_OF(lcut)//".gif -bar T")
  enddo


contains

  subroutine band_filter(l1, l2)
    COOP_INT,intent(IN)::l1, l2
    COOP_INT::l
    COOP_SINGLE::rms, mean
    COOP_INT,dimension(:),allocatable::listpix, ind_start
    COOP_INT::npix, i, j, ngroups, min_group_size
    call map%allocate_alms(min(l2, gmap%lmax))
    map%alm = 0.
    do l = l1, l2
       map%alm(l,0:l, 1) = gmap%alm(l, 0:l, 1)*sin((l-l1)/dble(l2-l1)*coop_pi)**2
    enddo
    call map%alm2map()
    call map%convert2nested()
    call mask%convert2nested()
    mean = sum(map%map(:,1)*mask%map(:,1))/summask
    rms = sqrt ( sum( (map%map(:,1)-mean)**2*mask%map(:,1))/summask)
    threshold = invmapping%eval(dble(mean + rms * nu))
    !$omp parallel do
    do i = 0, map%npix-1
       map%map(i, 1) = invmapping%eval(dble(map%map(i, 1)))
    enddo
    !$omp end parallel do

    npix = count(map%map(:,1).ge. threshold .and. mask%map(:,1).ge.0.5 )
    allocate(listpix(npix), ind_start(npix))
    j = 0
    do i = 0, map%npix-1
       if(map%map(i,1).ge. threshold .and. mask%map(i,1).ge.0.5)then
          j = j + 1
          listpix(j) = i
       endif
    enddo
    call supermask%convert2nested()
    call coop_healpix_group_connected_pixels(nside = map%nside, listpix = listpix, ngroups = ngroups, ind_start = ind_start)
    min_group_size = 400
    do  i = 1, ngroups
!!$       !!do some selection here
       if(ind_start(i+1)- ind_start(i)  .ge. min_group_size)then
          supermask%map(listpix(ind_start(i):ind_start(i+1)-1),1) = 0.
       endif
    enddo
    deallocate(listpix, ind_start)
  end subroutine band_filter

end program fd
