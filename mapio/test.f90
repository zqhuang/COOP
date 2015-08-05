program test
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_fitswrap_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::mask, dustmap
  COOP_INT, dimension(:), allocatable::listpix
  COOP_INT::nlist, i, j
  call mask%read("dust/lat30_mask_n1024.fits")
  call dustmap%read("dust/dust_i_n1024_15a.fits")
  nlist = count(mask%map(:,1).gt.0.5)
  allocate(listpix(nlist))
  j = 0
  do i=0, mask%npix-1
     if(mask%map(i, 1) .gt. 0.5)then
        j = j + 1
        listpix(j) = i
     endif
  enddo
  call coop_asy_histogram( x = dble(dustmap%map(listpix, 1)), nbins = 100, xlabel = "$\ln I$", ylabel = "$dP/d\ln I$", filename = "lnI_hist.txt", fit_gaussian = .true.) 
end program test
