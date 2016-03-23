program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT,parameter::lmin = 20, lmax = 1200, lmax_mask = 50
  COOP_REAL,dimension(:,:),allocatable::kernel
  COOP_REAL::Cl_pseudo(lmin:lmax), Cl(lmin:lmax)
  COOP_REAL::Cl_mask(0:lmax_mask)
  COOP_INT::l
  allocate(kernel(lmin:lmax, lmin:lmax))
  do l = 0, lmax_mask
     Cl_mask(l) = (4.d0*coop_pi)/(l+1.d0)**5
  enddo

  do l = lmin, lmax
     Cl_pseudo(l) = 1.d0/(l+1.d0)**2
  enddo
  call coop_prtsystime(.true.)
  call coop_pseudoCl_get_kernel(lmin = lmin, lmax = lmax, lmax_mask = lmax_mask, Cl_mask = Cl_mask, kernel = kernel)
  call coop_prtsystime()
  call coop_pseudoCl2Cl(lmin = lmin, lmax = lmax, kernel = kernel, Cl_pseudo = Cl_pseudo, Cl = Cl)
  call coop_prtsystime()
  do l = lmin, lmax
     print*, l, Cl(l)/Cl_pseudo(l)
  enddo
end program Test
