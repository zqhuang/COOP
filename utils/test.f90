program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT,parameter::lmin = 20, lmax = 30, lmax_mask = 10
  COOP_INT::l, i
  COOP_REAL::Cl(lmin:lmax, 6), Cl_pseudo(lmin:lmax, 6), kernel(lmin:lmax, lmin:lmax, 4), Cl_mask(0:lmax_mask)
  cl_mask(0) = coop_4pi
  call coop_pseudoCl_get_kernel_pol(lmax_mask = lmax_mask, cl_mask = cl_mask, lmin = lmin, lmax = lmax, kernel = kernel)
!!  call coop_pseudoCl_get_kernel(lmax_mask = lmax_mask, cl_mask = cl_mask, lmin = lmin, lmax = lmax, kernel = kernel(:,:,1))
  do l=lmin, lmax
     do i=1, 6
        cl_pseudo(l,i) = 1.d0*i/l/(l+1.d0)
     enddo
  enddo
!!  call coop_pseudocl2cl(lmin = lmin, lmax = lmax, cl_pseudo = cl_pseudo(:,1), cl = cl(:,1), kernel = kernel(:,:,1))
  call coop_pseudocl2cl_pol(lmin = lmin, lmax = lmax, cl_pseudo = cl_pseudo(:,:), cl = cl(:,:), kernel = kernel(:,:,:))
  do l = lmin, lmax
     print*, l, cl(l,1:6)*l*(l+1.d0)
  enddo
end program Test
