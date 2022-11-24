program test
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_fitswrap_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map
  COOP_INT,parameter::lmax = 200
  type(coop_dynamic_array_real)::mydata, model
  COOP_INT::l
  call map%read("lowl/commander_dx11d2_extdata_temp_cmb_n0256_60arc_v1_cr.fits")
  call mydata%load_txt("cls_commander_2014_rc_TT.dat")
  call model%load_txt("planck14best_lensedCls.dat")
  call map%map2alm(lmax = lmax)
  do l = 2, lmax/2
     write(*,*) l, map%cl(l, 1)*l*(l+1.d0)/coop_2pi, mydata%f(l-1, 4)*coop_Gaussian_filter(60.d0, l)**2, model%f(l-1, 2)*coop_Gaussian_filter(60.d0, l)**2
  enddo
  
end program test
