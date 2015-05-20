program test
  use coop_wrapper_utils
  use coop_healpix_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map, mask, m2
  type(coop_healpix_inpaint_space)::inp
  COOP_INT,parameter ::lmax = 64
  COOP_REAL::cls(0:lmax)
  COOP_INT::l, ell, i
  type(coop_file)::fp
  call map%read("planck14/dx11_v2_smica_int_cmb_010a_1024.fits")
  call mask%read("planck14/dx11_v2_smica_int_mask_010a_1024.fits")
  Cls(0:1) = 0.d0
  call fp%open("planck14best_lensedCls.dat", "r")
  do l=2, lmax
     read(fp%unit, *) ell, Cls(l)
     Cls(l) = Cls(l)/l/(l+1.d0)*coop_2pi
  enddo
  call inp%init(map, mask, lmax, cls)
  call inp%map%write("inp_mean.fits")
end program test
