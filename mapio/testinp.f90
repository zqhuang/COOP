program test
  use coop_wrapper_utils
  use coop_healpix_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map, mask, m2
  type(coop_healpix_inpaint_space)::inp
  COOP_INT,parameter ::lmax = 5
  COOP_REAL::cls(0:lmax), c0, theta
  COOP_INT::l, ell, i
  type(coop_file)::fp
  call coop_MPI_init()
  call coop_random_init()
  call map%read("planck14/dx11_v2_smica_int_cmb_010a_1024.fits")
  call mask%read("planck14/dx11_v2_smica_int_mask_010a_1024.fits")
  Cls(0:1) = 0.d0
  call fp%open("planck14best_lensedCls.dat", "r")
  c0 = 0.d0
  theta = sqrt(coop_4pi/(12*4**2))/4.d0
  do l=2, lmax
     read(fp%unit, *) ell, Cls(l)
     Cls(l) = Cls(l)/l/(l+1.d0)*coop_2pi
     c0 = c0 + cls(l)*(2.d0*l+1.d0)*exp(-l*(l+1.d0)*theta**2)
  enddo
  print*, c0/coop_4pi
  call inp%init(map, mask, lmax, cls)
  call inp%map%write("inp_mean.fits")
  inp%map%map(inp%unknown,1) = map%bad_data
  call inp%map%write("inpmasked.fits")  
  call coop_MPI_finalize()  
end program test
