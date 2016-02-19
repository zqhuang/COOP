program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::hm, mask
  type(coop_fits_image_cea)::fm
  logical::pol = .true.
  call fm%read("act16/deep56_coadd_I.fits")
  if(pol)then
     call hm%read("planck14/COM_CMB_IQU-smica_1024_R2.02_full.fits", nmaps_wanted = 3)
     call mask%read("planck14/COM_CMB_IQU-smica_1024_R2.02_full.fits", nmaps_wanted = 5)  
     hm%map(:,1) = mask%map(:,5)
     call mask%free()
     call hm%rotate_coor(from = "G", to = "C")
     call fm%from_healpix(hm, 2)
     call fm%write("act16/planck_smica_Q.fits")
     call fm%from_healpix(hm, 3)
     call fm%write("act16/planck_smica_U.fits")
     call fm%from_healpix(hm, 1)
     call fm%write("act16/planck_smica_polmask.fits")
  else
     call hm%read("planck14/COM_CMB_IQU-smica-field-Int_2048_R2.01_full.fits", nmaps_wanted = 1)
     call mask%read("planck14/COM_CMB_IQU-smica-field-Int_2048_R2.01_full.fits", nmaps_wanted = 2)
     call hm%rotate_coor(from = "G", to = "C")
     call fm%from_healpix(hm, 1)
     call fm%write("act16/planck_smica_I.fits")
     call hm%read("planck14/COM_CMB_IQU-smica-field-Int_2048_R2.01_full.fits", nmaps_wanted = 1)
     hm%map(:,1) = mask%map(:,2)
     call mask%free()
     call hm%rotate_coor(from = "G", to = "C")
     call fm%from_healpix(hm, 1)
     call fm%write("act16/planck_smica_imask.fits")
  endif
end program test
