program shells
  use coop_hnn_mod
  implicit none
#include "constants.h"
  COOP_INT,parameter::fsky = 30
  COOP_REAL::fsky_eff
  type(coop_healpix_maps)::map, mask
  COOP_REAL,parameter::cut(3) = (/ 1.d4, 500.d0, 500.d0 /)
  COOP_INT::i

  call map%read("plmaps/DustFilaments_TQU_NS2048_Nfil180p5M_CMBS4BICEP_353p0GHz.fits", nmaps_wanted=3)
!  call map%read("plmaps/HFI_SkyMap_353-psb_2048_R3.01_halfmission-1.fits", nmaps_wanted=3)  
  call mask%read("plmaps/mask_apo5deg_LR80_n2048.fits")
  call map%apply_mask(mask)
  do i=1, map%nmaps
     map%map(:, i) = cut(i)*tanh(map%map(:, i)/cut(i))
  enddo
  call map%smooth(fwhm=4.818d0*coop_SI_arcmin, l_lower=100, l_upper=600, delta_l = 10)
  ! call map%write("plmaps/HFI_353_HM1_L100H600_n2048.fits")
   call map%write("plmaps/DF_L100H600_n2048.fits")
end program shells
