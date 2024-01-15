program main
  use coop_hnn_mod
  use tmp
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map1, map2, mask
  call map1%read("sims/noBB_conv_TEB_submap002.fits", nmaps_wanted=1)
  call map2%read("sims/Qr_noBB.fits", nmaps_wanted=1)  
  call mask%read("sims/mask256.fits", nmaps_wanted=1)  
  print*, sum(map1%map(:,1)*map2%map(:, 1)*mask%map(:, 1))/sqrt(sum( (map1%map(:,1)*mask%map(:,1))**2))/ sqrt(sum( (map2%map(:,1)*mask%map(:,1))**2 ))
end program main
