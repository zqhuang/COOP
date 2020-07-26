program test
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_fitswrap_mod
  use coop_stacking_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map, mask
  type(coop_function)::mapping
  COOP_INT::pix
  call coop_MPI_init()
  call map%read("dust/dust_i_n1024_30a.fits")  
  call mask%read("dust/lat30_mask_n1024.fits")
  call map%Gaussianize(mask, 1, mapping)
  call coop_asy_plot_function(mapping, "mapping.txt", xlabel = "$\ln I$", ylabel = "$x_g$")
  call map%apply_mask(mask, bad_data = .true.)
  call map%write("dust/dust_xg_n1024_30a.fits")
  call coop_MPI_Finalize()
end program test
