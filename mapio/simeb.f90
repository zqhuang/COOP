program simeb
  use coop_wrapper_utils
  use coop_healpix_mod
  implicit none
#include "constants.h"
  !!working directory, put all maps there
  character(LEN=*),parameter::dir = "inps/"
  !!the input IQU maps
  character(LEN=*),parameter::map_file = "predx11_iqu_nside256.fits"
  !!temperature mask
  character(LEN=*),parameter::imask_file = "predx11_imask_nside256.fits"
  !!Polarization mask
  character(LEN=*),parameter::polmask_file = "predx11_polmask_nside256.fits"

  !!mask smoothing scale
  real, parameter:: mask_smooth_scale = 2. * coop_SI_degree
  
  call coop_healpix_inpainting(mode = "IQU", map_file = dir//map_file, mask_file = dir//imask_file, maskpol_file = dir//polmask_file, output_freq = 10, output_types= "TEB")

  

end program simeb
