program simeb
  use coop_wrapper_utils
  use coop_healpix_mod
  implicit none
#include "constants.h"
  !!working directory, put all maps there
  character(LEN=*),parameter::dir = "inps/"
  !!the input IQU maps
  character(LEN=*),parameter::map_file = "sim2_iqu_nside512.fits"
  !!temperature mask
  character(LEN=*),parameter::imask_file = "predx11_imask_nside512.fits"
  !!Polarization mask
  character(LEN=*),parameter::polmask_file = "predx11_polmask_nside512.fits"

  !!mask smoothing scale
  real*8, parameter:: mask_smooth_scale = 2.5 * coop_SI_degree
  
  call coop_healpix_inpainting(mode = "IQU", map_file = dir//map_file, mask_file = dir//imask_file, maskpol_file = dir//polmask_file, output_freq = 20, output_types= "IQU", mask_smooth_scale = mask_smooth_scale)

  

end program simeb
