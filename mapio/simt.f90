program simt
  use coop_wrapper_utils
  use coop_healpix_mod
  implicit none
#include "constants.h"
  !!working directory, put all maps there
  character(LEN=*),parameter::dir = "inps/"
  !!the input IQU maps
  character(LEN=*),parameter::map_file = "predx11_iqu_nside256_submap001.fits"
  !!temperature mask
  character(LEN=*),parameter::imask_file = "predx11_imask_nside256.fits"
  !!Polarization mask
  !!mask smoothing scale
  real*8, parameter:: mask_smooth_scale = 2. * coop_SI_degree
  
  call coop_healpix_inpainting(mode = "I", map_file = dir//map_file, mask_file = dir//imask_file, output_freq = 20, mask_smooth_scale = mask_smooth_scale)

  

end program simt
