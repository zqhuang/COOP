module coop_constants
  implicit none
#include "constants.h"
!! universal
  integer,parameter::coop_integer_length = kind(1)  
  integer,parameter::coop_real_length = kind(1.d0)  ! double precision
  integer,parameter::coop_string_length = 256
  integer,parameter::coop_short_string_length = 32
  integer,parameter::coop_long_string_length = 8192
  integer,parameter::coop_default_array_size = 1024  ! default array size for interpolation

  COOP_REAL, parameter:: coop_pi = asin(1.d0)*2.d0
  COOP_STRING, parameter:: coop_version = "0.0"

!!maximum number of species
  COOP_INT, parameter::coop_max_num_species = 16

!!for dark energy background/background.f90
  COOP_INT, parameter:: coop_max_num_DE_params  = 16
 
!!for primordial power spectra firstorder/primordialpower.f90
  COOP_INT, parameter:: coop_max_num_PP_params = 32

 


end module coop_constants
