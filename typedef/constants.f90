module coop_constants
  implicit none
#include "constants.h"
!! universal
  integer,parameter::coop_integer_length = kind(1)  
  integer,parameter::coop_real_length = kind(1.d0)  ! double precision
  integer,parameter::coop_string_length = 256
  integer,parameter::coop_short_string_length = 32
  integer,parameter::coop_long_string_length = 8192
  integer,parameter::coop_default_array_size = 4096  ! default array size for interpolation
  integer,parameter::coop_default_chebyshev_fit_order  = 15 !!for chebyshev fit
  COOP_REAL, parameter:: coop_infinity = 1.e30_coop_integer_length

  !!math constants
  COOP_REAL, parameter:: coop_pi = asin(1.d0)*2.d0
  COOP_REAL,parameter:: coop_Riemannzeta3 = 1.2020569031595942853997d0
  COOP_REAL,parameter:: coop_Riemannzeta5  = 1.0369277551433699263313d0
  COOP_REAL,parameter:: coop_Riemannzeta7  = 1.0083492773819228268397d0
  COOP_REAL,parameter:: coop_Riemannzeta9 = 1.00200839282608221441785d0



  COOP_STRING, parameter:: coop_version = "0.0"

  !!minimal a
  COOP_REAL, parameter:: coop_min_scale_factor = 1.d-9


  !!maximum number of species
  COOP_INT, parameter::coop_max_num_species = 16

  !!for dark energy
  COOP_INT, parameter:: coop_max_num_DE_params  = 16
 
  !!for primordial power spectra 
  COOP_INT, parameter:: coop_max_num_PP_params = 32



end module coop_constants
