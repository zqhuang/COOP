#include "default_cosmology_parameters.h"

#define COOP_INT integer(coop_integer_length)
#define COOP_INT_ARRAY integer(coop_integer_length),dimension(coop_default_array_size)
#define COOP_REAL  real(coop_real_length)
#define COOP_REAL_ARRAY  real(coop_real_length),dimension(coop_default_array_size)
#define COOP_STRING character(len=coop_string_length)
#define COOP_SHORT_STRING character(len=coop_short_string_length)
#define COOP_LONG_STRING character(len=coop_long_string_length)
#define COOP_UNKNOWN_STRING character(len=*)
#define COOP_REAL_OF(x) real(x, coop_real_length)
#define COOP_INT_OF(x) int(x, coop_int_length)

#define COOP_INTERPOLATE_LINEAR 1
#define COOP_INTERPOLATE_QUDRATIC 2
#define COOP_INTERPOLATE_SPLINE 3
#define COOP_INTERPOLATE_CHEBYSHEV 4


#define COOP_SPECIES_FLUID 0
#define COOP_SPECIES_CDM 1
#define COOP_SPECIES_MASSIVE_FERMION 2
#define COOP_SPECIES_MASSIVE_BOSON 3
#define COOP_SPECIES_MASSLESS 4
#define COOP_SPECIES_SCALAR_FIELD 5
#define COOP_SPECIES_COSMOLOGICAL_CONSTANT 6
#define COOP_SPECIES_LAMBDA COOP_SPECIES_COSMOLOGICAL_CONSTANT 

