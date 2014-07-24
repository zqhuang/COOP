#include "default_cosmology_parameters.h"
#define COOP_YES 1
#define COOP_NO 0

#define COOP_INT integer(coop_integer_length)
#define COOP_INT_ARRAY integer(coop_integer_length),dimension(coop_default_array_size)
#define COOP_LONG_INT integer(coop_long_integer_length)
#define COOP_REAL  real(coop_real_length)
#define COOP_REAL_ARRAY  real(coop_real_length),dimension(coop_default_array_size)
#define COOP_COMPLEX complex(coop_complex_length)
#define COOP_COMPLEX_ARRAY  complex(coop_complex_length),dimension(coop_default_array_size)
#define COOP_STRING character(len=coop_string_length)
#define COOP_SHORT_STRING character(len=coop_short_string_length)
#define COOP_LONG_STRING character(len=coop_long_string_length)
#define COOP_UNKNOWN_STRING character(len=*)
#define COOP_CONSTANT_STRING character(len=*),parameter
#define COOP_REAL_OF(x) real(x, coop_real_length)
#define COOP_INT_OF(x) int(x, coop_int_length)

#define COOP_INTERPOLATE_LINEAR 1
#define COOP_INTERPOLATE_QUDRATIC 2
#define COOP_INTERPOLATE_SPLINE 3
#define COOP_INTERPOLATE_CHEBYSHEV 4

#define COOP_ODE_DVERK 1
#define COOP_ODE_RK4 2
#define COOP_ODE_RK6 3
#define COOP_ODE_RK8 4
#define COOP_ODE_GL4 5
#define COOP_ODE_GL6 6
#define COOP_ODE_GL8 7


#define COOP_SPECIES_FLUID 0
#define COOP_SPECIES_CDM 1
#define COOP_SPECIES_MASSIVE_FERMION 2
#define COOP_SPECIES_MASSIVE_BOSON 3
#define COOP_SPECIES_MASSLESS 4
#define COOP_SPECIES_SCALAR_FIELD 5
#define COOP_SPECIES_COSMOLOGICAL_CONSTANT 6
#define COOP_SPECIES_LAMBDA COOP_SPECIES_COSMOLOGICAL_CONSTANT 


#define COOP_COSMO  coop_global_cosmology
#define COOP_COSMO_PARAMS  coop_global_cosmological_parameters

#define COOP_OMEGABH2  COOP_COSMO_PARAMS%r(1)
#define COOP_OMEGACH2  COOP_COSMO_PARAMS%r(2)
#define COOP_100THETA COOP_COSMO_PARAMS%r(3)
#define COOP_TAU COOP_COSMO_PARAMS%r(4)
#define COOP_OMEGAK COOP_COSMO_PARAMS%r(5)
#define COOP_MNU    COOP_COSMO_PARAMS%r(6)

#define COOP_DE_MODEL  COOP_COSMO_PARAMS%i(1)
#define COOP_INDEX_DE      COOP_COSMO_PARAMS%i(2)
#define COOP_NUM_DE     COOP_COSMO_PARAMS%i(3)
#define COOP_PP_MODEL     COOP_COSMO_PARAMS%i(4)
#define COOP_INDEX_PP COOP_COSMO_PARAMS%i(5)
#define COOP_NUM_PP    COOP_COSMO_PARAMS%i(6)
#define COOP_PP_PARAMS(i) COOP_COSMO_PARAMS%r(COOP_INDEX_PP + i - 1)
#define COOP_LN10TO10AS COOP_COSMO_PARAMS%r(COOP_INDEX_PP)
#define COOP_NS COOP_COSMO_PARAMS%r(COOP_INDEX_PP+1)
#define COOP_NRUN COOP_COSMO_PARAMS%r(COOP_INDEX_PP+2)
#define COOP_NRUNRUN COOP_COSMO_PARAMS%r(COOP_INDEX_PP+3)
#define COOP_AMP_RATIO COOP_COSMO_PARAMS%r(COOP_INDEX_PP+4)
#define COOP_NT COOP_COSMO_PARAMS%r(COOP_INDEX_PP+5)
#define COOP_NTRUN COOP_COSMO_PARAMS%r(COOP_INDEX_PP+6)

#define COOP_DE_COSMOLOGICAL_CONSTANT 0
#define COOP_DE_W0 1
#define COOP_DE_W0WA 2
#define COOP_DE_QUINTESSENCE 3
#define COOP_DE_COUPLED_QUINTESSENCE 4

#define COOP_PP_STANDARD 0
#define COOP_PP_SCAN_SPLINE 1
#define COOP_PP_SCAN_LINEAR 2
#define COOP_PP_GENERAL_SINGLE_FIELD 3

#define COOP_RING 1
#define COOP_NESTED 2


#define COOP_MODULAS_SQUARE(z)  (real(z)**2 + aimag(z)**2)
#define COOP_MULT_REAL(z1, z2)  (real(z1)*real(z2) + aimag(z1)*aimag(z2))
#define COOP_POLAR_ANGLE(x, y)  (atan2(y, x+sign(tiny(x), x)))

#define COOP_FORMAT_STRING 0
#define COOP_FORMAT_INTEGER 1
#define COOP_FORMAT_FLOAT 2
#define COOP_FORMAT_LOGICAL -1
