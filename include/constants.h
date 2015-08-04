#include "default_cosmology_parameters.h"
#define COOP_YES 1
#define COOP_NO 0

#define DO_ZETA_TRANS COOP_YES
#define DO_EFT_DE COOP_YES

#define COOP_INT integer(coop_integer_length)
#define COOP_INT_ARRAY integer(coop_integer_length),dimension(coop_default_array_size)
#define COOP_LONG_INT integer(coop_long_integer_length)
#define COOP_REAL  real(coop_real_length)
#define COOP_SINGLE  real(coop_single_real_length)
#define COOP_REAL_ARRAY  real(coop_real_length),dimension(coop_default_array_size)
#define COOP_COMPLEX complex(coop_complex_length)
#define COOP_SINGLE_COMPLEX complex(coop_single_complex_length)
#define COOP_COMPLEX_ARRAY  complex(coop_complex_length),dimension(coop_default_array_size)
#define COOP_STRING character(len=coop_string_length)
#define COOP_SHORT_STRING character(len=coop_short_string_length)
#define COOP_LONG_STRING character(len=coop_long_string_length)
#define COOP_UNKNOWN_STRING character(len=*)
#define COOP_CONSTANT_STRING character(len=*),parameter
#define COOP_REAL_OF(x) real(x, coop_real_length)
#define COOP_INT_OF(x) int(x, coop_int_length)
#define COOP_STR_OF(x) trim(coop_num2str(x))
#define COOP_NICESTR_OF(x) trim(coop_num2goodstr(x))
#define COOP_FILESTR_OF(x) trim(coop_num2goodstr(x, "-", "pt"))

#define COOP_INTERPOLATE_LINEAR 1
#define COOP_INTERPOLATE_QUADRATIC 2
#define COOP_INTERPOLATE_SPLINE 3
#define COOP_INTERPOLATE_CHEBYSHEV 4
#define COOP_INTERPOLATE_POLYNOMIAL 5

#define COOP_ODE_DVERK 1
#define COOP_ODE_RK4 2
#define COOP_ODE_RK6 3
#define COOP_ODE_RK8 4
#define COOP_ODE_GL4 5
#define COOP_ODE_GL6 6
#define COOP_ODE_GL8 7

#define COOP_PROPER_SCALE_FACTOR(a) max(min(a,coop_scale_factor_today),coop_min_scale_factor)

#define COOP_SPECIES_FLUID 0
#define COOP_SPECIES_CDM 1
#define COOP_SPECIES_MASSIVE_FERMION 2
#define COOP_SPECIES_MASSIVE_BOSON 3
#define COOP_SPECIES_MASSLESS 4
#define COOP_SPECIES_SCALAR_FIELD 5
#define COOP_SPECIES_COSMOLOGICAL_CONSTANT 6
#define COOP_SPECIES_LAMBDA COOP_SPECIES_COSMOLOGICAL_CONSTANT
#define COOP_SPECIES_COUPLED 7
#define COOP_SPECIES_EFT 8


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

#define COOP_INFLATION_CONSISTENCY  COOP_COSMO_PARAMS%l(1)

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
#define COOP_PP_BUMP  4  

#define COOP_UNKNOWN_ORDERING 0
#define COOP_RING 1
#define COOP_NESTED 2


#define COOP_MODULAS_SQUARE(z)  (real(z)**2 + aimag(z)**2)
#define COOP_MULT_REAL(z1, z2)  (real(z1)*real(z2) + aimag(z1)*aimag(z2))
#define COOP_POLAR_ANGLE(x, y)  (atan2(y, x+sign(tiny(x), x)))
#define COOP_LN(x)  (log(max(x, tiny(x))))

#define COOP_FORMAT_STRING 0
#define COOP_FORMAT_INTEGER 1
#define COOP_FORMAT_FLOAT 2
#define COOP_FORMAT_LOGICAL -1


#define O0_BARYON(x) x%species(x%index_baryon)
#define O0_CDM(x) x%species(x%index_cdm)
#define O0_RADIATION(x) x%species(x%index_radiation)
#define O0_NU(x)  x%species(x%index_nu)
#define O0_MASSIVENU(x) x%species(x%index_massivenu)
#define O0_DE(x)  x%species(x%index_de)

#define O1_PSI y(pert%metric%i(0))
#define O1_PSIPR y(pert%metric%i(1))
#define O1_VEC_V y(pert%metric%i(1))
#define O1_TEN_H y(pert%metric%i(2))
#define O1_TEN_HPR y(pert%metric%i(3))
#define O1_DELTA_B y(pert%baryon%i(0))
#define O1_V_B y(pert%baryon%i(1)) 
#define O1_DELTA_C y(pert%cdm%i(0))
#define O1_V_C y(pert%cdm%i(1))
#define O1_T(l)  y(pert%T%i(l))
#define O1_E(l)  y(pert%E%i(l))
#define O1_B(l)  y(pert%B%i(l))
#define O1_NU(l) y(pert%nu%i(l))
#define O1_MASSIVENU(l, iq) y(pert%massivenu(iq)%i(l))
#define O1_DE(l) y(pert%de%i(l))
#define O1_DELTA_DE O1_DE(0)
#define O1_V_DE  O1_DE(1)
#define O1_DELTA_PHI O1_DE(0)
#define O1_DELTA_PHIPR O1_DE(1)

#define O1_DE_HPI O1_DE(0)    
#define O1_DE_HPIPR O1_DE(1)

#define O1_PHI  pert%O1_Phi
#define O1_PHI_PRIME  pert%O1_phipr
#define O1_PSI_PRIME yp(pert%metric%i(0))
#define O1_PSIPR_PRIME yp(pert%metric%i(1))
#define O1_VEC_V_PRIME yp(pert%metric%i(1))
#define O1_TEN_H_PRIME yp(pert%metric%i(2))
#define O1_TEN_HPR_PRIME yp(pert%metric%i(3))
#define O1_DELTA_B_PRIME yp(pert%baryon%i(0))
#define O1_V_B_PRIME yp(pert%baryon%i(1)) 
#define O1_DELTA_C_PRIME yp(pert%cdm%i(0))
#define O1_V_C_PRIME yp(pert%cdm%i(1))
#define O1_T_PRIME(l)  yp(pert%T%i(l))
#define O1_E_PRIME(l)  yp(pert%E%i(l))
#define O1_B_PRIME(l)  yp(pert%B%i(l))
#define O1_NU_PRIME(l) yp(pert%nu%i(l))
#define O1_MASSIVENU_PRIME(l, iq) yp(pert%massivenu(iq)%i(l))
#define O1_DE_PRIME(l) yp(pert%de%i(l))
#define O1_DELTA_DE_PRIME O1_DE_PRIME(0)
#define O1_V_DE_PRIME  O1_DE_PRIME(1)
#define O1_DELTA_PHI_PRIME O1_DE_PRIME(0)
#define O1_DELTA_PHIPR_PRIME O1_DE_PRIME(1)
  
#define O1_DE_HPI_PRIME  O1_DE_PRIME(0)
#define O1_DE_HPIPR_PRIME O1_DE_PRIME(1)

#define COOP_PERT_NONE 0
#define COOP_PERT_METRIC 1
#define COOP_PERT_PERFECT_FLUID 2
#define COOP_PERT_HIERARCHY 3
#define COOP_PERT_SCALAR_FIELD 4
#define COOP_PERT_EFT 5

#define COOP_INTERP_SOURCE(source, ind, idense, ik, itau) (source%s(ind, ik, itau)*source%a_dense(idense) + source%s(ind, ik-1, itau)*source%b_dense(idense) + source%s2(ind, ik, itau)*source%a2_dense(idense) + source%s2(ind, ik-1, itau)*source%b2_dense(idense))
