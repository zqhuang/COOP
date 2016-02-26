#define USER_PARAM(i)  args%r(i)
!! you are allowed to use up to 10 parameters for primordial power 

!!Input:
!!kbykpiv =  k / k_{pivot}
!!cosmology
!!args (contains the user defined parameters)
!!Output:
!!ps = dimensionless scalar  primordial  power spectrum
!!pt = dimensionless tensor primordial power spectrum

!!100 is a sample code (not used anywhere); do not change it; the user-defined subroutines are from 101 to 110.
subroutine coop_user_defined_primordial_power_100(kbykpiv, ps, pt, cosmology, args) 
  COOP_REAL::kbykpiv, ps, pt
  type(coop_cosmology_firstorder)::cosmology
  type(coop_arguments)::args
  COOP_REAL::nt
#define SCALAR_AS USER_PARAM(1)
#define SCALAR_NS USER_PARAM(2)
#define TENSOR_R  USER_PARAM(3)
  if(cosmology%inflation_consistency)then
     nt = - TENSOR_R/8.d0
  else
     nt = 0.d0
  endif
  ps =  SCALAR_AS * kbykpiv ** (SCALAR_NS - 1.d0)
  pt =  (SCALAR_AS * TENSOR_R) * kbykpiv ** nt
#undef SCALAR_AS
#undef SCALAR_NS
#undef TENSOR_R
end subroutine coop_user_defined_primordial_power_100


subroutine coop_user_defined_primordial_power_101(kbykpiv, ps, pt, cosmology, args) 
  COOP_REAL::kbykpiv, ps, pt
  type(coop_cosmology_firstorder)::cosmology
  type(coop_arguments)::args
!!write your code here
!!ps = ...
!!pt = ...
!!to use this model; put "pp_genre = 101" in your ini file
  stop "primordial power 101 has not been defined."
end subroutine coop_user_defined_primordial_power_101


subroutine coop_user_defined_primordial_power_102(kbykpiv, ps, pt, cosmology, args) 
  COOP_REAL::kbykpiv, ps, pt
  type(coop_cosmology_firstorder)::cosmology
  type(coop_arguments)::args
!!write your code here
!!ps = ...
!!pt = ...
!!to use this model; put "pp_genre = 102" in your ini file
  stop "primordial power 102 has not been defined."
end subroutine coop_user_defined_primordial_power_102


subroutine coop_user_defined_primordial_power_103(kbykpiv, ps, pt, cosmology, args) 
  COOP_REAL::kbykpiv, ps, pt
  type(coop_cosmology_firstorder)::cosmology
  type(coop_arguments)::args
!!write your code here
!!ps = ...
!!pt = ...
!!to use this model; put "pp_genre = 103" in your ini file
  stop "primordial power 103 has not been defined."
end subroutine coop_user_defined_primordial_power_103


subroutine coop_user_defined_primordial_power_104(kbykpiv, ps, pt, cosmology, args) 
  COOP_REAL::kbykpiv, ps, pt
  type(coop_cosmology_firstorder)::cosmology
  type(coop_arguments)::args
!!write your code here
!!ps = ...
!!pt = ...
!!to use this model; put "pp_genre = 104" in your ini file
  stop "primordial power 104 has not been defined."
end subroutine coop_user_defined_primordial_power_104


subroutine coop_user_defined_primordial_power_105(kbykpiv, ps, pt, cosmology, args) 
  COOP_REAL::kbykpiv, ps, pt
  type(coop_cosmology_firstorder)::cosmology
  type(coop_arguments)::args
!!write your code here
!!ps = ...
!!pt = ...
!!to use this model; put "pp_genre = 105" in your ini file
  stop "primordial power 105 has not been defined."
end subroutine coop_user_defined_primordial_power_105

subroutine coop_user_defined_primordial_power_106(kbykpiv, ps, pt, cosmology, args) 
  COOP_REAL::kbykpiv, ps, pt
  type(coop_cosmology_firstorder)::cosmology
  type(coop_arguments)::args
!!write your code here
!!ps = ...
!!pt = ...
!!to use this model; put "pp_genre = 106" in your ini file
  stop "primordial power 106 has not been defined."
end subroutine coop_user_defined_primordial_power_106

subroutine coop_user_defined_primordial_power_107(kbykpiv, ps, pt, cosmology, args) 
  COOP_REAL::kbykpiv, ps, pt
  type(coop_cosmology_firstorder)::cosmology
  type(coop_arguments)::args
!!write your code here
!!ps = ...
!!pt = ...
!!to use this model; put "pp_genre = 107" in your ini file
  stop "primordial power 107 has not been defined."
end subroutine coop_user_defined_primordial_power_107


subroutine coop_user_defined_primordial_power_108(kbykpiv, ps, pt, cosmology, args) 
  COOP_REAL::kbykpiv, ps, pt
  type(coop_cosmology_firstorder)::cosmology
  type(coop_arguments)::args
!!write your code here
!!ps = ...
!!pt = ...
!!to use this model; put "pp_genre = 108" in your ini file
  stop "primordial power 108 has not been defined."
end subroutine coop_user_defined_primordial_power_108


subroutine coop_user_defined_primordial_power_109(kbykpiv, ps, pt, cosmology, args) 
  COOP_REAL::kbykpiv, ps, pt
  type(coop_cosmology_firstorder)::cosmology
  type(coop_arguments)::args
!!write your code here
!!ps = ...
!!pt = ...
!!to use this model; put "pp_genre = 109" in your ini file
  stop "primordial power 109 has not been defined."
end subroutine coop_user_defined_primordial_power_109

subroutine coop_user_defined_primordial_power_110(kbykpiv, ps, pt, cosmology, args) 
  COOP_REAL::kbykpiv, ps, pt
  type(coop_cosmology_firstorder)::cosmology
  type(coop_arguments)::args
!!write your code here
!!ps = ...
!!pt = ...
!!to use this model; put "pp_genre = 110" in your ini file
  stop "primordial power 110 has not been defined."
end subroutine coop_user_defined_primordial_power_110





#undef USER_PARAM
