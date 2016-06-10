  COOP_INT, parameter::coop_default_array_size = 8192  ! default array size for interpolation
  COOP_INT, parameter::coop_default_chebyshev_fit_order  = 30 !!for chebyshev fit
  COOP_REAL, parameter:: coop_infinity = 1.e30
  COOP_REAL, parameter:: coop_tiny = 1.d0/coop_infinity
  COOP_REAL, parameter:: coop_logInfinity = log(coop_Infinity)
  COOP_REAL, parameter:: coop_logTiny = - coop_logInfinity

  COOP_INT, parameter::coop_tmp_file_unit = 9


  COOP_REAL, parameter:: coop_primordial_zeta_norm = 1.d0

  !!minimal a
  COOP_REAL, parameter:: coop_min_scale_factor = 1.d-12

  !!this should not be changed
  COOP_REAL, parameter:: coop_scale_factor_today = 1.d0

  !!maximum number of species
  COOP_INT, parameter::coop_max_num_species = 16

  !!for dark energy
  COOP_INT, parameter:: coop_max_num_DE_params  = 16
 
  !!for primordial power spectra 
  COOP_INT, parameter:: coop_max_num_PP_params = 32


  
