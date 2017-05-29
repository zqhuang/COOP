module coop_cl_indices_mod
  !!this makes the code faster and more accurate
  use coop_wrapper_utils
#include "constants.h"

  COOP_INT, parameter::coop_init_level_set_species = 0
  COOP_INT, parameter::coop_init_level_set_background = 1
  COOP_INT, parameter::coop_init_level_set_xe = 2
  COOP_INT, parameter::coop_init_level_set_pp = 3
  COOP_INT, parameter::coop_init_level_set_pert = 4
  COOP_INT, parameter::coop_init_level_set_Cls = 5
  COOP_INT, parameter::coop_init_level_set_tens = 6



  logical,parameter :: coop_firstorder_optimize = .true.
  COOP_INT, parameter :: coop_limber_ell = 550

  
  COOP_INT :: coop_Cls_lmax(0:2) = (/ 3100, 2000, 1500 /)

  COOP_REAL, parameter :: coop_power_kmin = 0.02d0 
  COOP_REAL, parameter :: coop_power_kmax = 1.d4
  COOP_REAL, parameter :: coop_visibility_amin = 1.8d-4
  COOP_REAL, parameter :: coop_initial_condition_epsilon = 2.d-7
  COOP_REAL, parameter :: coop_cosmology_firstorder_ode_accuracy = 2.d-7
  COOP_REAL, parameter :: coop_cosmology_firstorder_tc_cutoff = 0.01d0
  COOP_REAL, parameter :: coop_power_lnk_min = log(coop_power_kmin)
  COOP_REAL, parameter :: coop_power_lnk_max = log(coop_power_kmax)


  COOP_REAL, dimension(0:2), parameter :: coop_source_tau_step_factor = (/ 1.d0, 1.d0, 1.d0 /)
  COOP_REAL, dimension(0:2), parameter :: coop_source_k_weight = (/ 0.15d0, 0.15d0, 0.1d0 /)
  COOP_INT, dimension(0:2), parameter :: coop_source_k_n = (/ 420, 160, 160 /)
  COOP_REAL :: coop_source_k_index = 0.85d0


  COOP_INT, parameter :: coop_index_source_T = 1
  COOP_INT, parameter :: coop_index_source_E = 2
  COOP_INT, parameter :: coop_index_source_B = 3
#if DO_ZETA_TRANS
  COOP_INT, parameter :: coop_k_dense_fac = 50
  COOP_INT, parameter :: coop_index_source_zeta = 3
  COOP_INT, parameter :: coop_index_source_Len = 4
  COOP_REAL::coop_zeta_single_slice_chi = -1.d0  !!if set to negative, weight = visibility function
  type(coop_function)::coop_zeta_user_specified_weight
#else
  COOP_INT, parameter :: coop_k_dense_fac = 30
  COOP_INT, parameter :: coop_index_source_Len = 3
#endif

  COOP_INT, parameter::coop_index_ClTT = 1
  COOP_INT, parameter::coop_index_ClEE = 2
  COOP_INT, parameter::coop_index_ClBB = 3
  COOP_INT, parameter::coop_index_ClTE = 4
  COOP_INT, parameter::coop_index_ClLenLen = 5
  COOP_INT, parameter::coop_index_ClTLen = 6


  !!how many source terms you want to extract & save

#if DO_ZETA_TRANS
  COOP_INT, parameter::coop_index_ClTzeta = 7
  COOP_INT, parameter::coop_index_ClEzeta = 8
  COOP_INT, parameter::coop_index_Clzetazeta = 9
  COOP_INT, parameter::coop_num_Cls =  coop_index_Clzetazeta
  COOP_INT, dimension(0:2), parameter::coop_num_sources = (/ 4,  3,  3 /)
#else
  COOP_INT, parameter::coop_num_Cls =  coop_index_ClTLen
  COOP_INT, dimension(0:2), parameter::coop_num_sources = (/ 3,  3,  3 /)
#endif
  COOP_INT, dimension(0:2), parameter::coop_num_saux = (/ 7,  1,  1 /)
  !!odd numbers are used for splining  
  COOP_INT, parameter::coop_aux_index_Weyl = 2 
  COOP_INT, parameter::coop_aux_index_Psi = 4
  COOP_INT, parameter::coop_aux_index_delta_sync = 6  
  !!for scalar the auxiliary variables are
  !!   \Pi, Phi, Phi2, Psi, Psi2, delta_sync/k^2, delta_sync2/k^2
  
end module coop_cl_indices_mod
