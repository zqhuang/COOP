module coop_HSTlike_mod
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"

  !!from Riess et al 2019
  type coop_HST_object
     COOP_REAL:: H0 = 74.03
     COOP_REAL:: H0_err = 1.42
   contains
     procedure::LogLike => coop_HST_object_LogLike
  end type coop_HST_object

contains

  function coop_HST_object_LogLike(this, cosmology) result(logLike)
    class(coop_HST_object)::this
    type(coop_cosmology_firstorder)::cosmology
    COOP_REAL Heff, loglike, dlzeff
    Loglike = ((cosmology%h_value*100.d0 - this%H0)/this%H0_err)**2/2.d0
  end function coop_HST_object_LogLike
  
  
end module coop_HSTlike_mod
