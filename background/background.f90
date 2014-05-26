module coop_background
  use coop_constants
  use coop_types
  implicit none
#include "background.h"

private

public:: coop_de_lambda, coop_de_quintessence, coop_de_scalar_tensor

  type, extends(coop_species_dynamic)::coop_de_lambda
     contains 
       procedure:: wofa => coop_de_lambda_wofa
       procedure:: cs2ofa => coop_de_lambda_cs2ofa
  end type coop_de_lambda

contains

  function coop_de_lambda_wofa(this, a) result(w)
    class(coop_de_lambda)::this
    COOP_REAL a, w
    w = -1.
  end function coop_de_lambda_wofa
  
  function coop_de_lambda_cs2ofa(this, a) result(cs2)
    class(coop_de_lambda)::this
    COOP_REAL a, cs2
    cs2 = 1.
  end function coop_de_lambda_cs2ofa


end module coop_background
