module coop_firstorder_mod
  use coop_wrapper_background
  use coop_recfast_mod
  implicit none
#include "constants.h"
  
  type, extends(coop_cosmology_background) :: coop_cosmology_firstorder
     type(coop_function)::Ps, Pt, Xe
     COOP_REAL::kpiv_scalar, kpiv_tensor
   contains
     procedure:: setpower => coop_cosmology_firstorder_setpower
     procedure:: setxe => coop_cosmology_firstorder_setxe
  end type coop_cosmology_firstorder

contains

  subroutine coop_cosmology_firstorder_setpower(this, fps, fpt)
    class(coop_cosmology_firstorder)::this
    type(coop_function)::fps, fpt
    this%ps = fps
    this%pt = fpt
  end subroutine coop_cosmology_firstorder_setpower

  subroutine coop_cosmology_firstorder_setxe(this, fxe)
    class(coop_cosmology_firstorder)::this
    type(coop_function)::fxe
    this%xe = fxe

  end subroutine coop_cosmology_firstorder_setxe


end module coop_firstorder_mod
