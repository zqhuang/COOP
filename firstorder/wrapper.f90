module coop_wrapper_firstorder
  use coop_wrapper_background
  use coop_cl_indices_mod
  use coop_firstorder_mod
  use coop_pertobj_mod
  use coop_lensing_mod
  use coop_cls_postprocess_mod
  implicit none

#include "constants.h"
  
contains

  subroutine coop_wrapper_firstorder_print()
    write(*,*) "This is COOP VERSION "//trim(coop_version)
    write(*,*) "Wrapper for firstorder"
  end subroutine coop_wrapper_firstorder_print

end module coop_wrapper_firstorder
