module coop_lattice_mod
  use coop_wrapper_utils
  implicit none
#include "constants.h"

  type coop_lattice
     COOP_INT::n = 0
     COOP_INT::nfld = 1
     logical::metric_backreaction = .false.
     COOP_REAL,dimension(:,:,:,:),allocatable::f
     COOP_REAL,dimension(:,:,:),allocatable::Phi, Psi, hij
  end type coop_lattice
  
contains

  
end module coop_lattice_mod
