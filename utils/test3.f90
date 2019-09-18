program Daubechies
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT,parameter::n=50
  COOP_REAL,dimension(n,n,n)::f
  COOP_REAL,parameter::mean = 10.d0
  COOP_REAL,parameter::norm = 0.1d0
  COOP_INT::i
  COOP_REAL::x
  print*, coop_SI_planckmass*coop_SI_c**2/coop_SI_GeV,coop_SI_reduced_planckmass*coop_SI_c**2/coop_SI_GeV
end program Daubechies
