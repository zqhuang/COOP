module coop_Chronometer_mod
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"

  COOP_INT,parameter::coop_cc_npt = 31
  
  type coop_cc_object
     COOP_REAL,dimension(coop_cc_npt) ::z_list= (/ 0.07, 0.09, 0.12, 0.17, 0.179,  &
          0.199, 0.2, 0.27, 0.28, 0.352, &
          0.3802, 0.4, 0.4004, 0.4247, 0.4497, &
          0.4783, 0.48, 0.593, 0.68, 0.781, &
          0.875, 0.88, 0.9, 1.037, 1.3, &
          1.363, 1.430, 1.53, 1.75, 1.965, 0.47 /)
     COOP_REAL,dimension(coop_cc_npt)::H_list = (/ 69., 69., 68.6, 83., 75., &
          75., 72.9, 77., 88.8, 83., &
          83., 95., 77., 87.1, 92.8, &
          80.9, 97., 104., 92., 105., &
          125., 90., 117., 154., 168., &
          160., 177., 140., 202., 186.5, 89. /)
     COOP_REAL,dimension(coop_cc_npt)::dH_list = (/ 19.6, 12., 26.2, 8., 4., &
          5., 29.6, 14., 36.6,  14., &
          13.5, 17., 10.2, 11.2, 12.9, &
          9., 62., 13., 8., 12., &
          17., 40., 23., 20., 17., &
          33.6, 18., 14., 40., 50.4, 23. /)
   contains
     procedure::LogLike => coop_cc_object_LogLike
  end type coop_cc_object

contains

  function coop_cc_object_LogLike(this, cosmology) result(logLike)
    class(coop_cc_object)::this
    type(coop_cosmology_firstorder)::cosmology
    COOP_REAL:: loglike
    COOP_INT::i
    loglike = 0.d0
    do i = 1, coop_cc_npt
       loglike = loglike + ((cosmology%H_of_z(this%z_list(i))*100.*cosmology%h_value - this%H_list(i))/this%dH_list(i))**2
    enddo
    loglike = loglike/2.d0
  end function coop_cc_object_LogLike
  
  
end module coop_Chronometer_mod
