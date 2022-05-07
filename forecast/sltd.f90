!!strong lensing time delay likelihood module
module coop_sltd_mod
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"

  COOP_INT,parameter::coop_sltd_npt = 6
  
  type coop_sltd_object
     COOP_REAL,dimension(coop_sltd_npt) ::zd_list= (/ &
          0.6304, 0.295, 0.4546, 0.745, 0.6575, 0.311 /)
     COOP_REAL,dimension(coop_sltd_npt)::zs_list = (/ &
          1.394, 0.654, 1.693, 1.789, 1.662, 1.722 /)
     COOP_REAL,dimension(coop_sltd_npt)::ln_Ddt_list  = (/ &
          5156., 2096., 2707., 5769., 4784., 1470.  /)
     COOP_REAL,dimension(coop_sltd_npt)::ln_Dd_list = (/ &
          1228., 804., 1000., 1805., 1200., 697. /)
     COOP_REAL,dimension(coop_sltd_npt)::delta_ln_Ddt  = (/ &
          266., 90.5, 175.5, 530., 323.5, 132. /)
     COOP_REAL,dimension(coop_sltd_npt)::delta_ln_Dd  = (/ &
          266., 90.5, 175.5, 530., 323.5, 132. /)
     COOP_REAL,dimension(coop_sltd_npt):: corr  = (/ &
          0., 0., 0., 0., 0., 0. /)
   contains
     procedure::LogLike => coop_sltd_object_LogLike
  end type coop_sltd_object

contains

  function coop_sltd_object_LogLike(this, cosmology) result(logLike)
    class(coop_sltd_object)::this
    type(coop_cosmology_firstorder)::cosmology
    COOP_REAL:: loglike
    COOP_INT::i
    loglike = 0.d0
    do i = 1, coop_sltd_npt

    enddo
    loglike = loglike/2.d0
  end function coop_sltd_object_LogLike
  
  
end module coop_sltd_mod
