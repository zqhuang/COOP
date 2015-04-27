program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  COOP_REAL::Q0 = 0.d0
  COOP_REAL::tracking_n = 0.2d0
  COOP_INT, parameter::nk = 256
  type(coop_cosmology_firstorder)::cosmology
  COOP_REAL k(nk), matterPk(nk), khMpc(nk)
  COOP_INT iQ
  type(coop_asy)::fig
  !!initialize cosmology
  call cosmology%set_standard_cosmology(Omega_b=0.05d0, Omega_c=0.25d0, h = 0.72d0, tau_re = 0.078d0, As = 2.195d-9, ns = 0.9655d0)
  call cosmology%compute_source(0)! , de_Q = Q0, de_tracking_n = tracking_n, de_dlnQdphi = 0.d0, de_dUdphi = 0.d0, de_d2Udphi2 = 0.d0 )
  print*, cosmology%growth_of_z(0.d0), cosmology%growth_of_z(20.d0)
end program test
