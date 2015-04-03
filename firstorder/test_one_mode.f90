program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  COOP_REAL::Q0 = 0.d0
  COOP_REAL::tracking_n = 0.01d0
  COOP_INT, parameter::ntau = 1000
  type(coop_cosmology_firstorder)::fod
  COOP_REAL k
  COOP_REAL tau(ntau)
  type(coop_file)::fp
  type(coop_asy)::fig
  !!initialize cosmology
  call fod%set_standard_cosmology(Omega_b=0.0485374d0, Omega_c=0.2585497252d0, h = 0.67766d0, tau_re = 0.08193d0, As = 2.2098d-9, ns = 0.968d0, nrun = 0.d0, r = 0.d0, nt = -0.01d0, YHe = 0.248d0, Nnu = 3.d0, de_Q = Q0, de_tracking_n = tracking_n, de_dlnQdphi = 0.d0, de_dUdphi = 0.d0, de_d2Udphi2 = 0.d0 )
  !!***************************************************
  !! V = V0 / phi^n exp( C1  * phi + 1/2 C2 phi**2)
  !! n = de_tracking_n;  C1 = de_dUdphi; C2 = de_d2Udphi2
  !! V0 is determined by the condition Omega_k =0
  !!***************************************************
  !! Q = Q0 exp( A * phi)
  !! Q0 = de_Q,  A = de_dlnQdphi
  !!***************************************************


  !!set k/H0  
  k = 10.d0

  !!set tau
  call coop_set_uniform(ntau, tau, 5.d-3, fod%TAU0)
  
  
!!test energy conservation
  call fod%init_source(0, k = (/ k /), tau = tau)
  print*, fod%source(0)%ntau
  call fod%compute_source_k(fod%source(0), 1, do_test_energy_conservation = .true.)
  print*, fod%source(0)%saux(3, 1,  fod%source(0)%ntau),  fod%source(0)%saux(2, 1,  fod%source(0)%ntau)- fod%source(0)%saux(3, 1,  fod%source(0)%ntau)
  
  
  


end program test
