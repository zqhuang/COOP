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
  call fig%open("matterpower.txt")
  call fig%init(xlabel = "$k\, [h/{\rm Mpc}]$", ylabel = "$P(k)\, [(h^{-1}{\rm Mpc})^3]$", xlog = .true., ylog = .true., ymax = 3.2e4, ymin = 1.e2)
  call coop_set_uniform(nk, k, 0.4d0, 2.d3, logscale = .true.)
  do iQ = 0, 2
     Q0 = iQ * 0.1d0
     !!initialize cosmology; Planck 2015 bestfit +  DE models
     call cosmology%set_standard_cosmology(Omega_b=0.04904d0, Omega_c=0.2642d0, h = 0.6731d0, tau_re = 0.078d0, As = 2.195d-9, ns = 0.9655d0, de_Q = Q0, de_epsilon_s = 0.25d0, de_epsilon_inf = 0.2d0 )!, de_tracking_n = tracking_n, de_dlnQdphi = 0.d0, de_dUdphi = 0.d0, de_d2Udphi2 = 0.d0 )
     !!***************************************************
     !! V = V0 / phi^n exp( C1  * phi + 1/2 C2 phi**2)
     !! n = de_tracking_n;  C1 = de_dUdphi; C2 = de_d2Udphi2
     !! V0 is determined by the condition Omega_k =0
     !!***************************************************
     !! Q = Q0 exp( A * phi)
     !! Q0 = de_Q,  A = de_dlnQdphi
     !!***************************************************

     khMpc = k * cosmology%H0Mpc()/cosmology%h()  !!k/H0 * (H0 * Mpc) / h = k in unit of h Mpc^{-1}
     call cosmology%compute_source(0)
     call cosmology%get_MPhi_power(z=0.d0, nk = nk, k = k, Pk = matterPk)  !!this returns k^3 |\delta_k|^2 /(2pi^2)
     matterPk = matterPk * (2.d0*coop_pi**2)/khMpc**3  !!obtain |\delta_k|^2 in unit of (Mpc/h)^3
     call fig%curve(khMpc, matterPk, color=fig%color(iQ+1), linetype=fig%linetype(iQ+1), legend="$Q="//COOP_STR_OF(Q0)//"$")
  enddo
  call fig%legend(0.3, 0.5)
  call fig%close()
  


end program test
