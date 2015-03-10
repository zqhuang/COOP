program bgtest
  use coop_wrapper_background
  implicit none
#include "constants.h"  
  type(coop_cosmology_background)::bg
  COOP_INT, parameter::nvars = 3
  COOP_INT, parameter::n = 512
  COOP_REAL:: lna(n),a(n), Qcpl, tracking_n, phi, phidot, V, Vp, Vpp, weff, Omegac, OmegaL, omegab, dlnQdphi, dUdphi, d2Udphi2, hub, rhocdm, rhode, epsilon_inf, epsilon_s, dlnVdphi_eq, lnVpp_eq, aeq, zeta_s, wp1, at, qwp1, qnow, apiv
  COOP_INT::i, index_CDM , index_DE
  type(coop_ode)::ode
  type(coop_species)::quint_de
  type(coop_file)::fp, fpw
  !!DE  parameters
  Qcpl = 0.5d0 !!coupling between DE and CDM
  tracking_n = 0.d0 !!  V \propto 1 / phi^n e^{U(phi)}   (n > 0)
  dlnQdphi = 0.d0
  dUdphi = 0.d0
  d2Udphi2 = 0.d0
  apiv = 0.1d0
  !!=============  set up background ===============
  omegab = 0.046d0
  omegac = 0.25d0

 
  call bg%init(h=0.7d0)
  call bg%add_species(coop_baryon(omegab))
  call bg%add_species(coop_radiation(bg%Omega_radiation()))
  call bg%add_species(coop_neutrinos_massless(bg%Omega_massless_neutrinos_per_species()*(bg%Nnu())))
  call coop_background_add_coupled_DE(bg, Omega_c = omegac, Q = Qcpl, tracking_n = tracking_n, dlnQdphi = dlnQdphi, dUdphi = dUdphi, d2Udphi2 = d2Udphi2)
  
  index_CDM = bg%index_of("CDM")
  index_DE = bg%index_of("Dark Energy")


  !!
  call bg%setup_background()

  aeq = 1.d0
  do while(bg%species(index_DE)%density(aeq) + bg%species(index_cdm)%density(aeq) .gt.  3.d0*(2.d0*omegac+omegab)/aeq**3)
     aeq = aeq * 0.99
  enddo
  at = aeq/4.d0
  phi=bg%species(index_DE)%DE_phi(aeq)
  dlnVdphi_eq = bg%species(index_DE)%DE_dlnVdphi( phi )
  lnVpp_eq = bg%species(index_DE)%DE_d2lnVdphi2( phi )
  epsilon_s = (dlnVdphi_eq + bg%species(index_DE)%DE_Q(phi) )**2/2.d0
  epsilon_inf = bg%species(index_DE)%wp1ofa(apiv)*1.5d0

  zeta_s = 0.d0
  print*, "a_eq = ", aeq
  print*, "epsilon_s = ", epsilon_s
  print*, "epsilon_infty = ", epsilon_inf
  print*, "zeta_s = ", zeta_s
  quint_de = coop_de_quintessence(1.d0-omegab-omegac, epsilon_s, epsilon_inf, zeta_s)
  print*, quint_de%wp1ofa(apiv), bg%species(index_DE)%wp1ofa(apiv)
  !!=============== test energy conservation ===========
#define LNHUBBLE y(1)
#define LNRHODE y(2)
#define LNRHOCDM y(3)    
#define DLNHUBBLE_DLNA yp(1)
#define DLNRHODE_DLNA yp(2)
#define DLNRHOCDM_DLNA yp(3)  
  
  call coop_set_uniform(n, lna, -15.d0, 0.d0)
  a = exp(lna)
  call ode%init(nvars)
  call ode%set_initial_conditions(lna(1), (/ log(bg%Hratio(a(1))),  log(bg%species(index_DE)%density(a(1))), log(bg%species(index_cdm)%density(a(1))) /) )
  call fpw%open("w_output.dat","w")
  call fp%open("conservation.dat", "w")
  write(fpw%unit, "(20A19)") "# a", " 1 + w_phi ", " 1 + w_phi_approx ", " 1 + w_DE_eff ", "1 + w_DE_eff_approx"
  write(fp%unit, "(20A19)") "#ln a", "ln H", "ln sqrt(rho/3)",  "ln rho_DE_int", "ln rho_DE_theory", "ln rho_c_int", "ln rho_c_theory"
  do i=1, n
     if(i.gt.1)call ode%evolve(fcn, lna(i))
     phi = bg%species(index_DE)%DE_phi(a(i))
     phidot = bg%species(index_DE)%DE_phidot(a(i))
     V =  bg%species(index_DE)%DE_V(phi)
     hub = exp(ode%LNHUBBLE)
     rhocdm = bg%species(index_CDM)%density(a(i))
     rhode =bg%species(index_DE)%density(a(i))
     wp1 = bg%species(index_DE)%wp1ofa(a(i))
     qwp1 = quint_de%wp1ofa(a(i))
     qnow = bg%species(index_DE)%DE_Q(phi)
     if(a(i).gt. 0.2)then
        weff = - ( (bg%species(index_CDM)%dlnrhodlna(a(i))+3.d0)*rhocdm + (bg%species(index_DE)%dlnrhodlna(a(i))+3.d0)*rhode)/( Rhocdm + rhode - 3.d0*Omegac/a(i)**3 )/3.d0

        write(fpw%unit, "(20E19.9)") a(i), wp1, qwp1, weff+1.d0, qwp1 +qnow*sqrt(qwp1*omegac/(1.d0-omegab-omegac))/coop_sqrt3*sqrt(omegac/((omegac+omegab)*(a(i)**3+at**3) + (1.d0-omegac-omegab)*a(i)**6))
     endif
     write(fp%unit, "(20E19.9)")  lna(i), ode%LNHUBBLE, log(bg%Hratio(a(i))), ode%LNRHODE, log(rhode), ode%LNRHOCDM, log(rhocdm)  !!check energy conservation etc.
  enddo
  call fp%close()
  call fpw%close()


  !!===============  use gnuplot  to plot the results ============
  !!test energy conservation:  
  !!    plot "background_output.dat" u 1:($2-$3) w l


contains

  
  subroutine fcn(n, lna, y, yp)
    COOP_INT n
    COOP_REAL lna, y(n), yp(n), a
    a = exp(lna)
    DLNHUBBLE_DLNA = bg%HdotbyHsq(a)
    DLNRHODE_DLNA = -3.d0*bg%species(index_DE)%wp1effofa(a)
    DLNRHOCDM_DLNA = -3.d0*bg%species(index_cdm)%wp1effofa(a)    
  end subroutine fcn
  
end program bgtest
