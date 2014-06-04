program test
  use coop_wrapper
  implicit none
#include "constants.h"
#define PHI y(1)
#define PHIDOT y(2)
#define RHOM  y(3)
#define DPHI_DLNA yp(1)    
#define DPHIDOT_DLNA yp(2)    
#define DRHOM_DLNA yp(3)
#define Q_COUPLING args%r(1)
#define C_RUN args%r(2)
#define C_FLAT args%r(3)
#define N_POWER args%i(1)

  COOP_INT, parameter::nvars = 3  !!number of variables to evolve
  COOP_INT, parameter::nsteps = 1000  !!number of steps for output
  COOP_REAL, parameter:: MPl = 1.d0 !!define mass unit
  COOP_REAL, parameter:: H0 = 1.d0  !!define time unit

  COOP_INT, parameter::npower = 1
  COOP_REAL, parameter:: Q = 1./coop_sqrt6
  COOP_REAL,parameter:: rhom_ini = 1.e4 * (Mpl**2*H0**2)
  COOP_REAL,parameter:: phieq_ini = 0.01 * Mpl  
  COOP_REAL, parameter::c_flat = rhom_ini* 1.e-4
  COOP_REAL, parameter::c_run = Q * rhom_ini * phieq_ini ** (npower + 1)/npower
  COOP_REAL, parameter::lna_start = 0.
  COOP_REAL, parameter::lna_end = log(rhom_ini/c_flat)/3.d0 + 3.d0  !!add  3 more efolds to make sure rho_m < rho_phi happens

  type(coop_ode)::co
  type(coop_arguments)::args
  COOP_INT i
  COOP_REAL::phieq_dot_ini, hubble_ini
  COOP_REAL:: phi(nsteps), lna(nsteps), dotphi(nsteps), w(nsteps), Ev, Ek
  type(coop_asy)::fp

  args = coop_arguments( i = (/ npower /),  r = (/ Q, c_run, c_flat /) )
  hubble_ini = sqrt((potential(phieq_ini, args) + rhom_ini)/(3.d0*Mpl**2))
  phieq_dot_ini = 3.d0*hubble_ini/(npower+1.d0)*phieq_ini
  call co%init(n = nvars, method = COOP_ODE_DVERK)  !!initialize the ode solver
  call co%set_arguments(args = args)
  call co%set_initial_conditions( xini = lna_start, yini = (/ phieq_ini, phieq_dot_ini, rhom_ini /) )
  call coop_set_uniform(nsteps, lna, lna_start, lna_end)
  do  i = 1, nsteps
     if(i.gt.1)call co%evolve(get_yprime, lna(i))
     phi(i) = co%PHI / Mpl
     dotphi(i) = co%PHIDOT
     Ek = dotphi(i)**2/2.d0
     Ev = potential(phi(i), args)
     w(i) = (Ek -Ev)/(Ek+Ev)
  enddo
  call fp%open("phi.txt")
  call fp%init(xlabel = "$\ln a$", ylabel = "$\phi / M_p$")
  call coop_asy_curve(fp, x=lna, y=phi, color = "black", linewidth = 1.5, linetype="solid")
  call fp%close()

  call fp%open("w.txt")
  call fp%init(xlabel = "$\ln a$", ylabel = "$w_\phi$")
  call coop_asy_curve(fp, x=lna, y=w, color = "black", linewidth = 1.5, linetype="solid")
  call fp%close()

contains


  function potential(phi, args) result(V)
    COOP_REAL phi, V
    type(coop_arguments) args
    V = C_RUN / phi**N_POWER + C_FLAT
  end function potential

  function dVdphi(phi, args) 
    COOP_REAL phi, dVdphi
    type(coop_arguments) args
    dVdphi = - C_RUN * N_POWER/phi**(N_POWER +1)
  end function dVdphi

  subroutine get_yprime(n, x, y, yp, args)
    COOP_INT n
    COOP_REAL x, y(n), yp(n), hubble, Ek, Ep, Vp
    type(coop_arguments) args
    Ek = PHIDOT ** 2 /2.d0
    Ep = potential(PHI, args)
    Vp = dVdphi(PHI, args)
    hubble = sqrt((Ek + Ep + RHOM)/(3.d0*MPl**2))
    DPHI_DLNA = PHIDOT / hubble
    DPHIDOT_DLNA = (-3.d0*hubble*PHIDOT - Vp - Q_COUPLING/Mpl * RHOM)/hubble
    DRHOM_DLNA = (-3.d0 +  Q_COUPLING/Mpl * PHIDOT / hubble)*RHOM
  end subroutine get_yprime

end program test
