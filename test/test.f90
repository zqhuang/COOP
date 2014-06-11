program combined2
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
#define N_POWER args%r(4)

  !!example for coupled quintessence model 
  !!V(phi) = c0 / phi^n + c1

  COOP_INT, parameter::nvars = 3  !!number of variables to evolve
  COOP_INT, parameter::nsteps = 10000  !!number of steps for output
  COOP_INT, parameter::nminimize = 10000 !!number of elements in input const_integration arrays for minimization procedure
  COOP_REAL, parameter:: MPl = 1.d0 !!define mass unit
  COOP_REAL, parameter:: H0 = 1.d0  !!define time unit
  COOP_REAL, parameter:: rho_total_today = 3.d0*Mpl**2*H0**2
  COOP_REAL, parameter:: rhom_today = 0.3*rho_total_today

  COOP_REAL, parameter::npower = 0.8d0
  COOP_REAL:: Q , c_run
  COOP_REAL,parameter:: rhom_ini = 1.e8 * rhom_today
  COOP_REAL,parameter:: phieq_ini = 1.d-4 * Mpl  
  COOP_REAL, parameter::c_flat = 0.7 * rho_total_today

  COOP_REAL, parameter::lna_start = 0.
  COOP_REAL, parameter::lna_end = log(rhom_ini/c_flat)/3.d0 + 3.d0  !!add  3 more efolds to make sure rho_m < rho_phi happens
  COOP_INT iQ, itoday, ieq, iplot, imiddle 
  COOP_REAL::plot_a_power = 1.
  COOP_REAL eps_inf

  integer,parameter::iQ_steps = 3
  type(coop_ode)::co
  type(coop_arguments)::args
  COOP_INT i, j, it
  COOP_REAL::phieq_dot_ini, phieq_dot_dot_ini, hubble_ini, hubble_dot_ini, discriminant, deltaphi
  COOP_REAL:: phi(nsteps), phieq(nsteps), lna(nsteps), dotphi(nsteps), w(nsteps), e_s(nsteps), dotphieq(nsteps), dotdotphi(nsteps), dotdotphieq(nsteps), Ev, Ek, Vpp, F(nsteps)
  COOP_REAL:: e_s_eq, beta
  COOP_REAL  w_slow(nsteps), delta_w(nsteps), delta_w_ana(nsteps), w_ana(nsteps), a(nsteps)

  type(coop_asy):: phiplot, dot, dotdot, es, wplot, dwplot

  call wplot%open("w.txt")
  call dwplot%open("dw.txt")
  call wplot%init(xlabel = "$a$", ylabel = "$w_\phi$", caption = "$V=C_0+C_1\phi^{-"//trim(coop_num2str(npower))//"}$",width=8., height=6., ymin  = -1., ymax = -0.5)
  call dwplot%init(xlabel = "$F(a/a_{\rm eq})$", ylabel = "$\delta w_\phi$", caption = "$V=C_0+C_1\phi^{-"//trim(coop_num2str(npower))//"}$",width=8., height=6.)

  do iQ=1,iQ_steps
     Q= 0.15d0*iQ
     c_run = Q * rhom_ini * phieq_ini ** (npower + 1.d0)/npower
     beta = 3.d0/(npower + 1.d0)

     args = coop_arguments( r = (/ Q, c_run, c_flat, npower /) )
     hubble_ini = sqrt((potential(phieq_ini, args) + rhom_ini)/(3.d0*Mpl**2))
     phieq_dot_ini = 3.d0*hubble_ini/(npower+1.d0)*phieq_ini

     hubble_dot_ini = -phieq_dot_ini**2.d0/2.d0 - rhom_ini/2.d0
     phieq_dot_dot_ini = 3.d0/(npower+1.d0)*(hubble_ini*phieq_dot_ini+hubble_dot_ini*phieq_ini)

     discriminant = d2Vdphi2(phieq_ini, args)**2.d0-2.d0*d3Vdphi3(phieq_ini, args)*(phieq_dot_dot_ini+3.d0*hubble_ini*phieq_dot_ini)
     if (discriminant .lt. 0) then
        stop "Discriminant is negative."

     else
        print *, "Discriminant is not negative."
        deltaphi = (-d2Vdphi2(phieq_ini, args)+sqrt(discriminant))/(d3Vdphi3(phieq_ini, args))
        print *, "phi_0= ", phieq_ini+deltaphi
     endif

     Ek = phieq_dot_ini**2.d0/2.d0
     Ev = potential(phieq_ini+deltaphi, args)
     print *, "Initially, Ek/Ev = ", Ek/Ev


     if(abs(deltaphi/phieq_ini).gt. 0.1d0)then
        write(*,*) "V''=", Vpp
        write(*,*) deltaphi/Mpl, phieq_ini/Mpl
        stop "deltaphi too big"
     endif
     call co%init(n = nvars, method = COOP_ODE_DVERK, tol = 1.d-8)  !!initialize the ode solver
     call co%set_arguments(args = args)
     call co%set_initial_conditions( xini = lna_start, yini = (/ phieq_ini+deltaphi, phieq_dot_ini, rhom_ini /) )
     call coop_set_uniform(nsteps, lna, lna_start, lna_end)
     itoday = 0
     ieq = 0
     do  i = 1, nsteps
        if(i.gt.1)call co%evolve(get_yprime, lna(i))
        phi(i) = co%PHI / Mpl
        dotphi(i) = co%PHIDOT
        Ek = dotphi(i)**2.d0/2.d0
        Ev = potential(phi(i), args)
        w(i) = (Ek -Ev)/(Ek+Ev)
        phieq(i) = get_phi_eq(co%y, args)
        dotphieq(i) = get_phi_eq_dot(co%y, args)
        dotdotphieq(i) = get_phi_eq_dot_dot(co%y, args)
        dotdotphi(i) = get_phi_dot_dot(co%y, args)

        !! Early time solution
        e_s(i) = 0.5d0/potential(phi(i), args)**2.d0*(dVdphi(phi(i), args)+Q*co%RHOM)**2.d0
        if((Ek+Ev)/co%RHOM .gt. (0.7/0.3) .and. itoday .eq. 0)then
           itoday = i
           print*, lna(itoday)
        endif
        if((Ek+Ev)/co%RHOM .gt. (0.5/0.5) .and. ieq .eq. 0) then
           ieq = i
           e_s_eq = e_s(i)
           print *, "epsilon_s_tilde = ", e_s_eq
        endif
     enddo
     lna  = lna - lna(itoday)  !!renormalize a=1 today
     a = exp(lna)
     F = exp(lna - lna(ieq))
     F = sqrt(1.d0+F**3)/F**1.5d0 - log(F**1.5d0+sqrt(1.d0+F**3))/F**3
     !! Compute w_slow(a) and delta_w(a)
     w_slow = 2./3.*e_s_eq*F**2.-1.
     delta_w = sqrt(1.+w) - sqrt(1.+w_slow)
     it = 1
     do while(a(it)*2.5 .lt. a(ieq))
        it = it + 1
     enddo
     eps_inf = (coop_sqrt3/coop_sqrt2 * delta_w(ieq)/(1.d0-coop_sqrt2 * F(ieq)))**2
     delta_w_ana = coop_sqrt2/coop_sqrt3 * sqrt(eps_inf)*(1.d0- coop_sqrt2 * F) 
     print*, "epsilon_s  = ", e_s_eq
     call coop_asy_interpolate_curve(wplot, xraw=exp(lna(1:itoday)), yraw=w(1:itoday),interpolate="LinearLinear", color = wplot%color(iQ), linewidth = 1.2, linetype="solid", legend = "$Q="//trim(coop_num2str(Q))//"$ numerical")  

     w_ana = -1 + 2./3.*(F*(sqrt(e_s_eq)-sqrt(eps_inf*2.))+sqrt(eps_inf))**2 &
          * tanh(a/a(it))

     call coop_asy_interpolate_curve(wplot, xraw=exp(lna(1:itoday)), yraw=w_ana(1:itoday),interpolate="LinearLinear", color = wplot%color(iQ), linewidth = 0.8, linetype="dotted", legend = "$Q="//trim(coop_num2str(Q))//"$ analytical")  

     call coop_asy_interpolate_curve(dwplot, xraw=F(1:itoday), yraw=delta_w(1:itoday),interpolate="LinearLinear", color = dwplot%color(iQ), linewidth = dwplot%linewidth(iQ), linetype=dwplot%linetype(iQ), legend = "$Q="//trim(coop_num2str(Q))//"$ numerical")  
     call coop_asy_interpolate_curve(dwplot, xraw=F(1:itoday), yraw=delta_w_ana(1:itoday),interpolate="LinearLinear", color = dwplot%color(iQ+iQ_steps), linewidth = dwplot%linewidth(iQ+iQ_steps), linetype=dwplot%linetype(iQ+iQ_steps), legend = "$Q="//trim(coop_num2str(Q))//"$ analytical")  

  enddo

  call coop_asy_legend(dwplot) !, x = (0.9*dwplot%xmin + 0.1*dwplot%xmax), y = dwplot%ymin*0.8+dwplot%ymax*0.2, cols=2)
  call coop_asy_legend(wplot) !, x = (0.9*wplot%xmin + 0.1*wplot%xmax), y = wplot%ymin*0.8+wplot%ymax*0.2, cols=2)

  call dwplot%close()
  call wplot%close()

contains


  function potential(phi, args) result(V)
    COOP_REAL phi, V
    type(coop_arguments) args
    V = C_RUN / phi**N_POWER + C_FLAT
  end function potential

  function dVdphi(phi, args) 
    COOP_REAL phi, dVdphi
    type(coop_arguments) args
    dVdphi = - C_RUN * N_POWER/phi**(N_POWER +1.d0)
  end function dVdphi

  function d2Vdphi2(phi, args)
    COOP_REAL phi, d2Vdphi2
    type(coop_arguments) args
    d2Vdphi2 = C_RUN * N_POWER * (N_POWER + 1.d0)/phi**(N_POWER+2.d0)
  end function d2Vdphi2

  function d3Vdphi3(phi, args)
    COOP_REAL phi, d3Vdphi3
    type(coop_arguments) args
    d3Vdphi3 = - C_RUN * N_POWER * (N_POWER + 1.d0) * (N_POWER + 2.d0)/phi**(N_POWER+3.d0)
  end function d3Vdphi3

  function get_H(y, args)
    COOP_REAL y(nvars), Ek, Ep, get_H
    type(coop_arguments) args
    Ek = PHIDOT ** 2 /2.d0
    Ep = potential(PHI, args)
    get_H = sqrt((Ek + Ep + RHOM)/(3.d0*MPl**2))
  end function get_H

  function get_H_dot(y, args)
    COOP_REAL y(nvars), Ek, get_H_dot
    type(coop_arguments) args
    Ek = PHIDOT ** 2/2.d0
    get_H_dot = -Ek - RHOM/2.d0
  end function get_H_dot

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

  function get_phi_eq(y, args)
    COOP_REAL y(nvars), get_phi_eq
    type(coop_arguments) args
    get_phi_eq = (N_POWER * C_RUN/ (Q_COUPLING*RHOM))**(1.d0/(N_POWER+1))
  end function get_phi_eq

  function get_phi_eq_dot(y, args)
    COOP_REAL y(nvars), get_phi_eq_dot
    type(coop_arguments) args
    get_phi_eq_dot = 3.d0/(N_POWER+1.d0)*get_H(y, args)*get_phi_eq(y, args)
  end function get_phi_eq_dot

  function get_phi_dot_dot(y, args)
    COOP_REAL y(nvars), get_phi_dot_dot
    type(coop_arguments) args
    get_phi_dot_dot = -3.d0*get_H(y, args)*PHIDOT-dVdphi(PHI, args)-Q_COUPLING*RHOM
  end function get_phi_dot_dot

  function get_phi_eq_dot_dot(y, args)
    COOP_REAL y(nvars), get_phi_eq_dot_dot
    type(coop_arguments) args
    get_phi_eq_dot_dot = 3.d0/(N_POWER+1.d0)*(get_H(y, args)*get_phi_eq_dot(y, args)+get_H_dot(y, args)*get_phi_eq(y, args))
  end function get_phi_eq_dot_dot

end program combined2
