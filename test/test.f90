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
#define N_POWER args%r(4)

!!example for coupled quintessence model 
!!V(phi) = c0 / phi^n + c1
!!I am not doing perfect initial conditions, you might get some initial small oscillations

  COOP_INT, parameter::nvars = 3  !!number of variables to evolve
  COOP_INT, parameter::nsteps = 1000  !!number of steps for output
  COOP_REAL, parameter:: MPl = 1.d0 !!define mass unit
  COOP_REAL, parameter:: H0 = 1.d0  !!define time unit
  COOP_REAL,parameter::rho_total_today = 3.d0*Mpl**2*H0**2
  COOP_REAL,parameter::rhom_today = 0.3*rho_total_today
  COOP_REAL:: Q, c_run
  COOP_REAL, parameter::npower = 1.
  COOP_REAL,parameter:: rhom_ini = 1.e9 * rhom_today
  COOP_REAL,parameter:: phieq_ini = 0.0001 * Mpl  
  COOP_REAL, parameter::c_flat = 0.7* rho_total_today
  COOP_REAL, parameter::lna_start = 0.
  COOP_REAL, parameter::lna_end = log(rhom_ini/rhom_today)/3.d0 + 2.d0  !!add  3 more efolds to make sure rho_m < rho_phi happens
  COOP_INT:: itoday, iplotto

  COOP_INT:: iQ

  type(coop_ode)::co
  type(coop_arguments)::args
  COOP_INT i
  COOP_REAL::phieq_dot_ini, hubble_ini
  COOP_REAL:: phi(nsteps), lna(nsteps), phidot(nsteps), w(nsteps), phieq(nsteps), phieq_dot(nsteps), hubble(nsteps), Ev, Ek, Vpp, deltaphi, deltaphi_dot
  type(coop_asy)::plot_w, plot_phi, plot_phidot

  call plot_w%open("w.txt")
  call plot_phi%open("phi.txt")
  call plot_phidot%open("phidot.txt")
  
  call plot_w%init(xlabel = "$\ln a$", ylabel = "$1+w_\phi$", caption = "$\phi_{\rm ini} = "//trim(coop_num2str(phieq_ini/Mpl))//"M_p$")  
  call plot_phi%init(xlabel = "$\ln a$", ylabel = "$\phi/M_p$", caption = "$\phi_{\rm ini} = "//trim(coop_num2str(phieq_ini/Mpl))//"M_p$") 
  call plot_phidot%init(xlabel = "$\ln a$", ylabel = "$\dot\phi/(H_0M_p)$", caption = "$\phi_{\rm ini} = "//trim(coop_num2str(phieq_ini/Mpl))//"M_p$" )
     

  do iQ = 1, 3
     Q = iQ * 0.15
     c_run = Q * rhom_ini * phieq_ini ** (npower + 1)/npower
     args = coop_arguments( r = (/ Q, c_run, c_flat, npower /) )
     hubble_ini = sqrt((potential(phieq_ini, args) + rhom_ini)/(3.d0*Mpl**2))
     call get_equil(hubble_ini, phieq_ini, phieq_dot_ini, deltaphi, args)
     if(abs(deltaphi/phieq_ini).gt. 0.1d0)then
        write(*,*) "V''=", Vpp
        write(*,*) deltaphi/Mpl, phieq_ini/Mpl
        stop "deltaphi too big"
     else 
        write(*,*) "Correction: ", deltaphi/phieq_ini
        phieq_dot_ini = phieq_dot_ini * (1.d0 + deltaphi/phieq_ini)  !!assuming the ratio is varying slowly
     endif
     call co%init(n = nvars, method = COOP_ODE_DVERK)  !!initialize the ode solver
     call co%set_arguments(args = args)
     call co%set_initial_conditions( xini = lna_start, yini = (/ phieq_ini+deltaphi, phieq_dot_ini, rhom_ini /) )
     call coop_set_uniform(nsteps, lna, lna_start, lna_end)
     itoday = 0
     do  i = 1, nsteps
        if(i.gt.1)call co%evolve(get_yprime, lna(i))
        phi(i) = co%PHI / Mpl
        phidot(i) = co%PHIDOT
        Ek = phidot(i)**2/2.d0
        Ev = potential(phi(i), args)
        w(i) = (Ek -Ev)/(Ek+Ev)
        hubble(i) = sqrt((Ek+Ev+co%RHOM)/(3.d0*Mpl**2))
        phieq(i) = (co%C_RUN * co%N_POWER / co%Q_COUPLING/co%RHOM)**(1.d0/(CO%N_POWER + 1.d0))
        phieq_dot(i) = phieq(i)*(3.*hubble(i)/(co%N_POWER+1.d0))
        if((Ek+Ev)/co%RHOM .ge. (0.7/0.3) .and. itoday.eq.0)then
           itoday = i
           exit !!you don't really need later solutions
        endif
     enddo
     if(itoday .eq. 0)then
        write(*,*) "you might want to increase lna_end"
        stop
     endif
     lna = lna - lna(itoday)  !!renormalize a=1 today
     iplotto = itoday
     select case(iQ)
     case(1)
        call coop_asy_curve(plot_w, x=lna(1:iplotto), y=(1.d0+w(1:iplotto)), &
             color = "black", linewidth = 1.5, linetype="solid", &
             legend = "$Q = "//trim(coop_num2str(co%Q_COUPLING))//"$" )
        call coop_asy_curve(plot_phi, x=lna(1:iplotto), y=phi(1:iplotto)/Mpl, &
             color = "black", linewidth = 1.5, linetype="solid", &
             legend = "$\phi, Q = "//trim(coop_num2str(co%Q_COUPLING))//"$" )
        call coop_asy_curve(plot_phi, x=lna(1:iplotto), y=phieq(1:iplotto)/Mpl, &
             color = "black", linewidth = 0.5, linetype="solid", &
             legend = "$\phi_{\rm eq}, Q = "//trim(coop_num2str(co%Q_COUPLING)) //"$")
        call coop_asy_curve(plot_phidot, x=lna(1:iplotto), y=phidot(1:iplotto)/Mpl/H0, &
             color = "black", linewidth = 1.5, linetype="solid", &
             legend = "$\dot\phi, Q = "//trim(coop_num2str(co%Q_COUPLING))//"$" )
        call coop_asy_curve(plot_phidot, x=lna(1:iplotto), y=phieq_dot(1:iplotto)/Mpl/H0, &
             color = "black", linewidth = 0.5, linetype="solid", &
             legend = "$\dot\phi_{\rm eq}, Q = "//trim(coop_num2str(co%Q_COUPLING))//"$" )
     case(2)
        call coop_asy_curve(plot_w, x=lna(1:iplotto), y=(1.d0+w(1:iplotto)), &
             color = "red", linewidth = 1.5, linetype="dotted", &
             legend = "$Q = "//trim(coop_num2str(co%Q_COUPLING))//"$" )
        call coop_asy_curve(plot_phi, x=lna(1:iplotto), y=phi(1:iplotto)/Mpl, &
             color = "red", linewidth = 1.5, linetype="dotted", &
             legend = "$\phi, Q = "//trim(coop_num2str(co%Q_COUPLING))//"$" )
        call coop_asy_curve(plot_phi, x=lna(1:iplotto), y=phieq(1:iplotto)/Mpl, &
             color = "red", linewidth = 0.5, linetype="dotted", &
             legend = "$\phi_{\rm eq}, Q = "//trim(coop_num2str(co%Q_COUPLING))//"$" )
        call coop_asy_curve(plot_phidot, x=lna(1:iplotto), y=phidot(1:iplotto)/Mpl, &
             color = "red", linewidth = 1.5, linetype="dotted", &
             legend = "$\dot\phi, Q = "//trim(coop_num2str(co%Q_COUPLING))//"$" )
        call coop_asy_curve(plot_phidot, x=lna(1:iplotto), y=phieq_dot(1:iplotto)/Mpl, &
             color = "red", linewidth = 0.5, linetype="dotted", &
             legend = "$\dot\phi_{\rm eq}, Q = "//trim(coop_num2str(co%Q_COUPLING))//"$" )

     case(3)
        call coop_asy_curve(plot_w, x=lna(1:iplotto), y=(1.d0+w(1:iplotto)), &
             color = "blue", linewidth = 1.5, linetype="dashed", &
             legend = "$Q = "//trim(coop_num2str(co%Q_COUPLING))//"$" )
        call coop_asy_curve(plot_phi, x=lna(1:iplotto), y=phi(1:iplotto), &
             color = "blue", linewidth = 1.5, linetype="dashed", &
             legend = "$\phi, Q = "//trim(coop_num2str(co%Q_COUPLING))//"$" )
        call coop_asy_curve(plot_phi, x=lna(1:iplotto), y=phieq(1:iplotto), &
             color = "blue", linewidth = 0.5, linetype="dashed", &
             legend = "$\phi_{\rm eq}, Q = "//trim(coop_num2str(co%Q_COUPLING))//"$" )
        call coop_asy_curve(plot_phidot, x=lna(1:iplotto), y=phidot(1:iplotto), &
             color = "blue", linewidth = 1.5, linetype="dashed", &
             legend = "$\dot\phi, Q = "//trim(coop_num2str(co%Q_COUPLING))//"$" )
        call coop_asy_curve(plot_phidot, x=lna(1:iplotto), y=phieq_dot(1:iplotto), &
             color = "blue", linewidth = 0.5, linetype="dashed", &
             legend = "$\dot\phi_{\rm eq}, Q = "//trim(coop_num2str(co%Q_COUPLING))//"$" )


     end select

  enddo
  call coop_asy_legend(plot_w, x = plot_w%xmin*0.6+plot_w%xmax*0.4, y = plot_w%ymin*0.8 + plot_w%ymax*0.2)
  call coop_asy_legend(plot_phi, x = plot_phi%xmin*0.8+plot_phi%xmax*0.2, y =plot_phi%ymin*0.2 + plot_phi%ymax*0.8 )
  call coop_asy_legend(plot_phidot, x = plot_phidot%xmin*0.8+plot_phidot%xmax*0.2, y = plot_phidot%ymin*0.1 +  plot_phidot%ymax*0.9 )

  call coop_asy_label(plot_w, label = "$V(\phi) = C_0 + C_1 \phi^{-"//trim(coop_num2str(co%N_POWER))//"}$", x = plot_w%xmin*0.8+plot_w%xmax*0.2, y = plot_w%ymax*0.9 + plot_w%ymin*0.1)

  call plot_w%close()
  call plot_phi%close()
  call plot_phidot%close()

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

  function d2Vdphi2(phi, args)
    COOP_REAL phi, d2Vdphi2, dphi
    type(coop_arguments) args
    dphi = phi/100.d0
    d2Vdphi2 = (dVdphi(phi+dphi, args) - dVdphi(phi-dphi, args))/(2.d0*dphi)
  end function d2Vdphi2

  subroutine get_equil(H, phieq, phieq_dot, deltaphi, args)
    COOP_REAL phieq, phieq_dot, deltaphi, deltaphi_dot, H, Vpp
    type(coop_arguments) args
    phieq_dot = 3.d0*H/(N_POWER + 1.d0)*phieq
    Vpp = d2Vdphi2(phieq, args)
    deltaphi = -3.d0*H*phieq_dot/Vpp
  end subroutine get_equil

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
