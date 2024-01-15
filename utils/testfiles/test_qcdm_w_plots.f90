module tmp
#include "constants.h"      
  use coop_wrapper_utils
  implicit none  
  COOP_REAL,parameter::omegam = 0.31d0
  COOP_REAL,parameter::h = 0.68
  COOP_REAL,parameter::Mpl = 1.d0
  COOP_REAL,parameter::H0 = 1.d0
  COOP_REAL,parameter::rho_crit0 = 3.d0*Mpl**2*H0**2  
  COOP_REAL,parameter::omegar = 4.17e-5/h**2 !!calc "2.726**4*4*3.14159**4/45*1.381e-3**4/6.626e-4**2*3.0857e7**2/2.9979e8**4/4.341e-9**2*(1+7./8.*(4./11.)**(4./3.)*3.)"
  COOP_REAL,parameter::omegak = 0.1
  COOP_REAL,parameter::omegaphi = 1.d0-omegam-omegak-omegar
  COOP_REAL,parameter::aeq = (omegam/omegaphi)**(1.d0/3.d0)
  COOP_REAL::a_eq = aeq
  COOP_REAL::lnaeq = log(aeq)
  
  COOP_REAL::epsilon_s=0. !!TBD
  COOP_REAL::epsilon_inf=0. !!TBD
  COOP_REAL::zeta_s=0. !!TBD
  COOP_REAL::alpha  !!TBD
  COOP_REAL,parameter::beta = -1.d0
  COOP_REAL,parameter::lambda = 1.d0

#define EXPONENTIAL 0
#define POWERLAW 1
#define NPOW  2  

#define POTENTIAL_TYPE NPOW
  
contains

  function qwf(x, lam)
    COOP_REAL,parameter::xmax = 7.d0
    COOP_REAL::x, dx, lam, qwf, xx, twodx
    COOP_INT::i, n
    if(x.lt.0.02d0)then
       qwf = 3.d0*sqrt(x**3) * ( &
            2.d0/9.d0 + x*( &
            -lam/11.d0+x*(3.d0/52.d0*lam**2)))
       return
    endif
    if(x.gt.xmax)then
       n = ceiling(xmax/0.3d0) + 8
       dx = xmax/(n*2)
       xx = dx
       twodx = dx*2.d0
       qwf = integrand(xx)*2.d0
       do i = 1, n-1
          xx = xx + twodx       
          qwf = qwf + integrand(xx)*2.d0 + integrand(xx-dx)
       enddo
       qwf = ((qwf*2.d0 + integrand(xmax))* dx - (lam*(x-xmax)+log(x/xmax))*1.5d0)/x**3 + (1.d0 - (xmax/x)**3)
    else
       n = ceiling(x/0.3d0) + 8
       dx = x/(n*2)
       xx = dx
       twodx = dx*2.d0
       qwf = integrand(xx)*2.d0
       do i = 1, n-1
          xx = xx + twodx       
          qwf = qwf + integrand(xx)*2.d0 + integrand(xx-dx)
       enddo
       qwf = (qwf*2.d0 + integrand(x))*(dx/x**3)
    endif
    
  contains
    
    function integrand(x)
      COOP_REAL::integrand, x
      integrand = x**3*sqrt(x/(1+x**3+lam*x))
    end function integrand
    
  end function qwf

  function wp1_simple(a) 
    COOP_REAL::a, wp1_simple
    wp1_simple = (2.d0/3.d0)*epsilon_s*qwf(a/aeq, omegak/omegam)**2
  end function wp1_simple


  function wp1(a) 
    COOP_REAL::a, wp1
    COOP_REAL:: epss
    COOP_REAL::mu, mu3, sqrtepss, sqrtepsinf, diff, delta, f, s1, f2, s0,  qpsign, lam
    if(epsilon_s .ge. 0.d0)then
       qpsign = 1.d0
       epss = epsilon_s
    else
       qpsign = -1.d0
       epss = - epsilon_s
    endif
    sqrtepss = sqrt(epss)
    sqrtepsinf = sqrt(epsilon_inf)
    diff = sqrtepss - sqrt(2./(1.-omegak)) * sqrtepsinf
    delta = (sqrtepsinf + (0.91-0.78*omegam/(1.-omegak)+(0.236-0.76*omegam/(1.-omegak))*zeta_s)*diff)**2 &
         + (sqrtepsinf + (0.533-0.1*zeta_s)*diff)**2 
    a_eq = (omegam/omegaphi)**(1.d0/(3.d0-qpsign*delta))
    mu = a/a_eq
    mu3 = mu**3    
    s0 = sqrt(mu3)
    if(mu .gt. 0.05d0)then
       s1 = sqrt(1.d0+mu3)
       f = s1/s0 - log(s0+s1)/mu3       
       f2 = coop_sqrt2*(1.d0-log(1.d0+mu3)/mu3) - f
    else  !!asymptotic
       f = s0*((2.d0/3.d0)-0.2d0*mu3)       
       f2 = mu3*(1.d0/coop_sqrt2 - (coop_sqrt2/3.d0)*mu3)  - f
    endif
    s0 = sqrtepsinf*sqrt(((4.d0/3.d0*omegar) + omegam*a)/(omegar+omegam*a))
    lam = omegak/omegam*a_eq
    wp1 = (2.d0/3.d0)*qpsign*(s0 + (sqrtepss  -sqrt(2./(1.-omegak))*s0)*(qwf(mu, lam) + zeta_s * f2))**2
   ! wp1 = (2.d0/3.d0)*qpsign*(s0 + (sqrtepss  - coop_sqrt2*s0)*(f + zeta_s * f2))**2
  end function wp1


  function potential(phi) result(V)
    COOP_REAL::V, phi
#if POTENTIAL_TYPE == EXPONENTIAL
    V = alpha*exp(lambda*phi/Mpl)
#elif POTENTIAL_TYPE == NPOW || POTENTIAL_TYPE == POWERLAW
    V = alpha*(phi/Mpl)**beta
#endif    
  end function potential

  function dVdphi(phi)
    COOP_REAL::dVdphi, phi
#if POTENTIAL_TYPE == EXPONENTIAL
    dVdphi = (alpha*lambda/Mpl)*exp(lambda*phi/Mpl)
#elif POTENTIAL_TYPE == NPOW || POTENTIAL_TYPE == POWERLAW
    dVdphi = (alpha*beta/Mpl)*(phi/Mpl)**(beta-1.d0)
#endif
  end function dVdphi

  !!y(1) \ln a 
  !!y(3) \phi
  !!y(4) \dot\phi
#define LNA y(1)  
#define PHI y(2)
#define DOTPHI y(3)
#define HUBBLE yp(1)
#define PHI_PRIME yp(2)
#define DOTPHI_PRIME yp(3)  
  subroutine scalar_field_evolve(n, t, y, yp)
    COOP_INT::n
    COOP_REAL::t, y(n), yp(n), a
    a = exp(LNA)
    HUBBLE = sqrt(rho_crit0*((omegar/a + omegam)/a+omegak)/a**2 + DOTPHI**2/2.d0 + potential(PHI))/(coop_sqrt3*Mpl)
    PHI_PRIME = DOTPHI
    DOTPHI_PRIME=  - 3.d0 * HUBBLE * DOTPHI - dVdphi(PHI)
  end subroutine scalar_field_evolve



end module tmp



program Test
#include "constants.h"    
  use coop_wrapper_utils
  use tmp
  implicit none
  type(coop_ode)::ode
  COOP_REAL::tini, yini(3), aini, dt, tend, phi_ini, Hini, dotphi_ini
  COOP_REAL::up, down, mid, Hup, Hdown, Hmid
  type(coop_file)::output
  COOP_INT::i
  COOP_STRING::fname
  fname = trim(coop_InputArgs(1))
  tini = 0.d0
  aini = 5.d-6
  alpha = Mpl**2*H0**2
#if POTENTIAL_TYPE == EXPONENTIAL
  phi_ini = 0.d0
#elif  POTENTIAL_TYPE == POWERLAW
  phi_ini = 9.d0**(1.d0/beta)*Mpl  
#elif POTENTIAL_TYPE == NPOW 
  phi_ini = 1.d-10*Mpl
#endif
  Hini = sqrt((rho_crit0*((omegar/aini+omegam)/aini+omegak)/aini**2 + potential(phi_ini))/(3.d0*Mpl**2))
  dotphi_ini = - dVdphi(phi_ini)/(3*Hini)  
  do while(dotphi_ini**2 .gt. potential(phi_ini)*0.1)  !!adjust phi_ini to a reasonable range
     phi_ini = phi_ini * 2.d0     
     Hini = sqrt((rho_crit0*((omegar/aini+omegam)/aini+omegak)/aini**2 + potential(phi_ini)+dotphi_ini**2/2.d0)/(3.d0*Mpl**2))
     dotphi_ini = - dVdphi(phi_ini)/(3*Hini)
  enddo

  call setnorm()
  print*, "w(z=0) = ", wp1(1.d0)-1.d0 !!do not delete this line; it is used to set up a_eq
  print*, "a_eq = ", a_eq 
  lnaeq = log(a_eq) !!now use the more accurate a_eq to set epsilon_s etc.
  call getH(mid)
  if(trim(fname).ne.'')then
     call output%open(fname)
  else
     call output%open("woflna.txt")
  endif
  write(output%unit, "(A4, I6, 4G14.5)") "#   ", POTENTIAL_TYPE , alpha/Mpl**2/H0**2/3.d0, beta, lambda, omegak
  write(output%unit, "(A4, 2G14.5)") "#   ", phi_ini/Mpl, dotphi_ini/Mpl**2
  write(output%unit, "(A4, 4G14.5)") "#   ", epsilon_s, epsilon_inf, zeta_s, a_eq  
  call getH(mid, output)
  call output%close()
contains


  subroutine setnorm()
    yini = (/ log(aini), phi_ini, dotphi_ini /)
    up = 5.d0
    down = 0.2d0
    Hup = 0.9*H0
    do while(Hup .lt. H0)
       up = up*2.d0
       call getH(up)
       Hup = ode%HUBBLE
    enddo
    Hdown = 1.1*H0
    do while(Hdown .gt. H0)
       down = down/2.d0
       call getH(down)
       Hdown = ode%HUBBLE
    enddo
    do while(Hup - Hdown .gt. 5.d-2*H0)
       mid = (up+down)/2.d0
       call getH(mid)
       Hmid = ode%HUBBLE
       write(*,*) mid, Hmid/H0-1.d0, epsilon_s
       if(Hmid .gt. H0)then
          up = mid
          Hup = Hmid
       else
          down = mid
          Hdown = Hmid
       endif
    enddo
    do while(abs(Hmid/H0-1.d0) .gt. 1.d-8)
       mid = (up*(H0-Hdown) + down*(Hup-H0))/(Hup - Hdown)
       call getH(mid)
       Hmid = ode%HUBBLE
       write(*,*) mid, Hmid/H0-1.d0, epsilon_s
       if(Hmid .gt. H0)then
          up = mid
          Hup = Hmid
       else
          down = mid
          Hdown = Hmid
       endif
    end do
  end subroutine setnorm

  subroutine getH(norm, fp)
    COOP_REAL,parameter::am1 = 1.d0/51.d0
    COOP_INT::i
    COOP_REAL::norm
    type(coop_file),optional::fp
    COOP_REAL::phi_up, phi_down, lna_up, lna_down, phieq, wm1, lastdm1, dm1, wp1_up, wp1_down, wp1eq, goodzs, minwdiff, wdiff
    alpha = norm*Mpl**2*H0**2
    call ode%init(n=3, method=COOP_ODE_DVERK, tol = 1.d-6)
    call ode%set_initial_conditions(tini, yini)
    call scalar_field_evolve(3, ode%x, ode%y, ode%yp)  
    if(present(fp))write(fp%unit,"(3G14.5)") exp(ode%LNA), (ode%DOTPHI**2/2.d0-potential(ode%PHI))/(ode%DOTPHI**2/2.d0+potential(ode%PHI)), wp1(exp(ode%LNA))-1.d0
    lna_down = 100.
    lastdm1 = -100.
    do while(ode%y(1) .lt. -0.1d0)
       dt = 0.1/ode%HUBBLE
       tend = ode%x + dt
       call ode%evolve(scalar_field_evolve, tend)
       call scalar_field_evolve(3, ode%x, ode%y, ode%yp)
       if(present(fp))then
          write(fp%unit,"(3G14.5)") exp(ode%LNA), (ode%DOTPHI**2/2.d0-potential(ode%PHI))/(ode%DOTPHI**2/2.d0+potential(ode%PHI)), wp1(exp(ode%LNA))-1.d0
       else
#if   POTENTIAL_TYPE == NPOW           
          if(lastdm1 .lt. 0.)then
             dm1 = ode%LNA - log(am1)
             if(dm1 .ge. 0.)then
                epsilon_inf = 1.5d0* ode%DOTPHI**2/(ode%DOTPHI**2/2.d0+potential(ode%PHI))
             endif
             lastdm1 = dm1             
          endif
#endif          
          if(lna_down .gt. 99.)then
             if(ode%LNA .lt. lnaeq)then
                lna_up = ode%LNA
                phi_up = ode%PHI
                wp1_up = ode%DOTPHI**2 /(ode%DOTPHI**2/2.d0+potential(ode%PHI))                
             else
                lna_down = ode%LNA
                phi_down = ode%PHI
                wp1_down = ode%DOTPHI**2 /(ode%DOTPHI**2/2.d0+potential(ode%PHI))                             
             endif
          endif
       endif
    enddo
    do i=1, 3
       dt = -ode%LNA/ode%HUBBLE
       tend = ode%x + dt
       call ode%evolve(scalar_field_evolve, tend)
       call scalar_field_evolve(3, ode%x, ode%y, ode%yp)
       if(present(fp))then
          write(fp%unit,"(3G14.5)") exp(ode%LNA), (ode%DOTPHI**2/2.d0-potential(ode%PHI))/(ode%DOTPHI**2/2.d0+potential(ode%PHI)), wp1(exp(ode%LNA))-1.d0
       elseif(lna_down .gt. -99.)then       
          if(ode%LNA .gt. lnaeq)then
             lna_up = ode%LNA
             phi_up = ode%PHI
             wp1_up = ode%DOTPHI**2 /(ode%DOTPHI**2/2.d0+potential(ode%PHI))
          else
             lna_down = ode%LNA
             phi_down = ode%PHI
             wp1_down = ode%DOTPHI**2 /(ode%DOTPHI**2/2.d0+potential(ode%PHI))             
          endif
       endif       
    enddo    
    if(.not. present(fp))then
       phieq = (phi_up * (lnaeq - lna_down) + phi_down * (lna_up - lnaeq))/(lna_up - lna_down)
       wp1eq = (wp1_up * (lnaeq - lna_down) + wp1_down * (lna_up - lnaeq))/(lna_up - lna_down)       
       epsilon_s = (dVdphi(phieq)/potential(phieq)*Mpl)**2/2.d0       
#if POTENTIAL_TYPE == EXPONENTIAL || POTENTIAL_TYPE == POWERLAW 
       epsilon_inf = 0.d0
#endif
       minwdiff = 1.d99
       do i=-100, 100
          zeta_s = i*0.01d0
          wdiff = abs(wp1(a_eq)-wp1eq)
          if(wdiff .lt. minwdiff)then
             goodzs = zeta_s
             minwdiff = wdiff
          endif
       enddo
       zeta_s = goodzs
    endif
  end subroutine getH
    
end program Test



