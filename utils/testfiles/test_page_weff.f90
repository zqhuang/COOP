program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"

  COOP_REAL::eta, t0
  COOP_INT,parameter::n=128
  COOP_REAL::weff(n), z(n), a(n), oml(n)
  COOP_INT::i
  type(coop_asy)::fig
  
  t0 =  0.9 !wcdm_tofa(omegam, omegak, w, 1.d0)
  eta = -0.4 !1.d0 -  1.5 * (omegam/2 + (1.+3.*w)/2. * (1.-omegam-omegak) + 1.) * t0**2
  call fig%open("wde.txt")
  call fig%init(xlabel = "$z$", ylabel = "$w_{\rm de}$")  
  call coop_set_uniform(n, z, 0.d0, 5.d0)
  a = 1.d0/(1.d0+z)
  do i=1, n
     weff(i) = page_effective_wdeofa(t0, eta, a(i))
     oml(i) = 1.d0 -page_effective_ommofa(t0, eta, a(i))
     write(*,*) z(i), oml(i), weff(i), page_Hofa(t0,eta,a(i))**2*page_effective_ommofa(t0, eta, a(i))*a(i)**3
  enddo
  call fig%plot(z, weff)
  call fig%close()

contains

  function ScalarPower(k)
    COOP_INT,parameter::dl = kind(1.d0)
    real(dl),parameter::As = 2.4e-9, ns = 0.967
    real(dl), intent(in) :: k
    real(dl) ScalarPower
    real(dl) lnrat
    !!added for features
    real(dl),dimension(3), parameter::curl_A = (/ 0.01d0 , 0.d0, 0.d0 /)
    real(dl),parameter::curl_xs = 1.0
    real(dl),parameter::curl_ks = 0.05 !!Mpc^{-1}
    real(dl)::krat, curl_I0, curl_I1, lnP0, curl_D
    real(dl),dimension(3)::curl_W, curl_V
    !ScalarPower = const for scale invariant spectrum
    !The normalization is defined so that for adiabatic perturbations the gradient of the 3-Ricci
    !scalar on co-moving hypersurfaces receives power
    ! < |D_a R^{(3)}|^2 > = int dk/k 16 k^6/S^6 (1-3K/k^2)^2 ScalarPower(k)
    !In other words ScalarPower is the power spectrum of the conserved curvature perturbation given by
    !-chi = \Phi + 2/3*\Omega^{-1} \frac{H^{-1}\Phi' - \Psi}{1+w}
    !(w=p/\rho), so < |\chi(x)|^2 > = \int dk/k ScalarPower(k).
    !Near the end of inflation chi is equal to 3/2 Psi.
    !Here nu^2 = (k^2 + curv)/|curv|

    !This power spectrum is also used for isocurvature modes where
    !< |\Delta(x)|^2 > = \int dk/k ScalarPower(k)
    !For the isocurvture velocity mode ScalarPower is the power in the neutrino heat flux.

    
    lnrat = log(k/0.05d0)
    lnP0 = lnrat * (this%ns - 1)
    krat = k/curl_ks
    curl_D = damp(krat/curl_xs)
    curl_W = (/ W1ofx(krat), W2ofx(krat), W3ofx(krat) /)
    curl_V =(/  V1ofx(krat), V2ofx(krat), V3ofx(krat) /)
    curl_I0 = sum(curl_A * curl_W) * curl_D
    curl_I1 = (1.5708*(1.0-this%ns)+sum(curl_A * curl_V)*curl_D)/sqrt(2.d0)
    ScalarPower = As * exp(lnP0 + curl_I0 + log(1.0+curl_I1 ** 2))

    
  contains

    function damp(x)
      real(dl)::x, damp
      if(x.lt. 1.d-3)then
         damp = 1.d0 - x*x/6.d0
      else
         damp = x/sinh(x)
      endif
    end function damp

    function W1ofx(x) result(W1)
      real(dl)::x, W1
      if(x .lt. 1.d-3)then
         W1 = ((4.d0/5.d0) - (24.d0/35.d0)*x**2) * x**2
      else
         W1 = (x*(18.d0-6.d0*x**2)*cos(2.d0*x)+(15.d0*x**2-9.d0)*sin(2.d0*x))/(2.d0*x**3)
      endif
    end function W1ofx

    function W2ofx(x) result(W2)
      real(dl)::x, W2
      if(x .lt. 1.d-3)then
         W2 = 1.d0 + (2.d0/5.d0) * x**2
      else
         W2 = 3.d0/(2.d0*x**3)*((1.d0-x**2)*sin(2.d0*x)-2.d0*x*cos(2.d0*x))
      endif
    end function W2ofx

    function W3ofx(x) result(W3)
      real(dl)::x, W3
      if(x .lt. 1.d-3)then
         W3 = x**2*(-32.d0/15.d0 + x**2*(64.d0/105.d0))
      else
         W3 = (6*x*cos(2*x)+(4*x**2-3)*sin(2*x))/x**3
      endif
    end function W3ofx

    function V1ofx(x) result(V1)
      real(dl)::x, V1
      if(x .lt. 1.d-3)then
         V1 = x**3*(1.d0-x**2/3.d0)
      else
         V1 = (-3.d0/x**3)*(x*cos(x)-sin(x))*(3.d0*x*cos(x)+(2.d0*x**2-3.d0)*sin(x))
      endif
    end function V1ofx

    function V2ofx(x) result(V2)
      real(dl)::x, V2
      if(x .lt. 1.d-3)then
         V2 = x**3*(1.d0/3.d0 - x**2/15.d0)
      else
         V2 = (sin(x)-x*cos(x))**2*(3/x**3)
      endif
    end function V2ofx

    function V3ofx(x) result(V3)
      real(dl)::x, V3
      if(x .lt. 1.d-3)then
         V3 = x*(2.d0-(4.d0/3.d0)*x**2)
      else
         V3 = ((3-4*x**2)*cos(2*x)+6*x*sin(2*x)-3-2*x**2)/x**3
      endif
    end function V3ofx
    
  end function TInitialPowerLaw_ScalarPower

  
end program Test  
