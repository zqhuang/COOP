Module coop_page_mod
  implicit none
  !!a: scale factor, normalized to 1 today
  !!t: H0 * age_then
  !!z: redshift, 1/a-1
  !!H: H(z)/H_0
  !!dL: H_0 * luminosity distance
  !!dA: H_0 * angular diameter distance
  !!dc: H0 * comoving distance
  !!t0: H0* age_now

  !! We work with the H0 = 1 units.
  !!model: H = 2/(3t) + H_0 - 2/(3t0) +  2 eta/(3t0)  (t/t0-1) 
  !!The advantage of this parametrization is that we more or less know that
  !!      t0 ~ 1  (in LCDM it is (0.2629/Omega_m)**0.277
  !!   d H / dt < 0 => eta < 1
  !!      q0 = 2(1-eta)/(3t0^2) - 1
  !!      j0 =  4/(3t0^3) - 3q0 - 2
contains


  !!====================General Applications =============================
  function zerop_integrate(func,  a, b, step) result(integral)
    external func !! func(t) returns real*8
    real*8:: func, a, b, step, integral, sumodd
    real*8::t, dt
    integer::n, i
    n = floor((b-a)/step)+2
    dt  = (b-a)/n/2.d0
    t = a
    sumodd =  0.d0
    integral = (func( a)+func( b))/2.d0
    do i=1, n-1
       t  = t + dt
       sumodd = sumodd + func( t)
       t = t + dt
       integral = integral + func(t)
    enddo
    integral = ((sumodd + func(b-dt))*2.d0 + integral)*(dt/1.5d0)
  end function zerop_integrate


  
  function onep_integrate(func, param1, a, b, step) result(integral)
    external func !! func(param1, t) returns real*8
    real*8::param1, func, a, b, step, integral, sumodd
    real*8::t, dt
    integer::n, i
    n = floor((b-a)/step)+2
    dt  = (b-a)/n/2.d0
    t = a
    sumodd =  0.d0
    integral = (func(param1, a)+func(param1,  b))/2.d0
    do i=1, n-1
       t  = t + dt
       sumodd = sumodd + func(param1,  t)
       t = t + dt
       integral = integral + func(param1, t)
    enddo
    integral = ((sumodd + func(param1, b-dt))*2.d0 + integral)*(dt/1.5d0)
  end function onep_integrate


  
  function twop_integrate(func, param1, param2, a, b, step) result(integral)
    external func !! func(param1, param2, t) returns real*8
    real*8::param1,param2,func, a, b, step, integral, sumodd
    real*8::t, dt
    integer::n, i
    n = floor((b-a)/step)+2
    dt  = (b-a)/n/2.d0
    t = a
    sumodd =  0.d0
    integral = (func(param1, param2, a)+func(param1, param2, b))/2.d0
    do i=1, n-1
       t  = t + dt
       sumodd = sumodd + func(param1, param2, t)
       t = t + dt
       integral = integral + func(param1, param2, t)
    enddo
    integral = ((sumodd + func(param1, param2, b-dt))*2.d0 + integral)*(dt/1.5d0)
  end function twop_integrate

  function threep_integrate(func, param1, param2, param3, a, b, step) result(integral)
    external func !! func(param1, param2, param3, t) returns real*8
    real*8::param1,param2,param3, func, a, b, step, integral, sumodd
    real*8::t, dt
    integer::n, i
    n = floor((b-a)/step)+2
    dt  = (b-a)/n/2.d0
    t = a
    sumodd =  0.d0
    integral = (func(param1, param2, param3, a)+func(param1, param2, param3, b))/2.d0
    do i=1, n-1
       t  = t + dt
       sumodd = sumodd + func(param1, param2, param3, t)
       t = t + dt
       integral = integral + func(param1, param2, param3, t)
    enddo
    integral = ((sumodd + func(param1, param2, param3, b-dt))*2.d0 + integral)*(dt/1.5d0)
  end function threep_integrate


  function fourp_integrate(func, param1, param2, param3, param4, a, b, step) result(integral)
    external func !! func(param1, param2, param3, param4, t) returns real*8
    real*8::param1,param2,param3,param4, func, a, b, step, integral, sumodd
    real*8::t, dt
    integer::n, i
    n = floor((b-a)/step)+2
    dt  = (b-a)/n/2.d0
    t = a
    sumodd =  0.d0
    integral = (func(param1, param2, param3, param4,a)+func(param1, param2, param3,param4, b))/2.d0
    do i=1, n-1
       t  = t + dt
       sumodd = sumodd + func(param1, param2, param3,param4, t)
       t = t + dt
       integral = integral + func(param1, param2, param3,param4, t)
    enddo
    integral = ((sumodd + func(param1, param2, param3,param4, b-dt))*2.d0 + integral)*(dt/1.5d0)
  end function fourp_integrate
  
 
  function curved_rofchi(chi, omegak) result(r)
    real*8 chi, omegak, r, oc2
    oc2 = omegak*chi**2  
    if(abs(oc2).lt. 0.2d0)then
       r = chi * (1.d0 + oc2*(1.d0/6.d0 + oc2*(1.d0/120.d0+oc2*(1.d0/5040.d0+oc2/362880.d0))))
    elseif(oc2 .gt. 0.d0)then
       oc2 = sqrt(oc2)
       r = chi*(sinh(oc2)/oc2)
    else
       oc2 = sqrt(-oc2)
       r = chi*(sin(oc2)/oc2)
    endif
  end function curved_rofchi


  !!========================= cosmology applications for "PAge" model ====================
  

  function page_Hoft(t0, eta, t) result(H)
    real*8:: t0, eta, t, H
    H = (2.d0/3.d0)/t + 1.d0 - (2.d0/3.d0)/t0 * (1.d0 - eta*(t/t0-1.d0))
  end function page_Hoft

  function page_aoft(t0, eta, t) result(a)
    real*8::t0, eta, a, t
    a =(t/t0)**(2.d0/3.d0) * exp(eta/3.d0 * ((t/t0)**2 - 1)  + (t0-(2.d0/3.d0)*(1.d0+eta))*(t/t0-1.d0) )
  end function page_aoft

  function page_tofa(t0, eta, a) result(t)
    real*8::t0, eta, a, last, ap, t, s, m, u, eps, c0, c1, c2,c3, c4, c5, dap
    integer::i
    c0 = t0*1.5d0    
    u = c0 - 1.d0 - eta
    ap = a**1.5d0
    if(ap .gt. 0.7)then  !!expand around t/t_0 - 1
       dap = 1.d0-ap
       c1 = (eta**2/2 + eta*u + 3*eta/2 + u**2/2 + u)
       c2 = (-eta**3/6 - eta**2*u/2 - eta**2 - eta*u**2/2 - 3*eta*u/2 - eta/2 - u**3/6 - u**2/2)
       c3 = (eta**4/24 + eta**3*u/6 + 5*eta**3/12 + eta**2*u**2/4 + eta**2*u + 5*eta**2/8 + eta*u**3/6 + 3*eta*u**2/4 + eta*u/2 + u**4/24 + u**3/6)
       c4 = (-eta**5/120 - eta**4*u/24 - eta**4/8 - eta**3*u**2/12 - 5*eta**3*u/12 - 3*eta**3/8 - eta**2*u**3/12 - eta**2*u**2/2 - 5*eta**2*u/8 - eta**2/8 - eta*u**4/24 - eta*u**3/4 - eta*u**2/4 - u**5/120 - u**4/24)
       c5 = (eta**6/720 + eta**5*u/120 + 7*eta**5/240 + eta**4*u**2/48 + eta**4*u/8 + 7*eta**4/48 + eta**3*u**3/36 + 5*eta**3*u**2/24 + 3*eta**3*u/8 + 7*eta**3/48 + eta**2*u**4/48 + eta**2*u**3/6 + 5*eta**2*u**2/16 + eta**2*u/8 + eta*u**5/120 + eta*u**4/16 + eta*u**3/12 + u**6/720 + u**5/120)
       eps =  dap/c0
       last = 1.d0 - eps
       eps = dap/(c0 - eps*(c1+eps*(c2 + eps*(c3+eps*(c4+eps*c5)))))
       eps = dap/(c0 - eps*(c1+eps*(c2 + eps*(c3+eps*(c4+eps*c5)))))
       t = 1.d0-eps
    else
       s = eta/2.d0
       m = u + s
       t = ap / exp( (s*ap  + m)*(ap-1.d0))
       t = ap / exp( (s*t  + m)*(t-1.d0))
    endif
    t = t*t0
    do i=1, 30
       ap = a - page_aoft(t0,eta,t)
       if(abs(ap) .lt. 1.d-6) return
       t = max(min(t*1.05,min(t + ap/page_Hoft(t0,eta,t)/a, t0)), t*0.95)
    enddo
    write(*,*) "Error:  page_tofa does not converge", t0, eta, t, a, page_aoft(t0, eta, t)
    stop
  end function page_tofa
  
  function page_Hofa(t0, eta, a) result(H)
    real*8::a, H, t0, eta, t
    t = page_tofa(t0, eta, a)
    H = page_Hoft(t0, eta, t)
  end function page_Hofa

  function page_qoft(t0, eta, t) result(q) !! deceleration parameter q = - a \ddot a/ (\dot a)^2 =  - \dot H/H^2 - 1
    real*8::t0,eta, t, q
    q = (2.d0/3.d0)*(1.d0/t**2 -eta/t0**2)/page_Hoft(t0,eta,t)**2 - 1.d0
  end function page_qoft

  function page_qofa(t0, eta, a) result(q) !! deceleration parameter q = - a \ddot a/ (\dot a)^2 =  - \dot H/H^2 - 1
    real*8::t0,eta, a, t, q
    t = page_tofa(t0, eta, a)
    q = page_qoft(t0, eta, t)
  end function page_qofa

  function page_effective_wofa(t0, eta, a) result(w)  !!effective w  for total energy density
    real*8::t0, eta, a, w
    w = (2.d0*page_qofa(t0,eta,a)-1.d0)/3.d0
  end function page_effective_wofa


  function page_effective_ommofa(t0, eta, a) result(omm)  !!effective w for DE assuming Omega_k = 0
    real*8::t0, eta, a, omm
    real*8,parameter::atiny = 1.d-9
    omm = (page_Hofa(t0, eta, atiny)/page_Hofa(t0,eta,a))**2*(atiny/a)**3
  end function page_effective_ommofa

  
  function page_effective_wdeofa(t0, eta, a) result(wde)  !!effective w for DE assuming Omega_k = 0
    real*8::t0, eta, a, wde, oml
    real*8,parameter::atiny = 1.d-8
    oml = 1.d0 - page_effective_ommofa(t0, eta, a)
    wde = (2.d0*page_qofa(t0,eta,a)-1.d0)/(3.d0*oml)
  end function page_effective_wdeofa
  
  
  function page_dchidt(t0, eta, t) result(intchi)
    real*8::t0, eta, t, intchi
    intchi = 1.d0/page_aoft(t0, eta, t)
  end function page_dchidt
  
  function page_chioft(t0, eta, t) result(chi)
    real*8,parameter::step = 2.d-3
    real*8,parameter::pivot = step*10.d0
    real*8::t0, eta, t, chi
    if(t/t0 .gt. pivot)then
       chi = twop_integrate(page_dchidt, t0, eta, t, t0, step)
    else
       chi = twop_integrate(page_dchidt, t0, eta, pivot, t0, step) &
            + twop_integrate(page_dchidt, t0, eta, t, pivot, step/10.d0)
    endif
  end function page_chioft

  function page_chiofa(t0, eta, a) result(chi)
    real*8::t0, eta, t, a, chi
    t = page_tofa(t0, eta, a)
    chi = page_chioft(t0, eta, t)
  end function page_chiofa

  function page_dlofa(t0, eta, omegak, a) result(dl)
    real*8::t0,eta, omegak,a, dl
    dl = curved_rofchi(page_chiofa(t0,eta,a), omegak)/a       
  end function page_dlofa

  function page_distance_moduli(t0, eta, omegak, h, zhel, zcmb) result(mu)
    real*8::t0, eta ,zhel, zcmb, mu, omegak, h
    mu = 5.d0*log10((1.0+zhel)/(1.0+zcmb)*page_dlofa(t0, eta, omegak, 1.d0/(1.d0+zcmb))*(3000.d0/h))
  end function page_distance_moduli



  !!============================ lcdm applications ==============================
  function lcdm_dtda(omegam, omegak, a) result(dtda)
    real*8::omegam, omegak, a, dtda
    dtda = sqrt(a/(omegam+(1.d0-omegam-omegak)*a**3+omegak*a))
  end function lcdm_dtda

  
  function lcdm_tofa(omegam, omegak, a) result(t)
    real*8,parameter::step = 3.d-3
    real*8,parameter::pivot =step*5.d0
    real*8::omegam, omegak, a, t
    if(a .lt. pivot)then
       t = twop_integrate(lcdm_dtda, omegam, omegak, 0.d0, a, step/10.d0)
    else
       t = twop_integrate(lcdm_dtda, omegam, omegak, 0.d0, pivot, step/10.d0)  &
            + twop_integrate(lcdm_dtda, omegam, omegak, pivot, a, step)
    endif
  end function lcdm_tofa

  function lcdm_Hofa(omegam, omegak, a) result(H)
    real*8::omegam, omegak, a, H
    H = sqrt((omegam/a + omegak)/a**2+1.d0-omegam-omegak)
  end function lcdm_Hofa

  function lcdm_dchida(omegam, omegak, a) result(intchi)
    real*8::omegam, omegak, a, intchi
    intchi = 1.d0/sqrt((omegam+ (omegak + (1.d0-omegam-omegak)*a**2)*a)*a)
  end function lcdm_dchida

  function lcdm_chiofa(omegam, omegak, a) result(chi)
    real*8,parameter::step = 2.d-3
    real*8,parameter::pivot =step*5.d0
    real*8:: omegam, omegak, a, chi
    if(a .gt. pivot)then
       chi = twop_integrate(lcdm_dchida, omegam, omegak, a, 1.d0, step)
    else
       chi = twop_integrate(lcdm_dchida, omegam, omegak, pivot, 1.d0, step) &
            + twop_integrate(lcdm_dchida, omegam, omegak, a, pivot, step/10.d0)
    endif
  end function lcdm_chiofa

  function lcdm_dlofa(omegam, omegak, a) result(dl)
    real*8::omegam, omegak,a, dl
    dl = curved_rofchi(lcdm_chiofa(omegam, omegak,a), omegak)/a       
  end function lcdm_dlofa

  function lcdm_distance_moduli(omegam, omegak, h, zhel, zcmb) result(mu)
    real*8::omegam ,zhel, zcmb, mu, omegak, h
    mu = 5.d0*log10((1.0+zhel)/(1.0+zcmb)*lcdm_dlofa(omegam, omegak, 1.d0/(1.d0+zcmb))*(3000.d0/h))
  end function lcdm_distance_moduli

  !!=============================wcdm applications====================================

  function wcdm_dtda(omegam, omegak, w, a) result(dtda)
    real*8::omegam, omegak, w, a, dtda
    dtda = sqrt(a/(omegam+(1.d0-omegam-omegak)*a**(-3.d0*w)+omegak*a))
  end function wcdm_dtda

  
  function wcdm_tofa(omegam, omegak, w, a) result(t)
    real*8,parameter::step = 3.d-3
    real*8,parameter::pivot =step*5.d0    
    real*8::omegam, omegak, a, t, w
    if(a .lt. pivot)then
       t = threep_integrate(wcdm_dtda, omegam, omegak, w, 0.d0, a, step/10.d0)
    else
       t = threep_integrate(wcdm_dtda, omegam, omegak, w, 0.d0, pivot, step/10.d0) &
            + threep_integrate(wcdm_dtda, omegam, omegak, w, pivot, a, step)
    endif
  end function wcdm_tofa

  function wcdm_Hofa(omegam, omegak, w, a) result(H)
    real*8::omegam, omegak, a, H, w
    H = sqrt((omegam/a + omegak)/a**2+(1.d0-omegam-omegak)*a**(-3.d0*(1.d0+w)))
  end function wcdm_Hofa

  function wcdm_dchida(omegam, omegak, w, a) result(intchi)
    real*8::omegam, omegak, a, intchi, w
    intchi = 1.d0/sqrt((omegam+ (omegak + (1.d0-omegam-omegak)*a**(-1.d0-3.d0*w))*a)*a)
  end function wcdm_dchida

  function wcdm_chiofa(omegam, omegak, w, a) result(chi)
    real*8,parameter::step = 3.d-3
    real*8,parameter::pivot = step * 5.d0    
    real*8:: omegam, omegak, a, chi, w
    if(a.gt. pivot)then
       chi = threep_integrate(wcdm_dchida, omegam, omegak, w, a, 1.d0, step)
    else
       chi = threep_integrate(wcdm_dchida, omegam, omegak, w, pivot, 1.d0, step) &
            + threep_integrate(wcdm_dchida, omegam, omegak, w, a, pivot, step/10.d0)       
    endif
  end function wcdm_chiofa

  function wcdm_dlofa(omegam, omegak, w, a) result(dl)
    real*8::omegam, omegak,a, dl, w
    dl = curved_rofchi(wcdm_chiofa(omegam, omegak, w, a ), omegak)/a       
  end function wcdm_dlofa

  function wcdm_distance_moduli(omegam, omegak, w, h, zhel, zcmb) result(mu)
    real*8::omegam ,zhel, zcmb, mu, omegak, w, h
    mu = 5.d0*log10((1.0+zhel)/(1.0+zcmb)*wcdm_dlofa(omegam, omegak, w, 1.d0/(1.d0+zcmb))*(3000.d0/h))
  end function wcdm_distance_moduli


  !!=============================w0wa applications====================================
  function w0wa_dtda(omegam, omegak, w0, wa, a) result(dtda)
    real*8::omegam, omegak, w0, wa, a, dtda
    dtda = sqrt(a/(omegam+(1.d0-omegam-omegak)*exp(-3.d0*((w0+wa)*log(max(a, 1.d-20))+wa*(1.d0-a))) + omegak*a))
  end function w0wa_dtda

  
  function w0wa_tofa(omegam, omegak, w0, wa, a) result(t)
    real*8,parameter::step = 3.d-3
    real*8,parameter::pivot =step*5.d0    
    real*8::omegam, omegak, a, t, w0, wa
    if(a .lt. pivot)then
       t = fourp_integrate(w0wa_dtda, omegam, omegak, w0, wa, 0.d0, a, step/10.d0)
    else
       t = fourp_integrate(w0wa_dtda, omegam, omegak, w0, wa, 0.d0, pivot, step/10.d0) &
            + fourp_integrate(w0wa_dtda, omegam, omegak, w0, wa, pivot, a, step)
    endif
  end function w0wa_tofa

  function w0wa_aoft(omegam, omegak, w0, wa, t) result(a)
    real*8,parameter::small_step = 2.d-4
    real*8,parameter::large_step = 2.d-3
    real*8,parameter::asmall = small_step * 10.d0    
    real*8::omegam, omegak, a, t, w0, wa, tsum, adot1, adot2, adotm, tlarge, da
    tlarge = t  * 0.98
    a = 1.d-8
    tsum = 0.d0
    adot1 = w0wa_dtda(omegam, omegak, w0, wa, a)    
    do while(tsum .lt. t )
       if(t .gt. tlarge  .or.  a .lt. asmall)then
          adot2 = adot1
          adotm = w0wa_dtda(omegam, omegak, w0, wa, a+small_step/2.d0)          
          adot1 = w0wa_dtda(omegam, omegak, w0, wa, a+small_step)
          tsum = tsum + (adot1 + 4.d0*adotm + adot2)*(small_step/6.d0)
          a = a + small_step
       else
          adot2 = adot1
          adotm = w0wa_dtda(omegam, omegak, w0, wa, a+large_step/2.d0)                    
          adot1 = w0wa_dtda(omegam, omegak, w0, wa, a+large_step)
          tsum = tsum + (adot1 + 4.d0*adotm + adot2)*(large_step/2.d0)
          a = a + large_step
       endif
    enddo
    da = (t-tsum)/w0wa_dtda(omegam, omegak, w0, wa, a)    
    a = a + (da + (t-tsum)/w0wa_dtda(omegam, omegak, w0, wa, a+da))/2.d0    
  end function w0wa_aoft
  
  function w0wa_Hofa(omegam, omegak, w0, wa, a) result(H)
    real*8::omegam, omegak, a, H, w0, wa
    H = sqrt((omegam/a + omegak)/a**2+(1.d0-omegam-omegak)*exp(-3.d0*((1.d0+w0+wa)*log(max(a, 1.d-20))+wa*(1.d0-a))))
  end function w0wa_Hofa

  function w0wa_dchida(omegam, omegak, w0, wa, a) result(intchi)
    real*8::omegam, omegak, a, intchi, w0, wa
    intchi = 1.d0/sqrt((omegam+ omegak*a)*a + (1.d0-omegam-omegak) *  exp((1.d0-3.d0*(w0+wa))*log(max(a, 1.d-20))-3.d0*wa*(1.d0-a)) )
  end function w0wa_dchida

  function w0wa_chiofa(omegam, omegak, w0, wa, a) result(chi)
    real*8,parameter::step = 3.d-3
    real*8,parameter::pivot = step * 5.d0    
    real*8:: omegam, omegak, a, chi, w0, wa
    if(a.gt. pivot)then
       chi = fourp_integrate(w0wa_dchida, omegam, omegak, w0, wa, a, 1.d0, step)
    else
       chi = fourp_integrate(w0wa_dchida, omegam, omegak, w0, wa, pivot, 1.d0, step) &
            + fourp_integrate(w0wa_dchida, omegam, omegak, w0, wa, a, pivot, step/10.d0)       
    endif
  end function w0wa_chiofa

  function w0wa_dlofa(omegam, omegak, w0, wa, a) result(dl)
    real*8::omegam, omegak,a, dl, w0, wa
    dl = curved_rofchi(w0wa_chiofa(omegam, omegak, w0, wa, a ), omegak)/a       
  end function w0wa_dlofa

  function w0wa_distance_moduli(omegam, omegak, w0, wa, h, zhel, zcmb) result(mu)
    real*8::omegam ,zhel, zcmb, mu, omegak, w0, wa, h
    mu = 5.d0*log10((1.0+zhel)/(1.0+zcmb)*w0wa_dlofa(omegam, omegak, w0, wa, 1.d0/(1.d0+zcmb))*(3000.d0/h))
  end function w0wa_distance_moduli


  !!============================ rhct applications ==============================
  function rhct_dtda(a) result(dtda)
    real*8::a, dtda
    dtda = 1.d0
  end function rhct_dtda

  
  function rhct_tofa(a) result(t)
    real*8::t, a
    t = a
  end function rhct_tofa

  function rhct_Hofa(a) result(H)
    real*8::a, H
    H = 1.d0/a
  end function rhct_Hofa

  function rhct_dchida(a) result(intchi)
    real*8::a, intchi
    intchi = 1.d0/a
  end function rhct_dchida

  function rhct_chiofa(a) result(chi)
    real*8:: a, chi
    chi = -log(a)
  end function rhct_chiofa

  function rhct_dlofa(omegak, a) result(dl)
    real*8::a, dl, omegak
    dl = curved_rofchi(rhct_chiofa(a), omegak)/a       
  end function rhct_dlofa

  function rhct_distance_moduli(omegak, h, zhel, zcmb) result(mu)
    real*8::zhel, zcmb, mu, omegak, h
    mu = 5.d0*log10((1.0+zhel)/(1.0+zcmb)*rhct_dlofa(omegak, 1.d0/(1.d0+zcmb))*(3000.d0/h))
  end function rhct_distance_moduli
  


  !!=============================DGP applications====================================
  
  function DGP_dtda(omegam, a) result(dtda)
    real*8::omegam, a, dtda, sqrta
    sqrta = sqrt(a)
    dtda = 2.d0*sqrta/((1.d0-omegam)*a*sqrta+ sqrt((1.d0-omegam)**2*a**3 + 4.d0*omegam))
  end function DGP_dtda

  
  function DGP_tofa(omegam, a) result(t)
    real*8,parameter::step = 3.d-3
    real*8,parameter::pivot =step*5.d0
    real*8::omegam, a, t
    if(a .lt. pivot)then
       t = onep_integrate(DGP_dtda, omegam, 0.d0, a, step/10.d0)
    else
       t = onep_integrate(DGP_dtda, omegam, 0.d0, pivot, step/10.d0)  &
            + onep_integrate(DGP_dtda, omegam, pivot, a, step)
    endif
  end function DGP_tofa

  function DGP_Hofa(omegam, a) result(H)
    real*8::omegam, a, H
    H = ((1.d0-omegam)+ sqrt((1.d0-omegam)**2 + 4.d0*omegam/a**3))/2.d0
  end function DGP_Hofa

  function DGP_dchida(omegam, a) result(intchi)
    real*8::omegam, a, intchi
    intchi = DGP_dtda(omegam, a)/a
  end function DGP_dchida

  function DGP_chiofa(omegam, a) result(chi)
    real*8,parameter::step = 2.d-3
    real*8,parameter::pivot =step*5.d0
    real*8:: omegam, a, chi
    if(a .gt. pivot)then
       chi = onep_integrate(DGP_dchida, omegam, a, 1.d0, step)
    else
       chi = onep_integrate(DGP_dchida, omegam, pivot, 1.d0, step) &
            + onep_integrate(DGP_dchida, omegam, a, pivot, step/10.d0)
    endif
  end function DGP_chiofa

  function DGP_dlofa(omegam, omegak, a) result(dl)
    real*8::omegam, omegak, a, dl
    dl = curved_rofchi(DGP_chiofa(omegam, a), omegak)/a       
  end function DGP_dlofa

  function DGP_distance_moduli(omegam, omegak, h, zhel, zcmb) result(mu)
    real*8::omegam ,zhel, zcmb, mu, omegak, h
    mu = 5.d0*log10((1.0+zhel)/(1.0+zcmb)*DGP_dlofa(omegam, omegak, 1.d0/(1.d0+zcmb))*(3000.d0/h))
  end function DGP_distance_moduli


  !!=============================Generalized Chapylygin Gas applications====================================
  function GCG_Hofa(omegab, As, zeta, a) result(H)
    real*8::omegab, As, zeta, a, H
    H = sqrt(omegab/a**3+(1.d0-omegab)*(As + (1.d0-As)/a**(3.d0*(1.d0+zeta)))**(1.d0/(1.d0+zeta)))
  end function GCG_Hofa

  function GCG_dtda(omegab, As, zeta, a) result(dtda)
    real*8::omegab, As, zeta, a, dtda
    dtda = 1.d0/GCG_Hofa(omegab, As, zeta, a)/a
  end function GCG_dtda

  
  function GCG_tofa(omegab, As, zeta, a) result(t)
    real*8,parameter::step = 3.d-3
    real*8,parameter::pivot =step*5.d0    
    real*8::omegab, As, zeta, a, t
    if(a .lt. pivot)then
       t = threep_integrate(GCG_dtda, omegab, As, zeta, 1.d-10, a, step/10.d0)
    else
       t = threep_integrate(GCG_dtda, omegab, As, zeta, 1.d-10, pivot, step/10.d0) &
            + threep_integrate(GCG_dtda, omegab, As, zeta, pivot, a, step)
    endif
  end function GCG_tofa


  function GCG_dchida(omegab, As, zeta, a) result(intchi)
    real*8::omegab, As, a, intchi, zeta
    intchi = 1.d0/GCG_Hofa(omegab, As, zeta, a)/a**2
  end function GCG_dchida

  function GCG_chiofa(omegab, As, zeta, a) result(chi)
    real*8,parameter::step = 3.d-3
    real*8,parameter::pivot = step * 5.d0    
    real*8:: omegab, As, a, chi, zeta
    if(a.gt. pivot)then
       chi = threep_integrate(GCG_dchida, omegab, As, zeta, a, 1.d0, step)
    else
       chi = threep_integrate(GCG_dchida, omegab, As, zeta, pivot, 1.d0, step) &
            + threep_integrate(GCG_dchida, omegab, As, zeta, a, pivot, step/10.d0)       
    endif
  end function GCG_chiofa

  function GCG_dlofa(omegab, As, zeta, omegak, a) result(dl)
    real*8::omegab, omegak,a, dl, As, zeta
    dl = curved_rofchi(GCG_chiofa(omegab, As, zeta, a), omegak)/a       
  end function GCG_dlofa

  function GCG_distance_moduli(omegab, As, zeta, omegak, h, zhel, zcmb) result(mu)
    real*8::omegab, As ,zhel, zcmb, mu, omegak, zeta, h
    mu = 5.d0*log10((1.0+zhel)/(1.0+zcmb)*GCG_dlofa(omegab, As, zeta, omegak, 1.d0/(1.d0+zcmb))*(3000.d0/h))
  end function GCG_distance_moduli

  
  
End Module coop_page_mod



