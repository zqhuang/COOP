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
  !!model: H = 2/(3t) + 1-2/(3t0) +  2 eta/(3t0)  (t/t0-1) 
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
    real*8::t0, eta, a, last, ap, t, s, m
    integer::i
    ap = a**1.5d0
    s = eta/2.d0
    m = t0*1.5d0-1.d0-eta + s
    t = ap / exp( (s*ap  + m)*(ap-1.d0))
    t = ap / exp( (s*t  + m)*(t-1.d0))
    t = ap / exp( (s*t  + m)*(t-1.d0))
    last = t    
    t = ap / exp( (s*t  + m)*(t-1.d0))
    do i=1,25
       if(abs(last - t) .lt. 1.d-7) exit
       last = t              
       t = ap / exp( (s*t  + m)*(t-1.d0))
    enddo
    t = t*t0
  end function page_tofa
  
  function page_Hofa(t0, eta, a) result(H)
    real*8::a, H, t0, eta, t
    t = page_tofa(t0, eta, a)
    H = page_Hoft(t0, eta, t)
  end function page_Hofa

  function page_qoft(t0, eta, t) result(q) !! deceleration parameter q = - a \ddot a/ (\dot a)^2 =  - \dot H/H^2 - 1
    real*8::t0,eta, t, q
    q = (2.d0/3.d0)*(1.d0/t**2 -eta/t0**2) - 1.d0
  end function page_qoft

  function page_qofa(t0, eta, a) result(q) !! deceleration parameter q = - a \ddot a/ (\dot a)^2 =  - \dot H/H^2 - 1
    real*8::t0,eta, a, t, q
    t = page_tofa(t0, eta, a)
    q = page_qoft(t0, eta, t)
  end function page_qofa

  function page_effective_wofa(t0, eta, a) result(w)  !!effective w assuming \Omega_k = 0
    real*8::t0, eta, a, w, t
    
  end function page_effective_wofa
  
  function page_dchidt(t0, eta, t) result(intchi)
    real*8::t0, eta, t, intchi
    intchi = 1.d0/page_aoft(t0, eta, t)
  end function page_dchidt
  
  function page_chioft(t0, eta, t) result(chi)
    real*8::t0, eta, t, chi
    chi = twop_integrate(page_dchidt, t0, eta, t, t0, 1.d-2)
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
    real*8::omegam, omegak, a, t
    t = twop_integrate(lcdm_dtda, omegam, omegak, 0.d0, a, 0.01d0)
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
    real*8:: omegam, omegak, a, chi
    chi = twop_integrate(lcdm_dchida, omegam, omegak, a, 1.d0, 1.d-2)
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
    real*8::omegam, omegak, a, t, w
    t = threep_integrate(wcdm_dtda, omegam, omegak, w, 0.d0, a, 0.01d0)
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
    real*8:: omegam, omegak, a, chi, w
    chi = threep_integrate(wcdm_dchida, omegam, omegak, w, a, 1.d0, 1.d-2)
  end function wcdm_chiofa

  function wcdm_dlofa(omegam, omegak, w, a) result(dl)
    real*8::omegam, omegak,a, dl, w
    dl = curved_rofchi(wcdm_chiofa(omegam, omegak, w, a ), omegak)/a       
  end function wcdm_dlofa

  function wcdm_distance_moduli(omegam, omegak, w, h, zhel, zcmb) result(mu)
    real*8::omegam ,zhel, zcmb, mu, omegak, w, h
    mu = 5.d0*log10((1.0+zhel)/(1.0+zcmb)*wcdm_dlofa(omegam, omegak, w, 1.d0/(1.d0+zcmb))*(3000.d0/h))
  end function wcdm_distance_moduli


  
End Module coop_page_mod



