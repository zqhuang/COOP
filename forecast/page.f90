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

  function page_int_chi(t0, eta, t) result(intchi)
    real*8::t0, eta, t, intchi
    intchi = 1.d0/page_aoft(t0, eta, t)
  end function page_int_chi
  
  function page_chioft(t0, eta, t) result(chi)
    real*8::t0, eta, t, chi
    chi = twop_integrate(page_int_chi, t0, eta, t, t0, 1.d-2)
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



  function page_snpre_age_distr(t0, eta, t, tau) result(p)
    real*8::t0, eta, tau, t, P
    if(t-tau .lt. 1.d-2)then  !!exclude very high-redshift star formation (avoid numeric overflow in sfr)
       P = 0.d0
       return
    endif
    P = dtd(tau) * sfr(t-tau)
  contains

    function dtd(t)
      real*8::dtd, t
      real*8,parameter::tp=0.0215 !! equivalent to 0.3(h/0.7) Gyr
      integer,parameter::alpha = 20
      dtd = t/tp
      if(dtd .lt. 1.d0)then
         dtd = dtd**alpha/(1.d0+dtd**(alpha-1))
      else
         dtd = dtd/(1.d0+dtd**(1-alpha))
      endif
    end function dtd

    function sfr(t)
      real*8::sfr, t, z
      real*8,parameter::a=-0.997, b=0.241, z0 = 1.243
      z = 1.d0/page_aoft(t0, eta, t) - 1.d0 
      sfr = 1.d0/(10.d0**(a*(z-z0)) + 10.d0**(b*(z-z0)))
    end function sfr
  end function page_snpre_age_distr

  function page_snpre_aveage(t0, eta, t) result(aveage)
    real*8::t0, eta, t, aveage, dtau, sump_odd, sump_even, sumt_odd, sumt_even, p, tau
    integer::n , i
    sump_odd = 0.d0
    sump_even = 0.d0
    sumt_odd = 0.d0
    sumt_even = 0.d0
    n = floor(t/0.01)+1
    dtau = t/n/2.d0
    tau = 0.d0
    do i=1, n-1
       tau = tau + dtau
       p = page_snpre_age_distr(t0, eta, t, tau)
       sump_odd = sump_odd + p
       sumt_odd = sumt_odd + p * tau
       tau = tau + dtau
       p = page_snpre_age_distr(t0, eta, t, tau)
       sump_even = sump_even + p
       sumt_even = sumt_even + p * tau
    enddo
    tau = tau + dtau
    p = page_snpre_age_distr(t0, eta, t, tau)
    sump_odd = sump_odd + p
    sumt_odd = sumt_odd + p * tau
    aveage = (sumt_odd*4.d0 + sumt_even * 2.d0) &
         / (sumt_odd*4.d0 + sumt_even * 2.d0  &
         + page_snpre_age_distr(t0, eta, t, 0.d0))
  end function page_snpre_aveage


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

  function lcdm_int_chi(omegam, omegak, a) result(intchi)
    real*8::omegam, omegak, a, intchi
    intchi = 1.d0/sqrt((omegam+ (omegak + (1.d0-omegam-omegak)*a**2)*a)*a)
  end function lcdm_int_chi

  function lcdm_chiofa(omegam, omegak, a) result(chi)
    real*8:: omegam, omegak, a, chi
    chi = twop_integrate(lcdm_int_chi, omegam, omegak, a, 1.d0, 1.d-2)
  end function lcdm_chiofa

  function lcdm_dlofa(omegam, omegak, a) result(dl)
    real*8::omegam, omegak,a, dl
    dl = curved_rofchi(lcdm_chiofa(omegam, omegak,a), omegak)/a       
  end function lcdm_dlofa

  function lcdm_distance_moduli(omegam, omegak, h, zhel, zcmb) result(mu)
    real*8::omegam ,zhel, zcmb, mu, omegak, h
    mu = 5.d0*log10((1.0+zhel)/(1.0+zcmb)*lcdm_dlofa(omegam, omegak, 1.d0/(1.d0+zcmb))*(3000.d0/h))
  end function lcdm_distance_moduli
  
  
End Module coop_page_mod



