Module coop_jerk_mod
  implicit none
  !!a: scale factor, normalized to 1 today
  !!t: H0 * age_then
  !!z: redshift, 1/a-1
  !!H: H(z)/H_0
  !!dL: H_0 * luminosity distance
  !!dA: H_0 * angular diameter distance
  !!dc: H0 * comoving distance
  !!t0: H0* age_now
  
  !!model: H = lambda (t0/t) + 1-lambda +  eta lambda (t/t0-1)
  !!The advantage of this parametrization is that we more or less know that
  !!      lambda ~ 2/3, t0 ~ 1.4
  !!      dot H = -lambda/t (t0/t) + eta lambda /t0 < 0 => eta < 1
  !!      q0 =  (1-eta) lambda - 1
  
contains


  function jerk_aoft(q0, j0, t) result(a)
    real*8::q0, j0    
    real*8::a, t
    a = 1.d0 - t*(1.d0 + t*(- q0/2.d0 - t * j0/6.d0))
  end function jerk_aoft

  function jerk_tofa(q0, j0, a) result(t)
    real*8::q0, j0    
    real*8::t, a, s, g, m
    s = 1.d0 - a
    g = - q0/2.d0
    m = -j0/6.d0
    t = s/(1.d0 + s*(g + s * (q0**2/2.d0+m)))
    t = s/(1.d0 + t*(g + t * m))
    t = s/(1.d0 + t*(g + t * m))
    t = s/(1.d0 + t*(g + t * m))
    t = s/(1.d0 + t*(g + t * m))        
  end function jerk_tofa

  function jerk_adotoft(q0, j0, t) result(adot)
    real*8::q0, j0
    real*8::t, adot
    adot = 1.d0 + t*(q0+t*j0/2.d0)
  end function jerk_adotoft

  function jerk_Hofa(q0, j0, a) result(H)
    real*8::q0,j0
    real*8::a, H, t
    t = jerk_tofa(q0, j0, a)
    H = jerk_adotoft(q0,j0,t)/a
  end function jerk_Hofa

  function int_jerk_dc(q0, j0, a) result(intdc)
    real*8::q0,j0
    real*8::a, intdc, t
    t = jerk_tofa(q0, j0, a)
    intdc = 1.d0/jerk_adotoft(q0,j0,t)/a
  end function int_jerk_dc
  
  function jerk_dcofa(q0, j0, a) result(dc)
    integer::n, i
    real*8::q0,j0
    real*8::a, dc, t, da, sumodd, sumeven, anow
    n = floor((1.d0-a)/1.d-2)+1
    da = (1.d0-a)/n/2.d0
    sumodd = 0.d0
    sumeven = 0.d0
    anow = a
    do i=1, n-1
       anow = anow + da
       sumodd = sumodd + int_jerk_dc(q0, j0, anow)
       anow = anow + da
       sumeven = sumeven + int_jerk_dc(q0, j0, anow)       
    enddo
    sumodd = sumodd + int_jerk_dc(q0, j0, anow+da)    
    dc = (sumodd*4.d0 + sumeven*2.d0 + int_jerk_dc(q0, j0, a)+int_jerk_dc(q0, j0, 1.d0))*(da/3.d0)
  end function jerk_dcofa

  function jerk_r_of_chi(chi, Omega_k) result(r)
    real*8 chi, Omega_k, r, oc2
    oc2 = Omega_k*chi**2  
    if(abs(oc2).lt. 0.2d0)then
       r = chi * (1.d0 + oc2*(1.d0/6.d0 + oc2*(1.d0/120.d0+oc2*(1.d0/5040.d0+oc2/362880.d0))))
    elseif(oc2 .gt. 0.d0)then
       oc2 = sqrt(oc2)
       r = chi*(sinh(oc2)/oc2)
    else
       oc2 = sqrt(-oc2)
       r = chi*(sin(oc2)/oc2)
    endif
  end function jerk_r_of_chi
  
  function jerk_dlofa(q0, j0, a, omegak) result(dl)
    real*8::q0,j0,a, dl
    real*8,optional::omegak
    if(present(omegak))then
       dl = jerk_r_of_chi(jerk_dcofa(q0,j0,a), omegak)/a       
    else
       dl = jerk_dcofa(q0,j0,a)/a
    endif
  end function jerk_dlofa


  function jerk_dAofa(q0, j0, a, omegak) result(dA)
    real*8::q0,j0,a, dA
    real*8,optional::omegak
    if(present(omegak))then
       dA = jerk_r_of_chi(jerk_dcofa(q0,j0,a), omegak)*a
    else
       dA = jerk_dcofa(q0,j0,a)*a
    endif
  end function jerk_dAofa

  function jerk_pregenitor_age_distr(q0, j0, tau, t) result(P)
    real*8::q0, j0, tau, t, P
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
      z = 1.d0/jerk_aoft(q0, j0, t) - 1.d0 
      sfr = 1.d0/(10.d0**(a*(z-z0)) + 10.d0**(b*(z-z0)))
    end function sfr
    
  end function jerk_pregenitor_age_distr

  function jerk_pregenitor_age_ave(q0, j0, t) result(aveage)
    real*8::q0, j0, t, aveage, dtau, sump_odd, sump_even, sumt_odd, sumt_even, p, tau
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
       p = jerk_pregenitor_age_distr(q0, j0, tau, t)
       sump_odd = sump_odd + p
       sumt_odd = sumt_odd + p * tau
       tau = tau + dtau
       p = jerk_pregenitor_age_distr(q0, j0, tau, t)
       sump_even = sump_even + p
       sumt_even = sumt_even + p * tau
    enddo
    tau = tau + dtau
    p = jerk_pregenitor_age_distr(q0, j0, tau, t)
    sump_odd = sump_odd + p
    sumt_odd = sumt_odd + p * tau
    aveage = (sumt_odd*4.d0 + sumt_even * 2.d0) &
         / (sumt_odd*4.d0 + sumt_even * 2.d0  &
         + jerk_pregenitor_age_distr(q0, j0, 0.d0, t)) 
  end function jerk_pregenitor_age_ave

  function jerk_distance_moduli(q0, j0, omegak, h, zhel, zcmb) result(mu)
    real*8::q0,j0,a, zhel, zcmb, mu, omegak, h
    mu = 5.d0*log10((1.0+zhel)/(1.0+zcmb)*jerk_dlofa(q0,j0, 1.d0/(1.d0+zcmb), omegak)*(3000.d0/h))
  end function jerk_distance_moduli
  
End Module coop_jerk_mod
