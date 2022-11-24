module coop_expint_mod

  use coop_wrapper_utils
  implicit none  
#include "constants.h"
contains

  function coop_expint( f, a, b, accuracy) !!must have b>a>0 and f(x)>0 for all x
    external f
    COOP_REAL::f, a, b, fa, fb, fm, lna, lnb, lnf2, lnf1, dt, t, coop_expint, tmp, acc, dlnf, cut
    COOP_REAL, optional::accuracy
    COOP_INT::n, i
    if(present(accuracy))then
       acc = accuracy
    else
       acc = 1.d-5
    endif
    lna = dlog(a)
    lnb = dlog(b)    
    n = ceiling(abs(lnb-lna)/acc) + 9
    dt = (lnb - lna)/n
    coop_expint = 0.d0
    tmp = 0.d0
    fa = f(a)
    fb = f(b)
    fm = f((a+b)/2.d0)
    lnf2 = dlog(fa) + lna
    t = lna
    cut = 0.d0
    do i = 1, n-1
       lnf1 = lnf2
       t = t + dt
       lnf2 = dlog(f(dexp(t)))+t
       dlnf = lnf2 - lnf1
       if(abs(dlnf) .lt. 3.d-5)then
          tmp = tmp + dexp(lnf1)*(1.d0+dlnf*(0.5d0+dlnf/6.d0))
       else
          tmp = tmp + dexp(lnf1)*(dexp(dlnf)-1.d0)/dlnf
       endif
       if(tmp > cut)then
          coop_expint = coop_expint + tmp
          tmp = 0.d0
          cut  = coop_expint * 1.d-5
       endif
    enddo
    coop_expint = (coop_expint + tmp)*dt
  end function coop_expint


  
  !!integrate \int_a^b e^{zt} f(t) dt 
  function coop_period_int(f, a, b, z, a_cut, b_cut, accuracy) result(pint)
    external f
    COOP_INT, parameter::min_level = 3, max_level = 10    
    COOP_REAL::a, b, f
    COOP_REAL::x, y, absy
    COOP_COMPLEX :: z, expp, expm, pint, intres(min_level-1:max_level)
    COOP_REAL:: halfp, period
    logical,optional::a_cut, b_cut
    logical::acut, bcut
    COOP_REAL::tstart, tend, acc
    COOP_INT::level
    COOP_REAL,optional::accuracy
    if(present(accuracy))then
       acc = accuracy
    else
       acc = 1.d-5
    endif
    if(present(a_cut))then
       acut = a_cut
    else
       acut = .false.
    endif
    if(present(b_cut))then
       bcut = b_cut
    else
       bcut = .false.
    endif    
    x = real(z)
    y = aimag(z)
    absy = abs(y)
    if(absy .gt. coop_pi/abs(b-a))then
       halfp = coop_pi/absy
       tstart = halfp*nint(a/halfp)
       if(tstart .gt. a) tstart = tstart - halfp
       tend = tstart + (2.d0*halfp) * ceiling((b - tstart)/(2.d0*halfp)) 
    else
       halfp = 0.5d0
       tstart = a
       tend = b
    endif
    tend = tend + halfp/2**(max_level +1)
    period = 2.d0*halfp
    expp = exp(halfp*z)
    expm = exp(-halfp*z)
    level = min_level
    intres(level) = int_grid(tstart, tend, halfp/2**level) - (g_eval(tstart) + g_eval(tend))*halfp/2**(level+3)
    do 
       level = level + 1       
       intres(level) = (int_grid(tstart + halfp/2**level, tend, halfp/2**(level-1)) + intres(level-1))/2.d0
       if(abs(intres(level) - intres(level-1)) .lt. acc*abs(intres(level)) .or. level .eq. max_level)then
          pint = intres(level)
          return
       endif
    enddo
  contains
    
    function int_grid(ts, te, dt) result(intall)
      COOP_COMPLEX::intall, tmp
      COOP_REAL::fcut, ts, te, dt, t
      tmp = 0.d0
      intall = 0.d0
      t=ts
      if(acut)then
         if(bcut)then
            do while(t .le. te)
               tmp = tmp  + g_leftright_cuts(t)
               t = t + dt
               if(abs(tmp) > fcut)then
                  intall = intall + tmp
                  tmp = 0.d0
                  fcut = abs(intall)*1.d-4
               endif
            enddo
         else
            do while(t .le. te)
               tmp = tmp  + g_left_cut(t)
               t = t + dt
               if(abs(tmp) > fcut)then
                  intall = intall + tmp
                  tmp = 0.d0
                  fcut = abs(intall)*1.d-4
               endif
            enddo
         endif
      else
         if(bcut)then
            do while(t .le. te)
               tmp = tmp  + g_right_cut(t)
               t = t + dt
               if(abs(tmp) > fcut)then
                  intall = intall + tmp
                  tmp = 0.d0
                  fcut = abs(intall)*1.d-4
               endif
            enddo
         else
            do while(t .le. te)
               tmp = tmp  + g_no_cuts(t)
               t = t + dt
               if(abs(tmp) > fcut)then
                  intall = intall + tmp
                  tmp = 0.d0
                  fcut = abs(intall)*1.d-4
               endif
            enddo
         endif
      endif
      intall = (intall+tmp)*dt/4.d0
    end function int_grid


    function g_no_cuts(t) result(g)
      COOP_COMPLEX::g, gl, gr
      COOP_REAL::t
      g = (2*f(t) + f(t - halfp)*expm + f(t + halfp)*expp) * exp(z*t)
    end function g_no_cuts
    
    function g_leftright_cuts(t) result(g)
      COOP_COMPLEX::g, gl, gr
      COOP_REAL::t
      if(t .lt. a)then
         gl = 0.d0
      else
         gl = f(t - halfp)*expm
      endif
      if(t .gt. b)then
         gr = 0.d0
      else
         gr = f(t + halfp)*expp
      endif
      g = (2*f(t) + gl + gr) * exp(z*t)
    end function g_leftright_cuts
    
    function g_left_cut(t) result(g)
      COOP_COMPLEX::g, gl, gr
      COOP_REAL::t
      if(t .lt. a)then
         gl = 0.d0
      else
         gl = f(t - halfp)*expm
      endif
      gr = f(t + halfp)*expp
      g = (2*f(t) + gl + gr) * exp(z*t)
    end function g_left_cut

    function g_right_cut(t) result(g)
      COOP_COMPLEX::g, gl, gr
      COOP_REAL::t
      gl = f(t - halfp)*expm
      if(t .gt. b)then
         gr = 0.d0
      else
         gr = f(t + halfp)*expp
      endif
      g = (2*f(t) + gl + gr) * exp(z*t)
    end function g_right_cut

    function g_eval(t) result(g)
      COOP_COMPLEX::g
      COOP_REAL::t
      if(acut)then
         if(bcut)then
            g = g_leftright_cuts(t)
         else
            g = g_left_cut(t)
         endif
      else
         if(bcut)then
            g = g_right_cut(t)
         else
            g = g_no_cuts(t)
         endif         
      endif
    end function g_eval
    
  end function coop_period_int
  
end module coop_expint_mod
