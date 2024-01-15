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
    COOP_INT, parameter::min_level = 3, max_level = 12    
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
       if( abs(a+b) .lt. 1.d-5)then
          tend = -tstart
       else
          tend = tstart + (2.d0*halfp) * ceiling((b - tstart)/(2.d0*halfp))
       endif
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


    !!integrate \int_{-a}^a e^{zt} f(t) dt 
  function coop_periodsym_int(f, a, z, accuracy) result(pint)
    external f
    COOP_INT, parameter::min_level = 3, max_level = 11    
    COOP_REAL::a,  f
    COOP_REAL::x, y, absy
    COOP_COMPLEX :: z, expp, expm, pint, intres(min_level-1:max_level)
    COOP_REAL:: p1, p2, p3, ep1, em1, ep2, em2, ep3, em3, acc
    COOP_INT::level, np
    COOP_REAL,optional::accuracy
    COOP_INT, parameter::sumw = 16  !! 1+2+1 = 4 (2nd order algorithm) or 1+4+6+4+1 = 16 (4th order algorithm) 1 + 6 + 15 + 20 + 15 +  6 + 1 = 64
    if(present(accuracy))then
       acc = accuracy
    else
       acc = 1.d-5
    endif
    x = real(z)
    y = aimag(z)
    absy = abs(y)
    if(absy*a .lt. 2.d0)then
       write(*,*) "periodsym_int is designed for periodic function with the imaginary part of z not too small"
       stop
    endif
    p1 = coop_pi/absy
    np = nint(a/p1)
    ep1 = exp(p1*x)
    em1 = exp(-p1*x)
    if(sumw .gt. 4)then
       p2 = p1*2.d0
       ep2 = exp(p2*x)
       em2 = exp(-p2*x)
    endif
    if(sumw .gt. 16)then
       p3 = p1*3.d0
       ep3 = exp(p3*x)
       em3 = exp(-p3*x)
    endif
    level = min_level
    intres(level) = int_grid(p1/2**level,  np * 2**level) 
    do 
       level = level + 1       
       intres(level) = (int_grid_odd(p1/2**(level-1), np * 2**level) + intres(level-1))/2.d0
       if(abs(intres(level) - intres(level-1)) .lt. acc*abs(intres(level)) .or. level .eq. max_level)then
          pint = intres(level)
          return
       endif
    enddo
  contains

    function g_eval(t) result(g)
      COOP_COMPLEX::g, gl, gr
      COOP_REAL::t
      select case(sumw)
      case(4)
         g = (2*f(t) - (f(t + p1)*ep1 + f(t - p1)*em1)) * exp(z*t)
      case(16)
         g = ( 6*f(t) - 4*(f(t + p1)*ep1 + f(t - p1)*em1) + (f(t+p2)*ep2 + f(t-p2)*em2) ) * exp(z*t)
      case(64)
         g = ( 20*f(t) - 15*(f(t + p1)*ep1 + f(t - p1)*em1) + 6* (f(t+p2)*ep2 + f(t-p2)*em2) - (f(t+p3) * ep3 + f(t-p3) * em3) ) * exp(z*t)
      end select
    end function g_eval
    
    
    function int_grid(dt, nsteps) result(intall)
      COOP_INT,parameter::nthreads = 64
      COOP_COMPLEX::intall, tmp(nthreads), int_thread(nthreads)
      COOP_REAL::dt, fcut(nthreads)
      COOP_INT::nsteps, i, ithread
      tmp = 0.d0
      int_thread = 0.d0
      fcut = 0.d0
      !$omp parallel do
      do ithread = 1, nthreads
         do i=ithread, nsteps, nthreads
            tmp(ithread) = tmp(ithread) + (g_eval(i*dt) + g_eval(-i*dt))
            if( abs(real(tmp(ithread))) + abs(aimag(tmp(ithread))) .gt. fcut(ithread))then
               int_thread(ithread) = int_thread(ithread)  + tmp(ithread)
               tmp(ithread) = 0.d0
               fcut(ithread) = (abs(real(int_thread(ithread))) + abs(aimag(int_thread(ithread))))*1.d-7
            endif
         enddo
      enddo
      !$omp end parallel do
      intall = (g_eval(0.d0) + sum(tmp) + sum(int_thread)) * (dt/sumw)
    end function int_grid


    function int_grid_odd(dt, nsteps) result(intall)
      COOP_INT,parameter::nthreads = 64
      COOP_COMPLEX::intall, tmp(nthreads), int_thread(nthreads)
      COOP_REAL::dt, fcut(nthreads)
      COOP_INT::nsteps, i, ithread
      tmp = 0.d0
      int_thread = 0.d0
      fcut = 0.d0
      !$omp parallel do
      do ithread = 1, nthreads
         do i=ithread, nsteps, nthreads
            tmp(ithread) = tmp(ithread) + (g_eval((i-0.5d0)*dt) + g_eval(-(i-0.5d0)*dt))
            if( abs(real(tmp(ithread))) + abs(aimag(tmp(ithread))) .gt. fcut(ithread))then
               int_thread(ithread) = int_thread(ithread)  + tmp(ithread)
               tmp(ithread) = 0.d0
               fcut(ithread) = (abs(real(int_thread(ithread))) + abs(aimag(int_thread(ithread))))*1.d-7
            endif
         enddo
      enddo
      !$omp end parallel do
      intall = (sum(tmp) + sum(int_thread)) * (dt/sumw)
    end function int_grid_odd
    
  end function coop_periodsym_int




    !!integrate \int_{-a}^a e^{zt + f(t)} dt 
  function coop_periodsym_intln(f, a, z, accuracy) result(pint)
    external f
    COOP_INT, parameter::min_level = 3, max_level = 11    
    COOP_REAL::a,  f
    COOP_REAL::x, y, absy
    COOP_COMPLEX :: z,  pint, intres(min_level-1:max_level)
    COOP_REAL:: p1, p1x, p2, p2x, p3, p3x, acc
    COOP_INT::level, np
    COOP_REAL,optional::accuracy
    COOP_INT, parameter::sumw = 16  !! 1+2+1 = 4 (2nd order algorithm) or 1+4+6+4+1 = 16 (4th order algorithm) 1 + 6 + 15 + 20 + 15 +  6 + 1 = 64
    if(present(accuracy))then
       acc = accuracy
    else
       acc = 1.d-5
    endif
    x = real(z)
    y = aimag(z)
    absy = abs(y)
    if(absy*a .lt. 2.d0)then
       write(*,*) "periodsym_intln is designed for periodic function with the imaginary part of z not too small"
       stop
    endif
    np = nint(a/p1)
    p1x = p1*x
    if(sumw > 4)then
       p2 = p1 * 2
       p2x = p2*x
    endif
    if(sumw > 16)then
       p3 = p1 * 3
       p3x = p3 * x
    endif
    level = min_level
    intres(level) = int_grid(p1/2**level,  np * 2**level) 
    do 
       level = level + 1       
       intres(level) = (int_grid_odd(p1/2**(level-1), np * 2**level) + intres(level-1))/2.d0
       if(abs(intres(level) - intres(level-1)) .lt. acc*abs(intres(level)) .or. level .eq. max_level)then
          pint = intres(level)
          return
       endif
    enddo
  contains

    function g_eval(t) result(g)
      COOP_COMPLEX::g, gl, gr
      COOP_REAL::t
      select case(sumw)
      case(4)
         g = (2.d0 - exp(f(t - p1) - f(t) - p1x) - exp(f(t + p1)-f(t) + p1x)) * exp(f(t) + z*t)  !2nd order algorithm
      case(16)
         g = (6.d0 - 4*(exp(f(t - p1) - f(t) - p1x) + exp(f(t + p1)-f(t) + p1x)) + (exp(f(t - p2) - f(t) - p2x) + exp(f(t + p2)-f(t) + p2x)) ) * exp(f(t) + z*t)  !4th order algorithm
      case(64)
         g = (20.d0 - 15*(exp(f(t - p1) - f(t) - p1x) + exp(f(t + p1)-f(t) + p1x)) + 6*(exp(f(t - p2) - f(t) - p2x) + exp(f(t + p2)-f(t) + p2x)) - (exp(f(t - p3) - f(t) - p3x) + exp(f(t + p3)-f(t) + p3x))  ) * exp(f(t) + z*t)  !6th order algorithm
      end select
    end function g_eval
    
    
    function int_grid(dt, nsteps) result(intall)
      COOP_INT,parameter::nthreads = 64
      COOP_COMPLEX::intall, tmp(nthreads), int_thread(nthreads)
      COOP_REAL::dt, fcut(nthreads)
      COOP_INT::nsteps, i, ithread
      tmp = 0.d0
      int_thread = 0.d0
      fcut = 0.d0
      !$omp parallel do
      do ithread = 1, nthreads
         do i=ithread, nsteps, nthreads
            tmp(ithread) = tmp(ithread) + (g_eval(i*dt) + g_eval(-i*dt))
            if( abs(real(tmp(ithread))) + abs(aimag(tmp(ithread))) .gt. fcut(ithread))then
               int_thread(ithread) = int_thread(ithread)  + tmp(ithread)
               tmp(ithread) = 0.d0
               fcut(ithread) = (abs(real(int_thread(ithread))) + abs(aimag(int_thread(ithread))))*1.d-7
            endif
         enddo
      enddo
      !$omp end parallel do
      intall = (g_eval(0.d0) + sum(tmp) + sum(int_thread)) * (dt/sumw)
    end function int_grid


    function int_grid_odd(dt, nsteps) result(intall)
      COOP_INT,parameter::nthreads = 64
      COOP_COMPLEX::intall, tmp(nthreads), int_thread(nthreads)
      COOP_REAL::dt, fcut(nthreads)
      COOP_INT::nsteps, i, ithread
      tmp = 0.d0
      int_thread = 0.d0
      fcut = 0.d0
      !$omp parallel do
      do ithread = 1, nthreads
         do i=ithread, nsteps, nthreads
            tmp(ithread) = tmp(ithread) + (g_eval((i-0.5d0)*dt) + g_eval(-(i-0.5d0)*dt))
            if( abs(real(tmp(ithread))) + abs(aimag(tmp(ithread))) .gt. fcut(ithread))then
               int_thread(ithread) = int_thread(ithread)  + tmp(ithread)
               tmp(ithread) = 0.d0
               fcut(ithread) = (abs(real(int_thread(ithread))) + abs(aimag(int_thread(ithread))))*1.d-7
            endif
         enddo
      enddo
      !$omp end parallel do
      intall = (sum(tmp) + sum(int_thread)) * (dt/sumw)
    end function int_grid_odd
    
    
  end function coop_periodsym_intln
  
  
end module coop_expint_mod
