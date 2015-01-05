module coop_gaussian_peak_stat_mod
  use coop_wrapper_typedef
  use coop_special_function_mod
  use coop_integrate_mod
  use coop_file_mod
  use coop_random_mod
  implicit none
#include  "constants.h"

#define GF_DIM  args%i(1)
#define GF_SIGMA0  args%r(1)
#define GF_SIGMA1 args%r(2)
#define GF_SIGMA2  args%r(3)
#define GF_SIN_BETA args%r(4)
#define GF_COS_BETA args%r(5)
#define GF_CONST args%r(6)

  

contains

  subroutine coop_gaussian_npeak_set_args(args, dim, sigma0, sigma1, sigma2)
    COOP_INT dim
    COOP_REAL sigma0, sigma1, sigma2, gam
    type(coop_arguments)::args
    gam = sigma0 * sigma2
    if(gam .le. sigma1**2)call coop_return_error("coop_Gaussian_npeak_set_args", "gamma>=1", "stop")
    gam = sigma1**2/gam
    call args%init( i = (/ dim /), r = (/ sigma0, sigma1, sigma2, sqrt(1.d0-gam**2), gam, (sigma2/(sqrt(dble(dim))*sigma1))**dim/coop_2pi**(dim/2.d0+1) /) )
  end subroutine coop_gaussian_npeak_set_args

  function  coop_gaussian_peak_f(v, args) result(f)
    COOP_REAL v, f
    COOP_REAL,parameter::sqrt5h = sqrt(2.5d0)
    type(coop_arguments)::args
    select case(GF_DIM)
    case(1)
       f = abs(v)
    case(2)
       f = v**2-1.d0 + exp(-v**2)
    case(3)
       f = (3.d0 - v**2)*v*(erf(-sqrt5h * v) + erf(-sqrt5h/2.d0 * v))/2.d0 + (1.d0/sqrt5h/coop_sqrtpi)*((31.d0/4.d0 * v**2 + 8.d0/5.d0)*exp(-5.d0/8.d0*v**2) + (v**2/2.d0-8.d0/5.d0)*exp(-2.5d0*v**2))
    case default
       stop "coop_gaussian_peak_f only support dim = 1, 2, 3"
    end select
  end function coop_gaussian_peak_f
    
  !!u, v, dnpk
  function coop_gaussian_npeak_differential(u, v, args) result(dnpk)
    COOP_REAL u, v, dnpk
    type(coop_arguments)::args
    dnpk = GF_CONST/GF_SIN_BETA * coop_gaussian_peak_f(v, args) * exp(- (u**2 + ((v + u*GF_COS_BETA)/GF_SIN_BETA)**2)/2.d0)
  end function coop_gaussian_npeak_differential


  function coop_gaussian_nmax_differential(u, args) result(dnpk)
    COOP_REAL u, v, dnpk
    type(coop_arguments)::args
    dnpk = GF_CONST * coop_gaussian_peak_G(u, args) 
  end function coop_gaussian_nmax_differential

  function coop_gaussian_peak_G(u, args) result(g)
    COOP_REAL u, g, t, q
    type(coop_arguments)::args    
    select case(GF_DIM)
    case(1)
       t = u * GF_COS_BETA/GF_SIN_BETA/coop_sqrt2       
       g = exp(- t**2)*GF_SIN_BETA + sqrt(coop_pi/2.d0)*u * GF_COS_BETA*(1.d0 + erf(t))
    case(2)
       t = u * GF_COS_BETA/GF_SIN_BETA/coop_sqrt2       
       q = sqrt(3.d0-2.d0*GF_COS_BETA**2)
       g = sqrt(coop_pio2)/q*(1.d0+erf(t/q))*exp(-(u*GF_COS_BETA/q)**2) &
            + sqrt(coop_pio2)*(u**2-1.d0)*GF_COS_BETA**2*(1.d0+erf(t)) &
            + u*GF_SIN_BETA*GF_COS_BETA*exp(-t**2)
    case(3)
       t = -u*GF_COS_BETA
       q = 5.d0*GF_SIN_BETA
       if(t + q .lt. 0.d0 .and. t .lt. -3.d0)then
          g = sqrt(coop_2pi)*(3.d0*GF_COS_BETA**2 - t**2)*t
       else
          g = coop_integrate(dgdv, min(t, 0.d0)-q, min(0.d0, t+q), 1.d-7) / GF_SIN_BETA
       endif
    case default
       stop "gaussian_peak_G only support dim = 1, 2, 3"
    end select
    g= g*exp(-u**2/2.d0)
  contains
    function dGdv(v)
      COOP_REAL dGdv, v
      dGdv =  coop_gaussian_peak_f(v, args) * exp(- (((v + u*GF_COS_BETA)/GF_SIN_BETA)**2)/2.d0)
    end function dGdv
  end function coop_gaussian_peak_G

  function coop_gaussian_peak_G_slow(u, args) result(g)
    COOP_REAL t, q, g, u
    type(coop_arguments)::args    
    t = -u*GF_COS_BETA
    q = 6.d0*GF_SIN_BETA
    g = coop_integrate(dgdv, min(t, 0.d0)-q, min(0.d0, t+q), 1.d-8) / GF_SIN_BETA*exp(-u**2/2.d0)
  contains
    function dGdv(v)
      COOP_REAL dGdv, v
      dGdv =  coop_gaussian_peak_f(v, args) * exp(- (((v + u*GF_COS_BETA)/GF_SIN_BETA)**2)/2.d0)
    end function dGdv
  end function coop_gaussian_peak_G_slow


  function coop_gaussian_peak_Gv(u, args) result(g)
    COOP_REAL u, g, t, q
    type(coop_arguments)::args
    
    select case(GF_DIM)
    case(1)
       q = u * GF_COS_BETA       
       t = q/GF_SIN_BETA/coop_sqrt2
       g =  - sqrt(coop_pio2)*(q**2 + GF_SIN_BETA**2)*(1.d0 + Erf(t)) &
            - exp(-t**2)*q*GF_SIN_BETA
    case(2)
       t = u * GF_COS_BETA/GF_SIN_BETA/coop_sqrt2       
       q = sqrt(3.d0-2.d0*GF_COS_BETA**2)
       g =  - sqrt(coop_pio2)/q**3*u*GF_COS_BETA*(1.d0+erf(t/q))*exp(-(u*GF_COS_BETA/q)**2) &
            - sqrt(coop_pio2)*u*GF_COS_BETA*(3.d0*GF_SIN_BETA**2+(u*GF_COS_BETA)**2-1.d0)*(1.d0+erf(t)) &
            - GF_SIN_BETA*exp(-t**2)*(1.d0/q**2 - 1 + 2.d0*GF_SIN_BETA**2 + (u*GF_COS_BETA)**2)
    case(3)
       t = -u*GF_COS_BETA
       q = 5.d0*GF_SIN_BETA
       if(t + q .lt. 0.d0 .and. t .lt. -3.d0)then
          g = -sqrt(coop_2pi)*t**2*(t**2+(6.d0*GF_SIN_BETA**2-3.d0))
       else
          g = coop_integrate(dgdv, min(t, 0.d0)-q, min(0.d0, t+q), 1.d-7) / GF_SIN_BETA
       endif
    case default
       stop "gaussian_peak_G only support dim = 1, 2, 3"
    end select
    g= g*exp(-u**2/2.d0)
  contains
    function dGdv(v)
      COOP_REAL dGdv, v
      dGdv =  coop_gaussian_peak_f(v, args) * v * exp(- (((v + u*GF_COS_BETA)/GF_SIN_BETA)**2)/2.d0)
    end function dGdv
  end function coop_gaussian_peak_Gv


  function coop_gaussian_peak_Gv_slow(u, args) result(g)
    COOP_REAL t, q, g, u
    type(coop_arguments)::args    
    t = -u*GF_COS_BETA
    q = 6.d0*GF_SIN_BETA
    g = coop_integrate(dgdv, min(t, 0.d0)-q, min(0.d0, t+q), 1.d-8) / GF_SIN_BETA*exp(-u**2/2.d0)
  contains
    function dGdv(v)
      COOP_REAL dGdv, v
      dGdv =  coop_gaussian_peak_f(v, args) * v * exp(- (((v + u*GF_COS_BETA)/GF_SIN_BETA)**2)/2.d0)
    end function dGdv
  end function coop_gaussian_peak_Gv_slow  
  

  function coop_gaussian_peak_Gu(u, args) result(g)
    COOP_REAL u, g
    type(coop_arguments)::args
    g= u * coop_gaussian_peak_G(u, args)
  end function coop_gaussian_peak_Gu  
  

  !!before marginalization of e
  function coop_gaussian_npeak_differential_2D_full(u, v, e, args) result(dnpk)
    COOP_REAL u, v, e, dnpk
    type(coop_arguments)::args
    dnpk = (GF_CONST*2.d0) * v**4*e*(1.d0-e**2)*exp(-(v*e)**2 - (u**2 + ((v + u*GF_COS_BETA)/GF_SIN_BETA)**2)/2.d0) 
  end function coop_gaussian_npeak_differential_2D_full

  function coop_gaussian_nmax(nu, args) result(nmax)
    COOP_REAL nmax, nu
    type(coop_arguments)::args
    if(nu .gt. 3.d0)then
       nmax = GF_CONST*coop_integrate(coop_gaussian_peak_G, nu, nu + 3.d0, args, 1.d-7)       
    else
       nmax = GF_CONST*coop_integrate(coop_gaussian_peak_G, max(nu, -5.d0), 6.d0, args, 1.d-7)
    endif
  end function coop_gaussian_nmax

  function coop_gaussian_peak_intu(nu, args) result(intu)
    COOP_REAL intu, nu
    type(coop_arguments)::args
    if(nu .gt. 3.d0)then
       intu = GF_CONST*coop_integrate(coop_gaussian_peak_Gu, nu, nu + 3.d0, args, 1.d-7)       
    else
       intu = GF_CONST*coop_integrate(coop_gaussian_peak_Gu, max(nu, -5.d0), 6.d0, args, 1.d-7)
    endif

  end function coop_gaussian_peak_intu

  function coop_gaussian_peak_intv(nu, args) result(intv)
    COOP_REAL intv, nu
    type(coop_arguments)::args
    if(nu .gt. 3.d0)then
       intv = GF_CONST*coop_integrate(coop_gaussian_peak_Gv, nu, nu + 3.d0, args, 1.d-7)       
    else
       intv = GF_CONST*coop_integrate(coop_gaussian_peak_Gv, max(nu, -5.d0), 6.d0, args, 1.d-7)
    endif
  end function coop_gaussian_peak_intv
  
  subroutine coop_gaussian_get_nonoriented_stacking_weights(nu, args, weights)
    type(coop_arguments)::args
    COOP_REAL nu
    COOP_REAL:: A(2,2), nmax, weights(4)
    A(1,1) = 1.d0
    A(2,2) = 1.d0
    A(1,2) = GF_COS_BETA
    A(2,1) = GF_COS_BETA
    nmax = coop_gaussian_nmax(nu, args)
    weights(1) = coop_gaussian_peak_intu(nu, args)/nmax
    weights(2) = coop_gaussian_peak_intv(nu, args)/nmax
    weights(1:2) = matmul(A, weights(1:2))/GF_SIN_BETA**2
    weights(1) = weights(1)/GF_SIGMA0
    weights(2) = weights(2)/GF_SIGMA2
    weights(3:4) = 0.d0
  end subroutine coop_gaussian_get_nonoriented_stacking_weights


  subroutine coop_gaussian_get_oriented_stacking_weights(nu, args, weights)
    type(coop_arguments)::args
    COOP_REAL nu
    COOP_REAL:: A(2,2), nmax, weights(4)
    A(1,1) = 1.d0
    A(2,2) = 1.d0
    A(1,2) = GF_COS_BETA
    A(2,1) = GF_COS_BETA
    call coop_gaussian_peak_2Dmax_mean(4, weights, get_iq, args)
    weights(1:2) = matmul(A, weights(1:2))/GF_SIN_BETA**2
    weights(1) = weights(1)/GF_SIGMA0
    weights(2) = weights(2)/GF_SIGMA2
    weights(3:4) = matmul(A, weights(3:4))*(coop_sqrt2/GF_SIN_BETA**2)
    weights(3) = weights(3)/GF_SIGMA0
    weights(4) = weights(4)/GF_SIGMA2
  contains
    subroutine get_iq(i,q,u,args,s,want)
      COOP_REAL::i(0:1),q(0:1),u(0:1)
      type(coop_arguments)::args
      COOP_REAL::s(4)
      logical want
      want = ( i(0).gt.nu )
      if(want)then
         s(1:2) = i
         !!rotate by 1/2arctan(u(0)/q(0))         
         s(3) = sqrt(q(0)**2 + u(0)**2) 
         s(4) = (q(1)*q(0)+u(1)*u(0))/s(3)
      endif
    end subroutine get_iq
  end subroutine coop_gaussian_get_oriented_stacking_weights


    subroutine coop_gaussian_get_pmax_stacking_weights(nu, args, weights)
    type(coop_arguments)::args
    COOP_REAL nu
    COOP_REAL:: A(2,2), nmax, weights(4), twonu2
    A(1,1) = 1.d0
    A(2,2) = 1.d0
    A(1,2) = GF_COS_BETA
    A(2,1) = GF_COS_BETA
    twonu2 = 2.d0*nu**2
    call coop_gaussian_peak_pmax_mean(2, weights(3:4), getp, args)
    weights(1:2) = 0.d0
    weights(3:4) = matmul(A, weights(3:4))*(coop_sqrt2/GF_SIN_BETA**2)
    weights(3) = weights(3)/GF_SIGMA0
    weights(4) = weights(4)/GF_SIGMA2
  contains
    subroutine getp(q, u, qx, qy, qq, uu, qu, uq, args, s, want)
      COOP_REAL:: q(0:1), u(0:1), qx, qy, qq, uu, qu, uq
      type(coop_arguments)::args
      COOP_REAL::s(2), psq
      logical want
      psq = q(0)**2+u(0)**2
      want =  (psq .gt. twonu2)
      if(want)then
         s(1) = sqrt(psq)
         s(2) = (q(1)*q(0)+u(1)*u(0))/s(1)
      endif
    end subroutine getp
  end subroutine coop_gaussian_get_pmax_stacking_weights


  subroutine coop_gaussian_radial_modes_I( weights, cr, frI )
    COOP_REAL weights(4), cr(0:1, 0:2), frI(0:2)
    frI(2) = 0.d0
    frI(1) = weights(3)*cr(0,1) + weights(4)*cr(1,1)
    frI(0) = weights(1)*cr(0,0) + weights(2)*cr(1,0)
  end subroutine coop_gaussian_radial_modes_I

  subroutine coop_gaussian_radial_modes_QU( weights, cr, frQU )
    COOP_REAL weights(4), cr(0:1, 0:2), frQU(0:2)
    frQU(0) = (weights(3)*cr(0,0)+weights(4)*cr(1,0))/2.d0
    frQU(1) = weights(1)*cr(0,1) + weights(2)*cr(1,1)
    frQU(2) = (weights(3)*cr(0,2) + weights(4)*cr(1,2))/2.d0
  end subroutine coop_gaussian_radial_modes_QU

  !!compute the mean of some function of i,q,u,args
  !!==============================================  
  !!subroutine f(i, q, u, args, s, want)
  !!  COOP_REAL i(2),q(2),u(2), s(nf)
  !!  type(coop_arguments)::args
  !!  logical want
  !!-------------------
  !!  input i,qu,args
  !!  return s, want
  !!==============================================
  subroutine coop_gaussian_peak_2Dmax_mean(nf, mean, f, args)
    COOP_INT,parameter::num_sims = 80000
    COOP_INT,parameter::max_try = 100*num_sims
    type(coop_arguments)::args
    COOP_INT::nf
    COOP_REAL:: mean(nf)
    external f
    COOP_REAL::i(0:1), q(0:1), u(0:1), s(nf), vol, dvol
    logical want
    type(coop_file)::fp
    COOP_INT::iaccept, itry
    iaccept = 0
    itry = 0
    vol = 0.d0
    mean = 0.d0
    do while(iaccept.lt.num_sims)
       itry = itry+1
       if(itry.gt.max_try)then
          if(iaccept .eq. 0)then
             stop "gaussian_peak_2Dmax_mean: cannot find constrained samples"
          else
             if(iaccept .lt. 100)then
                write(*,*) "Warning: gaussian_peak_2Dmax_mean too few samples, the result might not accurate"
             endif
          endif
          exit
       endif
       call coop_gaussian_get_peak_samples_i(i, args) !!only take i(1) negative
       call coop_gaussian_get_peak_samples_general(q, args)
       call coop_gaussian_get_peak_samples_general(u, args)
       dvol = i(1)**2-(q(1)**2+u(1)**2)/2.d0       
       if(dvol .gt. 0.d0)then  !!maxima, accepted
          call f(i,q,u,args,s,want)
          if(want)then
             iaccept = iaccept + 1
             vol = vol + dvol
             mean = mean + s*dvol
          endif
       endif
    enddo
    mean = mean/vol
  end subroutine coop_gaussian_peak_2Dmax_mean



  !!compute the mean of some function of i,q,u,args
  !!==============================================  
  !!subroutine f(i, q, u, args, s, want)
  !!  COOP_REAL i(2),q(2),u(2), s(nf)
  !!  type(coop_arguments)::args
  !!  logical want
  !!-------------------
  !!  input i,qu,args
  !!  return s, want
  !!==============================================
  subroutine coop_gaussian_peak_Pmax_mean(nf, mean, f, args)
    COOP_INT,parameter::num_sims = 100000
    COOP_INT,parameter::max_try = 100*num_sims
    type(coop_arguments)::args
    COOP_INT::nf
    COOP_REAL:: mean(nf)
    external f
    COOP_REAL::q(0:1), u(0:1), Qx, Qy, qq, qu, uq, uu, s(nf), vol, dvol, delta1, delta2, scal, psq
    logical want
    type(coop_file)::fp
    COOP_INT::iaccept, itry
    iaccept = 0
    itry = 0
    vol = 0.d0
    mean = 0.d0
    do while(iaccept.lt.num_sims)
       itry = itry+1
       if(itry.gt.max_try)then
          if(iaccept .eq. 0)then
             stop "gaussian_peak_2Dmax_mean: cannot find constrained samples"
          else
             if(iaccept .lt. 100)then
                write(*,*) "Warning: gaussian_peak_2Dmax_mean too few samples, the result might not accurate"
             endif
          endif
          exit
       endif
       call coop_gaussian_get_peak_samples_general(q, args)
       call coop_gaussian_get_peak_samples_general(u, args)
       call coop_random_get_gaussian_pair(qq, uu)
       call coop_random_get_gaussian_pair(qu, uq)
       call coop_random_get_gaussian_pair(qx, qy)
       psq = q(0)**2+u(0)**2
       scal = u(0)/sqrt(psq)
       qx = qx*scal
       qy = qy*scal
       dvol = abs(u(1)**2 - (qu**2+3.d0*uu**2)/4.d0)*scal**2
       delta1 = qx**2*psq*GF_COS_BETA + u(0)**2*(q(0)*(q(1)+coop_sqrt3/2.d0*qq)+u(0)*(u(1)+qu/2.d0))
       if(delta1 .lt. 0.d0)then  
          delta2 = delta1*(qy**2*psq*GF_COS_BETA + u(0)**2*(q(0)*(q(1)-coop_sqrt3/2.d0*qq)+u(0)*(u(1)-qu/2.d0))) - (qx*qy*psq*GF_COS_BETA+u(0)**2/2.d0*(q(0)*uq+coop_sqrt3*u(0)*uu))**2
          if(delta2 .gt. 0.d0)then
             call f(q,u,qx, qy, qq, uu, qu, uq, args,s,want)
             if(want)then
                iaccept = iaccept + 1
                vol = vol + dvol
                mean = mean + s*dvol
             endif
          endif
       endif
    enddo
    mean = mean/vol
  end subroutine coop_gaussian_peak_Pmax_mean
  


  subroutine coop_gaussian_get_peak_samples_general(v, args)
    COOP_REAL r, v(2), s
    type(coop_arguments)::args
100 call random_number(v)
    v = 2.d0*v - 1.d0
    r=v(1)**2 + v(2)**2
    if(r.ge. 1.d0) goto 100
    v = v * dsqrt(-2.d0*dlog(r)/r)
    v(1) = v(1)*GF_SIN_BETA - v(2)*GF_COS_BETA
  end subroutine coop_gaussian_get_peak_samples_general

  subroutine coop_gaussian_get_peak_samples_i(v, args)  !!v(2) is always negative in this case
    COOP_REAL r, v(2), s
    type(coop_arguments)::args
100 call random_number(v)
    v(1) = 2.d0*v(1) - 1.d0
    v(2) = v(2)-1.d0
    r=v(1)**2 + v(2)**2
    if(r.ge. 1.d0) goto 100
    v = v * dsqrt(-2.d0*dlog(r)/r)
    v(1) = v(1)*GF_SIN_BETA - v(2)*GF_COS_BETA
  end subroutine coop_gaussian_get_peak_samples_i

  
  
end module coop_gaussian_peak_stat_mod
