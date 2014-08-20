module coop_ode_mod
  use coop_wrapper_typedef
  implicit none

#include "constants.h"

  private
  integer,parameter::dl = coop_real_length
  COOP_REAL,parameter::coop_ode_default_precision = 3.d-7
  
  public::coop_ode, coop_symplectic_2nd, coop_symplectic_4th, coop_symplectic_6th, coop_dverk, coop_dverk_with_arguments, coop_dverk_with_cosmology

  type coop_ode
     COOP_INT  method
     type(coop_arguments)::args
     integer n, ind
     logical has_args
     COOP_REAL c(24)
     COOP_REAL,dimension(:),allocatable:: y, yp
     COOP_REAL,dimension(:,:),allocatable:: w
     COOP_REAL  x, tol
   contains
     procedure::init => coop_ode_init
     procedure::set_arguments => coop_ode_set_arguments
     procedure::set_initial_conditions => coop_ode_set_initial_conditions
     procedure::evolve => coop_ode_evolve
  end type coop_ode


contains

  subroutine coop_ode_init(this, n, method, tol)
    class(coop_ode)::this
    COOP_INT n
    COOP_INT,optional::method
    COOP_REAL,optional::tol
    this%n = n
    if(present(method))then
       this%method = method
    else
       this%method = COOP_ODE_DVERK
    end if
    this%c = 0
    this%ind = 1
    if(allocated(this%y))then
       if(size(this%y).ne.n)then
          deallocate(this%w, this%yp, this%y)
          allocate(this%w(n, 9), this%y(n), this%yp(n))
       endif
    else
       allocate(this%w(n, 9), this%y(n), this%yp(n))
    endif
    this%w = 0
    this%has_args = .false.
    if(present(tol))then
       this%tol = tol
    else
       this%tol = coop_ode_default_precision
    endif
  end subroutine coop_ode_init

  subroutine coop_ode_set_arguments(this, args)
    class(coop_ode)::this
    type(coop_arguments)::args
    this%args = args
    this%has_args = .true.
  end subroutine coop_ode_set_arguments

  subroutine coop_ode_set_initial_conditions(this, xini, yini)
    class(coop_ode)::this
    COOP_REAL xini
    COOP_REAL,dimension(:)::yini
    this%x = xini
    this%ind = 1
    if(size(yini).ne. this%n)then
       stop "error in set_initial_condtions: you need to initialize the ode object with correct dimension"
    endif
    this%y = yini
  end subroutine coop_ode_set_initial_conditions


  subroutine coop_ode_evolve(this, fcn, xend)
    class(coop_ode)::this
    external fcn
    COOP_REAL xend
    select case(this%method)
    case(COOP_ODE_DVERK)
       if(this%has_args)then
          call coop_dverk_with_arguments(this%n, fcn, this%args, this%x, this%y, xend, this%tol, this%ind, this%c, this%n, this%w)
       else
          call coop_dverk(this%n, fcn, this%x, this%y, xend, this%tol, this%ind, this%c, this%n, this%w)
       endif
    case(COOP_ODE_RK4)
       if(this%has_args)then
          call coop_RungeKutta4th_with_arguments(this%n, fcn, this%args, this%x, this%y, xend - this%x)
       else
          call coop_RungeKutta4th(this%n, fcn, this%x, this%y, xend - this%x)
       end if
    case(COOP_ODE_RK6)
       if(this%has_args)then
          call coop_RungeKutta6th_with_arguments(this%n, fcn, this%args, this%x, this%y, xend - this%x)
       else
          call coop_RungeKutta6th(this%n, fcn, this%x, this%y, xend - this%x)
       end if
    case(COOP_ODE_RK8)
       if(this%has_args)then
          call coop_RungeKutta8th_with_arguments(this%n, fcn, this%args, this%x, this%y, xend - this%x)
       else
          call coop_RungeKutta8th(this%n, fcn, this%x, this%y, xend - this%x)
       end if
    case(COOP_ODE_GL4)
       if(this%has_args)then
          call gauss_legendre_4th_with_arguments(this%n, fcn, this%args, this%x, this%y, xend - this%x)
       else
          call gauss_legendre_4th(this%n, fcn, this%x, this%y, xend - this%x)
       end if       
    case(COOP_ODE_GL6)
       if(this%has_args)then
          call gauss_legendre_6th_with_arguments(this%n, fcn, this%args, this%x, this%y, xend - this%x)
       else
          call gauss_legendre_6th(this%n, fcn, this%x, this%y, xend - this%x)
       end if       
    case(COOP_ODE_GL8)
       if(this%has_args)then
          call gauss_legendre_8th_with_arguments(this%n, fcn, this%args, this%x, this%y, xend - this%x)
       else
          call gauss_legendre_8th(this%n, fcn, this%x, this%y, xend - this%x)
       end if
    case default
       write(*,*) "ODE_METHOD = ", this%method
       stop "This ode method is not supported yet"
    end select
  end subroutine coop_ode_evolve



  !!fixed step Runge-Kutta ODE integrator
  subroutine coop_RungeKutta4th(n, fcn, t, y, h)
    external fcn
    !! fcn(n, x, y(1:n), yp(1:n))

    COOP_INT  n
    COOP_REAL y(n)
    COOP_REAL t, h
    COOP_REAL k1(n), k2(n), k3(n), k4(n)
    call fcn(n, t, y, k1)
    call fcn(n, t+h/2._dl, y+k1*(h/2._dl), k2)
    call fcn(n, t+h/2._dl, y+k2*(h/2._dl), k3)
    t = t + h
    call fcn(n, t, y+k3*h, k4)
    y = y+(k1+2.*(k2+k3)+k4)*(h/6._dl)
  end subroutine coop_RungeKutta4th

  !!fixed step Runge-Kutta ODE integrator
  subroutine coop_RungeKutta6th(n, fcn, t, y, h)
    external fcn
    !! fcn(n, t, y(1:n), yp(1:n))
    COOP_REAL,parameter::sqrt21 = 4.5825756949558400066_dl
    COOP_INT  n
    COOP_REAL y(n)
    COOP_REAL t, h
    COOP_REAL k1(n), k2(n), k3(n), k4(n), k5(n), k6(n), k7(n)
    call fcn(n, t, y, k1)  !!0
    call fcn(n, t+h, y+k1*h, k2)  !! 1
    call fcn(n, t+h/2._dl, y+(3._dl*k1+k2)*(h/8._dl), k3) !!0.5
    call fcn(n, t+(2._dl/3._dl)*h, y+ (4._dl*(k1+k3)+k2)*((2._dl/27._dl)*h), k4) !! 0.67
    call fcn(n, t+((7._dl-sqrt21)/14._dl)*h, y+((3._dl*(3._dl*sqrt21 - 7._dl)/392._dl)*k1 - (8._dl*(7._dl-sqrt21)/392._dl)*k2 + (48._dl*(7._dl-sqrt21)/392._dl)*k3 - (3._dl*(21._dl - sqrt21)/392._dl)*k4)*h, k5)  !! 0.17
    call fcn(n, t+((7._dl+sqrt21)/14._dl)*h, y+ ((-5._dl*(231._dl+51._dl*sqrt21)/1960._dl)*k1 - (40._dl*(7._dl+sqrt21)/1960._dl)*k2 - (320._dl*sqrt21/1960._dl)*k3 + (3._dl*(21._dl + 121._dl*sqrt21)/1960._dl)*k4 + (392._dl*(6._dl+sqrt21)/1960._dl)*k5)*h, k6)  !!0.8
    t = t+h
    call fcn(n, t, y+(((22._dl+7._dl*sqrt21)/12._dl)*k1 +(2._dl/3._dl)*k2 + (2._dl*(7._dl*sqrt21-5._dl)/9._dl)*k3 - (7._dl*(3._dl*sqrt21-2._dl)/20._dl)*k4 - (7._dl*(49._dl+9._dl*sqrt21)/90._dl)*k5 + (7._dl*(7._dl-sqrt21)/18._dl)*k6)*h, k7) !! 1
    y = y+  ((k1+k7)/20._dl + (49._dl/180._dl)*(k5+k6) + (16._dl/45._dl)*k3)*h
  end subroutine coop_RungeKutta6th
  
  subroutine coop_RungeKutta8th (n, fcn, x0, y, h, err)
    COOP_REAL,parameter:: c1 = 103.d0 / 1680.d0
    COOP_REAL,parameter:: c8 = -27.d0 / 140.d0
    COOP_REAL,parameter:: c9 = 76.d0 / 105.d0
    COOP_REAL,parameter:: c10 = -201.d0 / 280.d0
    COOP_REAL,parameter:: c11 = 1024.d0 / 1365.d0
    COOP_REAL,parameter:: c12 = 3.d0 / 7280.d0
    COOP_REAL,parameter:: c13 = 12.d0 / 35.d0
    COOP_REAL,parameter:: c14 = 9.d0 / 280.d0
    COOP_REAL,parameter:: a2 = 1.d0 / 12.d0
    COOP_REAL,parameter:: a3 = 1.d0 / 9.d0
    COOP_REAL,parameter:: a4 = 1.d0 / 6.d0
    COOP_REAL,parameter:: a5 = 2.d0 * (1.d0 + coop_sqrt6) / 15.d0
    COOP_REAL,parameter:: a6 = (6.d0 + coop_sqrt6) / 15.d0
    COOP_REAL,parameter:: a7 = (6.d0 - coop_sqrt6) / 15.d0
    COOP_REAL,parameter:: a8 = 2.d0 / 3.d0
    COOP_REAL,parameter:: a9 = 1.d0 / 2.d0
    COOP_REAL,parameter:: a10 = 1.d0 / 3.d0
    COOP_REAL,parameter:: a11 = 1.d0 / 4.d0
    COOP_REAL,parameter:: a12 = 4.d0 / 3.d0
    COOP_REAL,parameter:: a13 = 5.d0 / 6.d0
    COOP_REAL,parameter:: a15 = 1.d0 / 6.d0
    COOP_REAL,parameter:: b31 = 1.d0 / 27.d0
    COOP_REAL,parameter:: b32 = 2.d0 / 27.d0
    COOP_REAL,parameter:: b41 = 1.d0 / 24.d0
    COOP_REAL,parameter:: b43 = 3.d0 / 24.d0
    COOP_REAL,parameter:: b51 = (4.d0 + 94.d0 * coop_sqrt6) / 375.d0
    COOP_REAL,parameter:: b53 = -(282.d0 + 252.d0 * coop_sqrt6) / 375.d0
    COOP_REAL,parameter:: b54 = (328.d0 + 208.d0 * coop_sqrt6) / 375.d0
    COOP_REAL,parameter:: b61 = (9.d0 - coop_sqrt6) / 150.d0
    COOP_REAL,parameter:: b64 = (312.d0 + 32.d0 * coop_sqrt6) / 1425.d0
    COOP_REAL,parameter:: b65 = (69.d0 + 29.d0 * coop_sqrt6) / 570.d0
    COOP_REAL,parameter:: b71 = (927.d0 - 347.d0 * coop_sqrt6) / 1250.d0
    COOP_REAL,parameter:: b74 = (-16248.d0 + 7328.d0 * coop_sqrt6) / 9375.d0
    COOP_REAL,parameter:: b75 = (-489.d0 + 179.d0 * coop_sqrt6) / 3750.d0
    COOP_REAL,parameter:: b76 = (14268.d0 - 5798.d0 * coop_sqrt6) / 9375.d0
    COOP_REAL,parameter:: b81 = 4.d0 / 54.d0
    COOP_REAL,parameter:: b86 = (16.d0 - coop_sqrt6) / 54.d0
    COOP_REAL,parameter:: b87 = (16.d0 + coop_sqrt6) / 54.d0
    COOP_REAL,parameter:: b91 = 38.d0 / 512.d0
    COOP_REAL,parameter:: b96 = (118.d0 - 23.d0 * coop_sqrt6) / 512.d0
    COOP_REAL,parameter:: b97 = (118.d0 + 23.d0 * coop_sqrt6) / 512.d0
    COOP_REAL,parameter:: b98 = - 18.d0 / 512.d0
    COOP_REAL,parameter:: b10_1 = 11.d0 / 144.d0
    COOP_REAL,parameter:: b10_6 = (266.d0 - coop_sqrt6) / 864.d0
    COOP_REAL,parameter:: b10_7 = (266.d0 + coop_sqrt6) / 864.d0
    COOP_REAL,parameter:: b10_8 = - 1.d0 / 16.d0
    COOP_REAL,parameter:: b10_9 = - 8.d0 / 27.d0
    COOP_REAL,parameter:: b11_1 = (5034.d0 - 271.d0 * coop_sqrt6) / 61440.d0
    COOP_REAL,parameter:: b11_7 = (7859.d0 - 1626.d0 * coop_sqrt6) / 10240.d0
    COOP_REAL,parameter:: b11_8 = (-2232.d0 + 813.d0 * coop_sqrt6) / 20480.d0
    COOP_REAL,parameter:: b11_9 = (-594.d0  + 271.d0 * coop_sqrt6) / 960.d0
    COOP_REAL,parameter:: b11_10 = (657.d0 - 813.d0 * coop_sqrt6) / 5120.d0
    COOP_REAL,parameter:: b12_1 = (5996.d0 - 3794.d0 * coop_sqrt6) / 405.d0
    COOP_REAL,parameter:: b12_6 = (-4342.d0 - 338.d0 * coop_sqrt6) / 9.d0
    COOP_REAL,parameter:: b12_7 = (154922.d0 - 40458.d0 * coop_sqrt6) / 135.d0
    COOP_REAL,parameter:: b12_8 = (-4176.d0 + 3794.d0 * coop_sqrt6) / 45.d0
    COOP_REAL,parameter:: b12_9 = (-340864.d0 + 242816.d0 * coop_sqrt6) / 405.d0
    COOP_REAL,parameter:: b12_10 = (26304.d0 - 15176.d0 * coop_sqrt6) / 45.d0
    COOP_REAL,parameter:: b12_11 = -26624.d0 / 81.d0
    COOP_REAL,parameter:: b13_1 = (3793.d0 + 2168.d0 * coop_sqrt6) / 103680.d0
    COOP_REAL,parameter:: b13_6 = (4042.d0 + 2263.d0 * coop_sqrt6) / 13824.d0
    COOP_REAL,parameter:: b13_7 = (-231278.d0 + 40717.d0 * coop_sqrt6) / 69120.d0
    COOP_REAL,parameter:: b13_8 = (7947.d0 - 2168.d0 * coop_sqrt6) / 11520.d0
    COOP_REAL,parameter:: b13_9 = (1048.d0 - 542.d0 * coop_sqrt6) / 405.d0
    COOP_REAL,parameter:: b13_10 = (-1383.d0 + 542.d0 * coop_sqrt6) / 720.d0
    COOP_REAL,parameter:: b13_11 = 2624.d0 / 1053.d0
    COOP_REAL,parameter:: b13_12 = 3.d0 / 1664.d0
    COOP_REAL,parameter:: b14_1 = -137.d0 / 1296.d0
    COOP_REAL,parameter:: b14_6 = (5642.d0 - 337.d0 * coop_sqrt6) / 864.d0
    COOP_REAL,parameter:: b14_7 = (5642.d0 + 337.d0 * coop_sqrt6) / 864.d0
    COOP_REAL,parameter:: b14_8 = -299.d0 / 48.d0
    COOP_REAL,parameter:: b14_9 = 184.d0 / 81.d0
    COOP_REAL,parameter:: b14_10 = -44.d0 / 9.d0
    COOP_REAL,parameter:: b14_11 = -5120.d0 / 1053.d0
    COOP_REAL,parameter:: b14_12 = -11.d0 / 468.d0
    COOP_REAL,parameter:: b14_13 = 16.d0 / 9.d0
    COOP_REAL,parameter:: b15_1 = (33617.d0 - 2168.d0 * coop_sqrt6) / 518400.d0
    COOP_REAL,parameter:: b15_6 = (-3846.d0 + 31.d0 * coop_sqrt6) / 13824.d0
    COOP_REAL,parameter:: b15_7 = (155338.d0 - 52807.d0 * coop_sqrt6) / 345600.d0
    COOP_REAL,parameter:: b15_8 = (-12537.d0 + 2168.d0 * coop_sqrt6) / 57600.d0
    COOP_REAL,parameter:: b15_9 = (92.d0 + 542.d0 * coop_sqrt6) / 2025.d0
    COOP_REAL,parameter:: b15_10 = (-1797.d0 - 542.d0 * coop_sqrt6) / 3600.d0
    COOP_REAL,parameter:: b15_11 = 320.d0 / 567.d0
    COOP_REAL,parameter:: b15_12 = -1.d0 / 1920.d0
    COOP_REAL,parameter:: b15_13 = 4.d0 / 105.d0
    COOP_REAL,parameter:: b16_1 = (-36487.d0 - 30352.d0 * coop_sqrt6) / 279600.d0
    COOP_REAL,parameter:: b16_6 = (-29666.d0 - 4499.d0 * coop_sqrt6) / 7456.d0
    COOP_REAL,parameter:: b16_7 = (2779182.d0 - 615973.d0 * coop_sqrt6) / 186400.d0
    COOP_REAL,parameter:: b16_8 = (-94329.d0 + 91056.d0 * coop_sqrt6) / 93200.d0
    COOP_REAL,parameter:: b16_9 = (-232192.d0 + 121408.d0 * coop_sqrt6) / 17475.d0
    COOP_REAL,parameter:: b16_10 = (101226.d0 - 22764.d0 * coop_sqrt6) / 5825.d0
    COOP_REAL,parameter:: b16_11 = - 169984.d0 / 9087.d0
    COOP_REAL,parameter:: b16_12 = - 87.d0 / 30290.d0
    COOP_REAL,parameter:: b16_13 =  492.d0 / 1165.d0
    COOP_REAL,parameter:: b16_15 =  1260.d0 / 233.d0
    COOP_REAL,parameter:: e1 = -1911.d0 / 109200.d0
    COOP_REAL,parameter:: e8 = 34398.d0 / 109200.d0
    COOP_REAL,parameter:: e9 = -61152.d0 / 109200.d0
    COOP_REAL,parameter:: e10 = 114660.d0 / 109200.d0
    COOP_REAL,parameter:: e11 = -114688.d0 / 109200.d0
    COOP_REAL,parameter:: e12 = -63.d0 / 109200.d0
    COOP_REAL,parameter:: e13 = -13104.d0 / 109200.d0
    COOP_REAL,parameter:: e14 = -3510.d0 / 109200.d0
    COOP_REAL,parameter:: e15 = 39312.d0 / 109200.d0
    COOP_REAL,parameter:: e16 = 6058.d0 / 109200.d0

    COOP_INT  n
    external fcn
    !! fcn(n, t, y(1:n), yp(1:n))
    COOP_REAL x0, y(n), h
    COOP_REAL,optional::err
    COOP_REAL,dimension(n)::k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16
    COOP_REAL h6, h12
    h12 = a2 * h
    h6 = a3 * h

    call fcn(n, x0, y, k1)
    call fcn(n, x0+h12, y +  h12 * k1, k2)
    call fcn(n, x0+a3*h, y +  h * ( b31*k1 + b32*k2), k3 )
    call fcn(n, x0+a4*h, y +  h * ( b41*k1 + b43*k3), k4 )
    call fcn(n, x0+a5*h, y +  h * ( b51*k1 + b53*k3 + b54*k4), k5 )
    call fcn(n, x0+a6*h, y +  h * ( b61*k1 + b64*k4 + b65*k5), k6 )
    call fcn(n, x0+a7*h, y +  h * ( b71*k1 + b74*k4 + b75*k5 + b76*k6), k7)
    call fcn(n, x0+a8*h, y +  h * ( b81*k1 + b86*k6 + b87*k7), k8 )
    call fcn(n, x0+a9*h, y +  h * ( b91*k1 + b96*k6 + b97*k7 + b98*k8), k9 )
    call fcn(n, x0+a10*h, y +  h * ( b10_1*k1 + b10_6*k6 + b10_7*k7 + b10_8*k8 &
         + b10_9*k9 ), k10 )
    call fcn(n, x0+a11*h, y +  h * ( b11_1*k1 + b11_7*k7 + b11_8*k8 + b11_9*k9 &
         + b11_10 * k10 ), k11 )
    call fcn(n, x0+a12*h, y +  h * ( b12_1*k1 + b12_6*k6 + b12_7*k7 + b12_8*k8 &
         + b12_9*k9 + b12_10 * k10 + b12_11 * k11 ), k12 )
    call fcn(n, x0+a13*h, y +  h * ( b13_1*k1 + b13_6*k6 + b13_7*k7 + b13_8*k8 &
         + b13_9*k9 + b13_10*k10 + b13_11*k11 + b13_12*k12 ), k13 )
    call fcn(n, x0+h, y +  h * ( b14_1*k1 + b14_6*k6 + b14_7*k7 + b14_8*k8 &
         + b14_9*k9 + b14_10*k10 + b14_11*k11 + b14_12*k12 + b14_13*k13 ), k14 )
   call fcn(n, x0+a15*h, y +  h * ( b15_1*k1 + b15_6*k6 + b15_7*k7 + b15_8*k8 &
        + b15_9*k9 + b15_10*k10 + b15_11*k11 + b15_12*k12 + b15_13*k13 ), k15 )
   x0 = x0 + h
   call fcn(n, x0, y +  h * ( b16_1*k1 + b16_6*k6 + b16_7*k7 + b16_8*k8 &
        + b16_9*k9 + b16_10*k10 + b16_11*k11 + b16_12*k12 + b16_13*k13 &
        + b16_15*k15), k16 )
   y = y +  h * ( c1 * k1 + c8 * k8 + c9 * k9 + c10 * k10 + c11 * k11 &
        + c12 * k12 + c13 * k13 + c14 * k14 )
   if(present(err)) &
        err = maxval(abs(e1*k1 + e8*k8 + e9*k9 + e10*k10 + e11*k11 + e12*k12 + e13*k13  &
        + e14*k14 + e15*k15 + e16*k16 ))
 end subroutine coop_RungeKutta8th



!C===============================================================================
  subroutine coop_Dverk(n, fcn, x, y, xend, tol, ind, c, nw, w)
    implicit COOP_REAL (a-h,o-z)
    implicit COOP_INT (i-n)
    dimension y(n),c(24), w(nw,9)
    external fcn
!c
!c     ******************************************************************
!c     * begin initialization, parameter checking, interrupt re-entries *
!c     ******************************************************************
!c
!c  ......abort if ind out of range 1 to 6
    if (ind.lt.1 .or. ind.gt.6) goto 500
!c
!c        cases - initial entry, normal re-entry, interrupt re-entries
    goto (5, 5, 45, 1111, 2222, 2222), ind
!c        case 1 - initial entry (ind .eq. 1 or 2)
!c  .........abort if n.gt.nw or tol.le.0
5   if (n.gt.nw .or. tol.le.0.d0) go to 500
    if (ind.eq. 2) go to 15
!c              initial entry without options (ind .eq. 1)
!c              set c(1) to c(9) equal to 0
    do k = 1, 9
       c(k) = 0.d0
    enddo
10 continue
    go to 35
15  continue
!c              initial entry with options (ind .eq. 2)
!c              make c(1) to c(9) non-negative
    do  k = 1, 9
       c(k) = dabs(c(k))
    enddo
20  continue
!c              make floor values non-negative if they are to be used
    if (c(1).ne.4.d0 .and. c(1).ne.5.d0) go to 30
    do  k = 1, n
       c(k+30) = dabs(c(k+30))
    enddo
25  continue
30  continue
35  continue
!c           initialize rreb, dwarf, prev xend, flag, counts
    c(10) = 2.d0**(-56)
    c(11) = 1.d-35
!c           set previous xend initially to initial value of x
    c(20) = x
    do  k = 21, 24
       c(k) = 0.d0
    enddo
40  continue
    go to 50
!c        case 2 - normal re-entry (ind .eq. 3)
!c  .........abort if xend reached, and either x changed or xend not
45  if (c(21).ne.0.d0 .and. & 
         (x.ne.c(20) .or. xend.eq.c(20))) go to 500
!c           re-initialize flag
    c(21) = 0.d0
    go to 50
!c        case 3 - re-entry following an interrupt (ind .eq. 4 to 6)
!c           transfer control to the appropriate re-entry point..........
!c           this has already been handled by the computed go to        .
!c        end cases                                                     v
50  continue
!c
!c     end initialization, etc.
!c
!c     ******************************************************************
!c     * loop through the following 4 stages, once for each trial  step *
!c     * until the occurrence of one of the following                   *
!c     *    (a) the normal return (with ind .eq. 3) on reaching xend in *
!c     *        stage 4                                                 *
!c     *    (b) an error return (with ind .lt. 0) in stage 1 or stage 4 *
!c     *    (c) an interrupt return (with ind  .eq.  4,  5  or  6),  if *
!c     *        requested, in stage 1 or stage 4                        *
!c     ******************************************************************
!c
99999 continue
!c
!c        ***************************************************************
!c        * stage 1 - prepare - do calculations of  hmin,  hmax,  etc., *
!c        * and some parameter  checking,  and  end  up  with  suitable *
!c        * values of hmag, xtrial and htrial in preparation for taking *
!c        * an integration step.                                        *
!c        ***************************************************************
!c
!c***********error return (with ind=-1) if no of fcn evals too great
    if (c(7).eq.0.d0 .or. c(24).lt.c(7)) go to 100
    ind = -1
    return
100 continue
!c
!c           calculate slope (adding 1 to no of fcn evals) if ind .ne. 6
    if (ind .eq. 6) go to 105
    call fcn(n, x, y, w(1,1))
    c(24) = c(24) + 1.d0
105 continue
!c
!c           calculate hmin - use default unless value prescribed
    c(13) = c(3)
    if (c(3) .ne. 0.d0) go to 165
!c              calculate default value of hmin
!c              first calculate weighted norm y - c(12) - as specified
!c              by the error control indicator c(1)
    temp = 0.d0
    if (c(1) .ne. 1.d0) go to 115
!c                 absolute error control - weights are 1
    do k = 1, n
       temp = dmax1(temp, dabs(y(k)))
    enddo
110 continue
    c(12) = temp
    go to 160
115 if (c(1) .ne. 2.d0) go to 120
!c                 relative error control - weights are 1/dabs(y(k)) so
!c                 weighted norm y is 1
    c(12) = 1.d0
    go to 160
120 if (c(1) .ne. 3.d0) go to 130
!c                 weights are 1/max(c(2),abs(y(k)))
    do  k = 1, n
       temp = dmax1(temp, dabs(y(k))/c(2))
    enddo
125 continue
    c(12) = dmin1(temp, 1.d0)
    go to 160
130 if (c(1) .ne. 4.d0) go to 140
!c                 weights are 1/max(c(k+30),abs(y(k)))
    do k = 1, n
       temp = dmax1(temp, dabs(y(k))/c(k+30))
    enddo
135 continue
    c(12) = dmin1(temp, 1.d0)
    go to 160
140 if (c(1) .ne. 5.d0) go to 150
!c                 weights are 1/c(k+30)
    do k = 1, n
       temp = dmax1(temp, dabs(y(k))/c(k+30))
    enddo
145 continue
    c(12) = temp
    go to 160
150 continue
!c                 default case - weights are 1/max(1,abs(y(k)))
    do  k = 1, n
       temp = dmax1(temp, dabs(y(k)))
    enddo
155 continue
    c(12) = dmin1(temp, 1.d0)
160 continue
    c(13) = 10.d0*dmax1(c(11),c(10)*dmax1(c(12)/tol,dabs(x)))
165 continue
 !c
!c           calculate scale - use default unless value prescribed
    c(15) = c(5)
    if (c(5) .eq. 0.d0) c(15) = 1.d0
!c
!c           calculate hmax - consider 4 cases
!c           case 1 both hmax and scale prescribed
    if (c(6).ne.0.d0 .and. c(5).ne.0.d0) & 
         c(16) = dmin1(c(6), 2.d0/c(5))
!c           case 2 - hmax prescribed, but scale not
    if (c(6).ne.0.d0 .and. c(5).eq.0.d0) c(16) = c(6)
!c           case 3 - hmax not prescribed, but scale is
    if (c(6).eq.0.d0 .and. c(5).ne.0.d0) c(16) = 2.d0/c(5)
!c           case 4 - neither hmax nor scale is provided
    if (c(6).eq.0.d0 .and. c(5).eq.0.d0) c(16) = 2.d0
!c
!c***********error return (with ind=-2) if hmin .gt. hmax
    if (c(13) .le. c(16)) go to 170
    ind = -2
    return
170 continue
!c
!c           calculate preliminary hmag - consider 3 cases
    if (ind .gt. 2) go to 175
!c           case 1 - initial entry - use prescribed value of hstart, if
!c              any, else default
    c(14) = c(4)
    if (c(4) .eq. 0.d0) c(14) = c(16)*tol**(1./6.)
    go to 185
175 if (c(23) .gt. 1.d0) go to 180
!c           case 2 - after a successful step, or at most  one  failure,
!c              use min(2, .9*(tol/est)**(1/6))*hmag, but avoid possible
!c              overflow. then avoid reduction by more than half.
    temp = 2.d0*c(14)
    if (tol .lt. (2.d0/.9d0)**6*c(19)) & 
         temp = .9d0*(tol/c(19))**(1./6.)*c(14)
    c(14) = dmax1(temp, .5d0*c(14))
    go to 185
180 continue
!c           case 3 - after two or more successive failures
    c(14) = .5d0*c(14)
185 continue
!c
!c           check against hmax
    c(14) = dmin1(c(14), c(16))
!c
!c           check against hmin
    c(14) = dmax1(c(14), c(13))
!c
!c***********interrupt no 1 (with ind=4) if requested
    if (c(8) .eq. 0.d0) go to 1111
    ind = 4
    return
!c           resume here on re-entry with ind .eq. 4   ........re-entry..
1111 continue
!c
!c           calculate hmag, xtrial - depending on preliminary hmag, xend
    if (c(14) .ge. dabs(xend - x)) go to 190
!c              do not step more than half way to xend
    c(14) = dmin1(c(14), .5d0*dabs(xend - x))
    c(17) = x + dsign(c(14), xend - x)
    go to 195
190 continue
!c              hit xend exactly
    c(14) = dabs(xend - x)
    c(17) = xend
195 continue
!c
!c           calculate htrial
    c(18) = c(17) - x
!c
!c        end stage 1
!c
!c        ***************************************************************
!c        * stage 2 - calculate ytrial (adding 7 to no of  fcn  evals). *
!c        * w(*,2), ... w(*,8)  hold  intermediate  results  needed  in *
!c        * stage 3. w(*,9) is temporary storage until finally it holds *
!c        * ytrial.                                                     *
!c        ***************************************************************
!c
    temp = c(18)/1398169080000.d0
!c
    do k = 1, n
       w(k,9) = y(k) + temp*w(k,1)*233028180000.d0
    enddo
200 continue
    call fcn(n, x + c(18)/6.d0, w(1,9), w(1,2))
!c
    do k = 1, n
       w(k,9) = y(k) + temp*(   w(k,1)*74569017600.d0 & 
            + w(k,2)*298276070400.d0  )
    enddo
205 continue
    call fcn(n, x + c(18)*(4.d0/15.d0), w(1,9), w(1,3))
!c
    do  k = 1, n
       w(k,9) = y(k) + temp*(   w(k,1)*1165140900000.d0 & 
            - w(k,2)*3728450880000.d0 & 
            + w(k,3)*3495422700000.d0 )
    enddo
210 continue
    call fcn(n, x + c(18)*(2.d0/3.d0), w(1,9), w(1,4))
!c
    do k = 1, n
       w(k,9) = y(k) + temp*( - w(k,1)*3604654659375.d0 & 
            + w(k,2)*12816549900000.d0 & 
            - w(k,3)*9284716546875.d0 & 
            + w(k,4)*1237962206250.d0 )
    enddo
    continue
    call fcn(n, x + c(18)*(5.d0/6.d0), w(1,9), w(1,5))
!c
    do k = 1, n
       w(k,9) = y(k) + temp*(   w(k,1)*3355605792000.d0 & 
            - w(k,2)*11185352640000.d0 & 
            + w(k,3)*9172628850000.d0 & 
            - w(k,4)*427218330000.d0 & 
            + w(k,5)*482505408000.d0  )
    enddo
220 continue
    call fcn(n, x + c(18), w(1,9), w(1,6))
!c
    do k = 1, n
       w(k,9) = y(k) + temp*( - w(k,1)*770204740536.d0 & 
            + w(k,2)*2311639545600.d0 & 
            - w(k,3)*1322092233000.d0 & 
            - w(k,4)*453006781920.d0 & 
            + w(k,5)*326875481856.d0  )
    enddo
225    continue
    call fcn(n, x + c(18)/15.d0, w(1,9), w(1,7))
!c
    do k = 1, n
       w(k,9) = y(k) + temp*(   w(k,1)*2845924389000.d0 & 
            - w(k,2)*9754668000000.d0 & 
            + w(k,3)*7897110375000.d0 & 
            - w(k,4)*192082660000.d0 & 
            + w(k,5)*400298976000.d0 & 
            + w(k,7)*201586000000.d0  )
    enddo
230 continue
    call fcn(n, x + c(18), w(1,9), w(1,8))
!c
!c           calculate ytrial, the extrapolated approximation and store
!c              in w(*,9)
    do k = 1, n
       w(k,9) = y(k) + temp*(   w(k,1)*104862681000.d0 & 
            + w(k,3)*545186250000.d0 & 
            + w(k,4)*446637345000.d0 & 
            + w(k,5)*188806464000.d0 & 
            + w(k,7)*15076875000.d0 & 
            + w(k,8)*97599465000.d0   )
    enddo
235    continue
!c
!c           add 7 to the no of fcn evals
    c(24) = c(24) + 7.d0
!c
!c        end stage 2
!c
!c        ***************************************************************
!c        * stage 3 - calculate the error estimate est. first calculate *
!c        * the  unweighted  absolute  error  estimate vector (per unit *
!c        * step) for the unextrapolated approximation and store it  in *
!c        * w(*,2).  then  calculate the weighted max norm of w(*,2) as *
!c        * specified by the error  control  indicator  c(1).  finally, *
!c        * modify  this result to produce est, the error estimate (per *
!c        * unit step) for the extrapolated approximation ytrial.       *
!c        ***************************************************************
!c
!c           calculate the unweighted absolute error estimate vector
    do k = 1, n
       w(k,2) = (   w(k,1)*8738556750.d0 & 
            + w(k,3)*9735468750.d0 & 
            - w(k,4)*9709507500.d0 & 
            + w(k,5)*8582112000.d0 & 
            + w(k,6)*95329710000.d0 & 
            - w(k,7)*15076875000.d0 & 
            - w(k,8)*97599465000.d0)/1398169080000.d0
    enddo
300 continue
!c
!c           calculate the weighted max norm of w(*,2) as specified by
!c           the error control indicator c(1)
    temp = 0.d0
    if (c(1) .ne. 1.d0) go to 310
!c              absolute error control
    do  k = 1, n
       temp = dmax1(temp,dabs(w(k,2)))
    enddo
305 continue
    go to 360
310 if (c(1) .ne. 2.d0) go to 320
!c              relative error control
    do k = 1, n
       temp = dmax1(temp, dabs(w(k,2)/y(k)))
    enddo
315 continue
    go to 360
320 if (c(1) .ne. 3.d0) go to 330
!c              weights are 1/max(c(2),abs(y(k)))
    do k = 1, n
       temp = dmax1(temp, dabs(w(k,2)) & 
            / dmax1(c(2), dabs(y(k))) )
    enddo
325 continue
    go to 360
330 if (c(1) .ne. 4.d0) go to 340
!c              weights are 1/max(c(k+30),abs(y(k)))
    do k = 1, n
       temp = dmax1(temp, dabs(w(k,2)) & 
            / dmax1(c(k+30), dabs(y(k))) )
    enddo
335 continue
    go to 360
340 if (c(1) .ne. 5.d0) go to 350
!c              weights are 1/c(k+30)
    do k = 1, n
       temp = dmax1(temp, dabs(w(k,2)/c(k+30)))
    enddo
345 continue
    go to 360
350 continue
!c              default case - weights are 1/max(1,abs(y(k)))
    do k = 1, n
       temp = dmax1(temp, dabs(w(k,2)) & 
            / dmax1(1.d0, dabs(y(k))) )
    enddo
355 continue
360 continue
!c
!c           calculate est - (the weighted max norm of w(*,2))*hmag*scale
!c              - est is intended to be a measure of the error  per  unit
!c              step in ytrial
    c(19) = temp*c(14)*c(15)
!c
!c        end stage 3
!c
!c        ***************************************************************
!c        * stage 4 - make decisions.                                   *
!c        ***************************************************************
!c
!c           set ind=5 if step acceptable, else set ind=6
    ind = 5
    if (c(19) .gt. tol) ind = 6
!c
!c***********interrupt no 2 if requested
    if (c(9) .eq. 0.d0) go to 2222
    return
!c           resume here on re-entry with ind .eq. 5 or 6   ...re-entry..
2222 continue
!c
    if (ind .eq. 6) go to 410
!c              step accepted (ind .eq. 5), so update x, y from xtrial,
!c                 ytrial, add 1 to the no of successful steps, and set
!c                 the no of successive failures to zero
    x = c(17)
    do k = 1, n
       y(k) = w(k,9)
    enddo
400 continue
    c(22) = c(22) + 1.d0
    c(23) = 0.d0
!c**************return(with ind=3, xend saved, flag set) if x .eq. xend
    if (x .ne. xend) go to 405
    ind = 3
    c(20) = xend
    c(21) = 1.d0
    return
405 continue
    go to 420
410 continue
!c              step not accepted (ind .eq. 6), so add 1 to the no of
!c                 successive failures
    c(23) = c(23) + 1.d0
!c**************error return (with ind=-3) if hmag .le. hmin
    if (c(14) .gt. c(13)) go to 415
    ind = -3
    return
415 continue
420 continue
!c
!c        end stage 4
!c
    go to 99999
!c     end loop
!c
!c  begin abort action
500 continue
!c
    write(*,*) 'dverk crashes with ind =',  ind
    write(*,*) 'x = ', x
    stop
!c
!c  end abort action
!c
  end subroutine coop_Dverk

!!fcn(n, x, y, yprime, args)
  subroutine coop_dverk_with_arguments(n, fcn, args, x, y, xend, tol, ind, c, nw, w)
    implicit COOP_REAL (a-h,o-z)
    implicit COOP_INT (i-n)
    type(coop_arguments) args
    dimension y(n),c(24), w(nw,9)
    external fcn
!c
!c     ******************************************************************
!c     * begin initialization, parameter checking, interrupt re-entries *
!c     ******************************************************************
!c
!c  ......abort if ind out of range 1 to 6
    if (ind.lt.1 .or. ind.gt.6) goto 500
!c
!c        cases - initial entry, normal re-entry, interrupt re-entries
    goto (5, 5, 45, 1111, 2222, 2222), ind
!c        case 1 - initial entry (ind .eq. 1 or 2)
!c  .........abort if n.gt.nw or tol.le.0
5   if (n.gt.nw .or. tol.le.0.d0) go to 500
    if (ind.eq. 2) go to 15
!c              initial entry without options (ind .eq. 1)
!c              set c(1) to c(9) equal to 0
    do k = 1, 9
       c(k) = 0.d0
    enddo
10 continue
    go to 35
15  continue
!c              initial entry with options (ind .eq. 2)
!c              make c(1) to c(9) non-negative
    do  k = 1, 9
       c(k) = dabs(c(k))
    enddo
20  continue
!c              make floor values non-negative if they are to be used
    if (c(1).ne.4.d0 .and. c(1).ne.5.d0) go to 30
    do  k = 1, n
       c(k+30) = dabs(c(k+30))
    enddo
25  continue
30  continue
35  continue
!c           initialize rreb, dwarf, prev xend, flag, counts
    c(10) = 2.d0**(-56)
    c(11) = 1.d-35
!c           set previous xend initially to initial value of x
    c(20) = x
    do  k = 21, 24
       c(k) = 0.d0
    enddo
40  continue
    go to 50
!c        case 2 - normal re-entry (ind .eq. 3)
!c  .........abort if xend reached, and either x changed or xend not
45  if (c(21).ne.0.d0 .and. & 
         (x.ne.c(20) .or. xend.eq.c(20))) go to 500
!c           re-initialize flag
    c(21) = 0.d0
    go to 50
!c        case 3 - re-entry following an interrupt (ind .eq. 4 to 6)
!c           transfer control to the appropriate re-entry point..........
!c           this has already been handled by the computed go to        .
!c        end cases                                                     v
50  continue
!c
!c     end initialization, etc.
!c
!c     ******************************************************************
!c     * loop through the following 4 stages, once for each trial  step *
!c     * until the occurrence of one of the following                   *
!c     *    (a) the normal return (with ind .eq. 3) on reaching xend in *
!c     *        stage 4                                                 *
!c     *    (b) an error return (with ind .lt. 0) in stage 1 or stage 4 *
!c     *    (c) an interrupt return (with ind  .eq.  4,  5  or  6),  if *
!c     *        requested, in stage 1 or stage 4                        *
!c     ******************************************************************
!c
99999 continue
!c
!c        ***************************************************************
!c        * stage 1 - prepare - do calculations of  hmin,  hmax,  etc., *
!c        * and some parameter  checking,  and  end  up  with  suitable *
!c        * values of hmag, xtrial and htrial in preparation for taking *
!c        * an integration step.                                        *
!c        ***************************************************************
!c
!c***********error return (with ind=-1) if no of fcn evals too great
    if (c(7).eq.0.d0 .or. c(24).lt.c(7)) go to 100
    ind = -1
    return
100 continue
!c
!c           calculate slope (adding 1 to no of fcn evals) if ind .ne. 6
    if (ind .eq. 6) go to 105
    call fcn(n, x, y, w(1,1), args)
    c(24) = c(24) + 1.d0
105 continue
!c
!c           calculate hmin - use default unless value prescribed
    c(13) = c(3)
    if (c(3) .ne. 0.d0) go to 165
!c              calculate default value of hmin
!c              first calculate weighted norm y - c(12) - as specified
!c              by the error control indicator c(1)
    temp = 0.d0
    if (c(1) .ne. 1.d0) go to 115
!c                 absolute error control - weights are 1
    do k = 1, n
       temp = dmax1(temp, dabs(y(k)))
    enddo
110 continue
    c(12) = temp
    go to 160
115 if (c(1) .ne. 2.d0) go to 120
!c                 relative error control - weights are 1/dabs(y(k)) so
!c                 weighted norm y is 1
    c(12) = 1.d0
    go to 160
120 if (c(1) .ne. 3.d0) go to 130
!c                 weights are 1/max(c(2),abs(y(k)))
    do  k = 1, n
       temp = dmax1(temp, dabs(y(k))/c(2))
    enddo
125 continue
    c(12) = dmin1(temp, 1.d0)
    go to 160
130 if (c(1) .ne. 4.d0) go to 140
!c                 weights are 1/max(c(k+30),abs(y(k)))
    do k = 1, n
       temp = dmax1(temp, dabs(y(k))/c(k+30))
    enddo
135 continue
    c(12) = dmin1(temp, 1.d0)
    go to 160
140 if (c(1) .ne. 5.d0) go to 150
!c                 weights are 1/c(k+30)
    do k = 1, n
       temp = dmax1(temp, dabs(y(k))/c(k+30))
    enddo
145 continue
    c(12) = temp
    go to 160
150 continue
!c                 default case - weights are 1/max(1,abs(y(k)))
    do  k = 1, n
       temp = dmax1(temp, dabs(y(k)))
    enddo
155 continue
    c(12) = dmin1(temp, 1.d0)
160 continue
    c(13) = 10.d0*dmax1(c(11),c(10)*dmax1(c(12)/tol,dabs(x)))
165 continue
 !c
!c           calculate scale - use default unless value prescribed
    c(15) = c(5)
    if (c(5) .eq. 0.d0) c(15) = 1.d0
!c
!c           calculate hmax - consider 4 cases
!c           case 1 both hmax and scale prescribed
    if (c(6).ne.0.d0 .and. c(5).ne.0.d0) & 
         c(16) = dmin1(c(6), 2.d0/c(5))
!c           case 2 - hmax prescribed, but scale not
    if (c(6).ne.0.d0 .and. c(5).eq.0.d0) c(16) = c(6)
!c           case 3 - hmax not prescribed, but scale is
    if (c(6).eq.0.d0 .and. c(5).ne.0.d0) c(16) = 2.d0/c(5)
!c           case 4 - neither hmax nor scale is provided
    if (c(6).eq.0.d0 .and. c(5).eq.0.d0) c(16) = 2.d0
!c
!c***********error return (with ind=-2) if hmin .gt. hmax
    if (c(13) .le. c(16)) go to 170
    ind = -2
    return
170 continue
!c
!c           calculate preliminary hmag - consider 3 cases
    if (ind .gt. 2) go to 175
!c           case 1 - initial entry - use prescribed value of hstart, if
!c              any, else default
    c(14) = c(4)
    if (c(4) .eq. 0.d0) c(14) = c(16)*tol**(1./6.)
    go to 185
175 if (c(23) .gt. 1.d0) go to 180
!c           case 2 - after a successful step, or at most  one  failure,
!c              use min(2, .9*(tol/est)**(1/6))*hmag, but avoid possible
!c              overflow. then avoid reduction by more than half.
    temp = 2.d0*c(14)
    if (tol .lt. (2.d0/.9d0)**6*c(19)) & 
         temp = .9d0*(tol/c(19))**(1./6.)*c(14)
    c(14) = dmax1(temp, .5d0*c(14))
    go to 185
180 continue
!c           case 3 - after two or more successive failures
    c(14) = .5d0*c(14)
185 continue
!c
!c           check against hmax
    c(14) = dmin1(c(14), c(16))
!c
!c           check against hmin
    c(14) = dmax1(c(14), c(13))
!c
!c***********interrupt no 1 (with ind=4) if requested
    if (c(8) .eq. 0.d0) go to 1111
    ind = 4
    return
!c           resume here on re-entry with ind .eq. 4   ........re-entry..
1111 continue
!c
!c           calculate hmag, xtrial - depending on preliminary hmag, xend
    if (c(14) .ge. dabs(xend - x)) go to 190
!c              do not step more than half way to xend
    c(14) = dmin1(c(14), .5d0*dabs(xend - x))
    c(17) = x + dsign(c(14), xend - x)
    go to 195
190 continue
!c              hit xend exactly
    c(14) = dabs(xend - x)
    c(17) = xend
195 continue
!c
!c           calculate htrial
    c(18) = c(17) - x
!c
!c        end stage 1
!c
!c        ***************************************************************
!c        * stage 2 - calculate ytrial (adding 7 to no of  fcn  evals). *
!c        * w(*,2), ... w(*,8)  hold  intermediate  results  needed  in *
!c        * stage 3. w(*,9) is temporary storage until finally it holds *
!c        * ytrial.                                                     *
!c        ***************************************************************
!c
    temp = c(18)/1398169080000.d0
!c
    do k = 1, n
       w(k,9) = y(k) + temp*w(k,1)*233028180000.d0
    enddo
200 continue
    call fcn(n, x + c(18)/6.d0, w(1,9), w(1,2), args)
!c
    do k = 1, n
       w(k,9) = y(k) + temp*(   w(k,1)*74569017600.d0 & 
            + w(k,2)*298276070400.d0  )
    enddo
205 continue
    call fcn(n, x + c(18)*(4.d0/15.d0), w(1,9), w(1,3), args)
!c
    do  k = 1, n
       w(k,9) = y(k) + temp*(   w(k,1)*1165140900000.d0 & 
            - w(k,2)*3728450880000.d0 & 
            + w(k,3)*3495422700000.d0 )
    enddo
210 continue
    call fcn(n, x + c(18)*(2.d0/3.d0), w(1,9), w(1,4), args)
!c
    do k = 1, n
       w(k,9) = y(k) + temp*( - w(k,1)*3604654659375.d0 & 
            + w(k,2)*12816549900000.d0 & 
            - w(k,3)*9284716546875.d0 & 
            + w(k,4)*1237962206250.d0 )
    enddo
    continue
    call fcn(n, x + c(18)*(5.d0/6.d0), w(1,9), w(1,5), args)
!c
    do k = 1, n
       w(k,9) = y(k) + temp*(   w(k,1)*3355605792000.d0 & 
            - w(k,2)*11185352640000.d0 & 
            + w(k,3)*9172628850000.d0 & 
            - w(k,4)*427218330000.d0 & 
            + w(k,5)*482505408000.d0  )
    enddo
220 continue
    call fcn(n, x + c(18), w(1,9), w(1,6), args)
!c
    do k = 1, n
       w(k,9) = y(k) + temp*( - w(k,1)*770204740536.d0 & 
            + w(k,2)*2311639545600.d0 & 
            - w(k,3)*1322092233000.d0 & 
            - w(k,4)*453006781920.d0 & 
            + w(k,5)*326875481856.d0  )
    enddo
225    continue
    call fcn(n, x + c(18)/15.d0, w(1,9), w(1,7), args)
!c
    do k = 1, n
       w(k,9) = y(k) + temp*(   w(k,1)*2845924389000.d0 & 
            - w(k,2)*9754668000000.d0 & 
            + w(k,3)*7897110375000.d0 & 
            - w(k,4)*192082660000.d0 & 
            + w(k,5)*400298976000.d0 & 
            + w(k,7)*201586000000.d0  )
    enddo
230 continue
    call fcn(n, x + c(18), w(1,9), w(1,8), args)
!c
!c           calculate ytrial, the extrapolated approximation and store
!c              in w(*,9)
    do k = 1, n
       w(k,9) = y(k) + temp*(   w(k,1)*104862681000.d0 & 
            + w(k,3)*545186250000.d0 & 
            + w(k,4)*446637345000.d0 & 
            + w(k,5)*188806464000.d0 & 
            + w(k,7)*15076875000.d0 & 
            + w(k,8)*97599465000.d0   )
    enddo
235    continue
!c
!c           add 7 to the no of fcn evals
    c(24) = c(24) + 7.d0
!c
!c        end stage 2
!c
!c        ***************************************************************
!c        * stage 3 - calculate the error estimate est. first calculate *
!c        * the  unweighted  absolute  error  estimate vector (per unit *
!c        * step) for the unextrapolated approximation and store it  in *
!c        * w(*,2).  then  calculate the weighted max norm of w(*,2) as *
!c        * specified by the error  control  indicator  c(1).  finally, *
!c        * modify  this result to produce est, the error estimate (per *
!c        * unit step) for the extrapolated approximation ytrial.       *
!c        ***************************************************************
!c
!c           calculate the unweighted absolute error estimate vector
    do k = 1, n
       w(k,2) = (   w(k,1)*8738556750.d0 & 
            + w(k,3)*9735468750.d0 & 
            - w(k,4)*9709507500.d0 & 
            + w(k,5)*8582112000.d0 & 
            + w(k,6)*95329710000.d0 & 
            - w(k,7)*15076875000.d0 & 
            - w(k,8)*97599465000.d0)/1398169080000.d0
    enddo
300 continue
!c
!c           calculate the weighted max norm of w(*,2) as specified by
!c           the error control indicator c(1)
    temp = 0.d0
    if (c(1) .ne. 1.d0) go to 310
!c              absolute error control
    do  k = 1, n
       temp = dmax1(temp,dabs(w(k,2)))
    enddo
305 continue
    go to 360
310 if (c(1) .ne. 2.d0) go to 320
!c              relative error control
    do k = 1, n
       temp = dmax1(temp, dabs(w(k,2)/y(k)))
    enddo
315 continue
    go to 360
320 if (c(1) .ne. 3.d0) go to 330
!c              weights are 1/max(c(2),abs(y(k)))
    do k = 1, n
       temp = dmax1(temp, dabs(w(k,2)) & 
            / dmax1(c(2), dabs(y(k))) )
    enddo
325 continue
    go to 360
330 if (c(1) .ne. 4.d0) go to 340
!c              weights are 1/max(c(k+30),abs(y(k)))
    do k = 1, n
       temp = dmax1(temp, dabs(w(k,2)) & 
            / dmax1(c(k+30), dabs(y(k))) )
    enddo
335 continue
    go to 360
340 if (c(1) .ne. 5.d0) go to 350
!c              weights are 1/c(k+30)
    do k = 1, n
       temp = dmax1(temp, dabs(w(k,2)/c(k+30)))
    enddo
345 continue
    go to 360
350 continue
!c              default case - weights are 1/max(1,abs(y(k)))
    do k = 1, n
       temp = dmax1(temp, dabs(w(k,2)) & 
            / dmax1(1.d0, dabs(y(k))) )
    enddo
355 continue
360 continue
!c
!c           calculate est - (the weighted max norm of w(*,2))*hmag*scale
!c              - est is intended to be a measure of the error  per  unit
!c              step in ytrial
    c(19) = temp*c(14)*c(15)
!c
!c        end stage 3
!c
!c        ***************************************************************
!c        * stage 4 - make decisions.                                   *
!c        ***************************************************************
!c
!c           set ind=5 if step acceptable, else set ind=6
    ind = 5
    if (c(19) .gt. tol) ind = 6
!c
!c***********interrupt no 2 if requested
    if (c(9) .eq. 0.d0) go to 2222
    return
!c           resume here on re-entry with ind .eq. 5 or 6   ...re-entry..
2222 continue
!c
    if (ind .eq. 6) go to 410
!c              step accepted (ind .eq. 5), so update x, y from xtrial,
!c                 ytrial, add 1 to the no of successful steps, and set
!c                 the no of successive failures to zero
    x = c(17)
    do k = 1, n
       y(k) = w(k,9)
    enddo
400 continue
    c(22) = c(22) + 1.d0
    c(23) = 0.d0
!c**************return(with ind=3, xend saved, flag set) if x .eq. xend
    if (x .ne. xend) go to 405
    ind = 3
    c(20) = xend
    c(21) = 1.d0
    return
405 continue
    go to 420
410 continue
!c              step not accepted (ind .eq. 6), so add 1 to the no of
!c                 successive failures
    c(23) = c(23) + 1.d0
!c**************error return (with ind=-3) if hmag .le. hmin
    if (c(14) .gt. c(13)) go to 415
    ind = -3
    return
415 continue
420 continue
!c
!c        end stage 4
!c
    go to 99999
!c     end loop
!c
!c  begin abort action
500 continue
!c
    write(*,*) 'dverk crashes with ind =',  ind
    write(*,*) 'x = ', x
    stop
!c
!c  end abort action
!c
  end subroutine coop_dverk_with_arguments






  !!sometimes you want to test dverk
  function coop_RungeKutta4th_test(n, fcn, t, y, h, illegal, update) result(good)
    external fcn, illegal
    logical illegal
    logical good
    COOP_INT  n
    COOP_REAL y(n)
    COOP_REAL t, h
    COOP_REAL k1(n), k2(n), k3(n), k4(n)
    COOP_REAL ytmp(n)
    logical,optional::update
    call fcn(n, t, y, k1)
    ytmp = y+k1* (h/2._dl)
    if(illegal(ytmp) .or. coop_isnan(ytmp))then
       good = .false.
       return
    endif
    call fcn(n, t+h/2._dl, ytmp, k2)
    ytmp = y+k2*(h/2._dl)
    if(illegal(ytmp) .or. coop_isnan(ytmp))then
       good = .false.
       return
    endif
    call fcn(n, t+h/2._dl, ytmp, k3)
    ytmp = y+k3*h
    if(illegal(ytmp) .or. coop_isnan(ytmp))then
       good = .false.
       return
    endif
    call fcn(n, t+h, ytmp, k4)
    ytmp = y+(k1+2.*(k2+k3)+k4)*(h/6._dl)
    good = .not. (illegal(ytmp) .or. coop_isnan(ytmp))
    if(present(update))then
       if(good .and. update)then
          y = ytmp
          t = t+h
       endif
    endif
    return
  end function coop_RungeKutta4th_test


  subroutine coop_RungeKutta4th_with_arguments(n, fcn, args, t, y, h)
    type(coop_arguments) args
    external fcn
    COOP_INT  n
    COOP_REAL y(n)
    COOP_REAL t, h
    COOP_REAL k1(n), k2(n), k3(n), k4(n)
    call fcn(n, t, y, k1, args)
    call fcn(n, t+h/2._dl, y+k1*(h/2._dl), k2, args)
    call fcn(n, t+h/2._dl, y+k2*(h/2._dl), k3, args)
    t = t + h
    call fcn(n, t, y+k3*h, k4, args)
    y = y+(k1+2.*(k2+k3)+k4)*(h/6._dl)
  end subroutine coop_RungeKutta4th_with_arguments



  !!fixed step Runge-Kutta ODE integrator
  subroutine coop_RungeKutta6th_with_arguments(n, fcn, args, t, y, h)
    type(coop_arguments)args
    external fcn
    COOP_REAL,parameter::sqrt21 = 4.5825756949558400066_dl
    COOP_INT  n
    COOP_REAL y(n)
    COOP_REAL t, h
    COOP_REAL k1(n), k2(n), k3(n), k4(n), k5(n), k6(n), k7(n)
    call fcn(n, t, y, k1, args)
    call fcn(n, t+h, y+k1*h, k2, args)
    call fcn(n, t+h/2._dl, y+(3._dl*k1+k2)*(h/8._dl), k3, args)
    call fcn(n, t+(2._dl/3._dl)*h, y+ (4._dl*(k1+k3)+k2)*((2._dl/27._dl)*h), k4, args)
    call fcn(n, t+((7._dl-sqrt21)/14._dl)*h, y+((3._dl*(3._dl*sqrt21 - 7._dl)/392._dl)*k1 - (8._dl*(7._dl-sqrt21)/392._dl)*k2 + (48._dl*(7._dl-sqrt21)/392._dl)*k3 - (3._dl*(21._dl - sqrt21)/392._dl)*k4)*h, k5, args)
    call fcn(n, t+((7._dl+sqrt21)/14._dl)*h, y+ ((-5._dl*(231._dl+51._dl*sqrt21)/1960._dl)*k1 - (40._dl*(7._dl+sqrt21)/1960._dl)*k2 - (320._dl*sqrt21/1960._dl)*k3 + (3._dl*(21._dl + 121._dl*sqrt21)/1960._dl)*k4 + (392._dl*(6._dl+sqrt21)/1960._dl)*k5)*h, k6, args)
    t = t+h
    call fcn(n, t, y+(((22._dl+7._dl*sqrt21)/12._dl)*k1 +(2._dl/3._dl)*k2 + (2._dl*(7._dl*sqrt21-5._dl)/9._dl)*k3 - (7._dl*(3._dl*sqrt21-2._dl)/20._dl)*k4 - (7._dl*(49._dl+9._dl*sqrt21)/90._dl)*k5 + (7._dl*(7._dl-sqrt21)/18._dl)*k6)*h, k7, args)
    y = y+  ((k1+k7)/20._dl + (49._dl/180._dl)*(k5+k6) + (16._dl/45._dl)*k3)*h
  end subroutine coop_RungeKutta6th_with_arguments
  
  subroutine coop_RungeKutta8th_with_arguments (n, fcn, args, x0, y, h, err)
    type(coop_arguments)args
    COOP_REAL,parameter:: c1 = 103.d0 / 1680.d0
    COOP_REAL,parameter:: c8 = -27.d0 / 140.d0
    COOP_REAL,parameter:: c9 = 76.d0 / 105.d0
    COOP_REAL,parameter:: c10 = -201.d0 / 280.d0
    COOP_REAL,parameter:: c11 = 1024.d0 / 1365.d0
    COOP_REAL,parameter:: c12 = 3.d0 / 7280.d0
    COOP_REAL,parameter:: c13 = 12.d0 / 35.d0
    COOP_REAL,parameter:: c14 = 9.d0 / 280.d0
    COOP_REAL,parameter:: a2 = 1.d0 / 12.d0
    COOP_REAL,parameter:: a3 = 1.d0 / 9.d0
    COOP_REAL,parameter:: a4 = 1.d0 / 6.d0
    COOP_REAL,parameter:: a5 = 2.d0 * (1.d0 + coop_sqrt6) / 15.d0
    COOP_REAL,parameter:: a6 = (6.d0 + coop_sqrt6) / 15.d0
    COOP_REAL,parameter:: a7 = (6.d0 - coop_sqrt6) / 15.d0
    COOP_REAL,parameter:: a8 = 2.d0 / 3.d0
    COOP_REAL,parameter:: a9 = 1.d0 / 2.d0
    COOP_REAL,parameter:: a10 = 1.d0 / 3.d0
    COOP_REAL,parameter:: a11 = 1.d0 / 4.d0
    COOP_REAL,parameter:: a12 = 4.d0 / 3.d0
    COOP_REAL,parameter:: a13 = 5.d0 / 6.d0
    COOP_REAL,parameter:: a15 = 1.d0 / 6.d0
    COOP_REAL,parameter:: b31 = 1.d0 / 27.d0
    COOP_REAL,parameter:: b32 = 2.d0 / 27.d0
    COOP_REAL,parameter:: b41 = 1.d0 / 24.d0
    COOP_REAL,parameter:: b43 = 3.d0 / 24.d0
    COOP_REAL,parameter:: b51 = (4.d0 + 94.d0 * coop_sqrt6) / 375.d0
    COOP_REAL,parameter:: b53 = -(282.d0 + 252.d0 * coop_sqrt6) / 375.d0
    COOP_REAL,parameter:: b54 = (328.d0 + 208.d0 * coop_sqrt6) / 375.d0
    COOP_REAL,parameter:: b61 = (9.d0 - coop_sqrt6) / 150.d0
    COOP_REAL,parameter:: b64 = (312.d0 + 32.d0 * coop_sqrt6) / 1425.d0
    COOP_REAL,parameter:: b65 = (69.d0 + 29.d0 * coop_sqrt6) / 570.d0
    COOP_REAL,parameter:: b71 = (927.d0 - 347.d0 * coop_sqrt6) / 1250.d0
    COOP_REAL,parameter:: b74 = (-16248.d0 + 7328.d0 * coop_sqrt6) / 9375.d0
    COOP_REAL,parameter:: b75 = (-489.d0 + 179.d0 * coop_sqrt6) / 3750.d0
    COOP_REAL,parameter:: b76 = (14268.d0 - 5798.d0 * coop_sqrt6) / 9375.d0
    COOP_REAL,parameter:: b81 = 4.d0 / 54.d0
    COOP_REAL,parameter:: b86 = (16.d0 - coop_sqrt6) / 54.d0
    COOP_REAL,parameter:: b87 = (16.d0 + coop_sqrt6) / 54.d0
    COOP_REAL,parameter:: b91 = 38.d0 / 512.d0
    COOP_REAL,parameter:: b96 = (118.d0 - 23.d0 * coop_sqrt6) / 512.d0
    COOP_REAL,parameter:: b97 = (118.d0 + 23.d0 * coop_sqrt6) / 512.d0
    COOP_REAL,parameter:: b98 = - 18.d0 / 512.d0
    COOP_REAL,parameter:: b10_1 = 11.d0 / 144.d0
    COOP_REAL,parameter:: b10_6 = (266.d0 - coop_sqrt6) / 864.d0
    COOP_REAL,parameter:: b10_7 = (266.d0 + coop_sqrt6) / 864.d0
    COOP_REAL,parameter:: b10_8 = - 1.d0 / 16.d0
    COOP_REAL,parameter:: b10_9 = - 8.d0 / 27.d0
    COOP_REAL,parameter:: b11_1 = (5034.d0 - 271.d0 * coop_sqrt6) / 61440.d0
    COOP_REAL,parameter:: b11_7 = (7859.d0 - 1626.d0 * coop_sqrt6) / 10240.d0
    COOP_REAL,parameter:: b11_8 = (-2232.d0 + 813.d0 * coop_sqrt6) / 20480.d0
    COOP_REAL,parameter:: b11_9 = (-594.d0  + 271.d0 * coop_sqrt6) / 960.d0
    COOP_REAL,parameter:: b11_10 = (657.d0 - 813.d0 * coop_sqrt6) / 5120.d0
    COOP_REAL,parameter:: b12_1 = (5996.d0 - 3794.d0 * coop_sqrt6) / 405.d0
    COOP_REAL,parameter:: b12_6 = (-4342.d0 - 338.d0 * coop_sqrt6) / 9.d0
    COOP_REAL,parameter:: b12_7 = (154922.d0 - 40458.d0 * coop_sqrt6) / 135.d0
    COOP_REAL,parameter:: b12_8 = (-4176.d0 + 3794.d0 * coop_sqrt6) / 45.d0
    COOP_REAL,parameter:: b12_9 = (-340864.d0 + 242816.d0 * coop_sqrt6) / 405.d0
    COOP_REAL,parameter:: b12_10 = (26304.d0 - 15176.d0 * coop_sqrt6) / 45.d0
    COOP_REAL,parameter:: b12_11 = -26624.d0 / 81.d0
    COOP_REAL,parameter:: b13_1 = (3793.d0 + 2168.d0 * coop_sqrt6) / 103680.d0
    COOP_REAL,parameter:: b13_6 = (4042.d0 + 2263.d0 * coop_sqrt6) / 13824.d0
    COOP_REAL,parameter:: b13_7 = (-231278.d0 + 40717.d0 * coop_sqrt6) / 69120.d0
    COOP_REAL,parameter:: b13_8 = (7947.d0 - 2168.d0 * coop_sqrt6) / 11520.d0
    COOP_REAL,parameter:: b13_9 = (1048.d0 - 542.d0 * coop_sqrt6) / 405.d0
    COOP_REAL,parameter:: b13_10 = (-1383.d0 + 542.d0 * coop_sqrt6) / 720.d0
    COOP_REAL,parameter:: b13_11 = 2624.d0 / 1053.d0
    COOP_REAL,parameter:: b13_12 = 3.d0 / 1664.d0
    COOP_REAL,parameter:: b14_1 = -137.d0 / 1296.d0
    COOP_REAL,parameter:: b14_6 = (5642.d0 - 337.d0 * coop_sqrt6) / 864.d0
    COOP_REAL,parameter:: b14_7 = (5642.d0 + 337.d0 * coop_sqrt6) / 864.d0
    COOP_REAL,parameter:: b14_8 = -299.d0 / 48.d0
    COOP_REAL,parameter:: b14_9 = 184.d0 / 81.d0
    COOP_REAL,parameter:: b14_10 = -44.d0 / 9.d0
    COOP_REAL,parameter:: b14_11 = -5120.d0 / 1053.d0
    COOP_REAL,parameter:: b14_12 = -11.d0 / 468.d0
    COOP_REAL,parameter:: b14_13 = 16.d0 / 9.d0
    COOP_REAL,parameter:: b15_1 = (33617.d0 - 2168.d0 * coop_sqrt6) / 518400.d0
    COOP_REAL,parameter:: b15_6 = (-3846.d0 + 31.d0 * coop_sqrt6) / 13824.d0
    COOP_REAL,parameter:: b15_7 = (155338.d0 - 52807.d0 * coop_sqrt6) / 345600.d0
    COOP_REAL,parameter:: b15_8 = (-12537.d0 + 2168.d0 * coop_sqrt6) / 57600.d0
    COOP_REAL,parameter:: b15_9 = (92.d0 + 542.d0 * coop_sqrt6) / 2025.d0
    COOP_REAL,parameter:: b15_10 = (-1797.d0 - 542.d0 * coop_sqrt6) / 3600.d0
    COOP_REAL,parameter:: b15_11 = 320.d0 / 567.d0
    COOP_REAL,parameter:: b15_12 = -1.d0 / 1920.d0
    COOP_REAL,parameter:: b15_13 = 4.d0 / 105.d0
    COOP_REAL,parameter:: b16_1 = (-36487.d0 - 30352.d0 * coop_sqrt6) / 279600.d0
    COOP_REAL,parameter:: b16_6 = (-29666.d0 - 4499.d0 * coop_sqrt6) / 7456.d0
    COOP_REAL,parameter:: b16_7 = (2779182.d0 - 615973.d0 * coop_sqrt6) / 186400.d0
    COOP_REAL,parameter:: b16_8 = (-94329.d0 + 91056.d0 * coop_sqrt6) / 93200.d0
    COOP_REAL,parameter:: b16_9 = (-232192.d0 + 121408.d0 * coop_sqrt6) / 17475.d0
    COOP_REAL,parameter:: b16_10 = (101226.d0 - 22764.d0 * coop_sqrt6) / 5825.d0
    COOP_REAL,parameter:: b16_11 = - 169984.d0 / 9087.d0
    COOP_REAL,parameter:: b16_12 = - 87.d0 / 30290.d0
    COOP_REAL,parameter:: b16_13 =  492.d0 / 1165.d0
    COOP_REAL,parameter:: b16_15 =  1260.d0 / 233.d0
    COOP_REAL,parameter:: e1 = -1911.d0 / 109200.d0
    COOP_REAL,parameter:: e8 = 34398.d0 / 109200.d0
    COOP_REAL,parameter:: e9 = -61152.d0 / 109200.d0
    COOP_REAL,parameter:: e10 = 114660.d0 / 109200.d0
    COOP_REAL,parameter:: e11 = -114688.d0 / 109200.d0
    COOP_REAL,parameter:: e12 = -63.d0 / 109200.d0
    COOP_REAL,parameter:: e13 = -13104.d0 / 109200.d0
    COOP_REAL,parameter:: e14 = -3510.d0 / 109200.d0
    COOP_REAL,parameter:: e15 = 39312.d0 / 109200.d0
    COOP_REAL,parameter:: e16 = 6058.d0 / 109200.d0

    COOP_INT  n
    external fcn
    COOP_REAL x0, y(n), h
    COOP_REAL,optional::err
    COOP_REAL,dimension(n)::k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16
    COOP_REAL h6, h12
    h12 = a2 * h
    h6 = a3 * h

    call fcn(n, x0, y, k1, args)
    call fcn(n, x0+h12, y +  h12 * k1, k2, args)
    call fcn( n, x0+a3*h, y +  h * ( b31*k1 + b32*k2), k3 , args)
    call fcn( n, x0+a4*h, y +  h * ( b41*k1 + b43*k3), k4 , args)
    call fcn(n, x0+a5*h, y +  h * ( b51*k1 + b53*k3 + b54*k4), k5 , args)
    call fcn( n, x0+a6*h, y +  h * ( b61*k1 + b64*k4 + b65*k5), k6 , args)
    call fcn(n, x0+a7*h, y +  h * ( b71*k1 + b74*k4 + b75*k5 + b76*k6), k7, args)
    call fcn(n, x0+a8*h, y +  h * ( b81*k1 + b86*k6 + b87*k7), k8 , args)
    call fcn(n, x0+a9*h, y +  h * ( b91*k1 + b96*k6 + b97*k7 + b98*k8), k9 , args)
    call fcn(n, x0+a10*h, y +  h * ( b10_1*k1 + b10_6*k6 + b10_7*k7 + b10_8*k8 &
         + b10_9*k9 ), k10 , args)
    call fcn( n, x0+a11*h, y +  h * ( b11_1*k1 + b11_7*k7 + b11_8*k8 + b11_9*k9 &
         + b11_10 * k10 ), k11 , args)
    call fcn( n, x0+a12*h, y +  h * ( b12_1*k1 + b12_6*k6 + b12_7*k7 + b12_8*k8 &
         + b12_9*k9 + b12_10 * k10 + b12_11 * k11 ), k12 , args)
    call fcn(n, x0+a13*h, y +  h * ( b13_1*k1 + b13_6*k6 + b13_7*k7 + b13_8*k8 &
         + b13_9*k9 + b13_10*k10 + b13_11*k11 + b13_12*k12 ), k13 , args)
    call fcn( n, x0+h, y +  h * ( b14_1*k1 + b14_6*k6 + b14_7*k7 + b14_8*k8 &
         + b14_9*k9 + b14_10*k10 + b14_11*k11 + b14_12*k12 + b14_13*k13 ), k14 , args)
    call fcn(n, x0+a15*h, y +  h * ( b15_1*k1 + b15_6*k6 + b15_7*k7 + b15_8*k8 &
         + b15_9*k9 + b15_10*k10 + b15_11*k11 + b15_12*k12 + b15_13*k13 ), k15 , args)
   x0 = x0 + h
   call fcn(n, x0, y +  h * ( b16_1*k1 + b16_6*k6 + b16_7*k7 + b16_8*k8 &
        + b16_9*k9 + b16_10*k10 + b16_11*k11 + b16_12*k12 + b16_13*k13 &
        + b16_15*k15), k16 , args)
   y = y +  h * ( c1 * k1 + c8 * k8 + c9 * k9 + c10 * k10 + c11 * k11 &
        + c12 * k12 + c13 * k13 + c14 * k14 )
   if(present(err)) &
        err = maxval(abs(e1*k1 + e8*k8 + e9*k9 + e10*k10 + e11*k11 + e12*k12 + e13*k13  &
        + e14*k14 + e15*k15 + e16*k16 ))
 end subroutine coop_RungeKutta8th_with_arguments


 !!-------------------------- coop_symplectic integrator ------------------- !!
!!coop_symplectic integrator is useful when a Hamiltonian is splitable, i.e., can be written as
!!  H = H_1 + H_2 + ... + H_n
!! where H_i and H_j are commutable for any i and j

  subroutine coop_symplectic_2nd(dt, nsteps, Hamiltonian, numterms)
    !! 2nd-order coop_symplectic integrator
    COOP_REAL dt
    COOP_INT  nsteps, j, numterms
    external Hamiltonian
    !!Hamiltonian(dt, i) evolves the system with operator exp( H_i dt), where H_i is the i-th term of Hamiltonian
    call Hamiltonian(dt/2._dl, 1)
    do j=1,nsteps-1
       call coop_symplectic_o2step(dt,1._dl, 1._dl, Hamiltonian, numterms)
    enddo
    call coop_symplectic_o2step(dt,1._dl, 0._dl, Hamiltonian, numterms)    
  end subroutine coop_symplectic_2nd


  subroutine coop_symplectic_4th(dt,nsteps, Hamiltonian, numterms)
    !! 4th-order coop_symplectic integrator
    COOP_REAL dt
    COOP_INT  nsteps, j, numterms
    external Hamiltonian
    !!Hamiltonian(dt, i) evolves the system with operator exp( H_i dt), where H_i is the i-th term of Hamiltonian
    COOP_REAL,parameter:: c1 = 1._dl/(2._dl - 2._dl**(1._dl/3._dl))
    COOP_REAL,parameter:: c0 = 1._dl - 2._dl*c1
    call Hamiltonian(c1*dt/2._dl,1)
    do j=1,nsteps - 1
       call coop_symplectic_o2step(dt, c1, c0, Hamiltonian, numterms)
       call coop_symplectic_o2step(dt, c0, c1, Hamiltonian, numterms)
       call coop_symplectic_o2step(dt, c1, c1, Hamiltonian, numterms)
    enddo
    call coop_symplectic_o2step(dt, c1, c0, Hamiltonian, numterms)
    call coop_symplectic_o2step(dt, c0, c1, Hamiltonian, numterms)
    call coop_symplectic_o2step(dt, c1, 0._dl, Hamiltonian, numterms)
  end subroutine coop_symplectic_4th

  subroutine coop_symplectic_6th(dt,nsteps, Hamiltonian, numterms)
    COOP_REAL dt
    COOP_INT  nsteps,j, numterms
    external Hamiltonian
    !!Hamiltonian(dt, i) evolves the system with operator exp( H_i dt), where H_i is the i-th term of Hamiltonian
    COOP_REAL,parameter:: c1 = -1.177679984178871007_dl, c2 = 0.235573213359358134_dl, c3 = 0.784513610477557264_dl
    COOP_REAL,parameter:: c0 = 1._dl-2._dl*(c1+c2+c3)
    call Hamiltonian(c3*dt/2._dl,1)
    do j=1,nsteps-1
       call coop_symplectic_o2step(dt, c3, c2, Hamiltonian, numterms)
       call coop_symplectic_o2step(dt, c2, c1, Hamiltonian, numterms)
       call coop_symplectic_o2step(dt, c1, c0, Hamiltonian, numterms)
       call coop_symplectic_o2step(dt, c0, c1, Hamiltonian, numterms)
       call coop_symplectic_o2step(dt, c1, c2, Hamiltonian, numterms)
       call coop_symplectic_o2step(dt, c2, c3, Hamiltonian, numterms)
       call coop_symplectic_o2step(dt, c3, c3, Hamiltonian, numterms)
    enddo
    call coop_symplectic_o2step(dt, c3, c2, Hamiltonian, numterms)
    call coop_symplectic_o2step(dt, c2, c1, Hamiltonian, numterms)
    call coop_symplectic_o2step(dt, c1, c0, Hamiltonian, numterms)
    call coop_symplectic_o2step(dt, c0, c1, Hamiltonian, numterms)
    call coop_symplectic_o2step(dt, c1, c2, Hamiltonian, numterms)
    call coop_symplectic_o2step(dt, c2, c3, Hamiltonian, numterms)
    call coop_symplectic_o2step(dt, c3, 0._dl, Hamiltonian, numterms)
  end subroutine coop_symplectic_6th


  !!auxilliary routine
  subroutine coop_symplectic_o2step(dt, c1, c2, Hamiltonian, numterms)
    COOP_REAL dt, c1, c2
    COOP_INT  numterms, i
    external Hamiltonian
    !!Hamiltonian(dt, i) evolves the system with operator exp( H_i dt), where H_i is the i-th term of Hamiltonian
    do i = 2, numterms-1
       call Hamiltonian(c1*dt/2.d0, i)
    end do
    call Hamiltonian(c1*dt, numterms)
    do i=numterms-1, 2, -1
       call Hamiltonian(c1*dt/2.d0, i)
    enddo
    call Hamiltonian((c1+c2)*dt/2.d0, 1)
  end subroutine coop_symplectic_o2step


!!========================================================================
!!the gauss_legendre integrator only works well when y' = f(y) , i.e. no explicit dependence on t
!!this can always be achieved by taking t as an extra variable and using dt/dt  = 1
!!for portability I still assume fcn has the form fcn(n, t, y(1:n), yp(1:n))
!! but you should not make use of t in the subroutine fcn (i.e. output yp should not depend on t)
!!-------------------------------------------------------------------------
  ! 4th order implicit Gauss-Legendre integrator
  subroutine gauss_legendre_4th(n, fcn, t, y, h)
    COOP_INT , parameter :: s = 2
    ! Butcher tableau for 4th order Gauss-Legendre method
    COOP_REAL, parameter :: a(s,s) = reshape( (/ 0.25d0, 0.25d0 - 0.5d0/coop_sqrt3, 0.25d0 + 0.5d0/coop_sqrt3, 0.25d0 /), (/ s, s /) )
    COOP_REAL, parameter ::   b(s) = (/ 0.5d0, 0.5d0 /)

    COOP_INT  n, i, k
    COOP_REAL y(n), g(n,s), h, t
    external fcn
    do i = 1,s
       call fcn(n, t, y, g(:,i))
    end do
    do k = 2,16
       g = matmul(g,a)
       do i = 1, s
          call fcn(n, t, y + g(:,i)*h, g(:,i))
       end do
    end do
    t = t + h
    y = y + matmul(g,b)*h
  end subroutine gauss_legendre_4th

! 6th order implicit Gauss-Legendre integrator
  subroutine gauss_legendre_6th(n, fcn, t, y, h)
    COOP_INT , parameter :: s = 3
    COOP_REAL, parameter ::   b(s) = (/ 5.d0/18.d0, 4.d0/9.d0, 5.d0/18.d0/)
    ! Butcher tableau for 6th order Gauss-Legendre method
    COOP_REAL, parameter :: a(s,s) =  reshape( (/ &
         5.d0/36.d0, 2.d0/9.d0 - 1.d0/sqrt(15.d0), 5.d0/36.d0 - 0.5d0/sqrt(15.d0), &
         5.d0/36.d0 + sqrt(15.d0)/24.d0, 2.d0/9.d0, 5.d0/36.d0 - sqrt(15.d0)/24.d0, &
         5.d0/36.d0 + 0.5d0/sqrt(15.d0), 2.d0/9.d0 + 1.d0/sqrt(15.d0), 5.d0/36.d0 /), (/ s, s /) )
    COOP_INT  n, i, k
    COOP_REAL y(n), g(n,s), h, t
    do i = 1,s
       call fcn(n, t, y, g(:,i))
    end do
    do k = 2,16
       g = matmul(g,a)
       do i = 1,s
          call fcn(n, t, y + g(:,i)*h, g(:,i))
       end do
    end do
    t = t + h
    y = y + matmul(g,b)*h
  end subroutine gauss_legendre_6th


  ! 8th order implicit Gauss-Legendre integrator
  subroutine gauss_legendre_8th(n, fcn, t, y, h)
    COOP_INT , parameter :: s = 4
    ! Butcher tableau for 8th order Gauss-Legendre method
    COOP_REAL, parameter :: a(s,s) = reshape( (/ &
         0.869637112843634643432659873054998518d-1, -0.266041800849987933133851304769531093d-1, &
         0.126274626894047245150568805746180936d-1, -0.355514968579568315691098184956958860d-2, &
         0.188118117499868071650685545087171160d0,   0.163036288715636535656734012694500148d0,  &
         -0.278804286024708952241511064189974107d-1,  0.673550059453815551539866908570375889d-2, &
         0.167191921974188773171133305525295945d0,   0.353953006033743966537619131807997707d0,  &
         0.163036288715636535656734012694500148d0,  -0.141906949311411429641535704761714564d-1, &
         0.177482572254522611843442956460569292d0,   0.313445114741868346798411144814382203d0,  &
         0.352676757516271864626853155865953406d0,   0.869637112843634643432659873054998518d-1 /) , (/ s, s /) )
    COOP_REAL, parameter ::   b(s) = (/ &
         0.173927422568726928686531974610999704d0,   0.326072577431273071313468025389000296d0,  &
         0.326072577431273071313468025389000296d0,   0.173927422568726928686531974610999704d0  /)
    COOP_INT  n, i, k
    COOP_REAL y(n), g(n,s)
    COOP_REAL t, h
    external fcn
    do i = 1,s
       call fcn(n, t, y, g(:,i))
    end do
    do k = 2,16
       g = matmul(g,a)
       do i = 1, s
          call fcn(n, t, y + g(:,i)*h, g(:,i))
       end do
    end do
    t = t + h
    y = y + matmul(g,b)*h
  end subroutine gauss_legendre_8th




  subroutine gauss_legendre_4th_with_arguments(n, fcn,args, t, y, h)
    type(coop_arguments)::args
    COOP_INT , parameter :: s = 2
    ! Butcher tableau for 4th order Gauss-Legendre method
    COOP_REAL, parameter :: a(s,s) = reshape( (/ 0.25d0, 0.25d0 - 0.5d0/coop_sqrt3, 0.25d0 + 0.5d0/coop_sqrt3, 0.25d0 /), (/ s, s /) )
    COOP_REAL, parameter ::   b(s) = (/ 0.5d0, 0.5d0 /)

    COOP_INT  n, i, k
    COOP_REAL y(n), g(n,s), h, t
    external fcn
    do i = 1,s
       call fcn(n, t, y, g(:,i), args)
    end do
    do k = 2,16
       g = matmul(g,a)
       do i = 1, s
          call fcn(n, t, y + g(:,i)*h, g(:,i), args)
       end do
    end do
    t = t + h
    y = y + matmul(g,b)*h
  end subroutine gauss_legendre_4th_with_arguments

! 6th order implicit Gauss-Legendre integrator
  subroutine gauss_legendre_6th_with_arguments(n, fcn, args, t, y, h)
    type(coop_arguments)::args
    COOP_INT , parameter :: s = 3
    COOP_REAL, parameter ::   b(s) = (/ 5.d0/18.d0, 4.d0/9.d0, 5.d0/18.d0/)
    ! Butcher tableau for 6th order Gauss-Legendre method
    COOP_REAL, parameter :: a(s,s) =  reshape( (/ &
         5.d0/36.d0, 2.d0/9.d0 - 1.d0/sqrt(15.d0), 5.d0/36.d0 - 0.5d0/sqrt(15.d0), &
         5.d0/36.d0 + sqrt(15.d0)/24.d0, 2.d0/9.d0, 5.d0/36.d0 - sqrt(15.d0)/24.d0, &
         5.d0/36.d0 + 0.5d0/sqrt(15.d0), 2.d0/9.d0 + 1.d0/sqrt(15.d0), 5.d0/36.d0 /), (/ s, s /) )
    COOP_INT  n, i, k
    COOP_REAL y(n), g(n,s), h, t
    do i = 1,s
       call fcn(n, t, y, g(:,i), args)
    end do
    do k = 2,16
       g = matmul(g,a)
       do i = 1,s
          call fcn(n, t, y + g(:,i)*h, g(:,i), args)
       end do
    end do
    t = t + h
    y = y + matmul(g,b)*h
  end subroutine gauss_legendre_6th_with_arguments


  ! 8th order implicit Gauss-Legendre integrator
  subroutine gauss_legendre_8th_with_arguments(n, fcn, args, t, y, h)
    type(coop_arguments)::args
    COOP_INT , parameter :: s = 4
    ! Butcher tableau for 8th order Gauss-Legendre method
    COOP_REAL, parameter :: a(s,s) = reshape( (/ &
         0.869637112843634643432659873054998518d-1, -0.266041800849987933133851304769531093d-1, &
         0.126274626894047245150568805746180936d-1, -0.355514968579568315691098184956958860d-2, &
         0.188118117499868071650685545087171160d0,   0.163036288715636535656734012694500148d0,  &
         -0.278804286024708952241511064189974107d-1,  0.673550059453815551539866908570375889d-2, &
         0.167191921974188773171133305525295945d0,   0.353953006033743966537619131807997707d0,  &
         0.163036288715636535656734012694500148d0,  -0.141906949311411429641535704761714564d-1, &
         0.177482572254522611843442956460569292d0,   0.313445114741868346798411144814382203d0,  &
         0.352676757516271864626853155865953406d0,   0.869637112843634643432659873054998518d-1 /) , (/ s, s /) )
    COOP_REAL, parameter ::   b(s) = (/ &
         0.173927422568726928686531974610999704d0,   0.326072577431273071313468025389000296d0,  &
         0.326072577431273071313468025389000296d0,   0.173927422568726928686531974610999704d0  /)
    COOP_INT  n, i, k
    COOP_REAL y(n), g(n,s)
    COOP_REAL t, h
    external fcn
    do i = 1,s
       call fcn(n, t, y, g(:,i), args)
    end do
    do k = 2,16
       g = matmul(g,a)
       do i = 1, s
          call fcn(n, t, y + g(:,i)*h, g(:,i), args)
       end do
    end do
    t = t + h
    y = y + matmul(g,b)*h
  end subroutine gauss_legendre_8th_with_arguments






!!fcn(n, x, y, yprime, args)
  subroutine coop_dverk_with_cosmology(n, fcn, cosmology, args, x, y, xend, tol, ind, c, nw, w)
    implicit COOP_REAL (a-h,o-z)
    implicit COOP_INT (i-n)
    class(coop_cosmology) cosmology
    type(coop_arguments)::args
    dimension y(n),c(24), w(nw,9)
    external fcn
!c
!c     ******************************************************************
!c     * begin initialization, parameter checking, interrupt re-entries *
!c     ******************************************************************
!c
!c  ......abort if ind out of range 1 to 6
    if (ind.lt.1 .or. ind.gt.6) goto 500
!c
!c        cases - initial entry, normal re-entry, interrupt re-entries
    goto (5, 5, 45, 1111, 2222, 2222), ind
!c        case 1 - initial entry (ind .eq. 1 or 2)
!c  .........abort if n.gt.nw or tol.le.0
5   if (n.gt.nw .or. tol.le.0.d0) go to 500
    if (ind.eq. 2) go to 15
!c              initial entry without options (ind .eq. 1)
!c              set c(1) to c(9) equal to 0
    do k = 1, 9
       c(k) = 0.d0
    enddo
10 continue
    go to 35
15  continue
!c              initial entry with options (ind .eq. 2)
!c              make c(1) to c(9) non-negative
    do  k = 1, 9
       c(k) = dabs(c(k))
    enddo
20  continue
!c              make floor values non-negative if they are to be used
    if (c(1).ne.4.d0 .and. c(1).ne.5.d0) go to 30
    do  k = 1, n
       c(k+30) = dabs(c(k+30))
    enddo
25  continue
30  continue
35  continue
!c           initialize rreb, dwarf, prev xend, flag, counts
    c(10) = 2.d0**(-56)
    c(11) = 1.d-35
!c           set previous xend initially to initial value of x
    c(20) = x
    do  k = 21, 24
       c(k) = 0.d0
    enddo
40  continue
    go to 50
!c        case 2 - normal re-entry (ind .eq. 3)
!c  .........abort if xend reached, and either x changed or xend not
45  if (c(21).ne.0.d0 .and. & 
         (x.ne.c(20) .or. xend.eq.c(20))) go to 500
!c           re-initialize flag
    c(21) = 0.d0
    go to 50
!c        case 3 - re-entry following an interrupt (ind .eq. 4 to 6)
!c           transfer control to the appropriate re-entry point..........
!c           this has already been handled by the computed go to        .
!c        end cases                                                     v
50  continue
!c
!c     end initialization, etc.
!c
!c     ******************************************************************
!c     * loop through the following 4 stages, once for each trial  step *
!c     * until the occurrence of one of the following                   *
!c     *    (a) the normal return (with ind .eq. 3) on reaching xend in *
!c     *        stage 4                                                 *
!c     *    (b) an error return (with ind .lt. 0) in stage 1 or stage 4 *
!c     *    (c) an interrupt return (with ind  .eq.  4,  5  or  6),  if *
!c     *        requested, in stage 1 or stage 4                        *
!c     ******************************************************************
!c
99999 continue
!c
!c        ***************************************************************
!c        * stage 1 - prepare - do calculations of  hmin,  hmax,  etc., *
!c        * and some parameter  checking,  and  end  up  with  suitable *
!c        * values of hmag, xtrial and htrial in preparation for taking *
!c        * an integration step.                                        *
!c        ***************************************************************
!c
!c***********error return (with ind=-1) if no of fcn evals too great
    if (c(7).eq.0.d0 .or. c(24).lt.c(7)) go to 100
    ind = -1
    return
100 continue
!c
!c           calculate slope (adding 1 to no of fcn evals) if ind .ne. 6
    if (ind .eq. 6) go to 105
    call fcn(n, x, y, w(1,1), cosmology, args)
    c(24) = c(24) + 1.d0
105 continue
!c
!c           calculate hmin - use default unless value prescribed
    c(13) = c(3)
    if (c(3) .ne. 0.d0) go to 165
!c              calculate default value of hmin
!c              first calculate weighted norm y - c(12) - as specified
!c              by the error control indicator c(1)
    temp = 0.d0
    if (c(1) .ne. 1.d0) go to 115
!c                 absolute error control - weights are 1
    do k = 1, n
       temp = dmax1(temp, dabs(y(k)))
    enddo
110 continue
    c(12) = temp
    go to 160
115 if (c(1) .ne. 2.d0) go to 120
!c                 relative error control - weights are 1/dabs(y(k)) so
!c                 weighted norm y is 1
    c(12) = 1.d0
    go to 160
120 if (c(1) .ne. 3.d0) go to 130
!c                 weights are 1/max(c(2),abs(y(k)))
    do  k = 1, n
       temp = dmax1(temp, dabs(y(k))/c(2))
    enddo
125 continue
    c(12) = dmin1(temp, 1.d0)
    go to 160
130 if (c(1) .ne. 4.d0) go to 140
!c                 weights are 1/max(c(k+30),abs(y(k)))
    do k = 1, n
       temp = dmax1(temp, dabs(y(k))/c(k+30))
    enddo
135 continue
    c(12) = dmin1(temp, 1.d0)
    go to 160
140 if (c(1) .ne. 5.d0) go to 150
!c                 weights are 1/c(k+30)
    do k = 1, n
       temp = dmax1(temp, dabs(y(k))/c(k+30))
    enddo
145 continue
    c(12) = temp
    go to 160
150 continue
!c                 default case - weights are 1/max(1,abs(y(k)))
    do  k = 1, n
       temp = dmax1(temp, dabs(y(k)))
    enddo
155 continue
    c(12) = dmin1(temp, 1.d0)
160 continue
    c(13) = 10.d0*dmax1(c(11),c(10)*dmax1(c(12)/tol,dabs(x)))
165 continue
 !c
!c           calculate scale - use default unless value prescribed
    c(15) = c(5)
    if (c(5) .eq. 0.d0) c(15) = 1.d0
!c
!c           calculate hmax - consider 4 cases
!c           case 1 both hmax and scale prescribed
    if (c(6).ne.0.d0 .and. c(5).ne.0.d0) & 
         c(16) = dmin1(c(6), 2.d0/c(5))
!c           case 2 - hmax prescribed, but scale not
    if (c(6).ne.0.d0 .and. c(5).eq.0.d0) c(16) = c(6)
!c           case 3 - hmax not prescribed, but scale is
    if (c(6).eq.0.d0 .and. c(5).ne.0.d0) c(16) = 2.d0/c(5)
!c           case 4 - neither hmax nor scale is provided
    if (c(6).eq.0.d0 .and. c(5).eq.0.d0) c(16) = 2.d0
!c
!c***********error return (with ind=-2) if hmin .gt. hmax
    if (c(13) .le. c(16)) go to 170
    ind = -2
    return
170 continue
!c
!c           calculate preliminary hmag - consider 3 cases
    if (ind .gt. 2) go to 175
!c           case 1 - initial entry - use prescribed value of hstart, if
!c              any, else default
    c(14) = c(4)
    if (c(4) .eq. 0.d0) c(14) = c(16)*tol**(1./6.)
    go to 185
175 if (c(23) .gt. 1.d0) go to 180
!c           case 2 - after a successful step, or at most  one  failure,
!c              use min(2, .9*(tol/est)**(1/6))*hmag, but avoid possible
!c              overflow. then avoid reduction by more than half.
    temp = 2.d0*c(14)
    if (tol .lt. (2.d0/.9d0)**6*c(19)) & 
         temp = .9d0*(tol/c(19))**(1./6.)*c(14)
    c(14) = dmax1(temp, .5d0*c(14))
    go to 185
180 continue
!c           case 3 - after two or more successive failures
    c(14) = .5d0*c(14)
185 continue
!c
!c           check against hmax
    c(14) = dmin1(c(14), c(16))
!c
!c           check against hmin
    c(14) = dmax1(c(14), c(13))
!c
!c***********interrupt no 1 (with ind=4) if requested
    if (c(8) .eq. 0.d0) go to 1111
    ind = 4
    return
!c           resume here on re-entry with ind .eq. 4   ........re-entry..
1111 continue
!c
!c           calculate hmag, xtrial - depending on preliminary hmag, xend
    if (c(14) .ge. dabs(xend - x)) go to 190
!c              do not step more than half way to xend
    c(14) = dmin1(c(14), .5d0*dabs(xend - x))
    c(17) = x + dsign(c(14), xend - x)
    go to 195
190 continue
!c              hit xend exactly
    c(14) = dabs(xend - x)
    c(17) = xend
195 continue
!c
!c           calculate htrial
    c(18) = c(17) - x
!c
!c        end stage 1
!c
!c        ***************************************************************
!c        * stage 2 - calculate ytrial (adding 7 to no of  fcn  evals). *
!c        * w(*,2), ... w(*,8)  hold  intermediate  results  needed  in *
!c        * stage 3. w(*,9) is temporary storage until finally it holds *
!c        * ytrial.                                                     *
!c        ***************************************************************
!c
    temp = c(18)/1398169080000.d0
!c
    do k = 1, n
       w(k,9) = y(k) + temp*w(k,1)*233028180000.d0
    enddo
200 continue
    call fcn(n, x + c(18)/6.d0, w(1,9), w(1,2), cosmology, args)
!c
    do k = 1, n
       w(k,9) = y(k) + temp*(   w(k,1)*74569017600.d0 & 
            + w(k,2)*298276070400.d0  )
    enddo
205 continue
    call fcn(n, x + c(18)*(4.d0/15.d0), w(1,9), w(1,3), cosmology, args)
!c
    do  k = 1, n
       w(k,9) = y(k) + temp*(   w(k,1)*1165140900000.d0 & 
            - w(k,2)*3728450880000.d0 & 
            + w(k,3)*3495422700000.d0 )
    enddo
210 continue
    call fcn(n, x + c(18)*(2.d0/3.d0), w(1,9), w(1,4), cosmology, args)
!c
    do k = 1, n
       w(k,9) = y(k) + temp*( - w(k,1)*3604654659375.d0 & 
            + w(k,2)*12816549900000.d0 & 
            - w(k,3)*9284716546875.d0 & 
            + w(k,4)*1237962206250.d0 )
    enddo
    continue
    call fcn(n, x + c(18)*(5.d0/6.d0), w(1,9), w(1,5), cosmology, args)
!c
    do k = 1, n
       w(k,9) = y(k) + temp*(   w(k,1)*3355605792000.d0 & 
            - w(k,2)*11185352640000.d0 & 
            + w(k,3)*9172628850000.d0 & 
            - w(k,4)*427218330000.d0 & 
            + w(k,5)*482505408000.d0  )
    enddo
220 continue
    call fcn(n, x + c(18), w(1,9), w(1,6), cosmology, args)
!c
    do k = 1, n
       w(k,9) = y(k) + temp*( - w(k,1)*770204740536.d0 & 
            + w(k,2)*2311639545600.d0 & 
            - w(k,3)*1322092233000.d0 & 
            - w(k,4)*453006781920.d0 & 
            + w(k,5)*326875481856.d0  )
    enddo
225    continue
    call fcn(n, x + c(18)/15.d0, w(1,9), w(1,7), cosmology, args)
!c
    do k = 1, n
       w(k,9) = y(k) + temp*(   w(k,1)*2845924389000.d0 & 
            - w(k,2)*9754668000000.d0 & 
            + w(k,3)*7897110375000.d0 & 
            - w(k,4)*192082660000.d0 & 
            + w(k,5)*400298976000.d0 & 
            + w(k,7)*201586000000.d0  )
    enddo
230 continue
    call fcn(n, x + c(18), w(1,9), w(1,8), cosmology, args)
!c
!c           calculate ytrial, the extrapolated approximation and store
!c              in w(*,9)
    do k = 1, n
       w(k,9) = y(k) + temp*(   w(k,1)*104862681000.d0 & 
            + w(k,3)*545186250000.d0 & 
            + w(k,4)*446637345000.d0 & 
            + w(k,5)*188806464000.d0 & 
            + w(k,7)*15076875000.d0 & 
            + w(k,8)*97599465000.d0   )
    enddo
235    continue
!c
!c           add 7 to the no of fcn evals
    c(24) = c(24) + 7.d0
!c
!c        end stage 2
!c
!c        ***************************************************************
!c        * stage 3 - calculate the error estimate est. first calculate *
!c        * the  unweighted  absolute  error  estimate vector (per unit *
!c        * step) for the unextrapolated approximation and store it  in *
!c        * w(*,2).  then  calculate the weighted max norm of w(*,2) as *
!c        * specified by the error  control  indicator  c(1).  finally, *
!c        * modify  this result to produce est, the error estimate (per *
!c        * unit step) for the extrapolated approximation ytrial.       *
!c        ***************************************************************
!c
!c           calculate the unweighted absolute error estimate vector
    do k = 1, n
       w(k,2) = (   w(k,1)*8738556750.d0 & 
            + w(k,3)*9735468750.d0 & 
            - w(k,4)*9709507500.d0 & 
            + w(k,5)*8582112000.d0 & 
            + w(k,6)*95329710000.d0 & 
            - w(k,7)*15076875000.d0 & 
            - w(k,8)*97599465000.d0)/1398169080000.d0
    enddo
300 continue
!c
!c           calculate the weighted max norm of w(*,2) as specified by
!c           the error control indicator c(1)
    temp = 0.d0
    if (c(1) .ne. 1.d0) go to 310
!c              absolute error control
    do  k = 1, n
       temp = dmax1(temp,dabs(w(k,2)))
    enddo
305 continue
    go to 360
310 if (c(1) .ne. 2.d0) go to 320
!c              relative error control
    do k = 1, n
       temp = dmax1(temp, dabs(w(k,2)/y(k)))
    enddo
315 continue
    go to 360
320 if (c(1) .ne. 3.d0) go to 330
!c              weights are 1/max(c(2),abs(y(k)))
    do k = 1, n
       temp = dmax1(temp, dabs(w(k,2)) & 
            / dmax1(c(2), dabs(y(k))) )
    enddo
325 continue
    go to 360
330 if (c(1) .ne. 4.d0) go to 340
!c              weights are 1/max(c(k+30),abs(y(k)))
    do k = 1, n
       temp = dmax1(temp, dabs(w(k,2)) & 
            / dmax1(c(k+30), dabs(y(k))) )
    enddo
335 continue
    go to 360
340 if (c(1) .ne. 5.d0) go to 350
!c              weights are 1/c(k+30)
    do k = 1, n
       temp = dmax1(temp, dabs(w(k,2)/c(k+30)))
    enddo
345 continue
    go to 360
350 continue
!c              default case - weights are 1/max(1,abs(y(k)))
    do k = 1, n
       temp = dmax1(temp, dabs(w(k,2)) & 
            / dmax1(1.d0, dabs(y(k))) )
    enddo
355 continue
360 continue
!c
!c           calculate est - (the weighted max norm of w(*,2))*hmag*scale
!c              - est is intended to be a measure of the error  per  unit
!c              step in ytrial
    c(19) = temp*c(14)*c(15)
!c
!c        end stage 3
!c
!c        ***************************************************************
!c        * stage 4 - make decisions.                                   *
!c        ***************************************************************
!c
!c           set ind=5 if step acceptable, else set ind=6
    ind = 5
    if (c(19) .gt. tol) ind = 6
!c
!c***********interrupt no 2 if requested
    if (c(9) .eq. 0.d0) go to 2222
    return
!c           resume here on re-entry with ind .eq. 5 or 6   ...re-entry..
2222 continue
!c
    if (ind .eq. 6) go to 410
!c              step accepted (ind .eq. 5), so update x, y from xtrial,
!c                 ytrial, add 1 to the no of successful steps, and set
!c                 the no of successive failures to zero
    x = c(17)
    do k = 1, n
       y(k) = w(k,9)
    enddo
400 continue
    c(22) = c(22) + 1.d0
    c(23) = 0.d0
!c**************return(with ind=3, xend saved, flag set) if x .eq. xend
    if (x .ne. xend) go to 405
    ind = 3
    c(20) = xend
    c(21) = 1.d0
    return
405 continue
    go to 420
410 continue
!c              step not accepted (ind .eq. 6), so add 1 to the no of
!c                 successive failures
    c(23) = c(23) + 1.d0
!c**************error return (with ind=-3) if hmag .le. hmin
    if (c(14) .gt. c(13)) go to 415
    ind = -3
    return
415 continue
420 continue
!c
!c        end stage 4
!c
    go to 99999
!c     end loop
!c
!c  begin abort action
500 continue
!c
    write(*,*) 'dverk crashes with ind =',  ind
    write(*,*) 'x = ', x
    stop
!c
!c  end abort action
!c
  end subroutine coop_dverk_with_cosmology



end module coop_ode_mod
