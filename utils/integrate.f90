module coop_integrate_mod
  use coop_wrapper_typedef
  implicit none


#include "constants.h"

  !!your function should be like
  !! COOP_REAL function(COOP_REAL x)
  !!   call coop_integrate(func, a, b, precision)
  !!
  !! COOP_REAL function(COOP_REAL x, type(coop_arguments) args)
  !!   call coop_integrate(func, a, b, args, precision) 
  !!
  !! COOP_REAL function(COOP_REAL x, COOP_REAL y)
  !!   call coop_integrate(func, xmin, xmax, ymin, ymax, precision)
  !!
  !! COOP_REAL function(COOP_REAL x, COOP_REAL y,  type(coop_arguments) args)    
  !!   call coop_integrate(func, xmin, xmax, ymin, ymax, args, precision)
  
  private

  public::coop_integrate, coop_polint

  interface coop_integrate
     module procedure coop_qrombc, coop_qromb_with_arguments, coop_2d_integrate, coop_2d_integrate_with_arguments
  end interface coop_integrate


contains


  function coop_qrombc(func,a, b, precision) result(integral)
#define QROMB_ARGUMENTS
#include "qromb.h"    
#undef QROMB_ARGUMENTS    
  end function coop_qrombc

!!func(x, args)
  function coop_qromb_with_arguments(func, a, b, args, precision) result(integral)
    type(coop_arguments)::args
#define QROMB_ARGUMENTS  ,args
#include "qromb.h"
#undef QROMB_ARGUMENTS
  end function coop_qromb_with_arguments

  function coop_2d_integrate(func, xmin, xmax, ymin, ymax, precision) result(integral)
    external func
    COOP_REAL integral
    COOP_REAL xmin, xmax, ymin, ymax, func
    COOP_REAL, optional::precision
    if(present(precision))then
       if(precision .le. 1.d-6)then
          integral = coop_2d_integrate_128(func, xmin, xmax, ymin, ymax)          
       elseif(precision .le. 1.d-5)then
          integral = coop_2d_integrate_64(func, xmin, xmax, ymin, ymax)                    
       elseif(precision .le. 1.d-4)then
          integral = coop_2d_integrate_32(func, xmin, xmax, ymin, ymax)                    
       elseif(precision .le. 1.d-3)then
          integral = coop_2d_integrate_16(func, xmin, xmax, ymin, ymax)                    
       else
          integral = coop_2d_integrate_8(func, xmin, xmax, ymin, ymax)                    
       endif
    else
       integral = coop_2d_integrate_32(func, xmin, xmax, ymin, ymax)
    endif
  end function coop_2d_integrate


  function coop_2d_integrate_with_arguments(func, xmin, xmax, ymin, ymax, args, precision) result(integral)
    external func
    COOP_REAL integral
    COOP_REAL xmin, xmax, ymin, ymax, func
    type(coop_arguments)::args    
    COOP_REAL, optional::precision
    if(present(precision))then
       if(precision .le. 1.d-6)then
          integral = coop_2d_integrate_with_arguments_128(func, xmin, xmax, ymin, ymax, args)          
       elseif(precision .le. 1.d-5)then
          integral = coop_2d_integrate_with_arguments_64(func, xmin, xmax, ymin, ymax, args)                    
       elseif(precision .le. 1.d-4)then
          integral = coop_2d_integrate_with_arguments_32(func, xmin, xmax, ymin, ymax, args)                    
       elseif(precision .le. 1.d-3)then
          integral = coop_2d_integrate_with_arguments_16(func, xmin, xmax, ymin, ymax, args)                    
       else
          integral = coop_2d_integrate_with_arguments_8(func, xmin, xmax, ymin, ymax, args)                    
       endif
    else
       integral = coop_2d_integrate_with_arguments_32(func, xmin, xmax, ymin, ymax, args)
    endif
  end function coop_2d_integrate_with_arguments


  function coop_2d_integrate_8(func, xmin, xmax, ymin, ymax) result(integral)
    type(coop_arguments)::args
#define INT2D_ARGUMENTS 
#define INT2D_N_BASE  8
#include "int2d.h"
#undef INT2D_ARGUMENTS
#undef INT2D_N_BASE    
  end function coop_2d_integrate_8
  
  
  function coop_2d_integrate_with_arguments_8(func, xmin, xmax, ymin, ymax, args) result(integral)
    type(coop_arguments)::args
#define INT2D_ARGUMENTS ,args
#define INT2D_N_BASE  8
#include "int2d.h"
#undef INT2D_ARGUMENTS
#undef INT2D_N_BASE
  end function coop_2d_integrate_with_arguments_8
  

  function coop_2d_integrate_16(func, xmin, xmax, ymin, ymax) result(integral)
    type(coop_arguments)::args
#define INT2D_ARGUMENTS 
#define INT2D_N_BASE  16
#include "int2d.h"
#undef INT2D_ARGUMENTS
#undef INT2D_N_BASE    
  end function coop_2d_integrate_16
  
  
  function coop_2d_integrate_with_arguments_16(func, xmin, xmax, ymin, ymax, args) result(integral)
    type(coop_arguments)::args
#define INT2D_ARGUMENTS ,args
#define INT2D_N_BASE  16
#include "int2d.h"
#undef INT2D_ARGUMENTS
#undef INT2D_N_BASE
  end function coop_2d_integrate_with_arguments_16

  
  function coop_2d_integrate_32(func, xmin, xmax, ymin, ymax) result(integral)
    type(coop_arguments)::args
#define INT2D_ARGUMENTS 
#define INT2D_N_BASE  32
#include "int2d.h"
#undef INT2D_ARGUMENTS
#undef INT2D_N_BASE    
  end function coop_2d_integrate_32
  
  
  function coop_2d_integrate_with_arguments_32(func, xmin, xmax, ymin, ymax, args) result(integral)
    type(coop_arguments)::args
#define INT2D_ARGUMENTS ,args
#define INT2D_N_BASE  32
#include "int2d.h"
#undef INT2D_ARGUMENTS
#undef INT2D_N_BASE
  end function coop_2d_integrate_with_arguments_32

  function coop_2d_integrate_64(func, xmin, xmax, ymin, ymax) result(integral)
    type(coop_arguments)::args
#define INT2D_ARGUMENTS 
#define INT2D_N_BASE  64
#include "int2d.h"
#undef INT2D_ARGUMENTS
#undef INT2D_N_BASE    
  end function coop_2d_integrate_64
  
  
  function coop_2d_integrate_with_arguments_64(func, xmin, xmax, ymin, ymax, args) result(integral)
    type(coop_arguments)::args
#define INT2D_ARGUMENTS ,args
#define INT2D_N_BASE  64
#include "int2d.h"
#undef INT2D_ARGUMENTS
#undef INT2D_N_BASE
  end function coop_2d_integrate_with_arguments_64


  function coop_2d_integrate_128(func, xmin, xmax, ymin, ymax) result(integral)
    type(coop_arguments)::args
#define INT2D_ARGUMENTS 
#define INT2D_N_BASE  128
#include "int2d.h"
#undef INT2D_ARGUMENTS
#undef INT2D_N_BASE    
  end function coop_2d_integrate_128
  
  
  function coop_2d_integrate_with_arguments_128(func, xmin, xmax, ymin, ymax, args) result(integral)
    type(coop_arguments)::args
#define INT2D_ARGUMENTS ,args
#define INT2D_N_BASE  128
#include "int2d.h"
#undef INT2D_ARGUMENTS
#undef INT2D_N_BASE
  end function coop_2d_integrate_with_arguments_128
  
  
 



  subroutine Coop_polint(xa,ya,n,x,y,dy)
    COOP_INT  :: n,NMAX
    COOP_REAL  :: dy,x,y,xa(n),ya(n)
    parameter (NMAX=10)
    COOP_INT  :: i,m,ns
    COOP_REAL  :: den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
    ns=1
    dif=ABS(x-xa(1))
    do i=1,n
       dift=ABS(x-xa(i))
       if (dift.LT.dif) then
          ns=i
          dif=dift
       endif
       c(i)=ya(i)
       d(i)=ya(i)
    enddo
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
       do i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) then
             write(*,*) 'failure in polint'
             STOP
          endif
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
       enddo
       if (2*ns.LT.n-m)then
          dy=c(ns+1)
       ELSE
          dy=d(ns)
          ns=ns-1
       endif
       y=y+dy
    enddo
    return
  end subroutine Coop_polint


end module coop_integrate_mod
