program potential
  use coop_wrapper
  implicit none
#include "constants.h"
#define C_1 args%r(1)
#define C_2 args%r(2)
#define N_POWER args%r(3)


  COOP_REAL max_R
  COOP_INT i
  COOP_INT, parameter::length=5000
  COOP_REAL R_val(length), der_f_val(length), R_left(length), R_right(length), V_left(length), V_right(length), phi_left(length), phi_right(length)

  COOP_REAL, parameter:: mu=1.d0
  COOP_REAL, parameter:: R_upper = 10000.d0*mu
  COOP_REAL c1,c2, f_R0, R_lower 
  COOP_REAL npower
  type(coop_function)::f, Vofphi
  type(coop_arguments)::args
  type(coop_asy)::plotV

  npower = 4.
  f_R0 = -1.d-3
  c2 = -18.*npower/(41.**(npower+1.)*f_R0)
  c1 = 18.*c2
  args = coop_arguments( r = (/ c1, c2, npower /) ) 
  R_lower = 30.*mu
  f = coop_function(f = broken, xmin = R_lower, xmax = R_upper, xlog=.true., ylog=.false., args=args)


!!$  do i = 1, length
!!$     R_val(i) = R_lower + dble(i-1)*(R_upper-R_lower)/dble(length-1)
!!$     der_f_val(i) = f%derivative(R_val(i))
!!$  enddo
!!$
!!$  max_R = R_val(minloc(der_f_val, 1))
!!$  print *, max_R
  
  max_R = 40.*mu

  do i = 1, length
    ! R_left(i) = R_lower + dble(i-1)*(max_R-R_lower)/(dble(length-1))
     R_right(i) = exp(log(max_R) + dble(i-1)*log(R_upper*0.99/max_R)/dble(length-1))
    ! V_left(i) = V(R_left(i))
     V_right(i) = V(R_right(i))
     print*, R_right(i)/mu, V_right(i)/mu
    ! phi_left(i) = sqrt(1.5d0)*log(f%derivative(R_left(i)))
     phi_right(i) = sqrt(1.5d0)*log(1.d0+f%derivative(R_right(i)))
  end do
     
  call plotV%open("plotV.txt")
  call plotV%init(xlabel = "$\phi$", ylabel = "$V(\phi)$",caption = "$f(R) = R -\mu \frac{c_1 (R/\mu)^n}{c_2 (R/\mu)^n+1}$" , height=6.)
 ! call coop_asy_interpolate_curve(plotV, xraw=phi_left(1:length), yraw=V_left(1:length), interpolate="LinearLinear", color = "red", linewidth = 1.2, linetype="solid", legend = "Late time")
  call coop_asy_interpolate_curve(plotV, xraw=phi_right(1:length), yraw=V_right(1:length)-9.d0,  color = "brown", interpolate="LinearLinear", linewidth = 1.2, linetype="solid", legend = "Early time")
  call coop_asy_legend(plotV)
  call plotV%close()
  

contains

  function V(R)
    implicit none
    COOP_REAL R, V

    V = (f%derivative(R)*R-f%eval(R))/(2.d0*(1.d0+f%derivative(R))**2.d0)
    
  end function V


  function broken(R, args)
    implicit none
    COOP_REAL R, broken
    type(coop_arguments) args

    broken =  - mu*(C_1*(R/mu)**N_POWER)/(C_2*(R/mu)**N_POWER + 1.)

  end function broken

end program potential
