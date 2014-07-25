program potential
  use coop_wrapper
  implicit none
#include "constants.h"
#define C_1 args%r(1)
#define C_2 args%r(2)
#define N_POWER args%r(3)


  COOP_REAL R_lower_used
  COOP_INT i
  COOP_INT, parameter::length=5000
  COOP_REAL R_val(length), der_f_val(length),  myR(length),  myV(length), myphi(length), myVp(length)

  COOP_REAL, parameter:: mu=1.d0
  COOP_REAL, parameter:: R_upper = 10000.d0*mu
  COOP_REAL c1,c2, f_R0, R_lower 
  COOP_REAL npower
  type(coop_function)::f, minusVofminusphi
  type(coop_arguments)::args
  type(coop_asy)::plotV

  npower = 4.
  f_R0 = -1.d-3
  c2 = -18.*npower/(41.**(npower+1.)*f_R0)
  c1 = 18.*c2
  args = coop_arguments( r = (/ c1, c2, npower /) ) 
  R_lower = 30.*mu
  f = coop_function(f = broken, xmin = R_lower, xmax = R_upper*1.01, xlog=.true., ylog=.false., args=args)

  
  R_lower_used = 40.*mu

  do i = 1, length
     myR(i) = exp(log(R_lower_used) + dble(i-1)*log(R_upper*0.99/R_lower_used)/dble(length-1))
     myV(i) = V(myR(i))
     myphi(i) = sqrt(1.5d0)*(f%derivative(myR(i)) - f%derivative(myR(i))**2/2)
     print*, myR(i)/mu, myphi(i), myV(i)/mu
  end do
  
  call minusVofminusphi%init_NonUniform(x=-myphi, f=-myV, xlog=.true., ylog = .true. )
  do i=3, length-2
     myVp(i) = -minusVofminusphi%derivative2(-myphi(i))
  enddo
  stop
  call plotV%open("plotVprime.txt")
  call plotV%init(xlabel = "$\phi$", ylabel = "$V'(\phi)$",caption = "$f(R) = R -\mu \frac{c_1 (R/\mu)^n}{c_2 (R/\mu)^n+1}$" , height=6., ylog=.true., xlog=.true.)

  call coop_asy_interpolate_curve(plotV, xraw=-myphi(3:length-2), yraw=myVp(3:length-2),  color = "brown", interpolate="LogLog", linewidth = 1.2, linetype="solid")
  call coop_asy_legend(plotV)
  call plotV%close()
  

contains

  function V(R)
    implicit none
    COOP_REAL R, V

    V = (f%derivative(R)*R-f%eval(R) + C_1/C_2*mu*(f%derivative(R))**2)/(2.d0*(1.d0+f%derivative(R))**2)  !! C_1/C_2*mu/2 removed
    
  end function V


  function broken(R, args)
    implicit none
    COOP_REAL R, broken
    type(coop_arguments) args

    broken =  - mu*(-C_1/C_2)/(C_2*(R/mu)**N_POWER + 1.) !!- (C_1/C_2)*mu removed

  end function broken

end program potential
