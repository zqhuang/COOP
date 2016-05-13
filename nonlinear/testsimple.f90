program test
  use coop_wrapper_firstorder
  use coop_ellipse_collapse_mod
  implicit none
#include "constants.h"
  COOP_REAL::delta, fr, theta_low, theta_high, theta_mid, r_low, r_high, r_mid
  call coop_get_command_line_argument(key="fpk", arg=delta)
  call coop_get_command_line_argument(key="fr", arg=fr)
  theta_low = coop_pi
  theta_high = coop_2pi
  r_low = radius(theta_low)
  r_high = radius(theta_high)
  if(r_low .lt. fr .or. r_high .gt. fr) stop "invalid fr (0<fr<0.56)"
  do while(theta_high - theta_low .gt. 1.d-10)
     theta_mid = (theta_high+theta_low)/2.d0
     r_mid = radius(theta_mid)
     if(r_mid .lt. fr)then
        theta_high = theta_mid
     else
        theta_low = theta_mid
     endif
  enddo
  theta_mid = (theta_high + theta_low)/2.d0
  print*, "z_collapse = ", redshift(theta_mid)
contains
  function radius(theta)
    COOP_REAL::radius, theta
    radius = (1.d0-cos(theta))/(theta-sin(theta))**(2.d0/3.d0)*2.d0/6.d0**(2.d0/3.d0) 
  end function radius

  function redshift(theta)
    COOP_REAL::redshift, theta
    redshift = 1.d0/((3.d0/4.d0)**(2.d0/3.d0)*(theta-sin(theta))**(2.d0/3.d0)/(5.d0/3.d0)/delta) - 1.d0
  end function redshift

end program test
