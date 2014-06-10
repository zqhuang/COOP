program derivative
  use coop_wrapper
  implicit none
#include "constants.h"
  COOP_REAL, parameter::Mpl = 1.d0
  COOP_REAL, parameter::H0 = 1.d0
  type(coop_function)::V
  type(coop_arguments)::args
  COOP_REAL phi

#define MASS args%r(1) 
#define LAMBDA args%r(2)

  args = coop_arguments( r = (/ 1.e-5*Mpl, 3.*Mpl**2*H0**2*0.7 /) )
  V = coop_function( f = potential, xmin = MASS/10., xmax =Mpl*10.d0, xlog = .true., ylog = .true., args = args, method = COOP_INTERPOLATE_SPLINE)

  phi = 0.2*MASS
  print*, V%eval(phi), potential(phi, args)
  print*, V%derivative(phi), dVdphi(phi, args)
  print*, V%derivative2(phi), d2Vdphi2(phi, args)

  print*
  phi = Mpl
  print*, V%eval(phi), potential(phi, args)
  print*, V%derivative(phi), dVdphi(phi, args)
  print*, V%derivative2(phi), d2Vdphi2(phi, args)

contains

  function potential(phi, args) result(V)
    type(coop_arguments) args
    COOP_REAL V, phi
    V = exp(-(MASS/phi)**2)*(MASS/phi)**2 * LAMBDA
  end function potential

  function dVdphi(phi, args) result(Vp)
    type(coop_arguments) args
    COOP_REAL Vp, phi
    Vp = (2.d0*MASS**2/phi**3 - 2.d0/phi)*potential(phi, args)
  end function dVdphi

  function d2Vdphi2(phi,args) result(Vpp)
    COOP_REAL Vpp, phi
    type(coop_arguments) args
    Vpp = ((2.d0*MASS**2/phi**3 - 2.d0/phi)**2 - 6.d0*MASS**2/phi**4 + 2.d0/phi**2)* potential(phi, args)

  end function d2Vdphi2

end program derivative
