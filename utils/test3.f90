<<<<<<< HEAD


program searchzeq
  use coop_wrapper_utils
  print*,Bessel_jn(10, 20.d0)
  print*,Bessel_jn(10, 100.d0)  
=======
program searchzeq
  use coop_wrapper_utils
#include "constants.h"
  COOP_REAL::omega_m, w, z
  read(*,*) omega_m, w, z
  print*, coop_fGrowth_fitting(Omega_m, w, z),  coop_fGrowth_naive(Omega_m, w, z) 

contains

    !!for quick tests; d\ln D/d\ln a fitting formula for wCDM cosmology
  function coop_fGrowth_fitting(Omega_m, w, z) result(f)
    COOP_REAL::Omega_m, w ,z ,f, f1, f2, f3
    f1 = 2.d0*(1.d0-2.d0*w)*(2.d0-3.d0*w)*(1.d0-Omega_m)
    f2 =  - w*(5.d0-6.d0*w)*(4.d0+w)*Omega_m
    f3 = (1.d0+z) ** (3.d0 * w)
    f = 1.d0 - 1.5*w*f3*(1.d0-Omega_m)/(Omega_m+(1.d0-Omega_m)*f3) + (3.d0+w*0.75d0)*w*f3*f1/(f1*f3+f2)
  end function coop_fGrowth_fitting


  function coop_fGrowth_naive(Omega_m, w, z) result(f)
    COOP_REAL::omm, f
    COOP_REAL::omega_m, w, z
    omm = Omega_m*(1.d0+z)**3
    omm = omm/(omm + (1.d0-Omega_m)*(1.d0+z)**(3.d0*(1.d0+w)))
    f = omm**0.56
  end function coop_fGrowth_naive

  
>>>>>>> b5d294d1d639acfc5d425cf9967c5415ac38dd4c
end program searchzeq

