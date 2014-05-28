program test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  type(coop_species):: baryon, cdm, radiation, de
  type(coop_cosmology_background) bg
  type(coop_function) wofa
  type(coop_arguments) w0wa
  w0wa = coop_arguments( i=(/ 0, 1, 2, 3 /), r =  (/ COOP_REAL_OF(0.), COOP_REAL_OF(0.) /))
  wofa = coop_function(my_w_function, xmin = coop_min_scale_factor, xmax = COOP_REAL_OF(1.), xlog = .true., arg = w0wa)
  call baryon%init(name = "Baryon", id = 1, Omega = COOP_REAL_OF(0.046), w = COOP_REAL_OF(0.), cs2= COOP_REAL_OF(0.))
  call cdm%init(name = "CDM", id = 2, Omega = COOP_REAL_OF(0.26), w = COOP_REAL_OF(0.), cs2= COOP_REAL_OF(0.))
  call radiation%init(name = "Radiation", id = 3, Omega = COOP_REAL_OF(0.00008), w = COOP_REAL_OF(1./3.), cs2= COOP_REAL_OF(1./3.))
  call de%init(name = "DE", id = 4, Omega = COOP_REAL_OF(0.693), w = my_w_function(COOP_REAL_OF(1.), w0wa), cs2 = COOP_REAL_OF(1.), fw = wofa )
  call bg%init(name = "wCDM", id = 1)
  call bg%add_species(baryon)
  call bg%add_species(cdm)
  call bg%add_species(radiation)
  call bg%add_species(de)
  call bg%print
contains

  function my_w_function(a, arg) result(w)
    use coop_wrapper_typedef
    implicit none
    COOP_REAL a, w
    type(coop_arguments) arg
    w = arg%r(1) + arg%r(2)*(1.-a)
  end function my_w_function


end program test
