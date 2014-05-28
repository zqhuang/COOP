program test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  type(coop_species):: baryon, cdm, radiation, de
  type(coop_cosmology_background) bg
  type(coop_function) wofa
  type(coop_arguments) w0wa
  w0wa = coop_arguments( r =  (/ COOP_REAL_OF(-1.), COOP_REAL_OF(0.2) /))
  wofa = coop_function(w0wa, xmin = COOP_REAL_OF(1.e-8), xmax = COOP_REAL_OF(1.), xlog = .true. )

  allocate(w0wa%r(w0wa%n_real))
  call baryon%init(name = "Baryon", id = 1, Omega = COOP_REAL_OF(0.046), w = COOP_REAL_OF(0.), cs2= COOP_REAL_OF(0.))
  call baryon%print
  call cdm%init(name = "CDM", id = 2, Omega = COOP_REAL_OF(0.26), w = COOP_REAL_OF(0.), cs2= COOP_REAL_OF(0.), fw = wofa)
  call cdm%print
  call radiation%init(name = "Radiation", id = 3, Omega = COOP_REAL_OF(0.00008), w = COOP_REAL_OF(1./3.), cs2= COOP_REAL_OF(1./3.))
  call radiation%print
  call de%init(name = "DE", id = 4, Omega = COOP_REAL_OF(0.693), w = wofa_w0wa(COOP_REAL_OF(1.)), cs2 = COOP_REAL_OF(1.), fw = wofa )
  call de%print
  call bg%init
  call bg%print
  call bg%add_species(baryon)
  call bg%print
  call bg%add_species(cdm)
  call bg%print
  call bg%add_species(radiation)
  call bg%print

contains

  function wofa_w0wa(a, arg) result(w)
    use coop_wrapper_typedef
    implicit none
    COOP_REAL a, w
    type(coop_arguments) arg
    w = arg%r(1) + arg%r(2)*(1.-a)
  end function wofa_w0wa


end program test
