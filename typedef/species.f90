module coop_type_species
  use coop_constants
  use coop_basicutils
  use coop_particle
  use coop_string
  use coop_type_function
  implicit none
#include "constants.h"

  private

  public:: coop_species

  type coop_species
     COOP_INT gengre
     logical w_dynamic, cs2_dynamic
     COOP_SHORT_STRING name
     COOP_INT id
     COOP_REAL Omega, w, cs2
     COOP_REAL Omega_massless, mbyT
     type(coop_function),pointer:: fw, fcs2
     type(coop_function)::frho
   contains
     procedure :: init => coop_species_initialize
     procedure :: print => coop_species_print
     procedure :: free  => coop_species_free
     procedure :: wofa => coop_species_wofa
     procedure :: cs2ofa => coop_species_cs2ofa
     procedure :: density_ratio => coop_species_density_ratio
  end type coop_species


  interface coop_species
     procedure coop_species_constructor
  end interface coop_species


contains

  function coop_species_constructor(gengre, name, id, Omega, w, cs2, Omega_massless, fw, fcs2) result(this)
    type(coop_species) :: this
#include "species_init.h"
  end function coop_species_constructor

  subroutine coop_species_initialize(this, gengre, name, id, Omega, w, cs2, Omega_massless, fw, fcs2)
    class(coop_species) :: this
#include "species_init.h"
  end subroutine coop_species_initialize

  function coop_species_wofa(this, a) result(w)
    class(coop_species)::this
    COOP_REAL w, a
    if(this%w_dynamic)then
       w = this%fw%eval(a)
    else
       w = this%w
    endif
  end function coop_species_wofa


  function coop_species_cs2ofa(this, a) result(cs2)
    class(coop_species)::this
    COOP_REAL cs2, a
    if(this%cs2_dynamic)then
       cs2 = this%fcs2%eval(a)
    else
       cs2 = this%cs2
    endif
  end function coop_species_cs2ofa


  subroutine coop_species_print(this)
    class(coop_species):: this
    write(*,"(A)") "Species Name: "//trim(this%name)
    write(*,"(A)") "Species ID: "//trim(coop_num2str(this%id))
    select case(this%gengre)
    case(COOP_SPECIES_MASSLESS)
       write(*,"(A)") "Species gengre: massless particles"
    case(COOP_SPECIES_MASSIVE_FERMION)
       write(*,"(A)") "Species gengre: massive fermion"
       write(*,"(A, G15.6)") "Species Omega_massless: ", this%Omega_massless
    case(COOP_SPECIES_MASSIVE_BOSON)
       write(*,"(A)") "Species gengre: massive boson"
       write(*,"(A, G15.6)") "Species Omega_massless: ", this%Omega_massless
    case(COOP_SPECIES_COSMOLOGICAL_CONSTANT)
       write(*,"(A)") "Species gengre: cosmological constant"
    case(COOP_SPECIES_FLUID)
       write(*,"(A)") "Species gengre: fluid"
    case(COOP_SPECIES_SCALAR_FIELD)
       write(*,"(A)") "Species gengre: scalar field"
    case(COOP_SPECIES_CDM)
       write(*,"(A)") "Species gengre: cold dark matter"
    case default
       write(*,"(A)") "Species gengre: unknown"
    end select
    write(*,"(A, G15.6)") "Species Omega: ", this%Omega
    if(this%w_dynamic)then
       write(*,"(A, G14.5)") "Species w(z=10000) = ", this%wofa(COOP_REAL_OF(1.d0/10001.d0))
       write(*,"(A, G14.5)") "Species w(z=1000) = ", this%wofa(COOP_REAL_OF(1.d0/1001.d0))
       write(*,"(A, G14.5)") "Species w(z=10) = ", this%wofa(COOP_REAL_OF(1.d0/11.d0))
       write(*,"(A, G14.5)") "Species w(z=1) = ", this%wofa(COOP_REAL_OF(0.5d0))
       write(*,"(A, G14.5)") "Species w(z=0.5) = ", this%wofa(COOP_REAL_OF(2.d0/3.d0))
       write(*,"(A, G14.5)") "Species w(z=0) = ", this%wofa(COOP_REAL_OF(1.))
       write(*,"(A, G14.5)") "Species rho_ratio(z=10000) = ", this%density_ratio(COOP_REAL_OF(1.d0/10001.d0))
       write(*,"(A, G14.5)") "Species rho_ratio(z=1000) = ", this%density_ratio(COOP_REAL_OF(1.d0/1001.d0))
       write(*,"(A, G14.5)") "Species rho_ratio(z=10) = ", this%density_ratio(COOP_REAL_OF(1.d0/11.d0))
       write(*,"(A, G14.5)") "Species rho_ratio(z=1) = ", this%density_ratio(COOP_REAL_OF(0.5d0))
       write(*,"(A, G14.5)") "Species rho_ratio(z=0.5) = ", this%density_ratio(COOP_REAL_OF(2.d0/3.d0))
       write(*,"(A, G14.5)") "Species rho_ratio(z=0) = ", this%density_ratio(COOP_REAL_OF(1.d0))
       
    else
       write(*,"(A, G14.5)") "Species w: ", this%w
    endif
    if(this%cs2_dynamic)then
       write(*,"(A, G14.5)") "Species cs^2(z=10000) = ", this%cs2ofa(COOP_REAL_OF(1.d0/10001.d0))
       write(*,"(A, G14.5)") "Species cs^2(z=1000) = ", this%cs2ofa(COOP_REAL_OF(1.d0/1001.d0))
       write(*,"(A, G14.5)") "Species cs^2(z=10) = ", this%cs2ofa(COOP_REAL_OF(1.d0/11.d0))
       write(*,"(A, G14.5)") "Species cs^2(z=1) = ", this%cs2ofa(COOP_REAL_OF(0.5d0))
       write(*,"(A, G14.5)") "Species cs^2(z=0.5) = ", this%cs2ofa(COOP_REAL_OF(2.d0/3.d0))
       write(*,"(A, G14.5)") "Species cs^2(z=0) = ", this%cs2ofa(COOP_REAL_OF(1.))
    else
       write(*,"(A, G14.5)") "Species cs^2: ", this%cs2
    endif
  end subroutine coop_species_print

  function coop_species_density_ratio(this, a) result(density)
    class(coop_species)::this
    COOP_REAL a
    COOP_REAL density
    if(this%w_dynamic)then
       density = this%frho%eval(a)
    else
       density = a**(-3.*(1.+this%w))
    endif
  end function coop_species_density_ratio

  subroutine coop_species_free(this)
    class(coop_species)::this
    call this%fw%free()
    call this%fcs2%free()
    call this%frho%free()
  end subroutine coop_species_free

end module coop_type_species
