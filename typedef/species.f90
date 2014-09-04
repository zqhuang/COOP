module coop_species_mod
  use coop_constants_mod
  use coop_basicutils_mod
  use coop_particle_mod
  use coop_string_mod
  use coop_function_mod
  implicit none
#include "constants.h"

  private
  public:: coop_species

  type coop_species
     COOP_INT genre
     logical w_dynamic, cs2_dynamic
     COOP_SHORT_STRING name
     COOP_INT id
     COOP_REAL Omega, wp1, cs2
     COOP_REAL Omega_massless, mbyT
     type(coop_function):: fwp1, fcs2, fwp1eff
     type(coop_function)::flnrho
   contains
     procedure :: init => coop_species_initialize
     procedure :: print => coop_species_print
     procedure :: free  => coop_species_free
     procedure :: wofa => coop_species_wofa
     procedure :: dwda => coop_species_dwda
     procedure :: dlnrhodlna => coop_species_dlnrhodlna
     procedure :: wp1ofa => coop_species_wp1ofa
     procedure :: weffofa => coop_species_weffofa
     procedure :: wp1effofa => coop_species_wp1effofa
     procedure :: cs2ofa => coop_species_cs2ofa
     procedure :: density_ratio => coop_species_density_ratio
     procedure :: rhoa2_ratio => coop_species_rhoa2_ratio
     procedure :: rhoa3_ratio => coop_species_rhoa3_ratio
     procedure :: rhoa4_ratio => coop_species_rhoa4_ratio
     procedure :: density => coop_species_density
     procedure :: pressure => coop_species_pressure
     procedure :: rhoa2 => coop_species_rhoa2
     procedure :: pa2 => coop_species_pa2
  end type coop_species


  interface coop_species
     module procedure coop_species_constructor
  end interface coop_species

contains

  function coop_species_constructor(genre, name, id, Omega, w, cs2, Omega_massless, fwp1, fcs2, fwp1eff) result(this)
    type(coop_species) :: this
#include "species_init.h"
  end function coop_species_constructor

  subroutine coop_species_initialize(this, genre, name, id, Omega, w, cs2, Omega_massless, fwp1, fcs2, fwp1eff)
    class(coop_species) :: this
#include "species_init.h"
  end subroutine coop_species_initialize

  function coop_species_wofa(this, a) result(w)
    class(coop_species) :: this
    COOP_REAL w, a
    w = coop_species_wp1ofa(this, a) - 1.d0
  end function coop_species_wofa


  function coop_species_dwda(this, a) result(dwda)
    class(coop_species) :: this
    COOP_REAL dwda, a
    if(this%w_dynamic)then
       dwda = this%fwp1%derivative(COOP_PROPER_SCALE_FACTOR(a))
    else
       dwda = 0.d0
    endif
  end function coop_species_dwda


  function coop_species_wp1ofa(this, a) result(wp1)
    class(coop_species) :: this
    COOP_REAL wp1, a
    if(this%w_dynamic)then
       wp1 = this%fwp1%eval(COOP_PROPER_SCALE_FACTOR(a))
    else
       wp1 = this%wp1
    endif
  end function coop_species_wp1ofa


  function coop_species_weffofa(this, a) result(w)
    class(coop_species) :: this
    COOP_REAL w, a
    w = coop_species_wp1effofa(this, a) - 1.d0
  end function coop_species_weffofa

  function coop_species_wp1effofa(this, a) result(wp1)
    class(coop_species) :: this
    COOP_REAL wp1, a
    if(this%w_dynamic)then
       wp1 = this%fwp1eff%eval(COOP_PROPER_SCALE_FACTOR(a))
    else
       wp1 = this%wp1
    endif
  end function coop_species_wp1effofa


  function coop_species_cs2ofa(this, a) result(cs2)
    class(coop_species) :: this
    COOP_REAL cs2, a
    if(this%cs2_dynamic)then
       cs2 = this%fcs2%eval(COOP_PROPER_SCALE_FACTOR(a))
    else
       cs2 = this%cs2
    endif
  end function coop_species_cs2ofa


  subroutine coop_species_print(this)
    class(coop_species) :: this
    write(*,"(A)") "Species Name: "//trim(this%name)
    write(*,"(A)") "Species ID: "//trim(coop_num2str(this%id))
    select case(this%genre)
    case(COOP_SPECIES_MASSLESS)
       write(*,"(A)") "Species genre: massless particles"
    case(COOP_SPECIES_MASSIVE_FERMION)
       write(*,"(A)") "Species genre: massive fermion"
       write(*,"(A, G15.6)") "Species Omega_massless: ", this%Omega_massless
    case(COOP_SPECIES_MASSIVE_BOSON)
       write(*,"(A)") "Species genre: massive boson"
       write(*,"(A, G15.6)") "Species Omega_massless: ", this%Omega_massless
    case(COOP_SPECIES_COSMOLOGICAL_CONSTANT)
       write(*,"(A)") "Species genre: cosmological constant"
    case(COOP_SPECIES_FLUID)
       write(*,"(A)") "Species genre: fluid"
    case(COOP_SPECIES_SCALAR_FIELD)
       write(*,"(A)") "Species genre: scalar field"
    case(COOP_SPECIES_CDM)
       write(*,"(A)") "Species genre: cold dark matter"
    case default
       write(*,"(A)") "Species genre: unknown"
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
       write(*,"(A, G14.5)") "Species w: ", this%wp1-1.d0
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

  function coop_species_rhoa4_ratio(this, a) result(rhoa4)
    class(coop_species)::this
    COOP_REAL a, an
    COOP_REAL rhoa4
    if(this%w_dynamic)then
       an = COOP_PROPER_SCALE_FACTOR(a)
       rhoa4 = exp(this%flnrho%eval(an)+4.d0*log(an))
    else
       rhoa4 = max(a, 1.d-99)**(4.d0-3.d0*this%wp1)
    endif
  end function coop_species_rhoa4_ratio

  function coop_species_rhoa3_ratio(this, a) result(rhoa3)
    class(coop_species)::this
    COOP_REAL a, an
    COOP_REAL rhoa3
    if(this%w_dynamic)then
       an = COOP_PROPER_SCALE_FACTOR(a)
       rhoa3 = exp(this%flnrho%eval(an)+3.d0*log(an))
    else
       rhoa3 = max(a, 1.d-99)**(3.d0-3.d0*this%wp1)
    endif
  end function coop_species_rhoa3_ratio

  function coop_species_rhoa2_ratio(this, a) result(rhoa2)
    class(coop_species)::this
    COOP_REAL a, an
    COOP_REAL rhoa2
    if(this%w_dynamic)then
       an = COOP_PROPER_SCALE_FACTOR(a)
       rhoa2 = exp(this%flnrho%eval(an)+2.d0*log(an))
    else
       rhoa2 = max(a, 1.d-99)**(2.d0-3.d0*this%wp1)
    endif
  end function coop_species_rhoa2_ratio



  function coop_species_density(this, a) result(density)  !!unit H_0^2M_p^2
    class(coop_species)::this
    COOP_REAL a
    COOP_REAL density
    density = 3.d0*this%Omega * this%density_ratio(a)
  end function coop_species_density


  function coop_species_rhoa2(this, a) result(density)  !!unit H_0^2M_p^2
    class(coop_species)::this
    COOP_REAL a
    COOP_REAL density
    density = 3.d0*this%Omega * this%rhoa2_ratio(a)
  end function coop_species_rhoa2


  function coop_species_dlnrhodlna(this, a) result(dlnrhodlna)  
    class(coop_species)::this
    COOP_REAL a
    COOP_REAL dlnrhodlna
    if(this%w_dynamic)then
       dlnrhodlna = this%flnrho%derivative_bare(log(a))
    else
       dlnrhodlna = -3.*this%wp1
    endif
  end function coop_species_dlnrhodlna

  function coop_species_drhoa2da(this, a) result(drhoa2da)
    class(coop_species)::this
    COOP_REAL a, drhoa2da
    drhoa2da = (this%dlnrhodlna(a) + 2.d0)/a*this%rhoa2(a)
  end function coop_species_drhoa2da


  function coop_species_pressure(this, a) result(p)  !!unit H_0^2M_p^2
    class(coop_species)::this
    COOP_REAL a
    COOP_REAL p
    p = this%wofa(a) * this%density(a)
  end function coop_species_pressure


  function coop_species_pa2(this, a) result(p)  !!unit H_0^2M_p^2
    class(coop_species)::this
    COOP_REAL a
    COOP_REAL p
    p = this%wofa(a) * this%rhoa2(a)
  end function coop_species_pa2

  function coop_species_dpa2da(this, a) result(dpa2da)
    class(coop_species)::this
    COOP_REAL a, dpa2da
    dpa2da = ((this%dlnrhodlna(a) + 2.d0)/a*this%wofa(a) + this%dwda(a))*this%rhoa2(a)
  end function coop_species_dpa2da

  function coop_species_density_ratio(this, a) result(density)
    class(coop_species)::this
    COOP_REAL a
    COOP_REAL density
    if(this%w_dynamic)then
       density = exp(this%flnrho%eval(COOP_PROPER_SCALE_FACTOR(a)))
    else
       density = a**(-3.*this%wp1)
    endif
  end function coop_species_density_ratio

  subroutine coop_species_free(this)
    class(coop_species)::this
    call this%flnrho%free()
    call this%fwp1%free()
    call this%fcs2%free()
    call this%fwp1eff%free()
  end subroutine coop_species_free



end module coop_species_mod
