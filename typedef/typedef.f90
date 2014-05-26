module coop_types
  use coop_constants
  use coop_basictype
  use coop_basicutils
  implicit none
#include "constants.h"

private

public:: coop_species_basic, coop_species_dynamic, coop_cosmology_null, coop_cosmology_background



  type coop_species_basic
     COOP_SHORT_STRING name
     COOP_INT id
     COOP_REAL Omega, w, cs2
     COOP_REAL params(coop_max_num_DE_params)
   contains
     procedure :: init => coop_species_initialize
     procedure :: print => coop_species_print
     procedure :: density_ratio => coop_species_density_ratio
  end type coop_species_basic


  type, extends(coop_species_basic),abstract  :: coop_species_dynamic
     COOP_REAL_ARRAY :: ln_rho_ratio
   contains
     procedure(coop_wofa), deferred:: wofa
     procedure(coop_cs2ofa), deferred:: cs2ofa
  end type coop_species_dynamic

  abstract interface 
     function coop_wofa(this, a) result(wofa)
       use coop_constants
       import coop_species_dynamic
       class(coop_species_dynamic):: this
       COOP_REAL a, wofa
     end function coop_wofa

     function coop_cs2ofa(this, a) result(cs2ofa)
       use coop_constants
       import coop_species_dynamic
       class(coop_species_dynamic):: this
       COOP_REAL a, cs2ofa
     end function coop_cs2ofa
  end interface

  type coop_cosmology_null
     COOP_SHORT_STRING  name
     COOP_INT id
   contains
     procedure :: init => coop_cosmology_initialize
     procedure :: print=> coop_cosmology_print
  end type coop_cosmology_null

  type, extends(coop_cosmology_null):: coop_cosmology_background
     COOP_INT num_species
     class(coop_species_basic),dimension(coop_max_num_species), pointer ::species
   contains
     procedure::add_species=>coop_cosmology_background_add_species
  end type coop_cosmology_background


  interface coop_species_basic
     procedure coop_species_basic_constructor
  end interface

  interface coop_cosmology_background
     procedure coop_cosmology_background_constructor
  end interface

  interface Coop_cosmology_null
     procedure coop_cosmology_null_construct
  end interface



contains

  function coop_cosmology_null_construct() result(cosmo)
    type(coop_cosmology_null) :: cosmo
    call cosmo%init()
  end function coop_cosmology_null_construct

  subroutine coop_cosmology_background_add_species(this, species)
    class(coop_cosmology_background)::this
    class(coop_species_basic):: species
    if(this%num_species .ge. coop_max_num_species) stop "coop_cosmology_background_add_species: too many species"
    this%num_species = this%num_species+1
    allocate(this%species(this%num_species:this%num_species), source=species)
  end subroutine coop_cosmology_background_add_species

  function coop_cosmology_background_constructor() result(cosmo)
    type(coop_cosmology_background) cosmo
    call cosmo%init()
  end function coop_cosmology_background_constructor

  subroutine coop_cosmology_initialize(this, name, id)
    class(coop_cosmology_null)::this
    COOP_UNKNOWN_STRING, optional::name
    COOP_INT, optional::id
    if(present(name))then
       this%name = name
    else
       this%name = "COOP_COSMOLOGY"
    endif
    if(present(id))then
       this%id =  id
    else
       this%id = 0
    endif
    select type (this)
    class is (coop_cosmology_background)
       deallocate(this%species)
       cosmo%num_species = 0
    end select
  end subroutine coop_cosmology_initialize

  subroutine coop_cosmology_print(this)
    class(coop_cosmology_null)::this
    integer i
    write(*,"(A)") "================================="
    write(*,"(A)") "Cosmology Class: background"
    write(*,"(A)") "Cosmology Name = "//trim(this%name)
    write(*,"(A, I8)") "Cosmology id = ", this%id
    select type(this)
    class is (coop_cosmology_background)
       do i=1, cosmo%num_species
          write(*,"(A)") "---------------------------------"
          write(*,"(A, I5)") "Species #:", i
          call cosmo%species(i)%print
       enddo
       write(*,"(A)") "================================="
    end select
  end subroutine coop_cosmology_print

  function coop_species_basic_constructor()  result(species)
    type(coop_species_basic) species
    call species%init()
  end function coop_species_basic_constructor

  subroutine coop_species_initialize(this, name, id, Omega, w, cs2, params)
    class(coop_species_basic) :: this
    COOP_UNKNOWN_STRING, optional::name
    COOP_INT, optional:: id
    COOP_REAL, optional::Omega
    COOP_REAL, optional::w
    COOP_REAL, optional::cs2
    COOP_REAL, dimension(coop_max_num_de_params):: params
    if(present(name))then
       this%name = name
    else
       this%name = ""
    endif
    if(present(id))then
       this%id = id
    else
       this%id = 0
    endif
    if(present(Omega))then
       this%Omega = Omega
    else
       this%Omega = 1
    endif
    if(present(w))then
       this%w = w
    else
       this%w = 0
    endif
    if(present(cs2))then
       this%cs2 = cs2
    else
       this%cs2 = 0.
    endif
    if(present(params))then
       this%params = params
    else
       this%params = 0.
    endif
  end subroutine coop_species_initialize

  subroutine coop_species_print(this)
    class(coop_species_basic):: this
    write(*,"(A)") "Species Name: "//trim(this%name)
    write(*,"(A, I8)") "Species ID: ", this%id
    write(*,"(A, G15.6)") "Species Omega: ", this%Omega
    select type(this)
    type is (coop_species_basic)
       write(*,"(A, G15.6)") "Species w: ", this%w
       write(*,"(A, G15.6)") "Species cs^2: ", this%cs2
    class is (coop_species_dynamic)
       write(*,"(A)") "Species w(a=0.00001) = ", this%wofa(COOP_REAL_OF(0.00001))
       write(*,"(A)") "Species w(a=0.01) = ", this%wofa(COOP_REAL_OF(0.01))
       write(*,"(A)") "Species w(a=0.1) = ", this%wofa(COOP_REAL_OF(0.1))
       write(*,"(A)") "Species w(a=0.5) = ", this%wofa(COOP_REAL_OF(0.5))
       write(*,"(A)") "Species w(a=1) = ", this%wofa(COOP_REAL_OF(1.))
       write(*,"(A)") "Species cs^2(a=0.00001) = ", this%cs2ofa(COOP_REAL_OF(0.00001))
       write(*,"(A)") "Species cs^2(a=0.01) = ", this%cs2ofa(COOP_REAL_OF(0.01))
       write(*,"(A)") "Species cs^2(a=0.1) = ", this%cs2ofa(COOP_REAL_OF(0.1))
       write(*,"(A)") "Species cs^2(a=0.5) = ", this%cs2ofa(COOP_REAL_OF(0.5))
       write(*,"(A)") "Species cs^2(a=1) = ", this%cs2ofa(COOP_REAL_OF(1.))
    end select
  end subroutine coop_species_print

  function coop_species_density_ratio(this, a) result(density)
    class(coop_species_basic)::this
    COOP_REAL density
    select type(this)
    type is (coop_species_basic)
       density = this%Omega * a**(-3.*(1.+this%w))
    class is(coop_species_dynamic)

    end select
  end function coop_species_density_ratio





end module coop_types
