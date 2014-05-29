module coop_type_cosmology
  use coop_constants
  use coop_basicutils
  use coop_string
  use coop_type_arguments
  use coop_type_function
  use coop_type_species
  implicit none
#include "constants.h"

  private

  public:: coop_cosmology, coop_cosmology_background



  type coop_cosmology
     COOP_SHORT_STRING  name
     COOP_INT id
   contains
     procedure :: init => coop_cosmology_initialize
     procedure :: print=> coop_cosmology_print
     procedure :: free => coop_cosmology_free
  end type coop_cosmology

  type, extends(coop_cosmology):: coop_cosmology_background
     COOP_REAL Omega_k
     COOP_INT num_species
     type(coop_species), dimension(coop_max_num_species)::species
   contains
     procedure::add_species=>coop_cosmology_background_add_species
  end type coop_cosmology_background

  interface coop_cosmology
     procedure coop_cosmology_constructor
  end interface coop_cosmology

  interface coop_cosmology_background
     procedure coop_cosmology_background_constructor
  end interface coop_cosmology_background

contains


  subroutine coop_cosmology_free(this)
    class(coop_cosmology):: this
    integer i
    select type (this)
    class is(coop_cosmology_background)
       do i= 1, this%num_species
          call this%species(i)%free
       enddo
       this%num_species = 0
       this%Omega_k = 1.
    end select
  end subroutine coop_cosmology_free


  subroutine coop_cosmology_initialize(this, name, id)
    class(coop_cosmology)::this
    COOP_UNKNOWN_STRING, optional::name
    COOP_INT, optional::id
    call this%free()
    if(present(name))then
       this%name = name
    else
       select type (this)
       type is (coop_cosmology)
          this%name = "COOP_COSMOLOGY"
       type is(coop_cosmology_background)
          this%name = "COOP_COSMOLOGY_BACKGROUND"
          this%Omega_k = 1.
       class default
          this%name = "COOP_COSMOLOGY_UNKNOWN"
       end select
    endif
    if(present(id))then
       this%id =  id
    else
       this%id = 0
    endif
  end subroutine coop_cosmology_initialize


  function coop_cosmology_constructor(name, id) result(this)
    type(coop_cosmology)::this
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
  end function coop_cosmology_constructor


  function coop_cosmology_background_constructor(name, id) result(this)
    type(coop_cosmology_background)::this
    COOP_UNKNOWN_STRING, optional::name
    COOP_INT, optional::id
    this%num_species = 0
    this%Omega_k = 1.
    if(present(name))then
       this%name = name
    else
       this%name = "COOP_COSMOLOGY_BACKGROUND"
    endif
    if(present(id))then
       this%id =  id
    else
       this%id = 1
    endif
  end function coop_cosmology_background_constructor

  subroutine coop_cosmology_background_add_species(this, species)
    class(coop_cosmology_background)::this
    type(coop_species):: species
    if(this%num_species .ge. coop_max_num_species) stop "coop_cosmology_background_add_species: too many species"
    this%num_species = this%num_species+1
    this%species(this%num_species:this%num_species) = species
    this%Omega_k = this%Omega_k - species%Omega
  end subroutine coop_cosmology_background_add_species


  subroutine coop_cosmology_print(this)
    class(coop_cosmology)::this
    integer i
    write(*,"(A)") "================================="
    select type(this)
    type is(coop_cosmology)
       write(*,"(A)") "Cosmology Class: Null"
    type is(coop_cosmology_Background)
       write(*,"(A)") "Cosmology Class: Background"
    class default
       write(*,"(A)") "Cosmology Class: Unknown"
    end select
    write(*,"(A)") "Cosmology Name = "//trim(this%name)
    write(*,"(A)") "Cosmology id = "//trim(coop_num2str(this%id))
    select type(this)
    class is (coop_cosmology_background)
       do i=1, this%num_species
          write(*,"(A)") "---------------------------------"
          write(*,"(A)") "Species #: "//trim(coop_num2str(i))
          call this%species(i)%print
       enddo
       write(*,"(A)") "================================="
    end select
  end subroutine coop_cosmology_print



end module coop_type_cosmology
