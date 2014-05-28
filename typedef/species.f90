module coop_type_species
  use coop_constants
  use coop_basicutils
  use coop_type_function
  implicit none
#include "constants.h"

  private

  public:: coop_species

  type coop_species
     logical w_dynamic, cs2_dynamic
     COOP_SHORT_STRING name
     COOP_INT id
     COOP_REAL Omega, w, cs2
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


contains

  subroutine coop_species_initialize(this, name, id, Omega, w, cs2, fw, fcs2)
    class(coop_species) :: this
    COOP_UNKNOWN_STRING, optional::name
    COOP_INT, optional:: id
    COOP_REAL, optional::Omega
    COOP_REAL, optional::w
    COOP_REAL, optional::cs2
    type(coop_function),optional::fw
    type(coop_function),optional::fcs2
    COOP_INT i
    COOP_REAL_ARRAY::lnrat
    COOP_REAL amin, amax, lnamin, lnamax, w1, w2, dlna

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
    if(present(fw))then
       this%w_dynamic = .true.
       allocate(this%fw, source = fw)
    else
       this%w_dynamic = .false.
    endif
    if(present(fcs2))then
       this%cs2_dynamic = .true.
       allocate(this%fcs2, source = fcs2)
    else
       this%cs2_dynamic = .false.
    endif
    if(this%w_dynamic)then
       amin = COOP_min_scale_factor
       amax = 1.
       lnamin = log(amin)
       lnamax = log(amax)
       dlna = (lnamax - lnamin)/(coop_default_array_size  - 1)
       w1 = this%wofa(amax)
       lnrat(coop_default_array_size) = 0.
       do i= coop_default_array_size - 1, 1, -1
          w2 = this%wofa(exp((lnamin + dlna*(i-1))))
          lnrat(i) = lnrat(i+1) - (6.+w2 + w1 + 4.*this%wofa(this%wofa(exp((lnamin + dlna*(i-0.5))))))
          w1 = w2
       enddo
       lnrat = lnrat*(-dlna/2.)
       call this%frho%init(coop_default_array_size, amin, amax, exp(lnrat), method = COOP_INTERPOLATE_SPLINE, xlog = .true., ylog = .true.)
    endif
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
    write(*,"(A, I8)") "Species ID: ", this%id
    write(*,"(A, G15.6)") "Species Omega: ", this%Omega
    if(this%w_dynamic)then
       write(*,"(A, G14.5)") "Species w(a=0.00001) = ", this%wofa(COOP_REAL_OF(0.00001))
       write(*,"(A, G14.5)") "Species w(a=0.01) = ", this%wofa(COOP_REAL_OF(0.01))
       write(*,"(A, G14.5)") "Species w(a=0.1) = ", this%wofa(COOP_REAL_OF(0.1))
       write(*,"(A, G14.5)") "Species w(a=0.5) = ", this%wofa(COOP_REAL_OF(0.5))
       write(*,"(A, G14.5)") "Species w(a=1) = ", this%wofa(COOP_REAL_OF(1.))
       write(*,"(A, G14.5)") "Species rho_de(a=0.00001) = ", this%density_ratio(COOP_REAL_OF(0.00001))
       write(*,"(A, G14.5)") "Species rho_de(a=0.01) = ", this%density_ratio(COOP_REAL_OF(0.01))
       write(*,"(A, G14.5)") "Species rho_de(a=0.1) = ", this%density_ratio(COOP_REAL_OF(0.1))
       write(*,"(A, G14.5)") "Species rho_de(a=0.5) = ", this%density_ratio(COOP_REAL_OF(0.5))
       write(*,"(A, G14.5)") "Species rho_de(a=1) = ", this%density_ratio(COOP_REAL_OF(1.))
       
    else
       write(*,"(A, G14.5)") "Species w: ", this%w
    endif
    if(this%cs2_dynamic)then
       write(*,"(A, G14.5)") "Species cs^2(a=0.00001) = ", this%cs2ofa(COOP_REAL_OF(0.00001))
       write(*,"(A, G14.5)") "Species cs^2(a=0.01) = ", this%cs2ofa(COOP_REAL_OF(0.01))
       write(*,"(A, G14.5)") "Species cs^2(a=0.1) = ", this%cs2ofa(COOP_REAL_OF(0.1))
       write(*,"(A, G14.5)") "Species cs^2(a=0.5) = ", this%cs2ofa(COOP_REAL_OF(0.5))
       write(*,"(A, G14.5)") "Species cs^2(a=1) = ", this%cs2ofa(COOP_REAL_OF(1.))
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
