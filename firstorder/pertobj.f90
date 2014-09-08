module coop_pertobj_mod
  use coop_wrapper_background
  implicit none
#include "constants.h"

  COOP_INT, parameter:: coop_pert_default_lmax = 32  
  COOP_INT, parameter:: coop_pert_default_mmax = 2
  COOP_INT, parameter:: coop_pert_default_smax = 2


  COOP_INT, parameter:: coop_pert_default_nq = 5
  COOP_REAL, dimension(coop_pert_default_nq),parameter:: coop_pert_default_q = coop_fermion_int_q5
  COOP_REAL, dimension(coop_pert_default_nq),parameter:: coop_pert_default_q_kernel = coop_fermion_int_kernel5 


  type coop_pert_species
     COOP_INT::genre = COOP_PERT_NONE  !!NONE = no perturbations, other options are PERFECT_FLUID, HIERARCHY
     COOP_INT::m = 0
     COOP_INT::s = 0
     COOP_INT::index = -1  !!index of this variable in the ODE array
     COOP_INT::lmin = 0
     COOP_INT::lmax = -1
     COOP_REAL,dimension(0:coop_pert_default_lmax)::F = 0.d0
     COOP_REAL::q, mass
  end type coop_pert_species

  type coop_pert_object
     logical::tight_coupling = .true.
     logical::has_massivenu = .false.
     logical::has_de = .false.
     COOP_INT::m = 0
     COOP_INT::ny = 0
     type(coop_pert_species)::baryon, cdm, T, E, B, nu, de
     type(coop_pert_species),dimension(coop_pert_default_nq)::massivenu !!massive neutrinos

     COOP_REAL,dimension(:),allocatable::y
   contains
     procedure::init =>  coop_pert_object_initialize
     procedure::free =>  coop_pert_object_free
  end type coop_pert_object


contains

  subroutine coop_pert_object_initialize(this, m,  nu_mass, de_genre)
    !!nu_mass = neutrino mass to temprature ratio
    class(coop_pert_object)::this
    COOP_REAL,optional::nu_mass
    COOP_INT, optional::de_genre
    COOP_INT::m, i, iq
    this%m = m

    !!Cold Dark Matter
    this%cdm%genre = COOP_PERT_PERFECT_FLUID
    this%cdm%m = m
    this%cdm%s = 0
    this%cdm%index = 1
    this%cdm%lmin = max(this%cdm%m, this%cdm%s)
    this%cdm%lmax = max(this%cdm%lmin - 1, 1)
    this%cdm%F = 0.d0
    this%cdm%q = 1.d-30
    this%cdm%mass = 1.d30
    if(this%cdm%lmax .lt. this%cdm%lmin) this%cdm%genre = COOP_PERT_NONE

    !!baryon
    this%baryon%genre = COOP_PERT_PERFECT_FLUID
    this%baryon%m = m
    this%baryon%s = 0
    this%baryon%index = this%cdm%index + this%cdm%lmax - this%cdm%lmin + 1
    this%baryon%lmin = max(this%baryon%m, this%baryon%s)
    this%baryon%lmax = max(this%baryon%lmin - 1, 1)
    this%baryon%F = 0.d0
    this%baryon%q = 1.d-30
    this%baryon%mass = 1.d30
    if(this%baryon%lmax .lt. this%baryon%lmin) this%baryon%genre = COOP_PERT_NONE

    !!T
    this%T%genre = COOP_PERT_HIERARCHY
    this%T%m = m
    this%T%s = 0
    this%T%index = this%baryon%index + this%baryon%lmax - this%baryon%lmin + 1
    this%T%lmin = max(this%T%m, this%T%s)
    this%T%lmax = max(this%T%lmin - 1, 16)
    this%T%F = 0.d0
    this%T%q = 1.d0
    this%T%mass = 0.d0

    !!E
    this%E%genre = COOP_PERT_HIERARCHY
    this%E%m = m
    this%E%s = 2
    this%E%index = this%T%index + this%T%lmax - this%T%lmin + 1
    this%E%lmin = max(this%E%m, this%E%s)
    this%E%lmax = max(this%E%lmin - 1, 12)
    this%E%F = 0.d0
    this%E%q = 1.d0
    this%E%mass = 0.d0
    
    !!B
    this%B%index = this%E%index + this%E%lmax - this%E%lmin + 1
    if(m.eq.0)then
       this%B%genre = COOP_PERT_NONE   
       this%B%lmin = 0
       this%B%lmax = -1
    else
       this%B%genre = COOP_PERT_HIERARCHY
       this%B%m = m
       this%B%s = 2
       this%B%lmin = max(this%B%m, this%B%s)
       this%B%lmax = max(this%B%lmin - 1, 12)
       this%B%F = 0.d0
       this%B%q = 1.d0
       this%B%mass = 0.d0
    endif

    !!Neutrinos
    this%nu%genre = COOP_PERT_HIERARCHY
    this%nu%m = m
    this%nu%s = 0
    this%nu%index = this%B%index + this%B%lmax - this%B%lmin + 1
    this%nu%lmin = max(this%nu%m, this%nu%s)
    this%nu%lmax = max(this%nu%lmin - 1, 12)
    this%nu%F = 0.d0
    this%nu%q = 1.d0
    this%nu%mass = 0.d0

    if(present(nu_mass))then
       this%has_massivenu = (nu_mass .ne. 0)
    else
       this%has_massivenu = .false.
    endif
    if(this%has_massivenu)then
       this%massivenu%genre = COOP_PERT_HIERARCHY
       this%massivenu%m = m
       this%massivenu%s = 0
       this%massivenu%lmin = max(this%nu%m, this%nu%s)
       this%massivenu%lmax = max(this%nu%lmin - 1, 12)
       this%massivenu%mass = nu_mass
       do iq = 1, coop_pert_default_nq
          this%massivenu(iq)%F = 0.d0
          this%massivenu(iq)%q = coop_pert_default_q(iq)
          if(iq .eq. 1)then
             this%massivenu(iq)%index = this%nu%index + this%nu%lmax - this%nu%lmin + 1
          else
             this%massivenu(iq)%index = this%massivenu(iq-1)%index + this%massivenu(iq-1)%lmax - this%massivenu(iq-1)%lmin + 1
          endif
       enddo
    else
       this%massivenu%genre = COOP_PERT_NONE
       this%massivenu%lmin = 0
       this%massivenu%lmax = -1
       this%massivenu%index = this%nu%index + this%nu%lmax - this%nu%lmin + 1
    endif

    !!massive neutrinos


    !!dark energy
    if(present(de_genre))then
       this%de%genre = de_genre
       this%de%m = m
       this%de%s = 0
       this%de%index = this%massivenu(coop_pert_default_nq)%index + this%massivenu(coop_pert_default_nq)%lmax - this%massivenu(coop_pert_default_nq)%lmin + 1
       this%de%lmin = max(this%de%m, this%de%s)
       select case(this%de%genre)
       case(COOP_PERT_NONE)
          this%de%lmax = max(this%de%lmin - 1, -1)
       case(COOP_PERT_PERFECT_FLUID)
          this%de%lmax = max(this%de%lmin - 1, 1)
       case(COOP_PERT_SCALAR_FIELD) 
          if(m.eq.0)then
             this%de%lmax = 1
          else
             this%de%lmax = this%de%lmin - 1
          endif
       case default
          call coop_return_error("pert_initialize", "unknown genre "//trim(coop_num2str(this%de%genre)), "stop")
       end select
       if(this%de%lmax .lt. this%de%lmin) this%de%genre = COOP_PERT_NONE
    endif
   
  end subroutine coop_pert_object_initialize


  subroutine coop_pert_object_free(this)
    class(coop_pert_object)::this
    if(allocated(this%y))deallocate(this%y)
    this%ny = 0
  end subroutine coop_pert_object_free


  subroutine coop_pert_object_set_ode_index(this)
    class(coop_pert_object)::this
    
  end subroutine coop_pert_object_set_ode_index
  
  subroutine coop_pert_object_read_var(this, var)
    class(coop_pert_object)::this
    COOP_REAL,dimension(:),optional::var

  end subroutine coop_pert_object_read_var

  subroutine coop_pert_object_write_var(this, var)
    class(coop_pert_object)::this
    COOP_REAL,dimension(:),optional::var

  end subroutine coop_pert_object_write_var
  


end module coop_pertobj_mod
