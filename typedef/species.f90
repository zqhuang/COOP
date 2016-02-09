module coop_species_mod
  use coop_constants_mod
  use coop_basicutils_mod
  use coop_particle_mod
  use coop_string_mod
  use coop_function_mod
  implicit none
#include "constants.h"


  private
  public:: coop_species, coop_species_constructor


  !!Omega is defined as rho / (3 M_p^2 H^2);
  !!The unit of energy density is H_0^2; of length/time is 1/H_0; of momentum is H_0.
  
  type coop_species
     COOP_INT::genre
     COOP_SHORT_STRING::name = ''
     COOP_INT::id = 0
     COOP_REAL::Omega = 0.
     COOP_REAL::wp1 = 1.d0
     COOP_REAL::cs2 = 0.d0
     COOP_REAL::Omega_massless = 0.d0
     COOP_REAL::mbyT = 0.d0
     COOP_REAL::Mpsq0 = 1.d0
     type(coop_function):: fwp1, fcs2, fwp1eff
     type(coop_function):: flnrho
#if DO_COUPLED_DE
     type(coop_function)::cplde_wp1
     type(coop_function)::cplde_Q
     type(coop_function)::cplde_dQdphi_lna
     type(coop_function)::cplde_lnV_lna
     type(coop_function)::cplde_dVdphibyH2_lna     !!log(- d ln V / d phi) as a function of ln a
     type(coop_function)::cplde_phi_lna
     type(coop_function)::cplde_phi_prime_lna     
     type(coop_function)::cplde_m2byH2_lna
#endif     
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
     procedure :: rhoa3 => coop_species_rhoa3
     procedure :: rhoa4 => coop_species_rhoa4     
     procedure :: pa2 => coop_species_pa2
     procedure :: dpa2da => coop_species_dpa2da
     procedure :: drhoa2da => coop_species_drhoa2da
  end type coop_species
  



contains



!!$  function coop_Mpsq(a) result(Mpsq)
!!$    COOP_REAL::Mpsq, a
!!$#if DO_EFT_DE
!!$    Mpsq = coop_EFT_DE_Mpsq%eval(a)
!!$#else
!!$    Mpsq = 1.d0
!!$#endif    
!!$  end function coop_Mpsq
!!$
!!$  function coop_alphaM(a) result(alphaM)
!!$    COOP_REAL::alphaM, a
!!$#if DO_EFT_DE
!!$    alphaM = coop_EFT_DE_alphaM%eval(a)
!!$#else
!!$    alphaM = 0.d0
!!$#endif        
!!$  end function coop_alphaM
!!$
!!$  function coop_alphaM_prime(a) result(alphaM_prime)
!!$    COOP_REAL::alphaM_prime, a
!!$#if DO_EFT_DE
!!$    alphaM_prime = coop_EFT_DE_alphaM%derivative(a)*a
!!$#else
!!$    alphaM_prime = 0.d0
!!$#endif        
!!$  end function coop_alphaM_prime
  

  function coop_species_constructor(genre, name, id, Omega, w, cs2, Omega_massless, fwp1, fcs2, fwp1eff, Mpsq0) result(this)
    type(coop_species) :: this
    COOP_INT, optional::genre
    COOP_UNKNOWN_STRING, optional::name
    COOP_INT, optional:: id
    COOP_REAL, optional::Omega
    COOP_REAL, optional::w
    COOP_REAL, optional::cs2
    COOP_REAL, optional::Omega_massless
    type(coop_function),optional::fwp1, fcs2, fwp1eff
    COOP_REAL,optional::Mpsq0
    COOP_INT i
    COOP_REAL_ARRAY::lnrat, warr, cs2arr
    COOP_REAL amin, amax, lnamin, lnamax, w1, w2, dlna, lnma, lnm, lnrho, lnp, dlnrho, dlnp
    if(present(Mpsq0))then
       this%Mpsq0 = Mpsq0
    else
       this%Mpsq0 = 1.d0
    endif
    
    if(present(name))then
       this%name = name
    else
       this%name = "COOP_SPECIES"
    endif
    if(present(id))then
       this%id = id
    else
       this%id = 0
    endif

    if(present(Omega))then
       this%Omega = Omega
    else
       this%Omega = 1.d0
    endif

    if(present(Omega_massless))then
       this%Omega_massless = max(min(Omega_massless, this%Omega), COOP_REAL_OF(0.))
    endif

    if(present(genre))then
       this%genre =genre
       if(this%genre .eq. COOP_SPECIES_MASSIVE_FERMION .or. this%genre .eq. COOP_SPECIES_MASSIVE_BOSON)then
          if(this%Omega_massless .ge. this%Omega * 0.9999)then
             this%genre = COOP_SPECIES_MASSLESS
          endif
       endif
    else
       if(present(Omega_massless))then
          if(Omega_massless .le. 0.)then
             this%genre = COOP_SPECIES_CDM
          elseif(Omega_massless/this%Omega .le. 0.9999)then
             this%genre = COOP_SPECIES_MASSIVE_FERMION
          else
             this%genre = COOP_SPECIES_MASSLESS
          endif
       else
          this%genre = COOP_SPECIES_FLUID
       endif
    endif

    select case(this%genre)
    case(COOP_SPECIES_MASSLESS)
       this%wp1 = 1.d0+1.d0/3.
       this%cs2 = 1.d0/3.d0
    case(COOP_SPECIES_CDM)
       this%wp1 = 1.d0
       this%cs2 = 0.d0
    case(COOP_SPECIES_COSMOLOGICAL_CONSTANT)
       this%wp1 = 0.d0
       this%cs2 = 1.d0
    case default
       if(present(w))then
          this%wp1 = w+1.d0
       else
          this%wp1 = 1.d0
       endif
       if(present(cs2))then
          this%cs2 = cs2
       else
          this%cs2 = 0.d0
       endif
    end select
    select case(this%genre)
    case(COOP_SPECIES_MASSLESS, COOP_SPECIES_COSMOLOGICAL_CONSTANT, COOP_SPECIES_CDM)
       return
    case(COOP_SPECIES_MASSIVE_BOSON)
       amin = coop_min_scale_factor
       amax = coop_scale_factor_today
       lnamin = log(amin)
       lnamax = log(amax)
       dlna = (lnamax - lnamin)/(coop_default_array_size  - 1)
       w1 = log(this%Omega/this%Omega_massless)
       call coop_boson_get_lnam(w1, lnm)
       this%mbyT = exp(lnm)
       !$omp parallel do private(lnma, lnrho, lnp, dlnrho, dlnp)
       do i=1, coop_default_array_size
          lnma = lnm + lnamin+(i-1)*dlna
          call coop_boson_get_lnrho(lnma, lnrho)
          call coop_boson_get_lnp(lnma, lnp)
          call coop_boson_get_dlnrho(lnma, dlnrho)
          call coop_boson_get_dlnp(lnma, dlnp)
          warr(i) = exp(lnp-lnrho)
          cs2arr(i) = exp(dlnp-dlnrho)*warr(i)
          lnrat(i) = lnrho- 4.d0*(lnamin+(i-1)*dlna) - w1
       enddo
       !$omp end parallel do
       call this%fwp1%init(coop_default_array_size, amin, amax, 1.d0+warr, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false., name = trim(this%name)//"1+w(a)")
       call this%fcs2%init(coop_default_array_size, amin, amax, cs2arr, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false.,name= trim(this%name)//" c_s^2(a)")
       call this%flnrho%init(coop_default_array_size, amin, amax, lnrat, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false., name= trim(this%name)//"\ln\rho(a) ratio")
    case(COOP_SPECIES_MASSIVE_FERMION)
       amin = coop_min_scale_factor
       amax = coop_scale_factor_today
       lnamin = log(amin)
       lnamax = log(amax)
       dlna = (lnamax - lnamin)/(coop_default_array_size  - 1)
       w1 = log(this%Omega/this%Omega_massless)
       call coop_fermion_get_lnam( w1 , lnm)
       this%mbyT = exp(lnm)
       !$omp parallel do private(lnma, lnrho, lnp, dlnrho, dlnp)
       do i=1, coop_default_array_size
          lnma = lnm + lnamin+(i-1)*dlna
          call coop_fermion_get_lnrho(lnma, lnrho)
          call coop_fermion_get_lnp(lnma, lnp)
          call coop_fermion_get_dlnrho(lnma, dlnrho)
          call coop_fermion_get_dlnp(lnma, dlnp)
          warr(i) = exp(lnp-lnrho)
          cs2arr(i) = exp(dlnp-dlnrho)*warr(i)
          lnrat(i) = lnrho- 4.d0*(lnamin+(i-1)*dlna) - w1
       enddo
       !$omp end parallel do
       call this%fwp1%init(coop_default_array_size, amin, amax, 1.d0+warr, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false., name = trim(this%name)//" 1+w(a)")
       call this%fcs2%init(coop_default_array_size, amin, amax, cs2arr, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false., name = trim(this%name)//" cs^2(a)")
       call this%flnrho%init(coop_default_array_size, amin, amax, lnrat, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false., name = trim(this%name)//" \ln\rho(a) ratio")
    case default
       if(present(fwp1))then
          this%fwp1 = fwp1
       else
          if(present(w))then
             this%wp1 = 1.d0 + w
          else
             stop "You have to pass either w or fwp1 to species_init"
          endif
       endif
       if(present(fwp1eff))  this%fwp1eff = fwp1eff
       if(present(fcs2))this%fcs2 = fcs2       
       if(this%fwp1eff%initialized .or. this%fwp1%initialized)then
          amin = coop_min_scale_factor
          amax = coop_scale_factor_today
          lnamin = log(amin)
          lnamax = log(amax)
          dlna = (lnamax - lnamin)/(coop_default_array_size  - 1)
          w1 = this%wp1effofa(amax)
          lnrat(coop_default_array_size) = 0.d0
          do i= coop_default_array_size - 1, 1, -1
             w2 = this%wp1effofa(exp((lnamin + dlna*(i-1))))
             lnrat(i) = lnrat(i+1) - (w2 + w1 + 4.*this%wp1effofa(exp((lnamin + dlna*(i-0.5)))))
             w1 = w2
          enddo
          lnrat = lnrat*(-dlna/2.)
          call this%flnrho%init(coop_default_array_size, amin, amax, lnrat, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false., name = trim(this%name)//"\ln\rho(a) ratio")
       endif
    end select
    if(this%fwp1%initialized)then
       this%wp1 = this%wp1ofa(coop_scale_factor_today)
       call this%flnrho%set_boundary(slopeleft = -3.*this%wp1effofa(amin), sloperight = -3.*this%wp1effofa(amax))
    endif
    if(this%fcs2%initialized) this%cs2 = this%cs2ofa(coop_scale_factor_today)
  end function coop_species_constructor
  
  
  subroutine coop_species_initialize(this, genre, name, id, Omega, w, cs2, Omega_massless, fwp1, fcs2, fwp1eff, Mpsq0)
    class(coop_species) :: this
    COOP_INT, optional::genre
    COOP_UNKNOWN_STRING, optional::name
    COOP_INT, optional:: id
    COOP_REAL, optional::Omega
    COOP_REAL, optional::w
    COOP_REAL, optional::cs2
    COOP_REAL, optional::Omega_massless
    type(coop_function),optional::fwp1, fcs2, fwp1eff
    COOP_REAL,optional::Mpsq0
    COOP_INT i
    COOP_REAL_ARRAY::lnrat, warr, cs2arr
    COOP_REAL amin, amax, lnamin, lnamax, w1, w2, dlna, lnma, lnm, lnrho, lnp, dlnrho, dlnp
    if(present(Mpsq0))then
       this%Mpsq0 = Mpsq0
    else
       this%Mpsq0 = 1.d0
    endif
    

    if(present(name))then
       this%name = name
    else
       this%name = "COOP_SPECIES"
    endif
    if(present(id))then
       this%id = id
    else
       this%id = 0
    endif

    if(present(Omega))then
       this%Omega = Omega
    else
       this%Omega = 1.d0
    endif

    if(present(Omega_massless))then
       this%Omega_massless = max(min(Omega_massless, this%Omega), COOP_REAL_OF(0.))
    endif

    if(present(genre))then
       this%genre =genre
       if(this%genre .eq. COOP_SPECIES_MASSIVE_FERMION .or. this%genre .eq. COOP_SPECIES_MASSIVE_BOSON)then
          if(this%Omega_massless .ge. this%Omega * 0.9999)then
             this%genre = COOP_SPECIES_MASSLESS
          endif
       endif
    else
       if(present(Omega_massless))then
          if(Omega_massless .le. 0.)then
             this%genre = COOP_SPECIES_CDM
          elseif(Omega_massless/this%Omega .le. 0.9999)then
             this%genre = COOP_SPECIES_MASSIVE_FERMION
          else
             this%genre = COOP_SPECIES_MASSLESS
          endif
       else
          this%genre = COOP_SPECIES_FLUID
       endif
    endif

    select case(this%genre)
    case(COOP_SPECIES_MASSLESS)
       this%wp1 = 1.d0+1.d0/3.
       this%cs2 = 1.d0/3.d0
    case(COOP_SPECIES_CDM)
       this%wp1 = 1.d0
       this%cs2 = 0.d0
    case(COOP_SPECIES_COSMOLOGICAL_CONSTANT)
       this%wp1 = 0.d0
       this%cs2 = 1.d0
    case default
       if(present(w))then
          this%wp1 = w+1.d0
       else
          this%wp1 = 1.d0
       endif
       if(present(cs2))then
          this%cs2 = cs2
       else
          this%cs2 = 0.d0
       endif
    end select
    select case(this%genre)
    case(COOP_SPECIES_MASSLESS, COOP_SPECIES_COSMOLOGICAL_CONSTANT, COOP_SPECIES_CDM)
       return
    case(COOP_SPECIES_MASSIVE_BOSON)
       amin = coop_min_scale_factor
       amax = coop_scale_factor_today
       lnamin = log(amin)
       lnamax = log(amax)
       dlna = (lnamax - lnamin)/(coop_default_array_size  - 1)
       w1 = log(this%Omega/this%Omega_massless)
       call coop_boson_get_lnam(w1, lnm)
       this%mbyT = exp(lnm)
       !$omp parallel do private(lnma, lnrho, lnp, dlnrho, dlnp)
       do i=1, coop_default_array_size
          lnma = lnm + lnamin+(i-1)*dlna
          call coop_boson_get_lnrho(lnma, lnrho)
          call coop_boson_get_lnp(lnma, lnp)
          call coop_boson_get_dlnrho(lnma, dlnrho)
          call coop_boson_get_dlnp(lnma, dlnp)
          warr(i) = exp(lnp-lnrho)
          cs2arr(i) = exp(dlnp-dlnrho)*warr(i)
          lnrat(i) = lnrho- 4.d0*(lnamin+(i-1)*dlna) - w1
       enddo
       !$omp end parallel do
       call this%fwp1%init(coop_default_array_size, amin, amax, 1.d0+warr, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false., name = trim(this%name)//"1+w(a)")
       call this%fcs2%init(coop_default_array_size, amin, amax, cs2arr, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false.,name= trim(this%name)//" c_s^2(a)")
       call this%flnrho%init(coop_default_array_size, amin, amax, lnrat, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false., name= trim(this%name)//"\ln\rho(a) ratio")
    case(COOP_SPECIES_MASSIVE_FERMION)
       amin = coop_min_scale_factor
       amax = coop_scale_factor_today
       lnamin = log(amin)
       lnamax = log(amax)
       dlna = (lnamax - lnamin)/(coop_default_array_size  - 1)
       w1 = log(this%Omega/this%Omega_massless)
       call coop_fermion_get_lnam( w1 , lnm)
       this%mbyT = exp(lnm)
       !$omp parallel do private(lnma, lnrho, lnp, dlnrho, dlnp)
       do i=1, coop_default_array_size
          lnma = lnm + lnamin+(i-1)*dlna
          call coop_fermion_get_lnrho(lnma, lnrho)
          call coop_fermion_get_lnp(lnma, lnp)
          call coop_fermion_get_dlnrho(lnma, dlnrho)
          call coop_fermion_get_dlnp(lnma, dlnp)
          warr(i) = exp(lnp-lnrho)
          cs2arr(i) = exp(dlnp-dlnrho)*warr(i)
          lnrat(i) = lnrho- 4.d0*(lnamin+(i-1)*dlna) - w1
       enddo
       !$omp end parallel do
       call this%fwp1%init(coop_default_array_size, amin, amax, 1.d0+warr, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false., name = trim(this%name)//" 1+w(a)")
       call this%fcs2%init(coop_default_array_size, amin, amax, cs2arr, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false., name = trim(this%name)//" cs^2(a)")
       call this%flnrho%init(coop_default_array_size, amin, amax, lnrat, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false., name = trim(this%name)//" \ln\rho(a) ratio")
    case default
       if(present(fwp1))then
          this%fwp1 = fwp1
       else
          if(present(w))then
             this%wp1 = 1.d0 + w
          else
             stop "You have to pass either w or fwp1 to species_init"
          endif
       endif
       if(present(fwp1eff))  this%fwp1eff = fwp1eff
       if(present(fcs2))this%fcs2 = fcs2       
       if(this%fwp1eff%initialized .or. this%fwp1%initialized)then
          amin = coop_min_scale_factor
          amax = coop_scale_factor_today
          lnamin = log(amin)
          lnamax = log(amax)
          dlna = (lnamax - lnamin)/(coop_default_array_size  - 1)
          w1 = this%wp1effofa(amax)
          lnrat(coop_default_array_size) = 0.d0
          do i= coop_default_array_size - 1, 1, -1
             w2 = this%wp1effofa(exp((lnamin + dlna*(i-1))))
             lnrat(i) = lnrat(i+1) - (w2 + w1 + 4.*this%wp1effofa(exp((lnamin + dlna*(i-0.5)))))
             w1 = w2
          enddo
          lnrat = lnrat*(-dlna/2.)
          call this%flnrho%init(coop_default_array_size, amin, amax, lnrat, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false., name = trim(this%name)//"\ln\rho(a) ratio")
       endif
    end select
    if(this%fwp1%initialized)then
       this%wp1 = this%wp1ofa(coop_scale_factor_today)
       call this%flnrho%set_boundary(slopeleft = -3.*this%wp1effofa(amin), sloperight = -3.*this%wp1effofa(amax))
    endif
    if(this%fcs2%initialized) this%cs2 = this%cs2ofa(coop_scale_factor_today)
  end subroutine coop_species_initialize

  function coop_species_wofa(this, a) result(w)
    class(coop_species) :: this
    COOP_REAL w, a
    w = coop_species_wp1ofa(this, a) - 1.d0
  end function coop_species_wofa


  function coop_species_dwda(this, a) result(dwda)
    class(coop_species) :: this
    COOP_REAL dwda, a
    if(this%fwp1%initialized)then
       dwda = this%fwp1%derivative(COOP_PROPER_SCALE_FACTOR(a))
    else
       dwda = 0.d0
    endif
  end function coop_species_dwda


  function coop_species_wp1ofa(this, a) result(wp1)
    class(coop_species) :: this
    COOP_REAL wp1, a
    if(this%fwp1%initialized)then
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
    if(this%fwp1eff%initialized)then
       wp1 = this%fwp1eff%eval(COOP_PROPER_SCALE_FACTOR(a))
    elseif(this%fwp1%initialized)then
       wp1 = this%fwp1%eval(COOP_PROPER_SCALE_FACTOR(a))
    else
       wp1 = this%wp1
    endif
  end function coop_species_wp1effofa


  function coop_species_cs2ofa(this, a) result(cs2)
    class(coop_species) :: this
    COOP_REAL cs2, a
    if(this%fcs2%initialized)then
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
    if(this%fwp1%initialized)then
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
       write(*,"(A, G14.5)") "Species rho_ratio(z=0) = ", this%density_ratio(COOP_REAL_OF(1.d0))
       
    else
       write(*,"(A, G14.5)") "Species w: ", this%wp1-1.d0
    endif
    if(this%fcs2%initialized)then
       write(*,"(A, G14.5)") "Species cs^2(z=10000) = ", this%cs2ofa(COOP_REAL_OF(1.d0/10001.d0))
       write(*,"(A, G14.5)") "Species cs^2(z=1000) = ", this%cs2ofa(COOP_REAL_OF(1.d0/1001.d0))
       write(*,"(A, G14.5)") "Species cs^2(z=10) = ", this%cs2ofa(COOP_REAL_OF(1.d0/11.d0))
       write(*,"(A, G14.5)") "Species cs^2(z=1) = ", this%cs2ofa(COOP_REAL_OF(0.5d0))
       write(*,"(A, G14.5)") "Species cs^2(z=0) = ", this%cs2ofa(COOP_REAL_OF(1.))
    else
       write(*,"(A, G14.5)") "Species cs^2: ", this%cs2
    endif
  end subroutine coop_species_print

  function coop_species_rhoa4_ratio(this, a) result(rhoa4)
    class(coop_species)::this
    COOP_REAL a, an
    COOP_REAL rhoa4
    if(this%flnrho%initialized)then
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
    if(this%flnrho%initialized)then
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
    if(this%flnrho%initialized)then
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
#if DO_EFT_DE
    density = 3.d0*this%Mpsq0*this%Omega * this%density_ratio(a)    
#else    
    density = 3.d0*this%Omega * this%density_ratio(a)
#endif    
  end function coop_species_density


  function coop_species_rhoa2(this, a) result(density)  !!unit H_0^2M_p^2
    class(coop_species)::this
    COOP_REAL a
    COOP_REAL density
#if DO_EFT_DE
    density = 3.d0*this%Omega *this%Mpsq0 * this%rhoa2_ratio(a)    
#else    
    density = 3.d0*this%Omega * this%rhoa2_ratio(a)
#endif    
  end function coop_species_rhoa2


  function coop_species_rhoa3(this, a) result(density)  !!unit H_0^2M_p^2
    class(coop_species)::this
    COOP_REAL a
    COOP_REAL density
#if DO_EFT_DE
    density = 3.d0*this%Omega * this%Mpsq0*this%rhoa3_ratio(a)    
#else    
    density = 3.d0*this%Omega * this%rhoa3_ratio(a)
#endif    
  end function coop_species_rhoa3


  function coop_species_rhoa4(this, a) result(density)  !!unit H_0^2M_p^2
    class(coop_species)::this
    COOP_REAL a
    COOP_REAL density
#if DO_EFT_DE
    density = 3.d0*this%Omega*this%Mpsq0 * this%rhoa4_ratio(a)    
#else    
    density = 3.d0*this%Omega * this%rhoa4_ratio(a)
#endif    
  end function coop_species_rhoa4
  

  function coop_species_dlnrhodlna(this, a) result(dlnrhodlna)  
    class(coop_species)::this
    COOP_REAL a
    COOP_REAL dlnrhodlna
    
    if(this%flnrho%initialized)then
       if(this%fwp1eff%initialized)then
          dlnrhodlna = -3.d0*this%wp1effofa(a)
       else
          dlnrhodlna = this%flnrho%derivative_bare(log(a))
       endif
    else
       dlnrhodlna = -3.d0*this%wp1
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


  function coop_species_dpa2da(this, a) result(p)  !!unit H_0^2M_p^2
    class(coop_species)::this
    COOP_REAL a
    COOP_REAL p
    p = (this%dwda(a) + this%wofa(a)*(this%dlnrhodlna(a) + 2.d0)/a) * this%rhoa2(a)
  end function coop_species_dpa2da

  function coop_species_density_ratio(this, a) result(density)
    class(coop_species)::this
    COOP_REAL a
    COOP_REAL density
    if(this%flnrho%initialized)then
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
    this%Mpsq0 = 1.d0
#if DO_COUPLED_DE    
    call this%cplde_wp1%free()
    call this%cplde_Q%free()    
    call this%cplde_lnV_lna%free()
    call this%cplde_phi_lna%free()    
    call this%cplde_phi_prime_lna%free()
    call this%cplde_m2byH2_lna%free()
    call this%cplde_dQdphi_lna%free()
    call this%cplde_dVdphibyH2_lna%free()    
#endif    
  end subroutine coop_species_free


end module coop_species_mod
