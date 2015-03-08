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
     COOP_INT::genre
     COOP_SHORT_STRING::name = ''
     COOP_INT::id = 0
     COOP_REAL::Omega = 0.
     COOP_REAL::wp1 = 1.d0
     COOP_REAL::cs2 = 0.d0
     COOP_REAL::Omega_massless = 0.d0
     COOP_REAL::mbyT = 0.d0
     type(coop_function):: fwp1, fcs2, fwp1eff
     type(coop_function):: flnrho
     !!for scalar field DE model
     type(coop_function)::fDE_U_of_phi, fDE_dUdphi, fDE_phi, fDE_dphidlna, fDE_Q_of_phi, fDE_phidot
     COOP_REAL::DE_tracking_n = 0.d0
     COOP_REAL::DE_a_start = 1.d-6
     COOP_REAL::DE_lnV0 = 0.d0
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
     procedure :: dpa2da => coop_species_dpa2da
     procedure :: drhoa2da => coop_species_drhoa2da
     !!for scalar field DE model
     procedure :: DE_V => coop_species_DE_V  !!V(phi)
     procedure :: DE_dlnVdphi => coop_species_DE_dlnVdphi  !!d ln V/d phi  (phi)
     procedure :: DE_d2lnVdphi2 => coop_species_DE_d2lnVdphi2     !!d^2ln V/d phi^2 (phi)
     procedure :: DE_get_VVpVpp => coop_species_DE_get_VVpVpp     !!get V, dV/dphi, d^2V/dphi^2
     procedure :: DE_phi => coop_species_DE_phi !!phi(a)
     procedure :: DE_phidot => coop_species_DE_phidot  !!phidot (a)
     procedure :: DE_Q => coop_species_DE_Q !!Q(a)
     procedure :: DE_dlnQdphi => coop_species_DE_dlnQdphi   !!d ln Q/d phi (phi)
  end type coop_species
  



contains

  subroutine coop_species_initialize(this, genre, name, id, Omega, w, cs2, Omega_massless, fwp1, fcs2, fwp1eff)
    class(coop_species) :: this
    COOP_INT, optional::genre
    COOP_UNKNOWN_STRING, optional::name
    COOP_INT, optional:: id
    COOP_REAL, optional::Omega
    COOP_REAL, optional::w
    COOP_REAL, optional::cs2
    COOP_REAL, optional::Omega_massless
    type(coop_function),optional::fwp1, fcs2, fwp1eff
    COOP_INT i
    COOP_REAL_ARRAY::lnrat, warr, cs2arr
    COOP_REAL amin, amax, lnamin, lnamax, w1, w2, dlna, lnma, lnm, lnrho, lnp, dlnrho, dlnp
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
    call this%fDE_U_of_phi%free()
    call this%fDE_dUdphi%free()
    call this%fDE_phi%free()
    call this%fDE_phidot%free()    
    call this%fDE_Q_of_phi%free()
  end subroutine coop_species_free


  function coop_species_DE_phi(this, a) result(phi)
    class(coop_species)::this
    COOP_REAL a, phi
    phi = this%fDE_phi%eval(a)
  end function coop_species_DE_phi

  
  function coop_species_DE_Q(this, phi) result(Q)
    class(coop_species)::this
    COOP_REAL Q, phi
    if(this%fDE_Q_of_phi%initialized)then    
       Q = this%fDE_Q_of_phi%eval(phi)
    else
       Q = 0.d0
    endif
  end function coop_species_DE_Q

  function coop_species_DE_dlnQdphi(this, phi) result(dlnQdphi)
    class(coop_species)::this
    COOP_REAL dlnQdphi, phi
    if(this%fDE_Q_of_phi%initialized)then    
       dlnQdphi = this%fDE_Q_of_phi%derivative(phi)/this%fDE_Q_of_phi%eval(phi)
    else
       dlnQdphi = 0.d0
    endif
  end function coop_species_DE_dlnQdphi


  function coop_species_DE_V(this, phi) result(V)
    class(coop_species)::this
    COOP_REAL V, phi
    V = exp(this%fDE_U_of_phi%eval(phi) + this%DE_lnV0 - this%DE_tracking_n * COOP_LN(phi))
  end function coop_species_DE_V

  function coop_species_DE_phidot(this, a) result(dphidlna)
    class(coop_species)::this
    COOP_REAL::a, dphidlna
    dphidlna = this%fDE_phidot%eval(a)
  end function coop_species_DE_phidot

  function coop_species_DE_dlnVdphi(this, phi) result(dlnVdphi)
    class(coop_species)::this
    COOP_REAL:: phi, dlnVdphi, V
    dlnVdphi = this%fDE_dUdphi%eval(phi) - this%DE_tracking_n/max(phi, tiny(phi))
  end function coop_species_DE_dlnVdphi

  function coop_species_DE_d2lnVdphi2(this, phi) result(d2lnVdphi2)
    class(coop_species)::this
    COOP_REAL::d2lnVdphi2, phi
    d2lnVdphi2 = this%fDE_dUdphi%derivative(phi) + this%DE_tracking_n/max(phi, tiny(phi))**2
  end function coop_species_DE_d2lnVdphi2

  subroutine coop_species_DE_get_VVpVpp(this, phi, V, Vp, Vpp)
    class(coop_species)::this
    COOP_REAL::phi, V, Vp, Vpp
    V = this%DE_V(phi)
    Vp = this%DE_dlnVdphi(phi)
    Vpp = this%DE_d2lnVdphi2(phi)
    Vpp = (Vpp + Vp**2)*V
    Vp = Vp * V
  end subroutine coop_species_DE_get_VVpVpp

end module coop_species_mod
