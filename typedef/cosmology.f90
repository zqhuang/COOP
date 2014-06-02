module coop_type_cosmology
  use coop_constants
  use coop_basicutils
  use coop_string
  use coop_particle
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
     private
     COOP_REAL Omega_k_value, h_value, Tcmb_value, YHe_value, Nnu_value
     COOP_INT num_species
     type(coop_species), dimension(coop_max_num_species)::species
     logical::need_setup_background
     type(coop_function):: fdis, ftime
   contains
     procedure::h => coop_cosmology_background_hubble  !!H_0  = 100 h km/s/Mpc
     procedure::Tcmb => coop_cosmology_background_Tcmb !!CMB temprature in Kelvin
     procedure::Tnu => coop_cosmology_background_Tnu !!CMB temprature in Kelvin
     procedure::YHe => coop_cosmology_background_YHe  !!Helium mass fraction
     procedure::Nnu => coop_cosmology_background_Nnu  !!Helium mass fraction

     procedure::Omega_k => coop_cosmology_background_Omega_k 
     procedure::AgeGyr => coop_cosmology_background_AgeGyr !!return age in Gyr
     procedure::Omega_radiation => coop_cosmology_background_Omega_radiation
     procedure::Omega_massless_neutrinos_per_species => coop_cosmology_background_Omega_massless_neutrinos_per_species
     procedure::Omega_massless_neutrinos => coop_cosmology_background_Omega_massless_neutrinos
     procedure::Omega_nu_from_mnu_eV => coop_cosmology_background_Omega_nu_from_mnu_eV
     procedure::mnu_eV_from_Omega_nu => coop_cosmology_background_mnu_eV_from_Omega_nu
     procedure::set_h => coop_cosmology_background_set_hubble
     procedure::set_Tcmb => coop_cosmology_background_set_Tcmb
     procedure::set_YHe => coop_cosmology_background_set_YHe
     procedure::set_Nnu => coop_cosmology_background_set_Nnu

     procedure::setup_background => coop_cosmology_background_setup_background
     procedure::H0Mpc => coop_cosmology_background_H0Mpc  !!H_0 * Mpc
     procedure::H0Gyr => coop_cosmology_background_H0Gyr  !!H_0 * Gyr
     procedure::Hasq => coop_cosmology_background_Hasq   !! input a , return H a^2 / H_0
     procedure::dadtau => coop_cosmology_background_Hasq   !! input a , return da/d tau /H_0
     procedure::Hratio => coop_cosmology_background_Hratio !! input a, return H/H_0
     procedure::time => coop_cosmology_background_time !!input a, return H_0 * time
     procedure::conformal_time => coop_cosmology_background_conformal_time !!input a, return H_0 * conformal time
     !!input  a1, a2 (optional, default = 1)
     !!output H_0 * distance(a1, a2)  
     !!you need to divide the output by this%H0Mpc() to get the distance in unit of Mpc
     procedure::comoving_distance => coop_cosmology_background_comoving_distance
     procedure::comoving_angular_diameter_distance => coop_cosmology_background_comoving_angular_diameter_distance  
     procedure::luminosity_distance => coop_cosmology_background_luminosity_distance   
     procedure::angular_diameter_distance => coop_cosmology_background_angular_diameter_distance 
     procedure::add_species=>coop_cosmology_background_add_species
  end type coop_cosmology_background

  interface coop_cosmology
     procedure coop_cosmology_constructor
  end interface coop_cosmology

  interface coop_cosmology_background
     procedure coop_cosmology_background_constructor
  end interface coop_cosmology_background

contains

  function coop_cosmology_background_hubble(this) result(h)
    class(coop_cosmology_background)::this
    COOP_REAL h
    h = this%h_value
  end function coop_cosmology_background_hubble

  function coop_cosmology_background_Omega_k(this) result(omega_k)
    class(coop_cosmology_background)::this
    COOP_REAL omega_k
    Omega_k = this%omega_k_value
  end function coop_cosmology_background_Omega_k

  function coop_cosmology_background_YHe(this) result(YHe)
    class(coop_cosmology_background)::this
    COOP_REAL YHe
    YHe = this%YHe_value
  end function coop_cosmology_background_YHe

  function coop_cosmology_background_Nnu(this) result(Nnu)
    class(coop_cosmology_background)::this
    COOP_REAL Nnu
    Nnu = this%Nnu_value
  end function coop_cosmology_background_Nnu

  function coop_cosmology_background_Tcmb(this, a) result(Tcmb)
    class(coop_cosmology_background)::this
    COOP_REAL, optional::a
    COOP_REAL Tcmb
    Tcmb = this%Tcmb_value
    if(present(a))Tcmb = Tcmb/a
  end function coop_cosmology_background_Tcmb

  subroutine coop_cosmology_background_set_hubble(this, h)
    class(coop_cosmology_background)::this
    COOP_REAL h
    this%h_value = h
  end subroutine coop_cosmology_background_set_hubble


  subroutine coop_cosmology_background_set_Tcmb(this, Tcmb)
    class(coop_cosmology_background)::this
    COOP_REAL Tcmb
    this%Tcmb_value = Tcmb
  end subroutine coop_cosmology_background_set_Tcmb

  subroutine coop_cosmology_background_set_YHe(this, YHe)
    class(coop_cosmology_background)::this
    COOP_REAL YHe
    this%YHe_value = YHe
  end subroutine coop_cosmology_background_set_YHe


  subroutine coop_cosmology_background_set_NNu(this, Nnu)
    class(coop_cosmology_background)::this
    COOP_REAL Nnu
    this%Nnu_value = Nnu
  end subroutine coop_cosmology_background_set_Nnu

  subroutine coop_cosmology_initialize(this, name, id, h, Tcmb, YHe, Nnu)
    class(coop_cosmology)::this
    COOP_UNKNOWN_STRING, optional::name
    COOP_INT, optional::id
    COOP_REAL, optional:: h 
    COOP_REAL, optional:: Tcmb
    COOP_REAL, optional:: YHe
    COOP_REAL, optional:: Nnu
    select type(this)
    type is (coop_cosmology)
#include "cosmology_init.h"
    type is (coop_cosmology_background)
#include "cosmology_background_init.h"
    class default
       Stop "coop_cosmology_initialize: Unknown class of cosmology."
    end select
  end subroutine coop_cosmology_initialize

  subroutine coop_cosmology_free(this)
    class(coop_cosmology):: this
    integer i
    select type (this)
    class is(coop_cosmology_background)
       call this%fdis%free
       call this%ftime%free
       do i= 1, this%num_species
          call this%species(i)%free
       enddo
       this%num_species = 0
       this%omega_k_value = 1.
       this%need_setup_background = .true.
    end select
  end subroutine coop_cosmology_free

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
       write(*,"(A)") "Omega_k =  "//trim(coop_num2str(this%omega_k_value))
       write(*,"(A)") "Hubble (100 km/s/Mpc) =  "//trim(coop_num2str(this%h_value))
       write(*,"(A)") "CMB temperature (K) =  "//trim(coop_num2str(this%Tcmb_value))
       write(*,"(A)") "Helium mass fraction =  "//trim(coop_num2str(this%YHe_value))
       write(*,"(A)") "number of neutrino species=  "//trim(coop_num2str(this%Nnu_value))
       do i=1, this%num_species
          write(*,"(A)") "---------------------------------"
          write(*,"(A)") "Species #: "//trim(coop_num2str(i))
          call this%species(i)%print
       enddo
       write(*,"(A)") "================================="
    end select
  end subroutine coop_cosmology_print



  function coop_cosmology_constructor(name, id) result(this)
    type(coop_cosmology)::this
    COOP_UNKNOWN_STRING, optional::name
    COOP_INT, optional::id
#include "cosmology_init.h"
  end function coop_cosmology_constructor


  function coop_cosmology_background_constructor(name, id, h, Tcmb, YHe, Nnu) result(this)
    type(coop_cosmology_background)::this
    COOP_UNKNOWN_STRING, optional::name
    COOP_INT, optional::id
    COOP_REAL, optional:: h 
    COOP_REAL, optional:: Tcmb
    COOP_REAL, optional:: YHe
    COOP_REAL, optional:: Nnu
#include "cosmology_background_init.h"
  end function coop_cosmology_background_constructor

  subroutine coop_cosmology_background_add_species(this, species)
    class(coop_cosmology_background)::this
    type(coop_species):: species
    if(this%num_species .ge. coop_max_num_species) stop "coop_cosmology_background_add_species: too many species"
    this%num_species = this%num_species+1
    this%species(this%num_species:this%num_species) = species
    this%omega_k_value = this%omega_k_value - species%Omega
    this%need_setup_background = .true.
  end subroutine coop_cosmology_background_add_species

  function coop_cosmology_background_H0Mpc(this)result(H0Mpc)
    class(coop_cosmology_background)::this
    COOP_REAL H0Mpc
    H0Mpc = this%h_value * (1.d5/coop_SI_c)
  end function coop_cosmology_background_H0Mpc

  function coop_cosmology_background_H0Gyr(this)result(H0Gyr)
    class(coop_cosmology_background)::this
    COOP_REAL H0Gyr
    H0Gyr = this%h_value * (1.d5*coop_SI_Gyr/coop_SI_Mpc)
  end function coop_cosmology_background_H0Gyr

  function coop_cosmology_background_Hasq(this, a)result(Hasq)
    class(coop_cosmology_background)::this
    COOP_REAL Hasq, a
    integer i
    Hasq =  this%Omega_k_value*a**2
    do i=1, this%num_species
       Hasq = Hasq + this%species(i)%Omega * this%species(i)%rhoa4_ratio(a)
    enddo
    Hasq = sqrt(Hasq)
  end function coop_cosmology_background_Hasq

  function coop_cosmology_background_Hratio(this, a)result(Hratio)
    class(coop_cosmology_background)::this
    COOP_REAL Hratio, a
    Hratio = this%Hasq(a)/a**2
  end function coop_cosmology_background_Hratio

  subroutine coop_cosmology_background_setup_background(this)
    class(coop_cosmology_background)::this
    integer,parameter::n = 16384
    COOP_REAL, parameter::amin = 0.d0, amax = coop_scale_factor_today
    COOP_REAL,dimension(n):: dis, a, t
    COOP_REAL  hasq1, hasq2, da, daby2
    integer i
    if(this%need_setup_background)then
       call coop_set_uniform(n, a, amin, amax)
       dis(1) = 0.d0
       t(1) = 0.d0
       da = (amax-amin)/(n-1.d0)
       daby2 = da/2.d0
       hasq1 = this%Hasq(amin)
       do i=2, n
          hasq2 = this%Hasq(a(i))
          dis(i) = dis(i-1) + (1.d0/Hasq1 + 1.d0/hasq2 + 4.d0/this%hasq(a(i)-daby2))
          t(i) = t(i-1) + (a(i-1)/Hasq1  + a(i)/hasq2 + 4.d0*(a(i)-daby2)/this%hasq(a(i)-daby2))
          Hasq1 = hasq2
       enddo
       dis = dis * (da/6.d0)
       t = t * (da/6.d0)
       call this%fdis%init(n, amin, amax, dis, COOP_INTERPOLATE_LINEAR)
       call this%ftime%init(n, amin, amax, t, COOP_INTERPOLATE_LINEAR)
       this%need_setup_background = .false. 
    endif
  end subroutine coop_cosmology_background_setup_background

  function coop_cosmology_background_conformal_time(this, a) result(tau)
    class(coop_cosmology_background)::this
    COOP_REAL a, tau
    if(this%need_setup_background) call this%setup_background()
    tau = this%fdis%eval(a)
  end function coop_cosmology_background_conformal_time

  function coop_cosmology_background_time(this, a) result(t)
    class(coop_cosmology_background)::this
    COOP_REAL a, t
    if(this%need_setup_background) call this%setup_background()
    t = this%ftime%eval(a)
  end function coop_cosmology_background_time

  function coop_cosmology_background_AgeGyr(this) result(AgeGyr)
    class(coop_cosmology_background)::this
    COOP_REAL AgeGyr 
    AgeGyr = this%time(coop_scale_factor_today)/this%H0Gyr()
  end function coop_cosmology_background_AgeGyr

  function coop_cosmology_background_comoving_distance(this, a1, a2) result(chi)
    class(coop_cosmology_background)::this
    COOP_REAL a1, chi
    COOP_REAL,optional::a2
    if(this%need_setup_background) call this%setup_background()
    if(present(a2))then
       chi = abs(this%fdis%eval(a2)  - this%fdis%eval(a1))
    else
       chi = abs(this%fdis%eval(coop_scale_factor_today)  - this%fdis%eval(a1))
    endif
  end function coop_cosmology_background_comoving_distance
  
  !!r, chi must be in unit of H_0^{-1}
  function coop_r_of_chi(chi, Omega_k) result(r)
    COOP_REAL chi, Omega_k, r, oc2
    oc2 = Omega_k*chi**2  
    if(abs(oc2).lt. 0.2d0)then
       r = chi * (1.d0 + oc2*(1.d0/6.d0 + oc2*(1.d0/120.d0+oc2*(1.d0/5040.d0+oc2/362880.d0))))
    elseif(oc2 .gt. 0.d0)then
       oc2 = sqrt(oc2)
       r = chi*(sinh(oc2)/oc2)
    else
       oc2 = sqrt(-oc2)
       r = chi*(sin(oc2)/oc2)
    endif
  end function coop_r_of_chi


  function coop_cosmology_background_comoving_angular_diameter_distance(this, a1, a2) result(r)
    class(coop_cosmology_background)::this
    COOP_REAL a1, r
    COOP_REAL,optional::a2
    if(present(a2))then
       r =coop_r_of_chi(abs(this%fdis%eval(a2)  - this%fdis%eval(a1)), this%Omega_k_value)
    else
       r = coop_r_of_chi(abs(this%fdis%eval(coop_scale_factor_today)  - this%fdis%eval(a1)), this%Omega_k_value)
    endif
  end function coop_cosmology_background_comoving_angular_diameter_distance


  function coop_cosmology_background_angular_diameter_distance(this, a) result(dA)
    class(coop_cosmology_background)::this
    COOP_REAL a, dA
    dA = coop_r_of_chi(abs(this%fdis%eval(coop_scale_factor_today)  - this%fdis%eval(a)), this%Omega_k_value)*a
  end function coop_cosmology_background_angular_diameter_distance

  function coop_cosmology_background_luminosity_distance(this, a) result(dlum)
    class(coop_cosmology_background)::this
    COOP_REAL a, dlum
    dlum = coop_r_of_chi(abs(this%fdis%eval(coop_scale_factor_today)  - this%fdis%eval(a)), this%Omega_k_value)/a
  end function coop_cosmology_background_luminosity_distance

  function coop_cosmology_background_Omega_radiation(this) result(Omega_r)
    class(coop_cosmology_background)::this
    COOP_REAL, parameter:: c =  coop_SI_blackbody_alpha*2.d0/coop_SI_rhocritbyh2
    COOP_REAL Omega_r
    Omega_r = c * (this%Tcmb_value **2 / this%h_value)**2
  end function coop_cosmology_background_Omega_radiation

  function coop_cosmology_background_Omega_massless_neutrinos_per_species(this) result(Omega_nu_massless)
    class(coop_cosmology_background)::this
    COOP_REAL, parameter:: c =  coop_SI_blackbody_alpha*2.d0/coop_SI_rhocritbyh2 * (7./8.) * (4.d0/11.)**(4.d0/3.d0)*coop_neutrinos_temperature_correction**4
    COOP_REAL Omega_nu_massless
    Omega_nu_massless =  c * (this%Tcmb_value **2 / this%h_value)**2
  end function coop_cosmology_background_Omega_massless_neutrinos_per_species


  function coop_cosmology_background_Omega_massless_neutrinos(this) result(Omega_nu_massless)
    class(coop_cosmology_background)::this
    COOP_REAL, parameter:: c =  coop_SI_blackbody_alpha*2.d0/coop_SI_rhocritbyh2 * (7./8.) * (4.d0/11.)**(4.d0/3.d0)*coop_neutrinos_temperature_correction**4
    COOP_REAL Omega_nu_massless
    Omega_nu_massless =  c * (this%Tcmb_value **2 / this%h_value)**2 * this%Nnu_value
  end function coop_cosmology_background_Omega_massless_neutrinos

  function coop_cosmology_background_Tnu(this, a) result(Tnu)
    class(coop_cosmology_background)::this
    COOP_REAL, parameter::c = (4.d0/11.)**(1.d0/3.d0)*coop_neutrinos_temperature_correction
    COOP_REAL, optional::a
    COOP_REAL Tnu
    Tnu = this%Tcmb_value * c
    if(present(a))Tnu = Tnu/a
  end function coop_cosmology_background_Tnu

  function coop_cosmology_background_Omega_nu_from_mnu_eV(this, mnu_eV) result(Omega_nu)
    class(coop_cosmology_background)::this
    COOP_REAL:: mnu_eV, Omega_nu, lnrho
    call coop_fermion_get_lnrho(log(mnu_eV/this%Tnu() * (coop_SI_eV/coop_SI_kB)), lnrho)
    Omega_nu = this%Omega_massless_neutrinos() * exp(lnrho)
  end function coop_cosmology_background_Omega_nu_from_mnu_eV

  function coop_cosmology_background_mnu_eV_from_Omega_nu(this, Omega_nu) result(mnu_eV)
    class(coop_cosmology_background)::this
    COOP_REAL:: mnu_eV, Omega_nu, lnam
    call coop_fermion_get_lnam(log(Omega_nu/this%Omega_massless_neutrinos()), lnam)
    mnu_eV = this%Tnu()* (coop_SI_kB/coop_SI_eV)*exp(lnam)

  end function coop_cosmology_background_mnu_eV_from_Omega_nu



end module coop_type_cosmology