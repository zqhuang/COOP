module coop_cosmology_mod
  use coop_constants_mod
  use coop_basicutils_mod
  use coop_string_mod
  use coop_particle_mod
  use coop_arguments_mod
  use coop_function_mod
  use coop_species_mod

  implicit none

#include "constants.h"

  private

  public:: coop_cosmology, coop_cosmology_background, coop_r_of_chi, coop_zrecomb_fitting, coop_cosmology_constructor, coop_cosmology_background_constructor

  type coop_cosmology
     COOP_SHORT_STRING  name
     COOP_INT id
   contains
     procedure :: init => coop_cosmology_initialize
     procedure :: print=> coop_cosmology_print
     procedure :: free => coop_cosmology_free
  end type coop_cosmology

  type, extends(coop_cosmology):: coop_cosmology_background
     COOP_REAL::Omega_k_value = 1.
     COOP_REAL::h_value =  COOP_DEFAULT_HUBBLE
     COOP_REAL::Tcmb_value = COOP_DEFAULT_TCMB
     COOP_REAL::YHe_value= COOP_DEFAULT_YHE
     COOP_REAL:: Nnu_value = COOP_DEFAULT_NUM_NEUTRINO_SPECIES
     logical::need_setup_background = .true.
     type(coop_function):: fdis, ftime, faoftau
#if DO_EFT_DE
     type(coop_function):: f_alpha_M, f_alpha_T, f_alpha_B, f_alpha_K, f_alpha_H, f_M2
#endif     
     COOP_REAL::dis_const, time_const, a_eq, Omega_r, Omega_m, a_switch, dis_switch
     COOP_INT :: num_species = 0
     type(coop_species), dimension(coop_max_num_species)::species
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
     procedure::Omega_nu_per_species_from_mnu_eV => coop_cosmology_background_Omega_nu_per_species_from_mnu_eV
     procedure::set_h => coop_cosmology_background_set_hubble
     procedure::set_Tcmb => coop_cosmology_background_set_Tcmb
     procedure::set_YHe => coop_cosmology_background_set_YHe
     procedure::set_Nnu => coop_cosmology_background_set_Nnu

     procedure::setup_background => coop_cosmology_background_setup_background
     procedure::H0Mpc => coop_cosmology_background_H0Mpc  !!H_0 * Mpc
     procedure::H0Gyr => coop_cosmology_background_H0Gyr  !!H_0 * Gyr
     
     procedure::H2a4 => coop_cosmology_background_H2a4   !! input a , return H a^2 / H_0
     procedure::rhoa4 => coop_cosmology_background_rhoa4 !!input a, return rho a^4 / (H_0^2 M_p^2)

     procedure::get_pa4_rhoa4 => coop_cosmology_background_get_pa4_rhoa4 !!input a, return rho a^4 / (H_0^2 M_p^2)
     
     procedure::Hasq => coop_cosmology_background_Hasq   !! input a , return H a^2 / H_0
     procedure::dadtau => coop_cosmology_background_Hasq   !! input a , return da/d tau /H_0
     procedure::Hratio => coop_cosmology_background_Hratio !! input a, return H/H_0
     procedure::HdotbyHsq => coop_cosmology_background_HdotbyHsq !! input a, return \dot H/H^2
     procedure::HddbyH3 => coop_cosmology_background_HddbyH3     !!input a, return \ddot H/ H^3
     procedure::time => coop_cosmology_background_time !!input a, return H_0 * time
     procedure::conformal_time => coop_cosmology_background_conformal_time !!input a, return H_0 * conformal time

     procedure::tauofa => coop_cosmology_background_conformal_time !!input a, return H_0 * conformal time

     procedure::aoftau => coop_cosmology_background_aoftau !!input  H_0 * conformal time, return a

     !!input  a1, a2 (optional, default = 1)
     !!output H_0 * distance(a1, a2)  
     !!you need to divide the output by this%H0Mpc() to get the distance in unit of Mpc
     procedure::comoving_distance => coop_cosmology_background_comoving_distance
     procedure::comoving_angular_diameter_distance => coop_cosmology_background_comoving_angular_diameter_distance  
     procedure::luminosity_distance => coop_cosmology_background_luminosity_distance   
     procedure::angular_diameter_distance => coop_cosmology_background_angular_diameter_distance
     procedure::add_species=>coop_cosmology_background_add_species
     procedure::delete_species=>coop_cosmology_background_delete_species     
     procedure::index_of => coop_cosmology_background_index_of

     !!functions of z
     procedure::Dv_of_z => coop_cosmology_background_Dv_of_z
     procedure::H_of_z => coop_cosmology_background_H_of_z
     procedure::dA_of_z => coop_cosmology_background_DA_of_z
     procedure::comoving_dA_of_z => coop_cosmology_background_comoving_DA_of_z
     procedure::dL_of_z => coop_cosmology_background_DL_of_z
#if DO_EFT_DE
     procedure::alpha_M => coop_cosmology_background_alpha_M
     procedure::alpha_T => coop_cosmology_background_alpha_T
     procedure::alpha_K => coop_cosmology_background_alpha_K
     procedure::alpha_B => coop_cosmology_background_alpha_B
     procedure::alpha_H => coop_cosmology_background_alpha_H
     procedure::alpha_M_prime => coop_cosmology_background_alpha_M_prime
     procedure::alpha_T_prime => coop_cosmology_background_alpha_T_prime
     procedure::alpha_K_prime => coop_cosmology_background_alpha_K_prime
     procedure::alpha_B_prime => coop_cosmology_background_alpha_B_prime
     procedure::alpha_H_prime => coop_cosmology_background_alpha_H_prime     
     procedure::M2 => coop_cosmology_background_M2     
#endif     
  end type coop_cosmology_background

!!$  interface coop_cosmology
!!$     module procedure coop_cosmology_constructor
!!$  end interface coop_cosmology
!!$
!!$  interface coop_cosmology_background
!!$     module procedure coop_cosmology_background_constructor
!!$  end interface coop_cosmology_background

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
    class is (coop_cosmology_background)
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
       call this%faoftau%free
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

  !!delete the last n species
  subroutine coop_cosmology_background_delete_species(this, ind)
    class(coop_cosmology_background)::this
    COOP_INT:: i, ind
    if(ind .le. 0 .or. ind .gt. this%num_species)return
    this%omega_k_value = this%omega_k_value + this%species(ind)%Omega
    do i = ind, this%num_species-1
       this%species(i) = this%species(i+1)
    enddo
    call this%species(this%num_species)%free()
    this%num_species = this%num_species - 1
    this%need_setup_background = .true.
  end subroutine coop_cosmology_background_delete_species
  

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


  function coop_cosmology_background_H2a4(this, a)result(H2a4)
    class(coop_cosmology_background)::this
    COOP_REAL H2a4, a
    COOP_INT i
    H2a4 =  this%Omega_k_value*a**2
    do i=1, this%num_species
       H2a4 = H2a4 + this%species(i)%Omega * this%species(i)%rhoa4_ratio(a)
    enddo
#if DO_EFT_DE    
    if(this%f_M2%initialized) H2a4 = H2a4/this%f_M2%eval(a)
#endif    
  end function coop_cosmology_background_H2a4


  function coop_cosmology_background_rhoa4(this, a)result(rhoa4)
    class(coop_cosmology_background)::this
    COOP_REAL rhoa4, a
    COOP_INT i
    rhoa4 = 0.d0
    do i=1, this%num_species
       rhoa4 = rhoa4 + this%species(i)%Omega * this%species(i)%rhoa4_ratio(a)
    enddo
    rhoa4 = rhoa4*3.d0
  end function coop_cosmology_background_rhoa4


  subroutine coop_cosmology_background_get_pa4_rhoa4(this, a, pa4, rhoa4)
    class(coop_cosmology_background)::this
    COOP_REAL pa4, rhoa4, a, r
    COOP_INT i
    rhoa4 = 0.d0
    pa4 = 0.d0
    do i=1, this%num_species
       r = this%species(i)%Omega * this%species(i)%rhoa4_ratio(a)
       pa4 = pa4 + r*this%species(i)%wofa(a)
       rhoa4 = rhoa4 +r
    enddo
    rhoa4 = rhoa4 * 3.d0
    pa4 = pa4 * 3.d0
  end subroutine coop_cosmology_background_get_pa4_rhoa4
  

  
  function coop_cosmology_background_Hasq(this, a)result(Hasq)
    class(coop_cosmology_background)::this
    COOP_REAL Hasq, a
    Hasq =  sqrt(this%H2a4(a))
  end function coop_cosmology_background_Hasq

  function coop_cosmology_background_Hratio(this, a)result(Hratio)
    class(coop_cosmology_background)::this
    COOP_REAL Hratio, a
    Hratio = this%Hasq(a)/a**2
  end function coop_cosmology_background_Hratio


  function coop_cosmology_background_HddbyH3(this, a) result(HddbyH3) 
    class(coop_cosmology_background)::this
    COOP_INT i
    COOP_REAL HddbyH3, a, pa4, rhoa4, w,  rhoa4_total, pa4_total, ddterm, H2term, addbyaHsq
    rhoa4_total = 0.d0
    pa4_total = 0.d0
    ddterm = 0.d0
    do i=1, this%num_species
       rhoa4 = this%species(i)%Omega * this%species(i)%rhoa4_ratio(a)
       w = this%species(i)%wofa(a)
       pa4 = w *rhoa4
       ddterm = ddterm + rhoa4*(this%species(i)%wp1effofa(a)*(1.d0+3.d0*w) - this%species(i)%dwda(a) * a)
       rhoa4_total = rhoa4_total + rhoa4
       pa4_total = pa4_total + pa4
    enddo
    H2term = (rhoa4_total + this%Omega_k_value*a**2)
    addbyaHsq =  -(rhoa4_total + 3.d0*pa4_total)/2.d0/H2term
#if DO_EFT_DE
    HddbyH3 = 1.5d0 * ddterm / H2term - 2.d0 * (addbyaHsq - 1.d0) - this%alpha_M(a) *addbyaHsq
#else    
    HddbyH3 = 1.5d0 * ddterm / H2term - 2.d0 * (addbyaHsq - 1.d0)
#endif    
  end function coop_cosmology_background_HddbyH3


  function coop_cosmology_background_HdotbyHsq(this, a)result(HdotbyHsq)
    class(coop_cosmology_background)::this
    COOP_INT i
    COOP_REAL HdotbyHsq, a, pa4, rhoa4,  rhoa4_total, pa4_total
    rhoa4_total = 0.d0
    pa4_total = 0.d0
    do i=1, this%num_species
       rhoa4 = this%species(i)%Omega * this%species(i)%rhoa4_ratio(a)
       pa4 = this%species(i)%wofa(a)*rhoa4
       rhoa4_total = rhoa4_total + rhoa4
       pa4_total = pa4_total + pa4
    enddo
    HdotbyHsq = -(rhoa4_total + 3.d0*pa4_total)/2.d0/(rhoa4_total + this%Omega_k_value*a**2) - 1.d0
  end function coop_cosmology_background_HdotbyHsq


  subroutine coop_cosmology_background_setup_background(this)
    class(coop_cosmology_background)::this
    integer,parameter::n = 32768
    COOP_REAL, parameter::amin = 5.d-4, amax = coop_scale_factor_today
    COOP_REAL,dimension(n):: dis, a, t
    COOP_REAL  hasq1, hasq2, da, daby2, Hasqmin, Hasqmax, M2
    integer i
    if(this%need_setup_background)then
       call coop_set_uniform(n, a, amin, amax)
       this%a_switch = amin
       this%Omega_r = this%H2a4(coop_min_scale_factor)
       this%Omega_m = ( this%H2a4(amin) - this%Omega_r ) / amin
       this%Omega_r = this%H2a4(coop_min_scale_factor) - this%Omega_m * coop_min_scale_factor
       this%Omega_m = ( this%H2a4(amin) - this%Omega_r ) / amin
       this%a_eq = this%Omega_r / this%Omega_m
#if DO_EFT_DE
       M2 = this%M2(this%a_eq)
       this%dis_const = 2.d0*this%a_eq/sqrt(this%Omega_r/M2)
       this%time_const = 4.d0/3.d0*this%a_eq**2/sqrt(this%Omega_r/M2)       
#else       
       this%dis_const = 2.d0*this%a_eq/sqrt(this%Omega_r)
       this%time_const = 4.d0/3.d0*this%a_eq**2/sqrt(this%Omega_r)       
#endif       

       this%dis_switch = this%dis_const * (sqrt(1.d0+a(1)/this%a_eq) - 1.d0)
       da = (amax-amin)/(n-1.d0)
       daby2 = da/2.d0
       dis(1) = this%dis_switch*(6.d0/da)
       t(1) = this%time_const *(1.d0 -  (1.d0- a(1)/2.d0/this%a_eq)*sqrt(1.d0+a(1)/this%a_eq)) *(6.d0/da)
       Hasqmin = this%Hasq(amin)
       Hasqmax = this%Hasq(amax)
       hasq1 = Hasqmin
       do i=2, n
          hasq2 = this%Hasq(a(i))
          dis(i) = dis(i-1) + (1.d0/Hasq1 + 1.d0/hasq2 + 4.d0/this%hasq(a(i)-daby2))
          t(i) = t(i-1) + (a(i-1)/Hasq1  + a(i)/hasq2 + 4.d0*(a(i)-daby2)/this%hasq(a(i)-daby2))
          Hasq1 = hasq2
       enddo
       dis = dis * (da/6.d0)
       t = t * (da/6.d0)
       call this%fdis%init(n, amin, amax, dis, COOP_INTERPOLATE_QUADRATIC, check_boundary = .false., slopeleft = 1.d0/Hasqmin, sloperight = 1.d0/Hasqmax, name = "comoving distance")
       call this%faoftau%init_NonUniform(dis, a, name="a(tau)")
       call this%faoftau%set_boundary(fleft = amin, fright=amax, slopeleft = Hasqmin, sloperight = Hasqmax)
       call this%ftime%init(n, amin, amax, t, COOP_INTERPOLATE_QUADRATIC, check_boundary = .false., name = "t(a)")
       this%need_setup_background = .false. 
    endif
  end subroutine coop_cosmology_background_setup_background

  function coop_cosmology_background_conformal_time(this, a) result(tau)
    class(coop_cosmology_background)::this
    COOP_REAL a, tau, eps
    if(this%need_setup_background) call this%setup_background()
    if(a .gt. this%a_switch)then
       tau = this%fdis%eval(a)
    else
       eps = a/this%a_eq
       if(eps .lt. 1.d-4)then
          tau = this%dis_const * (0.5d0 - eps/8.d0)*eps
       else
          tau = this%dis_const * ( sqrt(1.d0 + eps) - 1.d0)
       endif
    endif
  end function coop_cosmology_background_conformal_time

  function coop_cosmology_background_aoftau(this, tau) result(a)
    class(coop_cosmology_background)::this
    COOP_REAL a, tau, eps
    if(this%need_setup_background) call this%setup_background()
    if(tau .gt. this%dis_switch)then
       a = this%faoftau%eval(tau)
    else
       eps = tau/this%dis_const
       a = this%a_eq * (2.d0+eps)*eps
    endif
  end function coop_cosmology_background_aoftau

  function coop_cosmology_background_time(this, a) result(t)
    class(coop_cosmology_background)::this
    COOP_REAL a, t, eps
    if(this%need_setup_background) call this%setup_background()
    if(a .gt. this%a_switch)then
       t = this%ftime%eval(a)
    else
       eps = a/this%a_eq
       if(eps .lt. 1.d-3)then
          t = this%time_const * (3.d0/8.d0 - eps/8.d0)*eps**2
       else
          t = this%time_const * ( 1.d0 - (1.d0 - eps/2.d0)*sqrt(1.d0+eps))
       endif
    endif
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

  function coop_cosmology_background_Omega_nu_from_mnu_eV(this, mnu_eV) result(Omega_nu)  !!assuming equal mass for each species
    class(coop_cosmology_background)::this
    COOP_REAL:: mnu_eV, Omega_nu, lnrho
    call coop_fermion_get_lnrho(log(mnu_eV/this%Tnu() * (coop_SI_eV/coop_SI_kB)), lnrho)
    Omega_nu = this%Omega_massless_neutrinos() * exp(lnrho)
  end function coop_cosmology_background_Omega_nu_from_mnu_eV

  function coop_cosmology_background_mnu_eV_from_Omega_nu(this, Omega_nu) result(mnu_eV)  !!assuming equal mass for each species
    class(coop_cosmology_background)::this
    COOP_REAL:: mnu_eV, Omega_nu, lnam
    call coop_fermion_get_lnam(log(Omega_nu/this%Omega_massless_neutrinos()), lnam)
    mnu_eV = this%Tnu()* (coop_SI_kB/coop_SI_eV)*(exp(lnam) - 1.d0)

  end function coop_cosmology_background_mnu_eV_from_Omega_nu


  function coop_cosmology_background_Omega_nu_per_species_from_mnu_eV(this, mnu_eV) result(Omega_nu)  
    class(coop_cosmology_background)::this
    COOP_REAL:: mnu_eV, Omega_nu, lnrho
    call coop_fermion_get_lnrho(log(mnu_eV/this%Tnu() * (coop_SI_eV/coop_SI_kB)), lnrho)
    Omega_nu = this%Omega_massless_neutrinos_per_species() * exp(lnrho)
  end function coop_cosmology_background_Omega_nu_per_species_from_mnu_eV


  function coop_cosmology_background_index_of(this, name, error_not_found) result(ind)
    class(coop_cosmology_background)::this
    COOP_UNKNOWN_STRING::name
    logical,optional::error_not_found
    COOP_INT ind, i
    ind = 0
    do i = 1, this%num_species
       if(trim(this%species(i)%name).eq.trim(name))then
          ind = i
          return
       endif
    enddo
    if(present(error_not_found))then
       if(error_not_found)call coop_return_error("index_of", trim(name)//" not found", "stop")
    endif
  end function coop_cosmology_background_index_of

  !!useful fitting functions
  function coop_zrecomb_fitting(ombh2, omch2) result(zstar)
    COOP_REAL zstar, ombh2, omch2
    !!From Hu & Sugiyama
    zstar =  1048 * (1 + 0.00124 * ombh2**(-0.738))*(1+ &
         (0.0783 * ombh2 **(-0.238)/(1+39.5* ombh2 ** 0.763)) * &
         (omch2 + ombh2)**(0.560/(1+21.1* ombh2 **1.81)))
  end function coop_zrecomb_fitting
  

  function coop_cosmology_background_Dv_of_z(this, z) result(Dv)
    class(coop_cosmology_background)::this    
    COOP_REAL::z, Dv, a
    a = 1.d0/(1.d0+z)
    Dv = (z*this%comoving_angular_diameter_distance(a)**2/this%HRatio(a))**(1.d0/3.d0)
  end function coop_cosmology_background_Dv_of_z

  function coop_cosmology_background_H_of_z(this, z) result(Hz)
    class(coop_cosmology_background)::this
    COOP_REAL::z, Hz
    Hz = this%Hratio(1.d0/(1.d0+z))
  end function coop_cosmology_background_H_of_z

  function coop_cosmology_background_dA_of_z(this, z) result(dA)
    class(coop_cosmology_background)::this
    COOP_REAL::z, dA
    dA = this%angular_diameter_distance(1.d0/(1.d0+z))
  end function coop_cosmology_background_dA_of_z

  function coop_cosmology_background_dL_of_z(this, z) result(dL)
    class(coop_cosmology_background)::this
    COOP_REAL::z, dL
    dL = this%luminosity_distance(1.d0/(1.d0+z))
  end function coop_cosmology_background_dL_of_z

  function coop_cosmology_background_comoving_dA_of_z(this, z) result(dA)
    class(coop_cosmology_background)::this
    COOP_REAL::z, dA
    dA = this%comoving_angular_diameter_distance(1.d0/(1.d0+z))
  end function coop_cosmology_background_comoving_dA_of_z

#if DO_EFT_DE
  function coop_cosmology_background_M2(this,a) result(M2)
    class(coop_cosmology_background)::this
    COOP_REAL::M2,a
    if(this%f_M2%initialized)then
       M2 = this%f_M2%eval(a)
    else
       M2 = 1.d0
    endif
  end function coop_cosmology_background_M2
  
  function coop_cosmology_background_alpha_M(this,a) result(alpha_M)
    class(coop_cosmology_background)::this
    COOP_REAL::alpha_M,a
    if(this%f_alpha_M%initialized)then
       alpha_M = this%f_alpha_M%eval(a)
    else
       alpha_M = 0.d0
    endif
  end function coop_cosmology_background_alpha_M

  function coop_cosmology_background_alpha_M_prime(this,a) result(alpha_M_prime)
    class(coop_cosmology_background)::this
    COOP_REAL::alpha_M_prime,a
    if(this%f_alpha_M%initialized)then
       alpha_M_prime = this%f_alpha_M%derivative(a)*a
    else
       alpha_M_prime = 0.d0
    endif
  end function coop_cosmology_background_alpha_M_prime

  function coop_cosmology_background_alpha_T(this,a) result(alpha_T)
    class(coop_cosmology_background)::this
    COOP_REAL::alpha_T,a
    if(this%f_alpha_T%initialized)then
       alpha_T = this%f_alpha_T%eval(a)
    else
       alpha_T = 0.d0
    endif
  end function coop_cosmology_background_alpha_T

  function coop_cosmology_background_alpha_T_prime(this,a) result(alpha_T_prime)
    class(coop_cosmology_background)::this
    COOP_REAL::alpha_T_prime,a
    if(this%f_alpha_T%initialized)then
       alpha_T_prime = this%f_alpha_T%derivative(a)*a
    else
       alpha_T_prime = 0.d0
    endif
  end function coop_cosmology_background_alpha_T_prime


  function coop_cosmology_background_alpha_B(this,a) result(alpha_B)
    class(coop_cosmology_background)::this
    COOP_REAL::alpha_B,a
    if(this%f_alpha_B%initialized)then
       alpha_B = this%f_alpha_B%eval(a)
    else
       alpha_B = 0.d0
    endif
  end function coop_cosmology_background_alpha_B

  function coop_cosmology_background_alpha_B_prime(this,a) result(alpha_B_prime)
    class(coop_cosmology_background)::this
    COOP_REAL::alpha_B_prime,a
    if(this%f_alpha_B%initialized)then
       alpha_B_prime = this%f_alpha_B%derivative(a)*a
    else
       alpha_B_prime = 0.d0
    endif
  end function coop_cosmology_background_alpha_B_prime



  function coop_cosmology_background_alpha_H(this,a) result(alpha_H)
    class(coop_cosmology_background)::this
    COOP_REAL::alpha_H,a
    if(this%f_alpha_H%initialized)then
       alpha_H = this%f_alpha_H%eval(a)
    else
       alpha_H = 0.d0
    endif
  end function coop_cosmology_background_alpha_H

  function coop_cosmology_background_alpha_H_prime(this,a) result(alpha_H_prime)
    class(coop_cosmology_background)::this
    COOP_REAL::alpha_H_prime,a
    if(this%f_alpha_H%initialized)then
       alpha_H_prime = this%f_alpha_H%derivative(a)*a
    else
       alpha_H_prime = 0.d0
    endif
  end function coop_cosmology_background_alpha_H_prime

  
  function coop_cosmology_background_alpha_K(this,a) result(alpha_K)
    class(coop_cosmology_background)::this
    COOP_REAL::alpha_K,a
    if(this%f_alpha_K%initialized)then
       alpha_K = this%f_alpha_K%eval(a)
    else
       alpha_K = 0.d0
    endif
  end function coop_cosmology_background_alpha_K

  function coop_cosmology_background_alpha_K_prime(this,a) result(alpha_K_prime)
    class(coop_cosmology_background)::this
    COOP_REAL::alpha_K_prime,a
    if(this%f_alpha_K%initialized)then
       alpha_K_prime = this%f_alpha_K%derivative(a)*a
    else
       alpha_K_prime = 0.d0
    endif
  end function coop_cosmology_background_alpha_K_prime
#endif  

  
end module coop_cosmology_mod
