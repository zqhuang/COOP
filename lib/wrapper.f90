module coop_wrapper
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"

  type(coop_cosmology_background):: coop_global_cosmology
  type(coop_species):: coop_global_baryon
  type(coop_species):: coop_global_cdm
  type(coop_species):: coop_global_radiation
  type(coop_species):: coop_global_massless_neutrinos
  type(coop_species):: coop_global_massive_neutrinos
  type(coop_species):: coop_global_de
  type(coop_arguments)::coop_global_cosmological_parameters
  COOP_REAL::coop_scalar_pivot_k = 0.05d0
  COOP_REAL::coop_tensor_pivot_k = 0.002d0

  
contains

  subroutine coop_print_info()
    write(*,*) "This COOP Version "//trim(coop_version)
    write(*,*) "Author: Zhiqi Huang"
  end subroutine coop_print_info

  subroutine coop_setup_global_cosmology_with_h(h)
    COOP_REAL h
    call COOP_COSMO%init(name = "COOP_GLOBAL_COSMOLOGY",  id = 0, h = h)
    if(h.le.0.d0)return  !!return for bad h
    coop_global_baryon = coop_baryon(COOP_OMEGABH2/h**2)
    call COOP_COSMO%add_species(coop_global_baryon)
    coop_global_cdm = coop_cdm(COOP_OMEGACH2/h**2)
    call COOP_COSMO%add_species(coop_global_cdm)
    coop_global_radiation = coop_radiation(COOP_COSMO%Omega_radiation())
    call COOP_COSMO%add_species(coop_global_radiation)
    if(COOP_MNU .eq. 0.d0)then
       coop_global_massless_neutrinos = coop_neutrinos_massless(COOP_COSMO%Omega_massless_neutrinos())
       call COOP_COSMO%add_species(coop_global_massless_neutrinos)
    else
       coop_global_massless_neutrinos = coop_neutrinos_massless(COOP_COSMO%Omega_massless_neutrinos_per_species()*(COOP_COSMO%NNu()-1))
       call COOP_COSMO%add_species(coop_global_massless_neutrinos)  !!assuming one species
       coop_global_massive_neutrinos = coop_neutrinos_massive( &
            COOP_COSMO%Omega_nu_per_species_from_mnu_eV(COOP_MNU) ,&
            COOP_COSMO%Omega_massless_neutrinos_per_species())
       call COOP_COSMO%add_species(coop_global_massive_neutrinos)
    endif
    select case(COOP_DE_MODEL)
    case(COOP_DE_COSMOLOGICAL_CONSTANT)
       coop_global_de = coop_de_lambda( &
            COOP_COSMO%Omega_k()-COOP_OMEGAK &
            )
    case(COOP_DE_W0)
       coop_global_de = coop_de_w0( &
            COOP_COSMO%Omega_k()-COOP_OMEGAK, &
            COOP_COSMO_PARAMS%r(COOP_INDEX_DE) &
            )
    case(COOP_DE_W0WA)
       coop_global_de = coop_de_w0wa( &
            COOP_COSMO%Omega_k()-COOP_OMEGAK, &
            COOP_COSMO_PARAMS%r(COOP_INDEX_DE), &
            COOP_COSMO_PARAMS%r(COOP_INDEX_DE+1) &
            )
    case(COOP_DE_QUINTESSENCE)
      coop_global_de = coop_de_quintessence( &
            COOP_COSMO%Omega_k()-COOP_OMEGAK, &
            COOP_COSMO_PARAMS%r(COOP_INDEX_DE), &
            COOP_COSMO_PARAMS%r(COOP_INDEX_DE+1), &
            COOP_COSMO_PARAMS%r(COOP_INDEX_DE+2) &
            )
    case(COOP_DE_COUPLED_QUINTESSENCE)
       coop_global_de = coop_de_coupled_quintessence( &
            COOP_COSMO%Omega_k() - COOP_OMEGAK, &
            COOP_COSMO_PARAMS%r(COOP_INDEX_DE), &
            COOP_COSMO_PARAMS%r(COOP_INDEX_DE+1), &
            COOP_COSMO_PARAMS%r(COOP_INDEX_DE+2), &
            COOP_COSMO_PARAMS%r(COOP_INDEX_DE+3) &
            )
    case default
       stop "UNKNOWN DARK ENERGY MODEL"
    end select
    call COOP_COSMO%add_species(coop_global_de)
  end subroutine coop_setup_global_cosmology_with_h

  subroutine coop_setup_global_cosmology()
    COOP_REAL,parameter::hmin = 0.5, hmax = 0.9
    COOP_REAL hl, hr, hm, tl, tr, tm
    hl = hmin
    call coop_setup_global_cosmology_with_h(hl)
    tl = 100.d0*coop_global_cosmology_cosmomc_theta()
    hr = hmax
    call coop_setup_global_cosmology_with_h(hr)
    tr = 100.d0*coop_global_cosmology_cosmomc_theta()
    if(tl .lt. COOP_100THETA-1.d-8 .and. tr .gt. COOP_100THETA+1.d-8)then
       do while(hr - hl .gt. 0.0001)
          hm = (hl + hr)/2.
          call coop_setup_global_cosmology_with_h(hm)
          tm = 100.d0*coop_global_cosmology_cosmomc_theta()
          if(tm .gt. COOP_100THETA+1.d-8)then
             hr  = hm
          elseif(tm .lt. COOP_100THETA - 1.d-8)then
             hl = hm
          else
             return
          endif
       enddo
    else
       write(*,*) hl, tl
       write(*,*) hr, tr
       write(*,*) "Bad cosmological parameters. Setting h = 0."
       call coop_setup_global_cosmology_with_h(COOP_REAL_OF(0.))
    endif
  end subroutine coop_setup_global_cosmology


  function coop_global_cosmology_cosmomc_theta() result(theta)
    COOP_REAL theta, zstar, astar, rs, da
    zstar = coop_zrecomb_fitting(COOP_OMEGABH2, COOP_OMEGACH2)
    astar = 1.d0/(1.d0 + zstar)
    rs = coop_integrate( coop_global_cosmology_dsoundda, COOP_REAL_OF(0.d0), astar, COOP_REAL_OF(1.d-7) )
    da = coop_r_of_chi(coop_integrate( coop_global_cosmology_dtauda, astar, COOP_REAL_OF(1.d0), COOP_REAL_OF(1.d-7) ), COOP_OMEGAK)
    theta = rs/da
  end function coop_global_cosmology_cosmomc_theta

  function coop_global_cosmology_dsoundda(a) result(dsda)
    COOP_REAL a, dsda
    dsda = (1.d0/coop_sqrt3) / sqrt(1.d0 + 0.75d0/ COOP_COSMO%Omega_radiation()/COOP_COSMO%h()**2 * a * COOP_OMEGABH2 ) / COOP_COSMO%Hasq(a)
!    dsda = (1.d0/coop_sqrt3) / sqrt(1.d0 + 3.d4 * a * COOP_OMEGABH2 ) / COOP_COSMO%Hasq(a)  !!note that old cosmomc uses this
  end function coop_global_cosmology_dsoundda

  function coop_global_cosmology_dtauda(a) result(dtauda)
    COOP_REAL a, dtauda
    dtauda = 1.d0/COOP_COSMO%Hasq(a)
  end function coop_global_cosmology_dtauda


  function coop_zrecomb_fitting(ombh2, omch2) result(zstar)
    COOP_REAL zstar, ombh2, omch2
  !!From Hu & Sugiyama
    zstar =  1048 * (1 + 0.00124 * ombh2**(-0.738))*(1+ &
         (0.0783 * ombh2 **(-0.238)/(1+39.5* ombh2 ** 0.763)) * &
         (omch2 + ombh2)**(0.560/(1+21.1* ombh2 **1.81)))
  end function coop_zrecomb_fitting


  function coop_primordial_ps(k)  result(ps)
    COOP_REAL k, ps, lnk
    select case(COOP_PP_MODEL)
    case(COOP_PP_STANDARD)
       lnk = log(k/coop_scalar_pivot_k)
       ps = exp(COOP_LN10TO10AS + (COOP_NS - 1.d0)*lnk + (COOP_NRUN/2.d0)*lnk**2)
    case(COOP_PP_SPLINE)

    case(COOP_PP_LINEAR)

    case default

    end select
  end function coop_primordial_ps


  function coop_primordial_pt(k)  result(pt)
    COOP_REAL k, pt
    select case(COOP_PP_MODEL)
    end select    
  end function coop_primordial_pt


end module coop_wrapper
