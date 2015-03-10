module coop_wrapper
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"

  type(coop_cosmology_firstorder):: coop_global_cosmology
  type(coop_species):: coop_global_baryon
  type(coop_species):: coop_global_cdm
  type(coop_species):: coop_global_radiation
  type(coop_species):: coop_global_massless_neutrinos
  type(coop_species):: coop_global_massive_neutrinos
  type(coop_species):: coop_global_de

  type(coop_arguments)::coop_global_cosmological_parameters, coop_global_param_index

  COOP_REAL,parameter::coop_pp_lnkmin = -9.22d0
  COOP_REAL::coop_pp_lnkmax = -0.3d0



  COOP_INT:: cosmomc_de_index = 8
  COOP_INT:: cosmomc_de2pp_num_params = 7

  COOP_INT:: cosmomc_de_model = COOP_DE_COSMOLOGICAL_CONSTANT
  COOP_INT:: cosmomc_pp_model = COOP_PP_STANDARD
  COOP_INT:: cosmomc_de_num_params = 2
  COOP_INT:: cosmomc_pp_num_params = 8
  COOP_INT:: cosmomc_pp_num_origin = 8
  COOP_INT::  coop_pp_nleft, coop_pp_nright
  COOP_REAL:: coop_pp_lnk_per_knot
  COOP_REAL:: coop_pp_scalar_lnkpivot = log(0.05d0)
  COOP_REAL:: coop_pp_tensor_lnkpivot =  log(0.05d0)
  logical ::cosmomc_pp_inflation_consistency = .true.
  logical ::coop_global_cosmology_do_firstorder = .false.
  logical ::coop_global_cosmology_do_background = .false.
  COOP_INT, parameter::coop_pp_lmin = 2
  COOP_INT, parameter::coop_pp_lmax = 2850
  COOP_REAL::coop_pp_ells(coop_pp_lmin:coop_pp_lmax)  
  COOP_REAL::coop_pp_scalar_Cls(coop_num_cls, coop_pp_lmin:coop_pp_lmax)
  COOP_REAL::coop_pp_lensed_Cls(coop_num_cls, coop_pp_lmin:coop_pp_lmax)
  COOP_REAL::coop_pp_tensor_Cls(coop_num_cls, coop_pp_lmin:coop_pp_lmax)
  COOP_REAL::coop_pp_total_Cls(coop_num_cls, coop_pp_lmin:coop_pp_lmax)
  COOP_INT,parameter::coop_pp_n = 1024  
  COOP_REAL,dimension(coop_pp_n)::coop_pp_lnkMpc, coop_pp_lnps, coop_pp_lnpt, coop_pp_lnps2, coop_pp_lnpt2, coop_pp_lneps, coop_pp_lnV, coop_pp_phi, coop_pp_lnH
  COOP_INT::coop_pp_ipivot
 
  interface coop_setup_cosmology_from_cosmomc
     module procedure coop_setup_cosmology_from_cosmomc_s, coop_setup_cosmology_from_cosmomc_d
  end interface coop_setup_cosmology_from_cosmomc

contains

  subroutine coop_setup_cosmology_from_cosmomc_s(params, h)
    real params(:)
    double precision, optional::h
    call COOP_COSMO_PARAMS%init(r = COOP_REAL_OF(params), i = (/ cosmomc_de_model, cosmomc_de_index, cosmomc_de_num_params, cosmomc_pp_model, cosmomc_de_index + cosmomc_de_num_params + cosmomc_de2pp_num_params, cosmomc_pp_num_params /), l = (/ cosmomc_pp_inflation_consistency /)  )

    if(COOP_INFLATION_CONSISTENCY)then
       COOP_NT = - COOP_AMP_RATIO / 8.d0
    endif

    if(present(h))then
       call coop_setup_global_cosmology_with_h(COOP_REAL_OF(h))
    endif
  end subroutine coop_setup_cosmology_from_cosmomc_s

  subroutine coop_setup_cosmology_from_cosmomc_d(params, h)
    doubleprecision params(:)
    double precision, optional::h
    call COOP_COSMO_PARAMS%init(r = COOP_REAL_OF(params), i = (/ cosmomc_de_model, cosmomc_de_index, cosmomc_de_num_params, cosmomc_pp_model, cosmomc_de_index + cosmomc_de_num_params + cosmomc_de2pp_num_params, cosmomc_pp_num_params /), l = (/ cosmomc_pp_inflation_consistency /) )
    if(COOP_INFLATION_CONSISTENCY)then
       COOP_NT = - COOP_AMP_RATIO / 8.d0
    endif
    if(present(h))then
       call coop_setup_global_cosmology_with_h(COOP_REAL_OF(h))
    endif
  end subroutine coop_setup_cosmology_from_cosmomc_d


  subroutine coop_print_info()
    write(*,*) "This COOP Version "//trim(coop_version)
    write(*,*) "Author: Zhiqi Huang"
  end subroutine coop_print_info

  subroutine coop_setup_global_cosmology_with_h(h)
    COOP_REAL h
    type(coop_arguments) args
    COOP_INT l
    call COOP_COSMO%free()
    call COOP_COSMO%init(name = "COOP_GLOBAL_COSMOLOGY",  id = 0, h = h)
    if(h.le.0.d0)return  !!return for bad h
    coop_global_baryon = coop_baryon(COOP_OMEGABH2/h**2)
    call COOP_COSMO%add_species(coop_global_baryon)
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
       coop_global_cdm = coop_cdm(COOP_OMEGACH2/h**2)
       call COOP_COSMO%add_species(coop_global_cdm)
       coop_global_de = coop_de_lambda( &
            COOP_COSMO%Omega_k()-COOP_OMEGAK &
            )
       call COOP_COSMO%add_species(coop_global_de)
       COOP_COSMO%de_genre = COOP_PERT_NONE                     
    case(COOP_DE_W0)
       coop_global_cdm = coop_cdm(COOP_OMEGACH2/h**2)
       call COOP_COSMO%add_species(coop_global_cdm)             
       coop_global_de = coop_de_w0( &
            COOP_COSMO%Omega_k()-COOP_OMEGAK, &
            COOP_COSMO_PARAMS%r(COOP_INDEX_DE) &
            )
       call COOP_COSMO%add_species(coop_global_de)
       COOP_COSMO%de_genre = COOP_PERT_NONE                     
    case(COOP_DE_W0WA)
       coop_global_cdm = coop_cdm(COOP_OMEGACH2/h**2)
       call COOP_COSMO%add_species(coop_global_cdm)                    
       coop_global_de = coop_de_w0wa( &
            COOP_COSMO%Omega_k()-COOP_OMEGAK, &
            COOP_COSMO_PARAMS%r(COOP_INDEX_DE), &
            COOP_COSMO_PARAMS%r(COOP_INDEX_DE+1) &
            )
       call COOP_COSMO%add_species(coop_global_de)
       COOP_COSMO%de_genre = COOP_PERT_NONE                     
    case(COOP_DE_QUINTESSENCE)
       coop_global_cdm = coop_cdm(COOP_OMEGACH2/h**2)
       call COOP_COSMO%add_species(coop_global_cdm)                    
       coop_global_de = coop_de_quintessence( &
            COOP_COSMO%Omega_k()-COOP_OMEGAK, &
            COOP_COSMO_PARAMS%r(COOP_INDEX_DE), &
            COOP_COSMO_PARAMS%r(COOP_INDEX_DE+1), &
            COOP_COSMO_PARAMS%r(COOP_INDEX_DE+2) &
            )
       call COOP_COSMO%add_species(coop_global_de)
       COOP_COSMO%de_genre = COOP_PERT_NONE              
    case(COOP_DE_COUPLED_QUINTESSENCE)
       call coop_background_add_coupled_DE(COOP_COSMO, Omega_c = COOP_OMEGACH2/h**2, Q =  COOP_COSMO_PARAMS%r(COOP_INDEX_DE), tracking_n =  COOP_COSMO_PARAMS%r(COOP_INDEX_DE+1), dlnQdphi =  COOP_COSMO_PARAMS%r(COOP_INDEX_DE+2), dUdphi =  COOP_COSMO_PARAMS%r(COOP_INDEX_DE+3), d2Udphi2 =  COOP_COSMO_PARAMS%r(COOP_INDEX_DE+4))
       COOP_COSMO%de_genre = COOP_PERT_SCALAR_FIELD       
       coop_global_cdm = COOP_COSMO%species(COOP_COSMO%index_of("CDM"))
       coop_global_de = COOP_COSMO%species(COOP_COSMO%index_of("Dark Energy"))       
    case default
       stop "UNKNOWN DARK ENERGY MODEL"
    end select
    if(COOP_COSMO%h().le.0.d0)return  !!rejected model
    if(coop_global_cosmology_do_firstorder .or. coop_global_cosmology_do_background)then
       call COOP_COSMO%setup_background()
    endif
    if(coop_global_cosmology_do_firstorder)then
       COOP_COSMO%optre = COOP_TAU
       call COOP_COSMO%set_xe()
       COOP_COSMO%As = 1.d-10 * exp(COOP_LN10TO10AS)
       COOP_COSMO%ns = COOP_NS
       COOP_COSMO%nrun = COOP_NRUN
       COOP_COSMO%r =  COOP_AMP_RATIO
       COOP_COSMO%nt = COOP_NT
       COOP_COSMO%inflation_consistency = COOP_INFLATION_CONSISTENCY
       call coop_setup_pp()
       call COOP_COSMO%set_power(coop_pp_get_power, args)
       call COOP_COSMO%compute_source(0)
       do l = coop_pp_lmin, coop_pp_lmax
          coop_pp_ells(l) = dble(l)
       enddo

       call COOP_COSMO%source(0)%get_all_Cls(coop_pp_lmin, coop_pp_lmax, coop_pp_scalar_Cls)       
       if(COOP_COSMO%has_tensor)then
          call COOP_COSMO%compute_source(2)
          call COOP_COSMO%source(2)%get_all_Cls(coop_pp_lmin, coop_pp_lmax, coop_pp_tensor_Cls)
       else
          coop_pp_tensor_cls = 0.d0
       endif
       call coop_get_lensing_Cls(coop_pp_lmin, coop_pp_lmax, coop_pp_scalar_Cls, coop_pp_lensed_Cls)
       coop_pp_total_cls = coop_pp_scalar_cls + coop_pp_lensed_cls + coop_pp_tensor_cls
    endif
  end subroutine coop_setup_global_cosmology_with_h


  subroutine coop_setup_global_cosmology()
    COOP_REAL,parameter::hmin = 0.5d0, hmax = 0.9d0
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


  subroutine coop_setup_pp()
    COOP_REAL,dimension(:),allocatable::lnk, lnps, lnps2
    COOP_REAL  dlnk
    COOP_INT  coop_pp_nleft, coop_pp_nright, nknots, i
    select case(COOP_PP_MODEL)
    case(COOP_PP_STANDARD)
       call coop_set_uniform(coop_pp_n, coop_pp_lnkMpc, coop_pp_lnkmin, coop_pp_lnkmax)
       coop_pp_lnps = COOP_LN10TO10AS - 10.d0*coop_ln10 + ( COOP_NS - 1.d0 ) * (coop_pp_lnkMpc - coop_pp_scalar_lnkpivot) + (COOP_NRUN/2.d0) *   (coop_pp_lnkMpc - coop_pp_scalar_lnkpivot) ** 2 + (COOP_NRUNRUN/6.d0) *  (coop_pp_lnkMpc - coop_pp_scalar_lnkpivot) ** 3
    case(COOP_PP_BUMP)
       call coop_set_uniform(coop_pp_n, coop_pp_lnkMpc, coop_pp_lnkmin, coop_pp_lnkmax)
       coop_pp_lnps = COOP_LN10TO10AS - 10.d0*coop_ln10 + ( COOP_NS - 1.d0 ) * (coop_pp_lnkMpc - coop_pp_scalar_lnkpivot) + (COOP_NRUN/2.d0) *   (coop_pp_lnkMpc - coop_pp_scalar_lnkpivot) ** 2 + (COOP_NRUNRUN/6.d0) *  (coop_pp_lnkMpc - coop_pp_scalar_lnkpivot) ** 3 + COOP_COSMO_PARAMS%r(COOP_INDEX_PP+cosmomc_pp_num_origin) * exp(-((coop_pp_lnkMpc - COOP_COSMO_PARAMS%r(COOP_INDEX_PP+cosmomc_pp_num_origin+1))/COOP_COSMO_PARAMS%r(COOP_INDEX_PP+cosmomc_pp_num_origin+2))**2/2)

    case(COOP_PP_SCAN_SPLINE)
       nknots =  COOP_NUM_PP - cosmomc_pp_num_origin + 1
       if(nknots .lt. 5) stop "You need at least 5 knots for scan_spline mode"
       coop_pp_nleft = nint((nknots - 1)* (coop_pp_scalar_lnkpivot-coop_pp_lnkmin) / (-coop_pp_lnkmin))
       coop_pp_nright = nknots - 1 - coop_pp_nleft
       dlnk = (coop_pp_scalar_lnkpivot-coop_pp_lnkmin)/coop_pp_nleft
       coop_pp_lnk_per_knot = dlnk
       coop_pp_lnkmax = coop_pp_lnkmin + (nknots-1)*dlnk
       call coop_set_uniform(coop_pp_n, coop_pp_lnkMpc, coop_pp_lnkmin, coop_pp_lnkmax)       
       allocate(lnk(nknots), lnps(nknots), lnps2(nknots))
       call coop_set_uniform(nknots, lnk, -dlnk*coop_pp_nleft, dlnk*coop_pp_nright)
       lnps(1:coop_pp_nleft) = COOP_COSMO_PARAMS%r(COOP_INDEX_PP+cosmomc_pp_num_origin:COOP_INDEX_PP+cosmomc_pp_num_origin+coop_pp_nleft -1)
       lnps(coop_pp_nleft+1) = 0.d0
       lnps(coop_pp_nleft+2:nknots) = COOP_COSMO_PARAMS%r(COOP_INDEX_PP+cosmomc_pp_num_origin+coop_pp_nleft:COOP_INDEX_PP+COOP_NUM_PP-1)
       call coop_spline(nknots, lnk, lnps, lnps2)
       do i=1, coop_pp_n
          call coop_splint(nknots, lnk, lnps, lnps2, coop_pp_lnkMpc(i) - coop_pp_scalar_lnkpivot, coop_pp_lnps(i))
       enddo
       !!modified to resolve the bump issue
       coop_pp_lnps = coop_pp_lnps + COOP_LN10TO10AS -10.d0*coop_ln10&
            + ( COOP_NS - 1.d0 ) * (coop_pp_lnkMpc - coop_pp_scalar_lnkpivot)
       deallocate(lnk, lnps,  lnps2)
    case(COOP_PP_SCAN_LINEAR)
       nknots =  COOP_NUM_PP - cosmomc_pp_num_origin + 1
       if(nknots .lt. 5) stop "You need at least 5 knots for scan_linear mode"
       coop_pp_nleft = nint((nknots - 1)* (coop_pp_scalar_lnkpivot-coop_pp_lnkmin) / (-coop_pp_lnkmin))
       coop_pp_nright = nknots - 1 - coop_pp_nleft
       dlnk = (coop_pp_scalar_lnkpivot-coop_pp_lnkmin)/coop_pp_nleft
       coop_pp_lnk_per_knot = dlnk       
       coop_pp_lnkmax = coop_pp_lnkmin + (nknots-1)*dlnk
       call coop_set_uniform(coop_pp_n, coop_pp_lnkMpc, coop_pp_lnkmin, coop_pp_lnkmax)       
       allocate(lnk(nknots), lnps(nknots), lnps2(nknots))
       call coop_set_uniform(nknots, lnk, -dlnk*coop_pp_nleft, dlnk*coop_pp_nright)
       lnps(1:coop_pp_nleft) = COOP_COSMO_PARAMS%r(COOP_INDEX_PP+cosmomc_pp_num_origin:COOP_INDEX_PP+cosmomc_pp_num_origin+coop_pp_nleft -1)
       lnps(coop_pp_nleft+1) = 0.d0
       lnps(coop_pp_nleft+2:nknots) = COOP_COSMO_PARAMS%r(COOP_INDEX_PP+cosmomc_pp_num_origin+coop_pp_nleft:COOP_INDEX_PP+COOP_NUM_PP-1)
       do i=1, coop_pp_n
          call coop_linear_interp(nknots, lnk, lnps, coop_pp_lnkMpc(i) - coop_pp_scalar_lnkpivot, coop_pp_lnps(i))
       enddo
       !!modified to resolve the bump issue
       coop_pp_lnps = coop_pp_lnps + COOP_LN10TO10AS - 10.d0*coop_ln10 &
            + ( COOP_NS - 1.d0 ) * (coop_pp_lnkMpc - coop_pp_scalar_lnkpivot)
       deallocate(lnk, lnps, lnps2)
    case(COOP_PP_GENERAL_SINGLE_FIELD)
       stop "not applied yet"
    end select
    call coop_spline(coop_pp_n, coop_pp_lnkMpc, coop_pp_lnps, coop_pp_lnps2)

    if(COOP_AMP_RATIO .gt. 0.d0)then
       select case(COOP_PP_MODEL)
       case(COOP_PP_STANDARD, COOP_PP_SCAN_SPLINE, COOP_PP_SCAN_LINEAR, COOP_PP_BUMP)
          coop_pp_lnpt = log(COOP_AMP_RATIO * coop_primordial_ps(exp(coop_pp_tensor_lnkpivot)))+(COOP_NT)*(coop_pp_lnkMpc - coop_pp_tensor_lnkpivot) + COOP_NTRUN*(coop_pp_lnkMpc - coop_pp_tensor_lnkpivot)**2
       case (COOP_PP_GENERAL_SINGLE_FIELD)
          write(*,*) "COOP_PP_GENERAL_SINGLE_FIELD not applied yet"
          stop          
       case default
          write(*,*) "COOP_PP_MODEL"//COOP_STR_OF(COOP_PP_MODEL)//"not applied yet"
          stop
       end select
       call coop_spline(coop_pp_n, coop_pp_lnkMpc, coop_pp_lnpt, coop_pp_lnpt2)
    else
       coop_pp_lnpt = -50.d0
       coop_pp_lnpt2 = 0.d0
    endif
  end subroutine coop_setup_pp


  function coop_primordial_lnps(kMpc)  result(lnps)
    COOP_REAL kMpc, lnps
    call coop_splint(coop_pp_n, coop_pp_lnkMpc, coop_pp_lnps, coop_pp_lnps2, log(kMpc), lnps)
  end function coop_primordial_lnps


  function coop_primordial_ps(kMpc)  result(ps)
    COOP_REAL kMpc, ps
    ps = exp(coop_primordial_lnps(kMpc))
  end function coop_primordial_ps

  function coop_primordial_lnpt(kMpc)  result(lnpt)
    COOP_REAL kMpc, lnpt
    if(COOP_AMP_RATIO .eq. 0.d0)then
       lnpt = -50.d0
    else
       call coop_splint(coop_pp_n, coop_pp_lnkMpc, coop_pp_lnpt, coop_pp_lnpt2, log(kMpc), lnpt)
    endif
  end function coop_primordial_lnpt

  function coop_primordial_pt(kMpc)  result(pt)
    COOP_REAL kMpc, pt
    if(COOP_AMP_RATIO .eq. 0.d0)then
       pt = 0.d0
    else
       pt = exp(coop_primordial_lnpt(kMpc))
    endif
  end function coop_primordial_pt

  subroutine coop_pp_get_potential()
    integer i
    COOP_REAL, parameter::max_delta = 0.4
    COOP_REAL, parameter::max_lneps = log(0.4)
    COOP_REAL dlnk, fourdlnk, dphiby2(coop_pp_n), eps(coop_pp_n), delta(coop_pp_n)
    coop_pp_ipivot = coop_minloc(abs(coop_pp_lnkMpc - coop_pp_scalar_lnkpivot))
    dlnk = coop_pp_lnkMpc(2)-coop_pp_lnkMpc(1)
    fourdlnk = dlnk*4.d0
    do i=2, coop_pp_n - 1
       delta(i) = (coop_pp_lnps(i-1)-coop_pp_lnps(i+1))/fourdlnk
       if(abs(delta(i)) .gt. max_delta)then
          delta(i) = sign(max_delta, delta(i))
       endif
    enddo
    delta(1) = delta(2)
    delta(coop_pp_n) = delta(coop_pp_n  - 1)
    coop_pp_lneps = min(coop_pp_lnpt - coop_pp_lnps - log(16.d0), max_lneps) 
    eps = exp(coop_pp_lneps)
    coop_pp_lneps = min(max_lneps,  coop_pp_lneps - log( (1.d0 - (2.d0*(coop_ln2+coop_EulerC-1.d0))*eps ) /(1.d0-2.d0*eps+(2.d0*(2.d0-coop_EulerC-coop_ln2))*delta))) !!slow-roll correction
    eps = exp(coop_pp_lneps)
    coop_pp_lnH = log(coop_pi2/2.d0*exp(coop_pp_lnpt)/( 1.d0 - (2.d0*(coop_ln2+coop_EulerC-1.d0))*eps ))/2.d0
    coop_pp_lnV = 2.d0*coop_pp_lnH + log(3.d0*(1.d0-eps/3.d0))
    coop_pp_phi(1) = 0.d0
    dphiby2 = sqrt(2.d0*eps)/(1.d0-eps)*dlnk/2.d0
    do i=2, coop_pp_n 
       coop_pp_phi(i) = coop_pp_phi(i-1) + (dphiby2(i-1)+dphiby2(i))
    enddo
    coop_pp_phi = coop_pp_phi - coop_pp_phi(coop_pp_ipivot)
  end subroutine coop_pp_get_potential


  subroutine coop_pp_get_power(kbykpiv, ps, pt, cosmology, args)
    COOP_REAL kbykpiv, ps, pt, kMpc
    type(coop_cosmology_firstorder)::cosmology
    type(coop_arguments)::args
    kMpc = exp(coop_pp_scalar_lnkpivot)*kbykpiv
    ps = coop_primordial_ps(kMpc)
    pt = coop_primordial_pt(kMpc)
  end subroutine coop_pp_get_power
  

end module coop_wrapper
