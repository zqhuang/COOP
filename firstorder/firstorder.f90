module coop_firstorder_mod
  use coop_wrapper_background
  use coop_pertobj_mod
  implicit none
#include "constants.h"

  private



  public::coop_cosmology_firstorder, coop_cosmology_firstorder_source,  coop_recfast_get_xe, coop_power_lnk_min, coop_power_lnk_max,  coop_k_dense_fac, coop_index_ClTT, coop_index_ClTE, coop_index_ClEE, coop_index_ClBB, coop_index_ClLenLen, coop_index_ClTLen,  coop_num_Cls, coop_Cls_lmax, coop_bbks_trans, coop_index_source_T, coop_index_source_E, coop_index_source_B, coop_index_source_Len, coop_index_source_zeta


  !!this makes the code faster and more accurate
  logical,parameter:: coop_firstorder_optimize = .true.
  COOP_INT, parameter::coop_limber_ell = 600
  logical,parameter::coop_do_limber_separately = .true.      
  
  COOP_INT::coop_Cls_lmax(0:2) = (/ 2800, 2000, 1500 /)

  COOP_REAL, parameter :: coop_power_lnk_min = log(0.1d0) 
  COOP_REAL, parameter :: coop_power_lnk_max = log(5.d3) 
  COOP_REAL, parameter :: coop_visibility_amin = 1.8d-4
  COOP_REAL, parameter :: coop_initial_condition_epsilon = 1.d-7
  COOP_REAL, parameter :: coop_cosmology_firstorder_ode_accuracy = 1.d-8
  COOP_REAL, parameter :: coop_cosmology_firstorder_tc_cutoff = 0.001d0
  COOP_REAL, parameter ::  coop_cosmology_firstorder_omega_rad_cutoff = 1.d-2


  COOP_REAL, dimension(0:2), parameter::coop_source_tau_step_factor = (/ 1.d0, 1.d0, 1.d0 /)
  COOP_REAL, dimension(0:2), parameter::coop_source_k_weight = (/ 0.15d0, 0.15d0, 0.1d0 /)
  COOP_INT, dimension(0:2), parameter::coop_source_k_n = (/ 170, 120, 100 /)
  COOP_REAL, parameter::coop_source_k_index = 0.55d0
  COOP_INT, parameter:: coop_k_dense_fac = 40

  COOP_INT, parameter::coop_index_source_T = 1
  COOP_INT, parameter::coop_index_source_E = 2
  COOP_INT, parameter::coop_index_source_B = 3
  COOP_INT, parameter::coop_index_source_Len = 3
  COOP_INT, parameter::coop_index_source_zeta = 4  

  COOP_INT, parameter::coop_index_ClTT = 1
  COOP_INT, parameter::coop_index_ClEE = 2
  COOP_INT, parameter::coop_index_ClBB = 3
  COOP_INT, parameter::coop_index_ClTE = 4
  COOP_INT, parameter::coop_index_ClLenLen = 5
  COOP_INT, parameter::coop_index_ClTLen = 6


  !!how many source terms you want to extract & save

#if DO_ZETA_TRANS
  COOP_INT, parameter::coop_index_ClTzeta = 7
  COOP_INT, parameter::coop_index_ClEzeta = 8
  COOP_INT, parameter::coop_index_Clzetazeta = 9
  COOP_INT, parameter::coop_num_Cls =  coop_index_Clzetazeta
  public::  coop_index_ClTzeta,  coop_index_ClEzeta, coop_index_clzetazeta
  COOP_INT, dimension(0:2), parameter::coop_num_sources = (/ 4,  3,  3 /)
#else
  COOP_INT, parameter::coop_num_Cls =  coop_index_ClTLen
  COOP_INT, dimension(0:2), parameter::coop_num_sources = (/ 3,  3,  3 /)
#endif
  COOP_INT, dimension(0:2), parameter::coop_num_saux = (/ 7,  1,  1 /)

!!recfast head file
#include "recfast_head.h"


  type coop_cosmology_firstorder_source
     COOP_INT::m = 0
     COOP_INT::ntau = 0
     COOP_INT::nk = 0
     COOP_INT::nsrc = 0
     COOP_INT::nsaux = 0
     COOP_INT, dimension(:),allocatable::index_tc_off, index_rad_off
     COOP_INT, dimension(coop_pert_default_nq)::index_massivenu_on
     COOP_INT::index_massivenu_cold, index_rad_small, index_vis_end, index_vis_max, index_vis_start
     COOP_REAL::dkop, kopmin, kopmax, kmin, kmax, kweight, tauweight, bbks_keq, bbks_trans_kmax, distlss, dkop_dense
     COOP_REAL,dimension(coop_k_dense_fac)::a_dense, b_dense, a2_dense, b2_dense
     COOP_REAL, dimension(:),allocatable::k, kop, dk !!tau is conformal time, chi is comoving distance; in flat case, chi + tau = tau_0
     COOP_REAL, dimension(:),allocatable::tau, a, tauc, lna, dtau, chi, omega_rad, vis
     COOP_REAL, dimension(:,:),allocatable::k_dense, ws_dense, wt_dense, ps_dense
     COOP_REAL, dimension(:,:,:),allocatable::s, s2, saux
   contains
     procedure::free => coop_cosmology_firstorder_source_free
     procedure::get_transfer => coop_cosmology_firstorder_source_get_transfer
     procedure::get_Cls => coop_cosmology_firstorder_source_get_Cls
     procedure::get_Cls_limber => coop_cosmology_firstorder_source_get_Cls_limber     
     procedure::get_All_Cls => coop_cosmology_firstorder_source_get_All_Cls
     procedure::kop2k => coop_cosmology_firstorder_source_kop2k
     procedure::k2kop => coop_cosmology_firstorder_source_k2kop
     procedure::interpolate => coop_cosmology_firstorder_source_interpolate
     procedure::intbypart => coop_cosmology_firstorder_source_intbypart
     procedure::interpolate_one => coop_cosmology_firstorder_source_interpolate_one
     procedure::get_Psi_trans => coop_cosmology_firstorder_source_get_psi_trans
     procedure::get_PhiPlusPsi_trans => coop_cosmology_firstorder_source_get_PhiPlusPsi_trans
  end type coop_cosmology_firstorder_source
  
  type, extends(coop_cosmology_background) :: coop_cosmology_firstorder
     logical::do_reionization = .true.
     COOP_REAL::zrecomb, distlss, tau0, zrecomb_start, maxvis, taurecomb, arecomb, zrecomb_end, arecomb_start
     COOP_REAL::optre = 0.07d0
     COOP_REAL::zre = 8.d0
     COOP_REAL::deltaz = 1.5d0
     COOP_REAL::kMpc_pivot = 0.05d0
     COOP_INT ::de_genre = COOP_PERT_NONE
     COOP_REAL::k_pivot 
     COOP_REAL::dkappadtau_coef, ReionFrac, Omega_b, Omega_c, Omega_nu, Omega_g, tau_eq, mnu_by_Tnu, As, ns, nrun, r, nt, Omega_massivenu, bbks_keq, sigma_8
     !!these omega's are defined at a=1 (today)
     COOP_REAL::ombh2, omch2 !!these two parameters are defined in the a->0 limit (only matters when baryon or CDM is coupled to DE in some modified gravity models)
     logical::inflation_consistency
     type(coop_function)::Ps, Pt, Xe, ekappa, vis, Tb
     type(coop_cosmology_firstorder_source),dimension(0:2)::source
     COOP_INT::index_baryon, index_cdm, index_radiation, index_nu, index_massiveNu, index_de
     COOP_REAL, dimension(0:coop_pert_default_lmax, 0:coop_pert_default_mmax, 0:coop_pert_default_smax)::klms, klms_by_2lm1, klms_by_2lp1
     COOP_REAL,dimension(coop_pert_default_lmax)::fourbyllp1
     logical::klms_done = .false.
     logical::has_tensor = .false.
   contains
     procedure:: set_standard_cosmology =>  coop_cosmology_firstorder_set_standard_cosmology
     procedure:: set_standard_power => coop_cosmology_firstorder_set_standard_power
     procedure:: set_Planck_bestfit =>coop_cosmology_firstorder_set_Planck_bestfit
     procedure:: set_Planck_bestfit_with_r =>coop_cosmology_firstorder_set_Planck_bestfit_with_r
     procedure:: set_klms => coop_cosmology_firstorder_set_klms
     procedure:: set_source_tau => coop_cosmology_firstorder_set_source_tau
     procedure:: set_source_k => coop_cosmology_firstorder_set_source_k
     procedure:: set_power => coop_cosmology_firstorder_set_power
     procedure:: set_xe => coop_cosmology_firstorder_set_xe
     procedure:: set_zre_from_optre => coop_cosmology_firstorder_set_zre_from_optre
     procedure:: set_optre_from_zre => coop_cosmology_firstorder_set_optre_from_zre
     procedure:: set_initial_conditions => coop_cosmology_firstorder_set_initial_conditions
     procedure:: cosmomc_theta => coop_cosmology_firstorder_cosmomc_theta
     procedure:: xeofa => coop_cosmology_firstorder_xeofa
     procedure:: cs2bofa => coop_cosmology_firstorder_cs2bofa
     procedure:: Tbofa => coop_cosmology_firstorder_Tbofa
     procedure:: dxeda => coop_cosmology_firstorder_dxeda
     procedure:: dlnxedlna => coop_cosmology_firstorder_dlnxedlna
     procedure:: dkappada => coop_cosmology_firstorder_dkappada
     procedure:: dkappadtau => coop_cosmology_firstorder_dkappadtau
     procedure:: taucofa => coop_cosmology_firstorder_taucofa
     procedure:: dot_tauc => coop_cosmology_firstorder_dot_tauc
     procedure:: visofa => coop_cosmology_firstorder_visofa
     procedure:: ekappaofa => coop_cosmology_firstorder_ekappaofa
     procedure:: psofk => coop_cosmology_firstorder_psofk
     procedure:: ptofk => coop_cosmology_firstorder_ptofk
     procedure:: Clzetazeta => coop_cosmology_firstorder_Clzetazeta
     procedure:: Clzetazeta_at_r => coop_cosmology_firstorder_Clzetazeta_at_r
     procedure:: free => coop_cosmology_firstorder_free
     procedure:: pert2source => coop_cosmology_firstorder_pert2source
     procedure:: init_source => coop_cosmology_firstorder_init_source
     procedure:: compute_source =>  coop_cosmology_firstorder_compute_source
     procedure:: compute_source_k =>  coop_cosmology_firstorder_compute_source_k
     procedure:: get_matter_power => coop_cosmology_firstorder_get_matter_power
     procedure:: sigma_Tophat_R => coop_cosmology_firstorder_sigma_Tophat_R
     procedure:: sigma_Gaussian_R => coop_cosmology_firstorder_sigma_Gaussian_R
     procedure:: sigma_Gaussian_R_quick => coop_cosmology_firstorder_sigma_Gaussian_R_quick
     procedure::sigma_Gaussian_R_with_dervs => coop_cosmology_firstorder_sigma_Gaussian_R_with_dervs
     !!
  end type coop_cosmology_firstorder


contains


!!recfast code
#include "recfast_source.h"

!!this head file contains the evolution equations of the firstorder ODE system
#include "firstorder_equations.h"

!!this head file set the outputs (saved in this%source%s)
#include "firstorder_source.h"

!!this head file sets the initial conditions
#include "firstorder_ic.h"

  subroutine coop_cosmology_firstorder_set_Planck_Bestfit(this, Omega_nu)
    class(coop_cosmology_firstorder)::this
    COOP_REAL, optional::Omega_nu
    if(present(Omega_nu))then
       call this%set_standard_cosmology(Omega_b=0.04904d0, Omega_c=0.2642d0, h = 0.6731d0, tau_re = 0.078d0, As = 2.195d-9, ns = 0.9655d0,  Omega_nu = Omega_nu)
    else
       call this%set_standard_cosmology(Omega_b=0.04904d0, Omega_c=0.2642d0, h = 0.6731d0, tau_re = 0.078d0, As = 2.195d-9, ns = 0.9655d0)       
     
    endif
  end subroutine coop_cosmology_firstorder_set_Planck_Bestfit


  subroutine coop_cosmology_firstorder_set_Planck_Bestfit_with_r(this, r, Omega_nu)
    class(coop_cosmology_firstorder)::this
    COOP_REAL r
    COOP_REAL, optional::Omega_nu
    if(present(Omega_nu))then
       call this%set_standard_cosmology(Omega_b=0.04904d0, Omega_c=0.2642d0, h = 0.6731d0, tau_re = 0.078d0, As = 2.195d-9, ns = 0.9655d0,  Omega_nu = Omega_nu, r = r)
    else
       call this%set_standard_cosmology(Omega_b=0.04904d0, Omega_c=0.2642d0, h = 0.6731d0, tau_re = 0.078d0, As = 2.195d-9, ns = 0.9655d0, r = r)       
     
    endif
  end subroutine coop_cosmology_firstorder_set_Planck_Bestfit_with_r



  subroutine coop_cosmology_firstorder_set_standard_cosmology(this, h, omega_b, omega_c, tau_re, nu_mass_eV, Omega_nu, As, ns, nrun, r, nt, inflation_consistency, Nnu, YHe, de_Q, de_tracking_n, de_dlnQdphi, de_dUdphi, de_d2Udphi2)
    class(coop_cosmology_firstorder)::this
    COOP_REAL:: h, tau_re, Omega_b, Omega_c
    COOP_REAL, optional::nu_mass_eV, Omega_nu, As, ns, nrun, r, nt, Nnu, YHe, de_Q, de_tracking_n, de_dlnQdphi, de_dUdphi, de_d2Udphi2
    COOP_REAL::Q0, Q1, U1, U2, tracking_n
    logical::scalar_de
    logical,optional::inflation_consistency
    if(present(YHe))then
       if(present(Nnu))then
          call this%init(h=h, YHe = YHe, Nnu = NNu)
       else
          call this%init(h=h, YHe = YHe)
       endif
    else
       if(present(Nnu))then
          call this%init(h=h, Nnu = NNu)
       else
          call this%init(h=h)
       endif
    endif
    call this%add_species(coop_baryon(COOP_REAL_OF(Omega_b)))
    call this%add_species(coop_radiation(this%Omega_radiation()))
    if(present(nu_mass_eV))then
       if(present(Omega_nu))then
          call coop_return_error("set_standard_cosmology", "you cannot have both Omega_nu and nu_mass_eV in the arguments", "stop")
       endif
       if(nu_mass_eV .gt. 1.d-3)then
          call this%add_species(coop_neutrinos_massless(this%Omega_massless_neutrinos_per_species()*(this%Nnu()-1)))
          call this%add_species(coop_neutrinos_massive(this%Omega_nu_per_species_from_mnu_eV(nu_mass_eV), this%Omega_massless_neutrinos_per_species()))
       else
          call this%add_species(coop_neutrinos_massless(this%Omega_massless_neutrinos_per_species()*(this%Nnu())))
       endif
    elseif(present(Omega_nu))then
       if(Omega_nu .gt. this%Omega_massless_neutrinos_per_species()*1.01d0)then
          call this%add_species(coop_neutrinos_massless(this%Omega_massless_neutrinos_per_species()*(this%Nnu()-1)))
          call this%add_species(coop_neutrinos_massive(Omega_nu, this%Omega_massless_neutrinos_per_species()))
       else
          call this%add_species(coop_neutrinos_massless(this%Omega_massless_neutrinos_per_species()*(this%Nnu())))
       endif
    else
       call this%add_species(coop_neutrinos_massless(this%Omega_massless_neutrinos_per_species()*(this%Nnu())))
    endif
    scalar_de = .false.
    if(present(de_dUdphi))then
       U1 = de_dUdphi
       scalar_de = .true.
    else
       U1 = 0.d0
    endif
    if(present(de_d2Udphi2))then
       U2 = de_d2Udphi2
       scalar_de = .true.       
    else
       U2 = 0.d0
    endif
    if(present(de_Q))then
       Q0 = de_Q
       scalar_de = .true.       
       if(present(de_dlnQdphi))then
          Q1 = de_dlnQdphi
       else
          Q1 = 0.d0
       endif
    else
       Q0 = 0.d0
       Q1 = 0.d0
    endif
    if(present(de_tracking_n))then
       scalar_de = .true.       
       tracking_n = de_tracking_n
    else
       tracking_n = 0.d0
    endif
    
    if(scalar_de)then
       call coop_background_add_coupled_DE(this, Omega_c = Omega_c, Q = Q0, tracking_n = tracking_n, dlnQdphi = Q1, dUdphi = U1, d2Udphi2 = U2)
       this%de_genre = COOP_PERT_SCALAR_FIELD
    else
       call this%add_species(coop_cdm(omega_c))
       call this%add_species(coop_de_lambda(this%Omega_k()))
       this%de_genre = COOP_PERT_NONE       
    endif
    call this%setup_background()
    this%optre = tau_re
    call this%set_xe()
    if(present(As).or.present(ns) .or. present(nrun) .or. present(r) .or. present(nt) .or. present(inflation_consistency))then
       if(present(As))then
          this%As = As
       else
          this%As = 1.d0
       endif
       if(present(ns))then
          this%ns = ns
       else
          this%ns = 1.d0
       endif
       if(present(nrun))then
          this%nrun = nrun
       else
          this%nrun = 0.d0
       endif
       if(present(r))then
          this%r = r
       else
          this%r = 0.d0
       endif
       if(present(nt))then
          this%nt = nt
       else
          this%nt = 0.d0
       endif
       if(present(inflation_consistency))then
          this%inflation_consistency = inflation_consistency 
       else
          this%inflation_consistency = .true.
       endif
       call this%set_standard_power(this%As, this%ns, this%nrun, this%r, this%nt, this%inflation_consistency)
    endif
  end subroutine coop_cosmology_firstorder_set_standard_cosmology


  

  subroutine coop_cosmology_firstorder_source_get_transfer(source, l, trans)
    class(coop_cosmology_firstorder_source)::source
    COOP_INT::l
    COOP_REAL,dimension(:,:,:)::trans
    COOP_INT::n, ik, idense, itau, ii, nchis, ikcut, kchicut, itau_cut
    COOP_REAL::jl, xmin, xmax, dchi, chis
    COOP_REAL,dimension(:),allocatable::kmin, kmax, kfine
    if(size(trans,2).ne. coop_k_dense_fac .or. size(trans, 3).ne. source%nk .or. size(trans, 1) .ne. source%nsrc) call coop_return_error("get_transfer", "wrong size", "stop")
    if(l.ge.coop_limber_ell)then
       itau_cut = source%index_vis_end
    else
       itau_cut = source%ntau
    endif
    trans = 0.d0
    call coop_jl_startpoint(l, xmin)
    call coop_jl_check_init(l)
    xmax = (coop_jl_zero(l, 8)+coop_jl_zero(l, 9))/2.d0
    allocate(kmin(source%ntau), kmax(source%ntau), kfine(source%ntau))
    do itau = 1, itau_cut
       kmin(itau) = xmin/source%chi(itau)
       kmax(itau) = xmax/source%chi(itau)
       kfine(itau) = 1.8d0/source%dtau(itau) !!k > kfine use fine grid
    enddo
    ikcut = source%nk
    kchicut = min(max(l*2.8d0, 4000.d0), l*100.d0)
    do while(source%k(ikcut) * source%chi(1) .gt. kchicut .and. ikcut .gt. 2)
       ikcut = ikcut - 1
    enddo
    !$omp parallel do private(ik, itau, idense, jl, ii, dchi, nchis)
    do ik=2, ikcut
       do itau = source%index_vis_start, itau_cut
          if(source%k(ik)*source%chi(itau) .lt. xmin) exit
          if(source%k(ik) .lt. kfine(itau))then
             do idense = 1, coop_k_dense_fac
                jl = coop_jl(l, source%k_dense(idense, ik)*source%chi(itau))
                trans(:, idense, ik) = trans(:, idense, ik) + jl* source%dtau(itau)* COOP_INTERP_SOURCE(source, :, idense, ik, itau) 
             enddo
          elseif(l.lt.coop_limber_ell)then
             if(source%k(ik) .gt. kmin(itau) .and. source%k_dense(1,ik) .lt. kmax(itau))then            
                do idense = 1, coop_k_dense_fac
                   nchis = 1+ceiling(source%k_dense(idense, ik)*source%dtau(itau))
                   dchi = source%dtau(itau)/(nchis*2+1)
                   jl = 0.d0
                   do ii = -nchis, nchis
                      jl = jl + coop_jl(l, source%k_dense(idense, ik)*(source%chi(itau)+dchi*ii))
                   enddo
                   jl = jl/(2*nchis+1.d0)
                   trans(:,idense, ik) = trans(:,idense, ik) + jl* source%dtau(itau)*COOP_INTERP_SOURCE(source, :, idense, ik, itau)
                enddo
             endif
          endif
       enddo
    enddo
    !$omp end parallel do
    deallocate(kmin, kmax, kfine)

  end subroutine coop_cosmology_firstorder_source_get_transfer

  subroutine coop_cosmology_firstorder_source_get_Cls_limber(source, l, Cls)
    class(coop_cosmology_firstorder_source)::source
    COOP_INT l, ik, idense, itau
    COOP_REAL,dimension(coop_num_cls),intent(OUT)::Cls
    COOP_REAL::src(source%nsrc), kop, rk
    Cls = 0.d0
    do itau = source%index_vis_end + 1, source%ntau
       call source%k2kop(l/source%chi(itau), kop)
       rk = (kop-source%kopmin)/source%dkop+1.d0
       ik = ceiling(rk)
       if(ik.lt.2 .or. ik .gt. source%nk)cycle
       idense = coop_k_dense_fac - nint((ik - rk)*coop_k_dense_fac)
       if(idense .eq. 0)then
          ik = ik - 1
          if(ik.eq.1)cycle
          idense = coop_k_dense_fac
       endif
       src = COOP_INTERP_SOURCE(source, :, idense, ik, itau)
       rk = source%ps_dense(idense,ik)*source%chi(itau)*source%dtau(itau)
       Cls(coop_index_ClTT) = Cls(coop_index_ClTT) + rk*src(coop_index_source_T)**2
       Cls(coop_index_ClEE) = Cls(coop_index_ClEE) + rk*src(coop_index_source_E)**2
       Cls(coop_index_ClTE) = Cls(coop_index_ClTE) + rk*src(coop_index_source_T)*src(coop_index_source_E)
       if(source%m .eq. 0)then
          if(source%nsrc .ge. 3)then
             Cls(coop_index_ClLenLen) = Cls(coop_index_ClLenLen) + rk*src(coop_index_source_Len)**2
             Cls(coop_index_ClTLen) = Cls(coop_index_ClTLen) + rk*src(coop_index_source_T)*src(coop_index_source_Len)
#if DO_ZETA_TRANS             
             if(source%nsrc .ge. 4)then
                Cls(coop_index_ClTzeta) = Cls(coop_index_ClTzeta) + rk*src(coop_index_source_zeta)*src(coop_index_source_T)
                Cls(coop_index_ClEzeta) = Cls(coop_index_ClTzeta) + rk*src(coop_index_source_zeta)*src(coop_index_source_E)
                Cls(coop_index_Clzetazeta) = Cls(coop_index_ClTzeta) + rk*src(coop_index_source_zeta)**2
             endif
#endif             
          endif
       else
          Cls(coop_index_ClBB) = Cls(coop_index_ClBB) + rk*src(coop_index_source_B)**2
       endif
    enddo
    if(source%m .eq. 0 )then
       do itau = 1, source%index_vis_end
          call source%k2kop(l/source%chi(itau), kop)
          rk = (kop-source%kopmin)/source%dkop+1.d0
          ik = ceiling(rk)
          if(ik.lt.2 .or. ik .gt. source%nk)cycle
          idense = coop_k_dense_fac - nint((ik - rk)*coop_k_dense_fac)
          if(idense .eq. 0)then
             ik = ik - 1
             if(ik.eq.1)cycle
             idense = coop_k_dense_fac
          endif
          src = COOP_INTERP_SOURCE(source, :, idense, ik, itau)
          rk = source%ps_dense(idense,ik)*source%chi(itau)*source%dtau(itau)

          if(source%nsrc .ge. 3)then
             Cls(coop_index_ClLenLen) = Cls(coop_index_ClLenLen) + rk*src(coop_index_source_Len)**2
             Cls(coop_index_ClTLen) = Cls(coop_index_ClTLen) + rk*src(coop_index_source_T)*src(coop_index_source_Len)
#if DO_ZETA_TRANS             
             if(source%nsrc .ge. 4)then
                Cls(coop_index_ClTzeta) = Cls(coop_index_ClTzeta) + rk*src(coop_index_source_zeta)*src(coop_index_source_T)
                Cls(coop_index_ClEzeta) = Cls(coop_index_ClTzeta) + rk*src(coop_index_source_zeta)*src(coop_index_source_E)
                Cls(coop_index_Clzetazeta) = Cls(coop_index_ClTzeta) + rk*src(coop_index_source_zeta)**2
             endif
#endif             
          endif
       enddo
    endif
    Cls = Cls*((2.d0*coop_pi**2)/dble(l)**3)
  end subroutine coop_cosmology_firstorder_source_get_Cls_limber



  subroutine coop_cosmology_firstorder_source_get_Cls(source, l, Cls)
    class(coop_cosmology_firstorder_source)::source
    COOP_INT l
    COOP_REAL,dimension(coop_num_cls),intent(OUT)::Cls
    COOP_REAL,dimension(coop_num_cls)::Cls_limber
    COOP_REAL,dimension(:,:,:),allocatable::trans
    COOP_INT:: i, ik, idense
    COOP_REAL::tmp
    Cls = 0.d0
    allocate(trans(source%nsrc, coop_k_dense_fac, source%nk))
    call source%get_transfer(l, trans)
    select case(source%m)
    case(0)
       Cls(coop_index_ClTT) = sum(source%ws_dense * trans(coop_index_source_T, :, :)**2)*coop_4pi
       if(source%nsrc .ge. 2)then
          Cls(coop_index_ClTE) = sqrt((l+2.d0)*(l+1.d0)*l*(l-1.d0))*sum(source%ws_dense * trans(coop_index_source_T, :, :)*trans(coop_index_source_E, :, :))*coop_4pi
          Cls(coop_index_ClEE) = (l+2.d0)*(l+1.d0)*l*(l-1.d0)*sum(source%ws_dense * trans(coop_index_source_E,:,:)**2)*coop_4pi
       endif
       if(source%nsrc.ge.3 .and. ( l .lt. coop_limber_ell) )then
          Cls(coop_index_ClLenLen) =  sum(source%ws_dense * trans(coop_index_source_Len, :, :)**2)*coop_4pi
          Cls(coop_index_ClTLen) =  sum(source%ws_dense * trans(coop_index_source_T, :, :) * trans(coop_index_source_Len,:,:))*coop_4pi
       endif
#if DO_ZETA_TRANS        
       if(source%nsrc.ge.4 .and. l .lt. coop_limber_ell)then
          Cls(coop_index_ClTzeta) = sum(source%ws_dense * trans(coop_index_source_T, :, :)*trans(coop_index_source_zeta,:,:))*coop_4pi
          Cls(coop_index_ClEzeta) = sqrt((l+2.d0)*(l+1.d0)*l*(l-1.d0))*sum(source%ws_dense * trans(coop_index_source_E, :, :) * trans(coop_index_source_zeta,:,:))*coop_4pi
          Cls(coop_index_Clzetazeta) = sum(source%ws_dense * trans(coop_index_source_zeta, :, :)**2)*coop_4pi
       endif
#endif       
    case(1)
       call coop_tbw("get_Cls: vector")
    case(2)
       Cls(coop_index_ClTT) =(l+2.d0)*(l+1.d0)*l*(l-1.d0)*sum(source%wt_dense * trans(coop_index_source_T, :, :)**2)*coop_4pi
       if(source%nsrc .ge. 2)then
          Cls(coop_index_ClTE) = sqrt((l+2.d0)*(l+1.d0)*l*(l-1.d0))*sum(source%wt_dense * trans(coop_index_source_T, :, :)*trans(coop_index_source_E,:,:))*coop_4pi
          Cls(coop_index_ClEE) = sum(source%wt_dense * trans(coop_index_source_E,:,:)**2)*coop_4pi
          if(source%nsrc .ge. 3)then
             Cls(coop_index_ClBB) = sum(source%wt_dense * trans(coop_index_source_B, :, :)**2)*coop_4pi
          end if
       endif       
    case default
       call coop_return_error("get_Cls", "unknown m = "//trim(coop_num2str(source%m)), "stop")
    end select
    deallocate(trans)
    if(l .ge. coop_limber_ell)then
       call source%get_cls_limber(l, cls_limber)
       Cls = Cls + Cls_limber
    endif
    Cls(coop_index_ClLenLen) =  Cls(coop_index_ClLenLen)*l*(l+1.d0)
    Cls(coop_index_ClTLen) =  Cls(coop_index_ClLenLen)*l
  end subroutine coop_cosmology_firstorder_source_get_Cls


  subroutine coop_cosmology_firstorder_source_get_All_Cls(source, lmin, lmax, Cls)
    class(coop_cosmology_firstorder_source)::source
    COOP_INT lmin, lmax, lmax_compute, l, nc, i
    COOP_REAL,dimension(coop_num_cls, lmin:lmax),intent(OUT)::Cls
    COOP_REAL, dimension(:),allocatable::ls_computed
    COOP_REAL, dimension(:,:),allocatable::Cls_computed, Cls2_computed
    COOP_REAL::Cls_tmp(coop_num_cls)
    COOP_REAL, parameter::norm = 1.d10
    lmax_compute = min(coop_Cls_lmax(source%m), lmax)
    l = lmin 
    nc = 1
    do
       call next_l()
       nc = nc + 1
       if(l.ge.lmax_compute)then
          l = lmax
          exit
       endif
    enddo
    allocate(ls_computed(nc), Cls_computed(nc, coop_num_Cls), Cls2_computed(nc, coop_num_Cls))
    i  = 1
    ls_computed(i) = lmin
    l = lmin
    do
       call next_l()
       i = i + 1
       if(l.ge.lmax_compute)then
          ls_computed(i) = dble(lmax_compute)
          exit
       else
          ls_computed(i) = dble(l)
       endif
    enddo

    do i=1, nc
       call source%Get_Cls(nint(ls_computed(i)), Cls_tmp)
       Cls_Computed(i, :) = Cls_tmp*(ls_computed(i)*(ls_computed(i) + 1)*norm)
    enddo
    !$omp parallel do private(i, l)
    do i=1, coop_num_Cls
       if(all(Cls_Computed(:,i).eq.0.d0))then
          Cls(i, lmin:lmax) = 0.d0
       else
          call coop_spline(nc, ls_computed, Cls_Computed(:, i), Cls2_computed(:, i))
          do l = lmin, lmax_compute
             call coop_splint(nc, ls_computed, Cls_Computed(:, i), Cls2_computed(:, i), dble(l), Cls(i, l))
          enddo
          do l = lmax_compute+1, lmax
             Cls(i, l) = Cls(i, lmax_compute)*exp(-(dble(l)**2-dble(lmax_compute)**2)/1.2d6)
          enddo
       endif
    enddo
    !$omp end parallel do
    !$omp parallel do
    do l= lmin, lmax
       Cls(:, l) = Cls(:, l)/(l*(l+1.d0)*norm)
       Cls(coop_index_ClLenLen, l) =  Cls(coop_index_ClLenLen, l)/(l*(l+1.d0))
       Cls(coop_index_ClTLen, l) =  Cls(coop_index_ClLenLen, l)/l
    enddo
    !$omp end parallel do
    deallocate(ls_computed, Cls_computed,Cls2_computed)

  contains

    subroutine next_l()
      l = l + min(38, 12 + l/70, max(1, l/4))
    end subroutine next_l

  end subroutine coop_cosmology_firstorder_source_get_All_Cls


  subroutine coop_cosmology_firstorder_compute_source(this, m)
    class(coop_cosmology_firstorder)::this
    COOP_INT m, ik, itau, is
    call this%init_source(m)
    this%source(m)%bbks_keq = this%bbks_keq
    this%source(m)%bbks_trans_kmax = coop_bbks_trans(this%source(m)%kmax/this%bbks_keq)
    !$omp parallel do
    do ik = 1, this%source(m)%nk
       call this%compute_source_k(this%source(m), ik)
    enddo
    !$omp end parallel do


    do itau = 1, this%source(m)%ntau
       do is = 1, this%source(m)%nsrc
          call coop_naturalspline_uniform(this%source(m)%nk, this%source(m)%s(is, :, itau), this%source(m)%s2(is, :, itau))
       enddo
    enddo

    if(m.eq.0)then !!interpolate Phi, Psi
       do itau = 1, this%source(m)%ntau
          call coop_naturalspline_uniform(this%source(m)%nk, this%source(m)%saux(2, :, itau), this%source(m)%saux(4, :, itau))
          call coop_naturalspline_uniform(this%source(m)%nk, this%source(m)%saux(3, :, itau), this%source(m)%saux(5, :, itau))
       enddo
       this%sigma_8 = this%sigma_Tophat_R(z = 0.d0, r = 8.d0/this%h()*this%H0Mpc())
    endif
  end subroutine coop_cosmology_firstorder_compute_source


  subroutine coop_cosmology_firstorder_compute_source_k(this, source, ik, do_test_energy_conservation, transfer_only)
    class(coop_cosmology_firstorder)::this
    type(coop_cosmology_firstorder_source)::source
    COOP_INT ik, nvars, itau, iq
    type(coop_pert_object) pert
    COOP_REAL, dimension(:,:),allocatable::w
    logical,optional::do_test_energy_conservation, transfer_only
    COOP_REAL c(24)
    COOP_INT ind, i
    COOP_REAL tau_ini, lna, mnu_deltarho, mnu_deltav, T0i, T00
    tau_ini = min(coop_initial_condition_epsilon/source%k(ik), this%conformal_time(this%a_eq*coop_initial_condition_epsilon), source%tau(1)*0.999d0)
    call this%set_initial_conditions(pert, m = source%m, k = source%k(ik), tau = tau_ini)
    lna = log(this%aoftau(tau_ini))
    call coop_cosmology_firstorder_equations(pert%ny+1, lna, pert%y, pert%yp, this, pert)
    ind = 1
    c = 0.d0
    nvars = pert%ny + 1
    allocate(w(nvars, 9))
    w = 0.d0
    iq = 1
    do itau = 1, source%ntau
       call coop_dverk_firstorder(nvars, coop_cosmology_firstorder_equations, this, pert, lna,   pert%y, source%lna(itau),  coop_cosmology_firstorder_ode_accuracy, ind, c, nvars, w)
       pert%want_source = .true.
       call coop_cosmology_firstorder_equations(pert%ny+1, lna, pert%y, pert%yp, this, pert)
       pert%want_source  = .false.              
       !!for energy conservation test:
       if(present(do_test_energy_conservation))then
          if(do_test_energy_conservation)print*, log(pert%a), pert%delta_T00a2(), pert%delta_G00a2()       
       endif
       !!------------------------------------------------------------
       !!forcing v_b - v_g to tight coupling approximations
       !!you can remove this, no big impact
       if(source%m .eq. 0)then
          if(pert%tight_coupling)then
             pert%O1_V_B = (pert%O1_V_B * pert%rhoa2_b + (pert%O1_T(1)/4.d0+pert%slip)* pert%rhoa2_g)/(pert%rhoa2_b+pert%rhoa2_g) 
             pert%O1_T(1) = (pert%O1_V_B - pert%slip)*4.d0
          endif
       endif
       !!------------------------------------------------------------
       call this%pert2source(pert, source, itau, ik)

       if(itau .eq. source%index_tc_off(ik))then
          call pert%save_ode()
          call coop_cosmology_firstorder_equations(pert%ny+1, lna, pert%y, pert%yp, this, pert)   !!set tight coupling apprixmations
          pert%tight_coupling = .false.
          call pert%init(m = source%m, nu_mass = this%mnu_by_Tnu, de_genre = this%de_genre, a = source%a(itau))
          call pert%restore_ode()
          ind = 1
          c = 0.d0
          deallocate(w)
          nvars = pert%ny + 1
          allocate(w(nvars, 9))
       endif
       if(itau .eq. source%index_rad_off(ik))then
          call pert%save_ode()
          call coop_cosmology_firstorder_equations(pert%ny+1, lna, pert%y, pert%yp, this, pert)   !!set tight coupling apprixmations
          pert%has_rad_pert  = .false.
          call pert%init(m = source%m, nu_mass = this%mnu_by_Tnu, de_genre = this%de_genre, a = source%a(itau))
          call pert%restore_ode()
          ind = 1
          c = 0.d0
          deallocate(w)
          nvars = pert%ny + 1
          allocate(w(nvars, 9))
       endif
       if(iq.le.coop_pert_default_nq)then
          do while(itau .eq. source%index_massivenu_on(iq))
             call pert%save_ode()
             call pert%init(m = source%m, nu_mass = this%mnu_by_Tnu, de_genre = this%de_genre, a = source%a(itau))             
             call pert%restore_ode()
             ind = 1
             c = 0.d0
             deallocate(w)
             nvars = pert%ny + 1
             allocate(w(nvars, 9))
             iq = iq + 1
             if(iq .gt.coop_pert_default_nq) exit
          enddo
       endif
       if(itau .eq. source%index_massivenu_cold)then
          call pert%save_ode()
          if(source%m .le. 0)then
             mnu_deltarho = pert%rhoa2_nu * (pert%deltatr_mnu + pert%deltap_mnu)/pert%rhoa2_mnu
          endif
          if(source%m .le. 1)then
             mnu_deltav = 0.d0
             do iq = 1, coop_pert_default_nq
                mnu_deltav = mnu_deltav + pert%O1_MASSIVENU(1, iq)*coop_pert_default_q_kernel(iq)
             enddo
             mnu_deltav =(pert%rhoa2_nu + pert%pa2_nu) * pert%num_mnu_ratio/4.d0 *  mnu_deltav/(pert%rhoa2_mnu + pert%pa2_mnu)
          endif
          call pert%init(m = source%m, nu_mass  = this%mnu_by_Tnu, de_genre = this%de_genre, a = source%a(itau))
          if(source%m .le. 0) pert%massivenu(1)%F(0) = mnu_deltarho
          if(source%m .le. 1) pert%massivenu(1)%F(1) = mnu_deltav
          call pert%restore_ode()
          ind = 1
          c = 0.d0
          deallocate(w)
          nvars = pert%ny + 1
          allocate(w(nvars, 9))          
       endif
    enddo
    deallocate(w)
    call pert%free()
    if(present(transfer_only))then
       if(transfer_only)return
    endif
    call source%intbypart(ik)
  end subroutine coop_cosmology_firstorder_compute_source_k


  function coop_cosmology_firstorder_xeofa(this, a) result(xe)
    COOP_REAL a, xe
    class(coop_cosmology_firstorder)::this
    xe = this%xe%eval(max(a, 1.d-99))
  end function coop_cosmology_firstorder_xeofa

  function coop_cosmology_firstorder_dxeda(this, a) result(dxeda)
    COOP_REAL a, dxeda
    class(coop_cosmology_firstorder)::this
    if(a .lt. coop_recfast_a_initial)then
       dxeda = 0.d0 
    else
       dxeda = this%xe%derivative(a)
    endif
  end function coop_cosmology_firstorder_dxeda

  function coop_cosmology_firstorder_dlnxedlna(this, a) result(dlnxedlna)
    COOP_REAL a, dlnxedlna
    class(coop_cosmology_firstorder)::this
    if(a .lt. coop_recfast_a_initial)then
       dlnxedlna = 0
    else
       dlnxedlna = this%xe%derivative_bare(log(a))
    endif
  end function coop_cosmology_firstorder_dlnxedlna



  function coop_cosmology_firstorder_ekappaofa(this, a) result(ekappa)
    COOP_REAL a, ekappa
    class(coop_cosmology_firstorder)::this
    if(a .le. coop_visibility_amin)then
       ekappa = 0.d0
    else
       ekappa = this%ekappa%eval(a)
    endif    
  end function coop_cosmology_firstorder_ekappaofa

  function coop_cosmology_firstorder_visofa(this, a) result(vis)
    COOP_REAL a, vis
    class(coop_cosmology_firstorder)::this
    if(a .le. coop_visibility_amin)then
       vis = 0.d0
    else
       vis = this%vis%eval(a)
    endif    
  end function coop_cosmology_firstorder_visofa


  function coop_cosmology_firstorder_dkappadtau(this, a) result(dkappadtau)
    class(coop_cosmology_firstorder)::this    
    COOP_REAL a, dkappadtau
    dkappadtau =  this%dkappadtau_coef * this%xeofa(a)/a**2
  end function coop_cosmology_firstorder_dkappadtau


  function coop_cosmology_firstorder_taucofa(this, a) result(tauc)
    class(coop_cosmology_firstorder)::this    
    COOP_REAL a, tauc
    tauc =   a**2/(this%dkappadtau_coef * this%xeofa(a))
  end function coop_cosmology_firstorder_taucofa


  function coop_cosmology_firstorder_dot_tauc(this, a) result(dtauc)
    class(coop_cosmology_firstorder)::this    
    COOP_REAL a, dtauc
    dtauc =  a*(2.d0 - this%dlnxedlna(a))/this%xeofa(a)/this%dkappadtau_coef*this%Hasq(a)
  end function coop_cosmology_firstorder_dot_tauc


  function coop_cosmology_firstorder_dkappada(this, a) result(dkappada)
    class(coop_cosmology_firstorder)::this    
    COOP_REAL a, dkappada
    dkappada = this%dkappadtau(a) / this%dadtau(a)
  end function coop_cosmology_firstorder_dkappada

  subroutine coop_cosmology_firstorder_set_source_tau(this, source, step_factor)
    class(coop_cosmology_firstorder_source)::source
    class(coop_cosmology_firstorder)::this
    COOP_REAL step_factor, viscut
    !!basic stepsize
    COOP_REAL,parameter::step_ini = 8.2d-4
    COOP_REAL,parameter::step_min = 6.3d-4
    COOP_REAL,parameter::step_recend = 9.6d-4
    COOP_REAL,parameter::step_reion = 1.3d-3
    COOP_REAL,parameter::step_early = 1.5d-3
    COOP_REAL,parameter::step_late = 2.1d-2
    COOP_REAL,parameter::step_isw = 1.4d-2
    COOP_INT,parameter::nmax = 16384
    COOP_REAL, parameter::a_factor = 1.05d0
    COOP_REAL, parameter::vary = 1.05d0
    COOP_REAL, parameter::slowvary = 1.02d0
    COOP_REAL step
    COOP_REAL a(nmax), tau(nmax), aend, tautmp
    COOP_INT i, n, nw, iw, j
    logical::check_input, do_inds
    n = 1
    aend = 1.d0/(1.d0+this%zrecomb_start)
    a(n) = aend/2.5d0
    tau(n) = this%tauofa(a(n))
    
    step = step_ini

    do while(a(n).lt. aend)
       n =n+1
       if(n.gt.nmax) call coop_return_error("set_source_tau", "nmax<n", "stop")
       a(n) = a(n-1) * a_factor
       tau(n) = tau(n-1) + step*step_factor
       tautmp = this%tauofa(a(n))
       if(tautmp .lt. tau(n))then
          tau(n) = tautmp
       else
          a(n) = this%aoftau(tau(n))
       endif
    enddo

    aend = 1.d0/(1.d0+this%zrecomb)
    do while(a(n).lt. aend)
       n =n+1
       if(n.gt.nmax) call coop_return_error("set_source_tau", "nmax<n", "stop")
       a(n) = a(n-1) * a_factor
       tau(n) = tau(n-1) + step*step_factor
       tautmp = this%tauofa(a(n))
       if(tautmp .lt. tau(n))then
          tau(n) = tautmp
       else
          a(n) = this%aoftau(tau(n))
       endif
       step = max(step_min, step/slowvary)
    enddo

    aend = 1.d0/(1.d0+this%zrecomb_end - 100.d0)
    do while(a(n) .lt. aend)
       n =n+1
       if(n.gt.nmax) call coop_return_error("set_source_tau", "nmax<n", "stop")
       a(n) = a(n-1) * a_factor
       tau(n) = tau(n-1) + step*step_factor
       tautmp = this%tauofa(a(n))
       if(tautmp .lt. tau(n))then
          tau(n) = tautmp
       else
          a(n) = this%aoftau(tau(n))
       endif
       step = min(step_recend, step*slowvary)
    enddo
    if(this%optre.gt.0.01d0)then
       aend = min(1.d0/(1.d0 + this%zre + this%deltaz), 0.5d0)
       do while(a(n) .lt. aend)
          n =n+1
          if(n.gt.nmax) call coop_return_error("set_source_tau", "nmax<n", "stop")
          a(n) = a(n-1) * a_factor
          tau(n) = tau(n-1) + step*step_factor
          tautmp = this%tauofa(a(n))
          if(tautmp .lt. tau(n))then
             tau(n) = tautmp
          else
             a(n) = this%aoftau(tau(n))
          endif
          step = min(step_early, step*slowvary)
       enddo

       aend = min(1.d0/(1.d0 + this%zre), 0.5d0)
       do while(a(n) .lt. aend)
          n =n+1
          if(n.gt.nmax) call coop_return_error("set_source_tau", "nmax<n", "stop")
          a(n) = a(n-1) * a_factor
          tau(n) = tau(n-1) + step*step_factor
          tautmp = this%tauofa(a(n))
          if(tautmp .lt. tau(n))then
             tau(n) = tautmp
          else
             a(n) = this%aoftau(tau(n))
          endif
          step = max(step_reion, step/vary)
       enddo

       aend = min(1.d0/(1.d0 + this%zre-this%deltaz), 0.5d0)
       do while(a(n) .lt. aend)
          n =n+1
          if(n.gt.nmax) call coop_return_error("set_source_tau", "nmax<n", "stop")
          a(n) = a(n-1) * a_factor
          tau(n) = tau(n-1) + step*step_factor
          tautmp = this%tauofa(a(n))
          if(tautmp .lt. tau(n))then
             tau(n) = tautmp
          else
             a(n) = this%aoftau(tau(n))
          endif
          step = min(step_early, step*vary)
       enddo



       do while((1.d0-this%Omega_m) * a(n)**3 .lt. 0.1d0) 
          n =n+1
          if(n.gt.nmax) call coop_return_error("set_source_tau", "nmax<n", "stop")
          tau(n) = tau(n-1) + step*step_factor
          if(tau(n) .lt. 0.989*this%tau0)then
             a(n) = this%aoftau(tau(n))
             step = min(step_late, step*vary)
          else
             tau(n) = 0.99d0*this%tau0
             a(n) = this%aoftau(tau(n))
             exit
          endif
       enddo

       do 
          n =n+1
          if(n.gt.nmax) call coop_return_error("set_source_tau", "nmax<n", "stop")
          tau(n) = tau(n-1) + step*step_factor
          if(tau(n) .lt. 0.99*this%tau0)then
             a(n) = this%aoftau(tau(n))
             step = max(step_isw, step/vary)
          else
             tau(n) = 0.9995d0*this%tau0
             a(n) = this%aoftau(tau(n))
             exit
          endif
       enddo
    endif
    
    source%ntau = n
    if(allocated(source%a))then
       if(size(source%a).ne.n)then
          deallocate(source%a, source%tau, source%chi, source%dtau, source%tauc, source%lna, source%omega_rad, source%vis)
          allocate(source%a(n), source%tau(n), source%dtau(n), source%chi(n), source%tauc(n), source%lna(n), source%omega_rad(n), source%vis(n))
       endif
    else
       allocate(source%a(n),source%tau(n), source%dtau(n), source%chi(n),  source%tauc(n), source%lna(n), source%omega_rad(n), source%vis(n))
    endif
    source%a = a(1:n)
    source%tau = tau(1:n)
    source%dtau(1) = source%tau(2) - source%tau(1)
    do i = 2, n-1
       source%dtau(i) = (source%tau(i+1) - source%tau(i-1))/2.d0
    enddo
    source%dtau(n) = source%tau(n) - source%tau(n-1)
    source%chi = this%tau0 - source%tau
    source%lna = log(source%a)
    !$omp parallel do
    do i=1, n
       source%tauc(i) = this%taucofa(source%a(i))
       source%vis(i) = this%visofa(source%a(i))
       source%omega_rad(i) = this%Omega_r/( this%Omega_r + source%a(i)*(this%Omega_m+(1.d0-this%Omega_m)*source%a(i)**3))  !!this is just an approximation, only used for optimization (drop radiation when omega_rad is very small)       
    enddo
    !$omp end parallel do
    source%index_rad_small = 0
    do i = 1,  n 
       if(source%omega_rad(i) .lt. coop_cosmology_firstorder_omega_rad_cutoff)then
          source%index_rad_small = i
          exit
       endif
    enddo
    source%index_vis_max = coop_maxloc(source%vis)
    source%index_vis_end = source%index_vis_max
    viscut = source%vis(source%index_vis_max)*0.003d0
    do while(source%index_vis_end .lt. n .and. source%vis(source%index_vis_end) .gt. viscut )
       source%index_vis_end = source%index_vis_end + 1
    enddo
    viscut = source%vis(source%index_vis_max)*1.d-4
    source%index_vis_start = source%index_vis_max    
    do while(source%index_vis_start .gt. 1 .and. source%vis(source%index_vis_start) .gt. viscut )
       source%index_vis_start = source%index_vis_start - 1
    enddo

  end subroutine coop_cosmology_firstorder_set_source_tau


  subroutine coop_cosmology_firstorder_set_source_k(this, source, n, weight)
    class(coop_cosmology_firstorder)::this
    class(coop_cosmology_firstorder_source)::source
    COOP_INT n, i, iq, j
    COOP_REAL weight, dlnk
    this%k_pivot = this%kMpc_pivot/this%H0Mpc()
    source%kweight = weight
    source%nk = n
    if(allocated(source%k))then
       if(size(source%k).ne.n)then
          deallocate(source%k, source%dk, source%index_tc_off, source%kop, source%index_rad_off)
          allocate(source%k(n), source%dk(n), source%index_tc_off(n), source%index_rad_off(n), source%kop(n))
       endif
       if(size(source%k_dense, 2).ne.n .or. size(source%k_dense, 1).ne. coop_k_dense_fac )then
          deallocate(source%k_dense, source%ws_dense, source%wt_dense, source%ps_dense)
          allocate(source%k_dense(coop_k_dense_fac, n), source%ws_dense(coop_k_dense_fac, n), source%wt_dense(coop_k_dense_fac, n), source%ps_dense(coop_k_dense_fac, n) )
       endif
    else
       allocate(source%k(n), source%dk(n), source%index_tc_off(n), source%index_rad_off(n), source%kop(n))
       allocate(source%k_dense(coop_k_dense_fac, n), source%ws_dense(coop_k_dense_fac, n), source%wt_dense(coop_k_dense_fac, n), source%ps_dense(coop_k_dense_fac, n) )
    endif
    do i=1, coop_k_dense_fac
       source%a_dense(i) = dble(i)/coop_k_dense_fac
       source%b_dense(i) =  1.d0 - source%a_dense(i)
       source%a2_dense(i) = source%a_dense(i)*(source%a_dense(i)**2-1.d0)
       source%b2_dense(i) = source%b_dense(i)*(source%b_dense(i)**2-1.d0)
    enddo


    source%kmin = 0.2d0/this%distlss
    source%kmax = (min(max(1500, coop_Cls_lmax(source%m)), 3000)*1.7d0)/this%distlss
    call source%k2kop(source%kmin, source%kopmin)
    call source%k2kop(source%kmax, source%kopmax)
    source%dkop = (source%kopmax-source%kopmin)/(n-1)
    source%dkop_dense = source%dkop/coop_k_dense_fac    
    call coop_set_uniform(n, source%kop, source%kopmin, source%kopmax)

    !$omp parallel do private(i, j, dlnk)
    do i=1, n
       call source%kop2k(source%kop(i), source%k(i), source%dkop,  source%dk(i))
       do j=1, coop_k_dense_fac
          call source%kop2k(source%kop(i)+(j-coop_k_dense_fac)*source%dkop_dense, source%k_dense(j,i), source%dkop_dense,  dlnk)
          dlnk = dlnk/source%k_dense(j, i)
          source%ps_dense(j, i) = this%psofk(source%k_dense(j, i))
          source%ws_dense(j, i) = dlnk*source%ps_dense(j, i)
          source%wt_dense(j, i) = this%ptofk(source%k_dense(j, i))*dlnk
       enddo
    enddo
    !$omp end parallel do

    !!set index_tc_off and index_rad_off
    if(source%ntau .gt. 0 .and. allocated(source%tauc))then
       source%index_tc_off = 1
       !$omp parallel do
       do i = 1, n
          do while(source%index_tc_off(i) .lt. source%ntau-1 .and. source%tauc(source%index_tc_off(i)+1)*source%k(i) .le. coop_cosmology_firstorder_tc_cutoff)
             source%index_tc_off(i)= source%index_tc_off(i)+1
          enddo
       enddo
       !$omp end parallel do

       if(coop_firstorder_optimize)then
          source%index_rad_off = max(source%index_rad_small, source%index_tc_off + 1)          
          !$omp parallel do
          do i = 1, n
             do while(source%index_rad_off(i) .le. source%ntau)
                if(source%omega_rad(source%index_rad_off(i)) .gt.  source%k(i)/1.d4)then
                   source%index_rad_off(i)= source%index_rad_off(i) + 1
                else
                   exit
                endif
             enddo
          enddo
          !$omp end parallel do
       else
          source%index_rad_off = source%ntau + 1
       endif
    else
       call coop_return_error("set_source_k", "You need to call set_source_tau before calling set_source_k", "stop")
    endif
    
    !!set index_massivenu_on
    if(this%mnu_by_Tnu*source%a(source%ntau) .lt. coop_pert_massivenu_threshold(1))then
       source%index_massivenu_on = source%ntau + 1
       source%index_massivenu_cold = source%ntau + 1
       return
    endif
    source%index_massivenu_cold = source%ntau + 1
    do while(this%mnu_by_Tnu*source%a(source%index_massivenu_cold-1) .gt. coop_pert_massivenu_cold_threshold)
       source%index_massivenu_cold = source%index_massivenu_cold - 1
       if(source%index_massivenu_cold .le. 1) exit
    enddo

    source%index_massivenu_on(coop_pert_default_nq) = source%index_massivenu_cold
    iq = coop_pert_default_nq
    if(source%index_massivenu_on(iq) .gt. 1)then
       do while(this%mnu_by_Tnu*source%a(source%index_massivenu_on(iq)-1) .gt. coop_pert_massivenu_threshold(iq))
          source%index_massivenu_on(iq) = source%index_massivenu_on(iq)-1
          if(source%index_massivenu_on(iq) .le. 1)exit
       enddo
    endif
    do iq = coop_pert_default_nq - 1, 1, -1
       source%index_massivenu_on(iq) =  source%index_massivenu_on(iq+1)
       if(source%index_massivenu_on(iq) .le. 1)cycle
       do while(this%mnu_by_Tnu*source%a(source%index_massivenu_on(iq)-1) .gt. coop_pert_massivenu_threshold(iq))
          source%index_massivenu_on(iq) = source%index_massivenu_on(iq)-1
          if(source%index_massivenu_on(iq) .le. 1)exit
       enddo
    enddo
  end subroutine coop_cosmology_firstorder_set_source_k



  

  function coop_cosmology_firstorder_source_interpolate(source, k, chi) result(s)
    COOP_INT ik, ichi
    COOP_INT ilow, iup
    class(coop_cosmology_firstorder_source)::source
    COOP_REAL k, chi, kop, rk, a, b, rchi
    COOP_REAL s(source%nsrc)
    if(chi .le. source%chi(source%ntau) .or. chi .ge. source%chi(1))then
       s = 0.d0
       return
    endif
    call source%k2kop(k, kop)
    rk = (kop - source%kopmin)/source%dkop + 1.d0
    ik = floor(rk)
    if(ik.le. 0 .or. ik.ge.source%nk)then
       s = 0.d0
       return
    endif
    a = rk - ik
    b = 1.d0 - a
    ilow = 1
    iup = source%ntau
    do while(iup - ilow .gt. 1)
       ichi = (iup + ilow)/2
       if(source%chi(ichi) .gt. chi)then
          ilow = ichi
       else
          iup = ichi
       endif
    enddo
    rchi = source%chi(ilow) - chi
    s = ((source%s(:, ik, ilow)+source%s2(:, ik, ilow)*(b**2-1.d0))*b + (source%s(:, ik+1, ilow) + source%s2(:, ik+1, ilow)*(a**2-1.d0))*a)*(1.d0-rchi) &
         +  ((source%s(:, ik, iup)+source%s2(:, ik, iup)*(b**2-1.d0))*b + (source%s(:, ik+1, iup) + source%s2(:, ik+1, iup)*(a**2-1.d0))*a)*rchi
  end function coop_cosmology_firstorder_source_interpolate


  function coop_cosmology_firstorder_source_interpolate_one(source, k, chi, isrc) result(s)
    COOP_INT ik, ichi, isrc
    COOP_INT ilow, iup
    class(coop_cosmology_firstorder_source)::source
    COOP_REAL k, chi, kop, rk, a, b, rchi
    COOP_REAL s
    if(chi .le. source%chi(source%ntau) .or. chi .ge. source%chi(1))then
       s = 0.d0
       return
    endif
    call source%k2kop(k, kop)
    rk = (kop - source%kopmin)/source%dkop + 1.d0
    ik = floor(rk)
    if(ik.le. 0 .or. ik.ge.source%nk)then
       s = 0.d0
       return
    endif
    a = rk - ik
    b = 1.d0 - a
    ilow = 1
    iup = source%ntau
    do while(iup - ilow .gt. 1)
       ichi = (iup + ilow)/2
       if(source%chi(ichi) .gt. chi)then
          ilow = ichi
       else
          iup = ichi
       endif
    enddo
    rchi = source%chi(ilow) - chi
    s = ((source%s(isrc, ik, ilow)+source%s2(isrc, ik, ilow)*(b**2-1.d0))*b + (source%s(isrc, ik+1, ilow) + source%s2(isrc, ik+1, ilow)*(a**2-1.d0))*a)*(1.d0-rchi) &
         +  ((source%s(isrc, ik, iup)+source%s2(isrc, ik, iup)*(b**2-1.d0))*b + (source%s(isrc, ik+1, iup) + source%s2(isrc, ik+1, iup)*(a**2-1.d0))*a)*rchi
  end function coop_cosmology_firstorder_source_interpolate_one

  subroutine coop_cosmology_firstorder_source_get_Psi_trans(source, tau, nk, k, psi)
    class(coop_cosmology_firstorder_source) source
    COOP_INT nk
    COOP_REAL k(nk), tau, psi(nk)
    COOP_REAL kop, a, b, atau, btau
    COOP_INT ik, itau, im, it, ikop
    if(tau .le. source%tau(1))then
       itau = 1
       atau = 0.d0
    elseif(tau .ge. source%tau(source%ntau))then
       itau = source%ntau - 1
       atau = (tau - source%tau(itau))/(source%tau(itau+1)-source%tau(itau))
    else
       itau = 1
       it = source%ntau
       do while(it - itau .gt. 1)
          im = (itau+it)/2
          if(source%tau(im) .gt. tau)then
             it = im
          else
             itau = im
          endif
       enddo
       atau = (tau - source%tau(itau))/(source%tau(itau+1)-source%tau(itau))
    endif
    btau = 1.d0-atau
    !$omp parallel do private(ikop, kop, a, b)
    do ik=1, nk
       call source%k2kop(k(ik), kop)
       a = (kop - source%kopmin)/source%dkop + 1.d0
       ikop = floor(a)
       if(ikop .lt. 1)then
          psi(ik) = source%saux(3, 1, itau)*btau +  source%saux(3, 1, itau+1)*atau
       elseif(ikop .ge. source%nk)then
          psi(ik) = (source%saux(3, source%nk, itau)*btau +  source%saux(3, source%nk, itau+1)*atau) * coop_bbks_trans(k(ik)/source%bbks_keq)/source%bbks_trans_kmax

       else
          a = a - ikop
          b = 1.d0 - a
          psi(ik) = ((source%saux(3, ikop, itau)+source%saux(5, ikop, itau)*(1.d0-b**2))*b + (source%saux(3, ikop+1, itau) + source%saux(5, ikop+1, itau)*(1.d0-a**2) ) * a ) * btau +  ((source%saux(3, ikop, itau+1)+source%saux(5, ikop, itau+1)*(1.d0-b**2))*b + (source%saux(3, ikop+1, itau+1) + source%saux(5, ikop+1, itau+1)*(1.d0-a**2) ) * a ) * atau
       endif
    enddo
    !$omp end parallel do
    
  end subroutine coop_cosmology_firstorder_source_get_Psi_trans


  function coop_bbks_trans(x) result(bbkstrans)
    !!BBKS fitting formula of adiabatic CDM transfer function, x = k/k_eq
    COOP_REAL::c1 = 0.284d0
    COOP_REAL::c2 = 1.18d0**2
    COOP_REAL::c3 = 0.399d0**3
    COOP_REAL::c4 = 0.490d0**4
    COOP_REAL BBKSTrans, x
    BBKSTrans = log(1.+0.171*x)/(0.171*x)/ &
         ( 1.d0+ x*(c1+x*(c2+x*(c3+x*c4)))) ** 0.25d0
  end function coop_bbks_trans

 subroutine coop_cosmology_firstorder_source_get_PhiPlusPsi_trans(source, tau, nk, k, PhiPlusPsi)
    class(coop_cosmology_firstorder_source) source
    COOP_INT nk
    COOP_REAL k(nk), tau, PhiPlusPsi(nk)
    COOP_REAL kop, a, b, atau, btau
    COOP_INT ik, itau, im, it, ikop
    if(tau .le. source%tau(1))then
       itau = 1
       atau = 0.d0
    elseif(tau .ge. source%tau(source%ntau))then
       itau = source%ntau - 1
       atau = (tau - source%tau(itau))/(source%tau(itau+1)-source%tau(itau))
    else
       itau = 1
       it = source%ntau
       do while(it - itau .gt. 1)
          im = (itau+it)/2
          if(source%tau(im) .gt. tau)then
             it = im
          else
             itau = im
          endif
       enddo
       atau = (tau - source%tau(itau))/(source%tau(itau+1)-source%tau(itau))
    endif
    btau = 1.d0-atau
    !$omp parallel do private(ikop, kop, a, b, ik)
    do ik=1, nk
       call source%k2kop(k(ik), kop)
       a = (kop - source%kopmin)/source%dkop + 1.d0
       ikop = floor(a)
       if(ikop .lt. 1)then
          PhiPlusPsi(ik) = source%saux(2, 1, itau)*btau +  source%saux(2, 1, itau+1)*atau
       elseif(ikop .ge. source%nk)then
          PhiPlusPsi(ik) = (source%saux(2, source%nk, itau)*btau +  source%saux(2, source%nk, itau+1)*atau) * coop_bbks_trans(k(ik)/source%bbks_keq)/source%bbks_trans_kmax
       else
          a = a - ikop
          b = 1.d0 - a
          PhiPlusPsi(ik) = ((source%saux(2, ikop, itau)+source%saux(4, ikop, itau)*(1.d0-b**2))*b + (source%saux(2, ikop+1, itau) + source%saux(4, ikop+1, itau)*(1.d0-a**2) ) * a ) * btau +  ((source%saux(2, ikop, itau+1)+source%saux(4, ikop, itau+1)*(1.d0-b**2))*b + (source%saux(2, ikop+1, itau+1) + source%saux(4, ikop+1, itau+1)*(1.d0-a**2) ) * a ) * atau
       endif
    enddo
    !$omp end parallel do
  end subroutine coop_cosmology_firstorder_source_get_PhiPlusPsi_trans


!!return the dimensionless k^3 |\delta_k|^2/(2pi^2)
  subroutine coop_cosmology_firstorder_get_matter_power(this, z, nk, k, Pk)
    class(coop_cosmology_firstorder)::this
    COOP_INT nk, ik
    COOP_REAL z, a, k(nk), Pk(nk), tau, Psi(nk), Ps(nk), PhiPlusPsi(nk)
    a = 1.d0/(1.d0+z)
    tau = this%tauofa(a)
    call this%source(0)%get_Psi_trans(tau, nk, k, Psi)
    call this%source(0)%get_PhiPlusPsi_trans(tau, nk, k, PhiPlusPsi)
    !$omp parallel do
    do ik = 1, nk
       ps(ik) = this%psofk(k(ik))
    enddo
    !$omp end parallel do
    pk = (PhiPlusPsi-Psi) **2 * ps * (2.d0*k**2/(O0_BARYON(this)%density(a) + O0_CDM(this)%density(a)))**2
  end subroutine coop_cosmology_firstorder_get_matter_power

#include "firstorder_basic_utils.h"
#include "firstorder_compute_sigma.h"  

  
end module coop_firstorder_mod

