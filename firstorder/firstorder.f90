module coop_firstorder_mod
  use coop_wrapper_background
  use coop_pertobj_mod
  implicit none
#include "constants.h"

  private

  public::coop_cosmology_firstorder, coop_cosmology_firstorder_source,  coop_recfast_get_xe, coop_power_lnk_min, coop_power_lnk_max,  coop_k_dense_fac, coop_index_ClTT, coop_index_ClTE, coop_index_ClEE, coop_index_ClBB, coop_index_ClLenLen, coop_index_ClTLen,  coop_num_Cls, coop_Cls_lmax, coop_bbks_trans, coop_index_source_T, coop_index_source_E, coop_index_source_B, coop_index_source_Len, coop_next_l, coop_nl_to_lmax, coop_set_ells
#if DO_ZETA_TRANS  
  public::coop_index_source_zeta, coop_index_ClZetaZeta, coop_index_ClTZeta, coop_index_ClEZeta, coop_zeta_single_slice_chi, coop_zeta_user_specified_weight
#endif  

  !!this makes the code faster and more accurate
  logical,parameter :: coop_firstorder_optimize = .true.
  COOP_INT, parameter :: coop_limber_ell = 550

  
  COOP_INT :: coop_Cls_lmax(0:2) = (/ 2800, 2000, 1500 /)

  COOP_REAL, parameter :: coop_power_lnk_min = log(0.1d0) 
  COOP_REAL, parameter :: coop_power_lnk_max = log(5.d3) 
  COOP_REAL, parameter :: coop_visibility_amin = 1.8d-4
  COOP_REAL, parameter :: coop_initial_condition_epsilon = 2.d-7
  COOP_REAL, parameter :: coop_cosmology_firstorder_ode_accuracy = 2.d-7
  COOP_REAL, parameter :: coop_cosmology_firstorder_tc_cutoff = 0.01d0


  COOP_REAL, dimension(0:2), parameter :: coop_source_tau_step_factor = (/ 1.d0, 1.d0, 1.d0 /)
  COOP_REAL, dimension(0:2), parameter :: coop_source_k_weight = (/ 0.15d0, 0.15d0, 0.1d0 /)
  COOP_INT, dimension(0:2), parameter :: coop_source_k_n = (/ 300, 160, 160 /)
  COOP_REAL, parameter :: coop_source_k_index = 0.88d0


  COOP_INT, parameter :: coop_index_source_T = 1
  COOP_INT, parameter :: coop_index_source_E = 2
  COOP_INT, parameter :: coop_index_source_B = 3
#if DO_ZETA_TRANS
  COOP_INT, parameter :: coop_k_dense_fac = 25
  COOP_INT, parameter :: coop_index_source_zeta = 3
  COOP_INT, parameter :: coop_index_source_Len = 4
  COOP_REAL::coop_zeta_single_slice_chi = -1.d0  !!if set to negative, weight = visibility function
  type(coop_function)::coop_zeta_user_specified_weight
#else
  COOP_INT, parameter :: coop_k_dense_fac = 15
  COOP_INT, parameter :: coop_index_source_Len = 3
#endif

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
  COOP_INT, dimension(0:2), parameter::coop_num_sources = (/ 4,  3,  3 /)
#else
  COOP_INT, parameter::coop_num_Cls =  coop_index_ClTLen
  COOP_INT, dimension(0:2), parameter::coop_num_sources = (/ 3,  3,  3 /)
#endif
  COOP_INT, dimension(0:2), parameter::coop_num_saux = (/ 5,  1,  1 /)

!!recfast head file
#include "recfast_head.h"


  type coop_cosmology_firstorder_source
     COOP_INT::m = 0
     COOP_INT::ntau = 0
     COOP_INT::nk = 0
     COOP_INT::nsrc = 0
     COOP_INT::nsaux = 0
     COOP_INT, dimension(:),allocatable::index_tc_off
     COOP_INT::index_de_perturb_on
     COOP_INT, dimension(coop_pert_default_nq)::index_massivenu_on
     COOP_INT::index_massivenu_cold, index_vis_end, index_vis_max, index_vis_start
     COOP_REAL::dkop, kopmin, kopmax, kmin, kmax, kweight, tauweight, bbks_keq, bbks_trans_kmax, distlss, dkop_dense
     COOP_REAL,dimension(coop_k_dense_fac)::a_dense, b_dense, a2_dense, b2_dense
     COOP_REAL, dimension(:),allocatable::k, kop, dk !!tau is conformal time, chi is comoving distance; in flat case, chi + tau = tau_0
     COOP_REAL, dimension(:),allocatable::tau, a, tauc, lna, dtau, chi, omega_rad, vis, omega_de
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
     COOP_REAL::zrecomb, distlss, tau0, zrecomb_start, maxvis, taurecomb, arecomb, zrecomb_end, arecomb_start, z_drag, z_star, r_drag, r_star, tau_late
     COOP_REAL::optre = 0.07d0
     COOP_REAL::zre = 8.d0
     COOP_REAL::deltaz = 1.5d0
     COOP_REAL::sigma_8 = 0.8d0
     COOP_REAL::r = 0.d0
     COOP_REAL::As = 2.2d-9
     COOP_REAL::ns = 0.96d0
     COOP_REAL::nrun = 0.d0
     COOP_REAL::nt = 0.d0
     COOP_REAL::kMpc_pivot = 0.05d0
     COOP_INT ::de_genre = COOP_PERT_NONE
     COOP_REAL::k_pivot
     
     !!these omega's are defined at a=1 (today)     
     COOP_REAL::dkappadtau_coef, ReionFrac, Omega_b, Rbya, Omega_c, Omega_nu, Omega_g, tau_eq, mnu_by_Tnu,  Omega_massivenu, bbks_keq

     !!these two parameters are defined in the a->0 limit (only matters when baryon or CDM is coupled to DE in some modified gravity models)     
     COOP_REAL::ombh2 = 0.022d0
     COOP_REAL::omch2 = 0.12d0
     
     logical::inflation_consistency = .true.
     
     type(coop_function)::Ps, Pt, Xe, ekappa, vis, Tb
     type(coop_cosmology_firstorder_source),dimension(0:2)::source
     COOP_INT::index_baryon, index_cdm, index_radiation, index_nu, index_massiveNu, index_de
     COOP_REAL, dimension(0:coop_pert_default_lmax, 0:coop_pert_default_mmax, 0:coop_pert_default_smax)::klms, klms_by_2lm1, klms_by_2lp1
     COOP_REAL,dimension(coop_pert_default_lmax)::fourbyllp1
     logical::klms_done = .false.
     logical::has_tensor = .false.
   contains
     procedure:: set_standard_cosmology =>  coop_cosmology_firstorder_set_standard_cosmology
#if DO_EFT_DE     
     procedure:: set_EFT_cosmology =>  coop_cosmology_firstorder_set_EFT_cosmology
#endif
#if DO_COUPLED_DE
     procedure:: set_coupled_DE_cosmology =>  coop_cosmology_firstorder_set_coupled_DE_cosmology
#endif     
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
     procedure:: get_z_drag => coop_cosmology_firstorder_get_z_drag
     procedure:: get_z_star => coop_cosmology_firstorder_get_z_star     
     procedure:: dlnxedlna => coop_cosmology_firstorder_dlnxedlna
     procedure:: dkappadz => coop_cosmology_firstorder_dkappadz    !!(z)
     procedure:: doptdragdz => coop_cosmology_firstorder_doptdragdz    !!(z)
     procedure:: dsoundda => coop_cosmology_firstorder_dsoundda !!(a)
     procedure:: dkappada => coop_cosmology_firstorder_dkappada    !!(a)
     procedure:: dkappadtau => coop_cosmology_firstorder_dkappadtau  !!(a)
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
     procedure:: get_Mphi_power => coop_cosmology_firstorder_get_Mphi_power     
     procedure:: sigma_Tophat_R => coop_cosmology_firstorder_sigma_Tophat_R
     procedure:: sigma_Gaussian_R => coop_cosmology_firstorder_sigma_Gaussian_R
     procedure:: sigma_Gaussian_R_quick => coop_cosmology_firstorder_sigma_Gaussian_R_quick
     procedure::sigma_Gaussian_R_with_dervs => coop_cosmology_firstorder_sigma_Gaussian_R_with_dervs

     procedure::growth_of_z => coop_cosmology_firstorder_growth_of_z
     procedure::late_damp_factor => coop_cosmology_firstorder_late_damp_factor
     !!
#if DO_ZETA_TRANS
     procedure::set_zeta_weight_latevis => coop_cosmology_firstorder_set_zeta_weight_latevis
     procedure::set_zeta_weight_earlyvis => coop_cosmology_firstorder_set_zeta_weight_earlyvis
     procedure::set_zeta_weight_vis => coop_cosmology_firstorder_set_zeta_weight_vis
     procedure::set_zeta_weight_single_slice => coop_cosmology_firstorder_set_zeta_weight_single_slice          
#endif     
     
  end type coop_cosmology_firstorder


contains


!!recfast code
#include "recfast_source.h"

  !!this head file contains the evolution equations of the firstorder ODE system
#if DO_EFT_DE
#include "firstorder_equations_EFT.h"  
#elif DO_COUPLED_DE
#include "firstorder_equations_CPLDE.h"  
#else   
#include "firstorder_equations.h"
#endif
  
!!this head file set the outputs (saved in this%source%s)
#include "firstorder_source.h"

!!this head file sets the initial conditions
#if DO_EFT_DE
#include "firstorder_ic_EFT.h"
#elif DO_COUPLE_DE
#include "firstorder_ic_CPLDE.h"    
#else  
#include "firstorder_ic.h"
#endif
  
  subroutine coop_cosmology_firstorder_set_Planck_Bestfit(this, Omega_nu)
    class(coop_cosmology_firstorder)::this
    COOP_REAL, optional::Omega_nu
    if(present(Omega_nu))then
       call this%set_standard_cosmology(Omega_b=0.04917d0, Omega_c=0.2647d0, h = 0.6727d0, tau_re = 0.079d0, As = 2.206d-9, ns = 0.9645d0,  Omega_nu = Omega_nu)
    else
       call this%set_standard_cosmology(Omega_b=0.04917d0, Omega_c=0.2647d0, h = 0.6727d0, tau_re = 0.079d0, As = 2.206d-9, ns = 0.9645d0)
     
    endif
  end subroutine coop_cosmology_firstorder_set_Planck_Bestfit


  subroutine coop_cosmology_firstorder_set_Planck_Bestfit_with_r(this, r, Omega_nu)
    class(coop_cosmology_firstorder)::this
    COOP_REAL r
    COOP_REAL, optional::Omega_nu
    if(present(Omega_nu))then
       call this%set_standard_cosmology(Omega_b=0.04917d0, Omega_c=0.2647d0, h = 0.6727d0, tau_re = 0.079d0, As = 2.206d-9, ns = 0.9645d0,  Omega_nu = Omega_nu, r = r)
    else
       call this%set_standard_cosmology(Omega_b=0.04917d0, Omega_c=0.2647d0, h = 0.6727d0, tau_re = 0.079d0, As = 2.206d-9, ns = 0.9645d0, r = r)       
     
    endif
  end subroutine coop_cosmology_firstorder_set_Planck_Bestfit_with_r



  subroutine coop_cosmology_firstorder_source_get_Cls(source, l, Cls)
    class(coop_cosmology_firstorder_source)::source
    COOP_INT l
    COOP_REAL,dimension(coop_num_cls),intent(OUT)::Cls
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
       Cls(coop_index_ClLenLen) =  sum(source%ws_dense * trans(coop_index_source_Len, :, :)**2)*coop_4pi
       Cls(coop_index_ClTLen) =  sum(source%ws_dense * trans(coop_index_source_T, :, :) * trans(coop_index_source_Len,:,:))*coop_4pi
#if DO_ZETA_TRANS        
       Cls(coop_index_ClTzeta) = sum(source%ws_dense * trans(coop_index_source_T, :, :)*trans(coop_index_source_zeta,:,:))*coop_4pi
       Cls(coop_index_ClEzeta) = sqrt((l+2.d0)*(l+1.d0)*l*(l-1.d0))*sum(source%ws_dense * trans(coop_index_source_E, :, :) * trans(coop_index_source_zeta,:,:))*coop_4pi
       Cls(coop_index_Clzetazeta) = sum(source%ws_dense * trans(coop_index_source_zeta, :, :)**2)*coop_4pi
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
       call coop_next_l(l)
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
       call coop_next_l(l)
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
    enddo
    !$omp end parallel do
100 deallocate(ls_computed, Cls_computed,Cls2_computed)

  end subroutine coop_cosmology_firstorder_source_get_All_Cls

  subroutine coop_next_l(l)
    COOP_INT::l
    l = l + min(60, 12 + l/30, max(1, l/3))
  end subroutine coop_next_l

  function coop_nl_to_lmax(lmax) result(nl)
    COOP_INT::l, lmax, nl
    l = 2
    nl = 1
    do while(l .lt. lmax)
       nl = nl + 1
       call coop_next_l(l)
    enddo
  end function coop_nl_to_lmax

  subroutine coop_set_ells(ells, lmax)
    COOP_INT::nl, lmax, l
    COOP_INT::ells(:)
    l = 2
    nl = 0
    do while(l .lt. lmax)
       nl = nl + 1
       ells(nl) = l
       call coop_next_l(l)
    enddo
    nl = nl + 1
    ells(nl) = lmax
  end subroutine coop_set_ells


  subroutine coop_cosmology_firstorder_compute_source(this, m, success)
    class(coop_cosmology_firstorder)::this
    COOP_INT m, ik, itau, is
    logical, optional::success
    logical,dimension(:),allocatable::suc
    call this%init_source(m)
    this%source(m)%bbks_keq = this%bbks_keq
    this%source(m)%bbks_trans_kmax = coop_bbks_trans(this%source(m)%kmax/this%bbks_keq)
    if(present(success))then
       allocate(suc(this%source(m)%nk))
       suc = .true.
       !$omp parallel do
       do ik = 1, this%source(m)%nk
          if(all(suc))call this%compute_source_k(this%source(m), ik, success = suc(ik))
       enddo
       !$omp end parallel do
       success = all(suc)
       deallocate(suc)
       if(.not. success) return
    else
       !$omp parallel do
       do ik = 1, this%source(m)%nk
          call this%compute_source_k(this%source(m), ik)
       enddo
       !$omp end parallel do
    endif

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


  subroutine coop_cosmology_firstorder_compute_source_k(this, source, ik, do_output, transfer_only, success)
    COOP_REAL, parameter::eps = 1.d-20
    class(coop_cosmology_firstorder)::this
    type(coop_cosmology_firstorder_source)::source
    COOP_INT ik, nvars, itau, iq, scheme
    type(coop_pert_object) pert
    COOP_REAL, dimension(:,:),allocatable::w
    logical,optional::do_output, transfer_only, success
    COOP_REAL c(24)
    COOP_INT ind, i
    COOP_REAL tau_ini, lna, mnu_deltarho, mnu_deltav
    tau_ini = min(coop_initial_condition_epsilon/source%k(ik), this%conformal_time(this%a_eq*coop_initial_condition_epsilon), source%tau(1)*0.999d0)
    call this%set_initial_conditions(pert, m = source%m, k = source%k(ik), tau = tau_ini)
    lna = log(this%aoftau(tau_ini))
#if DO_EFT_DE    
    pert%de_scheme = 0
#endif
    select type(this)
    type is(coop_cosmology_firstorder)
       call coop_cosmology_firstorder_equations(pert%ny+1, lna, pert%y, pert%yp, this, pert)
    class default
       stop "For compatibility with lower versions of gfortran, firstorder equations only works with type coop_cosmology_firstorder"
    end select
    
    !!for energy conservation test:
    if(present(do_output))then
       if(do_output)then
          call pert%print()
       endif
    endif
    
    ind = 1
    c = 0.d0
    nvars = pert%ny + 1
    allocate(w(nvars, 9))
    w = 0.d0
    iq = 1
    if(present(success))then
       success = .true.
    endif
    do itau = 1, source%ntau
       select type(this)
       type is(coop_cosmology_firstorder)
       
          call coop_dverk_firstorder(nvars, coop_cosmology_firstorder_equations, this, pert, lna,   pert%y, source%lna(itau),  coop_cosmology_firstorder_ode_accuracy, ind, c, nvars, w)
          if(present(success) .and. ind .le. 0 .or. ind .gt.6)then
             success = .false.
             return
          endif
       class default
          stop "For compatibility with lower versions of gfortran, dverk_firstorder  only works with type coop_cosmology_firstorder"
       end select
       
       pert%want_source = .true.
       select type(this)
       type is(coop_cosmology_firstorder)
          call coop_cosmology_firstorder_equations(pert%ny+1, lna, pert%y, pert%yp, this, pert)
       class default
          stop "For compatibility with lower versions of gfortran, firstorder equations only works with type coop_cosmology_firstorder"
       end select
       if(present(success))then
          if(.not. (all(abs(pert%y).lt. 1.d20) .and. all(abs(pert%yp) .lt. 1.d20)))then
             success = .false.
             return
          endif
       endif
       pert%want_source  = .false.              
       !!for energy conservation test:
       if(present(do_output))then
          if(do_output)then
             call pert%print()
          endif
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
#if DO_EFT_DE
       if(itau .ge. source%index_de_perturb_on)then
          if(abs(pert%deMat(2, 4)).gt. eps)then
             scheme = 3
          else
             if(abs(pert%deMat(5, 4)) .gt. eps)then
                scheme = 2
             else
                if(abs(pert%deMat(6, 4)) .gt. eps)then
                   scheme = 1
                else
                   scheme = 0
                endif
             endif
          endif
          if(scheme .ne. pert%de_scheme)then
             ind = 1
             pert%de_scheme = scheme
          endif
       endif
#endif       

       if(itau .eq. source%index_tc_off(ik))then
          if(coop_feedback_level .ge. 5) write(*,*) "switching tight coupling off"
          call pert%save_ode()
       
          select type(this)
          type is(coop_cosmology_firstorder)          
             call coop_cosmology_firstorder_equations(pert%ny+1, lna, pert%y, pert%yp, this, pert)   !!set tight coupling apprixmations
          class default
             stop "For compatibility with lower versions of gfortran, firstorder equations only works with type coop_cosmology_firstorder"
          end select
             
          pert%tight_coupling = .false.
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


  function coop_cosmology_firstorder_doptdragdz(this, z) result(ddragdz)
    class(coop_cosmology_firstorder)::this
    COOP_REAL z, ddragdz, a
    ddragdz = this%dkappadz(z)*(1.d0+z)/this%Rbya
  end function coop_cosmology_firstorder_doptdragdz

  
  function coop_cosmology_firstorder_dkappadz(this, z) result(dkappadz)
    class(coop_cosmology_firstorder)::this
    COOP_REAL z, dkappadz, a
    a = 1.d0/(1.d0+z)
    dkappadz = this%dkappada(a)*a**2
  end function coop_cosmology_firstorder_dkappadz


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


  !!input:
  !!source: computed source
  !!tau: conformal time
  !!k(1:nk): k array
  !!output:
  !!phiLen(1:nk)
  subroutine coop_cosmology_firstorder_source_get_Len_trans(source, tau, nk, k, phiLen)
    class(coop_cosmology_firstorder_source) source
    COOP_INT nk
    COOP_REAL k(nk), tau, phiLen(nk)
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
          phiLen(ik) = source%s(coop_index_source_Len, 1, itau)*btau +  source%s(coop_index_source_Len, 1, itau+1)*atau
       elseif(ikop .ge. source%nk)then
          phiLen(ik) = source%s(coop_index_source_Len, source%nk, itau)*btau +  source%s(coop_index_source_Len, source%nk, itau+1)*atau
       else
          a = a - ikop
          b = 1.d0 - a
          phiLen(ik) = ((source%s(coop_index_source_Len, ikop, itau)+source%s2(coop_index_source_Len, ikop, itau)*(1.d0-b**2))*b + (source%s(coop_index_source_Len, ikop+1, itau) + source%s2(coop_index_source_Len, ikop+1, itau)*(1.d0-a**2) ) * a ) * btau +  ((source%s(coop_index_source_Len, ikop, itau+1)+source%s2(coop_index_source_Len, ikop, itau+1)*(1.d0-b**2))*b + (source%s(coop_index_source_Len, ikop+1, itau+1) + source%s2(coop_index_source_Len, ikop+1, itau+1)*(1.d0-a**2) ) * a ) * atau
       endif
    enddo
    !$omp end parallel do
  end subroutine coop_cosmology_firstorder_source_get_Len_trans
  

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
#if DO_EFT_DE
    pk = (PhiPlusPsi-Psi) **2 * ps * (2.d0*k**2/(O0_BARYON(this)%density(a) + O0_CDM(this)%density(a))/this%M2(a) )**2    
#else    
    pk = (PhiPlusPsi-Psi) **2 * ps * (2.d0*k**2/(O0_BARYON(this)%density(a) + O0_CDM(this)%density(a)))**2
#endif
  end subroutine coop_cosmology_firstorder_get_matter_power


  subroutine coop_cosmology_firstorder_get_MPhi_power(this, z, nk, k, Pk)
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
    pk = (PhiPlusPsi-Psi) **2 * ps * (2.d0*k**2/(3.d0*(this%omch2+this%ombh2)/this%h()**2/a**3))**2
  end subroutine coop_cosmology_firstorder_get_MPhi_power


  function coop_cosmology_firstorder_dsoundda(this, a) result(dsda)  !!use approximations, to be consistent with CosmoMC
    class(coop_cosmology_firstorder)::this
    COOP_REAL a, dsda, R
    dsda = (1.d0/coop_sqrt3) / sqrt(1.d0 +  this%Rbya *  a )/ this%Hasq(a)
  end function coop_cosmology_firstorder_dsoundda


  !! z_drag 
  subroutine coop_cosmology_firstorder_get_z_drag(this, z_drag)
    class(coop_cosmology_firstorder)::this
    COOP_REAL z_drag, intopt, dz, z, dilow, diup, dimid
    intopt = 0.d0
    z = this%zre + this%deltaz * 5.d0
    dilow = this%doptdragdz(z)
    dz = 1.d-2
    do while(intopt .lt. 6.d0)
       dimid = this%doptdragdz(z+dz/2.d0)
       z = z + dz
       diup = this%doptdragdz(z)       
       intopt = intopt + (diup + dimid*4.d0 + dilow)*dz
       dilow = diup
       if(z .gt. 1.d4) stop "z_drag does not converge"
    enddo
    dz = dz/10.d0
    do while(intopt .gt. 6.d0)
       dimid = this%doptdragdz(z - dz/2.d0)
       z = z - dz
       diup = this%doptdragdz(z)       
       intopt = intopt - (diup + dimid*4.d0 + dilow)*dz
       dilow = diup
    enddo
    if(diup .gt. 1.d-10)then
       z_drag = min(z + (1.d0-intopt/6.d0)/diup, z+dz)
    else
       z_drag = z
    endif
  end subroutine coop_cosmology_firstorder_get_z_drag

  subroutine coop_cosmology_firstorder_get_z_star(this, z_star)
    class(coop_cosmology_firstorder)::this
    COOP_REAL z_star, intopt, dz, z, dilow, diup, dimid
    intopt = 0.d0
    z = this%zre + this%deltaz * 5.d0
    dilow = this%dkappadz(z)
    dz = 1.d-2
    do while(intopt .lt. 6.d0)
       dimid = this%dkappadz(z+dz/2.d0)
       z = z + dz
       diup = this%dkappadz(z)       
       intopt = intopt + (diup + dimid*4.d0 + dilow)*dz
       dilow = diup
       if(z .gt. 1.d4) stop "z_star does not converge"
    enddo
    dz = dz/10.d0
    do while(intopt .gt. 6.d0)
       dimid = this%dkappadz(z - dz/2.d0)
       z = z - dz
       diup = this%dkappadz(z)       
       intopt = intopt - (diup + dimid*4.d0 + dilow)*dz
       dilow = diup
    enddo
    if(diup .gt. 1.d-10)then
       z_star = min(z + (1.d0-intopt/6.d0)/diup, z+dz)
    else
       z_star = z
    endif
  end subroutine coop_cosmology_firstorder_get_z_star
  
    
#include "firstorder_basic_utils.h"
#include "firstorder_compute_sigma.h"  


!!return phi(z)/phi_matter_dominate  
  function coop_cosmology_firstorder_growth_of_z(this, z, k) result(Dz)
    class(coop_cosmology_firstorder)::this
    COOP_REAL  atau, btau, tau
    COOP_INT itau, im, it, ik
    COOP_REAL rk, kop
    COOP_REAL::z, Dz
    COOP_REAL,optional::k
    COOP_REAL, parameter::omr_zero = 1.d-3
    COOP_INT::iz_ref
    iz_ref = this%source(0)%ntau
    do while( iz_ref .gt. 2 .and. this%source(0)%omega_rad(iz_ref-1).lt.omr_zero)
       iz_ref = iz_ref - 1
    enddo
    tau = this%tauofa(1.d0/(1.d0+z))
    if(tau .le. this%source(0)%tau(iz_ref))then
       itau = iz_ref
       atau = 0.d0
    elseif(tau .ge. this%source(0)%tau(this%source(0)%ntau))then
       itau = this%source(0)%ntau - 1
       atau = (tau - this%source(0)%tau(itau))/(this%source(0)%tau(itau+1)-this%source(0)%tau(itau))
    else
       itau = 1
       it = this%source(0)%ntau
       do while(it - itau .gt. 1)
          im = (itau+it)/2
          if(this%source(0)%tau(im) .gt. tau)then
             it = im
          else
             itau = im
          endif
       enddo
       atau = (tau - this%source(0)%tau(itau))/(this%source(0)%tau(itau+1)-this%source(0)%tau(itau))
    endif
    btau = 1.d0-atau
    if(present(k))then
       call this%source(0)%k2kop(k, kop)
       rk = (kop - this%source(0)%kopmin)/this%source(0)%dkop + 1.d0
       ik = floor(rk)
       rk = rk - ik
       if(ik .lt.1)then
          Dz = ((this%source(0)%saux(2, 1, itau) - this%source(0)%saux(3, 1, itau))*btau +  (this%source(0)%saux(2, 1, itau+1) - this%source(0)%saux(3, 1, itau+1))*atau)/  (this%source(0)%saux(2, 1, iz_ref) - this%source(0)%saux(3, 1, iz_ref))
          return
       elseif(ik .ge. this%source(0)%nk)then
          ik = this%source(0)%nk
          Dz = ((this%source(0)%saux(2, ik, itau) - this%source(0)%saux(3, ik, itau))*btau +  (this%source(0)%saux(2, ik, itau+1) - this%source(0)%saux(3, ik, itau+1))*atau)/ (this%source(0)%saux(2, ik, iz_ref) - this%source(0)%saux(3, ik, iz_ref))
          return
       else
          Dz =  ( &
               ((this%source(0)%saux(2, ik, itau) - this%source(0)%saux(3, ik, itau))*btau +  (this%source(0)%saux(2, ik, itau+1) - this%source(0)%saux(3, ik, itau+1))*atau)  * (1.d0 -rk) &
               +((this%source(0)%saux(2, ik+1, itau) - this%source(0)%saux(3, ik+1, itau))*btau +  (this%source(0)%saux(2, ik+1, itau+1) - this%source(0)%saux(3, ik+1, itau+1))*atau)  * rk ) &
               / ( (this%source(0)%saux(2, ik, iz_ref) - this%source(0)%saux(3, ik, iz_ref))*(1.d0-rk) &
               +   (this%source(0)%saux(2, ik+1, iz_ref) - this%source(0)%saux(3, ik+1, iz_ref))*rk )
          
       endif
    else !!k -> 0 limit
       Dz = ((this%source(0)%saux(2, 1, itau) - this%source(0)%saux(3, 1, itau))*btau +  (this%source(0)%saux(2, 1, itau+1) - this%source(0)%saux(3, 1, itau+1))*atau)/ (this%source(0)%saux(2, 1, iz_ref) - this%source(0)%saux(3, 1, iz_ref))
    endif
  end function coop_cosmology_firstorder_growth_of_z


  subroutine coop_cosmology_firstorder_set_standard_cosmology(this, h, omega_b, omega_c, tau_re, nu_mass_eV, Omega_nu, As, ns, nrun, r, nt, inflation_consistency, Nnu, YHe)
    class(coop_cosmology_firstorder)::this
    COOP_REAL:: h, tau_re, Omega_b, Omega_c
    COOP_REAL, optional::nu_mass_eV, Omega_nu, As, ns, nrun, r, nt, Nnu, YHe
    logical::scalar_de
    COOP_INT::err
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
    call this%add_species(coop_cdm(omega_c))
    call this%add_species(coop_de_lambda(this%Omega_k()))
    this%de_genre = COOP_PERT_NONE
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


#if DO_COUPLED_DE
  subroutine coop_cosmology_firstorder_set_coupled_DE_cosmology(this, h, omega_b, omega_c, tau_re, nu_mass_eV, Omega_nu, As, ns, nrun, r, nt, inflation_consistency, Nnu, YHe, fwp1, fQ)
    class(coop_cosmology_firstorder)::this
    COOP_REAL:: h, tau_re, Omega_b, Omega_c
    COOP_REAL, optional::nu_mass_eV, Omega_nu, As, ns, nrun, r, nt, Nnu, YHe
    type(coop_function)::fwp1, fQ  !!1+w(a)  and Q(a)
    COOP_INT::err  !!
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
    call coop_background_add_coupled_DE(this, Omega_c = Omega_c, fQ = fQ, fwp1 =fwp1, err = err)
    this%de_genre = COOP_PERT_SCALAR_FIELD
    if(err .ne. 0)then
       if(coop_feedback_level .gt. 0)write(*,*) "add_coupled_DE_failed"
       call this%set_h(0.d0)
       return
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
  end subroutine coop_cosmology_firstorder_set_coupled_DE_cosmology
#endif  
  
#if DO_EFT_DE
  subroutine coop_cosmology_firstorder_set_EFT_cosmology(this, h, omega_b, omega_c, tau_re, nu_mass_eV, Omega_nu, As, ns, nrun, r, nt, inflation_consistency, Nnu, YHe, wp1, alphaM, alphaB, alphaK, alphaT, alphaH)
    class(coop_cosmology_firstorder)::this
    COOP_REAL:: h, tau_re, Omega_b, Omega_c
    COOP_REAL, optional::nu_mass_eV, Omega_nu, As, ns, nrun, r, nt, Nnu, YHe
    type(coop_function),optional::alphaM, alphaB, alphaK, alphaT, alphaH, wp1
    COOP_INT::err
    logical,optional::inflation_consistency
    type(coop_species)::de
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
    call this%add_species(coop_cdm(omega_c))
    call coop_background_add_EFT_DE(this, wp1, alphaM, err)
    if(err .ne. 0)then
       if(coop_feedback_level.gt.0)write(*,*) "Warning: EFT DE settings failed"
       call this%set_h(0.d0)
       return
    endif
    this%de_genre = COOP_PERT_EFT
    if(present(alphaB))then
       this%f_alpha_B =  alphaB
    endif
    if(present(alphaK))then
       this%f_alpha_K = alphaK
    endif
    if(present(alphaT))then
       this%f_alpha_T = alphaT
    endif
    if(present(alphaH))then
       this%f_alpha_H = alphaH
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
  end subroutine coop_cosmology_firstorder_set_EFT_cosmology
#endif

#if DO_ZETA_TRANS
  subroutine coop_cosmology_firstorder_set_zeta_weight_latevis(this)
    class(coop_cosmology_firstorder)::this    
    COOP_INT, parameter::n = 1024
    COOP_REAL::chi(n), vis(n), a, chiend
    COOP_INT::i
    call coop_set_uniform(n, chi, 0.d0, this%tau0)
    chiend = this%tau0 - this%tauofa(1.d0/(1.d0+this%zre+this%deltaz*5.d0))
    do i=1, n
       if(chi(i) .lt. chiend)then
          a = this%aoftau(this%tau0 - chi(i))
          vis(i) = this%visofa(a)
       else
          vis(i) = 0.d0
       endif
    enddo
    vis = vis/(sum(vis)*(chi(2)-chi(1)))
    call coop_zeta_user_specified_weight%init(n = n, xmin = chi(1), xmax = chi(n), f = vis, method = COOP_INTERPOLATE_LINEAR, name="latevis")
    coop_zeta_single_slice_chi = -1.d0    
  end subroutine coop_cosmology_firstorder_set_zeta_weight_latevis


  subroutine coop_cosmology_firstorder_set_zeta_weight_earlyvis(this)
    class(coop_cosmology_firstorder)::this
    COOP_INT, parameter::n = 1024
    COOP_REAL::chi(n), vis(n), a, chiend
    COOP_INT::i
    call coop_set_uniform(n, chi, 0.d0, this%tau0)
    chiend = this%tau0 - this%tauofa(1.d0/(1.d0+this%zre+this%deltaz*5.d0))
    do i=1, n
       if(chi(i) .ge. chiend)then
          a = this%aoftau(this%tau0 - chi(i))
          vis(i) = this%visofa(a)
       else
          vis(i) = 0.d0
       endif
    enddo
    vis = vis/(sum(vis)*(chi(2)-chi(1)))
    call coop_zeta_user_specified_weight%init(n = n, xmin = chi(1), xmax = chi(n), f = vis, method = COOP_INTERPOLATE_LINEAR, name="earlyvis")
    coop_zeta_single_slice_chi = -1.d0    
  end subroutine coop_cosmology_firstorder_set_zeta_weight_earlyvis

  subroutine coop_cosmology_firstorder_set_zeta_weight_vis(this)
    class(coop_cosmology_firstorder)::this    
    call coop_zeta_user_specified_weight%free()
    coop_zeta_single_slice_chi = -1.d0
  end subroutine coop_cosmology_firstorder_set_zeta_weight_vis

  subroutine coop_cosmology_firstorder_set_zeta_weight_single_slice(this, chi)
    class(coop_cosmology_firstorder)::this
    COOP_REAL,optional::chi
    call coop_zeta_user_specified_weight%free()
    if(present(chi))then
       coop_zeta_single_slice_chi = chi
    else
       coop_zeta_single_slice_chi = this%distlss       
    endif
  end subroutine coop_cosmology_firstorder_set_zeta_weight_single_slice

#endif  
  
  
end module coop_firstorder_mod

