module coop_firstorder_mod
  use coop_cl_indices_mod
  use coop_lensing_mod
  use coop_wrapper_background
  use coop_pertobj_mod
  implicit none
#include "constants.h"
  private

  public::coop_cosmology_firstorder, coop_cosmology_firstorder_source,  coop_recfast_get_xe,  coop_num_user_defined_params, coop_source_kop2k_noindex


!!recfast head file
#include "recfast_head.h"
  COOP_INT, parameter:: coop_num_user_defined_params = 10

  type coop_cosmology_firstorder_trans
     COOP_INT::lmin = 0
     COOP_INT::lmax = -1
     COOP_INT::num_l = 0
     COOP_REAL,dimension(:),allocatable::l
     COOP_REAL,dimension(:,:),allocatable::Cls, Cls2
     COOP_REAL,dimension(:,:,:,:),allocatable::trans
   contains
     procedure::free => coop_cosmology_firstorder_trans_free
     procedure::set_l => coop_cosmology_firstorder_trans_set_l
  end type coop_cosmology_firstorder_trans


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
     type(coop_cosmology_firstorder_trans)::trans
     COOP_REAL,dimension(:,:),allocatable::Cls, Cls_lensed
   contains
     procedure::free => coop_cosmology_firstorder_source_free
     procedure::get_transfer => coop_cosmology_firstorder_source_get_transfer
     procedure::get_Cls => coop_cosmology_firstorder_source_get_Cls
     procedure::set_trans => coop_cosmology_firstorder_source_set_trans
     procedure::set_Cls => coop_cosmology_firstorder_source_set_Cls
     procedure::get_Cls_limber => coop_cosmology_firstorder_source_get_Cls_limber     
     procedure::get_All_Cls => coop_cosmology_firstorder_source_get_All_Cls
     procedure::kop2k => coop_cosmology_firstorder_source_kop2k
     procedure::k2kop => coop_cosmology_firstorder_source_k2kop
     procedure::interpolate => coop_cosmology_firstorder_source_interpolate
     procedure::intbypart => coop_cosmology_firstorder_source_intbypart
     procedure::interpolate_one => coop_cosmology_firstorder_source_interpolate_one
     procedure::get_Phi_trans => coop_cosmology_firstorder_source_get_phi_trans     
     procedure::get_Psi_trans => coop_cosmology_firstorder_source_get_psi_trans
     procedure::get_PhiPlusPsi_trans => coop_cosmology_firstorder_source_get_PhiPlusPsi_trans
     procedure::get_Weyl_trans => coop_cosmology_firstorder_source_get_PhiPlusPsi_trans
     procedure::get_delta_sync_trans => coop_cosmology_firstorder_source_get_delta_sync_trans          
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
     COOP_INT::pp_genre = COOP_PP_STANDARD
     COOP_REAL::k_pivot
     
     !!these omega's are defined at a=1 (today)     
     COOP_REAL::dkappadtau_coef, ReionFrac, Omega_b, Rbya, Omega_c, Omega_nu, Omega_g, tau_eq, mnu_by_Tnu,  Omega_massivenu, bbks_keq

     !!these two parameters are defined in the a->0 limit (only matters when baryon or CDM is coupled to DE in some modified gravity models)     
     COOP_REAL::ombM2h2 = 0.022d0
     COOP_REAL::omcM2h2 = 0.12d0
     logical::inflation_consistency = .true.
     logical::w_is_background = .true.
     type(coop_function)::Ps, Pt, Xe, ekappa, vis, Tb
     type(coop_cosmology_firstorder_source),dimension(0:2)::source
     COOP_REAL, dimension(0:coop_pert_default_lmax, 0:coop_pert_default_mmax, 0:coop_pert_default_smax)::klms, klms_by_2lm1, klms_by_2lp1
     COOP_REAL,dimension(coop_pert_default_lmax)::fourbyllp1
     logical::klms_done = .false.
     logical::has_tensor = .false.
   contains
     !!standard way to initialize the cosmology (from a coop_real_table object)
     procedure::set_up =>  coop_cosmology_firstorder_set_up
     procedure::set_primordial_power => coop_cosmology_firstorder_set_primordial_power
     procedure::set_Cls => coop_cosmology_firstorder_set_Cls
     procedure::update_Cls => coop_cosmology_firstorder_update_Cls
     !!other ways to initialize the cosmology
     procedure::init_from_dictionary => coop_cosmology_firstorder_init_from_dictionary
     procedure:: set_standard_cosmology =>  coop_cosmology_firstorder_set_standard_cosmology
#if DO_EFT_DE     
     procedure:: set_EFT_cosmology =>  coop_cosmology_firstorder_set_EFT_cosmology
#endif
#if DO_COUPLED_DE
     procedure:: set_coupled_DE_cosmology =>  coop_cosmology_firstorder_set_coupled_DE_cosmology
#endif     
     procedure:: set_standard_power => coop_cosmology_firstorder_set_standard_power
     procedure:: set_user_defined_power => coop_cosmology_firstorder_set_user_defined_power
     procedure:: set_Planck_bestfit =>coop_cosmology_firstorder_set_Planck_bestfit
     procedure:: set_Planck_bestfit_with_r =>coop_cosmology_firstorder_set_Planck_bestfit_with_r
     procedure:: set_klms => coop_cosmology_firstorder_set_klms
     procedure:: set_source_tau => coop_cosmology_firstorder_set_source_tau
     procedure:: set_source_k => coop_cosmology_firstorder_set_source_k
     procedure:: set_source_ws => coop_cosmology_firstorder_set_source_ws
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
     procedure:: matter_power => coop_cosmology_firstorder_matter_power
     procedure:: smeared_matter_power => coop_cosmology_firstorder_smeared_matter_power
     procedure:: Gaussian_smeared_matter_power => coop_cosmology_firstorder_Gaussian_smeared_matter_power
     procedure:: get_Weyl_power => coop_cosmology_firstorder_get_Weyl_power
     procedure:: get_Psi_power => coop_cosmology_firstorder_get_Psi_power
     procedure:: get_Phi_power => coop_cosmology_firstorder_get_Phi_power
     procedure:: sigma_Tophat_R => coop_cosmology_firstorder_sigma_Tophat_R
     procedure:: sigma_Gaussian_R => coop_cosmology_firstorder_sigma_Gaussian_R
     procedure:: sigma_Gaussian_R_quick => coop_cosmology_firstorder_sigma_Gaussian_R_quick
     
     procedure::sigma_Gaussian_R_with_dervs => coop_cosmology_firstorder_sigma_Gaussian_R_with_dervs
     procedure::fsigma8_of_z => coop_cosmology_firstorder_fsigma8_of_z
     procedure::fgrowth_of_z => coop_cosmology_firstorder_fgrowth_of_z  !!d ln D/d lna
     procedure::sigma8_of_z => coop_cosmology_firstorder_sigma8_of_z
     procedure::growth_of_z => coop_cosmology_firstorder_growth_of_z  !! D(z)
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
          call coop_naturalspline_uniform(this%source(m)%nk, this%source(m)%saux(coop_aux_index_Weyl, :, itau), this%source(m)%saux(coop_aux_index_Weyl+1, :, itau))
          call coop_naturalspline_uniform(this%source(m)%nk, this%source(m)%saux(coop_aux_index_Psi, :, itau), this%source(m)%saux(coop_aux_index_Psi+1, :, itau))
          call coop_naturalspline_uniform(this%source(m)%nk, this%source(m)%saux(coop_aux_index_delta_sync, :, itau), this%source(m)%saux(coop_aux_index_delta_sync+1, :, itau))          
       enddo
       this%sigma_8 = this%sigma_Tophat_R(z = 0.d0, r = 8.d0/this%h()*this%H0Mpc())
    endif
  end subroutine coop_cosmology_firstorder_compute_source


  subroutine coop_cosmology_firstorder_compute_source_k(this, source, ik, output, names, output_itau, transfer_only, success)
    COOP_REAL, parameter::eps = 1.d-8
    class(coop_cosmology_firstorder)::this
    type(coop_cosmology_firstorder_source)::source
    COOP_INT ik, nvars, itau, iq, scheme
    type(coop_pert_object) pert
    COOP_REAL, dimension(:,:),allocatable::w
    logical,optional::transfer_only, success
    COOP_INT, optional::output_itau
    COOP_INT, optional::output
    type(coop_list_string), optional::names
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
    
    if(present(output))then
       if(present(output_itau))then
          if(output_itau .eq. 0) then
             if(present(names))then
                call pert%print(this, unit = output, names= names)
             else
                call pert%print(this, unit = output)
             endif
             return
          endif
       else
          if(present(names))then
             call pert%print(this, unit = output, names= names)
          else
             call pert%print(this, unit = output)
          endif
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
       if(present(output))then
          if(present(output_itau))then
             if(output_itau .eq. itau)then
                if(present(names))then
                   call pert%print(this, unit = output, names= names)
                else
                   call pert%print(this, unit = output)
                endif
                return
             endif
          else
             if(present(names))then
                call pert%print(this, unit = output, names= names)
             else
                call pert%print(this, unit = output)
             endif
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
          if(abs(pert%deMat(2, 4)).gt. 1.d-6)then
             scheme = 3
          else
             if((abs(pert%deMat(5, 4)) .gt. eps .and. pert%deMat(2, 4).eq.0.d0) .or. abs(pert%deMat(5, 4)) .gt. 1.d-4 )then
                scheme = 2
             else
                if((abs(pert%deMat(6, 4)) .gt. eps .and. pert%deMat(2, 4).eq.0.d0 .and. pert%deMat(5, 4) .eq. 0.d0 ).or. abs(pert%deMat(6, 4)) .gt. 1.d-3) then
                   scheme = 1
                else
                   scheme = 0
                endif
             endif
          endif
          if(scheme .ne. pert%de_scheme)then
             ind = 1
             pert%de_scheme = scheme
!!             if(present(output))print*, "de scheme switched to ", scheme
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
  

  subroutine coop_cosmology_firstorder_source_get_delta_sync_trans(source, tau, nk, k, delta_sync)
    class(coop_cosmology_firstorder_source) source
    COOP_INT nk
    COOP_REAL k(nk), tau, delta_sync(nk)
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
          delta_sync(ik) = (source%saux(coop_aux_index_delta_sync, 1, itau)*btau +  source%saux(coop_aux_index_delta_sync, 1, itau+1)*atau)
       elseif(ikop .ge. source%nk)then
          delta_sync(ik) = (source%saux(coop_aux_index_delta_sync, source%nk, itau)*btau +  source%saux(coop_aux_index_delta_sync, source%nk, itau+1)*atau)  * coop_bbks_trans(k(ik)/source%bbks_keq)/source%bbks_trans_kmax 
       else
          a = a - ikop
          b = 1.d0 - a
          delta_sync(ik) = ((source%saux(coop_aux_index_delta_sync, ikop, itau)+source%saux(coop_aux_index_delta_sync+1, ikop, itau)*(1.d0-b**2))*b + (source%saux(coop_aux_index_delta_sync, ikop+1, itau) + source%saux(coop_aux_index_delta_sync+1, ikop+1, itau)*(1.d0-a**2) ) * a ) * btau +  ((source%saux(coop_aux_index_delta_sync, ikop, itau+1)+source%saux(coop_aux_index_delta_sync+1, ikop, itau+1)*(1.d0-b**2))*b + (source%saux(coop_aux_index_delta_sync, ikop+1, itau+1) + source%saux(coop_aux_index_delta_sync+1, ikop+1, itau+1)*(1.d0-a**2) ) * a ) * atau
       endif
    enddo
    !$omp end parallel do
    delta_sync = delta_sync*k**2
  end subroutine coop_cosmology_firstorder_source_get_delta_sync_trans

  
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
          psi(ik) = source%saux(coop_aux_index_psi, 1, itau)*btau +  source%saux(coop_aux_index_psi, 1, itau+1)*atau
       elseif(ikop .ge. source%nk)then
          psi(ik) = (source%saux(coop_aux_index_psi, source%nk, itau)*btau +  source%saux(coop_aux_index_psi, source%nk, itau+1)*atau) * coop_bbks_trans(k(ik)/source%bbks_keq)/source%bbks_trans_kmax

       else
          a = a - ikop
          b = 1.d0 - a
          psi(ik) = ((source%saux(coop_aux_index_psi, ikop, itau)+source%saux(coop_aux_index_psi+1, ikop, itau)*(1.d0-b**2))*b + (source%saux(coop_aux_index_psi, ikop+1, itau) + source%saux(coop_aux_index_psi+1, ikop+1, itau)*(1.d0-a**2) ) * a ) * btau +  ((source%saux(coop_aux_index_psi, ikop, itau+1)+source%saux(coop_aux_index_psi+1, ikop, itau+1)*(1.d0-b**2))*b + (source%saux(coop_aux_index_psi, ikop+1, itau+1) + source%saux(coop_aux_index_psi+1, ikop+1, itau+1)*(1.d0-a**2) ) * a ) * atau
       endif
    enddo
    !$omp end parallel do
    
  end subroutine coop_cosmology_firstorder_source_get_Psi_trans

  subroutine coop_cosmology_firstorder_source_get_Phi_trans(source, tau, nk, k, phi)
    class(coop_cosmology_firstorder_source) source
    COOP_REAL::tau
    COOP_INT::nk
    COOP_REAL::k(nk), psi(nk), Phi(nk), PhiPlusPsi(nk)
    call source%get_Psi_trans(tau, nk, k, psi)
    call source%get_PhiPlusPsi_trans(tau, nk, k, PhiPlusPsi)
    Phi = PhiPlusPsi - Psi
  end subroutine coop_cosmology_firstorder_source_get_Phi_trans

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
          PhiPlusPsi(ik) = source%saux(coop_aux_index_Weyl, 1, itau)*btau +  source%saux(coop_aux_index_Weyl, 1, itau+1)*atau
       elseif(ikop .ge. source%nk)then
          PhiPlusPsi(ik) = (source%saux(coop_aux_index_Weyl, source%nk, itau)*btau +  source%saux(coop_aux_index_Weyl, source%nk, itau+1)*atau) * coop_bbks_trans(k(ik)/source%bbks_keq)/source%bbks_trans_kmax
       else
          a = a - ikop
          b = 1.d0 - a
          PhiPlusPsi(ik) = ((source%saux(coop_aux_index_Weyl, ikop, itau)+source%saux(coop_aux_index_Weyl+1, ikop, itau)*(1.d0-b**2))*b + (source%saux(coop_aux_index_Weyl, ikop+1, itau) + source%saux(coop_aux_index_Weyl+1, ikop+1, itau)*(1.d0-a**2) ) * a ) * btau +  ((source%saux(coop_aux_index_Weyl, ikop, itau+1)+source%saux(coop_aux_index_Weyl+1, ikop, itau+1)*(1.d0-b**2))*b + (source%saux(coop_aux_index_Weyl, ikop+1, itau+1) + source%saux(coop_aux_index_Weyl+1, ikop+1, itau+1)*(1.d0-a**2) ) * a ) * atau
       endif
    enddo
    !$omp end parallel do
  end subroutine coop_cosmology_firstorder_source_get_PhiPlusPsi_trans


  !!compute k^3 [ k^2/a^2 Psi_k /(3/2H^2\Omega_m) ]^2/(2pi^2); in GR it is the same as matter power
  subroutine coop_cosmology_firstorder_get_Psi_power(this, z, nk, k, Pk)
    class(coop_cosmology_firstorder)::this
    COOP_INT nk, ik
    COOP_REAL z, a, k(nk), Pk(nk), tau, Psi(nk), Ps(nk)
    a = 1.d0/(1.d0+z)
    tau = this%tauofa(a)
    call this%source(0)%get_Psi_trans(tau, nk, k, Psi)
    !$omp parallel do
    do ik = 1, nk
       ps(ik) = this%psofk(k(ik))
    enddo
    !$omp end parallel do
    pk = Psi **2 * ps * (2.d0*k**2/(O0_BARYON(this)%density(a) + O0_CDM(this)%density(a)) /a**2 * this%Mpsq(a) )**2    
  end subroutine coop_cosmology_firstorder_get_Psi_power


  !!compute k^3 [ k^2/a^2 Phi_k /(3/2H^2\Omega_m) ]^2/(2pi^2); in GR it is the same as matter power
  subroutine coop_cosmology_firstorder_get_Phi_power(this, z, nk, k, Pk)
    class(coop_cosmology_firstorder)::this
    COOP_INT nk, ik
    COOP_REAL z, a, k(nk), Pk(nk), tau, Phi(nk), Ps(nk)
    a = 1.d0/(1.d0+z)
    tau = this%tauofa(a)
    call this%source(0)%get_Phi_trans(tau, nk, k, Phi)
    !$omp parallel do
    do ik = 1, nk
       ps(ik) = this%psofk(k(ik))
    enddo
    !$omp end parallel do
    pk = Phi **2 * ps * (2.d0*k**2/(O0_BARYON(this)%density(a) + O0_CDM(this)%density(a)) /a**2 * this%Mpsq(a) )**2    
  end subroutine coop_cosmology_firstorder_get_Phi_power


!!return the dimensionless k^3 |\delta_k|^2/(2pi^2)
  subroutine coop_cosmology_firstorder_get_matter_power(this, z, nk, k, Pk)
    class(coop_cosmology_firstorder)::this
    COOP_INT nk, ik
    COOP_REAL z, a, k(nk), Pk(nk), tau, delta_sync(nk), Ps(nk)
    a = 1.d0/(1.d0+z)
    tau = this%tauofa(a)
    call this%source(0)%get_delta_sync_trans(tau, nk, k, delta_sync)
    !$omp parallel do
    do ik = 1, nk
       ps(ik) = this%psofk(k(ik))
    enddo
    !$omp end parallel do
    pk = delta_sync **2 * ps
  end subroutine coop_cosmology_firstorder_get_matter_power


  function coop_cosmology_firstorder_matter_power(this, z,  k) result(Pk)
    class(coop_cosmology_firstorder)::this
    COOP_REAL z, a, k, Pk, tau, delta_sync(1)
    a = 1.d0/(1.d0+z)
    tau = this%tauofa(a)
    call this%source(0)%get_delta_sync_trans(tau, 1, (/ k /), delta_sync)
    pk = delta_sync(1) **2 * this%psofk(k)
  end function coop_cosmology_firstorder_matter_power


  function coop_cosmology_firstorder_smeared_matter_power(this, z,  k, nw, kw, Wsq) result(Pk)
    class(coop_cosmology_firstorder)::this
    COOP_INT::nw
    COOP_REAL::pk, z, k, kw(nw), Wsq(nw)
    COOP_INT, parameter::ntheta = 72
    COOP_REAL::dcost, kp(nw),  Pkp(nw), dk3Wsq(nw),  kcost, ksq, mu
    COOP_INT::itheta, iw
    if(nw .eq. 1)then !!assume Gaussian window, sigma_W = kw(1)
       pk = this%Gaussian_smeared_matter_power(z, k, kw(1))
       return
    endif
    Pk = 0.d0
    dcost = 2.d0/ntheta
    ksq = k**2
    dk3wsq(1) = wsq(1)* ((kw(1)+kw(2))/2.d0)**3
    dk3Wsq(nw) = wsq(nw)*(3.d0*kw(nw)**2*(kw(nw)-kw(nw-1)))
    do iw = 2, nw - 1
       dk3Wsq(iw) = Wsq(iw)*( ((kw(iw+1)+kw(iw))/2.d0)**3 - ((kw(iw-1)+kw(iw))/2.d0)**3 ) 
    enddo
    do itheta = 1, ntheta
       mu = -1.d0+(itheta-0.5d0)*dcost
       kcost = 2.d0*k*mu
       do iw = 1, nw
          kp(iw) = sqrt(kw(iw)*(kw(iw)+kcost) + ksq)
       enddo
       call this%get_matter_power(z, nw, kp, Pkp)
       Pk = Pk + sum(Pkp/kp**3*dk3Wsq) !!need to divide extra k^3 because get_matter_power returns the dimensionless k^3P(k)/(2pi^2)
    enddo
    pk = pk*k**3/(sum(dk3Wsq)*ntheta)  !!still convert to the dimensionless k^3/(2pi^2) P(k)
  end function coop_cosmology_firstorder_smeared_matter_power


  !!W^2  = exp(-k^2/sigma_W^2)
  function coop_cosmology_firstorder_Gaussian_smeared_matter_power(this, z,  k, sigma_W) result(Pk)
    class(coop_cosmology_firstorder)::this
    COOP_INT,parameter::nkp = 500
    COOP_REAL::pk, z, k, sigma_W
    COOP_REAL::minkp, maxkp, dkp, kp(nkp), Pkp(nkp), w(nkp), twosig2
    COOP_INT::ikp
    maxkp = k + sigma_W*5.d0
    minkp = max(k - sigma_W*5.d0, k*1.d-4) 
    dkp = (maxkp - minkp)/(nkp-1)
    kp(1) = minkp
    do ikp = 2, nkp
       kp(ikp) = kp(ikp-1)+dkp
    enddo
    twosig2 = 2.d0*sigma_W**2
    w = kp/k * (exp(-(kp-k)**2/twosig2) -  exp(-(kp+k)**2/twosig2))
    call this%get_matter_power(z, nkp, kp, Pkp)
    pk = sum(pkp/kp**3*w)/sum(w)*k**3
  end function coop_cosmology_firstorder_Gaussian_smeared_matter_power

  
  !!compute k^3 [ k^2/a^2 (Psi_k+Phi_k)/2 /(3/2\Omega_m0/a^3) ]^2/(2pi^2); in GR it is the same as matter power
  subroutine coop_cosmology_firstorder_get_Weyl_power(this, z, nk, k, Pk)
    class(coop_cosmology_firstorder)::this
    COOP_INT nk, ik
    COOP_REAL z, a, k(nk), Pk(nk), tau, Psi(nk), Ps(nk), PhiPlusPsi(nk)
    a = 1.d0/(1.d0+z)
    tau = this%tauofa(a)
    call this%source(0)%get_PhiPlusPsi_trans(tau, nk, k, PhiPlusPsi)
    !$omp parallel do
    do ik = 1, nk
       ps(ik) = this%psofk(k(ik))
    enddo
    !$omp end parallel do
    pk = (PhiPlusPsi) **2 * ps * (k**2*a/this%Omega_m/3.d0)**2
  end subroutine coop_cosmology_firstorder_get_Weyl_power


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
  
#include "user_defined_primordial_power.h"    
#include "firstorder_basic_utils.h"
#include "firstorder_compute_sigma.h"  


  function coop_cosmology_firstorder_fsigma8_of_z(this, z) result(fsigma8)
    class(coop_cosmology_firstorder)::this
    COOP_REAL::z, fsigma8
    fsigma8 = this%fgrowth_of_z(z)*this%sigma8_of_z(z)
  end function coop_cosmology_firstorder_fsigma8_of_z

  function coop_cosmology_firstorder_fgrowth_of_z(this, z, k) result(f)
    class(coop_cosmology_firstorder)::this
    COOP_REAL::step = 0.03d0
    COOP_REAL::z, f, lna, dlna1, dlna2, lnD0, lnD1, lnD2, k_ref
    COOP_REAL,optional::k
    if(present(k))then
       k_ref= k
    else
       k_ref = this%k_pivot
    endif
    lna = -log(1.d0+z)
    if(lna + step .lt. 0.d0)then
       f = log(this%growth_of_z(exp(-lna-step)-1.d0, k_ref)/this%growth_of_z(exp(-lna+step)-1.d0, k_ref))/(2.d0*step)
       return
    endif
    if(lna + step/2.d0 .lt. 0.d0)then
       dlna1 = -step
       dlna2 = -lna
    else
       dlna1 = -step/2.d0
       dlna2 = -step
    endif
    lnD0 = log(this%growth_of_z(z, k_ref))
    lnD1 = log(this%growth_of_z(exp(-lna-dlna1)-1.d0, k_ref)) - lnD0
    lnD2 = log(this%growth_of_z(exp(-lna-dlna2)-1.d0, k_ref)) - lnD0
    f = (lnD1 * dlna2**2 - lnD2*dlna1**2)/(dlna1*dlna2*(dlna2 - dlna1))
  end function coop_cosmology_firstorder_fgrowth_of_z



  function coop_cosmology_firstorder_sigma8_of_z(this, z) result(sigma8)
    class(coop_cosmology_firstorder)::this
    COOP_REAL::z, sigma8
    sigma8 = this%sigma_tophat_R(z = z, r = 8.d0/this%h()*this%H0Mpc())
  end function coop_cosmology_firstorder_sigma8_of_z
  
!!growth function D(z) \propto delta_m (z) (!!normalized to 1/(1+z) at matter dominated regime)
  function coop_cosmology_firstorder_growth_of_z(this, z, k) result(Dz)
    class(coop_cosmology_firstorder)::this
    COOP_REAL::z, Dz
    COOP_REAL, optional::k
    COOP_REAL,parameter::a_ref = 0.03 !!matter dominate
    COOP_REAL::k_ref(1), delta(1), delta_ref(1), tau, tau_ref
    if(present(k))then
       k_ref = k
    else
       k_ref=this%k_pivot
    endif
    tau = this%tauofa(1.d0/(1.d0+z))
    tau_ref = this%tauofa(a_ref)  !!at redshift 30 where matter dominates
    call this%source(0)%get_delta_sync_trans(tau, 1, k_ref, delta)
    call this%source(0)%get_delta_sync_trans(tau_ref, 1, k_ref, delta_ref)    
    Dz = delta(1)/delta_ref(1)*a_ref
  end function coop_cosmology_firstorder_growth_of_z


  subroutine coop_cosmology_firstorder_set_standard_cosmology(this, h, omega_b, omega_c, tau_re, nu_mass_eV, Omega_nu, As, ns, nrun, r, nt, inflation_consistency, Nnu, YHe, tcmb, level)
    class(coop_cosmology_firstorder)::this
    COOP_REAL:: h, tau_re, Omega_b, Omega_c
    COOP_REAL, optional::nu_mass_eV, Omega_nu, As, ns, nrun, r, nt, Nnu, YHe, tcmb
    logical::scalar_de
    COOP_INT::err
    logical,optional::inflation_consistency
    COOP_INT,optional::level
    COOP_INT::init_level
    if(present(level))then
       init_level = level
    else
       init_level = coop_init_level_set_pp
    endif
    if(present(tcmb))then
       if(present(YHe))then
          if(present(Nnu))then
             call this%init(h=h, YHe = YHe, Nnu = NNu, Tcmb = tcmb)
          else
             call this%init(h=h, YHe = YHe, Tcmb = Tcmb)
          endif
       else
          if(present(Nnu))then
             call this%init(h=h, Nnu = NNu, Tcmb = Tcmb)
          else
             call this%init(h=h, Tcmb = tcmb)
          endif
       endif
    else
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
    if(init_level .le. coop_init_level_set_species)return
    call this%setup_background()
    if(init_level .le. coop_init_level_set_background)return
    this%optre = tau_re
    call this%set_xe()
    if(init_level .le. coop_init_level_set_xe) return    
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
  subroutine coop_cosmology_firstorder_set_coupled_DE_cosmology(this, tcmb, h, omega_b, omega_c, tau_re, nu_mass_eV, Omega_nu, As, ns, nrun, r, nt, inflation_consistency, Nnu, YHe, fwp1, fQ, level)
    class(coop_cosmology_firstorder)::this
    COOP_REAL:: h, tau_re, Omega_b, Omega_c
    COOP_REAL, optional::nu_mass_eV, Omega_nu, As, ns, nrun, r, nt, Nnu, YHe, tcmb
    type(coop_function)::fwp1, fQ  !!1+w(a)  and Q(a)
    COOP_INT::err  !!
    logical,optional::inflation_consistency
    COOP_INT,optional::level
    COOP_INT::init_level
    if(present(level))then
       init_level = level
    else
       init_level = coop_init_level_set_pp
    endif
    if(present(tcmb))then
       if(present(YHe))then
          if(present(Nnu))then
             call this%init(h=h, YHe = YHe, Nnu = NNu, Tcmb = tcmb)
          else
             call this%init(h=h, YHe = YHe, Tcmb = Tcmb)
          endif
       else
          if(present(Nnu))then
             call this%init(h=h, Nnu = NNu, Tcmb = Tcmb)
          else
             call this%init(h=h, Tcmb = tcmb)
          endif
       endif
    else
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
    if(init_level .le. coop_init_level_set_species)return
    call this%setup_background()
    if(init_level .le. coop_init_level_set_background) return
    this%optre = tau_re
    call this%set_xe()
    if(init_level .le. coop_init_level_set_xe)return
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
  subroutine coop_cosmology_firstorder_set_EFT_cosmology(this, h, tcmb, omega_b, omega_c, tau_re, nu_mass_eV, Omega_nu, As, ns, nrun, r, nt, inflation_consistency, Nnu, YHe, wp1, wp1_background, alphaM, alphaB, alphaK, alphaT, alphaH, level)
    class(coop_cosmology_firstorder)::this
    COOP_REAL:: h, tau_re, Omega_b, Omega_c
    COOP_REAL, optional::nu_mass_eV, Omega_nu, As, ns, nrun, r, nt, Nnu, YHe, tcmb
    type(coop_function),optional::alphaM, alphaB, alphaK, alphaT, alphaH, wp1, wp1_background
    COOP_INT::err
    logical,optional::inflation_consistency
    type(coop_species)::de
    COOP_INT,optional::level
    COOP_INT::init_level
    if(present(level))then
       init_level = level
    else
       init_level = coop_init_level_set_pp
    endif
    if(present(tcmb))then
       if(present(YHe))then
          if(present(Nnu))then
             call this%init(h=h, YHe = YHe, Nnu = NNu, Tcmb = tcmb)
          else
             call this%init(h=h, YHe = YHe, Tcmb = Tcmb)
          endif
       else
          if(present(Nnu))then
             call this%init(h=h, Nnu = NNu, Tcmb = Tcmb)
          else
             call this%init(h=h, Tcmb = tcmb)
          endif
       endif
    else
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
    endif
    if(present(alphaM))then
       call this%set_alphaM(alphaM)
    else
       call this%set_alphaM(coop_function_polynomial( (/ 0.d0 /) ) )
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

    if(present(wp1))then
       call coop_background_add_EFT_DE(this, wp1, err)
       if(present(wp1_background))then
          stop "you should not pass both wp1 and wp1_background as the arguments"
       endif
    elseif(present(wp1_background))then
       call coop_background_add_EFT_DE_with_effective_w(this, wp1_background, err)
    else
       stop  "you should pass either wp1 or wp1_background to EFT DE"
    endif

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
    if(init_level.le.coop_init_level_set_species)return
    call this%setup_background()
    if(init_level.le.coop_init_level_set_background)return
    this%optre = tau_re
    call this%set_xe()
    if(init_level .le. coop_init_level_set_xe)return
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

  subroutine coop_cosmology_firstorder_init_from_dictionary(this, paramtable)
    class(coop_cosmology_firstorder)::this
    type(coop_dictionary)::paramtable
    logical success
    type(coop_real_table)::params
    call coop_dictionary_lookup(paramtable, "w_is_background", this%w_is_background)
    call coop_dictionary_lookup(paramtable, "inflation_consistency", this%inflation_consistency, .true.)
    call coop_dictionary_lookup(paramtable, "pp_genre", this%pp_genre, COOP_PP_STANDARD)
    call params%load_dictionary(paramtable)
    call this%set_up(params, success)
  end subroutine coop_cosmology_firstorder_init_from_dictionary




  subroutine coop_cosmology_firstorder_set_up(this, paramtable, success, level)
    class(coop_cosmology_firstorder)::this
    type(coop_real_table)::paramtable
    logical::success
    COOP_INT,optional::level
    COOP_INT::init_level
  !!----------------------------------------
    COOP_REAL:: tau_re, h, H0, theta, omega_b, omega_c, w0, wa, tcmb

#if DO_EFT_DE  
    COOP_REAL::alpha_M0 = 0.d0
    COOP_REAL::alpha_T0 = 0.d0
    COOP_REAL::alpha_B0 = 0.d0
    COOP_REAL::alpha_K0 = 0.d0
    COOP_REAL::alpha_H0 = 0.d0
    COOP_REAL::de_cs2 = -1.d0
    COOP_REAL::r_B = 0.d0
    COOP_REAL::r_T = 0.d0
    COOP_REAL::r_M = 0.d0
    COOP_REAL::r_H = 0.d0    
    
#elif DO_COUPLED_DE
    COOP_REAL::Q0 = 0.d0
    COOP_REAL::Qa = 0.d0  
#endif
    COOP_REAL::h_low, h_high, h_mid, theta_high, theta_low, theta_mid
    COOP_REAL, parameter::min_omega_m = 0.001d0, max_omega_m = 0.999d0, h_min = 0.4d0, h_max = 1.d0, min_omega_m_conservative = 0.1d0, max_omega_m_conservative = 0.5d0
    type(coop_function)::fwp1, fQ, alphaM, alphaB, alphaK, alphaT, alphaH
    logical::w_predefined, alpha_predefined

    COOP_REAL::eps_inf, eps_s, zeta_s, beta_s

    !!0 background
    !!1 x_e
    !!2 primoridal power
    !!3 linear perturbations
    !!4 C_l
    !!5 tensor C_l
    if(present(level))then
       init_level = level
    else
       init_level = coop_init_level_set_pp
    endif

    call paramtable%lookup( "ombh2", this%ombm2h2, 0.d0)
    if(this%ombm2h2.eq.0.d0) call paramtable%lookup( "omegabh2", this%ombm2h2)
    call paramtable%lookup( "omch2", this%omcm2h2, 0.d0)
    if(this%omcm2h2.eq.0.d0) call paramtable%lookup( "omegach2", this%omcm2h2)       
    if(this%ombm2h2 .le. 0.d0 .or. this%omcm2h2.le.0.d0)then
       write(*,*) "cosmology_firstorder_init: you must have input parameters omegabh2 (>0) and omegach2 (>0) "
       stop
    endif

    
    call paramtable%lookup( "H0", H0, 0.d0)
    h = H0/100.d0
    if(h.eq.0.d0) call paramtable%lookup( "h", h, 0.d0)
    call paramtable%lookup("tcmb", tcmb, COOP_DEFAULT_TCMB)

    call paramtable%lookup( "tau", tau_re)    

    call paramtable%lookup( "de_w", w0, -1000.d0)
    if(w0 .gt. -999.d0)then
       call paramtable%lookup( "de_wa",wa, 0.d0)  
       call fwp1%init_polynomial( (/ 1.d0+w0+wa, -wa /) )
       w_predefined = .true.
    else
       call paramtable%lookup( "de_epss", eps_s, -1000.d0)
       if(eps_s .gt. -999.d0)then
          call paramtable%lookup( "de_epsinf",eps_inf, 0.d0)
          call paramtable%lookup( "de_zetas",zeta_s, 0.d0)
          call paramtable%lookup( "de_betas",beta_s, 6.d0)
          w_predefined = .false.
       else
          w0 = -1.d0
          wa = 0.d0
          call fwp1%init_polynomial( (/ 0.d0 /) )
          w_predefined = .true.          
       endif
    endif

    
#if DO_EFT_DE
    call paramtable%lookup( "de_cs2", de_cs2, -1.d30)
    if(de_cs2 .ge. -100.d0)then
       if(de_cs2 .le. 0.d0) stop "Error: dark energy c_s^2 <= 0"
       call paramtable%lookup( "de_r_B",  r_B, 0.d0)
       call paramtable%lookup( "de_r_H",  r_H, 0.d0)
       call paramtable%lookup( "de_r_M",  r_M, 0.d0)
       call paramtable%lookup( "de_r_T",  r_T, 0.d0)
       call paramtable%lookup( "de_alpha0_max", coop_de_alpha0_max, 1.d0)
       if(.not. w_predefined .or. wa .ne. 0.d0) stop "Error: Incompatible settings for de_w and de_cs2. For EFT DE with fixed c_s^2, w is assumed to be a constant, too."
       alpha_predefined = .false.
    else
       call paramtable%lookup( "de_alpha_M0", alpha_M0, 0.d0 )
       call paramtable%lookup( "de_alpha_T0", alpha_T0, 0.d0 )
       call paramtable%lookup( "de_alpha_B0", alpha_B0, 0.d0 )
       call paramtable%lookup( "de_alpha_K0", alpha_K0, 0.d0 )
       call paramtable%lookup( "de_alpha_H0", alpha_H0 , 0.d0 )
       alphaM = coop_de_alpha_constructor(alpha_M0, "omega")
       alphaT = coop_de_alpha_constructor(alpha_T0, "omega")
       alphaH = coop_de_alpha_constructor(alpha_H0, "omega")
       alphaB = coop_de_alpha_constructor(alpha_B0, "omega")
       alphaK = coop_de_alpha_constructor(alpha_K0, "omega")
       alpha_predefined = .true.
       call this%set_alpham(alphaM)       
    endif
#elif DO_COUPLED_DE
    call paramtable%lookup( "de_Q", Q0, 0.d0)
    call paramtable%lookup( "de_Qa", Qa, 0.d0)  
#endif
    if(h .ne. 0.d0) then
       call setforH(h, success, init_level)
       goto 120
    endif
    call paramtable%lookup( "theta", theta, 0.d0)
    theta = theta/100.d0
    if(theta .eq. 0.d0)then
       stop " you must either set theta or H0 in the parameter file"
    endif
    if(w_predefined .and. alpha_predefined)then
       h_low = max(sqrt((this%ombm2h2 + this%omcm2h2)/max_omega_m), h_min)
       h_high = min(sqrt((this%ombm2h2 + this%omcm2h2)/min_omega_m), h_max)
    else
       !!for dynamic w models we require a more conservative range 0.1<Omega_m < 0.5
       h_low = max(sqrt((this%ombm2h2 + this%omcm2h2)/max_omega_m_conservative), h_min) 
       h_high = min(sqrt((this%ombm2h2 + this%omcm2h2)/min_omega_m_conservative), h_max)
    endif
    if(h_low .ge. h_high) then
       success  = .false.
       return
    endif

    call setforH(h_low, success, level = 0)
    if(.not. success) return
    theta_low = this%cosmomc_theta()
    call setforH(h_high, success, level = 0)
    if(.not. success) return
    theta_high = this%cosmomc_theta()

    if(theta_low .gt. theta .or. theta_high .lt. theta)then
       success = .false.
       return
    endif
    do while(h_high - h_low .gt. 1.d-4 .and. theta_high - theta_low .gt. 1.d-5)
       h_mid = (h_high + h_low)/2.d0
       call setforH(h_mid, success, level = 0)
       if(.not. success) return
       theta_mid = this%cosmomc_theta()
       if(theta_mid  .gt. theta)then
          h_high = h_mid
          theta_high = theta_mid
       else
          h_low = h_mid
          theta_low = theta_mid
       endif
    enddo
    h_mid = (h_low*(theta_high - theta) + h_high*(theta - theta_low))/(theta_high - theta_low)
100 call setforH(h_mid, success, level = init_level )
120 if(.not. success .or. init_level .le. coop_init_level_set_xe) return
    call this%set_primordial_power(paramtable)
    if(init_level .le. coop_init_level_set_pp)return
    call this%compute_source(0, success)
    if(.not. success .or. init_level .le. coop_init_level_set_pert) return
    call this%set_cls(0, 2, coop_cls_lmax(0))
    if(init_level .le. coop_init_level_set_Cls) return
    call this%compute_source(2, success)
    call this%set_cls(2, 2, coop_cls_lmax(2))
    if(.not. success) return
  contains
    
    
    subroutine setforH(hubble, success, level)
      logical:: success
      COOP_INT::iloop
      COOP_REAL::hubble, omlast
      COOP_INT::level

#if DO_EFT_DE      
      if(alpha_predefined)then
#endif         
         Omega_b = this%ombm2h2/this%Mpsq0/hubble**2
         Omega_c = this%omcm2h2/this%Mpsq0/hubble**2

         if(.not. w_predefined)then
            fwp1 = coop_function_constructor(coop_de_wp1_coupled_quintessence, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., args = coop_arguments_constructor( r = (/ 1.d0-omega_b - omega_c, eps_s, eps_inf, zeta_s , beta_s /) ), name = "DE 1+w")
         endif
#if DO_EFT_DE               
      else
         Omega_b = this%ombm2h2/hubble**2
         Omega_c = this%omcm2h2/hubble**2         
         omlast = -1000.d0         
         iloop = 0
         do while(abs(omega_b + omega_c - omlast ) .gt. 3.d-5)
            omlast = omega_b + omega_c            
            call coop_de_construct_alpha_from_cs2(omlast, w0, de_cs2, r_B, r_H, r_M, r_T, alphaB, alphaH, alphaK, alphaM, alphaT, success)
            if(.not. success)return
            call this%set_alpham(alphaM)
            Omega_b = this%ombm2h2/this%Mpsq0/hubble**2
            Omega_c = this%omcm2h2/this%Mpsq0/hubble**2
            iloop = iloop + 1
            if(iloop .gt. 100)then
               success = .false.
               return
            endif
         enddo
      endif
#endif      
#if DO_EFT_DE  
      !!initialize this
      if(this%w_is_background .or. .not. alpha_predefined)then
         call this%set_EFT_cosmology(Omega_b=omega_b, Omega_c=omega_c, h = hubble, Tcmb = tcmb, tau_re = tau_re,  wp1_background = fwp1, alphaM = alphaM, alphaK = alphaK, alphaB= alphaB, alphaH = alphaH, alphaT = alphaT, level = level)     
      else
         call this%set_EFT_cosmology(Omega_b=Omega_b, Omega_c=Omega_c, h = hubble, Tcmb =tcmb, tau_re = tau_re, wp1 = fwp1, alphaM = alphaM, alphaK = alphaK, alphaB= alphaB, alphaH = alphaH, alphaT = alphaT, level = level )
      endif
#elif DO_COUPLED_DE
      call fQ%init_polynomial( (/ Q0+Qa, -Qa /) )
      call this%set_coupled_DE_cosmology(Omega_b=Omega_b, Omega_c=Omega_c, h = hubble, Tcmb =tcmb, tau_re = tau_re, fwp1 = fwp1, fQ = fQ, level = level )
#else
      call this%set_standard_cosmology(Omega_b=Omega_b, Omega_c=Omega_c, h = hubble, Tcmb = tcmb, tau_re = tau_re, level = level)
#endif
      if(this%h_value .eq. 0.d0)then
         success = .false.
      else 
         success = .true.
      endif
    end subroutine setforH
  end subroutine coop_cosmology_firstorder_set_up
  
  subroutine coop_cosmology_firstorder_set_primordial_power(this, paramtable)
    class(coop_cosmology_firstorder)::this
    type(coop_real_table)::paramtable
    COOP_REAL::tau_re
    COOP_REAL::upar(coop_num_user_defined_params)
    COOP_INT::i
    select case(this%pp_genre)
    case(COOP_PP_STANDARD)
       call paramtable%lookup( "tau", tau_re)    
       call paramtable%lookup( "As", this%As, 0.d0)
       if(this%As .eq. 0.d0)then
          call paramtable%lookup( "logA", this%As, 0.d0)
          if(this%As .eq. 0.d0)then
             call paramtable%lookup( "logAm2tau", this%As, 0.d0)
             if(this%As .eq. 0.d0)then
                stop "You need to set As, logA,  or logAm2tau in the parameter file"
             endif
             this%As = exp(this%As + 2.d0*tau_re)*1.d-10
          else
             this%As = exp(this%As)*1.d-10
          endif
       endif
       call paramtable%lookup( "ns", this%ns)
       call paramtable%lookup( "r", this%r, 0.d0)
       call paramtable%lookup( "nt", this%nt, 0.d0)
       call paramtable%lookup( "nrun", this%nrun, 0.d0)
       call this%set_standard_power(this%As, this%ns, this%nrun, this%r, this%nt, this%inflation_consistency)
    case(100:110) 
       do i=1, coop_num_user_defined_params
          call paramtable%lookup("user_pp"//COOP_STR_OF(i), upar(i), 0.d0)
       enddo
       call this%set_user_defined_power(upar)
    case default
       write(*,*) "pp_genre = "//COOP_STR_OF(this%pp_genre)
       Write(*,*) "Unknown pp genre. For standard power spectrum use pp_genre = 0; for user defined power spectrum use pp_genre = 101, 102, ..., 110. Modify include/user_defined_primordial_power.h to define your model."
       stop
    end select
  end subroutine coop_cosmology_firstorder_set_primordial_power


  subroutine coop_cosmology_firstorder_trans_free(this)
    class(coop_cosmology_firstorder_trans)::this
    COOP_DEALLOC(this%l)
    COOP_DEALLOC(this%Cls)
    COOP_DEALLOC(this%Cls2)
    COOP_DEALLOC(this%trans)
    this%num_l = 0
  end subroutine coop_cosmology_firstorder_trans_free

  subroutine coop_cosmology_firstorder_trans_set_l(this, lmin, lmax)
    class(coop_cosmology_firstorder_trans)::this
    COOP_INT::lmin, lmax, l, i
    if(allocated(this%l) .and. this%num_l.gt.0)then
       if(this%l(1) .eq. lmin .and. this%l(this%num_l).eq. lmax)return
    endif
    COOP_DEALLOC(this%l)
    this%lmin = lmin
    this%lmax = lmax
    l = lmin
    this%num_l = 1
    do while(l.lt.lmax)
      call coop_next_l(l)
      this%num_l = this%num_l + 1
    enddo
    allocate(this%l(this%num_l))
    l = lmin
    i = 1
    this%l(1) = dble(lmin)
    do 
       call coop_next_l(l)
       i = i + 1
       if(l .ge. lmax)then
          this%l(i) = dble(lmax)
          exit
       else
          this%l(i) = dble(l)
       endif
    enddo
  end subroutine coop_cosmology_firstorder_trans_set_l


  subroutine coop_cosmology_firstorder_source_set_trans(source, lmin, lmax)
    class(coop_cosmology_firstorder_source)::source
    COOP_INT::lmin, lmax, i
    call source%trans%set_l(lmin, lmax)
    if(allocated(source%trans%trans))then
       if(size(source%trans%trans, 1) .ne. source%nsrc .or. size(source%trans%trans, 2).ne. coop_k_dense_fac .or. size(source%trans%trans,3).ne. source%nk .or. size(source%trans%trans, 4).ne. source%trans%num_l)then
          deallocate(source%trans%trans)
          COOP_DEALLOC(source%trans%Cls)
          COOP_DEALLOC(source%trans%Cls2)
          allocate(source%trans%trans(source%nsrc, coop_k_dense_fac, source%nk, source%trans%num_l), source%trans%Cls(source%trans%num_l, coop_num_Cls), source%trans%Cls2(source%trans%num_l, coop_num_Cls))
       endif
    else
       COOP_DEALLOC(source%trans%Cls)
       COOP_DEALLOC(source%trans%Cls2)
       allocate(source%trans%trans(source%nsrc, coop_k_dense_fac, source%nk, source%trans%num_l), source%trans%Cls(source%trans%num_l, coop_num_Cls), source%trans%Cls2(source%trans%num_l, coop_num_Cls))
    endif
    if(allocated(source%cls))then
       if(lbound(source%cls, 2) .ne. lmin .or. ubound(source%cls, 2) .ne. lmax)then
          deallocate(source%cls)
          COOP_DEALLOC(source%cls_lensed)
          allocate(source%Cls(coop_num_cls, lmin:lmax), source%Cls_lensed(coop_num_cls, lmin:lmax))
       endif
    else
       COOP_DEALLOC(source%cls)
       COOP_DEALLOC(source%cls_lensed)
       allocate(source%Cls(coop_num_cls, lmin:lmax), source%Cls_lensed(coop_num_cls, lmin:lmax))
    endif
    !$omp parallel do
    do i = 1, source%trans%num_l
       call source%get_transfer(nint(source%trans%l(i)), source%trans%trans(:,:,:,i))
    enddo
    !$omp end parallel do
  end subroutine coop_cosmology_firstorder_source_set_trans


  subroutine coop_cosmology_firstorder_set_source_ws(this, source)
    class(coop_cosmology_firstorder)::this
    class(coop_cosmology_firstorder_source)::source
    COOP_INT::i,j
    select case(source%m)
    case(0)
       !$omp parallel do private(i, j)
       do i=1, source%nk
          do j=1, coop_k_dense_fac
             source%ws_dense(j, i) = this%psofk(source%k_dense(j, i))*(source%dkop_dense  / (1.d0 + source%kweight * coop_source_k_index * source%k_dense(j,i)**coop_source_k_index ) )
          enddo
       enddo
       !$omp end parallel do
    case(2)
       !$omp parallel do private(i, j)
       do i=1, source%nk
          do j=1, coop_k_dense_fac
             source%wt_dense(j, i) = this%ptofk(source%k_dense(j, i))*(source%dkop_dense  / (1.d0 + source%kweight * coop_source_k_index * source%k_dense(j,i)**coop_source_k_index ) )
          enddo
       enddo
       !$omp end parallel do
    case default
       write(*,*) "Error: m  = ", source%m
       stop "COOP only support scalar and tensor perturbations."
    end select
  end subroutine coop_cosmology_firstorder_set_source_ws


  subroutine coop_cosmology_firstorder_source_set_Cls(source)
    class(coop_cosmology_firstorder_source)::source
    COOP_INT::i, j, l
    select case(source%m)
    case(0)
       do i=1, source%trans%num_l
          source%trans%Cls(i, coop_index_ClTT) = sum(source%ws_dense * source%trans%trans(coop_index_source_T, :, :, i)**2)*coop_4pi
          if(source%nsrc .ge. 2)then
             source%trans%Cls(i, coop_index_ClTE) = sqrt((source%trans%l(i)+2.d0)*(source%trans%l(i)+1.d0)*source%trans%l(i)*(source%trans%l(i)-1.d0))*sum(source%ws_dense * source%trans%trans(coop_index_source_T, :, :, i)*source%trans%trans(coop_index_source_E, :, :, i))*coop_4pi
             source%trans%Cls(i, coop_index_ClEE) = (source%trans%l(i)+2.d0)*(source%trans%l(i)+1.d0)*source%trans%l(i)*(source%trans%l(i)-1.d0)*sum(source%ws_dense * source%trans%trans(coop_index_source_E,:,:,i)**2)*coop_4pi
          endif
          source%trans%Cls(i,coop_index_ClLenLen) =  sum(source%ws_dense * source%trans%trans(coop_index_source_Len, :, :,i)**2)*coop_4pi
          source%trans%Cls(i,coop_index_ClTLen) =  sum(source%ws_dense * source%trans%trans(coop_index_source_T, :, :,i) * source%trans%trans(coop_index_source_Len,:,:,i))*coop_4pi
#if DO_ZETA_TRANS
          source%trans%Cls(i,coop_index_ClTzeta) = sum(source%ws_dense * source%trans%trans(coop_index_source_T, :, :,i)*source%trans%trans(coop_index_source_zeta,:,:,i))*coop_4pi
          source%trans%Cls(i,coop_index_ClEzeta) = sqrt((source%trans%l(i)+2.d0)*(source%trans%l(i)+1.d0)*source%trans%l(i)*(source%trans%l(i)-1.d0))*sum(source%ws_dense * source%trans%trans(coop_index_source_E, :, :,i) * source%trans%trans(coop_index_source_zeta,:,:,i))*coop_4pi
          source%trans%Cls(i, coop_index_Clzetazeta) = sum(source%ws_dense * source%trans%trans(coop_index_source_zeta, :, :, i)**2)*coop_4pi
#endif
       enddo
    case(1)
       call coop_tbw("get_source%trans%Cls: vector")
    case(2)
       do i=1, source%trans%num_l
          source%trans%Cls(i,coop_index_ClTT) =(source%trans%l(i)+2.d0)*(source%trans%l(i)+1.d0)*source%trans%l(i)*(source%trans%l(i)-1.d0)*sum(source%wt_dense * source%trans%trans(coop_index_source_T, :, :,i)**2)*coop_4pi
          if(source%nsrc .ge. 2)then

             source%trans%Cls(i, coop_index_ClTE) = sqrt((source%trans%l(i)+2.d0)*(source%trans%l(i)+1.d0)*source%trans%l(i)*(source%trans%l(i)-1.d0))*sum(source%wt_dense * source%trans%trans(coop_index_source_T, :, :,i)*source%trans%trans(coop_index_source_E,:,:,i))*coop_4pi
             source%trans%Cls(i,coop_index_ClEE) = sum(source%wt_dense * source%trans%trans(coop_index_source_E,:,:,i)**2)*coop_4pi
             if(source%nsrc .ge. 3)then
                source%trans%Cls(i,coop_index_ClBB) = sum(source%wt_dense * source%trans%trans(coop_index_source_B, :, :,i)**2)*coop_4pi
             end if
          endif
       enddo
    case default
       call coop_return_error("set_Cls", "unknown m = "//trim(coop_num2str(source%m)), "stop")
    end select
    do i = 1, source%trans%num_l
       source%trans%Cls(i,:) = source%trans%Cls(i,:) *(source%trans%l(i)*(source%trans%l(i)+1.d0))
    enddo
    do i=1, coop_num_Cls
       call coop_spline(source%trans%num_l, source%trans%l, source%trans%Cls(:, i), source%trans%Cls2(:, i))
       !$omp parallel do
       do l = source%trans%lmin, source%trans%lmax
          call coop_splint(source%trans%num_l, source%trans%l, source%trans%Cls(:, i), source%trans%Cls2(:, i), dble(l), source%Cls(i, l))
          source%Cls(i,l) = source%Cls(i,l)/(l*(l+1.d0))
       enddo
       !$omp end parallel do
    enddo
    if(source%m.eq.0.d0)then
       call coop_get_lensing_Cls(source%trans%lmin, source%trans%lmax, source%Cls, source%Cls_lensed)
       source%Cls_lensed = source%Cls_lensed + source%Cls
    else
       source%Cls_lensed = source%Cls
    endif
    
  end subroutine coop_cosmology_firstorder_source_set_Cls

  subroutine coop_cosmology_firstorder_set_Cls(this, m, lmin, lmax)
    class(coop_cosmology_firstorder)::this
    COOP_INT::m, lmin, lmax
    call this%source(m)%set_trans(lmin, lmax)
    call this%source(m)%set_Cls()
  end subroutine coop_cosmology_firstorder_set_Cls


  subroutine coop_cosmology_firstorder_update_Cls(this, m)
    class(coop_cosmology_firstorder)::this
    COOP_INT::m
    call this%set_source_ws(this%source(m))
    call this%source(m)%set_Cls()
  end subroutine coop_cosmology_firstorder_update_Cls
  
end module coop_firstorder_mod

