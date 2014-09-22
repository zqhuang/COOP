module coop_firstorder_mod
  use coop_wrapper_background
  use coop_pertobj_mod
  implicit none
#include "constants.h"

private

  public::coop_cosmology_firstorder, coop_cosmology_firstorder_source,  coop_recfast_get_xe, coop_power_lnk_min, coop_power_lnk_max,  coop_k_dense_fac, coop_index_ClTT, coop_index_ClTE, coop_index_ClEE, coop_index_ClBB, coop_index_ClEB, coop_index_ClTB, coop_index_ClLenLen, coop_index_ClTLen, coop_num_Cls, coop_scalar_lmax, coop_vector_lmax, coop_tensor_lmax

  COOP_INT::coop_scalar_lmax = 2500
  COOP_INT::coop_vector_lmax = 1500
  COOP_INT::coop_tensor_lmax = 1000
  COOP_REAL, parameter :: coop_power_lnk_min = log(0.1d0) 
  COOP_REAL, parameter :: coop_power_lnk_max = log(5.d3) 
  COOP_REAL, parameter :: coop_visibility_amin = 1.8d-4
  COOP_REAL, parameter :: coop_initial_condition_epsilon = 1.d-6
  COOP_REAL, parameter :: coop_cosmology_firstorder_ode_accuracy = 1.d-8
  COOP_REAL, parameter :: coop_cosmology_firstorder_tc_cutoff = 0.005d0


  COOP_REAL, dimension(0:2), parameter::coop_source_tau_weight = (/ 0.3d0, 0.4d0, 0.3d0 /)
  COOP_INT, dimension(0:2), parameter::coop_source_tau_n = (/ 1200, 800, 750 /)
  COOP_REAL, dimension(0:2), parameter::coop_source_k_weight = (/ 0.15d0, 0.15d0, 0.1d0 /)
  COOP_INT, dimension(0:2), parameter::coop_source_k_n = (/ 150, 120, 80 /)
  COOP_REAL, parameter::coop_source_k_index = 0.4d0
  COOP_INT, parameter:: coop_k_dense_fac = 30


  COOP_INT, parameter::coop_index_ClTT = 1
  COOP_INT, parameter::coop_index_ClEE = 2
  COOP_INT, parameter::coop_index_ClBB = 3
  COOP_INT, parameter::coop_index_ClTE = 4
  COOP_INT, parameter::coop_index_ClEB = 5
  COOP_INT, parameter::coop_index_ClTB = 6
  COOP_INT, parameter::coop_index_ClLenLen = 7
  COOP_INT, parameter::coop_index_ClTLen = 8
  COOP_INT, parameter::coop_num_Cls =  coop_index_ClTLen

!!how many source terms you want to extract & save
  COOP_INT, dimension(0:2), parameter::coop_num_sources = (/ 2,  3,  3 /)

!!recfast head file
#include "recfast_head.h"


  type coop_cosmology_firstorder_source
     COOP_INT::m = 0
     COOP_INT::ntau = 0
     COOP_INT::nk = 0
     COOP_INT::nsrc = 0
     COOP_INT::index_tc_max = 1
     COOP_INT, dimension(:),allocatable::index_tc_off
     COOP_INT, dimension(coop_pert_default_nq)::index_massivenu_on
     COOP_REAL::dkop, kopmin, kopmax, kmin, kmax, kweight, tauweight
     COOP_REAL,dimension(coop_k_dense_fac)::a_dense, b_dense, a2_dense, b2_dense
     COOP_REAL, dimension(:),allocatable::k, kop, tau, a, tauc, lna,  dk, dtau, chi !!tau is conformal time, chi is comoving distance; in flat case, chi + tau = tau_0
     COOP_REAL, dimension(:,:),allocatable::k_dense, ws_dense, wt_dense, dk_dense
     COOP_REAL, dimension(:,:,:),allocatable::s, s2
   contains
     procedure::free => coop_cosmology_firstorder_source_free
     procedure::get_transfer => coop_cosmology_firstorder_source_get_transfer
     procedure::get_Cls => coop_cosmology_firstorder_source_get_Cls
     procedure::get_All_Cls => coop_cosmology_firstorder_source_get_All_Cls
     procedure::kop2k => coop_cosmology_firstorder_source_kop2k
     procedure::k2kop => coop_cosmology_firstorder_source_k2kop
  end type coop_cosmology_firstorder_source
  
  type, extends(coop_cosmology_background) :: coop_cosmology_firstorder
     logical::do_reionization = .true.
     COOP_REAL::zrecomb, distlss, tau0
     COOP_REAL::optre = 0.07d0
     COOP_REAL::zre = 8.d0
     COOP_REAL::deltaz = 1.5d0
     COOP_REAL::kMpc_pivot = 0.05d0
     COOP_INT ::de_genre = COOP_PERT_NONE
     COOP_REAL::k_pivot 
     COOP_REAL::dkappadtau_coef, ReionFrac, Omega_b, Omega_c, Omega_nu, Omega_g, tau_eq, mnu_by_Tnu, As, ns, nrun, r, nt
     logical::inflation_consistency
     type(coop_function)::Ps, Pt, Xe, ekappa, vis, Tb
     type(coop_cosmology_firstorder_source),dimension(0:2)::source
     COOP_INT::index_baryon, index_cdm, index_radiation, index_nu, index_massiveNu, index_de
     COOP_REAL, dimension(0:coop_pert_default_lmax, 0:coop_pert_default_mmax, 0:coop_pert_default_smax)::klms, klms_by_2lm1, klms_by_2lp1
     logical::klms_done = .false.
   contains
     procedure:: set_standard_cosmology =>  coop_cosmology_firstorder_set_standard_cosmology
     procedure:: set_standard_power => coop_cosmology_firstorder_set_standard_power
     procedure:: set_Planck_bestfit =>coop_cosmology_firstorder_set_Planck_bestfit
     procedure:: set_klms => coop_cosmology_firstorder_set_klms
     procedure:: set_source_tau => coop_cosmology_firstorder_set_source_tau
     procedure:: set_source_k => coop_cosmology_firstorder_set_source_k
     procedure:: set_power => coop_cosmology_firstorder_set_power
     procedure:: set_xe => coop_cosmology_firstorder_set_xe
     procedure:: set_zre_from_optre => coop_cosmology_firstorder_set_zre_from_optre
     procedure:: set_optre_from_zre => coop_cosmology_firstorder_set_optre_from_zre
     procedure::set_initial_conditions => coop_cosmology_firstorder_set_initial_conditions
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
     procedure::pert2source => coop_cosmology_firstorder_pert2source
     procedure::init_source => coop_cosmology_firstorder_init_source
     procedure::compute_source =>  coop_cosmology_firstorder_compute_source
     procedure::compute_source_k =>  coop_cosmology_firstorder_compute_source_k
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
       call this%set_standard_cosmology(Omega_b=0.04801228639964265020d0, Omega_c=0.26099082900159878590d0, h = 0.67779d0, tau_re = 0.092d0, Omega_nu = Omega_nu, As = 2.210168d-9, ns = 0.9608933d0)
    else
       call this%set_standard_cosmology(Omega_b=0.04801228639964265020d0, Omega_c=0.26099082900159878590d0, h = 0.67779d0, tau_re = 0.092d0, As = 2.210168d-9, ns = 0.9608933d0)
    endif
  end subroutine coop_cosmology_firstorder_set_Planck_Bestfit



  subroutine coop_cosmology_firstorder_set_standard_cosmology(this, h, omega_b, omega_c, tau_re, nu_mass_eV, Omega_nu, As, ns, nrun, r, nt, inflation_consistency, Nnu, YHe)
    class(coop_cosmology_firstorder)::this
    COOP_REAL:: h, Omega_b, Omega_c,  tau_re
    COOP_REAL, optional::nu_mass_eV, Omega_nu, As, ns, nrun, r, nt, Nnu, YHe
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
    call this%add_species(coop_cdm(COOP_REAL_OF(Omega_c)))
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
    call this%add_species(coop_de_lambda(this%Omega_k()))
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


  subroutine coop_cosmology_firstorder_set_standard_power(this, As, ns, nrun, r, nt, inflation_consistency)
    class(coop_cosmology_firstorder)::this
    COOP_REAL kpivot, As, ns, nrun, r, nt
    type(coop_arguments):: args
    logical,optional::inflation_consistency
    if(present(inflation_consistency))then
       if(inflation_consistency)then
          call args%init( r = (/ As, ns, nrun, r, -r/8.d0 /) )
       else
          call args%init( r = (/ As, ns, nrun, r, nt /) )
       endif
    else
       call args%init( r = (/ As, ns, nrun, r, nt /) )
    endif
    call this%set_power(coop_cosmology_firstorder_standard_power, args)
    call args%free()
  end subroutine coop_cosmology_firstorder_set_standard_power

  subroutine coop_cosmology_firstorder_standard_power(kbykpiv, ps, pt, cosmology, args)
    type(coop_cosmology_firstorder)::cosmology
    COOP_REAL kbykpiv, ps, pt, lnkbykpiv
    type(coop_arguments) args
    lnkbykpiv = log(kbykpiv)
    ps = args%r(1) * exp((args%r(2)-1.d0 + args%r(3)/2.d0*lnkbykpiv)*lnkbykpiv)
    pt = args%r(4)*args%r(1)*exp(lnkbykpiv*(args%r(5)))
  end subroutine coop_cosmology_firstorder_standard_power
  

  subroutine coop_cosmology_firstorder_source_get_transfer(source, l, trans)
    class(coop_cosmology_firstorder_source)::source
    COOP_INT::l
    COOP_REAL,dimension(:,:,:)::trans
    COOP_INT::n, ik, idense, itau, ii, nchis
    COOP_REAL::jl, xmin, xmax, dchi
    COOP_REAL,dimension(:),allocatable::kmin, kmax, kfine

    if(size(trans,2).ne. coop_k_dense_fac .or. size(trans, 3).ne. source%nk .or. size(trans, 1) .ne. source%nsrc) call coop_return_error("get_transfer", "wrong size", "stop")
    trans = 0
    call coop_jl_startpoint(l, xmin)
    call coop_jl_check_init(l)
    xmax = coop_jl_zero(l, 7)
    allocate(kmin(source%ntau), kmax(source%ntau), kfine(source%ntau))
    do itau = 1, source%ntau
       kmin(itau) = xmin/source%chi(itau)
       kmax(itau) = xmax/source%chi(itau)
       kfine(itau) = 1.d0/source%dtau(itau) !!k > kfine use fine grid
    enddo

    !$omp parallel do private(ik, itau, idense, jl, ii, dchi, nchis)
    do ik=2, source%nk
       do itau = 1, source%ntau
          if(source%k(ik)*source%chi(itau) .lt. xmin) exit
          if(source%k(ik) .lt. kfine(itau))then
             do idense = 1, coop_k_dense_fac
                jl = coop_jl(l, source%k_dense(idense, ik)*source%chi(itau))
                trans(:,idense, ik) = trans(:,idense, ik) + jl* source%dtau(itau)* COOP_INTERP_SOURCE(source, :, idense, ik, itau)
             enddo
          else
             if(source%k(ik) .gt. kmin(itau) .and. source%k_dense(1,ik) .lt. kmax(itau))then            
                do idense = 1, coop_k_dense_fac
                   nchis = ceiling(source%k_dense(idense, ik)*source%dtau(itau))
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


    


  subroutine coop_cosmology_firstorder_source_get_Cls(source, l, Cls)
    class(coop_cosmology_firstorder_source)::source
    COOP_INT l
    COOP_REAL,dimension(coop_num_cls),intent(OUT)::Cls
    COOP_REAL,dimension(:,:,:),allocatable::trans
    COOP_INT:: i
    COOP_REAL::tmp
    Cls = 0.d0
    select case(source%m)
    case(0)
       allocate(trans(source%nsrc, coop_k_dense_fac, source%nk))
       call source%get_transfer(l, trans)
       Cls(coop_index_ClTT) = sum(source%ws_dense * trans(1, :, :)**2)*coop_4pi
       if(source%nsrc .ge. 2)then
          Cls(coop_index_ClTE) = sqrt((l+2.d0)*(l+1.d0)*l*(l-1.d0))*sum(source%ws_dense * trans(1, :, :)*trans(2,:,:))*coop_4pi
          Cls(coop_index_ClEE) = (l+2.d0)*(l+1.d0)*l*(l-1.d0)*sum(source%ws_dense * trans(2,:,:)**2)*coop_4pi
       endif
       if(source%nsrc.ge.3)then
          Cls(coop_index_ClLenLen) = sum(source%ws_dense * trans(3, :, :)**2)*coop_4pi
          Cls(coop_index_ClTLen) = sum(source%ws_dense * trans(1, :, :) * trans(3,:,:))*coop_4pi
       endif
       deallocate(trans)
    case(1)
       call coop_tbw("get_Cls: vector")
    case(2)
       call coop_tbw("get_Cls: tensor")
    case default
       call coop_return_error("get_Cls", "unknown m = "//trim(coop_num2str(source%m)), "stop")
    end select
  end subroutine coop_cosmology_firstorder_source_get_Cls


  subroutine coop_cosmology_firstorder_source_get_All_Cls(source, lmin, lmax, Cls)
    class(coop_cosmology_firstorder_source)::source
    COOP_INT lmin, lmax, l, nc, i
    COOP_REAL,dimension(coop_num_cls, lmin:lmax),intent(OUT)::Cls
    COOP_REAL, dimension(:),allocatable::ls_computed
    COOP_REAL, dimension(:,:),allocatable::Cls_computed, Cls2_computed
    COOP_REAL::Cls_tmp(coop_num_cls)
    COOP_REAL, parameter::norm = 1.d10
    l = lmin 
    nc = 1
    do
       call next_l()
       nc = nc + 1
       if(l.ge.lmax)then
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
       if(l.ge.lmax)then
          ls_computed(i) = dble(lmax)
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
          do l = lmin, lmax
             call coop_splint(nc, ls_computed, Cls_Computed(:, i), Cls2_computed(:, i), dble(l), Cls(i, l))
             Cls(i, l) = Cls(i, l)/(l*(l+1.d0)*norm)
          enddo
       endif
    enddo
    !$omp end parallel do
    deallocate(ls_computed, Cls_computed,Cls2_computed)

  contains

    subroutine next_l()
      l = l + min(32, 12 + l/70, max(1, l/4))
    end subroutine next_l

  end subroutine coop_cosmology_firstorder_source_get_All_Cls



  subroutine coop_cosmology_firstorder_compute_source(this, m)
    class(coop_cosmology_firstorder)::this
    COOP_INT m, ik, itau, is

    call this%init_source(m)

    !$omp parallel do
    do ik = 1, this%source(m)%nk
       call this%compute_source_k(this%source(m), ik)
    enddo
    !$omp end parallel do

    !$omp parallel do private(is, itau)
    do itau = 1, this%source(m)%ntau
       do is = 1, this%source(m)%nsrc
          call coop_naturalspline_uniform(this%source(m)%nk, this%source(m)%s(is, :, itau), this%source(m)%s2(is, :, itau))
       enddo
    enddo
    !$omp end parallel do

  end subroutine coop_cosmology_firstorder_compute_source

  subroutine coop_cosmology_firstorder_compute_source_k(this, source, ik)
    class(coop_cosmology_firstorder)::this
    type(coop_cosmology_firstorder_source)::source
    COOP_INT ik, nvars, nw, itau, iq
    type(coop_pert_object) pert
    COOP_REAL, dimension(:,:),allocatable::w
    COOP_REAL c(24)
    COOP_INT ind, i
    COOP_REAL tau_ini, lna
    tau_ini = min(coop_initial_condition_epsilon/source%k(ik), this%conformal_time(this%a_eq*coop_initial_condition_epsilon), source%tau(1)*0.999d0)
    call this%set_initial_conditions(pert, m = source%m, k = source%k(ik), tau = tau_ini)
    lna = log(this%aoftau(tau_ini))
    call coop_cosmology_firstorder_equations(pert%ny+1, lna, pert%y, pert%yp, this, pert)

    ind = 1
    c = 0.d0
    nvars = pert%ny + 1
    nw = nvars
    allocate(w(nw, 9))
    w = 0.d0
    iq = 1
    do itau = 1, source%ntau
       call coop_dverk_firstorder(nvars, coop_cosmology_firstorder_equations, this, pert, lna,   pert%y, source%lna(itau),  coop_cosmology_firstorder_ode_accuracy, ind, c, nw, w)
       call coop_cosmology_firstorder_equations(pert%ny+1, lna, pert%y, pert%yp, this, pert)

       !!------------------------------------------------------------
       !!forcing v_b - v_g to tight coupling approximations
       !!you can remove this, no big impact
       if(pert%tight_coupling)then
          pert%O1_V_B = (pert%O1_V_B * pert%rhoa2_b + (pert%O1_T(1)/4.d0+pert%slip)* pert%rhoa2_g)/(pert%rhoa2_b+pert%rhoa2_g) 
          pert%O1_T(1) = (pert%O1_V_B - pert%slip)*4.d0
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
          nw = nvars
          allocate(w(nw, 9))
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
             nw = nvars
             allocate(w(nw, 9))
             iq = iq + 1
             if(iq .gt.coop_pert_default_nq) exit
          enddo
       endif
    enddo
    deallocate(w)
    call pert%free()
  end subroutine coop_cosmology_firstorder_compute_source_k



!!============================================================================

  subroutine coop_cosmology_firstorder_set_klms(this)
    class(coop_cosmology_firstorder)::this
    COOP_INT k, l, m, s
    if(this%klms_done) return
    do s = 0, coop_pert_default_smax
       do m = 0,coop_pert_default_mmax
          do l = 0, coop_pert_default_lmax
             if(m.ge.l .or. s.ge.l)then
                this%klms(l,m,s) = 0.d0
             else
                this%klms(l, m, s) = sqrt(dble(l**2-m**2)*dble(l**2-s**2))/l
             endif
             this%klms_by_2lp1(l, m, s) = this%klms(l, m, s)/(2*l+1.d0)
             this%klms_by_2lm1(l, m, s) = this%klms(l, m, s)/(2*l-1.d0)
          enddo
       enddo
    enddo
    this%klms_done = .true.
  end subroutine coop_cosmology_firstorder_set_klms


  function coop_cosmology_firstorder_cs2bofa(this, a) result(cs2bofa)
    class(coop_cosmology_firstorder)::this
    COOP_REAL a, cs2bofa
    cs2bofa = this%species(this%index_baryon)%fcs2%eval(a)
  end function coop_cosmology_firstorder_cs2bofa

  function coop_cosmology_firstorder_Tbofa(this, a) result(Tbofa)
    class(coop_cosmology_firstorder)::this
    COOP_REAL a, Tbofa
    Tbofa = this%Tb%eval(a)
  end function coop_cosmology_firstorder_Tbofa


  subroutine coop_cosmology_firstorder_source_free(this)
    class(coop_cosmology_firstorder_source)::this
    if(allocated(this%k))deallocate(this%k, this%dk, this%index_tc_off)
    this%nk = 0
    if(allocated(this%tau))deallocate(this%tau, this%chi, this%dtau, this%a, this%tauc, this%lna)
    this%ntau = 0
    if(allocated(this%s))deallocate(this%s)
  end subroutine coop_cosmology_firstorder_source_free

  subroutine coop_cosmology_firstorder_free(this)
    class(coop_cosmology_firstorder)::this
    COOP_INT m, i
    do m=0,2
       call this%source(m)%free()
    enddo
    call this%Ps%free()
    call this%Pt%free()
    call this%Xe%free()
    call this%eKappa%free()
    call this%vis%free()
    call this%Tb%free()
    call this%fdis%free
    call this%ftime%free
    call this%faoftau%free
    do i= 1, this%num_species
       call this%species(i)%free
    enddo
    this%num_species = 0
    this%omega_k_value = 1.
    this%need_setup_background = .true.
  end subroutine coop_cosmology_firstorder_free


  subroutine coop_cosmology_firstorder_set_power(this, power, args)
!!the format of primordial power
!!subroutine power(k/k_pivot, ps, pt, cosmology, args)
!!input kMpc and args, output ps and pt    
    integer,parameter::n = 4096
    class(coop_cosmology_firstorder)::this
    external power
    type(coop_arguments) args
    COOP_REAL k(n), ps(n), pt(n)
    COOP_INT i
    this%k_pivot = this%kMpc_pivot/this%H0Mpc()
    call coop_set_uniform(n, k, coop_power_lnk_min, coop_power_lnk_max)
    k = exp(k)
    do i=1, n
       call power(k(i)/this%k_pivot, ps(i), pt(i), this, args)
    end do
    call this%ps%init(n = n, xmin = k(1), xmax = k(n), f = ps, xlog = .true., ylog = .true., fleft = ps(1), fright = ps(n), check_boundary = .false.)
    if(any(pt.eq.0.d0))then
       call this%pt%init(n = n, xmin = k(1), xmax = k(n), f = pt, xlog = .true., ylog = .false., fleft = pt(1), fright = pt(n), check_boundary = .false.)
    else
       call this%pt%init(n = n, xmin = k(1), xmax = k(n), f = pt, xlog = .true., ylog = .true., fleft = pt(1), fright = pt(n), check_boundary = .false.)
    endif
    this%As = this%psofk(this%k_pivot)
    this%ns = this%ps%derivative_bare(log(this%k_pivot)) + 1.d0
    this%nrun = this%ps%derivative2_bare(log(this%k_pivot))
    this%r =  this%ptofk(this%k_pivot)/this%As
    if(this%r .eq. 0.d0)then
       this%nt = 0.d0
    else
       this%nt = this%pt%derivative(this%k_pivot)*this%k_pivot/(this%r*this%As)
    endif
  end subroutine coop_cosmology_firstorder_set_power

  subroutine coop_cosmology_firstorder_set_xe(this)
    class(coop_cosmology_firstorder)::this
    integer i, m
    integer, parameter::n = 4096
    COOP_REAL ekappa(n), a(n), dkappadinva(n), dkappadtau(n), vis(n)
    call this%setup_background()
    call this%set_klms()
    this%index_baryon = this%index_of("Baryon", .true.)
    this%index_radiation = this%index_of("Radiation", .true.)
    this%index_cdm = this%index_of("CDM", .true.)
    this%index_Nu = this%index_of("Massless Neutrinos", .true.)
    this%index_massiveNu = this%index_of("Massive Neutrinos")
    this%index_de = this%index_of("Dark Energy", .true.)

    this%Omega_b = O0_BARYON(this)%Omega 
    this%Omega_c = O0_CDM(this)%Omega 
    this%Omega_g = O0_RADIATION(this)%Omega 
    if(this%index_massiveNu .ne. 0)then
       call coop_fermion_get_lnam(log(O0_MASSIVENU(this)%Omega/O0_MASSIVENU(this)%Omega_massless), this%mnu_by_Tnu)
       this%mnu_by_Tnu = exp(this%mnu_by_Tnu)
       this%Omega_nu = O0_NU(this)%Omega + O0_MASSIVENU(this)%Omega_massless
    else
       this%Omega_nu = O0_NU(this)%Omega 
       this%mnu_by_Tnu = 0.d0
    endif

    if(abs(this%Omega_m/( this%Omega_b + this%Omega_c ) - 1.d0) .gt. 1.d-3)then
       call coop_feedback("warning: nonstandard matter component.", 1)
    endif


    if(abs(this%Omega_r/( this%Omega_g + this%Omega_nu ) - 1.d0) .gt. 1.d-3)then
       call coop_feedback("warning: cannot accurately determine early-time radiation fraction.", 1)
    endif
    this%tau_eq = this%conformal_time(this%a_eq)

    this%dkappadtau_coef = this%species(this%index_baryon)%Omega * this%h() * coop_SI_sigma_thomson * (coop_SI_rhocritbyh2/coop_SI_c**2) * coop_SI_hbyH0 * coop_SI_c/ coop_SI_m_H * (1.d0 - this%YHe())
    if(this%do_reionization)then
       this%reionFrac = 1.d0 + this%YHe()/(coop_m_He_by_m_H * (1.d0-  this%YHe()))
    else
       this%reionFrac = 0.d0
    endif

    call this%set_zre_from_optre()
    call coop_recfast_get_xe(this, this%xe, this%Tb, this%reionFrac, this%zre, this%deltaz)
    call coop_set_uniform(n, a, log(coop_visibility_amin), log(coop_scale_factor_today))

    a = exp(a)
    !$omp parallel do
    do i=1, n
       dkappadtau(i) = this%dkappadtau(a(i))
       dkappadinva(i)= - dkappadtau(i)/this%Hratio(a(i))
    enddo
    !$omp end parallel do
    ekappa(n) = 0.d0
    do i=n-1,1, -1
       ekappa(i) = ekappa(i+1) + (dkappadinva(i+1) + dkappadinva(i))*(0.5d0/a(i) - 0.5d0/a(i+1))
       if(ekappa(i) .lt. -200.d0)then
          ekappa(1:i) = -200.d0
          exit
       endif
    enddo
    ekappa = exp(ekappa)
    vis = dkappadtau*ekappa
    call this%ekappa%init(n = n, xmin = a(1), xmax = a(n), f = ekappa, xlog = .true., ylog = .true., fleft = ekappa(1), fright = ekappa(n), check_boundary = .false.)
    call this%vis%init(n = n, xmin = a(1), xmax = a(n), f = vis, xlog = .true., ylog = .true., fleft = vis(1), fright = vis(n), check_boundary = .false.)
    this%zrecomb = 1.d0/this%vis%maxloc() - 1.d0
    this%distlss = this%comoving_distance(1.d0/(1.d0+this%zrecomb))
    this%tau0 = this%conformal_time(coop_scale_factor_today)
    
  end subroutine coop_cosmology_firstorder_set_xe

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


  subroutine coop_cosmology_firstorder_set_source_tau(this, source, n, weight)
    COOP_REAL, parameter::incr = 1.05d0
    COOP_INT::nbuffer = 35
    COOP_INT n, i, j
    COOP_REAL  top(n), dtop, amax, amin, amid, topmax, topmin, topmid, weight, taucut, acut
    class(coop_cosmology_firstorder_source)::source
    class(coop_cosmology_firstorder)::this
    if(n.lt. nbuffer+10) stop "n is too small for set_source_tau"
    source%ntau = n
    if(allocated(source%a))then
       if(size(source%a).ne.n)then
          deallocate(source%a, source%tau, source%chi, source%dtau, source%tauc, source%lna)
          allocate(source%a(n), source%tau(n), source%dtau(n), source%chi(n), source%tauc(n), source%lna(n))
       endif
    else
       allocate(source%a(n), source%tau(n), source%dtau(n), source%chi(n),  source%tauc(n), source%lna(n))
    endif
    amin = 1.d0/(1.d0+this%zrecomb)
    do while(this%visofa(amin).gt. 1.)
       amin = amin/1.02d0
    enddo
    topmin = this%tau0*(1.d0-weight)*this%ekappaofa(amin)+weight*this%tauofa(amin)
    dtop = (this%tau0 -topmin) / (n - nbuffer)
    do i= nbuffer+1, n
       top(i) = topmin + dtop*(i-nbuffer-0.5d0)
    enddo
    do i=nbuffer+1, n
       amax = amin*incr
       topmax = this%ekappaofa(amax)*(1.d0-weight)*this%tau0 + weight*this%tauofa(amax)
       do while( topmax .lt. top(i))
          amin = amax
          topmin = topmax
          amax = amax*incr
          if(amax .ge. coop_scale_factor_today)then
             amax = coop_scale_factor_today
             topmax = this%tau0
             exit
          endif
          topmax = this%ekappaofa(amax)*(1.d0-weight)*this%tau0 + weight*this%tauofa(amax)
       enddo
       do while((amax-amin)/amin .gt. 1.d-7)
          amid = (amin+amax)/2.d0
          topmid = this%ekappaofa(amid)*(1.d0-weight)*this%tau0 + weight*this%tauofa(amid)
          if(topmid .gt. top(i))then
             amax = amid
             topmax = topmid
          else
             amin = amid
             topmin = topmid
          endif
       enddo      
       source%a(i) = amid + (top(i)-topmid)/((weight+(1.d0-weight)*this%tau0*this%visofa(amid))/this%dadtau(amid))
       source%tau(i) = this%tauofa(source%a(i))
       source%dtau(i) = dtop/(weight+(1.d0-weight)*this%tau0*this%visofa(source%a(i)))
       source%tauc(i) = this%taucofa(source%a(i))
    enddo

    do i=nbuffer, 1, -1
       source%a(i) = source%a(i+1)*0.976d0
       source%tau(i) = this%conformal_time(source%a(i))
       source%tauc(i) = this%taucofa(source%a(i))
    enddo
    source%dtau(1) = source%tau(2)/2.d0

    do i=2, nbuffer+1
       source%dtau(i) = (source%tau(i+1) - source%tau(i-1))/2.d0
    enddo
    source%chi = this%tau0 - source%tau
    source%lna = log(source%a)
    source%index_tc_max = 1
    do while( source%index_tc_max .lt. source%ntau - 1 .and. source%tauc(source%index_tc_max+1) .le. coop_cosmology_firstorder_tc_cutoff * source%tau(source%index_tc_max+1) )
       source%index_tc_max = source%index_tc_max+1
    enddo
  end subroutine coop_cosmology_firstorder_set_source_tau


  subroutine coop_cosmology_firstorder_set_source_k(this, source, n, weight)
    class(coop_cosmology_firstorder)::this
    class(coop_cosmology_firstorder_source)::source
    COOP_INT n, i, iq, j
    COOP_REAL weight, dkop_dense
    this%k_pivot = this%kMpc_pivot/this%H0Mpc()
    source%kweight = weight
    source%nk = n
    if(allocated(source%k))then
       if(size(source%k).ne.n)then
          deallocate(source%k, source%dk, source%index_tc_off, source%kop)
          allocate(source%k(n), source%dk(n), source%index_tc_off(n), source%kop(n))
       endif
       if(size(source%k_dense, 2).ne.n .or. size(source%k_dense, 1).ne. coop_k_dense_fac )then
          deallocate(source%k_dense, source%ws_dense, source%wt_dense, source%dk_dense)
          allocate(source%k_dense(coop_k_dense_fac, n), source%ws_dense(coop_k_dense_fac, n), source%wt_dense(coop_k_dense_fac, n), source%dk_dense(coop_k_dense_fac, n) )
       endif
    else
       allocate(source%k(n), source%dk(n), source%index_tc_off(n), source%kop(n))
       allocate(source%k_dense(coop_k_dense_fac, n), source%ws_dense(coop_k_dense_fac, n), source%wt_dense(coop_k_dense_fac, n), source%dk_dense(coop_k_dense_fac, n) )
    endif
    do i=1, coop_k_dense_fac
       source%a_dense(i) = dble(i)/coop_k_dense_fac
       source%b_dense(i) =  1.d0 - source%a_dense(i)
       source%a2_dense(i) = source%a_dense(i)*(source%a_dense(i)**2-1.d0)
       source%b2_dense(i) = source%b_dense(i)*(source%b_dense(i)**2-1.d0)
    enddo


    source%kmin = 0.25d0/this%distlss
    select case(source%m)
    case(0)
       source%kmax = (coop_scalar_lmax*2.d0)/this%distlss
    case(1)
       source%kmax = (coop_vector_lmax*2.d0)/this%distlss
    case(2)
       source%kmax = (coop_tensor_lmax*2.d0)/this%distlss
    case default
       write(*,*) "Error: m = ", source%m
       stop "source spin must be 0, 1, 2"
    end select

    call source%k2kop(source%kmin, source%kopmin)
    call source%k2kop(source%kmax, source%kopmax)
    source%dkop = (source%kopmax-source%kopmin)/(n-1)
    call coop_set_uniform(n, source%kop, source%kopmin, source%kopmax)
    dkop_dense = source%dkop/coop_k_dense_fac


    !$omp parallel do private(i, j)
    do i=1, n
       call source%kop2k(source%kop(i), source%k(i), source%dkop,  source%dk(i))
       do j=1, coop_k_dense_fac
          call source%kop2k(source%kop(i)+(j-coop_k_dense_fac)*dkop_dense, source%k_dense(j,i), dkop_dense,  source%dk_dense(j, i))
          source%ws_dense(j, i) = source%dk_dense(j, i)*this%psofk(source%k_dense(j, i))/source%k_dense(j, i)
          source%wt_dense(j, i) = source%dk_dense(j, i)*this%ptofk(source%k_dense(j, i))/source%k_dense(j, i)
       enddo
    enddo
    !$omp end parallel do

    !!set index_tc_off
    if(source%ntau .gt. 0 .and. allocated(source%tauc))then
       source%index_tc_off = 1
       !$omp parallel do
       do i = 1, n
          do while(source%index_tc_off(i) .lt. source%index_tc_max .and. source%tauc(source%index_tc_off(i)+1)*source%k(i) .le. coop_cosmology_firstorder_tc_cutoff)
             source%index_tc_off(i)= source%index_tc_off(i)+1
          enddo
       enddo
       !$omp end parallel do
    else
       call coop_return_error("set_source_k", "You need to call set_source_tau before calling set_source_k", "stop")
    endif
    
    !!set index_massivenu_on
    if(this%mnu_by_Tnu*source%a(source%ntau) .lt. coop_pert_massivenu_threshold(1))then
       source%index_massivenu_on = source%ntau + 1
       return
    endif
    source%index_massivenu_on(coop_pert_default_nq) = source%ntau + 1
    iq = coop_pert_default_nq
    do while(this%mnu_by_Tnu*source%a(source%index_massivenu_on(iq)-1) .gt. coop_pert_massivenu_threshold(iq))
       source%index_massivenu_on(iq) = source%index_massivenu_on(iq)-1
       if(source%index_massivenu_on(iq) .le. 1)exit
    enddo
    do iq = coop_pert_default_nq - 1, 1, -1
       source%index_massivenu_on(iq) =  source%index_massivenu_on(iq+1)
       do while(this%mnu_by_Tnu*source%a(source%index_massivenu_on(iq)-1) .gt. coop_pert_massivenu_threshold(iq))
          source%index_massivenu_on(iq) = source%index_massivenu_on(iq)-1
          if(source%index_massivenu_on(iq) .le. 1)exit
       enddo
    enddo
  end subroutine coop_cosmology_firstorder_set_source_k

  subroutine coop_cosmology_firstorder_source_k2kop(source, k, kop)
    class(coop_cosmology_firstorder_source)::source
    COOP_REAL::k, kop
    kop = k**coop_source_k_index*source%kweight + log(k)
  end subroutine coop_cosmology_firstorder_source_k2kop

  subroutine coop_cosmology_firstorder_source_kop2k(source, kop, k, dkop, dk)
    class(coop_cosmology_firstorder_source)::source
    COOP_REAL kop, k, kalpha
    COOP_REAL, optional::dkop, dk
    call coop_source_kop2k_noindex(kop*coop_source_k_index, source%kweight*coop_source_k_index, kalpha)
    k = kalpha**(1.d0/coop_source_k_index)
    if(present(dkop) .and. present(dk))then
       dk = dkop * k / (1.d0 + source%kweight * coop_source_k_index * kalpha )
    endif
  end subroutine coop_cosmology_firstorder_source_kop2k

  subroutine coop_source_kop2k_noindex(kop, weight, k)
    !!solve the equation k * weight + log(k) = kop
    COOP_REAL kop, weight, k, kmin, kmax, kmid
    if(kop .gt. weight)then
       kmax = kop/weight
       kmin = (kop+1.d0)/(weight+1.d0)
    else
       kmax = exp(kop)
       kmin = max(kop/weight, 1.d-3)
    endif
    do while((kmax-kmin)/kmin .gt. 1.d-8)
       kmid = (kmax+kmin)/2.d0
       if(kmid*weight + log(kmid) .gt. kop)then
          kmax = kmid
       else
          kmin = kmid
       endif
    enddo
    k = (kmax+kmin)/2.d0
  end subroutine coop_source_kop2k_noindex
  
  function coop_cosmology_firstorder_psofk(this, k) result(ps)
    COOP_REAL k, ps
    class(coop_cosmology_firstorder)::this
    ps = this%ps%eval(k)
  end function coop_cosmology_firstorder_psofk


  function coop_cosmology_firstorder_ptofk(this, k) result(pt)
    COOP_REAL k, pt
    class(coop_cosmology_firstorder)::this
    pt = this%pt%eval(k)
  end function coop_cosmology_firstorder_ptofk



  subroutine coop_cosmology_firstorder_init_source(this, m)
    class(coop_cosmology_firstorder)::this
    COOP_INT :: m
    this%source(m)%m = m
    this%source(m)%nsrc = coop_num_sources(m)
    call this%set_source_tau(this%source(m), coop_source_tau_n(m), coop_source_tau_weight(m))

    call this%set_source_k(this%source(m), coop_source_k_n(m), coop_source_k_weight(m))

    if(allocated(this%source(m)%s))then
       if(size(this%source(m)%s, 1) .ne. this%source(m)%ntau .or. size(this%source(m)%s, 2) .ne. this%source(m)%nk)then
          deallocate(this%source(m)%s, this%source(m)%s2)
          allocate(this%source(m)%s(this%source(m)%nsrc, this%source(m)%nk , this%source(m)%ntau), this%source(m)%s2(this%source(m)%nsrc, this%source(m)%nk, this%source(m)%ntau) )
       endif
    else
       allocate(this%source(m)%s(this%source(m)%nsrc, this%source(m)%nk , this%source(m)%ntau), this%source(m)%s2(this%source(m)%nsrc, this%source(m)%nk, this%source(m)%ntau) )

    endif
  end subroutine coop_cosmology_firstorder_init_source




  function coop_cosmology_firstorder_Clzetazeta_at_R(this, l, r) result(Cl)
    class(coop_cosmology_firstorder)::this
    !!this computes 4 \pi \int_0^\infty |j_l(kr)|^2 (k^3P(k)/(2\pi^2)) d\ln k
    COOP_REAL::Cl, r
    COOP_INT l
    Cl = (coop_pi**2/2.d0)* this%psofk(1.d0/r) * 2.d0 ** ( this%ns - 1.d0) * exp(log_gamma(3.d0 - this%ns )+log_gamma(l + (this%ns - 1.d0)/2.d0)-log_gamma(l+(5.d0-this%ns)/2.d0)-2.d0*log_gamma((4.d0-this%ns)/2.d0))
  end function Coop_cosmology_firstorder_Clzetazeta_at_R

  function coop_cosmology_firstorder_Clzetazeta(this, l, r1, r2) result(Cl)
    !!this computes 4 \pi \int_0^\infty j_l(kr_1) j_l(kr_2) (k^3P(k)/(2\pi^2)) d\ln k
    class(coop_cosmology_firstorder)::this
    COOP_REAL Cl, r1
    COOP_REAL, optional::r2
    COOP_INT l
    if(.not. present(r2))then
       Cl = this%Clzetazeta_at_R(l, r1)
    else
       if(r1.eq.r2)then
          Cl = this%Clzetazeta_at_R(l, r1)
       else
          cl = this%psofk(2.d0/(r1+r2))*coop_4pi*coop_sphericalBesselCross(l,l,r1,r2,this%ns)
       endif
    endif
  end function Coop_cosmology_firstorder_Clzetazeta


end module coop_firstorder_mod

