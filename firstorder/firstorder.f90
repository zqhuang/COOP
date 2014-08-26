module coop_firstorder_mod
  use coop_wrapper_background
  use coop_recfast_mod
  use coop_pertobj_mod
  implicit none
#include "constants.h"

!private

public:: coop_cosmology_firstorder, coop_cosmology_firstorder_source

  COOP_REAL, parameter :: coop_power_lnk_min = log(0.1d0) 
  COOP_REAL, parameter :: coop_power_lnk_max = log(5.d3) 
  COOP_REAL, parameter :: coop_visibility_amin = 3.d-4
  COOP_REAL, parameter :: coop_initial_condition_epsilon = 1.d-5
  COOP_REAL, parameter :: coop_cosmology_firstorder_ode_accuracy = 1.d-7

  COOP_REAL, dimension(0:2), parameter::coop_source_tau_weight = (/ 0.3d0, 0.2d0, 0.1d0 /)
  COOP_INT, dimension(0:2), parameter::coop_source_tau_n = (/ 1000, 850, 800 /)

  COOP_REAL, dimension(0:2), parameter::coop_source_k_weight = (/ 0.2d0, 0.15d0, 0.1d0 /)

  COOP_INT, dimension(0:2), parameter::coop_source_k_n = (/ 160, 120, 80 /)

  COOP_INT, dimension(0:2), parameter::coop_num_sources = (/ 4,  2,  2 /)

  COOP_REAL, parameter::coop_source_k_index = 0.55d0


  type coop_cosmology_firstorder_source
     COOP_INT::m = 0
     COOP_INT::ntau = 0
     COOP_INT::nk = 0
     COOP_REAL, dimension(:),allocatable::k, tau, a, dk, dtau, chi  !!tau is conformal time, chi is comoving distance; in flat case, chi + tau = tau_0
     COOP_REAL, dimension(:,:,:),allocatable::source
   contains
     procedure::free => coop_cosmology_firstorder_source_free
  end type coop_cosmology_firstorder_source
  
  type, extends(coop_cosmology_background) :: coop_cosmology_firstorder
     logical::do_reionization = .true.
     COOP_REAL::zrecomb, distlss, tau0
     COOP_REAL::optre = 0.07d0
     COOP_REAL::zre = 8.d0
     COOP_REAL::deltaz = 0.5d0
     COOP_REAL::kMpc_pivot = 0.05d0
     COOP_REAL::k_pivot
     COOP_REAL::dkappadtau_coef, ReionFrac, Omega_m, Omega_r, a_eq, tau_eq
     type(coop_function)::Ps, Pt, Xe, ekappa, vis
     type(coop_cosmology_firstorder_source),dimension(0:2)::saved_source
     COOP_INT::index_baryon, index_cdm, index_radiation, index_nu, index_massiveNu, index_de
   contains
     procedure:: set_source_tau => coop_cosmology_firstorder_set_source_tau
     procedure:: set_source_k => coop_cosmology_firstorder_set_source_k
     procedure:: set_power => coop_cosmology_firstorder_set_power
     procedure:: set_xe => coop_cosmology_firstorder_set_xe
     procedure:: set_zre_from_optre => coop_cosmology_firstorder_set_zre_from_optre
     procedure:: set_optre_from_zre => coop_cosmology_firstorder_set_optre_from_zre
     procedure::set_tight_coupling => coop_cosmology_firstorder_set_tight_coupling
     procedure::set_late_approx => coop_cosmology_firstorder_set_late_approx
     procedure::set_initial_conditions => coop_cosmology_firstorder_set_initial_conditions
     procedure:: xeofa => coop_cosmology_firstorder_xeofa
     procedure:: dxeda => coop_cosmology_firstorder_dxeda
     procedure:: dkappada => coop_cosmology_firstorder_dkappada
     procedure:: dkappadtau => coop_cosmology_firstorder_dkappadtau
     procedure:: tauc => coop_cosmology_firstorder_tauc
     procedure:: dot_tauc => coop_cosmology_firstorder_dot_tauc
     procedure:: visofa => coop_cosmology_firstorder_visofa
     procedure:: ekappaofa => coop_cosmology_firstorder_ekappaofa
     procedure:: psofk => coop_cosmology_firstorder_psofk
     procedure:: ptofk => coop_cosmology_firstorder_ptofk
     procedure:: free => coop_cosmology_firstorder_free
  end type coop_cosmology_firstorder


contains

  subroutine coop_cosmology_firstorder_source_free(this)
    class(coop_cosmology_firstorder_source)::this
    if(allocated(this%k))deallocate(this%k, this%dk)
    this%nk = 0
    if(allocated(this%tau))deallocate(this%tau, this%chi, this%dtau, this%a)
    this%ntau = 0
    if(allocated(this%source))deallocate(this%source)
  end subroutine coop_cosmology_firstorder_source_free

  subroutine coop_cosmology_firstorder_free(this)
    class(coop_cosmology_firstorder)::this
    COOP_INT m, i
    do m=0,2
       call this%saved_source(m)%free()
    enddo
    call this%Ps%free()
    call this%Pt%free()
    call this%Xe%free()
    call this%eKappa%free()
    call this%vis%free()
    call this%fdis%free
    call this%ftime%free
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
    call this%pt%init(n = n, xmin = k(1), xmax = k(n), f = pt, xlog = .true., ylog = .true., fleft = pt(1), fright = pt(n), check_boundary = .false.)

  end subroutine coop_cosmology_firstorder_set_power

  subroutine coop_cosmology_firstorder_set_xe(this)
    class(coop_cosmology_firstorder)::this
    integer i, m
    integer, parameter::n = 4096
    COOP_REAL ekappa(n), a(n), dkappadinva(n), dkappadtau(n), vis(n)
    this%index_baryon = this%index_of("Baryon")
    this%index_radiation = this%index_of("Radiation")
    if(this%index_baryon .eq. 0 .or. this%index_radiation .eq. 0)then
       stop "Set_Xe: for a model without baryon, you cannot set up Xe"
    endif
    this%Omega_m = this%species(this%index_baryon)%Omega
    this%Omega_r = this%species(this%index_radiation)%Omega
    this%index_cdm = this%index_of("CDM")
    if(this%index_cdm.eq.0)then
       call coop_feedback("set_xe: warning: no CDM")
    else
       this%Omega_m  = this%Omega_m  + this%species(this%index_cdm)%Omega
    endif
    this%index_Nu = this%index_of("Massless Neutrinos")
    if(this%index_Nu.eq.0)then
       call coop_feedback("set_xe: warning: no massless Neutrinos")
    else
       this%Omega_r = this%Omega_r + this%species(this%index_Nu)%Omega
    endif
    this%index_massiveNu = this%index_of("Massive Neutrinos")
    if(this%index_massiveNu.eq.0)then
       call coop_feedback("set_xe: warning: no massive Neutrinos", 2)
    else
       this%Omega_r = this%Omega_r + this%species(this%index_massiveNu)%Omega_massless
    endif
    this%index_de = this%index_of("Dark Energy")
    if(this%index_de.eq.0)then
       call coop_feedback("set_xe: warning: no dark energy")
    endif

    this%a_eq = this%Omega_r / this%Omega_m
    this%tau_eq = this%conformal_time(this%a_eq)

    this%dkappadtau_coef = this%species(this%index_baryon)%Omega * this%h() * coop_SI_sigma_thomson * (coop_SI_rhocritbyh2/coop_SI_c**2) * coop_SI_hbyH0 * coop_SI_c/ coop_SI_m_H * (1.d0 - this%YHe())
    if(this%do_reionization)then
       this%reionFrac = 1.d0 + this%YHe()/(coop_m_He_by_m_H * (1.d0-  this%YHe()))
    else
       this%reionFrac = 0.d0
    endif
    call this%set_zre_from_optre()

    call coop_recfast_get_xe(this, this%xe, this%reionFrac, this%zre, this%deltaz)
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
    if(a .lt. 1.1d-4)then
       xe = this%xe%fleft
    else
       xe = this%xe%eval(a)
    endif
  end function coop_cosmology_firstorder_xeofa

  function coop_cosmology_firstorder_dxeda(this, a) result(dxeda)
    COOP_REAL a, dxeda
    class(coop_cosmology_firstorder)::this
    if(a .lt. 1.1d-4)then
       dxeda = 0.d0 
    else
       dxeda = this%xe%derivative(a)
    endif
  end function coop_cosmology_firstorder_dxeda


  function coop_cosmology_firstorder_ekappaofa(this, a) result(ekappa)
    COOP_REAL a, ekappa
    class(coop_cosmology_firstorder)::this
    if(a .lt. 1.1d-4)then
       ekappa = 0.d0
    else
       ekappa = this%ekappa%eval(a)
    endif    
  end function coop_cosmology_firstorder_ekappaofa

  function coop_cosmology_firstorder_visofa(this, a) result(vis)
    COOP_REAL a, vis
    class(coop_cosmology_firstorder)::this
    if(a .lt. 1.1d-4)then
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


  function coop_cosmology_firstorder_tauc(this, a) result(tauc)
    class(coop_cosmology_firstorder)::this    
    COOP_REAL a, tauc
    tauc =   a**2/(this%dkappadtau_coef * this%xeofa(a))
  end function coop_cosmology_firstorder_tauc


  function coop_cosmology_firstorder_dot_tauc(this, a) result(dtauc)
    class(coop_cosmology_firstorder)::this    
    COOP_REAL a, dtauc
    dtauc =  a*(2.d0 - a*this%dxeda(a)/this%xeofa(a))/this%xeofa(a)
  end function coop_cosmology_firstorder_dot_tauc



  function coop_cosmology_firstorder_dkappada(this, a) result(dkappada)
    class(coop_cosmology_firstorder)::this    
    COOP_REAL a, dkappada
    dkappada = this%dkappadtau(a) / this%dadtau(a)
  end function coop_cosmology_firstorder_dkappada


  subroutine coop_cosmology_firstorder_set_source_tau(this, source, n, weight)
    COOP_REAL, parameter::incr = 1.05d0
    COOP_INT n, i, j
    COOP_REAL  top(n), dtop, amax, amin, amid, topmax, topmin, topmid, weight
    class(coop_cosmology_firstorder_source)::source
    class(coop_cosmology_firstorder)::this
    if(n.lt. 10) stop "n is too small for set_source_tau"
    source%ntau = n
    if(allocated(source%a))then
       if(size(source%a).ne.n)then
          deallocate(source%a, source%tau, source%chi, source%dtau)
          allocate(source%a(n), source%tau(n), source%dtau(n), source%chi(n))
       endif
    else
       allocate(source%a(n), source%tau(n), source%dtau(n), source%chi(n))
    endif
    dtop = this%tau0 / n
    !$omp parallel do
    do i=1, n
       top(i) = dtop * (i-0.5d0)
    enddo
    !$omp end parallel do
    amin = coop_visibility_amin
    topmin = this%ekappaofa(amin)**weight*this%conformal_time(amin)
    do i=1, n
       amax = amin*incr
       topmax = this%ekappaofa(amax)**weight*this%conformal_time(amax)
       do while( topmax .lt. top(i))
          amin = amax
          topmin = topmax
          amax = amax*incr
          if(amax .ge. coop_scale_factor_today)then
             amax = coop_scale_factor_today
             topmax = this%tau0
             exit
          endif
          topmax = this%ekappaofa(amax)**weight*this%conformal_time(amax)       
       enddo
       do while((amax-amin)/amin .gt. 1.d-6)
          amid = (amin+amax)/2.d0
          topmid = this%ekappaofa(amid)**weight*this%conformal_time(amid)       
          if(topmid .gt. top(i))then
             amax = amid
             topmax = topmid
          else
             amin = amid
             topmin = topmid
          endif
       enddo
       source%a(i) = (amin+amax)/2.d0
       source%tau(i) = this%conformal_time(source%a(i))
       source%dtau(i) = dtop/((weight*this%dkappadtau(source%a(i))*source%tau(i)+1.d0)*this%ekappaofa(source%a(i))**weight)
    enddo
    source%chi = this%tau0 - source%tau
  end subroutine coop_cosmology_firstorder_set_source_tau
  
  subroutine coop_cosmology_firstorder_set_source_k(this, source, n, weight)
    class(coop_cosmology_firstorder)::this
    class(coop_cosmology_firstorder_source)::source
    COOP_INT n, i
    COOP_REAL weight, kop(n), dkop, kopmin, kopmax, kmin, kmax, kalpha
    this%k_pivot = this%kMpc_pivot/this%H0Mpc()
    source%nk = n
    if(allocated(source%k))then
       if(size(source%k).ne.n)then
          deallocate(source%k, source%dk)
          allocate(source%k(n), source%dk(n))
       endif
    else
       allocate(source%k(n), source%dk(n))
    endif
    kmin = 0.3d0/this%distlss
    select case(source%m)
    case(0)
       kmax = 6.2d3/this%distlss
    case(1)
       kmax = 4.2d3/this%distlss
    case(2)
       kmax = 2.2d3/this%distlss
    case default
       write(*,*) "Error: m = ", source%m
       stop "source spin must be 0, 1, 2"
    end select
    kopmin = kmin**coop_source_k_index*weight + log(kmin)
    kopmax = kmax**coop_source_k_index*weight + log(kmax)
    dkop = (kopmax-kopmin)/(n-1)
    call coop_set_uniform(n, kop, kopmin, kopmax)
    !$omp parallel do private(kalpha)
    do i=1, n
       call coop_source_kop2k_noindex(kop(i)*coop_source_k_index, weight*coop_source_k_index, kalpha)
       source%k(i) = kalpha**(1.d0/coop_source_k_index)
       source%dk(i) = dkop * source%k(i) / (1.d0 + weight * coop_source_k_index * kalpha )
    enddo
    !$omp end parallel do
  end subroutine coop_cosmology_firstorder_set_source_k


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
    pt = this%ps%eval(k)
  end function coop_cosmology_firstorder_ptofk


  subroutine coop_cosmology_firstorder_set_zre_from_optre(this)
    class(coop_cosmology_firstorder)::this
    COOP_REAL zremin, zremax, zremid
    COOP_REAL optremin, optremax, optremid, optre_wanted
    integer iloop
    if(this%ReionFrac .le. 0.d0)then
       this%zre = -1.d0
    endif
    optre_wanted = this%optre

    zremin = 0.d0
    this%zre = zremin
    call this%set_optre_from_zre()
    optremin = this%optre

    zremax = 20.d0
    this%zre = zremax
    call this%set_optre_from_zre()    
    optremax = this%optre

    do while(optre_wanted .gt. optremax)
       zremin = zremax
       optremin = optremax
       zremax = zremax * 1.3d0
       this%zre = zremax
       call this%set_optre_from_zre()    
       optremax = this%optre
    enddo
    iloop = 0
    do while((optremax - optremin) .gt. 1.d-5 .and. iloop .lt. 25)
       zremid = (zremin + zremax)/2.d0
       this%zre = zremid
       call this%set_optre_from_zre()
       optremid = this%optre
       if(optremid .gt. optre_wanted)then
          zremax = zremid
          optremax = optremid
       else
          zremin = zremid
          optremin = optremid
       endif
       iloop = iloop + 1
    enddo
    this%optre = optre_wanted
  end subroutine coop_cosmology_firstorder_set_zre_from_optre

  subroutine coop_cosmology_firstorder_set_optre_from_zre(this)
    class(coop_cosmology_firstorder)::this
    COOP_REAL :: optre, ReionFrac, zre, delta
    type(coop_arguments):: args
    if(this%ReionFrac .le. 0.d0)then
       this%optre = 0.d0
    else
       this%optre = coop_integrate(coop_optre_int, 0.d0, this%zre + this%deltaz * 10.d0, this, args, COOP_REAL_OF(1.e-7))
    endif
  end subroutine coop_cosmology_firstorder_set_optre_from_zre

  function coop_optre_int(z, cosmology, args) result(dkappadz)
    COOP_REAL z, dkappadz
    class(coop_cosmology_firstorder) cosmology
    type(coop_arguments)::args
    dkappadz = cosmology%dkappadtau_coef * coop_reionization_xe(z, cosmology%reionFrac, cosmology%zre, cosmology%deltaz) / cosmology%Hasq(1.d0/(1.d0+z)) 
  end function coop_optre_int


  subroutine coop_cosmology_firstorder_equations(n, a, var, varp, cosmology, pert )
    COOP_INT n
    class(coop_cosmology_firstorder)::cosmology
    class(coop_standard_o1pert)::pert
    COOP_REAL a, var(n), varp(n)
    call pert%read_var(var)
    if(pert%tight_coupling) call cosmology%set_tight_coupling(pert)
    if(pert%late_approx) call cosmology%set_late_approx(pert)
    varp = 0.d0  
  end subroutine coop_cosmology_firstorder_equations

  subroutine coop_cosmology_firstorder_set_tight_coupling(this, pert)
    class(coop_cosmology_firstorder)::this
    class(coop_standard_o1pert)::pert

  end subroutine coop_cosmology_firstorder_set_tight_coupling


  subroutine coop_cosmology_firstorder_set_late_approx(this, pert)
    class(coop_cosmology_firstorder)::this
    class(coop_standard_o1pert)::pert

  end subroutine coop_cosmology_firstorder_set_late_approx


  subroutine coop_cosmology_firstorder_set_initial_conditions(this, pert, k, tau)
    class(coop_cosmology_firstorder)::this
    class(coop_standard_o1pert)::pert
    COOP_REAL tau, k
    pert%k = k
    select case(trim(pert%initial_conditions))
    case("adiabatic")
    case default
       stop "unknown initial conditions"
    end select
  end subroutine coop_cosmology_firstorder_set_initial_conditions


  subroutine coop_cosmology_firstorder_compute_source(this, m)
    class(coop_cosmology_firstorder)::this
    COOP_INT, parameter::n_threads = 8
    type(coop_standard_o1pert),dimension(n_threads):: pert
    COOP_INT m, ik, itau, ithread
    COOP_REAL tau
    this%saved_source(m)%m = m
    call this%set_source_tau(this%saved_source(m), coop_source_tau_n(m), coop_source_tau_weight(m))
    call this%set_source_k(this%saved_source(m), coop_source_k_n(m), coop_source_k_weight(m))
    if(allocated(this%saved_source(m)%source))then
       if(size(this%saved_source(m)%source, 1) .ne. this%saved_source(m)%ntau .or. size(this%saved_source(m)%source, 2) .ne. this%saved_source(m)%nk)then
          deallocate(this%saved_source(m)%source)
          allocate(this%saved_source(m)%source(this%saved_source(m)%ntau, this%saved_source(m)%nk, coop_num_sources(m)))
       endif
    else
       allocate(this%saved_source(m)%source(this%saved_source(m)%ntau, this%saved_source(m)%nk, coop_num_sources(m)))
    endif
    pert%tight_coupling = .true.
    pert%late_approx = .false.
    pert%has_massiveNu = (this%index_massiveNu .ne. 0)
    if(this%index_de .ne. 0)then
       pert%has_de = (this%species(this%index_de)%genre .ne. COOP_SPECIES_LAMBDA)
    else
       pert%has_de = .false.
    endif
    do ithread = 1, n_threads
       call pert(ithread)%init(m)
       call pert(ithread)%set_ode_index()
    enddo
    !$omp parallel do private(ik, itau, ithread, tau)
    do ithread = 1, n_threads
       do ik = ithread, this%saved_source(m)%nk, n_threads
          tau = min(coop_initial_condition_epsilon/this%saved_source(m)%k(ik), this%tau_eq*coop_initial_condition_epsilon)
          call this%set_initial_conditions( pert(ithread), this%saved_source(m)%k(ik), tau)
          call pert(ithread)%read_var()
          do itau = 1, this%saved_source(m)%ntau
             call coop_dverk_firstorder(pert(ithread)%ode_nvars, coop_cosmology_firstorder_equations, this, pert(ithread), tau, pert(ithread)%y, this%saved_source(m)%tau(itau), coop_cosmology_firstorder_ode_accuracy, pert(ithread)%ode_ind, pert(ithread)%ode_c, pert(ithread)%ode_nvars, pert(ithread)%ode_w)
          enddo
       enddo
    enddo
    !$omp end parallel do
    
  end subroutine coop_cosmology_firstorder_compute_source


  subroutine coop_dverk_firstorder(n, fcn, cosmology, pert, x, y, xend, tol, ind, c, nw, w)
    implicit COOP_REAL (a-h,o-z)
    implicit COOP_INT (i-n)
    class(coop_cosmology) cosmology
    class(coop_standard_o1pert) pert
#define DVERK_ARGUMENTS ,cosmology,pert
#include "dverk.h"    
#undef DVERK_ARGUMENTS
  end subroutine coop_dverk_firstorder



end module coop_firstorder_mod
