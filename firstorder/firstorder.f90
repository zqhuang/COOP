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
  COOP_REAL, parameter :: coop_visibility_amin = 1.8d-4
  COOP_REAL, parameter :: coop_initial_condition_epsilon = 1.d-5
  COOP_REAL, parameter :: coop_cosmology_firstorder_ode_accuracy = 1.d-7
  COOP_REAL, parameter :: coop_cosmology_firstorder_tc_cutoff = 0.06d0

  COOP_REAL, dimension(0:2), parameter::coop_source_tau_weight = (/ 0.15d0, 0.15d0, 0.1d0 /)
  COOP_INT, dimension(0:2), parameter::coop_source_tau_n = (/ 1000, 850, 800 /)

  COOP_REAL, dimension(0:2), parameter::coop_source_k_weight = (/ 0.2d0, 0.15d0, 0.1d0 /)

  COOP_INT, dimension(0:2), parameter::coop_source_k_n = (/ 160, 120, 80 /)

  COOP_INT, dimension(0:2), parameter::coop_num_sources = (/ 4,  2,  2 /)

  COOP_REAL, parameter::coop_source_k_index = 0.55d0


  type coop_cosmology_firstorder_source
     COOP_INT::m = 0
     COOP_INT::ntau = 0
     COOP_INT::nk = 0
     COOP_INT::index_tc_max = 1
     COOP_INT, dimension(:),allocatable::index_tc_off
     COOP_REAL, dimension(:),allocatable::k, tau, a, tauc, dk, dtau, chi  !!tau is conformal time, chi is comoving distance; in flat case, chi + tau = tau_0
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
     COOP_REAL::dkappadtau_coef, ReionFrac, Omega_m, Omega_r, Omega_b, Omega_c, Omega_nu, Omega_g, a_eq, tau_eq, mnu_by_Tnu
     type(coop_function)::Ps, Pt, Xe, ekappa, vis, Tb
     type(coop_cosmology_firstorder_source),dimension(0:2)::saved_source
     COOP_INT::index_baryon, index_cdm, index_radiation, index_nu, index_massiveNu, index_de
     COOP_REAL, dimension(0:coop_pert_default_lmax, 0:coop_pert_default_mmax, 0:coop_pert_default_smax)::klms, klms_by_2lm1, klms_by_2lp1
     logical::klms_done = .false.
   contains
     procedure:: setup_all => coop_cosmology_firstorder_setup_all
     procedure:: set_klms => coop_cosmology_firstorder_set_klms
     procedure:: set_source_tau => coop_cosmology_firstorder_set_source_tau
     procedure:: set_source_k => coop_cosmology_firstorder_set_source_k
     procedure:: set_power => coop_cosmology_firstorder_set_power
     procedure:: set_xe => coop_cosmology_firstorder_set_xe
     procedure:: set_zre_from_optre => coop_cosmology_firstorder_set_zre_from_optre
     procedure:: set_optre_from_zre => coop_cosmology_firstorder_set_optre_from_zre
     procedure::set_tight_coupling => coop_cosmology_firstorder_set_tight_coupling
     procedure::set_metric => coop_cosmology_firstorder_set_metric

     procedure::set_late_approx => coop_cosmology_firstorder_set_late_approx
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
     procedure:: free => coop_cosmology_firstorder_free
     procedure::pert2source => coop_cosmology_firstorder_pert2source
     procedure::init_source => coop_cosmology_firstorder_init_source
     procedure::compute_source =>  coop_cosmology_firstorder_compute_source
  end type coop_cosmology_firstorder


contains

  subroutine coop_cosmology_firstorder_setup_all(this)
    class(coop_cosmology_firstorder)::this
  end subroutine coop_cosmology_firstorder_setup_all

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
    if(allocated(this%tau))deallocate(this%tau, this%chi, this%dtau, this%a, this%tauc)
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
    call this%pt%init(n = n, xmin = k(1), xmax = k(n), f = pt, xlog = .true., ylog = .true., fleft = pt(1), fright = pt(n), check_boundary = .false.)

  end subroutine coop_cosmology_firstorder_set_power

  subroutine coop_cosmology_firstorder_set_xe(this)
    class(coop_cosmology_firstorder)::this
    integer i, m
    integer, parameter::n = 4096
    COOP_REAL ekappa(n), a(n), dkappadinva(n), dkappadtau(n), vis(n), acal
    call this%set_klms()
    this%index_baryon = this%index_of("Baryon", .true.)
    this%index_radiation = this%index_of("Radiation", .true.)
    this%index_cdm = this%index_of("CDM", .true.)
    this%index_Nu = this%index_of("Massless Neutrinos", .true.)
    this%index_massiveNu = this%index_of("Massive Neutrinos")
    this%index_de = this%index_of("Dark Energy", .true.)
    
    acal  = coop_visibility_amin

    !!asymptotic Omega_b and Omega_c
    this%Omega_b = O0_BARYON(this)%Omega *  O0_BARYON(this)%rhoa3_ratio(acal)
    this%Omega_c = O0_CDM(this)%Omega * O0_CDM(this)%rhoa3_ratio(acal)

    this%Omega_m = this%Omega_b + this%Omega_c
    this%Omega_g = O0_RADIATION(this)%Omega * O0_RADIATION(this)%rhoa4_ratio(acal)
    if(this%index_massiveNu .ne. 0)then
       call coop_fermion_get_lnam(log(O0_MASSIVENU(this)%Omega/O0_MASSIVENU(this)%Omega_massless), this%mnu_by_Tnu)
       this%mnu_by_Tnu = exp(this%mnu_by_Tnu)
       this%Omega_nu = O0_NU(this)%Omega * O0_NU(this)%rhoa4_ratio(acal) + O0_MASSIVENU(this)%Omega * O0_MASSIVENU(this)%rhoa4_ratio(acal)
    else
       this%Omega_nu = O0_NU(this)%Omega * O0_NU(this)%rhoa4_ratio(acal)
    endif

    this%Omega_r = this%Omega_g + this%Omega_nu 

    this%a_eq = this%Omega_r / this%Omega_m
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
    if(a .lt. coop_recfast_a_initial)then
       xe = this%xe%fleft
    else
       xe = this%xe%eval(a)
    endif
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
    COOP_INT::nbuffer = 10
    COOP_INT n, i, j
    COOP_REAL  top(n), dtop, amax, amin, amid, topmax, topmin, topmid, weight, taucut, acut
    class(coop_cosmology_firstorder_source)::source
    class(coop_cosmology_firstorder)::this
    if(n.lt. 30) stop "n is too small for set_source_tau"
    source%ntau = n
    if(allocated(source%a))then
       if(size(source%a).ne.n)then
          deallocate(source%a, source%tau, source%chi, source%dtau, source%tauc)
          allocate(source%a(n), source%tau(n), source%dtau(n), source%chi(n), source%tauc(n))
       endif
    else
       allocate(source%a(n), source%tau(n), source%dtau(n), source%chi(n),  source%tauc(n))
    endif
    dtop = this%tau0 / (n - nbuffer)
    do i= nbuffer+1, n
       top(i) = dtop*(i-nbuffer-0.5d0)
    enddo
    do i=nbuffer, 1, -1
       top(i) = top(i+1)*0.001d0
    enddo
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
       do while((amax-amin)/amin .gt. 1.d-7)
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
       if(i.gt.nbuffer+1)source%dtau(i) = dtop/((weight*this%dkappadtau(source%a(i))*source%tau(i)+1.d0)*this%ekappaofa(source%a(i))**weight)
       source%tauc(i) = this%taucofa(source%a(i))
    enddo
    source%dtau(1) = (source%tau(1)+source%tau(2))/2.d0
    do i=2, nbuffer+1
       source%dtau(i) = (source%tau(i+1) - source%tau(i-1))/2.d0
    enddo
    source%chi = this%tau0 - source%tau
    source%index_tc_max = 1
    do while( source%index_tc_max .lt. source%ntau - 1 .and. source%tauc(source%index_tc_max+1) .le. 0.03 * source%tau(source%index_tc_max+1) )
       source%index_tc_max = source%index_tc_max+1
    enddo
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
          deallocate(source%k, source%dk, source%index_tc_off)
          allocate(source%k(n), source%dk(n), source%index_tc_off(n))
       endif
    else
       allocate(source%k(n), source%dk(n), source%index_tc_off(n))
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
    if(source%ntau .gt. 0 .and. allocated(source%tauc))then
       source%index_tc_off = 1
       do i = 1, n
          do while(source%index_tc_off(i) .lt. source%index_tc_max .and. source%tauc(source%index_tc_off(i)+1)*source%k(i) .le. coop_cosmology_firstorder_tc_cutoff)
             source%index_tc_off(i)= source%index_tc_off(i)+1
          enddo
       enddo
    else
       stop "You need to call set_source_tau before calling set_source_k"
    endif
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


  subroutine coop_cosmology_firstorder_equations(n, tau, var, varp, cosmology, pert, pertp)
    COOP_INT n
    class(coop_cosmology_firstorder)::cosmology
    class(coop_pert_object)::pert, pertp
    COOP_INT i, l, iq
    COOP_REAL tau, var(n), varp(n), a
    COOP_REAL source(0:coop_pert_default_lmax), psi, pi
    a = cosmology%aoftau(tau)
  end subroutine coop_cosmology_firstorder_equations

  subroutine coop_cosmology_firstorder_set_tight_coupling(this, pert, pertp, a)
    class(coop_cosmology_firstorder)::this
    class(coop_pert_object)::pert, pertp
    COOP_REAL a
  end subroutine coop_cosmology_firstorder_set_tight_coupling



  subroutine coop_cosmology_firstorder_set_metric(this, pert, pertp, a)
    class(coop_cosmology_firstorder)::this
    class(coop_pert_object)::pert, pertp
    COOP_REAL a
!!$    select case(pert%m)
!!$    case(0)
!!$       if(this%index_massivenu .ne. 0)then
!!$          if(pert%has_massivenu)then
!!$             pert%O1_PI =  (12.d0/5.d0) *  (O0_RADIATION(this)%pa2(a) * pert%O1T(2) + O0_NU(this)%pa2(a) * pert%O1NU(2))
!!$          else
!!$             pert%O1_PI =  (12.d0/5.d0) *  (O0_RADIATION(this)%pa2(a) * pert%O1T(2) + (O0_NU(this)%pa2(a) + O0_MASSIVENU(this)%pa2(a)) * pert%O1NU(2))
!!$          endif
!!$       else
!!$          pert%O1_PI =  (12.d0/5.d0) *  (O0_RADIATION(this)%pa2(a) * pert%O1T(2) + O0_NU(this)%pa2(a) * pert%O1NU(2))
!!$       endif
!!$       pert%O1_PHI = pert%O1_PSI - pert%O1_PI/pert%k**2
!!$    case(1)
!!$       call coop_tbw("set_metric")
!!$    case(2)
!!$       call coop_tbw("set_metric")
!!$    end select
  end subroutine coop_cosmology_firstorder_set_metric


  subroutine coop_cosmology_firstorder_set_late_approx(this, pert, pertp, a)
    class(coop_cosmology_firstorder)::this
    class(coop_pert_object)::pert, pertp
    COOP_REAL a
    call coop_tbw("set_late_approx")
  end subroutine coop_cosmology_firstorder_set_late_approx


  subroutine coop_cosmology_firstorder_set_initial_conditions(this, pert, k, tau)
    class(coop_cosmology_firstorder)::this
    class(coop_pert_object)::pert
    COOP_REAL tau, k, Rnu, normC
!!$    pert%k = k
!!$    select case(trim(pert%initial_conditions))
!!$    case("adiabatic")
!!$       Rnu = this%Omega_nu / this%Omega_r
!!$       normC = coop_primordial_zeta_norm/(1.5d0 + 0.4d0*Rnu)
!!$       select case(pert%m)
!!$       case(0)
!!$          !!gravitational potential
!!$          pert%O1_PHI = normC
!!$          pert%O1_PHIDOT = 0.d0
!!$          pert%O1_PSI = normC * (1.d0 + 0.4d0 * Rnu)
!!$          pert%O1_PSIDOT = 0.d0
!!$          !!density perturbations (the 0-th multipole =  density fluctuation /[3(1+w)])
!!$          pert%O1BARYON(0) = - normC/2.d0
!!$          pert%O1CDM(0) =   pert%O1BARYON(0)
!!$          pert%O1T(0) =   pert%O1BARYON(0)
!!$          pert%O1NU(0) =   pert%O1BARYON(0)
!!$          if(pert%has_massivenu) pert%O1MASSIVENU(0,:) =   pert%O1BARYON(0)
!!$          !!velocity fluctuations
!!$          pert%O1BARYON(1) = 0.5d0*k*tau*normC
!!$          pert%O1CDM(1) = pert%O1BARYON(1)
!!$          pert%O1T(1) = pert%O1BARYON(1)
!!$          pert%O1NU(1) = pert%O1BARYON(1)
!!$          if(pert%has_massivenu) pert%O1MASSIVENU(1,:) =   pert%O1BARYON(1)
!!$          !!anisotropy stress, photon does not have anisotropy stress because it is coupled to baryon
!!$          pert%O1NU(2) = (1.d0/6.d0)*(k*tau)**2*normC
!!$          pert%PI = (12.d0/5.d0)*pert%O1NU(2)
!!$          !!dark energy
!!$          select case(pert%de_pert_model)
!!$          case (COOP_DE_PERT_NONE)
!!$          case(COOP_DE_PERT_FLUID)
!!$             stop "de_pert_fluid not implemented"
!!$          case(COOP_DE_PERT_PPF)
!!$             stop "de_pert_ppf not implemented"
!!$          case(COOP_DE_PERT_QUINTESSENCE)
!!$             stop "de_pert_quintessence not implemented"
!!$          case(COOP_DE_PERT_COUPLED_QUINTESSENCE)
!!$             stop "de_pert_coupled_quintessence not implemented"
!!$          case default
!!$             call coop_return_error("firstorder_set_initial_conditions", "Unknown de_pert_model", "stop")
!!$          end select
!!$       case(1)
!!$          stop "primordial vector is supposed to vanish in the standard model; modify the code here if you want to implement some non-standard model."
!!$       case(2)
!!$          pert%O1_H = coop_primordial_zeta_norm/coop_sqrt6
!!$       end select
!!$    case default
!!$       stop "unknown initial conditions"
!!$    end select
  end subroutine coop_cosmology_firstorder_set_initial_conditions


  subroutine coop_cosmology_firstorder_init_source(this, m)
    class(coop_cosmology_firstorder)::this
    COOP_INT :: m
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
  end subroutine coop_cosmology_firstorder_init_source


  subroutine coop_cosmology_firstorder_compute_source(this, m)
    class(coop_cosmology_firstorder)::this
    COOP_INT, parameter::n_threads = 8
    type(coop_pert_object),dimension(n_threads):: pert, pertp
    COOP_INT m, ik, itau, ithread, iq
    COOP_REAL tau, a
    call this%init_source(m)
!!$    if(this%index_de .ne. 0)then
!!$       if(this%species(this%index_de)%genre .eq. COOP_SPECIES_LAMBDA)then
!!$          pert%de_pert_model = COOP_DE_PERT_NONE
!!$       else
!!$          pert%de_pert_model = COOP_DE_PERT_FLUID
!!$       endif
!!$    else
!!$       pert%de_pert_model = COOP_DE_PERT_NONE
!!$    endif
!!$    do ithread = 1, n_threads
!!$       call pert(ithread)%init(m)
!!$       call pert(ithread)%set_ode_index()
!!$       pertp(ithread) = pert(ithread)
!!$    enddo
!!$
!!$    !$omp parallel do private(ik, itau, ithread, tau, iq)
!!$    do ithread = 1, n_threads
!!$       do ik = ithread, this%saved_source(m)%nk, n_threads
!!$          tau = min(coop_initial_condition_epsilon/this%saved_source(m)%k(ik), this%tau_eq*coop_initial_condition_epsilon)
!!$          pert(ithread)%tight_coupling = .true.
!!$          call pert(ithread)%set_ode_index()
!!$          call this%set_initial_conditions( pert(ithread), this%saved_source(m)%k(ik), tau)
!!$          call pert(ithread)%write_var()
!!$          call coop_cosmology_firstorder_equations(pert(ithread)%ode_nvars, tau, pert(ithread)%y, pertp(ithread)%y, this, pert(ithread), pertp(ithread))
!!$          do itau = 1, this%saved_source(m)%ntau
!!$             if(itau .eq. this%saved_source(m)%index_tc_off(ik))then
!!$                call pert(ithread)%read_var() 
!!$                pert(ithread)%tight_coupling = .false.
!!$                call pert(ithread)%set_ode_index()
!!$                call pert(ithread)%write_var()    
!!$             endif
!!$             call coop_dverk_firstorder(pert(ithread)%ode_nvars, coop_cosmology_firstorder_equations, this, pert(ithread), pertp(ithread), tau, pert(ithread)%y, this%saved_source(m)%tau(itau), coop_cosmology_firstorder_ode_accuracy, pert(ithread)%ode_ind, pert(ithread)%ode_c, pert(ithread)%ode_nvars, pert(ithread)%ode_w)
!!$             call this%pert2source(pert(ithread), this%saved_source(m)%source(itau, ik, :))
!!$          enddo
!!$       enddo
!!$    enddo
!!$    !$omp end parallel do
    
  end subroutine coop_cosmology_firstorder_compute_source

  subroutine coop_cosmology_firstorder_pert2source(this, pert, source)
    class(coop_cosmology_firstorder)::this
    class(coop_pert_object)::pert
    COOP_REAL source(:)
    stop "TBW"
  end subroutine coop_cosmology_firstorder_pert2source


  subroutine coop_dverk_firstorder(n, fcn, cosmology, pert, pertp, x, y, xend, tol, ind, c, nw, w)
    implicit COOP_REAL (a-h,o-z)
    implicit COOP_INT (i-n)
    class(coop_cosmology) cosmology
    class(coop_pert_object) pert, pertp
#define DVERK_ARGUMENTS ,cosmology,pert,pertp
#include "dverk.h"    
#undef DVERK_ARGUMENTS
  end subroutine coop_dverk_firstorder



end module coop_firstorder_mod
