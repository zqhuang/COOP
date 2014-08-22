module coop_firstorder_mod
  use coop_wrapper_background
  use coop_recfast_mod
  implicit none
#include "constants.h"


  COOP_REAL, parameter :: coop_power_lnk_min = log(3.d-5) 
  COOP_REAL, parameter :: coop_power_lnk_max = log(1.2d0) 
  COOP_REAL, parameter :: coop_visibility_amin = 3.d-4
  COOP_REAL, dimension(0:2), parameter::coop_source_times_weight = (/ 0.3d0, 0.2d0, 0.1d0 /)
  COOP_INT, dimension(0:2), parameter::coop_source_times_n = (/ 1000, 850, 800 /)

  type coop_cosmology_firstorder_source
     logical::used = .false.
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
     COOP_REAL dkappadtau_coef, ReionFrac
     type(coop_function)::Ps, Pt, Xe, ekappa, vis
     type(coop_cosmology_firstorder_source),dimension(0:2)::perturbation
   contains
     procedure:: set_source_times => coop_cosmology_firstorder_set_source_times
     procedure:: set_source_k => coop_cosmology_firstorder_set_source_k
     procedure:: set_power => coop_cosmology_firstorder_set_power
     procedure:: set_xe => coop_cosmology_firstorder_set_xe
     procedure:: set_zre_from_optre => coop_cosmology_firstorder_set_zre_from_optre
     procedure:: set_optre_from_zre => coop_cosmology_firstorder_set_optre_from_zre
     procedure:: xeofa => coop_cosmology_firstorder_xeofa
     procedure:: dkappada => coop_cosmology_firstorder_dkappada
     procedure:: dkappadtau => coop_cosmology_firstorder_dkappadtau
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
    this%used = .false.
  end subroutine coop_cosmology_firstorder_source_free

  subroutine coop_cosmology_firstorder_free(this)
    class(coop_cosmology_firstorder)::this
    COOP_INT m, i
    do m=0,2
       call this%perturbation(m)%free()
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
!!subroutine power(kMpc, ps, pt, cosmology, args)
!!input kMpc and args, output ps and pt    
    integer,parameter::n = 4096
    class(coop_cosmology_firstorder)::this
    external power
    type(coop_arguments) args
    COOP_REAL k(n), ps(n), pt(n)
    COOP_INT i
    call coop_set_uniform(n, k, coop_power_lnk_min, coop_power_lnk_max)
    k = exp(k)
    do i=1, n
       call power(k(i), ps(i), pt(i), this, args)
    end do
    call this%ps%init(n = n, xmin = k(1), xmax = k(n), f = ps, xlog = .true., ylog = .true., fleft = ps(1), fright = ps(n), check_boundary = .false.)
    call this%pt%init(n = n, xmin = k(1), xmax = k(n), f = pt, xlog = .true., ylog = .true., fleft = pt(1), fright = pt(n), check_boundary = .false.)

  end subroutine coop_cosmology_firstorder_set_power

  subroutine coop_cosmology_firstorder_set_xe(this)
    class(coop_cosmology_firstorder)::this
    integer index_baryon, i, m
    integer, parameter::n = 4096
    COOP_REAL ekappa(n), a(n), dkappadinva(n), dkappadtau(n), vis(n)
    index_baryon = this%index_of("Baryon")
    if(index_baryon .eq. 0)then
       stop "Set_Xe: for a model without baryon, you cannot set up Xe"
    endif
    this%dkappadtau_coef = this%species(index_baryon)%Omega * this%h() * coop_SI_sigma_thomson * (coop_SI_rhocritbyh2/coop_SI_c**2) * coop_SI_hbyH0 * coop_SI_c/ coop_SI_m_H * (1.d0 - this%YHe())
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
    do m=0, 2
       this%perturbation(m)%m = m
       call this%set_source_times(this%perturbation(m), coop_source_times_n(m), coop_source_times_weight(m))
    enddo
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

  function coop_cosmology_firstorder_dkappada(this, a) result(dkappada)
    class(coop_cosmology_firstorder)::this    
    COOP_REAL a, dkappada
    dkappada = this%dkappadtau(a) / this%Hasq(a)
  end function coop_cosmology_firstorder_dkappada


  subroutine coop_cosmology_firstorder_set_source_times(this, source, n, weight)
    COOP_REAL, parameter::incr = 1.05d0
    COOP_INT n, i, j
    COOP_REAL  top(n), dtop, amax, amin, amid, topmax, topmin, topmid, weight
    class(coop_cosmology_firstorder_source)::source
    class(coop_cosmology_firstorder)::this
    if(n.lt. 10) stop "n is too small for set_source_times"
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
  end subroutine coop_cosmology_firstorder_set_source_times
  
  subroutine coop_cosmology_firstorder_set_source_k(this, source, n, weight)
    class(coop_cosmology_firstorder)::this
    class(coop_cosmology_firstorder_source)::source
    COOP_INT n
    COOP_REAL weight
    source%nk = n
    if(allocated(source%k))then
       if(size(source%k).ne.n)then
          deallocate(source%k, source%dk)
          allocate(source%k(n), source%dk(n))
       endif
    else
       allocate(source%k(n), source%dk(n))
    endif
  end subroutine coop_cosmology_firstorder_set_source_k

  


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



end module coop_firstorder_mod
