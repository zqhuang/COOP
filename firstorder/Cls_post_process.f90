module coop_cls_postprocess_mod
  use coop_wrapper_background
  use coop_pertobj_mod
  use coop_cl_indices_mod
  use coop_firstorder_mod
#ifdef HAS_HEALPIX
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
  use udgrade_nr
#endif
  implicit none
#include "constants.h"

  private



  COOP_INT::coop_zeta_nr
  COOP_REAL,dimension(:),allocatable::coop_zeta_r, coop_zeta_dr
  
  type coop_zeta_shell
     COOP_INT::lmax
     COOP_REAL::r
     COOP_SINGLE,dimension(:),allocatable::alm_real
   contains
     procedure::free=>coop_zeta_shell_free
     procedure::map_project => coop_zeta_shell_map_project
     procedure::set_lmax => coop_zeta_shell_set_lmax
  end type coop_zeta_shell

     

  public::coop_generate_3dzeta, coop_get_zeta_shells_cov,  coop_get_zeta_trans_l, coop_zeta_r, coop_zeta_nr,  coop_set_default_zeta_r, coop_zeta_dr, coop_zeta_shell

contains

  subroutine coop_zeta_shell_free(shell)
    class(coop_zeta_shell)::shell
    if(allocated(shell%alm_real))deallocate(shell%alm_real)
    shell%lmax = -1
  end subroutine coop_zeta_shell_free

  subroutine coop_zeta_shell_map_project(shell, fnl, lmax, alm_total, weight)
    class(coop_zeta_shell)::shell
    external fnl
    COOP_INT lmax, l, m
    COOP_REAL weight(0:lmax)
    complex alm_total(0:lmax, 0:lmax)
#ifdef HAS_HEALPIX
    real, dimension(:),allocatable::map
    integer npix, nside, i
    complex, dimension(:,:,:),allocatable::alms
    nside = 32
    do while(nside*2 .lt. shell%lmax)
       nside = nside*2
    enddo
    npix = 12*nside**2
    allocate(alms(0:shell%lmax, 0:shell%lmax, 1), map(0:npix-1))
    alms = 0.
    do l=2, shell%lmax
       alms(l, 0, 1) = cmplx(shell%alm_real(l*l+1), 0.)
       do m=1, l
          alms(l, m, 1) = cmplx(shell%alm_real(l*l+m*2)/coop_sqrt2, shell%alm_real(l*l+m*2+1)/coop_sqrt2 )
       enddo
    enddo
    call alm2map(nside, shell%lmax, shell%lmax, alms, map)
    do i = 0, npix-1
       call fnl(map(i))
    enddo
    call map2alm(nside, shell%lmax, shell%lmax, map, alms)
    do l=2, shell%lmax
       alm_total(l, 0:l) =   alm_total(l, 0:l) + alms(l, 0:l, 1)*weight(l)
    enddo

    deallocate(alms, map)
#else
    call coop_return_error("zeta_shell_map_project", "Cannot find Healpix", "stop")
#endif
    
  end subroutine coop_zeta_shell_map_project


  subroutine coop_zeta_shell_set_lmax(shell, lmax)
    class(coop_zeta_shell)::shell
    COOP_INT lmax
    shell%lmax = lmax
    if(allocated(shell%alm_real))then
       if(  size(shell%alm_real) .ne. (lmax+1)**2 )then
          deallocate(shell%alm_real)
          allocate(shell%alm_real( (lmax+1)**2 ))
       endif
    else
       allocate(shell%alm_real( (lmax+1)**2 ))
    endif
    shell%alm_real = 0.
  end subroutine coop_zeta_shell_set_lmax

  subroutine coop_set_default_zeta_r(this, source)
    class(coop_cosmology_firstorder)::this
    type(coop_cosmology_firstorder_source)::source
    COOP_INT,parameter::N = 8192
    COOP_REAL::chi(N), chimax, step
    COOP_INT::i, nchi
    nchi = source%ntau
    if(nchi .gt. n) call coop_return_error("set_default_zeta_r", "nr overflow", "stop")
    chi(1:nchi) = source%chi(nchi:1:-1)
    chimax = this%tau0*1.04d0
    do while(chi(nchi) .lt. chimax)
       nchi = nchi + 1
       if(nchi .gt. n) call coop_return_error("set_default_zeta_r", "nr overflow", "stop")
       chi(nchi) = chi(nchi-1)+ min((chi(nchi-1)-chi(nchi-2))*1.03d0, 1.d-3)
    enddo
    coop_zeta_nr = nchi
    if(allocated(coop_zeta_dr))deallocate(coop_zeta_r, coop_zeta_dr)
    allocate(coop_zeta_r(nchi), coop_zeta_dr(nchi))
    coop_zeta_r = chi(nchi:1:-1)
    do i = 2, nchi-1
       coop_zeta_dr(i) = (coop_zeta_r(i-1)-coop_zeta_r(i+1))/2.d0
    enddo
    coop_zeta_dr(1) = coop_zeta_r(1) - coop_zeta_r(2)
    coop_zeta_dr(nchi) = coop_zeta_r(nchi-1) - coop_zeta_r(nchi)
   
  end subroutine coop_set_default_zeta_r


  subroutine coop_generate_3dzeta_partial(cosmology, lmax, nr, r, alm_real, l, ir_start, ir_end)
    class(coop_cosmology_firstorder)::cosmology
    COOP_INT::lmax, nr, l, ir_start, ir_end, nused, nm, i
    COOP_REAL::alm_real(-lmax:lmax, nr), r(nr)
    COOP_REAL::cov(ir_end - ir_start + 1, ir_end - ir_start + 1)
    if(l.le.1 .or. l.gt.lmax)then
       return
    endif
    nused = ir_end - ir_start + 1
    call coop_get_zeta_shells_cov(cosmology, l, nused, r(ir_start:ir_end), cov)
    call coop_matsym_sqrt(Cov)
    nm = 2*l+1
    do i= ir_start, ir_end
       call coop_random_get_Gaussian_vector(nm, alm_real(-l:l, i))
    end do
    alm_real(-l:l,ir_start:ir_end) = matmul(alm_real(-l:l,ir_start:ir_end), Cov)
  end subroutine coop_generate_3dzeta_partial


  subroutine coop_get_zeta_shells_cov(cosmology, l, nr, r, cov)
    class(coop_cosmology_firstorder)::cosmology
    COOP_INT l, nr, i, j
    COOP_REAL r(nr), cov(nr, nr)
    !$omp parallel do private(i, j)
    do i=1, nr
       do j=1, i-1
          cov(i, j) = cosmology%Clzetazeta(l, r(i), r(j))
          cov(j, i) = cov(i, j)
       enddo
       cov(i, i) = cosmology%Clzetazeta(l, r(i))+1.d-20  !!add small tiny number for stability
    enddo
    !$omp end parallel do

  end subroutine coop_get_zeta_shells_cov


  subroutine coop_generate_3Dzeta(cosmology, ns, shells)
    class(coop_cosmology_firstorder)::cosmology
    COOP_INT:: ns
    type(coop_zeta_shell)::shells(ns)
    COOP_REAL,dimension(:,:),allocatable::zlms
    COOP_REAL::r(ns)
    COOP_INT i, iend, lmm, istart, l
    do i=1, ns
       r(i) = shells(i)%r
    enddo
    lmm = shells(1)%lmax
    do i = 2, ns
       if(shells(i)%lmax .gt. shells(i-1)%lmax)then
          stop "lmax must be decreasing order"
       endif
    enddo
    allocate(zlms(-lmm:lmm, ns))
    istart = 1
    iend = ns
    do l=2, lmm
       do while(shells(iend)%lmax .lt. l)
          iend = iend - 1
       enddo
       call coop_generate_3Dzeta_partial(cosmology, lmm, ns, r, zlms, l, istart, iend)
       do i = istart, iend
          shells(i)%alm_real(l**2+1:(l+1)**2) = zlms(-l:l, i)
       enddo
    enddo
    deallocate(zlms)
  end subroutine coop_generate_3Dzeta

!!$mapping 3D zeta to CMB alms
!!The formula is
!!a_{lm} = 2/\pi  \int_0^\infty  j_l(kr) \zeta_{lm}(r)  \Delta_l(k) k^2 dk r^2 dr
!!where \Delta_l(k) = \int S(k, \chi) j_l(k \chi) d\chi is the transfer function computed in get_transfer function (see firstorder.f90)
  subroutine coop_get_zeta_trans_l(source, l, nr, r, trans)
    COOP_INT,parameter::n_threads = 16
    type(coop_cosmology_firstorder_source)::source
    COOP_INT l, ik, ir, ichi, idense, isrc, i, nr, lasti, ithread
    COOP_REAL::r(nr), trans(nr, source%nsrc)
    COOP_REAL, dimension(:,:,:),allocatable::trans_chi_r
    logical,dimension(:),allocatable::computed
    logical,dimension(:,:),allocatable::compute_chi
    COOP_REAL, dimension(:,:,:),allocatable::ampchi, phasechi, ampr, phaser
    COOP_REAL,dimension(:,:),allocatable::kwindow
    COOP_REAL chi1, dchimin, kmax
    logical:: chi_r_overlap
    COOP_INT::nbuffer

    allocate(trans_chi_r(source%nsrc, source%ntau, nr), ampchi(coop_k_dense_fac, source%nk, source%ntau), phasechi(coop_k_dense_fac, source%nk, source%ntau), ampr(coop_k_dense_fac, source%nk, nr), phaser(coop_k_dense_fac, source%nk, nr), computed(nr), kwindow(coop_k_dense_fac, source%nk), compute_chi(source%ntau, n_threads))

    kmax = 3.d0*max(l+0.5d0, 600.d0)/source%chi(1)
    !$omp parallel do private(ik, idense)
    do ik = 1, source%nk
       do idense = 1, coop_k_dense_fac
          kwindow(idense, ik) = (1.d0-tanh((source%k_dense(idense, ik)/kmax-0.5d0)*20.d0))/2.d0
       enddo
    enddo
    !$omp end parallel do
    ir = 1
    chi1 = source%chi(1)+1.d-9
    chi_r_overlap = .true.
    do while(r(ir).gt. chi1 )
       ir = ir + 1
       if(ir.eq.nr)then
          chi_r_overlap = .false.
          exit
       endif
    enddo
    if(chi_r_overlap)then
       nbuffer = ir-1
       if( source%ntau  .ne. nr-ir+1)then
          chi_r_overlap = .false.
       else
          if(any(abs(source%chi - r(ir:nr)).gt. 1.d-9))then
             chi_r_overlap = .false.
          endif
       endif
    endif
    call coop_jl_setup_amp_phase(l)
    !$omp parallel do private(ichi, ik ,idense)
    do ichi = 1, source%ntau
       do ik = 2, source%nk
          do idense = 1, coop_k_dense_fac
             call coop_jl_get_amp_phase(l, source%k_dense(idense, ik)*source%chi(ichi), ampchi(idense, ik, ichi), phasechi(idense, ik, ichi))
          enddo
       enddo
    enddo
    !$omp end parallel do
    if(chi_r_overlap)then
       !$omp parallel do private(ir, ik ,idense)
       do ir = 1, nbuffer
          do ik = 2, source%nk
             do idense = 1, coop_k_dense_fac
                call coop_jl_get_amp_phase(l, source%k_dense(idense, ik)*r(ir), ampr(idense, ik, ir), phaser(idense, ik, ir))
             enddo
          enddo
       enddo
       !$omp end parallel do
       ampr(:,:,nbuffer+1:nr) = ampchi(:, :, source%ntau:1:-1)
       phaser(:,:,nbuffer+1:nr) = phasechi(:, :, source%ntau:1:-1)
    else
       !$omp parallel do private(ir, ik ,idense)
       do ir = 1, nr
          do ik = 2, source%nk
             do idense = 1, coop_k_dense_fac
                call coop_jl_get_amp_phase(l, source%k_dense(idense, ik)*r(ir), ampr(idense, ik, ir), phaser(idense, ik, ir))
             enddo
          enddo
       enddo
       !$omp end parallel do
    endif
    computed = .false.
    computed(1) = .true.
    lasti = 1
    ir = 1
    do while(ir.lt.nr)
       if(r(lasti) - r(ir) .lt. 0.002)then
          ir = ir + 1
       else
          computed(ir) = .true.
          lasti = ir
       endif
    enddo
    computed(nr) = .true.
    trans_chi_r = 0.d0

    !$omp parallel do private(ir, ichi, ithread, lasti)
    do ithread = 1, n_threads
       do ir = ithread, nr, n_threads
          if(computed(ir))then
             compute_chi(1, ithread) = .true.
             lasti = 1
             do ichi=2, source%ntau
                if(source%chi(lasti)- source%chi(ichi) .lt. 0.8d-3)then
                   compute_chi(ichi, ithread) = .false.
                else
                   compute_chi(ichi, ithread) = .true.
                   lasti = ichi
                endif
             enddo
             compute_chi(source%ntau, ithread) = .true.
             do ichi=1, source%ntau
                if(abs(source%chi(ichi) - r(nr)).lt. 0.05d0)then
                   compute_chi(ichi, ithread) = .true.
                endif
             enddo
             do ichi = 2, source%ntau-1
                if(compute_chi(ichi, ithread))call coop_get_zeta_trans_l_step(source, l, nr, r, ichi, ir, trans_chi_r(:, ichi, ir), ampchi, phasechi, ampr, phaser, kwindow)
             enddo
             call coop_spline_fill(source%nsrc, source%ntau, source%chi, trans_chi_r(:, :, ir), compute_chi(:, ithread), .false., .false.)
          endif
       enddo
    enddo
    !$omp end parallel do

    call coop_spline_fill(source%nsrc, source%ntau, nr, r, trans_chi_r, computed, .false., .false.)

    !$omp parallel do private(isrc, ir)
    do isrc = 1, source%nsrc
       do ir = 1, nr
          trans(ir, isrc) = (2.d0/coop_pi)*sum(trans_chi_r(isrc, :, ir)*source%dtau)*r(ir)**2
       enddo
    enddo
    !$omp end parallel do

    deallocate(trans_chi_r,ampchi, phasechi, ampr, phaser, computed, compute_chi, kwindow)
    call coop_jl_destroy_amp_phase(l)
  end subroutine coop_get_zeta_trans_l


  !!compute \int_0^\infty S(k, \chi) j_l(k \chi) j_l (k r) k^2 dk
  !!input trans must be initialized to zero 
  subroutine coop_get_zeta_trans_l_step(source,  l, nr, r, ichi, ir, trans, ampchi, phasechi, ampr, phaser, kwindow)
    COOP_REAL,parameter::dphase = 0.6d0
    type(coop_cosmology_firstorder_source)::source
    COOP_INT::l, ichi, ir, nr
    COOP_REAL::r(nr)
    COOP_REAL,intent(INOUT)::trans(source%nsrc)
    COOP_REAL::widthm, widthp
    COOP_REAL::ampchi(coop_k_dense_fac, source%nk, source%ntau), phasechi(coop_k_dense_fac, source%nk, source%ntau), ampr(coop_k_dense_fac, source%nk, nr), phaser(coop_k_dense_fac, source%nk, nr), kwindow(coop_k_dense_fac, source%nk)
    COOP_INT ik, idense
    COOP_REAL  kmin, xmin, phasediff, phasesum, last_phasediff, last_phasesum, sumfac, difffac
    logical::dosum
!!$    if(source%a(ichi) .gt. 0.3d0)then !! low redshift stuff, use approx S(k, chi)  = S(k->0, chi)
!!$       if(abs(source%chi(ichi) - r(ir)).lt. 1.d-6)then
!!$          kmin = (l+0.5d0)/r(ir)
!!$          ik = 1
!!$          do while(ik .lt. source%nk .and. source%k(ik).lt. kmin)
!!$             ik = ik + 1
!!$          enddo
!!$          if(ik.ge.source%nk)return
!!$          trans = coop_pio2 * source%s(:, ik, ichi) / source%chi(ichi)**2/source%dtau(ichi)
!!$       endif
!!$       return
!!$    endif
    call coop_jl_startpoint(l, xmin)
    kmin = xmin/ min(source%chi(ichi), r(ir))
    last_phasediff = abs(phasechi(1, 2, ichi)-phaser(1, 2, ir))
    last_phasesum = phasechi(1, 2, ichi)+phaser(1, 2, ir)
    dosum = .true.
    do ik = 2, source%nk
       if(kwindow(1, ik).lt.1.d-4)return
       if(source%k(ik) .le. kmin )then
          last_phasediff = abs(phasechi(coop_k_dense_fac, ik, ichi)-phaser(coop_k_dense_fac, ik, ir))
          last_phasesum = phasechi(coop_k_dense_fac, ik, ichi)+phaser(coop_k_dense_fac, ik, ir)
          cycle
       endif
       if(dosum)then
          do idense = 1, coop_k_dense_fac
             phasediff = abs(phasechi(idense, ik, ichi)-phaser(idense, ik, ir))
             phasesum = phasechi(idense, ik, ichi)+phaser(idense, ik, ir)
             if(phasesum - last_phasesum .gt. dphase)then
                dosum = .false.
                sumfac = 0.d0
             else
                sumfac = 1.d0+tanh(20.d0-(40.d0/dphase)*(phasesum - last_phasesum))
             endif
             difffac = 1.d0+tanh(20.d0-(40.d0/dphase)*(phasediff - last_phasediff))
             trans = trans + COOP_INTERP_SOURCE(source, :, idense, ik, ichi) * (ampchi(idense, ik, ichi) * ampr(idense, ik, ir)*( cos(phasediff)*difffac + cos(phasesum)*sumfac )*source%k_dense(idense, ik)*source%ws_dense(idense, ik)/source%ps_dense(idense, ik)/4.d0)
             last_phasediff = phasediff
             last_phasesum = phasesum
          enddo
       else
          do idense = 1, coop_k_dense_fac
             phasediff = abs(phasechi(idense, ik, ichi)-phaser(idense, ik, ir))
             difffac = 1.d0+tanh(20.d0-(40.d0/dphase)*(phasediff - last_phasediff))
             trans = trans + COOP_INTERP_SOURCE(source, :, idense, ik, ichi) * (ampchi(idense, ik, ichi) * ampr(idense, ik, ir)*cos(phasediff)*source%k_dense(idense, ik)*source%ws_dense(idense, ik)/source%ps_dense(idense, ik) * difffac *kwindow(idense, ik)/4.d0)
             if(difffac .lt. 1.d-4)return
             last_phasediff = phasediff
          enddo
       endif
    enddo
  end subroutine coop_get_zeta_trans_l_step

end module coop_cls_postprocess_mod
