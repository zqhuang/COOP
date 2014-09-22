module coop_cls_postprocess_mod
  use coop_wrapper_background
  use coop_pertobj_mod
  use coop_firstorder_mod
  implicit none
#include "constants.h"

  private

  COOP_INT::coop_zeta_nr
  COOP_REAL,dimension(:),allocatable::coop_zeta_r, coop_zeta_dr

  public::coop_generate_3dzeta, coop_load_2Dzeta, coop_get_zeta_shells_cov,  coop_get_zeta_trans_l, coop_zeta_r, coop_zeta_nr,  coop_set_default_zeta_r, coop_zeta_dr

contains

  subroutine coop_set_default_zeta_r(source, nbuffer)
    type(coop_cosmology_firstorder_source)::source
    COOP_INT nbuffer, i
    COOP_REAL dchi
    coop_zeta_nr = source%ntau + nbuffer
    if(allocated(coop_zeta_r))then
       if(size(coop_zeta_r).ne. coop_zeta_nr)then
          deallocate(coop_zeta_r, coop_zeta_dr)
          allocate(coop_zeta_r(coop_zeta_nr), coop_zeta_dr(coop_zeta_nr))
       endif
    else
       allocate(coop_zeta_r(coop_zeta_nr), coop_zeta_dr(coop_zeta_nr))
    endif
    coop_zeta_r(nbuffer+1:coop_zeta_nr) = source%chi
    coop_zeta_dr(nbuffer+1:coop_zeta_nr) = source%dtau
    dchi = source%chi(1)-source%chi(2)
    do i = nbuffer, 1, -1
       if(dchi .lt. 0.05)dchi = dchi*1.06d0
       coop_zeta_r(i) =  coop_zeta_r(i+1)+dchi
    enddo
    do i = nbuffer+1, 2, -1
       coop_zeta_dr(i) = (coop_zeta_r(i-1)-coop_zeta_r(i+1))/2.d0
    enddo
    coop_zeta_dr(1) = coop_zeta_r(1) - coop_zeta_r(2)
  end subroutine coop_set_default_zeta_r

  subroutine coop_generate_3dzeta_l(cosmology, l, lmax, nr, r, zetalm)
    class(coop_cosmology_firstorder)::cosmology
    COOP_INT::l, nr, i, j, nm, lmax
    COOP_REAL::r(nr)
    COOP_REAL::zetalm(-lmax:lmax, nr)
    COOP_REAL::cov(nr, nr)
    if(l.le.1 .or. l.gt.lmax)then
       zetalm = 0.d0
       return
    endif
    call coop_get_zeta_shells_cov(cosmology, l, nr, r, cov)
    call coop_matsym_sqrt(Cov)
    nm = 2*l+1
    do i=1, nr
       call coop_random_get_Gaussian_vector(nm, zetalm(-l:l, i))
    end do
    zetalm(-l:l,:) = matmul(zetalm(-l:l,:), Cov)
  end subroutine coop_generate_3dzeta_l


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
       cov(i, i) = cosmology%Clzetazeta(l, r(i))
    enddo
    !$omp end parallel do

  end subroutine coop_get_zeta_shells_cov



  subroutine coop_generate_3Dzeta(cosmology, fname, lmax, nr, r, nprocs, proc_id)
    class(coop_cosmology_firstorder)::cosmology
    COOP_INT lmax, nr
    COOP_REAL zlms(-lmax:lmax,nr), r(nr)
    COOP_UNKNOWN_STRING fname
    type(coop_file) fp
    COOP_INT l, blocksize, i, m, start, step, base
    COOP_INT,optional::nprocs, proc_id
    blocksize = (lmax+1)*(lmax+1)
    call coop_feedback("Generating 3D zeta map "//trim(fname)//" (size = "//trim(coop_num2str(blocksize*4*nr))//"Bytes)")
    if(present(nprocs) .and. present(proc_id))then
       if(coop_MPI_Rank().eq.0)then
          if(coop_file_exists(trim(fname))) then
             call coop_return_error("zeta3D", "file already exists", "stop")
          else
             !!create the file if it does not exist
             call system("touch "//trim(fname))
             call fp%open(trim(coop_file_add_postfix(fname,"__info")), "w")
             write(fp%unit, "(2I5)") lmax, nr
             write(fp%unit, "("//trim(coop_num2str(nr))//"E16.7)") r
             call fp%close()
          endif
       endif
       call coop_MPI_Barrier()
       start = mod(proc_id, nprocs)
       step = nprocs
    else
       start = 0
       step = 1
    endif
    call fp%open(trim(fname), "b", recl = 4)  !!open with binary mode
    do l=start, lmax, step
       call coop_generate_3Dzeta_l(cosmology, l, lmax, nr, r , zlms)
       do i=1, nr
          base = blocksize*(i-1)+l*l
          write(fp%unit, ERR=100, REC=base+1) real(zlms(0,i), 4)
          do m = 1, l
             write(fp%unit, ERR=100, REC=base+2*m) real(zlms(-m,i)/coop_sqrt2, 4)
             write(fp%unit, ERR=100, REC=base+2*m+1) real(zlms(m,i)/coop_sqrt2, 4)
          enddo
       enddo
    enddo
    call fp%close()
    return
100  call coop_MPI_Abort( "IO error in coop_generate_3Dzeta")
  end subroutine coop_generate_3Dzeta


  subroutine coop_load_2Dzeta(fname, lmax, nr, ishell, alms)
    COOP_INT lmax, nr, ishell
    real(4) zlms((lmax+1)**2)
    complex alms(0:lmax, 0:lmax)
    COOP_UNKNOWN_STRING::fname
    COOP_INT blocksize
    type(coop_file) fp
    COOP_INT:: l, m
    alms = 0.
    if(ishell .gt. nr .or. ishell .le. 0) call coop_MPI_Abort( "coop_load_2Dmap wrong layer number")
    blocksize = (lmax+1)**2
    call fp%open(trim(fname), 'br', recl = blocksize*4)
    read(fp%unit, ERR=100, REC = ishell) zlms
    call fp%close()
    do l=0, lmax
       alms(l, 0) = zlms(l*l+1)
       do m=1, l
          alms(l,m) = cmplx(zlms(l*l+m*2), zlms(l*l+m*2+1))
       enddo
    enddo
    return
100 call coop_MPI_Abort("coop_load_2Dmap: io error")
  end subroutine coop_load_2Dzeta

!!$mapping 3D zeta to CMB alms
!!The formula is
!!a_{lm} = 2/\pi  \int_0^\infty  j_l(kr) \zeta_{lm}(r)  \Delta_l(k) k^2 dk r^2 dr
!!where \Delta_l(k) = \int S(k, \chi) j_l(k \chi) d\chi is the transfer function computed in get_transfer function (see firstorder.f90)


  subroutine coop_get_zeta_trans_l(source, l, nr, r, trans)
    type(coop_cosmology_firstorder_source)::source
    COOP_INT l, ik, ir, ichi, idense, isrc, i, nr
    COOP_REAL::r(nr), trans(nr, source%nsrc)
    COOP_REAL, dimension(:,:,:),allocatable::trans_chi_r
    COOP_REAL, dimension(:,:,:),allocatable::ampchi, phasechi, ampr, phaser
    allocate(trans_chi_r(source%nsrc, source%ntau, nr), ampchi(coop_k_dense_fac, source%nk, source%ntau), phasechi(coop_k_dense_fac, source%nk, source%ntau), ampr(coop_k_dense_fac, source%nk, nr), phaser(coop_k_dense_fac, source%nk, nr))
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

    !$omp parallel do private(ir, ik ,idense)
    do ir = 1, nr
       do ik = 2, source%nk
          do idense = 1, coop_k_dense_fac
             call coop_jl_get_amp_phase(l, source%k_dense(idense, ik)*r(ir), ampr(idense, ik, ir), phaser(idense, ik, ir))
          enddo
       enddo
    enddo
    !$omp end parallel do

    !$omp parallel do private(ir, ichi)
    do ir = 1, nr
       do ichi = 1, source%ntau
          call coop_get_zeta_trans_l_step(source, l, nr, r, ichi, ir, trans_chi_r(:, ichi, ir), ampchi, phasechi, ampr, phaser)
       enddo
    enddo
    !$omp end parallel do

    !$omp parallel do private(isrc, ir)
    do isrc = 1, source%nsrc
       do ir = 1, nr
          trans(ir, isrc) = (2.d0/coop_pi)*sum(trans_chi_r(isrc, :, ir)*source%dtau)*source%chi(ir)**2
       enddo
    enddo
    !$omp end parallel do

    deallocate(trans_chi_r,ampchi, phasechi, ampr, phaser)

  end subroutine coop_get_zeta_trans_l


  !!compute \int_0^\infty S(k, \chi) j_l(k \chi) j_l (k r) k^2 dk
  subroutine coop_get_zeta_trans_l_step(source,  l, nr, r, ichi, ir, trans, ampchi, phasechi, ampr, phaser)
    type(coop_cosmology_firstorder_source)::source
    COOP_INT::l, ichi, ir, nr
    COOP_REAL::r(nr)
    COOP_REAL::trans(source%nsrc)
    COOP_REAL,parameter::width = coop_pi*8.d0
    COOP_REAL::ampchi(coop_k_dense_fac, source%nk, source%ntau), phasechi(coop_k_dense_fac, source%nk, source%ntau), ampr(coop_k_dense_fac, source%nk, nr), phaser(coop_k_dense_fac, source%nk, nr)
    COOP_INT ik, idense
    COOP_REAL chidiff,  dk_diff, kmin, xmin, kmax
    chidiff = abs(source%chi(ichi)-r(ir))
    dk_diff = 1.d0/max(chidiff, 1.d-5)
    trans = 0.d0
    call coop_jl_startpoint(l, xmin)
    kmin = xmin/ max(source%chi(ichi), r(ir))
    kmax = max((l+0.5d0)/max(min(source%chi(ichi), r(ir)), 0.1d0)*6.d0, 1.d3)
    do ik = 2, source%nk
          if(source%k(ik) .le. kmin )cycle
          if(source%k_dense(1,ik) .gt. kmax)exit
          if(abs(phasechi(idense, ik, ichi)-phaser(idense, ik, ir)).gt. (18.d0*width))exit
          if(source%dk_dense(1, ik) .gt. dk_diff) exit
          do idense = 1, coop_k_dense_fac
             trans = trans + COOP_INTERP_SOURCE(source, :, idense, ik, ichi) * ampchi(idense, ik, ichi) * ampr(idense, ik, ir)*( cos(phasechi(idense, ik, ichi) -phaser(idense, ik, ir))*(1.d0+ tanh(10.d0-abs(phasechi(idense, ik, ichi)-phaser(idense, ik, ir))/width)) + cos(phasechi(idense, ik, ichi)+phaser(idense, ik, ir))*(1.d0+ tanh(10.d0-(phasechi(idense, ik, ichi)+phaser(idense, ik, ir))/width)) )*source%k_dense(idense, ik)**2*source%dk_dense(idense, ik)*(1.d0+tanh(20.d0*(0.5d0-source%k_dense(idense, ik)/kmax)))/8.d0
          enddo
       enddo
 !!$   endif
  end subroutine coop_get_zeta_trans_l_step

end module coop_cls_postprocess_mod
