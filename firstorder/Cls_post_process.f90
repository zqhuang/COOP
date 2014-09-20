module coop_cls_postprocess_mod
  use coop_wrapper_background
  use coop_pertobj_mod
  use coop_firstorder_mod
  implicit none
#include "constants.h"

  private

  COOP_INT,parameter::sp = kind(1.)

  public::coop_generate_3dzeta, coop_load_2Dzeta

contains

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
    !$omp parallel do private(i, j)
    do i=1, nr
       do j=1, i-1
          cov(i, j) = cosmology%Clzetazeta(l, r(i), r(j))
          cov(j, i) = cov(i, j)
       enddo
       cov(i, i) = cosmology%Clzetazeta(l, r(i))
    enddo
    !$omp end parallel do
    call coop_matsym_sqrt(Cov)
    nm = 2*l+1
    do i=1, nr
       call coop_random_get_Gaussian_vector(nm, zetalm(-l:l, i))
    end do
    zetalm(-l:l,:) = matmul(zetalm(-l:l,:), Cov)
  end subroutine coop_generate_3dzeta_l




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


  subroutine coop_get_zeta_trans_l(source, ind, l, trans)
    type(coop_cosmology_firstorder_source)::source
    COOP_REAL::trans(:)
    COOP_INT l, ind, ik, ir, ichi, idense
    COOP_REAL, dimension(:,:),allocatable::trans_chi_r
    COOP_REAL, dimension(:,:,:),allocatable::jls
    allocate(trans_chi_r(source%ntau, source%ntau), jls(coop_k_dense_fac, source%nk, source%ntau))
    jls(:, 1, :) = 0.d0
    call coop_jl_check_init(l=l, save_memory = .true.)
    !$omp parallel do private(ir, ik ,idense)
    do ir = 1, source%ntau
       do ik = 2, source%nk
          do idense = 1, coop_k_dense_fac
             jls(idense, ik, ir) = coop_jl(l, source%k_dense(idense, ik)*source%chi(ir))
          enddo
       enddo
    enddo
    !$omp end parallel do

    !$omp parallel do private(ir, ichi)
    do ir = 1, source%ntau
       do ichi = 1, source%ntau
          call coop_get_zeta_trans_l_r1r2(source, ind, l, ichi, ir, trans_chi_r(ichi, ir))
       enddo
    enddo
    !$omp end parallel do

    !$omp parallel do
    do ir = 1, source%ntau
       trans(ir) = sum(trans_chi_r(:, ir)*source%dtau)*source%chi(ir)**2
    enddo
    !$omp end parallel do

    deallocate(trans_chi_r, jls)

  end subroutine coop_get_zeta_trans_l


  subroutine coop_get_zeta_trans_l_r1r2(source, ind, l, ichi, ir, trans)
    type(coop_cosmology_firstorder_source)::source
    COOP_INT::l, ind, ichi, ir
    COOP_REAL::trans
    trans = 0.d0
    
  end subroutine coop_get_zeta_trans_l_r1r2

!!$    real(dl),parameter::buffer = 80._dl
!!$    real(dl),parameter::buffer2 = 80._dl
!!$    integer l
!!$    integer nk
!!$    real(dl) kmax
!!$    real(dl),dimension(:),allocatable:: k, rres
!!$    real(dl),dimension(:,:),allocatable:: src, amp, phase, jl, map
!!$    integer,dimension(:),allocatable::rloc
!!$    integer ir, ik, ichi,kloc, irmax
!!$    real(dl) lna, kalpha, dk, a,b, aa, bb, x, z2, trans(zu_nchi, 0:zu_lmax)
!!$    integer trans_type
!!$    z2 = sphericalbesselj_zero(l, 2)
!!$    kmax = min(o1cls_kmax*0.999, o1wr_l2kmax(l)+buffer)
!!$    nk = min(floor(kmax/o1wr_kmin)-1, nint((kmax*DISTLSS)*(2. + 2./l))) !!for l = 2, 3 need a bit more coverage
!!$    allocate(k(nk), src(nk, zu_nchi), amp(nk, zu_nchi), phase(nk, zu_nchi), jl(nk, zu_nchi), rloc(zu_nchi), rres(zu_nchi), map(zu_nchi, zu_nchi))
!!$
!!$    dk = kmax/nk
!!$    irmax = 0
!!$    do ir = 1, zu_nchi
!!$       if(zu_chi(ir) .le. H0TAU0 - H0TAU_RECOMB_START)then
!!$          call o1cls_get_rloc(zu_chi(ir), rloc(ir), rres(ir))
!!$       else
!!$          irmax = ir - 1
!!$          exit
!!$       endif
!!$    enddo
!!$    src(:,irmax+1:zu_nchi) = 0
!!$    do ik = 1, nk
!!$       k(ik) = dk*ik
!!$       call o1cls_get_kloc(k(ik), kloc, a, b, aa, bb)
!!$       do ir = 1, irmax
!!$          call o1cls_source_interp(trans_type, rloc(ir), rres(ir), kloc, a, b, aa, bb, src(ik, ir))
!!$       enddo
!!$       do ir = 1, zu_nchi
!!$          x = k(ik)*zu_chi(ir)
!!$          call jlb_get_amp_phase(l, x, amp(ik, ir), phase(ik, ir))
!!$          if(x.ge.z2)then
!!$             jl(ik, ir) = amp(ik, ir)*cos(phase(ik, ir))
!!$          else
!!$             jl(ik, ir) = sphericalbesselJ(l, x) 
!!$          endif
!!$       enddo
!!$       x = (kmax-k(ik)-buffer)*(10.d0/buffer)
!!$       if(x.lt.10.d0) src(ik,:) = src(ik,:)*(1.d0+tanh(x))/2.d0 
!!$       !!smooth the tail to avoid truncation glitches
!!$    enddo
!!$    map = 0.d0
!!$    do ir = 1, irmax
!!$       do ichi = 1, zu_nchi
!!$          map(ir, ichi) = sum(src(:, ir)*k**2*(jl(:, ir)*jl(:,ichi)-amp(:, ir)*amp(:, ichi)*cos(phase(:, ir)+phase(:, ichi)) * (1.d0 - tanh((buffer2 - min(phase(:, ir), phase(:, ichi)))*(10.d0/buffer2)))/4.d0))
!!$       enddo
!!$    enddo
!!$    deallocate(rloc, rres)
!!$    deallocate(k, src, amp, phase, jl) 
!!$    x = (2.d0/const_pi)*O1cls_lfactor(l, trans_type)*dk
!!$    do ichi = 1, zu_nchi
!!$       trans(ichi, l) = x*sum(map(:, ichi)*zu_dchi)*zu_chi(ichi)**2
!!$    enddo
!!$    deallocate(map) 
!!$  end subroutine coop_get_general_trans_l

end module coop_cls_postprocess_mod
