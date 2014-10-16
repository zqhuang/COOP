module coop_lensing_mod
  use coop_firstorder_mod
  use coop_pertobj_mod
  use coop_wrapper_background
  implicit none
#include "constants.h"

!!adapted from CAMB/lensing.f90

  private

  type coop_bessfuncs
     type(coop_function) b(4)
   contains
     procedure::init => coop_bessfuncs_init
  end type coop_bessfuncs


  public coop_get_lensing_Cls, coop_bessfuncs


contains



  subroutine coop_get_lensing_Cls(lmin, lmax, Cls, Cls_lensed)
    !Do flat sky approx partially non-perturbative coop_lensing, coop_lensing_method=2

    !!flat sky approx partially non-perturbative lensing
    !!adapted from CAMB lensing_method=2
    COOP_INT::lmin, lmax
    COOP_INT, parameter::l_buffer = 300
    COOP_REAL::Cls(coop_num_Cls,lmin:lmax), Cls_lensed(coop_num_Cls, lmin:lmax)
    type(coop_bessfuncs)::bf
    COOP_REAL, parameter::coop_lensing_reduce_fac = 1.d0
    integer l, i
    integer :: npoints 
    COOP_REAL Cgl2,  sigmasq, theta
    COOP_REAL dtheta
    COOP_REAL  C2term, expsig, corr1, corr2, corr3, corr4, fac, fac1, fac2
    COOP_REAL,dimension(lmin:lmax+l_buffer) :: Cphil3, CTT, CTE,CEE, lfacs
    COOP_REAL bess(4, lmin:lmax+l_buffer)
    integer maxind
    Cls_lensed = 0.d0
    npoints = min(max(lmax*3, 1000), 8000)
    maxind = nint(npoints*coop_lensing_reduce_fac) - 1
    dtheta = coop_pi / npoints
    call bf%init(coop_pi * (lmax+l_buffer))
    !$omp parallel do private(l, fac)
    do l = lmin, lmax
       Cphil3(l) = Cls(coop_index_ClLenLen, l) * l * (l+1.d0) * (l+0.5d0) /(coop_2pi)
       fac = l/coop_2pi
       CTT(l) = Cls(coop_index_ClTT, l)*fac
       CEE(l) = Cls(coop_index_ClEE, l)*fac
       CTE(l) = Cls(coop_index_ClTE, l)*fac
       lfacs(l) = l/2.d0 * l
    end do
    !$omp end parallel do

    fac = exp(-6.d0/l_buffer)

    do l = lmax+1, lmax+l_buffer
       Cphil3(l) = Cphil3(l-1)* fac
       CTT(l) = CTT(l-1)*fac
       CEE(l) = CEE(l-1)*fac
       CTE(l) = CTE(l-1)*fac
       lfacs(l) = l/2.d0*l
    end do


    do i=1, maxind
       theta = i * dtheta 
       sigmasq =0
       Cgl2=0

       !$omp parallel do private(l) reduction(+:sigmasq, Cgl2)
       do l = lmin, lmax+l_buffer
          bess(:, l) =  coop_function_multeval_bare(4, bf%b, l*theta)
          sigmasq = sigmasq + (1-bess(1, l))*Cphil3(l) 
          Cgl2 =  Cgl2 + bess(2, l)*Cphil3(l)
       end do
       !$omp end parallel do

       !Get difference between lensed and unlensed correlation function
       corr1 = 0
       corr2 = 0
       corr3 = 0
       corr4 = 0.
       !$omp parallel do private(expsig, C2term, fac1, fac2, fac) reduction(+:corr1, corr2, corr3, corr4)
       do l = lmin, lmax 
          !For 2nd order perturbative result use 
          !         expsig = 1 -sigmasq*l**2/2._dl
          !         C2term = l**2*Cgl2/2._dl
          fac = sigmasq*lfacs(l)
          expsig = exp(-fac) 
          C2term = Cgl2*lfacs(l)
          !Put theta factor later  in here
          fac1 = expsig*theta
          fac2 = C2term*fac1
          fac1 = fac1 - theta  !we want expsig-1 to get lensing difference

          fac = fac1*bess(1, l) + fac2*bess(2, l) 

          !TT
          corr1 = corr1 + CTT(l) * fac                              

          !Q + U
          corr2 = corr2 + CEE(l) * fac                              
          fac2 = fac2/2.d0
          !Q-U
          corr3 = corr3 + CEE(l) * &
               (fac1*bess(3, l) + fac2*(bess(2, l)+bess(4, l)))                               
          !Cross
          corr4 = corr4 + CTE(l) * &
               (fac1*bess(2, l) + fac2*(bess(1, l)+bess(3, l)))                               
       end do
       !$omp end parallel do
       corr1 = corr1 * dtheta*coop_2pi
       corr2 = corr2 * dtheta*coop_pi
       corr3 = corr3 * dtheta*coop_pi
       corr4 = corr4 * dtheta*coop_2pi


       Cls_lensed(coop_index_ClTT, :) =  Cls_lensed(coop_index_ClTT, :) + corr1 * bess(1, :)
       Cls_lensed(coop_index_ClEE, :) =  Cls_lensed(coop_index_ClEE, :) + (corr2 *Bess(1, :) + corr3 *Bess(3, :))
       Cls_lensed(coop_index_ClBB, :) =  Cls_lensed(coop_index_ClBB, :) + (corr2*Bess(1, :) -corr3*Bess(3, :))
       Cls_lensed(coop_index_ClTE, :) =  Cls_lensed(coop_index_ClTE, :) + corr4*Bess(2, :)
    end do
    !$omp parallel do
    do l=lmin, lmax
       Cls_lensed(coop_index_ClBB, l) =  Cls_lensed(coop_index_ClBB, l) * coop_lensing_BB_corr(l)
    enddo
    !$omp end parallel do
  end subroutine Coop_get_lensing_Cls


  subroutine coop_bessfuncs_init(bf, MaxArg)
    class(coop_bessfuncs)::bf
    COOP_REAL, intent(in):: MaxArg
    integer i, max_bes_ix
    COOP_REAL::x
    COOP_REAL, parameter::dbessel = 0.05d0
    COOP_REAL, allocatable, dimension(:) ::  bess0, bess2, bess4, bess6
    max_bes_ix = nint(MaxArg / dbessel) + 3
    allocate(Bess0(max_bes_ix), Bess2(max_bes_ix),Bess4(max_bes_ix), Bess6(max_bes_ix))
    Bess0(1)=1.d0
    Bess2(1)=0.d0; Bess4(1)=0.d0; Bess6(1)=0.d0
    !$omp parallel do private(x)
    do i=2, max_bes_ix
       x = (i-1)*dbessel
       Bess0(i) = bessel_j0(x)
       Bess2(i) = bessel_jn(2,x)
       Bess4(i) = bessel_jn(4,x)
       Bess6(i) = bessel_jn(6,x)
    end do
    !$omp end parallel do
    x = (max_bes_ix-1)*dbessel
    call bf%b(1)%init(n = max_bes_ix, xmin = 0.d0, xmax =x, f = bess0, check_boundary = .false.)
    call bf%b(2)%init(n = max_bes_ix, xmin = 0.d0, xmax =x, f = bess2, check_boundary = .false.)
    call bf%b(3)%init(n = max_bes_ix, xmin = 0.d0, xmax =x, f = bess4, check_boundary = .false.)
    call bf%b(4)%init(n = max_bes_ix, xmin = 0.d0, xmax =x, f = bess6, check_boundary = .false.)
    deallocate(bess0, bess2, bess4, bess6)

  end subroutine Coop_bessfuncs_init

  !!empirical correction for full sky
  function coop_lensing_BB_corr(l) result(rat)
    COOP_INT l
    COOP_REAL rat
    rat = l/1000.d0
    rat =  1.0761451542474956d0 +        3.1707237388772097d-3 * log(rat)+( - 9.2811615578650142d-2 +  (0.15467835418873549d0 +(- 0.10064923553139785d0 +  2.9276943032231332d-2*rat)*rat)*rat)*rat
  end function coop_lensing_BB_corr


end module coop_lensing_mod

