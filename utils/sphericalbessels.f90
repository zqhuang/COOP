module SphericalBessel_utils
  use basic_utils
  use file_io_utils
  use specfunc_utils
  use general_utils
  use jlb_utils
  implicit none

#include "utils.h"

#define JL_SUPPRESS_TAIL NO

  interface SphericalBesselJ
     module procedure SphericalBesselJ_s, SphericalBesselJ_v
  End interface SphericalBesselJ


  interface SphericalBesselJ_Slow
     module procedure SphericalBesselJ_Slow_s, SphericalBesselJ_Slow_v
  End interface SphericalBesselJ_Slow


  interface SphericalBesselJ_fast
     module procedure SphericalBesselJ_fast_s, SphericalBesselJ_fast_v
  End interface SphericalBesselJ_fast



  type jltmptab
     integer l
     real(dl) dx, xmax, xmin, xminmdx, dxsqby6
     integer n
     real(dl),dimension(:),allocatable::jl, jl2
  end type jltmptab

!!*************************************************************************
!!------------------defining the spherical Bessel function table ----------
!!------------------you do not need to know what they are -----------------
!!*************************************************************************
  character(LEN=*),parameter::jl_file = MAINPATH // "/data/jlb_table.dat"

  type(jltmptab),dimension(:),allocatable::sphbess_tmp_tab
  
  Integer::sphbess_lmax = 0
  integer::sphbess_lmin = 0


!!*************************************************************************
!!-------------------------------end of definition ------------------------
!!*************************************************************************

contains

!!*************************************************************************
!!============================Main Routines that are useful ==============
!!*************************************************************************
  function sphericalbesselJ_lowercut(l) result(cut)
    integer l
    real(dl) cut
    call jlb_startpoint(l, cut)
  end function sphericalbesselJ_lowercut


  function sphericalbesselJ_uppercut(l) result(cut)
    integer l
    real(dl) cut
#if JL_SUPPRESS_TAIL
    cut = l * 1.5d0 + (const_2pi*250.d0)
#else
    cut = max(500.d0+200.d0*log10(l+0.5d0), 2.8d0*(l+0.5d0))
#endif
  end function sphericalbesselJ_uppercut


  function SphericalBesselJ_Slow_s(l, x) result(jl)
    integer l
    real(dl) x , jl
    call sphericalBessJ(l, x, jl)
  end function SphericalBesselJ_Slow_s

  function SphericalBesselJ_Slow_v(l, x) result(jl)
    integer l
    real(dl),dimension(:),intent(in):: x 
    real(dl) jl(size(x))
    integer i
    do i=1, size(x)
       call sphericalBessJ(l, x(i), jl(i))
    enddo
  end function SphericalBesselJ_Slow_v


  subroutine sphbess_initialize()  
    !!************************************************************
    !!Initialization; call this subroutine before doing anything.
    !!************************************************************
    integer l, l1, l2
    type(file_pointer) fp
    logical file_loaded
    if(sphbess_lmax .eq. 0 )then
       call specfunc_initialize()
       file_loaded = .false.
       if(file_exists(trim(jl_file)))then
          fp = open_file(trim(jl_file), "ur")
          read(fp%unit,ERR=100) l1, l2
          if(l1 .eq. jlb_lmin .and. l2 .eq. jlb_lmax)then
             do l=jlb_lmin, jlb_lmax
                read(fp%unit, ERR=100)  jlb_table(l)%n
                if(allocated(jlb_table(l)%data))deallocate(jlb_table(l)%data)
                allocate(jlb_table(l)%data(0:2,0:jlb_table(l)%n))
                read(fp%unit, ERR=100) jlb_table(l)%dx,jlb_table(l)%xmin, jlb_table(l)%xmax,jlb_table(l)%data
                jlb_tabulated(l) = .true.
             enddo
             file_loaded = .true.
          endif
          call close_file(fp)       
       endif
100    if(file_loaded)then
          if(global_feedback) write(*,"(A, I5)") "Spherical Bessel table loaded from l=0 to l=", jlb_lmax
       else
          if(global_feedback) write(*,"(A, I5)") "Generating spherical Bessel table from l=0 to l=", jlb_lmax
          !$omp parallel do
          do l=jlb_lmin, jlb_lmax
             call jlb_tabulate(l)
          enddo
          !$omp end parallel do
          if(MPIComm_Rank().eq.0)then
             fp = open_file(trim(jl_file), "u")
             write(fp%unit) jlb_lmin, jlb_lmax
             do l=jlb_lmin, jlb_lmax
                write(fp%unit)  jlb_table(l)%n
                write(fp%unit) jlb_table(l)%dx,jlb_table(l)%xmin, jlb_table(l)%xmax,jlb_table(l)%data
             enddo
             call close_file(fp)
             call feedback_print( "Spherical Bessel table generated in file: "//trim(jl_file))
          endif
       endif
       sphbess_lmin = jlb_lmin
       sphbess_lmax = jlb_lmax
       allocate(sphbess_tmp_tab(0:sphbess_lmax))
       sphbess_tmp_tab%n = 0
       sphbess_tmp_tab%l = -1
    endif
  end subroutine sphbess_initialize


  Function SphericalBesselJ_s(l, x) result(jl)
    !!*************************************************************************
    !!return j_l(x) 
    !!for x>=0; 0<=l<= lmax 
    !!    (lmax is defined when you call sphbess_initialize; default is 10000)
    !!accuracy 10^{-8}
    !!*************************************************************************
    integer l
    real(dl) x, jl
    call jlb_get_jl(l, x, jl)
  End Function SphericalBesselJ_s

  Function SphericalBesselJ_V(l, x) result(jl)
    !!*************************************************************************
    !!return j_l(x) 
    !!for x>=0; 0<=l<= lmax 
    !!    (lmax is defined when you call sphbess_initialize; default is 10000)
    !!accuracy 10^{-8}
    !!*************************************************************************
    integer l
    real(dl),dimension(:),intent(in):: x
    real(dl) jl(size(x))
    integer i
    !$omp parallel do
    do i=1, size(x)
       call jlb_get_jl(l, x(i), jl(i))
    enddo
    !$omp end parallel do
  End Function SphericalBesselJ_V


  function sphericalbesselJPrime(l, x) result(jlp)
    integer l
    real(dl) x, jlp, jl
    call sphericalbessel_jl_jlp(l, x, jl, jlp)
  end function sphericalbesselJPrime

  subroutine SphericalBessel_Jl_Jlp(l, x, jl, jlp) 
    integer l
    real(dl) x, jl, jlp
    call jlb_get_jl_and_jlp(l, x, jl, jlp)
  end subroutine SphericalBessel_Jl_Jlp


  subroutine get_all_jls(x, lmax, jls)
    integer,parameter::larr_size = 133
    integer,dimension(larr_size),parameter::lstart = &
         (/ 11 , 17, 25, 34, 44, 55, 67, 80, 94, 109, 124 &
         , 140, 157, 174, 192, 211, 230, 250, 270, 291, 313, 335 &
         , 358, 381, 405, 429, 454, 479, 505, 531, 557, 584, 611 &
         , 639, 667, 696, 725, 754, 784, 814, 845, 876, 907, 939 &
         , 971, 1003, 1036, 1069, 1103, 1137, 1171, 1206, 1241, 1276, 1312 &
         , 1348, 1384, 1421, 1458, 1495, 1533, 1571, 1609, 1648, 1687, 1726 &
         , 1765, 1805, 1845, 1885, 1926, 1967, 2008, 2050, 2092, 2134, 2176 &
         , 2219, 2262, 2305, 2349, 2393, 2437, 2481, 2526, 2571, 2616, 2661 &
         , 2707, 2753, 2799, 2845, 2892, 2939, 2986, 3033, 3081, 3129, 3177 &
         , 3225, 3274, 3323, 3372, 3421, 3471, 3521, 3571, 3621, 3672, 3723 &
         , 3774, 3825, 3877, 3929, 3981, 4033, 4086, 4139, 4192, 4245, 4298 &
         , 4352, 4406, 4460, 4514, 4569, 4624, 4679, 4734, 4789, 4845, 4901 &
         , 4957 /)
    integer,dimension(larr_size),parameter::lend = &
         (/ 16, 24, 33, 43, 54, 66, 79, 93, 108, 123,  &
         139, 156, 173, 191, 210, 229, 249, 269, 290, 312, 334,  &
         357, 380, 404, 428, 453, 478, 504, 530, 556, 583, 610,  &
         638, 666, 695, 724, 753, 783, 813, 844, 875, 906, 938,  &
         970, 1002, 1035, 1068, 1102, 1136, 1170, 1205, 1240, 1275, 1311,  &
         1347, 1383, 1420, 1457, 1494, 1532, 1570, 1608, 1647, 1686, 1725,  &
         1764, 1804, 1844, 1884, 1925, 1966, 2007, 2049, 2091, 2133, 2175,  &
         2218, 2261, 2304, 2348, 2392, 2436, 2480, 2525, 2570, 2615, 2660,  &
         2706, 2752, 2798, 2844, 2891, 2938, 2985, 3032, 3080, 3128, 3176,  &
         3224, 3273, 3322, 3371, 3420, 3470, 3520, 3570, 3620, 3671, 3722,  &
         3773, 3824, 3876, 3928, 3980, 4032, 4085, 4138, 4191, 4244, 4297,  &
         4351, 4405, 4459, 4513, 4568, 4623, 4678, 4733, 4788, 4844, 4900,  &
         4956, 5000  /)
    real(dl) x
    integer lmax, l, i, nl, nr, nm
    real(dl) jls(0:lmax)
    if(lmax .gt. 5000) stop "lmax > 5000, no jl data"
    !$omp parallel do
    do l=0, min(lmax, 10)
       jls(l) = SphericalBesselJ(l,x)
    enddo
    !$omp end parallel do
    if(x.lt.1.d0 .or. lmax .le. 10)then
       if(lmax.gt.10)jls(11:lmax) = 0.d0
       return
    endif
    nl = 1
    nr = larr_size
    do while(nr - nl > 1)
       nm = (nl+nr)/2
       if(lend(nm) .gt. lmax)then
          nr = nm
       else
          nl = nm
       endif
    enddo
    !$omp parallel do private(i, l)
    do i = 1, nl
       l = lend(i) - 1
       jls(l) = SphericalBesselJ(l, x)
       jls(l-1) = SphericalBesselJ(l-1, x)
       jls(l+1) = (2*l+1)/x*jls(l) - jls(l-1)
       do l = lend(i) - 1, lstart(i) + 1, -1
          jls(l-1) = (2*l+1)/x*jls(l) - jls(l+1)
       enddo
    enddo
    !$omp end parallel do
    jls(lmax) = SphericalBesselJ(lmax, x)
    jls(lmax-1) = SphericalBesselJ(lmax-1, x)
    do l = lmax-1, lstart(nr)+1, -1
       jls(l-1) = (2*l+1)/x*jls(l) - jls(l+1)
    enddo
  end subroutine get_all_jls

  
  !!you need to call sphbess_init_fast(l, xmax) before using this function
  !!xmax must be larger than x
  function SphericalBesselJ_fast_s(l, x) result(jl)
    integer l
    real(dl) x, jl
    call jltmptab_eval(sphbess_tmp_tab(l), x, jl)
  end function SphericalBesselJ_fast_s


  function SphericalBesselJ_fast_v(l, x) result(jl)
    integer l
    real(dl),dimension(:),intent(IN):: x
    real(dl) jl(size(x))
    integer i
    !$omp parallel do
    do i = 1, size(x)
       call jltmptab_eval(sphbess_tmp_tab(l), x(i), jl(i))
    enddo
    !$omp end parallel do
  end function SphericalBesselJ_fast_v


  subroutine jltmptab_init(tab, l, xmax)
    type(jltmptab) tab
    integer l, i
    real(dl) xmax, amp, phase, x
#if JL_SUPPRESS_TAIL
    if(tab%l .eq. l  .and. tab%n .gt. 0 )return
#else
    if(tab%l .eq. l  .and. tab%n .gt. 0 .and. tab%xmax .ge. xmax )return
#endif
    if(allocated(tab%jl))deallocate(tab%jl)
    if(allocated(tab%jl2))deallocate(tab%jl2)
    tab%l = l
    tab%xmin = SphericalBesselJ_lowercut(l)
#if JL_SUPPRESS_TAIL
    tab%xmax = sphericalbesselJ_uppercut(l)
    tab%n = ceiling((tab%xmax - tab%xmin)*10)  
#else
    tab%xmax = max(xmax, l+1.d0)
    tab%n = ceiling((tab%xmax - tab%xmin)*8)
#endif
    allocate(tab%jl(tab%n),tab%jl2(tab%n))
    tab%dx = (tab%xmax-tab%xmin)/(tab%n-1)
    tab%dxsqby6 = tab%dx**2/6.d0
    tab%xminmdx = tab%xmin-tab%dx
#if JL_SUPPRESS_TAIL
    !$omp parallel do private (x, phase, amp)
    do i=1,tab%n
       x = tab%xmin + tab%dx * (i-1)
       call jlb_get_amp_phase(l, x, amp, phase)
       tab%jl(i) = SphericalBesselJ(l, x) - amp*cos(phase)*(1.d0+tanh((x-tab%xmax+200.d0)*0.05d0))/2.d0
    enddo
    !$omp end parallel do
#else
    !$omp parallel do
    do i=1,tab%n
       tab%jl(i) = SphericalBesselJ(l, tab%xmin + tab%dx * (i-1))
    enddo
    !$omp end parallel do
#endif
    call spline_uniform(tab%xmin, tab%dx, tab%n, tab%jl, tab%jl2)
  end subroutine jltmptab_init

  subroutine jltmptab_free(tab)
    type(jltmptab) tab
    if(tab%n.gt.0 .and. tab%l .ge. 0)then
       if(allocated(tab%jl))deallocate(tab%jl)
       if(allocated(tab%jl2))deallocate(tab%jl2)
       tab%n = 0
       tab%l = -1
    endif
  end subroutine jltmptab_free

  function jltmptab_val(tab, x) result(jl)
    real(dl) jl,x
    type(jltmptab) tab
    call splint_uniform_bd(tab%xminmdx, tab%dx, tab%dxsqby6, tab%n, tab%jl, tab%jl2, x, jl, 0.d0, 0.d0)
  end function jltmptab_val

  subroutine jltmptab_eval(tab, x,jl)
    real(dl) jl,x
    type(jltmptab) tab
    call splint_uniform_bd(tab%xminmdx, tab%dx, tab%dxsqby6, tab%n, tab%jl, tab%jl2, x, jl, 0.d0, 0.d0)
  end subroutine jltmptab_eval


  subroutine sphbess_init_fast(l, xmax)
    integer l
    real(dl) xmax
    call jltmptab_init(sphbess_tmp_tab(l), l, xmax)
  end subroutine sphbess_init_fast

  subroutine sphbess_end_fast(l)
    integer l
    call jltmptab_free(sphbess_tmp_tab(l))
  end subroutine sphbess_end_fast


end module SphericalBessel_utils
