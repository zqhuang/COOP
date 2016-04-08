module coop_special_function_mod
  use coop_wrapper_typedef
  implicit none
#include "constants.h"

  private


  public:: coop_log2, coop_sinc, coop_sinhc, coop_asinh, coop_acosh, coop_atanh, coop_InverseErf, coop_InverseErfc, coop_Gaussian_nu_of_P, coop_bessj, coop_sphericalbesselJ, coop_sphericalBesselCross, Coop_Hypergeometric2F1, coop_gamma_product, coop_sqrtceiling, coop_sqrtfloor, coop_bessI, coop_legendreP, coop_cisia, coop_Ylm, coop_normalized_Plm, coop_IncompleteGamma, coop_threej000, coop_ThreeJSymbol, coop_ThreeJ_Array, coop_bessI0, coop_bessi1, coop_bessJ0, coop_bessJ1, coop_sphere_correlation, coop_sphere_correlation_init, coop_get_normalized_Plm_array, FT_Gaussian3D_window, FT_spherical_tophat, coop_pseudoCl_kernel_index_TT, coop_pseudoCl_kernel_index_TE, coop_pseudoCl_kernel_index_TB, coop_pseudoCl_kernel_index_EB, coop_pseudoCl_kernel_index_EE_plus_BB, coop_pseudoCl_kernel_index_EE_minus_BB, coop_TEB_index_T, coop_TEB_index_E, coop_TEB_index_B, coop_TEB_index_TT, coop_TEB_index_EE, coop_TEB_index_BB, coop_TEB_index_TE, coop_TEB_index_TB, coop_TEB_index_EB, coop_int3j,  coop_pseudoCl_matrix, coop_pseudoCl2Cl, coop_pseudoCl_get_kernel,  coop_pseudoCl_matrix_pol, coop_pseudoCl2Cl_pol, coop_pseudoCl_get_kernel_pol,coop_next_l, coop_nl_range, coop_set_ells


  !!define the index of kernels
  COOP_INT, parameter::coop_pseudoCl_kernel_index_TT = 1
  COOP_INT, parameter::coop_pseudoCl_kernel_index_TE = 2
  COOP_INT, parameter::coop_pseudoCl_kernel_index_TB = 2
  COOP_INT, parameter::coop_pseudoCl_kernel_index_EB = 3
  COOP_INT, parameter::coop_pseudoCl_kernel_index_EE_minus_BB = 3
  COOP_INT, parameter::coop_pseudoCl_kernel_index_EE_plus_BB = 4
  COOP_INT, parameter::coop_TEB_index_T = 1
  COOP_INT, parameter::coop_TEB_index_E = 2
  COOP_INT, parameter::coop_TEB_index_B = 3
  COOP_INT, parameter::coop_TEB_index_TT = COOP_MATSYM_INDEX(3, coop_TEB_index_T, coop_TEB_index_T)
  COOP_INT, parameter::coop_TEB_index_EE = COOP_MATSYM_INDEX(3, coop_TEB_index_E, coop_TEB_index_E)
  COOP_INT, parameter::coop_TEB_index_BB = COOP_MATSYM_INDEX(3, coop_TEB_index_B, coop_TEB_index_B)
  COOP_INT, parameter::coop_TEB_index_TE = COOP_MATSYM_INDEX(3, coop_TEB_index_T, coop_TEB_index_E)
  COOP_INT, parameter::coop_TEB_index_TB = COOP_MATSYM_INDEX(3, coop_TEB_index_T, coop_TEB_index_B)
  COOP_INT, parameter::coop_TEB_index_EB = COOP_MATSYM_INDEX(3, coop_TEB_index_E, coop_TEB_index_B)


  interface coop_InverseErf
     module procedure coop_InverseErf_s, coop_InverseErf_v
  end interface coop_InverseErf

  interface coop_InverseErfc
     module procedure coop_InverseErfc_s, coop_InverseErfc_v
  end interface coop_InverseErfc


  interface coop_log2
     module procedure coop_log2_s, coop_log2_v
  end interface coop_log2

  interface coop_sinc
     module procedure coop_sinc_s, coop_sinc_v
  end interface coop_sinc

  interface coop_sinhc
     module procedure coop_sinhc_s, coop_sinhc_v
  end interface coop_sinhc

  interface coop_asinh
     module procedure coop_asinh_s, coop_asinh_v
  end interface coop_asinh

  interface coop_acosh
     module procedure coop_acosh_s, coop_acosh_v
  end interface coop_acosh

  interface coop_atanh
     module procedure coop_atanh_s, coop_atanh_v
  end interface coop_atanh

  interface coop_sphericalbesselJ
     module procedure coop_sphericalBesselJ_s, coop_sphericalBesselJ_v
  end interface coop_sphericalbesselJ

  
  
contains

  function coop_InverseErfc_s(x) result(ierfc)
    COOP_REAL::x, ierfc, y, piysq, dx
    if(x .ge. 2.d0 .or. x .le. 0.d0) stop "InverseErfc: argument overflow"
    if( x .lt. 0.5d0)then
       ierfc = sqrt(-log(coop_sqrtpi/2.d0 * x))
       ierfc = sqrt(-log(coop_sqrtpi * ierfc  * x ))
       ierfc = sqrt(-log(coop_sqrtpi * ierfc  * x ))
       do
          dx =  (erfc(ierfc) - x)
          if(dx .eq. 0.d0)return
          y = exp(ierfc**2+log(abs((coop_sqrtpi/2.d0)*dx)))
          ierfc = ierfc + sign( y, dx)
          if(abs(y) .lt. 1.d-12) return
       enddo
       return
    endif
    if( x .gt. 1.5d0)then
       ierfc = -sqrt(-log(coop_sqrtpi/2.d0 * (2.d0-x)))
       ierfc = -sqrt(-log(-coop_sqrtpi * ierfc * (2.d0-x) ))
       ierfc = -sqrt(-log(-coop_sqrtpi * ierfc * (2.d0-x) ))
       if(ierfc .gt. 20.d0)return
       do
          dx =  (erfc(ierfc) - x)
          if(dx .eq. 0.d0)return
          y = exp(ierfc**2+log(abs((coop_sqrtpi/2.d0)*dx)))
          ierfc = ierfc + sign( y, dx)
          if(abs(y) .lt. 1.d-12) return
       enddo
       return
    endif
    y = 1.d0 - x
    piysq = coop_pi*y**2
    ierfc = (coop_sqrtpi/2.d0)*y*( &
         1.d0 + piysq * ( 1.d0/12.d0 + &
         piysq * ( 7.d0/480.d0 + &
         piysq * ( 127.d0/40320.d0 + &
         piysq * (4369.d0/5806080.d0 + &
         piysq * ( 34807.d0/182476800.d0 + &
         piysq*( 20036983.d0/3985293312.d2 + &
         piysq*(2280356863.d0/167382319104.d3 + &
         piysq*(3.7684666922417d-6 + &
         piysq*(1.0591652893189d-6 + &
         piysq*(3.0149710040390d-7 + &
         piysq*(8.6717030296252d-8 + &
         piysq*(2.5157456448260d-8 + &
         piysq*(7.3516855472913d-9 + &
         piysq/(1.d0-piysq*0.29564d0)*(2.1617630349901d-9  &
         )))))))))))))))
    return
  end function coop_InverseErfc_s

  function coop_InverseErf_s(x) result(ierf)
    COOP_REAL::x
    COOP_REAL::ierf
    ierf = coop_InverseErfc(1.d0 - x)
  end function coop_InverseErf_s

  function coop_InverseErfc_v(x) 
    COOP_REAL, dimension(:)::x
    COOP_REAL::coop_InverseErfc_v(size(x))    
    COOP_INT i
    !$omp parallel do
    do i = 1, size(x)
       coop_InverseErfc_v(i)  = coop_InverseErfc_s(x(i))
    enddo
    !$omp end parallel do
  end function coop_InverseErfc_v

  function coop_InverseErf_v(x) 
    COOP_REAL, dimension(:)::x
    COOP_REAL::coop_InverseErf_v(size(x))    
    COOP_INT i
    !$omp parallel do
    do i = 1, size(x)
       coop_InverseErf_v(i)  = coop_InverseErfc_s(1.d0-x(i))
    enddo
    !$omp end parallel do
  end function coop_InverseErf_v

  function coop_Gaussian_nu_of_P(P) result(nu)
    !!given P = Probability(x<nu sigma), return nu
    COOP_REAL::P, nu
    if(P.lt. 0.5d0)then
       nu =  - coop_sqrt2 * coop_InverseErfc_s( 2.d0*P )
    else
       nu = coop_sqrt2 * coop_InverseErfc_s( 2.d0*(1.d0-P) )       
    endif
  end function coop_Gaussian_nu_of_P

  function coop_asinh_s(x)
    COOP_REAL x
    COOP_REAL coop_asinh_s
    if(x.gt. 1.d-6)then
       coop_asinh_s = log(x+sqrt(1.d0+x**2))
    else
       coop_asinh_s = x * (1.d0 - x*x/6.d0)
    endif
  end function coop_asinh_s

  function coop_asinh_v(x)
    COOP_REAL,dimension(:),intent(IN):: x
    COOP_REAL coop_asinh_v(size(x))
    COOP_INT i
    !$omp parallel do
    do i=1, size(x)
       coop_asinh_v(i) = coop_asinh_s(x(i))
    enddo
    !$omp end parallel do
  end function coop_asinh_v

  function coop_acosh_s(x)
    COOP_REAL x
    COOP_REAL coop_acosh_s
    if(x.gt. 1.000001d0)then
       coop_acosh_s = log(x+sqrt((x-1.d0)*(x+1.d0)))
    else
       coop_acosh_s = sqrt(2.d0*(x-1.d0))*(1.d0 - (x-1.d0)/12.d0)
    endif
  end function coop_acosh_s

  function coop_acosh_v(x)
    COOP_REAL,dimension(:),intent(IN):: x
    COOP_REAL coop_acosh_v(size(x))
    COOP_INT i
    !$omp parallel do
    do i=1, size(x)
       coop_acosh_v(i) = coop_acosh_s(x(i))
    enddo
    !$omp end parallel do
  end function coop_acosh_v

  function coop_atanh_s(x)
    COOP_REAL x
    COOP_REAL coop_atanh_s
    coop_atanh_s = 0.5d0*log((1.d0+x)/(1.d0-x))
  end function coop_atanh_s

  function coop_atanh_v(x)
    COOP_REAL,dimension(:),intent(IN):: x
    COOP_REAL coop_atanh_v(size(x))
    coop_atanh_v = 0.5d0*log((1.d0+x)/(1.d0-x))
  end function coop_atanh_v  

  subroutine continued_fraction_recurrence(an, bn, s)
    COOP_REAL an, bn, s(4)
    s = (/ s(3), s(4), bn*s(3)+an*s(1), bn*s(4)+an*s(2) /)
  end subroutine continued_fraction_recurrence

  Function coop_sinc_s(x) 
    COOP_REAL x,coop_sinc_s
    if(abs(x).lt.0.001)then
       coop_sinc_s=1.d0 - x*x/6.d0
    else
       coop_sinc_s=dsin(x)/x
    endif
  End Function coop_sinc_s

  function coop_sinc_v(x)
    COOP_REAL,dimension(:),intent(in)::x
    COOP_REAL coop_sinc_v(size(x))
    COOP_INT i
    !$omp parallel do
    do i=1, size(x)
       coop_sinc_v(i) = coop_sinc_s(x(i))
    end do
    !$omp end parallel do
  end function coop_sinc_v

  Function Coop_sinhc_s(x)
    COOP_REAL x,coop_sinhc_s
    if(abs(x).lt.0.001)then
       coop_sinhc_s=1.d0 + x*x/6.d0
    else
       coop_sinhc_s=dsinh(x)/x
    end if
  End Function Coop_sinhc_s


  function coop_sinhc_v(x)
    COOP_REAL,dimension(:),intent(in)::x
    COOP_REAL coop_sinhc_v(size(x))
    COOP_INT i
    !$omp parallel do
    do i=1, size(x)
       coop_sinhc_v(i) = coop_sinhc_s(x(i))
    end do
    !$omp end parallel do
  end function coop_sinhc_v

  function coop_digamma(x) result(psi)
    COOP_REAL x, psi, y
    if(abs(x-nint(x)).lt. coop_tiny .and. nint(x).le.0)then
       psi = coop_gamma_product( (/ x, 1-x /), (/ 1, 1 /) )
       if(mod(nint(x),2).eq.0) psi = -psi
       return
    endif
    psi = 0 
    y = x
    do while( y .lt. 10.d0)
       psi = psi - 1.d0/y
       y = y + 1.d0
    enddo
    y = y - 0.5d0
    psi = psi + dlog(y) 
    y = 1.d0/y**2
    psi = psi + y*((1.d0/24.d0) + y*(-(7.d0/960.d0) + y*((31.d0/8064.d0) + y*(-(127.d0/30720.d0) + y*((511.d0/67584.d0)+y*(-1414477.d0/67092480.d0)) ))))
  end function coop_digamma


  function coop_gamma_product(x, p, lnnorm) result(f)
    COOP_REAL f
    COOP_REAL,dimension(:)::x
    COOP_INT,dimension(:)::p
    COOP_INT i, n, m, sgn
    COOP_REAL eps, lnf, pifac
    COOP_REAL,optional::lnnorm
    n = size(x)
    if(n.ne. size(p)) stop "invalid input in coop_gamma_product"
    sgn = 1
    lnf = 0
    do i=1,n
       if(x(i).gt.0.d0)then
          lnf = lnf + log_gamma(x(i))*p(i)
       else
          m = nint(x(i))
          if(mod(m,2).eq.0)then
             eps = x(i) - m
          else
             eps = m - x(i)
          endif
          if(abs(eps).lt. coop_tiny)then
             if(mod(m,2).eq.0)then
                eps = coop_tiny
             else
                eps = -coop_tiny
             endif
          endif
          pifac = dsin(coop_pi*eps)/coop_pi
          lnf  = lnf + (log_gamma(1.d0-x(i)) + dlog(abs(pifac)))*(-p(i))
          if(pifac.lt.0.d0 .and. mod(p(i),2) .ne. 0) sgn = -sgn
       endif
    enddo
    if(present(lnnorm))then
       f = dexp(lnf+lnnorm)*sgn
    else
       f = dexp(lnf)*sgn
    endif
  end function coop_gamma_product


  !! InCompleteGamma(a,x)=\int_x^\infty t^{a-1}e^{-t}dt, for X>=0, A>0 
  !! Numerical Recipies
  function coop_InCompleteGamma(A,X) result(incompletegamma)
    COOP_REAL InCompleteGamma
    COOP_REAL A,X
    if(X.LT.(A+1.D0))then
       InCompleteGamma=exp(Log_gamma(A))
       if(X.ne.0.d0)InCompleteGamma=InCompleteGamma-D_Gamma1(A,X)
    else
       InCompleteGamma=D_Gamma2(A,X)
    endif
  contains
    function D_Gamma1(A,X)
      COOP_REAL D_Gamma1,A,X
      COOP_REAL TEMP,ssum,TEMPA
      COOP_REAL,PARAMETER::PREC=1.D-9
      if(X.eq.0.D0)THEN
         D_Gamma1=0.D0
         return
      endif
      TEMP=1.D0/A
      ssum=TEMP
      TEMPA=A
      DO WHILE(TEMP.gt.ssum*PREC)
         TEMPA=TEMPA+1
         TEMP=TEMP*X/TEMPA
         ssum=ssum+TEMP
      endDO
      D_Gamma1=ssum*dexp(-X+A*dlog(X))
    end function D_Gamma1

    function D_Gamma2(A,X)
      COOP_REAL D_Gamma2,A,X
      COOP_REAL AN,B,C,D,H,DEL
      COOP_REAL,PARAMETER::PREC=1.D-9
      COOP_REAL,PARAMETER::FPMIN=1.D-50
      COOP_INT i
      B=X+1.D0-A 
      C=1./FPMIN
      D=1./B
      H=D
      I=1
      DEL=0.D0
      DO WHILE(Dabs(DEL-1.D0).gt.PREC)
         AN=-I*(I-A)
         B=B+2.D0
         D=AN*D+B
         if(Dabs(D).LT.FPMIN)D=FPMIN
         C=B+AN/C
         if(Dabs(C).LT.FPMIN)C=FPMIN
         D=1.D0/D
         DEL=D*C
         H=H*DEL
         I=I+1
      endDO
      D_Gamma2=dexp(-X+A*dlog(X))*H
    end function D_Gamma2
  end function Coop_InCompleteGamma



  subroutine coop_cisia(X,CI,SI)
    implicit COOP_REAL (A-H,O-Z)
    implicit COOP_INT (I-N)
    !
    !       =============================================
    !       Purpose: Compute cosine and sine integrals
    !                Si(x) and Ci(x)  ( x ò 0 )
    !       Input :  x  --- Argument of Ci(x) and Si(x)
    !       Output:  CI --- Ci(x)
    !                SI --- Si(x)
    !       =============================================
    !
    DIMENSION BJ(101)
    P2=1.570796326794897D0
    EL=.5772156649015329D0
    EPS=1.0D-15
    X2=X*X
    if (X.eq.0.0D0) THEN
       CI=-1.0D+300
       SI=0.0D0
    else if (X.LE.16.0D0) THEN
       XR=-.25D0*X2
       CI=EL+dlog(X)+XR
       DO  K=2,40
          XR=-.5D0*XR*(K-1)/(K*K*(2*K-1))*X2
          CI=CI+XR
          if (Dabs(XR).LT.Dabs(CI)*EPS) GO TO 15
       endDO
15     XR=X
       SI=X
       DO  K=1,40
          XR=-.5D0*XR*(2*K-1)/K/(4*K*K+4*K+1)*X2
          SI=SI+XR
          if (Dabs(XR).LT.Dabs(SI)*EPS) return
       endDO
    else if (X.LE.32.0D0) THEN
       M=int(47.2+.82*X)
       XA1=0.0D0
       XA0=1.0D-100
       DO  K=M,1,-1
          XA=4.0D0*K*XA0/X-XA1
          BJ(K)=XA
          XA1=XA0
          XA0=XA
       endDO
       XS=BJ(1)
       DO  K=3,M,2
          XS=XS+2.0D0*BJ(K)
       Enddo
       BJ(1)=BJ(1)/XS
       DO  K=2,M
          BJ(K)=BJ(K)/XS
       Enddo
       XR=1.0D0
       XG1=BJ(1)
       DO  K=2,M
          XR=.25D0*XR*(2.0*K-3.0)**2/((K-1.0)*(2.0*K-1.0)**2)*X
          XG1=XG1+BJ(K)*XR
       Enddo
       XR=1.0D0
       XG2=BJ(1)
       DO  K=2,M
          XR=.25D0*XR*(2.0*K-5.0)**2/((K-1.0)*(2.0*K-3.0)**2)*X
          XG2=XG2+BJ(K)*XR
       Enddo
       XCS=DCOS(X/2.0D0)
       XSS=DSIN(X/2.0D0)
       CI=EL+dlog(X)-X*XSS*XG1+2*XCS*XG2-2*XCS*XCS
       SI=X*XCS*XG1+2*XSS*XG2-DSIN(X)
    else
       XR=1.0D0
       XF=1.0D0
       DO  K=1,9
          XR=-2.0D0*XR*K*(2*K-1)/X2
          XF=XF+XR
       Enddo
       XR=1.0D0/X
       XG=XR
       DO  K=1,8
          XR=-2.0D0*XR*(2*K+1)*K/X2
          XG=XG+XR
       Enddo
       CI=XF*DSIN(X)/X-XG*DCOS(X)/X
       SI=P2-XF*DCOS(X)/X-XG*DSIN(X)/X
    endif
    return
  end subroutine coop_cisia


  !!for x>0 
  function sphericalbesselj0(x) result(jl)
    COOP_REAL x, jl
    if(x.lt. 0.02d0)then
       jl =1.d0- x**2*((1.d0/6.d0)- x**2*(1.d0/120.d0))
    else
       jl = dsin(x)/x
    endif
  end function sphericalbesselj0

  function sphericalbesselj1(x) result( jl)
    COOP_REAL x, jl
    if(x.lt. 0.08d0)then
       jl =x*((1.d0/3.d0)- x**2*((1.d0/30.d0) - x**2*(1.d0/840.d0)))
    else
       jl = (dsin(x)/x-dcos(x))/x
    endif
  end function sphericalbesselj1

  function sphericalbesselj2(x) result( jl)
    COOP_REAL x, jl
    if(x.lt. 0.2d0)then
       jl = x**2*((1.d0/15.d0)-x**2*((1.d0/210.d0)-x**2*(1.d0/7560.d0)))
    else
       jl = (-3.0d0*dcos(x)/x-dsin(x)*(1.d0-3.d0/x**2))/x
    endif
  end function sphericalbesselj2

  function sphericalbesselj3(x) result( jl)
    COOP_REAL x, jl
    if(x.lt. 0.6d0)then
       jl = x**3*((1.d0/105.d0)-x**2*((1.d0/1890.d0)-x**2*(1.d0/83160.d0)))
    else
       jl = (dcos(x)*(1.d0-15.d0/x**2)-dsin(x)*(6.d0-15.d0/x**2)/x)/x
    endif
  end function sphericalbesselj3


  function sphericalbesselj4(x) result(jl)
    COOP_REAL x, jl
    if(x.lt. 0.9d0)then
       jl = x**4*((1.d0/945.d0)-x**2*((1.d0/20790.d0)-x**2*(1.d0/1081080.d0)))
    else
       jl = (dsin(x)*(1.d0-(45.d0-105.d0/x**2)/x**2)+dcos(x)*(10.D0-105.D0/x**2)/x)/x
    endif
  end function sphericalbesselj4


  subroutine coop_SphericalBesselJ_s(l, x, jl)
    !!== MODifIED subroutine FOR SPHERICAL bessEL functionS.                       ==!!
    !!== CORRECTED THE SMALL BUGS IN PACKAGE CMBFAST&CAMB(for l=4,5, x~0.001-0.002)==!! 
    !!== CORRECTED THE SIGN OF J_L(X) FOR X<0 CASE                                 ==!!
    !!== WORKS FASTER AND MORE ACCURATE FOR LOW L, X<<L, AND L<<X cases            ==!! 
    !!== zqhuang@cita.utoronto.ca                                                 ==!!
    COOP_INT l
    COOP_REAL X,JL
    COOP_REAL AX,AX2
    COOP_REAL,PARAMETER::LN2=0.6931471805599453094D0
    COOP_REAL,PARAMETER::ONEMLN2=0.30685281944005469058277D0
    COOP_REAL,parameter::ROOTPI12 = 21.269446210866192327578D0
    COOP_REAL,parameter::Gamma1 =   2.6789385347077476336556D0 !/* Gamma function of 1/3 */
    COOP_REAL,parameter::Gamma2 =   1.3541179394264004169452D0 !/* Gamma function of 2/3 */
    COOP_REAL NU,NU2,BETA,BETA2,COSB
    COOP_REAL sx,sx2
    COOP_REAL cotb,cot3b,cot6b,secb,sec2b
    COOP_REAL trigarg,expterm,L3

    if(L.LT.0)THEN
       call coop_return_error('coop_sphericalbesselJ', 'Can not evaluate Spherical Bessel Function with index l<0', 'stop')
    endif
    AX=Dabs(X)
    AX2=AX**2
    if(l.lt.7)then
       select case(l)
       case(0)
          if(AX.LT.1.D-1)THEN
             JL=1.D0-AX2/6.D0*(1.D0-AX2/20.D0)
          else
             JL=DSIN(AX)/AX
          endif

       case(1)
          if(AX.LT.2.D-1)THEN
             JL=AX/3.D0*(1.D0-AX2/10.D0*(1.D0-AX2/28.D0))
          else
             JL=(DSIN(AX)/AX-DCOS(AX))/AX
          endif
       case(2)
          if(AX.LT.3.D-1)THEN
             JL=AX2/15.D0*(1.D0-AX2/14.D0*(1.D0-AX2/36.D0))
          else
             JL=(-3.0D0*DCOS(AX)/AX-DSIN(AX)*(1.D0-3.D0/AX2))/AX
          endif
       case(3)
          if(AX.LT.4.D-1)THEN
             JL=AX*AX2/105.D0*(1.D0-AX2/18.D0*(1.D0-AX2/44.D0))
          else
             JL=(DCOS(AX)*(1.D0-15.D0/AX2)-DSIN(AX)*(6.D0-15.D0/AX2)/AX)/AX
          endif
       case(4)
          if(AX.LT.6.D-1)THEN
             JL=AX2**2/945.D0*(1.D0-AX2/22.D0*(1.D0-AX2/52.D0))
          else
             JL=(DSIN(AX)*(1.D0-(45.D0-105.D0/AX2)/AX2)+DCOS(AX)*(10.D0-105.D0/AX2)/AX)/AX
          endif
       case(5)
          if(AX.LT.1.D0)THEN
             JL=AX2**2*AX/10395.D0*(1.D0-AX2/26.D0*(1.D0-AX2/60.D0))
          else
             JL=(DSIN(AX)*(15.D0-(420.D0-945.D0/AX2)/AX2)/AX-DCOS(AX)*(1.D0-(105.D0-945.0d0/AX2)/AX2))/AX
          endif
       case(6)
          if(AX.LT.1.D0)THEN
             JL=AX2**3/135135.D0*(1.D0-AX2/30.D0*(1.D0-AX2/68.D0))
          else
             JL=(DSIN(AX)*(-1.D0+(210.D0-(4725.D0-10395.D0/AX2)/AX2)/AX2)+ &
                  DCOS(AX)*(-21.D0+(1260.D0-10395.D0/AX2)/AX2)/AX)/AX
          endif
       end select
    else
       NU=0.5D0+L
       NU2=NU**2
       if(AX.LT.1.D-40)THEN
          JL=0.D0
       elseif((AX2/L).LT.5.D-1)THEN
          JL=dexp(L*dlog(AX/NU)-LN2+NU*ONEMLN2-(1.D0-(1.D0-3.5D0/NU2)/NU2/30.D0)/12.D0/NU) &
               /NU*(1.D0-AX2/(4.D0*NU+4.D0)*(1.D0-AX2/(8.D0*NU+16.D0)*(1.D0-AX2/(12.D0*NU+36.D0))))
       elseif((dble(l)**2/AX).LT.5.D-1)THEN
          BETA = AX- coop_pio2 * (L+1)
          JL=(DCOS(BETA)*(1.D0-(NU2-0.25D0)*(NU2-2.25D0)/8.D0/AX2*(1.D0-(NU2-6.25)*(NU2-12.25D0)/48.D0/AX2)) &
               -DSIN(BETA)*(NU2-0.25D0)/2.D0/AX* (1.D0-(NU2-2.25D0)*(NU2-6.25D0)/24.D0/AX2*(1.D0-(NU2-12.25)* &
               (NU2-20.25)/80.D0/AX2)) )/AX   
       else
          L3=nu**0.325
          if(AX .LT. NU-(1.31 + 5./nu2)*L3) then
             COSB=NU/AX
             SX = Dsqrt(NU2-AX2)
             COTB=NU/SX
             SECB=AX/NU
             BETA=dlog(COSB+SX/AX)
             COT3B=COTB**3
             COT6B=COT3B**2
             SEC2B=SECB**2
             EXPTERM=( (2.D0+3.D0*SEC2B)*COT3B/24.D0 &
                  - ( (4.D0+SEC2B)*SEC2B*COT6B/16.D0 &
                  + ((16.D0-(1512.D0+(3654.D0+375.D0*SEC2B)*SEC2B)*SEC2B)*COT3B/5760.D0 &
                  + (32.D0+(288.D0+(232.D0+13.D0*SEC2B)*SEC2B)*SEC2B)*SEC2B*COT6B/128.D0/NU)*COT6B/NU) &
                  /NU)/NU
             JL =dexp(-NU*BETA+NU/COTB-EXPTERM) / (2.D0*dsqrt(sx*ax))

             !          /**************** Region 2: x >> l ****************/

          elseif (AX .gt. NU+(1.48+ 10./nu2)*L3) then
             COSB=NU/AX
             SX=Dsqrt(AX2-NU2)
             COTB=NU/SX
             SECB=AX/NU
             BETA=DACOS(COSB)
             COT3B=COTB**3
             COT6B=COT3B**2
             SEC2B=SECB**2
             TRIGARG=Sx-NU*BETA-coop_pio4 &
                  -((2.0+3.0*SEC2B)*COT3B/24.D0  &
                  +(16.D0-(1512.D0+(3654.D0+375.D0*SEC2B)*SEC2B)*SEC2B)*COT3B*COT6B/5760.D0/NU2)/NU
             EXPTERM=( (4.D0+sec2b)*sec2b*cot6b/16.D0 &
                  -(32.D0+(288.D0+(232.D0+13.D0*SEC2B)*SEC2B)*SEC2B)*SEC2B*COT6B**2/128.D0/NU2)/NU2
             JL=dexp(-EXPTERM)*dcos(TRIGARG)/dsqrt(ax*sx) 

             !          /***************** Region 3: x near l ****************/

          else
             BETA=AX-NU
             BETA2=BETA**2
             SX=6.D0/AX
             SX2=SX**2
             SECB=SX**0.3333333333333333d0
             SEC2B=SECB**2
             JL=( Gamma1*SECB + BETA*Gamma2*SEC2B &
                  -(BETA2/18.D0-1.D0/45.D0)*BETA*SX*SECB*Gamma1 &
                  -((BETA2-1.D0)*BETA2/36.D0+1.D0/420.D0)*SX*SEC2B*Gamma2   &
                  +(((BETA2/1620.D0-7.D0/3240.D0)*BETA2+1.D0/648.D0)*BETA2-1.D0/8100.D0)*SX2*SECB*Gamma1 &
                  +(((BETA2/4536.D0-1.D0/810.D0)*BETA2+19.D0/11340.D0)*BETA2-13.D0/28350.D0)*BETA*SX2*SEC2B*Gamma2 &
                  -((((BETA2/349920.D0-1.D0/29160.D0)*BETA2+71.D0/583200.D0)*BETA2-121.D0/874800.D0)* &
                  BETA2+7939.D0/224532000.D0)*BETA*SX2*SX*SECB*Gamma1)*Dsqrt(SX)/ROOTPI12
          endif
       endif
    endif
    if(X.LT.0.AND.MOD(L,2).ne.0)JL=-JL
  end subroutine Coop_SphericalBesselJ_s

  subroutine coop_sphericalbesselJ_v(l, xarr, jl)
    COOP_REAL,dimension(:)::xarr
    COOP_REAL jl(:)
    COOP_INT i, l
    !$omp parallel do
    do i=1,size(xarr)
       call coop_SphericalBesselJ_s(l, xarr(i), jl(i))
    end do
    !$omp end parallel do
  end subroutine coop_sphericalbesselJ_v

  function FT_spherical_tophat(kR) 
    !!Fourier transformation of a tophat function in a sphere with radius R
    COOP_REAL kR, FT_spherical_tophat
    if(kR .gt. 0.02)then
       FT_spherical_tophat = (dsin(kR)/kR -  dcos(kR))/ kR**2 * 3.d0
    else
       FT_spherical_tophat = (10.d0 - kR**2*(1.d0 - kR**2/28.d0))/10.d0 
    endif
  end function FT_spherical_tophat


  function FT_Gaussian3D_Window(sigma, k, kp) result(w)
    COOP_REAL w, k, kp, sigma
    w = kp/k/(coop_sqrt2*coop_sqrtpi)/sigma * (exp(-((kp-k)/sigma)**2/2.) -  exp(-((kp+k)/sigma)**2/2.))
  end function FT_Gaussian3D_Window

!!J_0(x)
  function coop_bessj0(x)
    COOP_REAL coop_bessj0,x
    COOP_REAL ax,xx,z
    COOP_REAL p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6, &
         s1,s2,s3,s4,s5,s6,y
    SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4, &
         s5,s6
    DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4, &
         -.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1, &
         .1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
    DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0, &
         651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/,s1,s2, &
         s3,s4,s5,s6/57568490411.d0,1029532985.d0,9494680.718d0, &
         59272.64853d0,267.8532712d0,1.d0/

    if(abs(x).lt.8.d0)then
       y=x**2
       coop_bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y* &
            (s4+y*(s5+y*s6)))))
    else
       ax=abs(x)
       z=8.d0/ax
       y=z**2
       xx=ax-.785398164d0
       coop_bessj0=sqrt(.636619772d0/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y* &
            p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
    endif
  end function coop_bessj0


!!J_1(x)
  function coop_bessj1(X)
    COOP_REAL, intent(in) :: x
    COOP_REAL coop_bessj1,ax,z,xx
    COOP_REAL Y,P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5, &
         R6,S1,S2,S3,S4,S5,S6
    DATA R1,R2,R3,R4,R5,R6/72362614232.D0,-7895059235.D0,242396853.1D0,&
         -2972611.439D0,15704.48260D0,-30.16036606D0/, &
         S1,S2,S3,S4,S5,S6/144725228442.D0,2300535178.D0, &
         18583304.74D0,99447.43394D0,376.9991397D0,1.D0/
    DATA P1,P2,P3,P4,P5/1.D0,.183105D-2,-.3516396496D-4,.2457520174D-5, & 
         -.240337019D-6/, Q1,Q2,Q3,Q4,Q5/.04687499995D0,-.2002690873D-3, &     
         .8449199096D-5,-.88228987D-6,.105787412D-6/
    if(abs(X).LT.8.)then
       Y=X**2
       coop_bessj1=X*(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))) &
            /(S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))
    else
       ax=abs(x)
       Z=8.0d0/ax
       Y=Z**2
       XX=AX-2.356194491d0
       coop_bessj1=sqrt(.636619772d0/AX)*(COS(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y &
            *P5))))-Z*SIN(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5))))) &
            *SIGN(1.d0,x)
    endif
  end function coop_bessj1



!!J_n
  function coop_bessj(n, x) result(bessj)
    COOP_REAL bessj
    COOP_REAL, intent(in) :: x
    COOP_INT, intent(in) :: n
    COOP_INT, parameter :: iacc = 40
    COOP_REAL, parameter :: bigno=1.d10,bigni=1.d-10
    COOP_INT jsum,j,m
    COOP_REAL bj, bjm, bjp, tox, ssum
    if(n .lt. 2)then
       if(n.eq.1)then
          bessj = coop_bessj1(x)
       elseif(n.eq.0)then
          bessj = coop_bessj0(x)
       else
          stop "wrong argument n in bessj"
       endif
       return
    endif
    if( x**2 .lt. n/2.)then
       bessj = (x/2.d0)**n*exp(-log_gamma(n+1.d0))*(1.d0-x**2/(4.d0*(n+1))*(1.d0-x**2/(8.d0*(n+2))*(1.d0-x**2/(12.d0*(n+3))*(1.d0-x**2/(16.d0*(n+4))))))
       return
    endif
    tox = 2/X
    if(x .gt. real(n))then
       BJM=coop_bessj0(X)
       BJ=coop_bessj1(X)
       DO J=1,N-1
          BJP=J*tox*BJ-BJM
          BJM=BJ
          BJ=BJP
       end DO
       bessJ=BJ
    else
       M=2*((N+int(sqrt(real(iacc*N))))/2)
       bessJ=0.0d0
       jsum=0
       ssum=0.d0
       BJP=0.d0
       BJ=1.d0
       DO J=M,1,-1
          BJM=J*tox*BJ-BJP
          BJP=BJ
          BJ=BJM
          if(abs(BJ).gt.bigno)THEN
             BJ=BJ*bigni
             BJP=BJP*bigni
             bessJ=bessJ*bigni
             ssum=ssum*bigni
          endif
          if(jsum.ne.0)ssum=ssum+BJ
          jsum=1-jsum
          if(J.eq.N)bessJ=BJP
       end do
       ssum=2.* ssum-BJ
       bessJ=bessJ/ssum
    endif
  end function Coop_bessj

  ! ----------------------------------------------------------------------
  ! Auxiliary Bessel functions for N=0, N=1
  function coop_bessi0(X)
    double precision X,coop_bessi0,Y,P1,P2,P3,P4,P5,P6,P7,  &
         Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
    DATA P1,P2,P3,P4,P5,P6,P7/1.D0,3.5156229D0,3.0899424D0,1.2067429D0,  &
         0.2659732D0,0.360768D-1,0.45813D-2/
    DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1, &
         0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,  &
         0.2635537D-1,-0.1647633D-1,0.392377D-2/
    if(abs(X).LT.3.75D0) THEN
       Y=(X/3.75D0)**2
       coop_bessi0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
    else
       AX=abs(X)
       Y=3.75D0/AX
       BX=EXP(AX)/sqrt(AX)
       AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
       coop_bessi0=AX*BX
    endif
    return
  end function coop_bessi0
  ! ----------------------------------------------------------------------
  function coop_bessi1(X)
    double precision X,coop_bessi1,Y,P1,P2,P3,P4,P5,P6,P7,  &
         Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
    DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,  &
         0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
    DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1, &
         -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1, &
         -0.2895312D-1,0.1787654D-1,-0.420059D-2/
    if(abs(X).LT.3.75D0) THEN
       Y=(X/3.75D0)**2
       coop_bessi1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
    else
       AX=abs(X)
       Y=3.75D0/AX
       BX=EXP(AX)/sqrt(AX)
       AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
       coop_bessi1=AX*BX
    endif
    return
  end function coop_bessi1


  function coop_bessI(N,X)
    !from http://perso.orange.fr/jean-pierre.moreau/Fortran/tbessi_f90.txt
    !
    !     This subroutine calculates the first kind modified Bessel function
    !     of COOP_INT order N, for any real X. We use here the classical
    !     recursion formula, when X > N. For X < N, the Miller's algorithm
    !     is used to avoid overflows. 
    !     REFERENCE:
    !     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL functionS,
    !     MATHEMATICAL TABLES, VOL.5, 1962.
    COOP_INT, intent(in) :: N
    COOP_INT, PARAMETER :: iacc = 40
    COOP_INT m,j
    double precision, parameter ::  bigno = 1.D10, bigni = 1.D-10
    double precision X,coop_bessI,tox,biM,bi,biP
    if(n.lt.2)then
       if (N.eq.0) then
          coop_bessI = coop_bessi0(X)
          return
       endif
       if (N.eq.1) then
          coop_bessI = coop_bessi1(X)
          return
       endif
    endif
    if(abs(x) .lt. 1.d-20) then
       coop_bessI=0.d0
       return
    endif
    tox = 2.D0/X
    bip = 0.D0
    bi = 1.D0
    coop_bessI = 0.D0
    M = 2*((N+int(sqrt(real(iacc*N)))))
    DO J = M,1,-1
       bim = bip+ J*tox*bi
       bip = bi
       bi  = bim
       if (abs(bi).gt.bigno) then
          bi  = bi*bigni
          biP = biP*bigni
          coop_bessI = coop_bessI*bigni
       endif
       if (J.eq.N) coop_bessI = bip
    end DO
    coop_bessI = coop_bessI*coop_bessi0(X)/bi
    return
  end function coop_bessI


  !!return Legendre polynomial P_l(x)
  function coop_legendreP(l, x) result(Pl)
    COOP_INT l, ell
    COOP_REAL x, Pl, x2, Plm1, Plm2
    select case(l)
    case(0)
       Pl = 1.d0
       return
    case(1)
       Pl = x
       return
    case(2)
       Pl = 1.5d0*x*x - 0.5d0
       return
    case(3)
       Pl = (2.5d0 * x * x - 1.5d0) * x
       return
    case(4)
       x2 = x*x
       Pl = 3.d0/8. + x2 * ( -15.d0/4.d0 + x2 * (35.d0/8.d0)  )
       return
    case(5)
       x2 = x*x
       Pl = x * ( 15.d0/8. + x2 * ( -35.d0/4.d0 + x2 * (63.d0/8.)))
       return
    case(6)
       x2 = x*x
       Pl = -5.d0/16. + x2*( 105.d0/16. + x2*( -315.d0/16.+ x2*(231.d0/16.)))
       return
    case(7)
       x2 = x*x
       Pl = x * ( -35.d0/16. + x2 * ( 315.d0/16.d0 + x2 *( -693.d0/16. +  x2 *(429.d0/16.))))
       return
    case(8)
       x2 = x*x
       Pl = 35.d0/128. + x2 * ( -315.d0/32.d0 + x2 * (3465.d0/64.d0 + x2 *(-3003.d0/32. + x2 * (6435.d0/128.))))
       return
    case(9)
       x2 = x*x
       Pl = x * (315.d0/128. + x2 * ( -1155.d0/32. + x2 * ( 9009.d0/64. + x2 * ( -6435.d0/32. +  x2* (12155.d0/128.)))))
       return
    case(10)
       x2 = x * x
       Pl = -63.d0/256. + x2 * (3465.d0/256. + x2 * (-15015.d0/128. + x2 * (45045.d0/128. + x2 * ( -109395.d0/256.  + x2 * (46189.d0/256.)))))
       return
    case(11)
       x2 = x*x
       Pl = x * (-693.d0/256. +  x2 * ( 15015.d0/256. + x2 * (-45045.d0/128.+ x2 * ( 109395.d0/128. + x2 * (-230945.d0/256. + x2 * (88179.d0/256.))))))
       return
    case(12)
       x2 = x*x
       Pl = 231.d0/1024. + x2 * ( -9009.d0/512.+ x2 * (225225.d0/1024.+ x2 * (-255255.d0/256 + x2 * (2078505.d0/1024 + x2 * (-969969.d0/512. + x2 * (676039.d0/1024.))))))
       return
    case(13)
       x2 = x*x
       Pl = x*(3003.d0 /1024 + x2*(- 45045.d0/512. + x2*(765765.d0/1024 + x2*(-692835.d0/256. + x2*(4849845.d0/1024 + x2*(-2028117.d0/512 + x2*(1300075.d0/1024.)))))))
       return
    case(14)
       x2 = x*x
       Pl = -429.d0/2048 + x2*(45045.d0/2048 + x2*(-765765.d0/2048 + x2* (4849845.d0/2048 + x2*(-14549535.d0/2048. + x2*(22309287.d0/2048. + x2*(-16900975.d0/2048. + x2*(5014575.d0/2048.)))))))
       return
    case(15)
       x2 = x*x
       Pl=x*(-6435.d0/2048. + x2*(255255.d0/2048 + x2*(-2909907.d0/2048. + x2*(14549535.d0/2048. + x2*(-37182145.d0/2048. + x2*(50702925.d0/2048 + x2*(-35102025.d0/2048. + x2*(9694845.d0/2048.))))))))
    case(16)
       x2 = x*x
       Pl = 6435.d0/32768. + x2*(-109395.d0/4096. + x2*(4849845.d0/8192. +x2*(-20369349.d0/4096. + x2*(334639305.d0/16384. + x2*(-185910725.d0/4096. + x2*(456326325.d0/8192. + x2*(-145422675.d0/4096. + x2*(300540195.d0/32768.))))))))
    case(17)
       x2=x*x
       Pl = x*(109395.d0/32768. + x2*(-692835.d0/4096. + x2*(20369349.d0/8192. + x2*(-66927861.d0/4096. + x2*(929553625.d0/16384. + x2*(-456326325.d0/4096. + x2*(1017958725.d0/8192. + x2*(-300540195.d0/4096. + x2*(583401555.d0/32768.)))))))))
    case(18)
      x2=x*x
      Pl = -12155.d0/65536. + x2*(2078505.d0/65536. + x2*(-14549535.d0/16384. + x2*(156165009.d0/16384.+ x2*( -1673196525.d0/32768 + x2*(5019589575.d0/32768. + x2*(-4411154475.d0/16384. + x2*(4508102925.d0/16384. + x2*(-9917826435.d0/65536. + x2*(2268783825.d0/65536.)))))))))
    case default
       Plm1 = 1.d0
       Pl = x
       do ell = 2, l
          Plm2 = Plm1
          Plm1 = Pl
          Pl = ((2*ell -1)*x*Plm1 - (ell -1)*Plm2)/ell
       enddo
       return
    end select
  end function coop_legendreP

  subroutine coop_sphere_correlation_init(lmax, als, bls)
    COOP_INT lmax, l
    COOP_REAL als(0:lmax), bls(0:lmax)
    als(0) = 1.d0/coop_4pi
    bls(0) = 0.d0
    do l = 1, lmax
       als(l) = (2*l+1.d0)/l
       bls(l) = -(l-1.d0)/(2*l-3.d0)*(2*l+1.d0)/l
    enddo
  end subroutine coop_sphere_correlation_init

  function coop_Sphere_Correlation(lmax, Cls, als, bls, x) result(s) !!x = cos(theta)  !!I assume lmax>=1
    !!als and bls are auxilliary arrays. You need to initialize them using Sphere_Correlation_init(lmax, als, bls)
    COOP_INT lmax
    COOP_REAL Cls(0:lmax), als(0:lmax), bls(0:lmax), x
    COOP_REAL s
    COOP_REAL Pls(2)
    COOP_INT l
    Pls = (/ 0.d0, als(0) /)
    s = Cls(0)*Pls(2)
    do l = 1, lmax
       Pls = (/ Pls(2), als(l) * x * Pls(2) + bls(l)*Pls(1) /)
       s = s + Cls(l)*Pls(2)
    enddo
  end function Coop_Sphere_Correlation

  !! return sqrt(4 pi / (2l + 1) ) * Y_l^m (arccos x, 0)
  recursive function Coop_normalized_plm(l, m, x) result(Plm)
    COOP_INT l, m
    COOP_REAL x, x2, Plm
    select case(m)
    case(0)
       Plm = coop_legendreP(l, x)
       return
    case(1)
       select case(l)
       case(0)
          Plm = 0.d0
       case(1)
          Plm = - sqrt((1.-x**2)/2.)
       case(2)
          Plm = - sqrt(1.5*(1.-x**2)) * x
       case(3)
          Plm = sqrt(3.*(1.-x**2))*(0.25 - 1.25*x**2)
       case(4)
          Plm = sqrt(5.*(1.- x**2)) * x * (0.75 - 1.75 * x**2)
       case(5)
          x2 = x**2
          Plm = sqrt(7.5*(1.-x2))*(-1.d0/8. + x2 * (14.d0/8.d0  + x2 * (-21.d0/8.)))
       case(6)
          x2 = x*x
          Plm = sqrt(10.5*(1.-x2)) * x * ( -5.d0/8.d0 + x2 * (30.d0/8.d0 + x2 * (-33.d0/8.d0)))
       case(7)
          x2 = x*x
          Plm = sqrt(3.5*(1.-x2)) * (5./32.d0 + x2*(-135.d0/32. + x2*(495./32.d0 + x2 * (-429./32.d0))))
       case(8)
          x2 = x*x
          Plm = sqrt(2.*(1.-x2))*x*(105.d0/64. + x2 * (-1155.d0/64. + x2 * (3003.d0/64.+ x2* (-2145.d0/64.))))
       case(9)
          x2 = x*x
          Plm = sqrt(2.5*(1.-x2))*(-21./128.d0 + x2 * (231.d0/32. + x2 * (-3003.d0/64. + x2 * (3003.d0/32. + x2 * (-7293.d0/128.)))))
       case(10)
          x2 = x*x
          Plm = sqrt(27.5*(1.-x2)) * x * (-63./128.d0 + x2 *(273.d0/32. + x2 * (-2457.d0/64. + x2 *(1989.d0/32. + x2 * (-4199.d0/128.)))))
       case(11)
          x2 = x*x
          Plm = sqrt(33.d0*(1.-x2))*(21./512.d0 + x2*(-1365.d0/512. + x2*(6825.d0/256. + x2 *(-23205.d0/256. + x2 * (62985.d0/512. + x2 * (-29393.d0/512.))))))
       case(12)
          x2 = x*x
          Plm = sqrt(39.d0 * (1.-x2))*x*(231./512.d0 + x2 *(-5775./512.d0 + x2 * ( 19635.d0/256.+ x2 * (-53295.d0/256. + x2 * (124355.d0/512. + x2 * (-52003.d0/512.))))))
       case(13)
          x2 = x*x
          Plm = sqrt(45.5*(1.-x2)) * ( -33.d0/1024. + x2 * (1485.d0/512. + x2 * (-42075.d0/1024. + x2 * (53295.d0/256. + x2 * (-479655.d0/1024. + x2*(245157.d0/512. + x2 * (-185725.d0/1024.)))))))
       case(14)
          x2 = x*x
          Plm = sqrt(52.5*(1.-x2)) * x * (-429.d0/1024. + x2 * (7293.d0/512. + x2 * (-138567.d0/1024. + x2 * (138567.d0/256. + x2 * (-1062347.d0/1024. + x2 * (482885.d0/512. + x2 * (-334305.d0/1024.)))))))
       case(15)
          x2 = x*x
          Plm = sqrt(15.*(1.-x2))*(429.d0/8192.+x2*(-51051.d0/8192. + x2*(969969.d0/8192. + x2*(-6789783.d0/8192 + x2*(22309287.d0/8192. + x2*(-37182145.d0/8192. + x2*(30421755.d0/8192. + x2*(-9694845.d0/8192))))))))
       case(16)
          x2 = x*x
          Plm = sqrt(17.*(1.-x2))*x*(6435.d0/8192.+x2*(-285285.d0/8192. + x2*(3594591.d0/8192.+x2*( -19684665.d0/8192. + x2*(54679625.d0/8192. + x2*(-80528175.d0/8192. + x2*(59879925.d0/8192. +x2*(-17678835.d0/8192.))))))))
       case(17)
          x2 = x*x
          Plm = sqrt(8.5d0*(1.-x2))*(-2145.d0/32768.d0 + x2*(40755.d0/4096.d0 +x2*(-1996995.d0/8192.d0 + x2*(9186177.d0/4096. +x2*(-164038875.d0/16384.d0 + x2*(98423325.d0/4096.+x2*( -259479675.d0/8192. + x2*(88394175.d0/4096.+x2*( -194467185.d0/32768.)))))))))
       case(18)
          x2 = x*x
          Plm  = sqrt(9.5*(1.-x2))*x*(-36465.d0/32768. + x2*(255255.d0/4096.+x2*(-8219211.d0/8192. + x2*(29354325.d0/4096. + x2*(-440314875.d0/16384. + x2*(232166025.d0/4096. +x2*(- 553626675.d0/8192. + x2* (173996955.d0/4096. + x2*(-358229025.d0/32768.)))))))))
       case default
          goto 100
       end select
       return
    case(2)
       select case(l)
       case(0:1)
          Plm = 0.d0
          return
       case(2)
          Plm = (0.5d0*sqrt(3./2.d0))* ( 1.d0-x*x)
       case(3)
          Plm = (0.5d0 * sqrt(15./2.d0))*x*(1.d0 -x *x)
       case(4)
          Plm = (-0.25d0 * sqrt(2.5d0))*(1.-x**2)*(1.- 7.*x**2)
       case(5)
          Plm = (-0.25d0 * sqrt(105.d0/2.d0)) * x * (1.- x**2) * (1. - 3.*x**2)
       case(6)
          x2 = x*x
          Plm = (sqrt(105.d0) /32.d0)*(1.-x2) * (1. + x2 *(-18.d0 + x2 * 33.d0))
       case(7)
          x2 = x * x
          Plm = (sqrt(21.d0)/32.)*x*(1. - x2)*(15. + x2*(-110. + x2 * 143.))
       case(8)
          x2 = x*x
          Plm = (3.d0/64.*sqrt(35.d0))*(1.-x2)*(-1.+ x2*(33. + x2 *143.* (-1.+ x2)))
       case(9)
          x2 = x*x
          Plm = (3.d0/64. * sqrt(55.d0))*x*(1.-x2)*(-7.+x2*(91.+x2*(-273.+x2*221.)))
       case(10)
          x2 = x*x
          Plm = (sqrt(165.d0/2.)/256.)*(1-x2)*(7. + x2 * 13. * (-28. + x2*(210. + x2 * (-476. + x2*323.))))
       case(11)
          x2 = x*x
          Plm = (sqrt(2145.d0/2.)/256.)* x*(1.-x2)*(21. + x2 * (-420. + x2 * (2142. + x2 * (-3876. + x2 * 2261.))))
       case(12)
          x2 = x*x
          Plm = (sqrt(3003.d0/2.)/512.)*(1.-x2)*(-3.+x2*(225. + 17.*x2 * (-150.+ x2*(570. + x2*(-855. + x2 * 437.)))))
       case(13)
          x2 = x*x
          Plm = (sqrt(455.d0/2.)/512.) * x * (1.- x2)*(-99.+x2*17.*(165. + x2*19.*(-66. + x2*(198. + x2*(-253. + x2*115.)))))
       case(14)   
          x2 = x*x
          Plm = (sqrt(1365.d0/2.)/4096.d0) * (1.-x2)*(33. + x2 * 17.*(-198.+ x2*19.*(165.+x2*(-924.+x2*23.*(99. + x2*(-110.+x2*45.))))))
       case(15)
          x2 = x*x
          Plm = (sqrt(1785.d0/2.)/4096.d0)*x*(1.-x2)*(429.d0 + x2*(-16302.d0 + x2*(171171.d0 + x2*(-749892.d0 + x2*(1562275.d0 + x2*(-1533870.d0 + 570285.d0*x2))))))
       case(16)
          x2 = x*x
          Plm = (sqrt(255.d0/2.))*(1.-x2)*(-143/8192.d0 + x2*(19019.d0/8192.d0 + x2*( -399399.d0/8192.d0 + x2*(3062059.d0/8192.d0 + x2*(-10935925.d0/8192.d0 + x2*(19684665.d0/8192.d0+ x2*(-17298645.d0/8192.d0 + x2*(5892945.d0/8192.d0))))))))
       case(17)
          x2 = x*x
          Plm = (sqrt(323.d0/2.))*(1.-x2)*x*(-2145.d0/8192 + x2*(105105.d0/8192. + x2*( -1450449.d0/8192. + x2*(8633625.d0/8192. + x2*(-25900875.d0/8192. + x2*(40970475.d0/8192. + x2*(-32566275.d0/8192. + x2*(10235115.d0/8192.))))))))
       case(18)
          x2 = x*x
          Plm = (sqrt(1615.d0/2.))*(1.-x2)*(429.d0/65536.d0 + x2*(-9009.d0/8192.d0 + x2*(483483.d0/16384.d0 + x2*(-2417415.d0/8192. + x2*(46621575.d0/32768.d0 + x2*(-30045015.d0/8192. + x2*(84672315.d0/16384 + x2*( -30705345.d0/8192. + x2*(71645805.d0/65536.d0)))))))))
       case default
          goto 100
       end select
       return
    case(3)
       select case(l)
       case(0:2)
          Plm = 0.d0
       case(3)
          Plm = (-coop_sqrt5/4.d0)*(1.d0-x**2)**1.5d0
       case(4)
          Plm = (-sqrt(35.d0)/4.d0)*x*(1.d0-x**2)**1.5d0
       case(5)
          Plm = (sqrt(35.d0)/16.d0)*(1.d0-9.d0*x**2)*(1.d0-x**2)**1.5d0
       case(6)
          Plm = (sqrt(105.d0)/16.d0)*x*(3.d0-11.d0*x**2)*(1.d0-x**2)**1.5d0
       case(7)
          x2 = x*x
          Plm = (-sqrt(10.5d0)/32.d0)*(3 + x2*(- 66 + x2*143))*(1.d0-x2)**1.5d0
       case(8)
          x2 =x*x
          Plm = (-sqrt(577.5d0)/32.d0)*x*(3 +x2*(- 26 + x2*39))*(1.d0-x2)**1.5d0
       case(9)
          x2 = x*x
          Plm = (sqrt(1155.d0)/128.d0)*(1 + x2*(- 39+ x2*(195 + x2*(- 221.d0))))*(1.d0-x2)**1.5d0
       case(10)
          x2 = x*x
          Plm = (Sqrt(2145.d0)/128.d0)*x*(7 + x2*(-105.d0 + x2* (357.d0 + x2*(- 323.d0))))*(1.d0-x2)**1.5d0
       case default
          goto 100
       end select
       return
    case(4)
       select case(l)
       case(0:3)
          Plm = 0.d0
       case(4)
          Plm = (sqrt(17.5d0)/8.d0)*(1.d0-x**2)**2
       case(5)
          Plm =  (3.d0*sqrt(17.5d0)/8.d0)*x*(1.d0-x**2)**2
       case(6)
          Plm =  (3*Sqrt(3.5d0)/16.d0)*(-1.d0 + 11*x**2) *(1.d0-x**2)**2
       case(7)
          Plm =  (sqrt(115.5d0)/16.d0)*x*(-3 + 13*x**2)*(1.d0-x**2)**2
       case(8)
          x2 = x*x
          Plm = (3*sqrt(38.5d0)/64.d0)*(1 + x2*(-26 + x2 * 65)) *(1.d0-x2)**2
       case(9)
          x2 = x*x
          Plm =  (3*Sqrt(2502.5d0)/64.d0)*x*(1.d0 + x2*(- 10 + x2*17))*(1.d0-x2)**2
       case(10)
          x2 = x*x
          Plm = (Sqrt(1072.5)/128.d0)*(-1.d0 + x2*(45 + x2*(- 255 + x2*323))) *(1.d0-x2)**2
       case default
          goto 100
       end select
       return
    case(5)
       select case(l)
       case(0:4)
          Plm = 0.d0
          return
       case(5)
          Plm = (-3*Sqrt(7.d0)/16.)*(1.d0-x**2)**2.5d0
       case(6)
          Plm =  (-3*Sqrt(77.d0)/16.)*x*(1.d0-x**2)**2.5d0
       case(7)
          Plm = (Sqrt(115.5d0)/32.d0)*(1 - 13*x**2)*(1.d0-x**2)**2.5d0
       case(8)
          Plm = (-3*Sqrt(500.5d0)/32.d0)*x*(-1 + 5*x**2)*(1.d0-x**2)**2.5d0
       case(9)
          x2 = x*x
          Plm = (-3*Sqrt(143.d0)/128.d0)*(1 + x2*(-30 + x2*85))*(1.d0-x2)**2.5d0
       case(10)
          x2 = x*x
          Plm = (-Sqrt(429.d0)/128.d0)*x*(15 + x2*(- 170 + 323*x2))*(1.d0-x2)**2.5d0
       case default
          goto 100
       end select
       return
    case(6)
       select case(l)
       case(0:5)
          Plm = 0.d0
       case(6)
          Plm = (sqrt(231.d0)/32.d0)*(1.d0-x**2)**3
       case(7)
          Plm = (Sqrt(3003.d0)/32.d0)*x*(1.d0-x**2)**3
       case(8)
          Plm = (Sqrt(429.d0)/64.d0)*(-1 + 15*x**2)*(1.d0-x**2)**3
       case(9)
          Plm = (Sqrt(2145.d0)/64.d0)*x*(-3 + 17*x**2)*(1.d0-x**2)**3
       case(10)
          x2 = x*x
          Plm = (Sqrt(2145.d0)/512.d0)*(3 + x2*( - 102 + x2*323))*(1.d0-x2)**3
       case default
          goto 100
       end select
       return
    case(7)
       select case(l)
       case(0:6)
          Plm = 0.d0
       case(7)
          Plm = (-sqrt(214.5d0)/32.d0)*(1.d0-x**2)**3.5d0
       case(8)
          Plm = (-3*Sqrt(357.5d0)/32.d0)*x*(1.d0-x**2)**3.5d0
       case(9)
          Plm = (-3*Sqrt(715.d0)/256.d0)*(-1 + 17*x**2)*(1.d0-x**2)**3.5d0
       case(10)
          Plm = (sqrt(36465.d0)/256.d0)*x*(3 - 19*x**2)*(1.d0-x**2)**3.5d0
       case default
          goto 100
       end select
       return
    case(8)
       select case(l)
       case(0:7)
          Plm = 0.d0
       case(8)
          Plm = (3*Sqrt(357.5d0)/128.d0)*(1.d0-x**2)**4
       case(9)
          Plm =(3*Sqrt(6077.5d0)/128.d0) *x*(1.d0-x**2)**4
       case(10)
          Plm = (sqrt(6077.5d0)/256.d0)*(-1.d0+19.d0*x**2)*(1.d0-x**2)**4
       case default
          goto 100
       end select
       return
    case(9)
       select case(l)
       case(0:8)
          Plm = 0.d0
       case(9)
          Plm = (-sqrt(12155.d0)/256.d0)*(1.d0-x**2)**4.5d0
       case(10)
          Plm = (-sqrt(230945.d0)/256.d0)*x*(1.d0-x**2)**4.5d0
       case default
          goto 100
       end select
       return
    case(10)
       select case(l)
       case(0:9)
          Plm = 0.d0
       case(10)
          Plm = (sqrt(46189.d0)/512.d0)*(1.d0-x**2)**5
       case default
          goto 100
       end select
       return
    case(11:)
       if(l.ge. m)then
          goto 100
       else
          Plm = 0.d0
       endif
       return
    case(:-1)
       Plm = Coop_normalized_plm(l, -m, x)
       if(mod(m,2).ne.0)Plm = -Plm
       return
    end select
    return
100 call get_coop_normalized_plm(l,m,x,Plm)    
  End function Coop_normalized_plm


  subroutine get_coop_normalized_plm(l,m,x,Plm)
    COOP_REAL x, Plm_prev, Plm, Plm_next
    COOP_INT l,m, i
    if(m.gt.l)then
       Plm = 0.d0
       return
    endif
    Plm_prev = (4.d0*(1.d0-x**2))**(m/2.d0)*dexp(log_gamma(m+0.5d0)-log_gamma(2*m+1.d0)/2.d0)/coop_sqrtpi
    if(mod(m,2).ne.0) Plm_prev = - Plm_prev
    if(m .eq. l)then
       Plm = Plm_prev
       return
    endif
    Plm =  dsqrt(2*m+1.d0)*x*Plm_prev
    do i=m+2,l
       Plm_next = ((2*i-1)*x*Plm-dsqrt((i-m-1.d0)*(i+m-1.d0))*Plm_prev)/(dsqrt((i-m)*dble(i+m)))
       Plm_prev = Plm
       Plm = Plm_next
    enddo    
  end subroutine get_coop_normalized_plm


  subroutine coop_get_normalized_Plm_array(m, lmax, x, Plms)
    COOP_INT m, lmax, l
    COOP_REAL Plms(0:lmax), x
    if(m.eq.0)then
       Plms(0) = 1.d0
       if(lmax.le.0)return
       Plms(1) = x
       if(lmax.le.1)return
       do l=2, lmax
          Plms(l)  = ((2*l-1)*x*Plms(l-1)- (l-1)*Plms(l-2))/l
       enddo
    else
       Plms(0:min(m-1, lmax)) = 0.d0
       if(lmax .lt. m) return
       Plms(m) = (4.d0*(1.d0-x**2))**(m/2.d0)*dexp(log_gamma(m+0.5d0)-log_gamma(2*m+1.d0)/2.d0)/coop_sqrtpi
       if(mod(m,2).ne.0) Plms(m) = - Plms(m)
       if(lmax .le. m) return    
       Plms(m+1) = sqrt(2*m+1.d0)*x*Plms(m)
       do l= m+2, lmax
          Plms(l)  = ((2*l-1)*x*Plms(l-1)-sqrt((l-m-1.d0)*(l+m-1.d0))*Plms(l-2))/(dsqrt((l-m)*dble(l+m)))
       enddo
    endif
  end subroutine coop_get_normalized_Plm_array

!! return sqrt((l-m)!/(l+m)!) P_l^m(x)
!!no check
  


!!return Y_l^m(\theta, \phi) 
!!input l>=m>=0, theta, phi
  Function coop_Ylm(l,m,theta, phi) result(Ylm)
    COOP_REAL theta,phi
    COOP_COMPLEX Ylm
    COOP_INT l, m
    Ylm = dsqrt((2*l+1)/coop_4pi) * Coop_normalized_plm(l, m, dcos(theta))*cmplx(dcos(m*phi), dsin(m*phi))
  End Function Coop_Ylm

  function coop_legendreP_Approx(l, x) result(Pl)
    COOP_INT l
    COOP_REAL x, Pl, theta, phi, sx
    theta = dacos(x)
    sx = dsin(theta)
    phi = (l+0.5d0)*theta-coop_pio4
    Pl = dsqrt((2.d0/coop_pi)/l/sx)*( &
         (1.d0 - 0.25d0/l)*dcos(phi) &
         + 0.125d0*x/sx/l * dsin(phi) &
         )
  end function coop_legendreP_Approx

  subroutine coop_get_Cn_array(n, Cn)
    COOP_INT n, i
    COOP_REAL Cn(0:n)
    Cn(0) = 1.d0
    do i = 1, n
       Cn(i) = (n-i+1.d0) * Cn(i-1)/i
    enddo
  end subroutine coop_get_Cn_array

  function coop_factor_lm_pm(l, m)  result (f)
    !! sqrt((l+m)!/(l-m)!)
    COOP_INT l, m
    COOP_REAL f
    select case(m)
    case(0)
       f = 1.d0
    case(1)
       f = sqrt(l*(l+1.d0))
    case(2)
       f = sqrt((l-1)*l*(l+1.d0)*(l+2.d0))
    case(3)
       f = sqrt((l-2)*(l-1.d0)*l*(l+1.d0)*(l+2.d0)*(l+3.d0))
    case(4)
       f = sqrt((l-3)*(l-2.d0)*(l-1.d0)*l*(l+1.d0)*(l+2.d0)*(l+3.d0)*(l+4.d0))
    case(5)
       f = sqrt((l-4)*(l-3.d0)*(l-2.d0)*(l-1.d0)*l*(l+1.d0)*(l+2.d0)*(l+3.d0)*(l+4.d0)*(l+5.d0))       
    case default
       f = exp((log_gamma(l+m+1.d0)-log_gamma(l-m+1.d0))/2.d0)
    end select
  end function coop_factor_lm_pm

  function coop_sqrtfloor(n) result(intsqrt)
    COOP_INT n, intsqrt
    intsqrt = floor(sqrt(dble(n)+coop_tiny))
    do while(intsqrt*intsqrt.gt.n)
       intsqrt = intsqrt - 1 !!just to avoid roundoff errors
    enddo
  end function coop_sqrtfloor

  function coop_sqrtceiling(n) result(intsqrt)
    COOP_INT n, intsqrt
    intsqrt = ceiling(sqrt(dble(n)-coop_tiny))
    do while(intsqrt*intsqrt.lt.n)
       intsqrt = intsqrt + 1 !!just to avoid roundoff errors
    enddo
  end function coop_sqrtceiling

  function coop_numbits(m) result(numbits) !!number of bits in binary representation
    COOP_INT n, m, numbits
    numbits = 0
    n = m
    do while(n.gt.0)
       n=n/2
       numbits = numbits+1
    enddo
  end function coop_numbits

  function coop_smoothstep(x, width) result(s)
    COOP_REAL x, width, s
    s = (coop_atanh(x/width)+1.d0)/2.d0
  end function coop_smoothstep

  
  recursive function coop_Hypergeometric2F1(a, b, c, x) result(f)
    !!input 0<=x<=1
    COOP_REAL,parameter::xbound = 0.95d0
    COOP_REAL,parameter::eps = 1.d-3
    COOP_REAL,parameter::tol = 1.d-12
    COOP_REAL a, b, c, x, f , cmab, s, y, cma, cmb, r1, r2, w1, w2,  psia, psib, psi1, lny
    COOP_REAL term
    COOP_INT i, n
    if(x.lt.0.d0 .or. x.gt. 1.d0 .or. c .le. max(0.d0, a+b) )then
       write(*,*) "Hypergemetric2F1: input a+b,c,x = ", a+b,c,x
       write(*,*) "The current version only supports 0<=x<=1 and c>max(0, a+b)"
       stop
    endif
    !!********  x = 1 ***********************
    if(x .gt. 1.d0-tol)then
       f = coop_gamma_product( (/ c, c-a-b, c-a, c-b /), (/ 1, 1, -1, -1 /) ) 
       return
    endif
    !!*********** x <= xbound ******************
    if(x .le. xbound)then
       term = 1.d0
       f = 0
       i = 1
       do while(abs(term) .gt. tol)
          f = f + term
          term = term*(a+i-1)*(b+i-1)/(c+i-1)*x/i
          i = i + 1
       enddo
       f = f + term/(1.d0-x)
       return
    endif
    !!********** xbound < x < 1 *****************
    cmab = c - a - b
    cma = c-a
    cmb = c-b
    y = 1.d0-x
    n = nint(cmab)
    !! -------- if c-a-b is not an COOP_INT ---------
    if(abs(cmab-n).ge. eps*0.999d0 .or. ( n .eq. 0 .and. abs(cmab-n).gt. 1.d-6) )then  
       term = 1.d0
       f = 0
       i = 1
       do while(abs(term).gt.tol)
          f = f + term
          term = term*(a+i-1)*(b+i-1)*y/(i-cmab)/i
          i = i + 1
       enddo
       term = 1.d0
       s = 0
       i = 1
       do while(abs(term) .gt. tol)
          s = s + term
          term = term*(cma+i-1)*(cmb+i-1)/(cmab+i)/i*y
          i = i + 1
       enddo
       f = coop_gamma_product( (/ c, cmab, cma, cmb /), (/ 1, 1, -1, -1 /) ) * (f + coop_gamma_product( (/ cma, cmb, -cmab, a, b, cmab /), (/ 1, 1, 1, -1, -1, -1 /), cmab*dlog(y) )*s)
       return
    endif
    if(n.eq.0)then
       psi1 = -coop_eulerC
       psia = coop_digamma(a)
       psib = coop_digamma(b)
       lny = dlog(y)
       term = 1.d0
       s = 0.d0
       i = 1
       do 
          s  = s + term*(-lny+2*psi1-psia-psib)
          term = term*(a+i-1)*(b+i-1)/i/(i-cmab)*y
          if(abs(term).lt. tol) exit
          psi1 = psi1 + 1.d0/i
          psia = psia + 1.d0/(i-1+a)
          psib = psib + 1.d0/(i-1+b)
          i = i + 1
       enddo
       f = coop_gamma_product( (/ c, a, b /), (/ 1, -1, -1 /) ) * s
       return
    endif
    !!---------- if c-a-b is almost a positive COOP_INT and x > xbound, we do interpolation ---------
    r1 = coop_gamma_product( (/ n+eps, b+eps, a+eps, cmab, cma, cmb /), (/ 1, -1, -1, -1, 1, 1 /) )
    r2 = coop_gamma_product( (/ n-eps, b-eps, a-eps, cmab, cma, cmb /), (/ 1, -1, -1, -1, 1, 1 /) )
    if(r1-r2.ne.0)then
       w1 = (cmab - n + eps)/(2.d0*eps)
       w2 = (n+eps-cmab)/(2.d0*eps)
    else
       w1 = (1-r2)/(r1-r2)
       w2 = (r1-1)/(r1-r2)
       if(abs(w1).gt. 1.5d0 .or. abs(w2).ge.1.5d0)then !!extrapolation not safe
          w1 = (cmab - n + eps)/(2.d0*eps)
          w2 = (n+eps-cmab)/(2.d0*eps)
       endif
    endif
    f = coop_Hypergeometric2F1(a, b, a + b + n + eps, x )*w1 +  coop_Hypergeometric2F1(a, b, a + b + n - eps, x )*w2
    return
  end function Coop_Hypergeometric2F1

  function coop_sphericalBesselCross(ell1, ell2, r1, r2, ns) result(cr)
    !!this function returns the integral
    !!\int_0^infinity j_{ell1}(k r1) j_{ell2}(k r2) (k(r1+r2)/2)^(ns-1) dk/k 
    COOP_INT l1, l2, ell1, ell2
    COOP_REAL r1, r2, ns, cr, R, ratio
    if(ell1+ell2+ns .le. 1.d0 .or. ns .ge. 2.d0)then
       write(*,*) "Invalid input in sphericalBesselCross"
       write(*,*) "l1, l2, r1, r2, n_s = ", ell1, ell2, r1, r2, ns
       stop
    endif
    if(r1 .le. 0.d0 .or. r2 .le. 0.d0)then
       cr = 0
       return
    endif
    if(r1.gt.r2)then
       ratio = r2/r1
       R = r1
       l1 = ell1
       l2 = ell2
    else
       ratio = r1/r2
       R = r2
       l1 = ell2
       l2 = ell1
    endif
    if(4+l1-l2-ns.gt.0.d0)then
       cr = 2.d0**(-4+ns)*coop_pi* (2.d0*R/(r1+r2))**(1-ns)*coop_Hypergeometric2F1((-2-l1+l2+ns)/2.d0,(-1+l1+l2+ns)/2.d0,1.5d0+l2,ratio**2)* dexp(l2*dlog(ratio)+log_gamma((-1+l1+l2+ns)/2.d0) -log_gamma(l2+1.5d0) -log_gamma((4+l1-l2-ns)/2.d0))
    else
       cr = 2.d0**(-4+ns)*dsin(coop_pio2*(4+l1-l2-ns)) * (2.d0*R/(r1+r2))**(1-ns)*coop_Hypergeometric2F1((-2-l1+l2+ns)/2.d0,(-1+l1+l2+ns)/2.d0,1.5d0+l2,ratio**2)* dexp(l2*dlog(ratio)+log_gamma((-1+l1+l2+ns)/2.d0) -log_gamma(l2+1.5d0) + log_gamma(1.d0 - (4+l1-l2-ns)/2.d0))
    endif
  end function coop_sphericalBesselCross


  Function coop_log2_s(x) 
    COOP_REAL x,coop_log2_s
    coop_log2_s= log(x)/coop_ln2
  End Function coop_log2_s

  function coop_log2_v(x)
    COOP_REAL,dimension(:),intent(in)::x
    COOP_REAL coop_log2_v(size(x))
    COOP_INT i
    !$omp parallel do
    do i=1, size(x)
       coop_log2_v(i) = coop_log2_s(x(i))
    end do
    !$omp end parallel do
  end function coop_log2_v


  function coop_int3j(l1, l2, l3, m1, m2, m3) result(thrj)
    COOP_INT::l1, l2, l3, m1, m2, m3
    COOP_REAL::thrj
    thrj = coop_ThreeJSymbol(l1*2, l2*2, l3*2, m1*2, m2*2, m3*2)
  end function coop_int3j

  !!this function is fast and accurate for small m1, m2, m3
  !!for large m's it can be slow and inaccurate.
  recursive function coop_ThreeJSymbol(twoj1, twoj2, twoj3, twom1, twom2, twom3) result(w3j)
    COOP_REAL::w3j
    COOP_INT::twoj1, twoj2, twoj3, twom1, twom2, twom3
    COOP_INT::J, g
    if(twoj1 .lt. 0 .or. twoj2 .lt. 0 .or. twoj3 .gt. twoj1 + twoj2 .or. twoj3 .lt. abs(twoj1 - twoj2) .or. abs(twom1).gt. twoj1 .or. abs(twom2).gt.twoj2 .or. abs(twom3).gt.twoj3 .or. mod(twoj1+twom1,2).ne.0 .or. mod(twoj2+twom2, 2).ne. 0 .or. mod(twoj3 + twom3, 2).ne.0 .or. twom1 + twom2 + twom3 .ne. 0 .or. mod(twoj1 + twoj2 + twoj3,2).ne.0)then
       w3j = 0.d0
       return
    endif
    J = (twoj1 + twoj2 + twoj3)/2
    if(twom1 .eq. 0)then
       if( twom2 .eq. 0) then
          if(mod(J,2).ne.0)then
             w3j = 0.d0
             return
          endif
          g = J/2
          w3j = dexp(&
               log_gamma(g+1.d0)-log_gamma(g-twoj1/2.d0+1.d0) &
               -log_gamma(g-twoj2/2.d0+1.d0) - log_gamma(g-twoj3/2.d0+1.d0) &
               +( &
               log_gamma(J-twoj1+1.d0) + log_gamma(J-twoj2+1.d0) &
               + log_gamma(J-twoj3+1.d0) - log_gamma(J+2.d0) &
               )/2.d0 &
               )
          if(mod(g,2).ne.0) w3j = -w3j
          return
       endif
       if(twom2 .lt. 0)then
          if(twom2 .eq. -2 .and. mod(J,2).eq.0)then
             w3j = -coop_ThreeJSymbol(twoj1, twoj2, twoj3, 0, 0, 0)*(twoj2*(twoj2+2.d0)+twoj3*(twoj3+2.d0)-twoj1*(twoj1+2.d0) )/(2.d0*sqrt(twoj2*(twoj2+2.d0)*twoj3*(twoj3+2.d0)))
             return
          endif
          if(twom2 .lt. -2 )then
             w3j = -(coop_ThreeJSymbol(twoj1, twoj2, twoj3, twom1, twom2+2, twom3-2) * (twoj2*(twoj2+2.d0)+ twoj3*(twoj3+2.d0)-twoj1*(twoj1+2.d0)+2.d0*(twom2+2.d0)*(twom3-2.d0)) &
                  + coop_ThreeJSymbol(twoj1, twoj2, twoj3, twom1, twom2 + 4, twom3 - 4) * sqrt((twoj2-twom2-2.d0)*(twoj2+twom2+4.d0)*(twoj3+twom3-2.d0)*(twoj3-twom3+4.d0)) ) &
                  /sqrt((twoj2-twom2)*(twoj2+twom2+2.d0)*(twoj3+twom3)*(twoj3-twom3+2.d0))
          else
             w3j = (-coop_ThreeJSymbol(twoj1, twoj2-1, twoj3+1, twom1, twom2+1, twom3-1)*sqrt((J-twoj3)*(J-twoj2+1.d0)*(twoj2-twom2)/2.d0) &
                  + coop_ThreeJSymbol(twoj1, twoj2+1, twoj3+1, twom1, twom2+1, twom3-1)*sqrt((J+2.d0)*(J-twoj1+1.d0)*(twoj2+twom2+2.d0)/2.d0))/(sqrt((twoj3-twom3)/2.d0+1.d0)*(twoj2+1.d0))
          endif
          return
       endif
       w3j = coop_ThreeJSymbol(twoj1, twoj3, twoj2, twom1, twom3, twom2)
       if(mod(J,2).ne.0) w3j = -w3j
       return
    endif
    if(twom2 .eq. 0)then
       w3j = coop_ThreeJSymbol(twoj2, twoj3, twoj1, twom2, twom3, twom1)
       return
    endif
    if(twom3 .eq. 0)then
       w3j = coop_ThreeJSymbol(twoj3, twoj1, twoj2, twom3, twom1, twom2)
       return
    endif
    !!all nonzero
    if(abs(twom2) .le. abs(twom1) .and. abs(twom2) .le. abs(twom3))then
       if(twom2 .lt. 0)then
          if(twom1 .lt. twom3)then
             if(twom2 .lt. -2 )then
                w3j = -(coop_ThreeJSymbol(twoj1, twoj2, twoj3, twom1, twom2+2, twom3-2) * (twoj2*(twoj2+2.d0)+ twoj3*(twoj3+2.d0)-twoj1*(twoj1+2.d0)+2.d0*(twom2+2.d0)*(twom3-2.d0)) &
                     + coop_ThreeJSymbol(twoj1, twoj2, twoj3, twom1, twom2 + 4, twom3 - 4) * sqrt((twoj2-twom2-2.d0)*(twoj2+twom2+4.d0)*(twoj3+twom3-2.d0)*(twoj3-twom3+4.d0)) ) &
                     /sqrt((twoj2-twom2)*(twoj2+twom2+2.d0)*(twoj3+twom3)*(twoj3-twom3+2.d0))
             else
                w3j = (-coop_ThreeJSymbol(twoj1, twoj2-1, twoj3+1, twom1, twom2+1, twom3-1)*sqrt((J-twoj3)*(J-twoj2+1.d0)*(twoj2-twom2)/2.d0) &
                     + coop_ThreeJSymbol(twoj1, twoj2+1, twoj3+1, twom1, twom2+1, twom3-1)*sqrt((J+2.d0)*(J-twoj1+1.d0)*(twoj2+twom2+2.d0)/2.d0))/(sqrt((twoj3-twom3)/2.d0+1.d0)*(twoj2+1.d0))
             endif
          else
             w3j = coop_ThreeJSymbol(twoj3, twoj2, twoj1, twom3, twom2, twom1)
             if(mod(J,2).ne.0) w3j = -w3j
          endif
          return
       endif
       w3j = coop_ThreeJSymbol(twoj1, twoj2, twoj3, -twom1, -twom2, -twom3)
       if(mod(J,2).ne.0) w3j = -w3j
       return
    endif
    if(abs(twom1).le. abs(twom2) .and. abs(twom1) .le. abs(twom3))then
       w3j = coop_ThreeJSymbol(twoj3, twoj1, twoj2, twom3, twom1, twom2)
       return
    endif
    w3j = coop_ThreeJSymbol(twoj2, twoj3, twoj1, twom2, twom3, twom1)
    return
  end function coop_ThreeJSymbol


!!ThreeJSymbol[{l1, 0}, {l2, 0}, {l3, 0}]
  function Coop_threej000(l1,l2,l3) result(w3j)
    COOP_REAL  w3j
    COOP_INT  l1,l2,l3, J, g
    J = (l1+l2+l3)
    if(l3 .gt. l1+l2 .or. l3 .lt. abs(l1-l2) .or. mod(J,2).ne.0)then
       w3j = 0.d0
       return
    endif
    g = J/2
    w3j = dexp(&
         log_gamma(g+1.d0)-log_gamma(g-l1+1.d0) &
         -log_gamma(g-l2+1.d0) - log_gamma(g-l3+1.d0) &
         +( &
         log_gamma(J-2*l1+1.d0) + log_gamma(J-2*l2+1.d0) &
         + log_gamma(J-2*l3+1.d0) - log_gamma(J+2.d0) &
         )/2.d0 &
         )
    if(mod(g,2).ne.0) w3j = -w3j
  end function Coop_threej000





  subroutine coop_ThreeJ_Array(threejs,twoj2, twoj3, twom2, twom3, j1min, num_j1s)
    !Recursive evaluation of 3j symbols. Does minimal error checking on input parameters.
    !return in threejs: ThreeJSymbol[{j1, m1},{j2, m2},{j3, m3}]
    !!where m1 = -m2 - m3
    !! j1 varies from min(abs(j2-j3), abs(m1)) to j2+j3 
    COOP_INT,intent(IN)::twoj2, twoj3, twom2, twom3
    COOP_INT,intent(OUT)::num_j1s
    COOP_REAL, dimension(:):: threejs
    COOP_REAL,intent(OUT)::j1min
    COOP_INT::twom1, twoj1
    COOP_REAL:: j1, m1, j1max,  j2,j3, m2,m3
    COOP_REAL newfac,lmatch
    COOP_REAL :: sum1, oldfac, a1, a2, c1, dv, denom, c1old, &
         x, sumuni, c2, sumfor, x1, x2, x3, sum2, &
         a1s, a2s, y, sumbac, y1, y2, y3, ratio, cnorm, sign1, &
         sign2, thresh
    COOP_INT :: j1cmin, j1cmax, lstep, i, nstep2, num_j1sp1, num_j1sp2, &
         num_j1sp3, index, nlim, n

    ! routine to generate set of 3j-coeffs (j1,j2,j3\\ m1,m2,m3)

    ! by recursion from j1min = max(abs(j2-j3),abs(m1)) 
    !                to j1max = j2+j3
    ! the resulting 3j-coeffs are stored as threejs(j1-j1min+1)

    ! to achieve the numerical stability, the recursion will proceed
    ! simultaneously forwards and backwards, starting from j1min and j1max
    ! respectively.
    !
    ! lmatch is the j1-value at which forward and backward recursion are matched.

    COOP_REAL,parameter:: tiny = 1.d-30, srtiny = 1.0d-15, huge = 1.d30,srhuge=1.d15

    twom1 = - twom2 - twom3
    j2 = twoj2/2.d0
    j3 = twoj3/2.d0
    m2 = twom2/2.d0
    m3 = twom3/2.d0
    m1 = -(m2+m3)
    j1min = max(abs(j2-j3),abs(m1))
    j1max = j2+j3
    if (twoj2.lt. 0 .or. twoj3 .lt. 0 .or. twoj2 .lt. abs(twom2) .or. twoj3 .lt. abs(twom3) .or. mod(twoj2+twom2,2).ne.0 .or. mod(twoj3+twom3,2).ne.0 )then
       num_j1s = 0
       return
    endif
    num_j1s = nint(j1max-j1min)+1
    if(num_j1s .gt. size(threejs))then
       write(*,*) "ThreeJ_array error: array size overflow"
       stop
    endif
    lmatch = 0.d0
    
    if (num_j1s .eq. 1)then
       ! reached if j1 can take only one value, i.e.j1min=j1max
       threejs(1) = evenodd_sign(j2+m2-j3+m3)/sqrt(j1min+j2+j3+1.d0)
       return
    end if

    ! starting forward recursion from j1min taking nstep1 steps
    j1 = j1min
    threejs(1) = srtiny
    sum1 = (2*j1 + 1)*tiny

    lstep = 1

30  lstep = lstep+1
    j1 = j1+1
    oldfac = newfac
    a1 = (j1+j2+j3+1)*(j1-j2+j3)*(j1+j2-j3)*(-j1+j2+j3+1)
    a2 = (j1+m1)*(j1-m1)
    newfac = sqrt(a1*a2)
    if (j1 < 1.0001d0) then
       !if J1 = 1  (J1-1) HAS TO BE FACTORED OUT OF DV, HENCE
       c1 = -(2*j1-1)*j1*(m3-m2)/newfac
    else

       dv = -j2*(j2+1)*m1 + j3*(j3+1)*m1 + j1*(j1-1)*(m3-m2)
       denom = (j1-1)*newfac

       if (lstep > 2) c1old = abs(c1)
       c1 = -(2*j1-1)*dv/denom

    end if

    if (lstep<= 2) then
       ! if j1=j1min+1 the third term in the recursion eqn vanishes, hence
       x = srtiny*c1
       threejs(2) = x
       sum1 = sum1+tiny*(2*j1+1)*c1*c1
       if(lstep==num_j1s) then
          sumuni=sum1
          go to 230
       end if
       goto 30

    end if

    c2 = -j1*oldfac/denom

    ! recursion to the next 3j-coeff x  
    x = c1*threejs(lstep-1) + c2*threejs(lstep-2)
    threejs(lstep) = x
    sumfor = sum1
    sum1 = sum1 + (2*j1+1)*x*x
    if (lstep/=num_j1s) then

       ! see if last unnormalised 3j-coeff exceeds srhuge
       if (abs(x) >= srhuge) then
          ! REACHED if LAST 3J-COEFFICIENT LARGER THAN SRHUGE
          ! SO THAT THE RECURSION SERIES THREEJS(1), ... , THREEJS(LSTEP)
          ! HAS TO BE RESCALED TO PREVENT OVERFLOW
          do i = 1, lstep
             if (abs(threejs(i)) < srtiny) threejs(i)= 0.d0
             threejs(i) = threejs(i)/srhuge
          end do
          sum1 = sum1/huge
          sumfor = sumfor/huge
          x = x/srhuge
       end if
       ! as long as abs(c1) is decreasing, the recursion proceeds towards increasing
       ! 3j-valuse and so is numerically stable. Once an increase of abs(c1) is 
       ! detected, the recursion direction is reversed.
       if (c1old > abs(c1)) goto 30
    end if 

    ! keep three 3j-coeffs around lmatch for comparison with backward recursion
    lmatch = j1-1
    x1 = x
    x2 = threejs(lstep-1)
    x3 = threejs(lstep-2)
    nstep2 = num_j1s-lstep+3

    ! --------------------------------------------------------------------------
    !
    ! starting backward recursion from j1max taking nstep2 stpes, so that
    ! forward and backward recursion overlap at 3 points 
    ! j1 = lmatch-1, lmatch, lmatch+1

    num_j1sp1 = num_j1s+1
    num_j1sp2 = num_j1s+2
    num_j1sp3 = num_j1s+3
    j1 = j1max
    threejs(num_j1s) = srtiny
    sum2 = tiny*(2*j1+1)

    j1 = j1+2
    lstep=1

    do
       lstep = lstep + 1
       j1= j1-1

       oldfac = newfac
       a1s = (j1+j2+j3)*(j1-j2+j3-1)*(j1+j2-j3-1)*(-j1+j2+j3+2)
       a2s = (j1+m1-1)*(j1-m1-1)
       newfac = sqrt(a1s*a2s)

       dv = -j2*(j2+1)*m1 + j3*(j3+1)*m1 +j1*(j1-1)*(m3-m2)

       denom = j1*newfac
       c1 = -(2*j1-1)*dv/denom
       if (lstep <= 2) then

          ! if j2=j2max+1, the third term in the recursion vanishes

          y = srtiny*c1
          threejs(num_j1s-1) = y
          sumbac = sum2
          sum2 = sum2 + tiny*(2*j1-3)*c1*c1

          cycle

       end if

       c2 = -(j1-1)*oldfac/denom

       ! recursion to the next 3j-coeff y
       y = c1*threejs(num_j1sp2-lstep)+c2*threejs(num_j1sp3-lstep)

       if (lstep==nstep2) exit

       threejs(num_j1sp1-lstep) = y
       sumbac = sum2
       sum2 = sum2+(2*j1-3)*y*y

       ! see if last unnormalised 3j-coeff exceeds srhuge
       if (abs(y) >= srhuge) then

          ! reached if 3j-coeff larger than srhuge so that the recursion series
          ! threejs(num_j1s),..., threejs(num_j1s-lstep+1) has to be rescaled to prevent overflow

          do i = 1, lstep
             index=num_j1s-i+1
             if (abs(threejs(index)) < srtiny) threejs(index)=0.d0
             threejs(index) = threejs(index)/srhuge
          end do

          sum2=sum2/huge
          sumbac=sumbac/huge

       end if

    end do

    ! the forward recursion 3j-coeffs x1, x2, x3 are to be matched with the 
    ! corresponding backward recursion vals y1, y2, y3

    y3 = y
    y2 = threejs(num_j1sp2-lstep)
    y1 = threejs(num_j1sp3-lstep)

    ! determine now ratio such that yi=ratio*xi (i=1,2,3) holds with minimal error

    ratio = (x1*y1+x2*y2+x3*y3)/(x1*x1+x2*x2+x3*x3)
    nlim = num_j1s-nstep2+1

    if (abs(ratio) >= 1) then

       threejs(1:nlim) = ratio*threejs(1:nlim) 
       sumuni = ratio*ratio*sumfor + sumbac

    else

       nlim = nlim+1
       ratio = 1/ratio
       do n = nlim, num_j1s
          threejs(n) = ratio*threejs(n)
       end do
       sumuni = sumfor + ratio*ratio*sumbac

    end if
    ! normalise 3j-coeffs

230 cnorm = 1/sqrt(sumuni)

    ! sign convention for last 3j-coeff determines overall phase

    sign1 = sign(1.d0,threejs(num_j1s))
    sign2 = evenodd_sign(j2+m2-j3+m3)
    if (sign1*sign2 <= 0) then
       cnorm = -cnorm
    end if
    if (abs(cnorm) >= 1.d0) then
       threejs(1:num_j1s) = cnorm*threejs(1:num_j1s)
       return
    end if

    thresh = tiny/abs(cnorm)

    do n = 1, num_j1s
       if (abs(threejs(n)) < thresh) threejs(n) = 0.d0
       threejs(n) = cnorm*threejs(n)
    end do
    return 

  contains

    function evenodd_sign(x)
      COOP_REAL::x, evenodd_sign
      if(mod(nint(x),2).eq.0)then
         evenodd_sign = 1.d0
      else
         evenodd_sign = -1.d0
      endif
    end function evenodd_sign

  end subroutine coop_ThreeJ_Array

  function coop_pseudoCl_matrix(l_pseudo, l, lmax, Cl_mask) result(m)
    COOP_INT l, l_pseudo, l2, lmax
    COOP_REAL m, Cl_mask(0:lmax)
    do l2 = abs(l-l_pseudo), min(lmax, abs(l+l_pseudo))
       m = m + cl_mask(l2)*(2.d0*l2+1.d0) * coop_ThreeJ000(l_pseudo, l, l2)**2
    enddo
    m = m*(2.d0*l + 1.d0)/coop_4pi/dble(l)/(l+1.d0)*dble(l_pseudo)*(l_pseudo+1.d0)
  end function coop_pseudoCl_matrix


  !!for full (TT, TE or TB, EB, EE+BB, EE-BB) kernels
  function coop_pseudoCl_matrix_pol(l_pseudo, l, lmax, Cl_mask) result(m)
    COOP_INT l, l_pseudo, l2, lmax, evenodd_sign, lmin
    COOP_REAL m(4), Cl_mask(0:lmax), thrj000, thrj220, maskterm
    lmin = abs(l-l_pseudo)
    m = 0.d0
    evenodd_sign = 1- 2*mod(lmin+l_pseudo+l,2)
    do l2 = lmin, min(lmax, abs(l+l_pseudo))
       thrj000 = coop_ThreeJ000(l_pseudo, l, l2)       
       thrj220 = coop_int3j(l_pseudo, l, l2, 2, -2, 0)
       maskterm = cl_mask(l2)*(2.d0*l2+1.d0)
       m(coop_pseudoCl_kernel_index_TT) = m(coop_pseudoCl_kernel_index_TT) + maskterm*thrj000 ** 2
       if(evenodd_sign .gt. 0) &
            m(coop_pseudoCl_kernel_index_TE) = m(coop_pseudoCl_kernel_index_TE) + maskterm*thrj000*thrj220  
       maskterm = maskterm*thrj220**2
       m(coop_pseudoCl_kernel_index_EB) = m(coop_pseudoCl_kernel_index_EB) + maskterm*evenodd_sign
       m(coop_pseudoCl_kernel_index_EE_plus_BB) = m(coop_pseudoCl_kernel_index_EE_plus_BB) + maskterm
       evenodd_sign = -evenodd_sign
    enddo
    m = m*((2.d0*l + 1.d0)/coop_4pi)/dble(l)/(l+1.d0)*dble(l_pseudo)*(l_pseudo+1.d0)
  end function coop_pseudoCl_matrix_pol


  !!temperature only
  subroutine coop_pseudoCl_get_kernel(lmax_mask, Cl_mask, lmin, lmax, kernel)
    COOP_INT::lmin, lmax, lmax_mask
    COOP_REAL  Cl_mask(0:lmax_mask)
    COOP_REAL::kernel(lmin:lmax, lmin:lmax)
    COOP_INT::l1, l2
    !$omp parallel do private(l1, l2)
    do l1 = lmin, lmax
       do l2 = lmin, lmax
          kernel(l1, l2) = coop_pseudoCl_matrix(l1, l2, lmax_mask, cl_mask)
       enddo
    enddo
    !$omp end parallel do
  end subroutine coop_pseudoCl_get_kernel

  !!temperature only
  subroutine coop_pseudoCl2Cl(lmin, lmax, Cl_pseudo, kernel, Cl, smooth)
    COOP_INT::lmin, lmax, l, iell, low, up
    COOP_REAL Cl_pseudo(lmin:lmax), Cl(lmin:lmax), Dl_pseudo(lmin:lmax)
    COOP_REAL::kernel(lmin:lmax, lmin:lmax)
    logical,optional::smooth
    COOP_INT,dimension(:),allocatable::ells
    COOP_REAL,dimension(:,:),allocatable::reduced_kernel, P, Q
    COOP_INT::num_ells
    COOP_REAL,dimension(:),allocatable::smoothed_Dl, smoothed_Dl_pseudo
    !$omp parallel do
    do l = lmin, lmax
       Dl_pseudo(l) = l*(l+1.d0)*Cl_pseudo(l)
    enddo
    !$omp end parallel do
    if(present(smooth))then
       if(smooth .and. lmax - lmin  .gt. 50)then
          num_ells = coop_nl_range(lmin, lmax)
          allocate(ells(num_ells), reduced_kernel(num_ells, num_ells), smoothed_Dl(num_ells), smoothed_Dl_pseudo(num_ells), P(num_ells, lmin:lmax), Q(lmin:lmax, num_ells))
          call coop_set_ells(ells, lmin, lmax)
          P = 0.d0
          do iell = 1, num_ells
             if(iell .ge. 2)then
                low = ells(iell-1)
             else
                low = lmin
             endif
             if(iell .lt. num_ells)then
                up = ells(iell+1)
             else
                up = lmax
             endif
             do l = low, ells(iell)-1                
                P(iell, l) = exp(-2.*((l - ells(iell))/dble(low - ells(iell)))**2)
             enddo
             P(iell, ells(iell)) = 1.d0
             do l = ells(iell) + 1, up
                P(iell, l) = exp(-2.*((l - ells(iell))/dble(up - ells(iell)))**2)
             enddo
             P(iell, low:up) = P(iell, low:up)/sum(P(iell,low:up))
          enddo
          Q = transpose(P)
          do l = lmin, lmax
             Q(l, :) = Q(l, :)/sum(Q(l,:))
          enddo
          reduced_kernel = matmul(P, matmul(kernel, Q))
          smoothed_Dl_pseudo = matmul(P, Dl_pseudo)
          call coop_fit_template(num_ells, num_ells, smoothed_Dl_pseudo, reduced_kernel, smoothed_Dl) 
          Cl = matmul( Q, smoothed_Dl )
          deallocate(ells, reduced_kernel, smoothed_Dl, smoothed_Dl_pseudo, P, Q)
          goto 200
       endif
    endif
100 continue
    call coop_fit_template(lmax-lmin+1, lmax-lmin+1, Dl_pseudo(lmin:lmax), kernel, Cl)
200 continue
    !$omp parallel do
    do l = lmin, lmax
       Cl(l) = Cl(l)/( l*(l+1.d0) )
    enddo
    !$omp end parallel do
  end subroutine coop_pseudoCl2Cl


  !!TT, TE, TB, EB, EE+BB, EE-BB
  subroutine coop_pseudoCl_get_kernel_pol(lmax_mask, Cl_mask, lmin, lmax, kernel)
    COOP_INT::lmin, lmax, lmax_mask
    COOP_REAL  Cl_mask(0:lmax_mask)
    COOP_REAL::kernel(lmin:lmax, lmin:lmax, 4)
    COOP_INT::l1, l2
    !$omp parallel do private(l1, l2)
    do l1 = lmin, lmax
       do l2 = lmin, lmax
          kernel(l1, l2, :) = coop_pseudoCl_matrix_pol(l1, l2, lmax_mask, cl_mask)
       enddo
    enddo
    !$omp end parallel do
  end subroutine coop_pseudoCl_get_kernel_pol

  !!T, E, B
  subroutine coop_pseudoCl2Cl_pol(lmin, lmax, Cl_pseudo, kernel, Cl, smooth)
    COOP_INT::lmin, lmax, iell, i, low, up
    !!use camb ordering
    !!TT, EE, BB, TE, TB, EB
    COOP_REAL Cl_pseudo(lmin:lmax, 6), Dl_pseudo(lmin:lmax, 6), Cl(lmin:lmax, 6), tmp(lmin:lmax)
    COOP_REAL::kernel(lmin:lmax, lmin:lmax,4)
    COOP_INT::l
    !!for smoothed C_l's
    logical,optional::smooth
    COOP_INT,dimension(:),allocatable::ells
    COOP_REAL,dimension(:,:),allocatable:: P, Q
    COOP_REAL,dimension(:,:,:),allocatable::reduced_kernel
    COOP_INT::num_ells
    COOP_REAL,dimension(:,:),allocatable::smoothed_Dl, smoothed_Dl_pseudo

    !$omp parallel do
    do l = lmin, lmax
       Dl_pseudo(l,:) = Cl_pseudo(l,:)*l*(l+1.d0)
    enddo
    !$omp end parallel do
    if(present(smooth))then
       if(smooth .and. lmax - lmin  .gt. 50)then
          num_ells = coop_nl_range(lmin, lmax)
          allocate(ells(num_ells), reduced_kernel(num_ells, num_ells,4), smoothed_Dl(num_ells,6), smoothed_Dl_pseudo(num_ells,6), P(num_ells, lmin:lmax), Q(lmin:lmax, num_ells))
          call coop_set_ells(ells, lmin, lmax)
          P = 0.d0
          do iell = 1, num_ells
             if(iell .ge. 2)then
                low = ells(iell-1)
             else
                low = lmin
             endif
             if(iell .lt. num_ells)then
                up = ells(iell+1)
             else
                up = lmax
             endif
             do l = low, ells(iell)-1                
                P(iell, l) = exp(-2.*((l - ells(iell))/dble(low - ells(iell)))**2)
             enddo
             P(iell, ells(iell)) = 1.d0
             do l = ells(iell) + 1, up
                P(iell, l) = exp(-2.*((l - ells(iell))/dble(up - ells(iell)))**2)
             enddo
             P(iell, low:up) = P(iell, low:up)/sum(P(iell,low:up))
          enddo
          Q = transpose(P)
          do l = lmin, lmax
             Q(l, :) = Q(l, :)/sum(Q(l,:))
          enddo
          do i = 1, 4
             reduced_kernel(:,:,i) = matmul(P, matmul(kernel(:,:,i), Q))
          enddo
          smoothed_Dl_pseudo = matmul(P, Dl_pseudo)
          !!TT
          call coop_fit_template(num_ells, num_ells, smoothed_Dl_pseudo(:, coop_TEB_index_TT), reduced_kernel(:,:, coop_pseudoCl_kernel_index_TT), smoothed_Dl(:, coop_TEB_index_TT))
          !!TE
          call coop_fit_template(num_ells, num_ells, smoothed_Dl_pseudo(:, coop_TEB_index_TE), reduced_kernel(:,:, coop_pseudoCl_kernel_index_TE), smoothed_Dl(:, coop_TEB_index_TE))
          !!TB
          call coop_fit_template(num_ells, num_ells, smoothed_Dl_pseudo(:, coop_TEB_index_TB), reduced_kernel(:,:, coop_pseudoCl_kernel_index_TB), smoothed_Dl(:, coop_TEB_index_TB))
          !!EB
          call coop_fit_template(num_ells, num_ells, smoothed_Dl_pseudo(:, coop_TEB_index_EB), reduced_kernel(:,:, coop_pseudoCl_kernel_index_EB), smoothed_Dl(:, coop_TEB_index_EB))
          !!EE + BB => saved in EE
          call coop_fit_template(num_ells, num_ells, smoothed_Dl_pseudo(:, coop_TEB_index_EE)+smoothed_Dl_pseudo(:, coop_TEB_index_BB), reduced_kernel(:,:, coop_pseudoCl_kernel_index_EE_plus_BB), smoothed_Dl(:, coop_TEB_index_EE))
          !!EE - BB => saved in BB
          call coop_fit_template(num_ells, num_ells, smoothed_Dl_pseudo(:, coop_TEB_index_EE)-smoothed_Dl_pseudo(:, coop_TEB_index_BB), reduced_kernel(:,:, coop_pseudoCl_kernel_index_EE_minus_BB), smoothed_Dl(:, coop_TEB_index_BB))
          Cl = matmul( Q, smoothed_Dl )
          tmp = Cl(lmin:lmax, coop_TEB_index_BB)
          deallocate(ells, reduced_kernel, smoothed_Dl, smoothed_Dl_pseudo, P, Q)
          goto 200
       endif
    endif
100 continue
    !!TT
    call coop_fit_template(lmax-lmin+1, lmax-lmin+1, Dl_pseudo(lmin:lmax, coop_TEB_index_TT), kernel(lmin:lmax,lmin:lmax, coop_pseudoCl_kernel_index_TT), cl(lmin:lmax, coop_TEB_index_TT))
    !!TE
    call coop_fit_template(lmax-lmin+1, lmax-lmin+1, Dl_pseudo(lmin:lmax,  coop_TEB_index_TE), kernel(lmin:lmax,lmin:lmax, coop_pseudoCl_kernel_index_TE), cl(lmin:lmax, coop_TEB_index_TE))
    !!TB
    call coop_fit_template(lmax-lmin+1, lmax-lmin+1, Dl_pseudo(lmin:lmax, coop_TEB_index_TB), kernel(lmin:lmax,lmin:lmax, coop_pseudoCl_kernel_index_TB), cl(lmin:lmax, coop_TEB_index_TB))
    !!EB
    call coop_fit_template(lmax-lmin+1, lmax-lmin+1, Dl_pseudo(lmin:lmax, coop_TEB_index_EB), kernel(lmin:lmax,lmin:lmax, coop_pseudoCl_kernel_index_EB), cl(lmin:lmax,  coop_TEB_index_EB))
    !!EE + BB
    call coop_fit_template(lmax-lmin+1, lmax-lmin+1, Dl_pseudo(lmin:lmax, coop_TEB_index_EE) + Dl_pseudo(lmin:lmax,   coop_TEB_index_BB), kernel(lmin:lmax,lmin:lmax, coop_pseudoCl_kernel_index_EE_plus_BB), cl(lmin:lmax, coop_TEB_index_EE))
    !!EE - BB
    call coop_fit_template(lmax-lmin+1, lmax-lmin+1, Dl_pseudo(lmin:lmax, coop_TEB_index_EE) - Dl_pseudo(lmin:lmax, coop_TEB_index_BB), kernel(lmin:lmax,lmin:lmax, coop_pseudoCl_kernel_index_EE_minus_BB), tmp(lmin:lmax))
200 continue
    cl(lmin:lmax, coop_TEB_index_BB) = (cl(lmin:lmax, coop_TEB_index_EE) - tmp(lmin:lmax))/2.d0
    cl(lmin:lmax, coop_TEB_index_EE) = (cl(lmin:lmax, coop_TEB_index_EE) + tmp(lmin:lmax))/2.d0
    !$omp parallel do
    do l = lmin, lmax
       Cl(l,:) = Cl(l,:)/(l*(l+1.d0))
    enddo
    !$omp end parallel do
    
  end subroutine coop_pseudoCl2Cl_pol


  subroutine coop_next_l(l)
    COOP_INT::l
    l = l + min(60, 12 + l/30, max(1, l/3))
  end subroutine coop_next_l

  function coop_nl_range(lmin, lmax) result(nl)
    COOP_INT::lmin, lmax
    COOP_INT::l, nl
    l = lmin
    nl = 1
    do while(l .lt. lmax)
       nl = nl + 1
       call coop_next_l(l)
    enddo
  end function coop_nl_range


  subroutine coop_set_ells(ells, lmin, lmax)
    COOP_INT::nl, lmax, l, lmin
    COOP_INT::ells(:)
    l = lmin
    nl = 0
    do while(l .lt. lmax)
       nl = nl + 1
       ells(nl) = l
       call coop_next_l(l)
    enddo
    nl = nl + 1
    ells(nl) = lmax
  end subroutine coop_set_ells


end module coop_special_function_mod


