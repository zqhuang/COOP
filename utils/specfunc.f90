module coop_special_function_mod
  use coop_wrapper_typedef
  implicit none
#include "constants.h"

  private

  integer,parameter::dl  = coop_real_length
  COOP_INT,parameter::special_function_maxl = 5000
  COOP_INT,parameter::lngamma_maxl = special_function_maxl*3+2
  COOP_REAL,parameter::special_function_mineps = 1.d-50


  public:: coop_special_function, log2, sinc, sinhc, lngamma, asinh, acosh, atanh

  type coop_special_function
     COOP_REAL,dimension(:),allocatable::lngamma_int_table, lngamma_hi_table, sqrt_table
   contains
     procedure::init => coop_special_function_initialize
     procedure::bessj => coop_special_function_bessj
     procedure::ln_factorial => coop_special_function_ln_factorial
     procedure::ln_factorial_odd => coop_special_function_ln_factorial_odd
  end type coop_special_function

  interface coop_special_function
     procedure coop_special_function_constructor
  end interface coop_special_function

  interface log2
     module procedure log2_s, log2_v
  end interface log2

  interface sinc
     module procedure sinc_s, sinc_v
  end interface sinc

  interface sinhc
     module procedure sinhc_s, sinhc_v
  end interface sinhc


  interface asinh
     module procedure asinh_s, asinh_v
  end interface asinh

  interface acosh
     module procedure acosh_s, acosh_v
  end interface acosh


  interface atanh
     module procedure atanh_s, atanh_v
  end interface atanh

  interface LnGamma
     module procedure lnGamma_s, LnGamma_v
  end interface LnGamma



contains

  function asinh_s(x)
    COOP_REAL x
    COOP_REAL asinh_s
    if(x.gt. 1.d-6)then
       asinh_s = log(x+sqrt(1.d0+x**2))
    else
       asinh_s = x * (1.d0 - x*x/6.d0)
    endif
  end function asinh_s

  function asinh_v(x)
    COOP_REAL,dimension(:),intent(IN):: x
    COOP_REAL asinh_v(size(x))
    COOP_INT i
    !$omp parallel do
    do i=1, size(x)
       asinh_v(i) = asinh_s(x(i))
    enddo
    !$omp end parallel do
  end function asinh_v

  function acosh_s(x)
    COOP_REAL x
    COOP_REAL acosh_s
    if(x.gt. 1.000001d0)then
       acosh_s = log(x+sqrt((x-1.d0)*(x+1.d0)))
    else
       acosh_s = sqrt(2.d0*(x-1.d0))*(1.d0 - (x-1.d0)/12.d0)
    endif
  end function acosh_s

  function acosh_v(x)
    COOP_REAL,dimension(:),intent(IN):: x
    COOP_REAL acosh_v(size(x))
    COOP_INT i
    !$omp parallel do
    do i=1, size(x)
       acosh_v(i) = acosh_s(x(i))
    enddo
    !$omp end parallel do
  end function acosh_v

  function atanh_s(x)
    COOP_REAL x
    COOP_REAL atanh_s
    atanh_s = 0.5d0*log((1.d0+x)/(1.d0-x))
  end function atanh_s

  function atanh_v(x)
    COOP_REAL,dimension(:),intent(IN):: x
    COOP_REAL atanh_v(size(x))
    atanh_v = 0.5d0*log((1.d0+x)/(1.d0-x))
  end function atanh_v  

  function coop_special_function_constructor() result(this)
    type(coop_special_function):: this
    call this%init()
  end function coop_special_function_constructor


  subroutine coop_special_function_initialize(this)
    class(coop_special_function):: this
    COOP_INT l
    if(.not. allocated(this%sqrt_table))then
       allocate(this%lngamma_int_table(0:lngamma_maxl))
       allocate(this%lngamma_hi_table(0:lngamma_maxl))
       allocate(this%sqrt_table(0:lngamma_maxl))
       this%lngamma_hi_table(0) = 1.e30_dl
       this%lngamma_int_table(0) = 0._dl
       this%lngamma_hi_table(1) = coop_lnpi/2.d0
       this%lngamma_int_table(1) = 0._dl
       do l=2, lngamma_maxl
          this%lngamma_hi_table(l) = this%lngamma_hi_table(l-1) + dlog(l-1.5_dl)
          this%lngamma_int_table(l) = this%lngamma_int_table(l-1) + dlog(l-1._dl)
       enddo
       !$omp parallel do
       do l=0, lngamma_maxl
          this%sqrt_table(l) = dsqrt(dble(l))
       enddo
       !$omp end parallel do
    endif
  end subroutine coop_special_function_initialize


  subroutine continued_fraction_recurrence(an, bn, s)
    COOP_REAL an, bn, s(4)
    s = (/ s(3), s(4), bn*s(3)+an*s(1), bn*s(4)+an*s(2) /)
  end subroutine continued_fraction_recurrence

!!input COOP_INT i
!!return i*(i-1)*(i-2)*...*2*1
  function coop_special_function_ln_factorial(this, i) result(frac)
    class(coop_special_function) this
    COOP_INT i
    COOP_REAL frac
    frac = this%lnGamma_int_table(i+1)
  end function coop_special_function_ln_factorial

!!input COOP_INT i
!!return (2*i-1)*(2*i-3)*...*3*1
  function coop_special_function_ln_factorial_odd(this, i) result(frac)
    class(coop_special_function)::this
    COOP_INT i
    COOP_REAL frac
    frac = this%lnGamma_hi_table(i+1) + coop_ln2 * i - 0.5*coop_lnpi
  end function coop_special_function_ln_factorial_odd

  Function sinc_s(x) 
    COOP_REAL x,sinc_s
    if(abs(x).lt.0.001)then
       sinc_s=1._dl - x*x/6._dl
    else
       sinc_s=dsin(x)/x
    endif
  End Function sinc_s

  function sinc_v(x)
    COOP_REAL,dimension(:),intent(in)::x
    COOP_REAL sinc_v(size(x))
    COOP_INT i
    !$omp parallel do
    do i=1, size(x)
       sinc_v(i) = sinc_s(x(i))
    end do
    !$omp end parallel do
  end function sinc_v

  Function Sinhc_s(x)
    COOP_REAL x,sinhc_s
    if(abs(x).lt.0.001)then
       sinhc_s=1._dl + x*x/6._dl
    else
       sinhc_s=dsinh(x)/x
    end if
  End Function Sinhc_s


  function sinhc_v(x)
    COOP_REAL,dimension(:),intent(in)::x
    COOP_REAL sinhc_v(size(x))
    COOP_INT i
    !$omp parallel do
    do i=1, size(x)
       sinhc_v(i) = sinhc_s(x(i))
    end do
    !$omp end parallel do
  end function sinhc_v


  function LnGamma_s(X) !!ln (\Gamma(x))  
    COOP_REAL LnGamma_s,X
    COOP_INT I
    COOP_REAL C(6),S,R,Y,TEMP
    data C,S/76.18009172947146d0,-86.50532032941677d0, &
         24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
         -.5395239384953d-5,2.5066282746310005d0/
    Y=X
    TEMP=X+5.5D0
    TEMP=(X+0.5D0)*dlog(TEMP)-TEMP
    R=1.000000000190015D0
    DO I=1,6
       Y=Y+1.D0
       R=R+C(I)/Y
    endDO
    LnGamma_s = temp + dlog(S*R/X)
  end function LnGamma_s


  function LnGamma_v(x)
    COOP_REAL,dimension(:),intent(in)::x
    COOP_REAL LnGamma_v(size(x))
    COOP_INT i
    !$omp parallel do
    do i=1, size(x)
       LnGamma_v(i) = LnGamma_s(x(i))
    end do
    !$omp end parallel do
  end function LnGamma_v

  function digamma(x) result(psi)
    COOP_REAL x, psi, y
    if(abs(x-nint(x)).lt. special_function_mineps .and. nint(x).le.0)then
       psi = gamma_product( (/ x, 1-x /), (/ 1, 1 /) )
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
  end function digamma


  function gamma_product(x, p, lnnorm) result(f)
    COOP_REAL f
    COOP_REAL,dimension(:)::x
    COOP_INT,dimension(:)::p
    COOP_INT i, n, m, sgn
    COOP_REAL eps, lnf, pifac
    COOP_REAL,optional::lnnorm
    n = size(x)
    if(n.ne. size(p)) stop "invalid input in gamma_product"
    sgn = 1
    lnf = 0
    do i=1,n
       if(x(i).gt.0.d0)then
          lnf = lnf + lngamma(x(i))*p(i)
       else
          m = nint(x(i))
          if(mod(m,2).eq.0)then
             eps = x(i) - m
          else
             eps = m - x(i)
          endif
          if(abs(eps).lt.special_function_mineps)then
             if(mod(m,2).eq.0)then
                eps = special_function_mineps
             else
                eps = -special_function_mineps
             endif
          endif
          pifac = dsin(coop_pi*eps)/coop_pi
          lnf  = lnf + (lngamma(1.d0-x(i)) + dlog(abs(pifac)))*(-p(i))
          if(pifac.lt.0.d0 .and. mod(p(i),2) .ne. 0) sgn = -sgn
       endif
    enddo
    if(present(lnnorm))then
       f = dexp(lnf+lnnorm)*sgn
    else
       f = dexp(lnf)*sgn
    endif
  end function gamma_product

  function special_function_erfc(x)
    COOP_REAL x,special_function_erfc
    COOP_REAL Z,T
    z=dabs(x)
    t=1./(1.+0.5*z)
    special_function_erfc=T*dexp(-Z*Z-1.26551223+T*(1.00002368+T*(0.37409196+ &
         T*(0.09678418+T*(-0.18628806+T*(0.27886807+T*(-1.13520398+ &
         T*(1.48851587+T*(-0.82215223+T*0.17087277)))))))))
    if(X.LT.0)special_function_erfc=2.d0-special_function_erfc
  end function special_function_erfc

  function special_function_erf(x) 
    COOP_REAL x, special_function_erf
    special_function_erf = 1.d0 - special_function_erfc(x)
  end function special_function_erf

  !! InCompleteGamma(a,x)=\int_x^\infty t^{a-1}e^{-t}dt, for X>=0, A>0 
  !! Numerical Recipies
  function InCompleteGamma(A,X)
    COOP_REAL InCompleteGamma
    COOP_REAL A,X
    if(X.LT.(A+1.D0))then
       InCompleteGamma=exp(LNGamma(A))
       if(X.ne.0.d0)InCompleteGamma=InCompleteGamma-D_Gamma1(A,X)
    ELSE
       InCompleteGamma=D_Gamma2(A,X)
    endif
  end function InCompleteGamma

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




  subroutine sf_CISIA(X,CI,SI)
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
    ELSE if (X.LE.16.0D0) THEN
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
    ELSE if (X.LE.32.0D0) THEN
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
    ELSE
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
  end subroutine sf_CISIA


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


  subroutine SphericalBessJ(L,X,JL)
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
       write(*,*) 'Can not evaluate Spherical Bessel Function with index l<0'
       STOP
    endif
    AX=Dabs(X)
    AX2=AX**2
    if(l.lt.7)then
       select case(l)
       case(0)
          if(AX.LT.1.D-1)THEN
             JL=1.D0-AX2/6.D0*(1.D0-AX2/20.D0)
          ELSE
             JL=DSIN(AX)/AX
          endif

       case(1)
          if(AX.LT.2.D-1)THEN
             JL=AX/3.D0*(1.D0-AX2/10.D0*(1.D0-AX2/28.D0))
          ELSE
             JL=(DSIN(AX)/AX-DCOS(AX))/AX
          endif
       case(2)
          if(AX.LT.3.D-1)THEN
             JL=AX2/15.D0*(1.D0-AX2/14.D0*(1.D0-AX2/36.D0))
          ELSE
             JL=(-3.0D0*DCOS(AX)/AX-DSIN(AX)*(1.D0-3.D0/AX2))/AX
          endif
       case(3)
          if(AX.LT.4.D-1)THEN
             JL=AX*AX2/105.D0*(1.D0-AX2/18.D0*(1.D0-AX2/44.D0))
          ELSE
             JL=(DCOS(AX)*(1.D0-15.D0/AX2)-DSIN(AX)*(6.D0-15.D0/AX2)/AX)/AX
          endif
       case(4)
          if(AX.LT.6.D-1)THEN
             JL=AX2**2/945.D0*(1.D0-AX2/22.D0*(1.D0-AX2/52.D0))
          ELSE
             JL=(DSIN(AX)*(1.D0-(45.D0-105.D0/AX2)/AX2)+DCOS(AX)*(10.D0-105.D0/AX2)/AX)/AX
          endif
       case(5)
          if(AX.LT.1.D0)THEN
             JL=AX2**2*AX/10395.D0*(1.D0-AX2/26.D0*(1.D0-AX2/60.D0))
          ELSE
             JL=(DSIN(AX)*(15.D0-(420.D0-945.D0/AX2)/AX2)/AX-DCOS(AX)*(1.D0-(105.D0-945.0d0/AX2)/AX2))/AX
          endif
       case(6)
          if(AX.LT.1.D0)THEN
             JL=AX2**3/135135.D0*(1.D0-AX2/30.D0*(1.D0-AX2/68.D0))
          ELSE
             JL=(DSIN(AX)*(-1.D0+(210.D0-(4725.D0-10395.D0/AX2)/AX2)/AX2)+ &
                  DCOS(AX)*(-21.D0+(1260.D0-10395.D0/AX2)/AX2)/AX)/AX
          endif
       end select
    ELSE
       NU=0.5D0+L
       NU2=NU**2
       if(AX.LT.1.D-40)THEN
          JL=0.D0
       ELSEif((AX2/L).LT.5.D-1)THEN
          JL=dexp(L*dlog(AX/NU)-LN2+NU*ONEMLN2-(1.D0-(1.D0-3.5D0/NU2)/NU2/30.D0)/12.D0/NU) &
               /NU*(1.D0-AX2/(4.D0*NU+4.D0)*(1.D0-AX2/(8.D0*NU+16.D0)*(1.D0-AX2/(12.D0*NU+36.D0))))
       ELSEif((real(l,dl)**2/AX).LT.5.D-1)THEN
          BETA = AX- coop_pio2 * (L+1)
          JL=(DCOS(BETA)*(1.D0-(NU2-0.25D0)*(NU2-2.25D0)/8.D0/AX2*(1.D0-(NU2-6.25)*(NU2-12.25D0)/48.D0/AX2)) &
               -DSIN(BETA)*(NU2-0.25D0)/2.D0/AX* (1.D0-(NU2-2.25D0)*(NU2-6.25D0)/24.D0/AX2*(1.D0-(NU2-12.25)* &
               (NU2-20.25)/80.D0/AX2)) )/AX   
       ELSE
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

          ELSEif (AX .gt. NU+(1.48+ 10./nu2)*L3) then
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

          ELSE
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
  end subroutine SphericalBessJ

  function spherical_bessel_j_vec(l, xarr) result(jl)
    COOP_REAL,dimension(:)::xarr
    COOP_REAL jl(size(xarr))
    COOP_INT i, l
    !$omp parallel do
    do i=1,size(xarr)
       call SphericalBessJ(l, xarr(i), jl(i))
    end do
    !$omp end parallel do
  end function spherical_bessel_j_vec

  function FT_spherical_tophat(kR) 
    !!Fourier transformation of a tophat function in a sphere with radius R
    COOP_REAL kR, FT_spherical_tophat
    if(kR .gt. 0.02)then
       FT_spherical_tophat = (dsin(kR)/kR -  dcos(kR))/ kR**2 * 3._dl
    else
       FT_spherical_tophat = (10._dl - kR**2*(1._dl - kR**2/28._dl))/10._dl 
    endif
  end function FT_spherical_tophat


  function FT_Gaussian3D_Window(sigma, k, kp) result(w)
    COOP_REAL w, k, kp, sigma
    w = kp/k/(coop_sqrt2*coop_sqrtpi)/sigma * (exp(-((kp-k)/sigma)**2/2.) -  exp(-((kp+k)/sigma)**2/2.))
  end function FT_Gaussian3D_Window

!!J_0(x)
  function bessj0(x)
    COOP_REAL bessj0,x
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
       bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y* &
            (s4+y*(s5+y*s6)))))
    else
       ax=abs(x)
       z=8.d0/ax
       y=z**2
       xx=ax-.785398164d0
       bessj0=sqrt(.636619772d0/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y* &
            p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
    endif
  end function bessj0


!!J_1(x)
  function bessJ1(X)
    COOP_REAL, intent(in) :: x
    COOP_REAL bessj1,ax,z,xx
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
       bessJ1=X*(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))) &
            /(S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))
    else
       ax=abs(x)
       Z=8.0d0/ax
       Y=Z**2
       XX=AX-2.356194491d0
       bessJ1=sqrt(.636619772d0/AX)*(COS(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y &
            *P5))))-Z*SIN(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5))))) &
            *SIGN(1._dl,x)
    endif
  end function bessJ1



!!J_n
  function coop_special_function_bessj(this, n, x) result(bessj)
    class(coop_special_function)::this
    COOP_REAL bessj
    COOP_REAL, intent(in) :: x
    COOP_INT, intent(in) :: n
    COOP_INT, parameter :: iacc = 40
    COOP_REAL, parameter :: bigno=1.d10,bigni=1.d-10
    COOP_INT jsum,j,m
    COOP_REAL bj, bjm, bjp, tox, ssum
    if(n .lt. 2)then
       if(n.eq.1)then
          bessj = bessj1(x)
       elseif(n.eq.0)then
          bessj = bessj0(x)
       else
          stop "wrong argument n in bessj"
       endif
       return
    endif
    if( x**2 .lt. n/2.)then
       bessj = (x/2.d0)**n*exp(-this%lngamma_int_table(n+1))*(1.d0-x**2/(4.d0*(n+1))*(1.d0-x**2/(8.d0*(n+2))*(1.d0-x**2/(12.d0*(n+3))*(1.d0-x**2/(16.d0*(n+4))))))
       return
    endif
    tox = 2/X
    if(x .gt. real(n))then
       BJM=bessJ0(X)
       BJ=bessJ1(X)
       DO J=1,N-1
          BJP=J*tox*BJ-BJM
          BJM=BJ
          BJ=BJP
       end DO
       bessJ=BJ
    ELSE
       M=2*((N+int(sqrt(real(iacc*N))))/2)
       bessJ=0.0d0
       jsum=0
       ssum=0._dl
       BJP=0._dl
       BJ=1._dl
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
  end function Coop_special_function_bessj


  subroutine rec3jj(thrcof,l2,l3,m2,m3)
    !Recursive evaluation of 3j symbols. Does minimal error checking on input parameters.
    implicit none
    COOP_INT, intent(in) :: l2,l3, m2,m3
    COOP_INT:: l1, m1, l1min,l1max
    COOP_REAL newfac,lmatch
    COOP_REAL, dimension(*) :: thrcof
    COOP_INT ier
    COOP_REAL, parameter :: zero=0._dl, eps=0.01_dl, one=1._dl
    COOP_REAL :: srtiny, sum1, tiny, oldfac, a1, a2, c1, dv, denom, c1old, &
         x, sumuni, c2, sumfor, srhuge, huge, x1, x2, x3, sum2, &
         a1s, a2s, y, sumbac, y1, y2, y3, ratio, cnorm, sign1, &
         sign2, thresh
    COOP_INT :: l1cmin, l1cmax, nfin, lstep, i, nstep2, nfinp1, nfinp2, &
         nfinp3, index, nlim, n


    ! routine to generate set of 3j-coeffs (l1,l2,l3\\ m1,m2,m3)

    ! by recursion from l1min = max(abs(l2-l3),abs(m1)) 
    !                to l1max = l2+l3
    ! the resulting 3j-coeffs are stored as thrcof(l1-l1min+1)

    ! to achieve the numerical stability, the recursion will proceed
    ! simultaneously forwards and backwards, starting from l1min and l1max
    ! respectively.
    !
    ! lmatch is the l1-value at which forward and backward recursion are matched.
    !
    ! ndim is the length of the array thrcof
    !
    ! ier = -1 for all 3j vanish(l2-abs(m2)<0, l3-abs(m3)<0 or not COOP_INT)
    ! ier = -2 if possible 3j's exceed ndim
    ! ier >= 0 otherwise
    !
    data tiny,srtiny /1.0d-30,1.0d-15/
    data huge,srhuge /1.0d30,1.0d15/

    lmatch = zero
    m1 = -(m2+m3)

    ! check relative magnitude of l and m values
    ier = 0

    if (l2 < abs(m2) .or. l3 < m3) then
       ier = -1
       return
    end if

    ! limits for l1
    l1min = max(abs(l2-l3),abs(m1))
    l1max = l2+l3

    if (l1min >= l1max-eps) then
       if (l1min/=l1max) then
          ier = -1
          return
       end if

       ! reached if l1 can take only one value, i.e.l1min=l1max
       thrcof(1) = (-1)**abs(l2+m2-l3+m3)/sqrt(l1min+l2+l3+one)
       l1cmin = l1min
       l1cmax = l1max
       return

    end if

    nfin = int(l1max-l1min+one)


    ! starting forward recursion from l1min taking nstep1 steps
    l1 = l1min
    thrcof(1) = srtiny
    sum1 = (2*l1 + 1)*tiny

    lstep = 1

30  lstep = lstep+1
    l1 = l1+1

    oldfac = newfac
    a1 = (l1+l2+l3+1)*(l1-l2+l3)*(l1+l2-l3)*(-l1+l2+l3+1)
    a2 = (l1+m1)*(l1-m1)
    newfac = sqrt(a1*a2)
    if (l1 < one+eps) then
       !if L1 = 1  (L1-1) HAS TO BE FACTORED OUT OF DV, HENCE
       c1 = -(2*l1-1)*l1*(m3-m2)/newfac
    else

       dv = -l2*(l2+1)*m1 + l3*(l3+1)*m1 + l1*(l1-1)*(m3-m2)
       denom = (l1-1)*newfac

       if (lstep > 2) c1old = abs(c1)
       c1 = -(2*l1-1)*dv/denom

    end if

    if (lstep<= 2) then

       ! if l1=l1min+1 the third term in the recursion eqn vanishes, hence
       x = srtiny*c1
       thrcof(2) = x
       sum1 = sum1+tiny*(2*l1+1)*c1*c1
       if(lstep==nfin) then
          sumuni=sum1
          go to 230
       end if
       goto 30

    end if

    c2 = -l1*oldfac/denom

    ! recursion to the next 3j-coeff x  
    x = c1*thrcof(lstep-1) + c2*thrcof(lstep-2)
    thrcof(lstep) = x
    sumfor = sum1
    sum1 = sum1 + (2*l1+1)*x*x
    if (lstep/=nfin) then

       ! see if last unnormalised 3j-coeff exceeds srhuge
       if (abs(x) >= srhuge) then

          ! REACHED if LAST 3J-COEFFICIENT LARGER THAN SRHUGE
          ! SO THAT THE RECURSION SERIES THRCOF(1), ... , THRCOF(LSTEP)
          ! HAS TO BE RESCALED TO PREVENT OVERFLOW

          ier = ier+1
          do i = 1, lstep
             if (abs(thrcof(i)) < srtiny) thrcof(i)= zero
             thrcof(i) = thrcof(i)/srhuge
          end do

          sum1 = sum1/huge
          sumfor = sumfor/huge
          x = x/srhuge

       end if

       ! as long as abs(c1) is decreasing, the recursion proceeds towards increasing
       ! 3j-valuse and so is numerically stable. Once an increase of abs(c1) is 
       ! detected, the recursion direction is reversed.

       if (c1old > abs(c1)) goto 30

    end if !lstep/=nfin

    ! keep three 3j-coeffs around lmatch for comparison with backward recursion

    lmatch = l1-1
    x1 = x
    x2 = thrcof(lstep-1)
    x3 = thrcof(lstep-2)
    nstep2 = nfin-lstep+3

    ! --------------------------------------------------------------------------
    !
    ! starting backward recursion from l1max taking nstep2 stpes, so that
    ! forward and backward recursion overlap at 3 points 
    ! l1 = lmatch-1, lmatch, lmatch+1

    nfinp1 = nfin+1
    nfinp2 = nfin+2
    nfinp3 = nfin+3
    l1 = l1max
    thrcof(nfin) = srtiny
    sum2 = tiny*(2*l1+1)

    l1 = l1+2
    lstep=1

    do
       lstep = lstep + 1
       l1= l1-1

       oldfac = newfac
       a1s = (l1+l2+l3)*(l1-l2+l3-1)*(l1+l2-l3-1)*(-l1+l2+l3+2)
       a2s = (l1+m1-1)*(l1-m1-1)
       newfac = sqrt(a1s*a2s)

       dv = -l2*(l2+1)*m1 + l3*(l3+1)*m1 +l1*(l1-1)*(m3-m2)

       denom = l1*newfac
       c1 = -(2*l1-1)*dv/denom
       if (lstep <= 2) then

          ! if l2=l2max+1, the third term in the recursion vanishes

          y = srtiny*c1
          thrcof(nfin-1) = y
          sumbac = sum2
          sum2 = sum2 + tiny*(2*l1-3)*c1*c1

          cycle

       end if

       c2 = -(l1-1)*oldfac/denom

       ! recursion to the next 3j-coeff y
       y = c1*thrcof(nfinp2-lstep)+c2*thrcof(nfinp3-lstep)

       if (lstep==nstep2) exit

       thrcof(nfinp1-lstep) = y
       sumbac = sum2
       sum2 = sum2+(2*l1-3)*y*y

       ! see if last unnormalised 3j-coeff exceeds srhuge
       if (abs(y) >= srhuge) then

          ! reached if 3j-coeff larger than srhuge so that the recursion series
          ! thrcof(nfin),..., thrcof(nfin-lstep+1) has to be rescaled to prevent overflow

          ier=ier+1
          do i = 1, lstep
             index=nfin-i+1
             if (abs(thrcof(index)) < srtiny) thrcof(index)=zero
             thrcof(index) = thrcof(index)/srhuge
          end do

          sum2=sum2/huge
          sumbac=sumbac/huge

       end if

    end do

    ! the forward recursion 3j-coeffs x1, x2, x3 are to be matched with the 
    ! corresponding backward recursion vals y1, y2, y3

    y3 = y
    y2 = thrcof(nfinp2-lstep)
    y1 = thrcof(nfinp3-lstep)

    ! determine now ratio such that yi=ratio*xi (i=1,2,3) holds with minimal error

    ratio = (x1*y1+x2*y2+x3*y3)/(x1*x1+x2*x2+x3*x3)
    nlim = nfin-nstep2+1

    if (abs(ratio) >= 1) then

       thrcof(1:nlim) = ratio*thrcof(1:nlim) 
       sumuni = ratio*ratio*sumfor + sumbac

    else

       nlim = nlim+1
       ratio = 1/ratio
       do n = nlim, nfin
          thrcof(n) = ratio*thrcof(n)
       end do
       sumuni = sumfor + ratio*ratio*sumbac

    end if
    ! normalise 3j-coeffs

230 cnorm = 1/sqrt(sumuni)

    ! sign convention for last 3j-coeff determines overall phase

    sign1 = sign(one,thrcof(nfin))
    sign2 = (-1)**abs(l2+m2-l3+m3)
    if (sign1*sign2 <= 0) then
       cnorm = -cnorm
    end if
    if (abs(cnorm) >= one) then
       thrcof(1:nfin) = cnorm*thrcof(1:nfin)
       return
    end if

    thresh = tiny/abs(cnorm)

    do n = 1, nfin
       if (abs(thrcof(n)) < thresh) thrcof(n) = zero
       thrcof(n) = cnorm*thrcof(n)
    end do
    return 

  end subroutine rec3jj


  ! ----------------------------------------------------------------------
  ! Auxiliary Bessel functions for N=0, N=1
  function bessI0(X)
    double precision X,bessI0,Y,P1,P2,P3,P4,P5,P6,P7,  &
         Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
    DATA P1,P2,P3,P4,P5,P6,P7/1.D0,3.5156229D0,3.0899424D0,1.2067429D0,  &
         0.2659732D0,0.360768D-1,0.45813D-2/
    DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1, &
         0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,  &
         0.2635537D-1,-0.1647633D-1,0.392377D-2/
    if(abs(X).LT.3.75D0) THEN
       Y=(X/3.75D0)**2
       bessI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
    ELSE
       AX=abs(X)
       Y=3.75D0/AX
       BX=EXP(AX)/sqrt(AX)
       AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
       bessI0=AX*BX
    endif
    return
  end function bessI0
  ! ----------------------------------------------------------------------
  function bessI1(X)
    double precision X,bessI1,Y,P1,P2,P3,P4,P5,P6,P7,  &
         Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
    DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,  &
         0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
    DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1, &
         -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1, &
         -0.2895312D-1,0.1787654D-1,-0.420059D-2/
    if(abs(X).LT.3.75D0) THEN
       Y=(X/3.75D0)**2
       bessI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
    ELSE
       AX=abs(X)
       Y=3.75D0/AX
       BX=EXP(AX)/sqrt(AX)
       AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
       bessI1=AX*BX
    endif
    return
  end function bessI1


  function bessI(N,X)
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
    double precision X,bessI,tox,biM,bi,biP
    if(n.lt.2)then
       if (N.eq.0) then
          bessI = bessI0(X)
          return
       endif
       if (N.eq.1) then
          bessI = bessI1(X)
          return
       endif
    endif
    if(abs(x) .lt. 1.d-20) then
       bessI=0.d0
       return
    endif
    tox = 2.D0/X
    bip = 0.D0
    bi = 1.D0
    bessI = 0.D0
    M = 2*((N+int(sqrt(real(iacc*N)))))
    DO J = M,1,-1
       bim = bip+ J*tox*bi
       bip = bi
       bi  = bim
       if (abs(bi).gt.bigno) then
          bi  = bi*bigni
          biP = biP*bigni
          bessI = bessI*bigni
       endif
       if (J.eq.N) bessI = bip
    end DO
    bessI = bessI*bessI0(X)/bi
    return
  end function bessI


  !!return Legendre polynomial P_l(x)
  function Legendre_P(l, x) result(Pl)
    COOP_INT l, ell
    COOP_REAL x, Pl, x2, Plm1, Plm2
    select case(l)
    case(0)
       Pl = 1._dl
       return
    case(1)
       Pl = x
       return
    case(2)
       Pl = 1.5_dl*x*x - 0.5_dl
       return
    case(3)
       Pl = (2.5_dl * x * x - 1.5_dl) * x
       return
    case(4)
       x2 = x*x
       Pl = 3._dl/8. + x2 * ( -15._dl/4._dl + x2 * (35._dl/8._dl)  )
       return
    case(5)
       x2 = x*x
       Pl = x * ( 15._dl/8. + x2 * ( -35._dl/4._dl + x2 * (63._dl/8.)))
       return
    case(6)
       x2 = x*x
       Pl = -5._dl/16. + x2*( 105._dl/16. + x2*( -315._dl/16.+ x2*(231._dl/16.)))
       return
    case(7)
       x2 = x*x
       Pl = x * ( -35._dl/16. + x2 * ( 315._dl/16._dl + x2 *( -693._dl/16. +  x2 *(429._dl/16.))))
       return
    case(8)
       x2 = x*x
       Pl = 35._dl/128. + x2 * ( -315._dl/32._dl + x2 * (3465._dl/64._dl + x2 *(-3003._dl/32. + x2 * (6435._dl/128.))))
       return
    case(9)
       x2 = x*x
       Pl = x * (315._dl/128. + x2 * ( -1155._dl/32. + x2 * ( 9009._dl/64. + x2 * ( -6435._dl/32. +  x2* (12155._dl/128.)))))
       return
    case(10)
       x2 = x * x
       Pl = -63._dl/256. + x2 * (3465._dl/256. + x2 * (-15015._dl/128. + x2 * (45045._dl/128. + x2 * ( -109395._dl/256.  + x2 * (46189._dl/256.)))))
       return
    case(11)
       x2 = x*x
       Pl = x * (-693._dl/256. +  x2 * ( 15015._dl/256. + x2 * (-45045._dl/128.+ x2 * ( 109395._dl/128. + x2 * (-230945._dl/256. + x2 * (88179._dl/256.))))))
       return
    case(12)
       x2 = x*x
       Pl = 231._dl/1024. + x2 * ( -9009._dl/512.+ x2 * (225225._dl/1024.+ x2 * (-255255._dl/256 + x2 * (2078505._dl/1024 + x2 * (-969969._dl/512. + x2 * (676039._dl/1024.))))))
       return
    case(13)
       x2 = x*x
       Pl = x*(3003._dl /1024 + x2*(- 45045._dl/512. + x2*(765765._dl/1024 + x2*(-692835._dl/256. + x2*(4849845._dl/1024 + x2*(-2028117._dl/512 + x2*(1300075._dl/1024.)))))))
       return
    case(14)
       x2 = x*x
       Pl = -429._dl/2048 + x2*(45045._dl/2048 + x2*(-765765._dl/2048 + x2* (4849845._dl/2048 + x2*(-14549535._dl/2048. + x2*(22309287._dl/2048. + x2*(-16900975._dl/2048. + x2*(5014575._dl/2048.)))))))
       return
    case(15)
       x2 = x*x
       Pl=x*(-6435._dl/2048. + x2*(255255._dl/2048 + x2*(-2909907._dl/2048. + x2*(14549535._dl/2048. + x2*(-37182145._dl/2048. + x2*(50702925._dl/2048 + x2*(-35102025._dl/2048. + x2*(9694845._dl/2048.))))))))
    case(16)
       x2 = x*x
       Pl = 6435._dl/32768. + x2*(-109395._dl/4096. + x2*(4849845._dl/8192. +x2*(-20369349._dl/4096. + x2*(334639305._dl/16384. + x2*(-185910725._dl/4096. + x2*(456326325._dl/8192. + x2*(-145422675._dl/4096. + x2*(300540195._dl/32768.))))))))
    case(17)
       x2=x*x
       Pl = x*(109395._dl/32768. + x2*(-692835._dl/4096. + x2*(20369349._dl/8192. + x2*(-66927861._dl/4096. + x2*(929553625._dl/16384. + x2*(-456326325._dl/4096. + x2*(1017958725._dl/8192. + x2*(-300540195._dl/4096. + x2*(583401555._dl/32768.)))))))))
    case(18)
      x2=x*x
      Pl = -12155._dl/65536. + x2*(2078505._dl/65536. + x2*(-14549535._dl/16384. + x2*(156165009._dl/16384.+ x2*( -1673196525._dl/32768 + x2*(5019589575._dl/32768. + x2*(-4411154475._dl/16384. + x2*(4508102925._dl/16384. + x2*(-9917826435._dl/65536. + x2*(2268783825._dl/65536.)))))))))
    case default
       Plm1 = 1._dl
       Pl = x
       do ell = 2, l
          Plm2 = Plm1
          Plm1 = Pl
          Pl = ((2*ell -1)*x*Plm1 - (ell -1)*Plm2)/ell
       enddo
       return
    end select
  end function Legendre_P

  subroutine sphere_correlation_init(lmax, als, bls)
    COOP_INT lmax, l
    COOP_REAL als(0:lmax), bls(0:lmax)
    als(0) = 1.d0/coop_4pi
    bls(0) = 0.d0
    do l = 1, lmax
       als(l) = (2*l+1.d0)/l
       bls(l) = -(l-1.d0)/(2*l-3.d0)*(2*l+1.d0)/l
    enddo
  end subroutine sphere_correlation_init

  function Sphere_Correlation(lmax, Cls, als, bls, x) result(s) !!x = cos(theta)  !!I assume lmax>=1
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
  end function Sphere_Correlation

  !! return sqrt(4 pi / (2l + 1) ) * Y_l^m (arccos x, 0)
  recursive function Normalized_Associated_Legendre(l, m, x) result(Plm)
    COOP_INT l, m
    COOP_REAL x, x2, Plm
    select case(m)
    case(0)
       Plm = Legendre_P(l, x)
       return
    case(1)
       select case(l)
       case(0)
          Plm = 0._dl
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
          Plm = sqrt(7.5*(1.-x2))*(-1._dl/8. + x2 * (14._dl/8._dl  + x2 * (-21._dl/8.)))
       case(6)
          x2 = x*x
          Plm = sqrt(10.5*(1.-x2)) * x * ( -5._dl/8._dl + x2 * (30._dl/8._dl + x2 * (-33._dl/8._dl)))
       case(7)
          x2 = x*x
          Plm = sqrt(3.5*(1.-x2)) * (5./32._dl + x2*(-135._dl/32. + x2*(495./32._dl + x2 * (-429./32._dl))))
       case(8)
          x2 = x*x
          Plm = sqrt(2.*(1.-x2))*x*(105._dl/64. + x2 * (-1155._dl/64. + x2 * (3003._dl/64.+ x2* (-2145._dl/64.))))
       case(9)
          x2 = x*x
          Plm = sqrt(2.5*(1.-x2))*(-21./128._dl + x2 * (231._dl/32. + x2 * (-3003._dl/64. + x2 * (3003._dl/32. + x2 * (-7293._dl/128.)))))
       case(10)
          x2 = x*x
          Plm = sqrt(27.5*(1.-x2)) * x * (-63./128._dl + x2 *(273._dl/32. + x2 * (-2457._dl/64. + x2 *(1989._dl/32. + x2 * (-4199._dl/128.)))))
       case(11)
          x2 = x*x
          Plm = sqrt(33._dl*(1.-x2))*(21./512._dl + x2*(-1365._dl/512. + x2*(6825._dl/256. + x2 *(-23205._dl/256. + x2 * (62985._dl/512. + x2 * (-29393._dl/512.))))))
       case(12)
          x2 = x*x
          Plm = sqrt(39._dl * (1.-x2))*x*(231./512._dl + x2 *(-5775./512._dl + x2 * ( 19635._dl/256.+ x2 * (-53295._dl/256. + x2 * (124355._dl/512. + x2 * (-52003._dl/512.))))))
       case(13)
          x2 = x*x
          Plm = sqrt(45.5*(1.-x2)) * ( -33._dl/1024. + x2 * (1485._dl/512. + x2 * (-42075._dl/1024. + x2 * (53295._dl/256. + x2 * (-479655._dl/1024. + x2*(245157._dl/512. + x2 * (-185725._dl/1024.)))))))
       case(14)
          x2 = x*x
          Plm = sqrt(52.5*(1.-x2)) * x * (-429._dl/1024. + x2 * (7293._dl/512. + x2 * (-138567._dl/1024. + x2 * (138567._dl/256. + x2 * (-1062347._dl/1024. + x2 * (482885._dl/512. + x2 * (-334305._dl/1024.)))))))
       case(15)
          x2 = x*x
          Plm = sqrt(15.*(1.-x2))*(429._dl/8192.+x2*(-51051._dl/8192. + x2*(969969._dl/8192. + x2*(-6789783._dl/8192 + x2*(22309287._dl/8192. + x2*(-37182145._dl/8192. + x2*(30421755._dl/8192. + x2*(-9694845._dl/8192))))))))
       case(16)
          x2 = x*x
          Plm = sqrt(17.*(1.-x2))*x*(6435._dl/8192.+x2*(-285285._dl/8192. + x2*(3594591._dl/8192.+x2*( -19684665._dl/8192. + x2*(54679625._dl/8192. + x2*(-80528175._dl/8192. + x2*(59879925._dl/8192. +x2*(-17678835._dl/8192.))))))))
       case(17)
          x2 = x*x
          Plm = sqrt(8.5_dl*(1.-x2))*(-2145._dl/32768._dl + x2*(40755._dl/4096._dl +x2*(-1996995._dl/8192._dl + x2*(9186177._dl/4096. +x2*(-164038875._dl/16384._dl + x2*(98423325._dl/4096.+x2*( -259479675._dl/8192. + x2*(88394175._dl/4096.+x2*( -194467185._dl/32768.)))))))))
       case(18)
          x2 = x*x
          Plm  = sqrt(9.5*(1.-x2))*x*(-36465._dl/32768. + x2*(255255._dl/4096.+x2*(-8219211._dl/8192. + x2*(29354325._dl/4096. + x2*(-440314875._dl/16384. + x2*(232166025._dl/4096. +x2*(- 553626675._dl/8192. + x2* (173996955._dl/4096. + x2*(-358229025._dl/32768.)))))))))
       case default
          goto 100
       end select
       return
    case(2)
       select case(l)
       case(0:1)
          Plm = 0._dl
          return
       case(2)
          Plm = (0.5_dl*sqrt(3./2._dl))* ( 1._dl-x*x)
       case(3)
          Plm = (0.5_dl * sqrt(15./2._dl))*x*(1._dl -x *x)
       case(4)
          Plm = (-0.25_dl * sqrt(2.5_dl))*(1.-x**2)*(1.- 7.*x**2)
       case(5)
          Plm = (-0.25_dl * sqrt(105._dl/2._dl)) * x * (1.- x**2) * (1. - 3.*x**2)
       case(6)
          x2 = x*x
          Plm = (sqrt(105._dl) /32._dl)*(1.-x2) * (1. + x2 *(-18._dl + x2 * 33._dl))
       case(7)
          x2 = x * x
          Plm = (sqrt(21._dl)/32.)*x*(1. - x2)*(15. + x2*(-110. + x2 * 143.))
       case(8)
          x2 = x*x
          Plm = (3._dl/64.*sqrt(35._dl))*(1.-x2)*(-1.+ x2*(33. + x2 *143.* (-1.+ x2)))
       case(9)
          x2 = x*x
          Plm = (3._dl/64. * sqrt(55._dl))*x*(1.-x2)*(-7.+x2*(91.+x2*(-273.+x2*221.)))
       case(10)
          x2 = x*x
          Plm = (sqrt(165._dl/2.)/256.)*(1-x2)*(7. + x2 * 13. * (-28. + x2*(210. + x2 * (-476. + x2*323.))))
       case(11)
          x2 = x*x
          Plm = (sqrt(2145._dl/2.)/256.)* x*(1.-x2)*(21. + x2 * (-420. + x2 * (2142. + x2 * (-3876. + x2 * 2261.))))
       case(12)
          x2 = x*x
          Plm = (sqrt(3003._dl/2.)/512.)*(1.-x2)*(-3.+x2*(225. + 17.*x2 * (-150.+ x2*(570. + x2*(-855. + x2 * 437.)))))
       case(13)
          x2 = x*x
          Plm = (sqrt(455._dl/2.)/512.) * x * (1.- x2)*(-99.+x2*17.*(165. + x2*19.*(-66. + x2*(198. + x2*(-253. + x2*115.)))))
       case(14)   
          x2 = x*x
          Plm = (sqrt(1365._dl/2.)/4096._dl) * (1.-x2)*(33. + x2 * 17.*(-198.+ x2*19.*(165.+x2*(-924.+x2*23.*(99. + x2*(-110.+x2*45.))))))
       case(15)
          x2 = x*x
          Plm = (sqrt(1785._dl/2.)/4096._dl)*x*(1.-x2)*(429._dl + x2*(-16302._dl + x2*(171171._dl + x2*(-749892._dl + x2*(1562275._dl + x2*(-1533870._dl + 570285._dl*x2))))))
       case(16)
          x2 = x*x
          Plm = (sqrt(255._dl/2.))*(1.-x2)*(-143/8192._dl + x2*(19019._dl/8192._dl + x2*( -399399._dl/8192._dl + x2*(3062059._dl/8192._dl + x2*(-10935925._dl/8192._dl + x2*(19684665._dl/8192._dl+ x2*(-17298645._dl/8192._dl + x2*(5892945._dl/8192._dl))))))))
       case(17)
          x2 = x*x
          Plm = (sqrt(323._dl/2.))*(1.-x2)*x*(-2145._dl/8192 + x2*(105105._dl/8192. + x2*( -1450449._dl/8192. + x2*(8633625._dl/8192. + x2*(-25900875._dl/8192. + x2*(40970475._dl/8192. + x2*(-32566275._dl/8192. + x2*(10235115._dl/8192.))))))))
       case(18)
          x2 = x*x
          Plm = (sqrt(1615._dl/2.))*(1.-x2)*(429._dl/65536._dl + x2*(-9009._dl/8192._dl + x2*(483483._dl/16384._dl + x2*(-2417415._dl/8192. + x2*(46621575._dl/32768._dl + x2*(-30045015._dl/8192. + x2*(84672315._dl/16384 + x2*( -30705345._dl/8192. + x2*(71645805._dl/65536._dl)))))))))
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
          Plm = 0._dl
       endif
       return
    case(:-1)
       Plm = Normalized_Associated_Legendre(l, -m, x)
       if(mod(m,2).ne.0)Plm = -Plm
       return
    end select
    return
100 call get_normalized_Plm(l,m,x,Plm)    
  End function Normalized_Associated_Legendre


  subroutine get_normalized_Plm(l,m,x,Plm)
    COOP_REAL x, Plm_prev, Plm, Plm_next
    COOP_INT l,m, i
    if(m.gt.l)then
       Plm = 0.d0
       return
    endif
    Plm_prev = (4.d0*(1.d0-x**2))**(m/2.d0)*dexp(lnGamma(m+0.5d0)-lnGamma(2*m+1.d0)/2.d0)/coop_sqrtpi
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
  end subroutine get_normalized_Plm

!! return sqrt((l-m)!/(l+m)!) P_l^m(x)
!!no check
  


!!return Y_l^m(\theta, \phi) 
!!input l>=m>=0, theta, phi
  Function SphericalHarmonicY(l,m,theta, phi) result(Ylm)
    COOP_REAL theta,phi
    COOP_COMPLEX Ylm
    COOP_INT l, m
    Ylm = dsqrt((2*l+1)/coop_4pi) * Normalized_Associated_Legendre(l, m, dcos(theta))*cmplx(dcos(m*phi), dsin(m*phi))
  End Function SphericalHarmonicY

  function Legendre_P_Approx(l, x) result(Pl)
    COOP_INT l
    COOP_REAL x, Pl, theta, phi, sx
    theta = dacos(x)
    sx = dsin(theta)
    phi = (l+0.5d0)*theta-coop_pio4
    Pl = dsqrt((2.d0/coop_pi)/l/sx)*( &
         (1.d0 - 0.25d0/l)*dcos(phi) &
         + 0.125d0*x/sx/l * dsin(phi) &
         )
  end function Legendre_P_Approx

  subroutine get_Cn_array(n, Cn)
    COOP_INT n, i
    COOP_REAL Cn(0:n)
    Cn(0) = 1.d0
    do i = 1, n
       Cn(i) = (n-i+1.d0) * Cn(i-1)/i
    enddo
  end subroutine get_Cn_array

  function factor_lm_pm(l, m)  result (f)
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
       f = exp((lngamma(l+m+1.d0)-lngamma(l-m+1.d0))/2.d0)
    end select
  end function factor_lm_pm

  function sqrtfloor(n) result(intsqrt)
    COOP_INT n, intsqrt
    intsqrt = floor(sqrt(dble(n)+1.e-3_dl))
    do while(intsqrt*intsqrt.gt.n)
       intsqrt = intsqrt - 1 !!just to avoid roundoff errors
    enddo
  end function sqrtfloor

  function sqrtceiling(n) result(intsqrt)
    COOP_INT n, intsqrt
    intsqrt = ceiling(sqrt(dble(n)-1.e-3_dl))
    do while(intsqrt*intsqrt.lt.n)
       intsqrt = intsqrt + 1 !!just to avoid roundoff errors
    enddo
  end function sqrtceiling

  function numbits(m) !!number of bits in binary representation
    COOP_INT n, m, numbits
    numbits = 0
    n = m
    do while(n.gt.0)
       n=n/2
       numbits = numbits+1
    enddo
  end function numbits

  function ellPower(l, x, r)
    COOP_INT l
    COOP_REAL x, r, ellPower
    if(x .le. r)then
       ellPower = (x/r)**l
    else
       ellPower = (r/x)**(l+1)
    endif
  end function ellPower


  !!the integral \int j_l(kD) j_l(kx) k^{-m} k^2 dk (m>=0)
  !!can be approximated with (when m=0 it is exactly equal to) 
  !!   N(l, m) delta(x-D) D^(m-2)
  function jljlintN(l, m) result(Nlm)
    !!if Limber approximation works, returns N(l,m)
    !!if Limber approximation fails, return -1.1e30
    COOP_INT l, m
    COOP_REAL Nlm
    COOP_REAL,parameter::bad_value = -1.1e30_dl
    select case(m)
    case(0)
       Nlm = coop_pio2
    case(1)
       select case(l)
       case(0)
          Nlm = bad_value
       case(1)
          Nlm =  coop_pi2/8.d0
       case(2)
          Nlm = 2.d0/3.d0
       case(3)
          Nlm = coop_pi2*(3.d0/64.d0)
       case(4)
          Nlm = 16.d0/45.d0
       case(5)
          Nlm = coop_pi2*15.d0/512.d0
       case(6)
          Nlm = 128.d0/525.d0
       case default
          Nlm = coop_pio2/(l+0.5d0)/(1.+ 0.375d0/(l+0.5d0)**2)
       end select
    case(2)
       if(l.gt.0)then
          Nlm = coop_pio2/l/(l+1.d0)
       else
          Nlm =bad_value
       endif
    case(3)
       select case(l)
       case(0)
          Nlm = bad_value
       case(1)
          Nlm = (coop_pi2/16.d0)
       case(2)
          Nlm = (1.d0/9.d0)
       case(3)
          Nlm = coop_pi2/256.d0
       case(4)
          Nlm = 4.d0/225.d0
       case(5)
          Nlm = coop_pi2/1024.d0
       case(6)
          Nlm = 64.d0/11025.d0
       case default
          Nlm = coop_pio2/(l-0.5d0)/(l+0.5d0)/(l+1.5d0)*(1.+ 0.37d0/l/l)
       end select
    case(4)
       if(l.gt.1)then
          Nlm = coop_pio2/l/(l+1.d0)/(l+2.d0)/(l-1.d0)
       else
          Nlm=bad_value
       endif
    case default
       if(2*l.gt.m-2)then
          Nlm = coop_pio2*dexp(lnGamma(l-m/2.d0+1.d0)-lnGamma(l+m/2.d0+1.d0))
       else
          Nlm = bad_value
       endif
    end select
  end function jljlintN


  subroutine jljlint_setw(l, m, d, x, w, n)
    COOP_INT l,m,n, nmid, i
    COOP_REAL d, x(n), w(n), lnxbydmin, lnxbydmax,x1,x2, fact1, fact2
    select case(m)
    case(0)
       x(1) = d
       w(1) = coop_pio2/d**2
       n = 1
    case(2)
       if(mod(n,2).eq.0) n = n - 1
       if(n.eq.1)then
          x(1) = d
          w(1) = coop_pio2/l/(l+1.d0)
          return
       endif
       nmid = (n+1)/2 
       lnxbydmin = - 8.d0/l
       lnxbydmax =  8.d0/(l+1)
       call coop_set_uniform(nmid, w(1:nmid), exp(lnxbydmin*(l+1)), 1.d0)
       call coop_set_uniform(nmid, w(nmid:n), 1.d0, exp(-lnxbydmax*l))
       x(1:nmid) = log(w(1:nmid))/(l+1) 
       x(nmid+1:n) = - log(w(nmid+1:n))/l
       x(1:n) = exp(x(1:n))*d
       w(1:nmid-1) = coop_pio4/((l+0.5d0)*(l+1))/(nmid-0.5d0)
       w(nmid+1:n) = coop_pio4/((l+0.5d0)*l)/(nmid-0.5d0)
       w(nmid) = (w(1)+w(n))/2.d0
    case(4)
       if(mod(n,2).eq.0) n = n - 1
       if(n.eq.1)then
          if(l.le.2) stop "divergent integral"
          x(1) = d
          w(1) = d**2*coop_pio2/((l-2.d0)*l*(l+1.d0)*(l+3.d0))
          return
       endif
       nmid = (n+1)/2
       lnxbydmin = -8.d0/l
       lnxbydmax = 8.d0/(l-1)
       call coop_set_uniform(nmid, w(1:nmid), exp(lnxbydmin*(l+1)), 1.d0)
       call coop_set_uniform(nmid, w(nmid:n), 1.d0, exp(-lnxbydmax*(l-1)))
       x(1:nmid) = log(w(1:nmid))/(l+1.d0)
       x(nmid+1:n) = -log(w(nmid+1:n))/(l-1.d0)
       x2 = (x(1)+x(2))/2.d0 
       fact1 =1.d0/(l+1.d0)
       fact2 = (0.5d0-l)/(1.5d0+l)/(l+3.d0)
       w(1) = exp(x2*(l+1))*(fact1+fact2*exp(x2*2))
       do i=2, nmid-1
          x1 = x2
          x2 = (x(i+1)+x(i))/2.d0
          w(i) = exp(x2*(l+1))*((fact1+fact2*exp(x2*2))-exp((x1-x2)*(l+1))*(fact1+fact2*exp(x1*2)))
       enddo
       x1 = x2
       w(nmid) = fact1+fact2 - exp(x1*(l+1))*(fact1+fact2*exp(x1*2))
       x2 = (x(nmid)+x(nmid+1))/2.d0
       x1 = 0.d0
       if(l.eq.2)then  
          !!in this case the integral actually diverges
          !!the hope is that it is not a global divergence and can be fixed by ignoring large x
          fact2 = -(0.5d0-l)/(1.5d0+l)/l
          w(nmid) = w(nmid) + x2 + fact2*(exp(-2.d0*x2)  - 1.d0)
          do i = nmid+1, n-1
             x1 = x2
             x2 = (x(i+1)+x(i))/2.d0
             w(i) = x2 - x1 + fact2 * exp(-2.d0*x2)*(1.d0-exp(-2.d0*(x1-x2)))
          enddo
          x1 = x2
          x2 = x(n)
          w(n) = x2 - x1 + fact2 * exp(-2.d0*x2)*(1.d0-exp(-2.d0*(x1-x2)))
       else
          fact1 = -1.d0/(l-2.d0)
          fact2 = -(0.5d0-l)/(1.5d0+l)/l
          w(nmid) = w(nmid) + fact1*( exp((2-l)*x2)-1.d0) &
               + fact2*(exp(-l*x2) - 1.d0)
          do i = nmid+1, n-1
             x1 = x2
             x2 = (x(i+1)+x(i))/2.d0
             w(i) = fact1*(exp((2-l)*x2)*(1.d0-exp((2-l)*(x1-x2)))) &
                  + fact2*(exp(-l*x2)*(1.d0-exp(-l*(x1-x2))))
          enddo
          x1 = x2
          x2 = x(n)
          w(n) = fact1*(exp((2-l)*x2)*(1.d0-exp((2-l)*(x1-x2)))) &
                  + fact2*(exp(-l*x2)*(1.d0-exp(-l*(x1-x2))))
       endif
       w(1:n) = w(1:n) * ( coop_pi/16.d0 * d**2 / (l-0.5d0) / (l+0.5d0) )
       x(1:n) = exp(x(1:n)) * d

    case default
       write(*,*) "m = ", m
       stop "jljlint_setw input : wrong argument"
    end select
  end subroutine jljlint_setw

  function smoothstep(x, width) result(s)
    COOP_REAL x, width, s
    s = (tanh(x/width)+1.d0)/2.d0
  end function smoothstep

  
  recursive function Hypergeometric2F1(a, b, c, x) result(f)
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
       f = gamma_product( (/ c, c-a-b, c-a, c-b /), (/ 1, 1, -1, -1 /) ) 
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
       f = gamma_product( (/ c, cmab, cma, cmb /), (/ 1, 1, -1, -1 /) ) * (f + gamma_product( (/ cma, cmb, -cmab, a, b, cmab /), (/ 1, 1, 1, -1, -1, -1 /), cmab*dlog(y) )*s)
       return
    endif
    if(n.eq.0)then
       psi1 = -coop_eulerC
       psia = digamma(a)
       psib = digamma(b)
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
       f = gamma_product( (/ c, a, b /), (/ 1, -1, -1 /) ) * s
       return
    endif
    !!---------- if c-a-b is almost a positive COOP_INT and x > xbound, we do interpolation ---------
    r1 = gamma_product( (/ n+eps, b+eps, a+eps, cmab, cma, cmb /), (/ 1, -1, -1, -1, 1, 1 /) )
    r2 = gamma_product( (/ n-eps, b-eps, a-eps, cmab, cma, cmb /), (/ 1, -1, -1, -1, 1, 1 /) )
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
    f = Hypergeometric2F1(a, b, a + b + n + eps, x )*w1 +  Hypergeometric2F1(a, b, a + b + n - eps, x )*w2
    return
  end function Hypergeometric2F1

  function sphericalBesselCross(ell1, ell2, r1, r2, ns) result(cr)
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
       cr = 2.d0**(-4+ns)*coop_pi* (2.d0*R/(r1+r2))**(1-ns)*Hypergeometric2F1((-2-l1+l2+ns)/2.d0,(-1+l1+l2+ns)/2.d0,1.5d0+l2,ratio**2)* dexp(l2*dlog(ratio)+lnGamma((-1+l1+l2+ns)/2.d0) -lnGamma(l2+1.5d0) -lnGamma((4+l1-l2-ns)/2.d0))
    else
       cr = 2.d0**(-4+ns)*dsin(coop_pio2*(4+l1-l2-ns)) * (2.d0*R/(r1+r2))**(1-ns)*Hypergeometric2F1((-2-l1+l2+ns)/2.d0,(-1+l1+l2+ns)/2.d0,1.5d0+l2,ratio**2)* dexp(l2*dlog(ratio)+lnGamma((-1+l1+l2+ns)/2.d0) -lnGamma(l2+1.5d0) + lnGamma(1.d0 - (4+l1-l2-ns)/2.d0))
    endif
  end function sphericalBesselCross


  Function log2_s(x) 
    COOP_REAL x,log2_s
    log2_s= log(x)/coop_ln2
  End Function log2_s

  function log2_v(x)
    COOP_REAL,dimension(:),intent(in)::x
    COOP_REAL log2_v(size(x))
    COOP_INT i
    !$omp parallel do
    do i=1, size(x)
       log2_v(i) = log2_s(x(i))
    end do
    !$omp end parallel do
  end function log2_v


  function is_integer(x) result(is)
    logical is
    COOP_REAL x
    is = (x - nint(x)) .lt. 1.d-6
  end function is_integer


end module coop_special_function_mod


