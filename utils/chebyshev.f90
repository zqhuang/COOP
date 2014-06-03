module cheb_utils
  use interpolation_utils
  implicit none

  type chebfit_vars
     integer n
     real(dl),dimension(:,:),allocatable::c
     real(dl)::lbd, rbd, a
  end type chebfit_vars

  interface chebfit_eval_vars
     module procedure chebfit_eval_vars_s, chebfit_eval_vars_v
  end interface chebfit_eval_vars

  INTERFACE Cheb_T
     MODULE PROCEDURE D_Cheb_T,D_Cheb_T_V
  end INTERFACE

  INTERFACE Cheb_U
     MODULE PROCEDURE D_Cheb_U,D_Cheb_U_V
  end INTERFACE

  interface get_cheb_value
     module procedure get_cheb_value_s, get_cheb_value_v
  end interface get_cheb_value

  INTERFACE Cheb_Value
     MODULE PROCEDURE Cheb_Value_s,Cheb_Value_V, cheb_value_ss, cheb_value_vv
  end INTERFACE

  Interface ChebFit
     MODULE PROCEDURE d_ChebFit,d_ChebFit2, chebfit_s
  End Interface

  Interface ChebFit_value
     MODULE PROCEDURE d_ChebFit_value,d_ChebFit_value_v, ChebFit_value_ss, chebfit_value_vv
  End Interface

  Interface ChebFit_derv
     MODULE PROCEDURE d_ChebFit_derv,d_ChebFit_derv_v
  End Interface

contains


!!=============== evaluating chebyshevs =========================

  !!Chebyshev Polynomials
  !!%%%%%%%%%%%%%%%%% POLYCheb %%%%%%%%%%%%%%%%%%%%%%%!!
  !! CALCULATE THE CHEBYSHEV POLYNOMIAL COEFFICIENTS
  !! T_n(X)=c_1+ c_2* x + c_3 * X^2+...+ c_{n+1} * X^N
  !! THE ARRAY c IS RETURNED IN POLYCheb
  function PolyCheb(n)
    integer(IB) N,R
    real(DL) POLYCheb(n + 1)
    IF(N.EQ.0)THEN
       POLYCheb(1)=1.D0
       RETURN
    endIF
    POLYCheb=0.D0
    POLYCheb(N+1)=2.D0**(N-1)
    DO R=1,N/2
       POLYCheb(N+1-2*R)=-POLYCheb(N+3-2*R)*(N-2*R+1)*(N-2*R+2)/R/(N-R)/4.D0
    endDO
  end function POLYCheb

  !!%%%%%%%%%% from nodal point values to chebyshev coefficients %%%%%%%%%%%%%%%%%%%%%%%
  !! GIVEN y(x_i), i=1,2,...,n where x_i is the i-th nodel point (T_N(x_i)=0 , x_1<x_2<...<x_n)
  !! one can perform c=ChebTRANS(n)*Y 
  !! to get chebyshev interpolation y(x)=c_1+c_2 T_1(x) +...+c_n T_{N-1}(x)
  subroutine Get_ChebTRANS(n, chebTrans)
    integer(IB) n,i,j
    real(dl) ChebTrans(n,n)
    ChebTrans(1,:)=0.5d0
    do i=2,n
       do j=1,n
          ChebTrans(i,n+1-j)=dcos(const_pi*(i-1)*(j-0.5d0)/N)
       enddo
    enddo
    ChebTrans = (2.d0/N)*ChebTrans
  End subroutine Get_ChebTRANS

  function Mat_ChebTrans(n) result(chebtrans)
    integer(IB) n
    real(dl) ChebTrans(n,n)
    call get_chebtrans(n,chebtrans)
  end function Mat_ChebTrans

  !! caluculate T_n(x)
  function D_Cheb_T(n,x)
    real(DL) x,D_Cheb_T,y1, y2, twox
    integer(IB) n ,i
    select case(n)
    case(0)
       D_Cheb_T = 1.d0
       return
    case(1)
       D_cheb_T = x
       return
    case(2)
       D_cheb_T = 2.d0*x**2 - 1.d0
       return
    case(3)
       d_cheb_T = x*(4.d0*x**2-3.d0)
       return
    case(4)
       d_cheb_T = 8.d0*x**2*(x**2-1.d0)+1.d0
       return
    case(5:25)
       twox=2.d0*x
       y1 = twox*x - 1.d0
       y2 = twox*y1 - x
       do i=4,N-2
          D_Cheb_T=Twox*y2-y1
          y1=y2
          y2=D_Cheb_T
       enddo
       D_cheb_T = twox*y2 - y1
       D_cheb_T = twox*d_cheb_T - y2
    case default
       D_Cheb_T=dcos(n*dacos(X))
       return
    end select
  end function D_Cheb_T

  Function D_Cheb_T_V(N,x)
    real(dl),dimension(:),intent(in)::x
    integer(IB) n,i
    real(dl) D_Cheb_T_V(size(x))
    do i=1,size(x)
       d_Cheb_T_V(i)=d_Cheb_T(N,x(i))
    enddo
  End Function D_Cheb_T_V


  function D_Cheb_U(N,X)
    real(DL) X,THETA,D_Cheb_U,Y1,Y2,TWOX
    integer(IB) N,I
    IF(N.EQ.0)THEN
       D_Cheb_U=1.D0
       RETURN
    endIF
    IF(N.GT.20)THEN
       THETA=DACOS(X)
       IF(ABS(THETA).LT.1.d-8)THEN
          D_Cheb_U=N+1
          RETURN
       ELSE
          D_Cheb_U=DSIN((N+1)*THETA)/DSIN(THETA)
          RETURN
       endIF
    endIF
    TWOX=2.D0*X
    Y1=1.D0
    Y2=TWOX
    D_Cheb_U=TWOX
    DO I=1,N-1
       D_Cheb_U=TWOX*Y2-Y1
       Y1=Y2
       Y2=D_Cheb_U
    endDO
  end function D_Cheb_U


  Function D_Cheb_U_V(N,x)
    real(dl),dimension(:),intent(in)::x
    integer(IB) n,i
    real(dl) D_Cheb_U_V(size(x))
    do i=1,size(x)
       d_Cheb_U_V(i)=d_Cheb_U(N,x(i))
    enddo
  End Function D_Cheb_U_V

  subroutine cheb_T_eval_all(n, x, y)
    integer,intent(IN):: n
    real(dl),intent(IN):: x
    real(dl),intent(OUT)::y(0:n)
    real(dl) twox
    integer i
    y(0) = 1.d0
    if(n.eq.0) return
    y(1) = x
    twox = 2.d0*x
    do i=2, n
       y(i) = twox * y(i-1) - y(i-2)
    enddo
  end subroutine cheb_T_eval_all


  subroutine cheb_U_eval_all(n, x, y)
    integer n
    real(dl) x
    real(dl),intent(OUT)::y(0:n)
    integer i
    y(0) = 1.d0
    if(n.eq.0) return
    y(1) = 2.d0*x
    do i=2, n
       y(i) = y(1) * y(i-1) - y(i-2)
    enddo
  end subroutine cheb_U_eval_all



  !!%%%%%%%%%%%%%%%%% ChebFit %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !! least square fit y=c(1)+c(2)*T_1(t)+... in interval [a,b]
  !!  where t=2(x-a)/(b-a)-1
  Subroutine D_ChebFit(x,y,a,b,C)
    real(dl),intent(in)::a,b
    real(dl),dimension(:),INTENT(IN)::X,Y
    real(DL),dimension(:),INTENT(OUT)::C
    real(dl) Fx(size(X),size(C)),t(size(x))
    integer(IB) n,m,i
    n=getdim("ChebFit",size(x),size(y))
    m=size(c)
    if(m.gt.n) call ReturnError("ChebFit")
    if(m.gt.50)print*,"ChebFit warning: may be not stable"
    t=2.d0*(x-a)/(b-a)-1.d0
    do i=1,n
       call Cheb_T_eval_all(m-1,t(i), Fx(i,1:m))
    enddo
    call fit_template(n, m, y, fx, c)
  End subroutine D_ChebFit

  subroutine ChebFit_s(n, t, y, m, c)
    integer,intent(IN):: n, m
    real(dl),intent(IN)::y(n)
    real(dl) t(n),  c(m)
    real(dl) Fx(n, m)
    integer i,j
    if(m.gt.n) call ReturnError("ChebFit: not enough data")
    do j=1,m
       do i=1,n
          call Cheb_T_eval_all(m-1, t(i), Fx(i,1:m))
       enddo
    enddo
    call fit_template(n, m, y, fx, c)
  End subroutine ChebFit_s


  subroutine D_CHEBFIT2(X,Y,Sigma,a,b,C)
    real(dl),intent(in):: a,b
    real(DL),dimension(:),INTENT(IN)::X,Y,Sigma
    real(DL),dimension(:),INTENT(OUT)::C
    real(dl) Fx(size(X),size(C)),t(size(x))
    integer(IB) n,m,i
    n=getdim("ChebFit2",size(x),size(y),size(Sigma))
    m=size(c)
    if(m.gt.n.or.any(sigma.le.0.d0)) call ReturnError("ChebFit2")
    if(m.gt.30)print*,"ChebFit2 warning: may be not stable"
    t=2.d0*(x-a)/(b-a)-1.d0
    do i=1,n
       call Cheb_T_eval_all(m-1,t(i), Fx(i,:))
       Fx(i,:) = fx(i,:)/Sigma(i)
    enddo
    call fit_template(n, m, y/sigma, fx, c)
  end subroutine D_CHEBFIT2


  function Cheb_value_s(c,x) result(cheb)
    !! return c(1)*T_0(x)+c(2)*T_1(x)+...+c(n)*T_{n-1}(x)
    real(DL),dimension(:),INTENT(IN)::C
    real(DL) X, cheb
    call get_cheb_value_s(c, size(c), x, cheb)
  end function Cheb_value_s


  function Cheb_Value_ss(C, n, x) result(cheb)
    !! return c(1)*T_0(x)+c(2)*T_1(x)+...+c(n)*T_{n-1}(x)
    !! using Clenshaw summation
    integer n
    real(dl) c(n), x, cheb
    call get_cheb_value_s(c, n, x, cheb)
  end function Cheb_Value_ss


  function Cheb_value_V(C, x) result(cheb)
    !! return c(1)*T_0(x)+c(2)*T_1(x)+...+c(n)*T_{n-1}(x)
    integer(IB) i,n, m
    real(dl),dimension(:),intent(in)::c,x
    real(dl) cheb(size(x))
    n = size(x)
    m = size(c)
    do i=1, n
       call get_cheb_value_s(c, m, x(i), cheb(i))
    enddo
  end function Cheb_value_V

  function Cheb_value_Vv(C, n, x, m) result(cheb)
    !! return c(1)*T_0(x)+c(2)*T_1(x)+...+c(n)*T_{n-1}(x)
    integer(IB) i,n, m
    real(dl) c(n), x(m)
    real(dl) cheb(m)
    do i = 1, m
       call get_cheb_value_s(c, m, x(i), cheb(i))
    enddo
  end function Cheb_value_Vv

  subroutine get_cheb_value_s(c, n, x, cheb)
    integer n
    real(dl) c(n), x, cheb
    real(dl) twox, y1, y2, y3
    integer i
    twox = x+x
    y3 = 0.d0
    y2 = 0.d0
    do i = n, 2, -1
       y1 = c(i) + twox * y2 - y3
       y3 = y2
       y2 = y1
    enddo
    cheb = c(1) + x* y2 - y3    
  end subroutine get_cheb_value_s

  subroutine get_cheb_value_v(c, n, x, m, cheb)
    integer(IB) i,n, m
    real(dl) c(n), x(m)
    real(dl) cheb(m)
    do i = 1, m
       call get_cheb_value_s(c, n, x(i), cheb(i))
    enddo
  end subroutine get_cheb_value_v

  Function D_CHEBFIT_Value(a,b,c,x) result(cheb)
    real(dl),dimension(:),intent(in)::c
    real(dl) X
    real(dl) a,b
    real(dl) cheb
    real(dl) t
    t=2.d0*(x-a)/(b-a)-1.d0
    call get_Cheb_value(c, size(c), t, cheb)
  end function D_CHEBFIT_Value



  function D_ChebFit_value_V(a,b,c,x) result(cheb)
    real(dl),dimension(:),intent(in)::c,x
    real(dl) a,b, s
    integer(IB) i, n
    real(dl) cheb(size(x))
    n = size(c)
    s = 2.d0/(b-a)
    do i=1,size(x)
       call get_cheb_value(c, n, s*(x(i)-a)-1.d0, cheb(i))
    enddo
  end function D_ChebFit_value_V


  function ChebFit_value_ss(a,b,c,n, x) result(cheb)
    integer n
    real(dl) a,b,c(n), x, cheb
    call get_cheb_value(c(1:n), n, 2.d0*(x-a)/(b-a)-1.d0, cheb)
  end function ChebFit_value_ss

  function ChebFit_value_vv(a,b,c,n, x, m) result(cheb)
    integer n, m
    real(dl) a,b,c(n), x(m), cheb(m)
    call get_cheb_value(c(1:n), n, 2.d0*(x-a)/(b-a)-1.d0, m, cheb)
  end function ChebFit_value_vv


  function D_ChebFit_derv(a,b,c,x)
    real(dl),dimension(:),intent(in)::c
    real(dl) x,a,b,d_ChebFit_derv,twox,y1,y2,y3
    integer(IB) i,n
    n=size(c)
    if(n.eq.1)then
       d_ChebFit_derv=0._dl
       return
    endif
    twox=4._dl * (x-a) / (b-a) - 2._dl
    Y1=1._dl
    Y2=TWOX
    y3=TWOX
    d_ChebFit_derv=y1*c(2)
    DO i = 3, n
       y3=twox * Y2-Y1
       Y1=Y2
       Y2=y3
       d_ChebFit_derv=d_ChebFit_derv+y1*c(i)*(i-1)
    endDO
    d_ChebFit_derv=d_ChebFit_derv*2./(b-a)
  end function D_ChebFit_derv


  function D_ChebFit_derv_V(a,b,c,x)
    real(dl),dimension(:),intent(in)::c,x
    real(dl) a,b
    integer(IB) i
    real(dl) D_ChebFit_derv_V(size(x))
    do i=1,size(x)
       D_ChebFit_derv_v(i)=D_ChebFit_derv(a,b,c,x(i))
    enddo
  end function D_ChebFit_derv_V

  function ChebFit_int(a, b, c, lbd, ubd)
    real(dl) a,b, xinf, xsup, ChebFit_int, lbd, ubd
    real(dl),dimension(:),intent(in)::c
    real(dl) twoxinf, twoxsup, y1inf, y2inf, y3inf, y1sup, y2sup, y3sup
    integer(IB) i, n 
    xinf = 2._dl*(lbd-a)/(b-a)-1._dl
    xsup = 2._dl*(ubd-a)/(b-a)-1._dl
    twoxinf = 2._dl*Xinf
    twoxsup = 2._dl*xsup
    N=SIZE(c)
    ChebFit_int = c(1)*(xsup - xinf)*2._dl
    if(n.ge.2)then
       Y3inf = twoxinf * xinf - 1._dl
       y3sup = twoxsup * xsup - 1._dl
       ChebFit_int = ChebFit_int + c(2)*(y3sup - y3inf)/2._dl 
       Y1inf = xinf
       Y2inf = y3inf
       y1sup = xsup
       y2sup = y3sup
       DO i = 3, n
          Y3inf = twoxinf * Y2inf - Y1inf
          y3sup = twoxsup * y2sup - y1sup
          ChebFit_int = ChebFit_int + c(i)*((y3sup - y3inf)/i - (y1sup - y1inf)/(i-2))
          Y1inf = Y2inf
          Y2inf = Y3inf
          y1sup = y2sup
          y2sup = y3sup
       endDO
    endif
    ChebFit_int = ChebFit_int * (b-a)/4._dl
  end function ChebFit_int

  subroutine cheb_derv(n, c, t, y)
    integer n
    real(dl) c(n)
    real(dl) t, y,twox,y1,y2,y3
    integer i 
    twox=4.d0 *t  - 2.d0
    y3 = twox*twox - 1.d0
    Y1 = twox
    y2 = y3
    y=c(2)
    do i=3,n-1
       y=y+y1*c(i)*(i-1)
       y3=twox * y2-Y1
       Y1=Y2
       Y2=y3
    enddo
    y= (y+y1*c(n)*(n-1))* 2.d0
  end subroutine cheb_derv


  subroutine cheb_derivative(c, n, cder)
    INTEGER n
    REAL(dl) c(n),cder(n-1)
    INTEGER j
    cder(n-1) = 2*(n-1)*c(n)
    do j = n-2,1,-1
       cder(j)=cder(j+2)+2*j*c(j+1)
    enddo
    cder(1) = cder(1)/2.d0
    return
  end subroutine cheb_derivative

  subroutine chebfit_derivative(a, b, c, n, cder)
    INTEGER n
    REAL(dl) c(n), cder(n-1), a, b
    INTEGER j
    cder(n-1) = 2*(n-1)*c(n)
    do j = n-2,1,-1
       cder(j)=cder(j+2)+2*j*c(j+1)
    enddo
    cder(1) = cder(1)/2.d0
    cder = cder*(2.d0/(b-a))
    return
  end subroutine chebfit_derivative



  subroutine chebfit_all(n, x, y, a, b, c, m, sigma)
    integer,intent(IN):: n,m
    real(dl),intent(IN)::x(n), y(n), a, b
    real(dl),intent(OUT)::c(0:m, 0:m)
    real(dl),optional::sigma(n)
    real(dl) fx(n, 0:m)
    integer i
    do i = 1, n
       call cheb_T_eval_all(m, 2.d0*(x(i)-a)/(b-a)-1.d0, Fx(i, 0:m))
       if(present(sigma)) fx(i,0:m) = fx(i,0:m)/sigma(i)
    enddo
    if(present(sigma))then
       call fit_template(n, m+1, y/sigma, fx, c(0:m, 0))
    else
       call fit_template(n, m+1, y, fx, c(0:m, 0))
    endif
    do i = 1, m
       call chebfit_derivative(a, b, c(0:m-i+1, i-1), m-i+2, c(0:m-i, i))
       c(m-i+1:m, i) = 0.d0
    enddo
  end subroutine chebfit_all

  subroutine chebfit_eval_any(a, b, c, m, x, y, np)
    integer,intent(IN)::m
    real(dl),intent(IN):: a, b, c(0:m, 0:m), x
    real(dl) ,intent(OUT)::y
    integer,optional::np
    if(present(np))then
       if(np.gt.m)then
          y = 0
          return
       else
          call get_cheb_value_s(c(0:m-np, np), m-np+1, 2.d0*(x-a)/(b-a)-1.d0, y)           
       endif
    else
       call get_cheb_value_s(c(0:m, 0), m+1, 2.d0*(x-a)/(b-a)-1.d0, y)    
    endif
  end subroutine chebfit_eval_any

  subroutine chebfit_all_vars(n, x, y, lbd, rbd, m, chvar, sigma)
    integer,intent(IN)::n, m
    real(dl),intent(IN)::x(n), y(n), lbd, rbd
    real(dl),optional::sigma(n)
    type(chebfit_vars) chvar
    chvar%n = m
    if(allocated(chvar%c))then
       if(chvar%n .ne. m)then
          deallocate(chvar%c)
          allocate(chvar%c(0:m,0:m))
       endif
    else
       allocate(chvar%c(0:m,0:m))
    endif
    chvar%lbd = lbd
    chvar%rbd = rbd
    chvar%a = 2.d0/(rbd-lbd)
    if(present(sigma))then
       call chebfit_all(n, x, y, lbd, rbd, chvar%c, m, sigma)    
    else
       call chebfit_all(n, x, y, lbd, rbd, chvar%c, m)    
    endif
  end subroutine chebfit_all_vars

  subroutine chebfit_eval_vars_s(chvar, x, y, np)
    type(chebfit_vars),intent(IN)::chvar
    real(dl),intent(IN):: x
    real(dl) ,intent(OUT)::y
    integer,optional::np
    if(present(np))then
       if(np.gt.chvar%n)then
          y = 0
          return
       else
          call get_cheb_value_s(chvar%c(0:chvar%n-np, np), chvar%n-np+1, chvar%a*(x-chvar%lbd)-1.d0, y)           
       endif
    else
       call get_cheb_value_s(chvar%c(0:chvar%n, 0), chvar%n+1, chvar%a*(x-chvar%lbd)-1.d0, y)    
    endif
  end subroutine chebfit_eval_vars_s


  subroutine chebfit_eval_vars_v(chvar, x, y, np)
    type(chebfit_vars),intent(IN)::chvar
    real(dl),dimension(:),intent(IN):: x
    real(dl),dimension(:),intent(OUT)::y
    integer,optional::np
    integer i
    if(present(np))then
       if(np.gt.chvar%n)then
          y = 0
          return
       else
          do i=1,size(x)
             call get_cheb_value_s(chvar%c(0:chvar%n-np, np), chvar%n-np+1, chvar%a*(x(i)-chvar%lbd)-1.d0, y(i))
          enddo
       endif
    else
       do i=1,size(x)
          call get_cheb_value_s(chvar%c(0:chvar%n, 0), chvar%n+1, chvar%a*(x(i)-chvar%lbd)-1.d0, y(i))    
       enddo
    endif
  end subroutine chebfit_eval_vars_v


end module cheb_utils
