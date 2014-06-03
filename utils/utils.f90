Module general_utils !!this module is thread safe
  !!use precision
  use basic_utils
  use matrix_utils
  use specfunc_utils
  use string_utils
  use random_utils
  use interpolation_utils
  use integrate_utils
  use cheb_utils
  implicit none
  !!==============================

  !!======================================
  Interface isNan
     Module Procedure isNan_r, isNan_d, isNan_i, isnan_v
  End Interface

  Interface PolyValue
     Module Procedure PolyValue_S, PolyValue_V
  End interface PolyValue

contains


  Function R_of_Chi(chi,omk) !!using unit c/H_0
    real(dl) chi,omk
    real(dl) R_of_Chi
    if(omk.eq.0.d0)then
       R_of_Chi=chi
       return
    elseif(omk.gt.0.)then !!open universe, k<0
       R_of_Chi=sinhc(sqrt(omk)*chi)*chi
       return
    else
       R_of_Chi=sinc(sqrt(-omk)*chi)*chi !!closed universe, k>0
       return
    Endif
  End Function R_of_Chi


  !!polynomials
  function PolyCoeff(xa,ya)
    real(Dl), DIMENSION(:), INTENT(IN) :: xa,ya
    real(Dl), DIMENSION(size(xa)) :: PolyCoeff
    integer(IB):: j,k,m,n
    real(Dl) :: dy
    real(Dl), DIMENSION(size(xa)) :: x,y
    n=getdim("polcof",size(xa),size(ya))
    x=xa
    y=ya
    do j=1,n
       m=n+1-j
       call d_polint(x(1:m),y(1:m),m,0.d0,PolyCoeff(j),dy)
       k=iminloc(abs(x(1:m)))
       where (x(1:m) /= 0.0) y(1:m)=(y(1:m)-PolyCoeff(j))/x(1:m)
       y(k:m-1)=y(k+1:m)
       x(k:m-1)=x(k+1:m)
    end do
  end function PolyCoeff


  !!just for testing, you should use chebyshev fitting for serious application
  subroutine polyfit(n, x, y, m, p)
    integer n, m
    real(dl),intent(IN):: x(n), y(n)
    real(dl),intent(OUT):: p(m)
    real(dl) fx(n, m), norm, xnorm(n)
    integer i
    norm = maxval(abs(x))
    if(norm .eq.0.d0)then
       p = 0.d0
       return
    endif
    xnorm = x/norm
    fx(:, 1) = 1.d0
    do i=2, m
       fx(:, i) = fx(:, i-1)*xnorm
    enddo
    call fit_template(n, m, y, fx, p)
  end subroutine polyfit


  subroutine get_polyvalue_arr(p, n, m, x, y)
    integer n, m
    real(dl) p(n)
    real(dl) x(m), y(m)
    integer j
    y = p(n)
    do j= n-1,1,-1
       y= y*x + p(j)
    enddo
  end subroutine get_polyvalue_arr


  subroutine get_poly_dervs(p, n, m)
    integer n, m
    real(dl) p(n,0:m)
    integer i, j
    do i= 1, m
       p(n-i+1:n,i) = 0.d0
       do j=1, n-i
          p(j, i) = p(j+1, i-1)*j
       enddo
    enddo
  end subroutine get_poly_dervs

  function Polyvalue_s(P,X)
    real(DL) X
    real(DL),dimension(:),INTENT(IN)::P
    real(DL) POLYValue_s
    integer(IB) M,J
    real(DL),dimension(:),ALLOCATABLE::VEC
    real(DL) POW
    M=SIZE(P)
    IF(M.LT.25)THEN
       POLYValue_s=P(M)
       DO J=M-1,1,-1
          POLYValue_s=POLYValue_s*X+P(J)
       endDO
       RETURN
    endIF
    allocate(VEC(M+1))
    POW=X
    VEC(1:M)=P
    DO
       VEC(M+1)=0.D0
       J=ISHFT(M+1,-1)
       VEC(1:J)=VEC(1:M:2)+POW*VEC(2:M+1:2)
       IF(J .LE. 1) EXIT
       POW=POW*POW
       M=J
    endDO
    POLYValue_s=VEC(1)
    DEallocate(VEC)
  end function POLYValue_S

  subroutine Polyvalue_eval(psize, P, X, val)
    integer(IB),intent(in):: psize
    real(dl), intent(in) :: x
    real(dl), intent(out) :: val
    real(DL),intent(in)::p(psize)
    integer j,M
    real(DL) vec(psize + 1)
    real(DL) pow
    M = psize
    pow = x
    vec(1:M) = P
    do
       vec(M+1)=0.D0
       J = ISHFT(M+1,-1)
       VEC(1:J)=VEC(1:M:2)+POW*VEC(2:M+1:2)
       IF(J .LE. 1) EXIT
       POW=POW*POW
       M=J
    enddo
    val = vec(1)
  end subroutine Polyvalue_eval

  function PolyValue_v(P,X)
    real(dl) x(:)
    real(DL),dimension(:),INTENT(IN)::P
    real(DL) POLYValue_v(SIZE(X))
    integer(IB) i
    do i=1, size(X)
       PolyValue_v(i)= Polyvalue_s(P, x(i))
    enddo
  end function POLYValue_v





  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Function IsNan_d(x)
    real(dl) x
    Logical IsNan_d
    IsNaN_d=.not.((x.gt.0.) .or. (x.le.0.))
  End Function IsNan_d

  Function IsNan_r(x)
    real x
    Logical IsNan_r
    IsNaN_r=.not.((x.gt.0.) .or. (x.le.0.))
  End Function IsNan_r

  Function IsNan_i(x)
    integer(IB) x
    Logical IsNan_i
    IsNaN_i=.not.((x.gt.0.) .or. (x.le.0.))
  End Function IsNan_i

  function isnan_v(x)
    real(dl),dimension(:),intent(in):: x
    logical isnan_v
    integer i
    do i=1, size(x)
       if(isnan_d(x(i)))then
          isnan_v = .true.
          return
       endif
    enddo
    isnan_v = .false.
  end function isnan_v



  !!%%%%%%%%%%%%%%%%%% for multi-frequency correlations and weak lensing tomography %%%%%%%%%%%%%%%

  function arr_ind(I,J)
    !!put (i,j) pairs into an array, arr_ind is the index of (i,j) in that array
    integer(IB) I,J,arr_IND
    IF(I.GE.J)THEN
       arr_IND=I*(I-1)/2+J
    ELSE
       arr_IND=J*(J-1)/2+I
    endIF
  end function arr_ind


  !!the inverse routine that return i 
  function arr_IND1(N) !!find i such that i(i-1)/2<N<=i(i+1)/2
    Integer(IB) arr_IND1,N
    select case(N)
    case(0)
       arr_Ind1 = 0
       return
    case(1)
       arr_Ind1 = 1
       return
    case(2:3)
       arr_Ind1 = 2
       return
    case(4:6)
       arr_Ind1 = 3
       return
    case(7:10)
       arr_Ind1 = 4
       return
    case(11:15)
       arr_Ind1 = 5
       return
    case(16:21)
       arr_Ind1 = 6
    case default
       arr_Ind1=nint(sqrt(2.d0*n))  
       return
    end select
  end function arr_IND1

  !!the inverse routine that return N-i(i-1)/2 
  function arr_IND2(N)
    integer(IB) i,N,arr_Ind2
    i=arr_ind1(N)
    arr_IND2=N-i*(i-1)/2
  end function arr_IND2







  Subroutine Export_1D_distr(sample, m, filename)
    real(dl),dimension(:),intent(in)::sample
    character(LEN=*) filename
    integer(IB) m
    integer(IB) n,i,j
    real(dl) xmin,dx, xmax, prob(m)
    n=size(sample)
    xmin = minval(sample)
    xmax = maxval(sample)
    dx = ( xmax - xmin)/m
    if(dx.le.0.d0) stop "delta distribution can not be exported"
    prob = 0.d0
    do i=1,n
       j=ceiling((sample(i)-xmin)/dx)
       j = min(max(1,j), m)
       prob(j) = prob(j)+1
    enddo
    prob = prob/sum(prob)/dx
    open(tmp_file_unit, file=trim(filename))
    do j=1,m
       write(tmp_file_unit,'(2E16.7)') xmin+(j-0.5)*dx, prob(j)
    enddo
    close(tmp_file_unit)
  End Subroutine Export_1D_distr

  subroutine DervFit(fcn, fcnp, x, maxshift)
    external fcn, fcnp
    real(dl),dimension(:),intent(in)::maxshift
    real(dl),dimension(:),intent(inout)::x
    real(dl) derv(size(x)), fnow, fupdate
    integer j, l, n
    integer,parameter::maxnloop=10
    n = GetDim("generalFit", size(x), size(maxshift))
    call fcnp(x, derv)
    call fcn(x, fnow)
    do j=1, n
       if(abs(derv(j))*maxshift(j) .lt. 1.)then
          derv(j) = sign(maxshift(j), derv(j))
       else
          derv(j) = abs(fnow)*0.01/derv(j)
       endif
    enddo
    do l=1, maxnloop
       call fcn(x + derv, fupdate)
       if(fupdate .lt. fnow)then
          x = x + derv
       else
          derv = derv /2.
       endif
    enddo
  end subroutine DervFit

  subroutine MCMCFit(fcn, x, propose, xmin, xmax)
    real(dl),dimension(:),intent(in)::propose, xmin, xmax
    real(dl),dimension(:),intent(inout)::x
    external fcn
    real(dl) xupdate(size(x)), xbest(size(x)), fbest, fnow, fupdate
    integer i, l, n
    integer,parameter::maxnloop=100000
    n = GetDim("generalFit", size(xmin), size(xmax))
    call fcn(x, fnow)
    fbest = fnow
    xbest = x
    do l = 1, maxnloop
       do i=1,n
          xupdate(i) = random_Gaussian()
       enddo
       xupdate = x + propose * xupdate
       if(any(xupdate .lt. xmin) .or. any(xupdate .gt. xmax)) cycle
       call fcn(xupdate, fupdate)
       if(fnow.ge.fupdate)then
          x=xupdate
          fnow=fupdate
          if(fbest.gt.fnow)then
             fbest=fnow
             xbest=x
          endif
          cycle
       endif
       if( fupdate - fnow .le. random_exp())then
          x=xupdate
          fnow=fupdate
       endif
    enddo
    x = xbest
  end subroutine MCMCFit

  !!return log(exp(x) + exp(y))
  function LogSum(x, y)
    real(dl) x, y
    real(dl) Logsum
    if(x.gt.y)then
       LogSum = x + dlog(1._dl + dexp(y-x))
    else
       LogSum = y + dlog(1._dl + dexp(x-y))
    endif
  end function LogSum



  subroutine numderv_multiscale(n, m, dx, y, yp)
    !!n = 2*(m+m/5+m/6)
    integer n, m
    real(dl) dx, y(-n:n), yp(0:m), dscale
    yp(0) = y(0)
    if(m.eq.0) return
    if(n.lt.1)then
       yp(1:m) = 0
       return
    endif
    yp(1) = (y(1)-y(-1))/(2.d0*dx)
    if(m.eq.1 )return
    if(n.lt.2)then
       yp(2:m) = 0
       return
    endif
    dscale = dx * dx
    yp(2) =(14*y(-2) - 3*y(-1) - 22*y(0) - 3*y(1) + 14*y(2))/(53.*dscale)
    if(m.eq.2)return
    if(n.lt.3)then
       yp(3:m) = 0
       return
    endif
    dscale = dscale * dx
    yp(3) =  (-16*y(-3) + 13*y(-2) + 22*y(-1) - 22*y(1) - 13*y(2) + 16*y(3))/(102.*dscale)
    if(m.eq.3) return
    if(n.lt.4)then
       yp(4:m) = 0
       return
    endif
    dscale = dscale * dx
    yp(4) =  (1994*y(-4) - 2509*y(-3) - 2457*y(-2) + 505*y(-1) + 4934*y(0) + &
         505*y(1) - 2457*y(2) - 2509*y(3) + 1994*y(4))/(22369.*dscale)
    if(m.eq.4) return
    if(n.lt.5)then
       yp(5:m) = 0
       return
    endif 
    dscale = dscale * dx
    yp(5) = (-20872*y(-5) + 37524*y(-4) + 16879*y(-3) - 25316*y(-2) - &
         45741*y(-1) + 45741*y(1) + 25316*y(2) - 16879*y(3) - &
         37524*y(4) + 20872*y(5))/(392578.*dscale)
    if(m.eq.5 ) return
    if(n.lt.6)then
       yp(6:m) = 0
       return
    endif
    dscale = dscale * dx
    yp(6) = (90619657.d0*y(-6) - 205461400.d0*y(-5) - &
         23134128.d0*y(-4) + 173805763.d0*y(-3) + 177853997.d0*y(-2) -  &
         31294459.d0*y(-1) - 364778860.d0*y(0) - 31294459.d0*y(1) + &
         177853997.d0*y(2) + 173805763.d0*y(3) - 23134128.d0*y(4) - &
         205461400.d0*y(5) + 90619657.d0*y(6))/(2.946985323e9_dl* dscale)
  end subroutine numderv_multiscale

  !!input: y(-n:n), where y(i) = f(i * dx)
  !!output: yp(0:m), where yp(i) is the i-th derivative at x=0
  subroutine numderv_polyfit(n, m, dx, y, yp)
    real(dl) dx
    integer n,m
    real(dl) y(-n:n), yp(0:m)
    if(m.gt. 2*n .or. m.gt.6) stop "numderv_polyfit argument out of range"
    select case(n)
    case(0)
       yp(0) = y(0)
    case(1)
       select case(m)
       case(0)
          yp(0) = sum(y(-1:1))/3.d0
       case(1)
          yp(0:1) = (/ sum(y(-1:1))/3.d0, (y(1)- y(-1))/(2.d0*dx) /)
       case(2)
          yp(0:2) = (/ y(0), &
               (y(1)- y(-1))/(2.d0*dx), &
               (y(1)+y(-1)-2.d0*y(0))/(dx**2) /)
       end select
    case(2)          
       select case(m)
       case(0)
          yp(0) = sum(y(-2:2))/5.d0
       case(1)
          yp(0:1) = (/ (y(-2) + y(-1) + y(0) + y(1) + y(2))/5., &
               (-2*y(-2) - y(-1) + y(1) + 2*y(2))/(10.*dx) /)
       case(2)
          yp(0:2) = (/ (-3*y(-2) + 12*y(-1) + 17*y(0) + 12*y(1) - 3*y(2))/35., &
               (-2*y(-2) - y(-1) + y(1) + 2*y(2))/(10.*dx), &
               -(-2*y(-2) + y(-1) + 2*y(0) + y(1) - 2*y(2))/(7.*dx**2) /)
       case(3)
          yp(0:3) = (/ (-3*y(-2) + 12*y(-1) + 17*y(0) + 12*y(1) - 3*y(2))/35., &
               (y(-2) - 8*y(-1) + 8*y(1) - y(2))/(12.*dx), &
               -(-2*y(-2) + y(-1) + 2*y(0) + y(1) - 2*y(2))/(7.*dx**2), &
               (-y(-2) + 2*y(-1) - 2*y(1) + y(2))/(2.*dx**3) /)
       case(4)
          yp(0:4) = (/ y(0), (y(-2) - 8*y(-1) + 8*y(1) - y(2))/(12.*dx), &
               -(y(-2) - 16*y(-1) + 30*y(0) - 16*y(1) + y(2))/(12.*dx**2), &
               (-y(-2) + 2*y(-1) - 2*y(1) + y(2))/(2.*dx**3), &
               (y(-2) - 4*y(-1) + 6*y(0) - 4*y(1) + y(2))/dx**4 /)
       end select
    case(3)
       select case(m)
       case(0) 
          yp(0) = sum(y(-3:3))/7.d0
       case(1)
          yp(0:1) = (/ sum(y(-3:3))/7.d0,  (-3*y(-3) - 2*y(-2) - y(-1) + y(1) + 2*y(2) + 3*y(3))/(28.*dx) /)
       case(2)
          yp(0:2) = (/ (-2*y(-3) + 3*y(-2) + 6*y(-1) + 7*y(0) + 6*y(1) + 3*y(2) - 2*y(3))/21., &
               (-3*y(-3) - 2*y(-2) - y(-1) + y(1) + 2*y(2) + 3*y(3))/(28.*dx), &
               (5*y(-3) - 3*y(-1) - 4*y(0) - 3*y(1) + 5*y(3))/(42.*dx**2) /)
       case(3)          
          yp(0:3) = (/ (-2*y(-3) + 3*y(-2) + 6*y(-1) + 7*y(0) + 6*y(1) + 3*y(2) - 2*y(3))/21., &
               (22*y(-3) - 67*y(-2) - 58*y(-1) + 58*y(1) + 67*y(2) - 22*y(3))/(252.*dx), &
               (5*y(-3) - 3*y(-1) - 4*y(0) - 3*y(1) + 5*y(3))/(42.*dx**2), &
               (-y(-3) + y(-2) + y(-1) - y(1) - y(2) + y(3))/(6.*dx**3) /)
       case(4)
          yp(0:4) = (/ (5*y(-3) - 30*y(-2) + 75*y(-1) + 131*y(0) + 75*y(1) - 30*y(2) + 5*y(3))/231., &
               (22*y(-3) - 67*y(-2) - 58*y(-1) + 58*y(1) + 67*y(2) - 22*y(3))/(252.*dx), &
               -(13*y(-3) - 67*y(-2) + 19*y(-1) + 70*y(0) + 19*y(1) - 67*y(2) + 13*y(3))/(132.*dx**2), &
               (-y(-3) + y(-2) + y(-1) - y(1) - y(2) + y(3))/(6.*dx**3), &
               (3*y(-3) - 7*y(-2) + y(-1) + 6*y(0) + y(1) - 7*y(2) + 3*y(3))/(11.*dx**4) /)
       case(5)
          yp(0:5) = (/ (5*y(-3) - 30*y(-2) + 75*y(-1) + 131*y(0) + 75*y(1) - 30*y(2) + 5*y(3))/231., &
               (-y(-3) + 9*y(-2) - 45*y(-1) + 45*y(1) - 9*y(2) + y(3))/(60.*dx), &
               -(13*y(-3) - 67*y(-2) + 19*y(-1) + 70*y(0) + 19*y(1) - 67*y(2) + 13*y(3))/(132.*dx**2), &
               (y(-3) - 8*y(-2) + 13*y(-1) - 13*y(1) + 8*y(2) - y(3))/(8.*dx**3), &
               (3*y(-3) - 7*y(-2) + y(-1) + 6*y(0) + y(1) - 7*y(2) + 3*y(3))/ (11.*dx**4), &
               (-y(-3) + 4*y(-2) - 5*y(-1) + 5*y(1) - 4*y(2) + y(3))/(2.*dx**5) /)
       case(6)
          yp(0:6) = (/ y(0), &
               (-y(-3) + 9*y(-2) - 45*y(-1) + 45*y(1) - 9*y(2) + y(3))/(60.*dx), &
               (2*y(-3) - 27*y(-2) + 270*y(-1) - 490*y(0) + 270*y(1) - 27*y(2) + 2*y(3))/(180.*dx**2), &
               (y(-3) - 8*y(-2) + 13*y(-1) - 13*y(1) + 8*y(2) - y(3))/(8.*dx**3), &
               -(y(-3) - 12*y(-2) + 39*y(-1) - 56*y(0) + 39*y(1) - 12*y(2) + y(3))/ (6.*dx**4), &
               (-y(-3) + 4*y(-2) - 5*y(-1) + 5*y(1) - 4*y(2) + y(3))/(2.*dx**5), & 
               (y(-3) - 6*y(-2) + 15*y(-1) - 20*y(0) + 15*y(1) -  6*y(2) + y(3))/dx**6 /)
       end select
    case(4)
       select case(m)
       case(0)
          yp(0) = sum( y(-4:4))/9.d0
       case(1)
          yp(0:1) = (/  sum( y(-4:4))/9.d0, (-4*y(-4) - 3*y(-3) - 2*y(-2) - y(-1) + y(1) + 2*y(2) + 3*y(3) + 4*y(4))/(60.*dx) /)
       case(2)
          yp(0:2) = (/ (-21*y(-4) + 14*y(-3) + 39*y(-2) + 54*y(-1) + 59*y(0) + 54*y(1) + 39*y(2) + 14*y(3) - 21*y(4))/231.,  &
               (-4*y(-4) - 3*y(-3) - 2*y(-2) - y(-1) + y(1) + 2*y(2) + 3*y(3) + 4*y(4))/(60.*dx), &
               (28*y(-4) + 7*y(-3) - 8*y(-2) - 17*y(-1) - 20*y(0) - 17*y(1) - 8*y(2) + 7*y(3) + 28*y(4))/(462.*dx**2) /)
       case(3)
          yp(0:3) = (/ (-21*y(-4) + 14*y(-3) + 39*y(-2) + 54*y(-1) + 59*y(0) +  54*y(1) + 39*y(2) + 14*y(3) - 21*y(4))/231., &
               (86*y(-4) - 142*y(-3) - 193*y(-2) - 126*y(-1) + 126*y(1) + 193*y(2) + 142*y(3) - 86*y(4))/(1188.*dx), &
               (28*y(-4) + 7*y(-3) - 8*y(-2) - 17*y(-1) - 20*y(0) - 17*y(1) - 8*y(2) + 7*y(3) + 28*y(4))/(462.*dx**2), &
               (-14*y(-4) + 7*y(-3) + 13*y(-2) + 9*y(-1) - 9*y(1) - 13*y(2) - 7*y(3) + 14*y(4))/(198.*dx**3) /)
       case(4)
          yp(0:4) = (/ (15*y(-4) - 55*y(-3) + 30*y(-2) + 135*y(-1) + 179*y(0) + 135*y(1) + 30*y(2) - 55*y(3) + 15*y(4))/429., &
               (86*y(-4) - 142*y(-3) - 193*y(-2) - 126*y(-1) + 126*y(1) + 193*y(2) + 142*y(3) - 86*y(4))/(1188.*dx), &
               -(126*y(-4) - 371*y(-3) - 151*y(-2) + 211*y(-1) + 370*y(0) + 211*y(1) - 151*y(2) - 371*y(3) + 126*y(4))/(1716.*dx**2), &
               (-14*y(-4) + 7*y(-3) + 13*y(-2) + 9*y(-1) - 9*y(1) - 13*y(2) - 7*y(3) + 14*y(4))/(198.*dx**3), &
               (14*y(-4) - 21*y(-3) - 11*y(-2) + 9*y(-1) + 18*y(0) + 9*y(1) - 11*y(2) - 21*y(3) + 14*y(4))/(143.*dx**4) /)
       case(5)
          yp(0:5) = (/ (15*y(-4) - 55*y(-3) + 30*y(-2) + 135*y(-1) + 179*y(0) + 135*y(1) + 30*y(2) - 55*y(3) + 15*y(4))/429., &
               (-254*y(-4) + 1381*y(-3) - 2269*y(-2) - 2879*y(-1) + 2879*y(1) +  2269*y(2) - 1381*y(3) + 254*y(4))/(8580.*dx), &
               -(126*y(-4) - 371*y(-3) - 151*y(-2) + 211*y(-1) + 370*y(0) +  211*y(1) - 151*y(2) - 371*y(3) + 126*y(4))/(1716.*dx**2), &
               (100*y(-4) - 457*y(-3) + 256*y(-2) + 459*y(-1) - 459*y(1) -  256*y(2) + 457*y(3) - 100*y(4))/(1144.*dx**3), &
               (14*y(-4) - 21*y(-3) - 11*y(-2) + 9*y(-1) + 18*y(0) + 9*y(1) -  11*y(2) - 21*y(3) + 14*y(4))/(143.*dx**4), &
               (-4*y(-4) + 11*y(-3) - 4*y(-2) - 9*y(-1) + 9*y(1) + 4*y(2) - 11*y(3) + 4*y(4))/(26.*dx**5) /)
       case(6)
          yp(0:6) = (/ (-7*y(-4) + 56*y(-3) - 196*y(-2) + 392*y(-1) + 797*y(0) +392*y(1) - 196*y(2) + 56*y(3) - 7*y(4))/1287., &
               (-254*y(-4) + 1381*y(-3) - 2269*y(-2) - 2879*y(-1) + 2879*y(1) +  2269*y(2) - 1381*y(3) + 254*y(4))/(8580.*dx), &
               (11014*y(-4) - 83822*y(-3) + 250477*y(-2) - 37634*y(-1) -  280070*y(0) - 37634*y(1) + 250477*y(2) - 83822*y(3) + 11014*y(4))/(386100.*dx**2), &
               (100*y(-4) - 457*y(-3) + 256*y(-2) + 459*y(-1) - 459*y(1) - 256*y(2) + 457*y(3) - 100*y(4))/(1144.*dx**3), &
               (-136*y(-4) + 893*y(-3) - 1468*y(-2) + 11*y(-1) + 1400*y(0) + 11*y(1) - 1468*y(2) + 893*y(3) - 136*y(4))/(1170.*dx**4), &
               (-4*y(-4) + 11*y(-3) - 4*y(-2) - 9*y(-1) + 9*y(1) + 4*y(2) - 11*y(3) + 4*y(4))/(26.*dx**5), &
               (4*y(-4) - 17*y(-3) + 22*y(-2) + y(-1) - 20*y(0) + y(1) + 22*y(2) - 17*y(3) + 4*y(4))/(15.*dx**6) /)
       end select
    case(5)
       select case(m)
       case(0)
          yp(0) = sum(y(-5:5))/11.d0
       case(1)
          yp(0:1) = (/ sum(y(-5:5))/11.d0, (-5*y(-5) - 4*y(-4) - 3*y(-3) - 2*y(-2) - y(-1) + y(1) + 2*y(2) +  3*y(3) + 4*y(4) + 5*y(5))/(110.*dx) /)
       case(2)
          yp(0:2) = (/ (-36*y(-5) + 9*y(-4) + 44*y(-3) + 69*y(-2) + 84*y(-1) +  89*y(0) + 84*y(1) + 69*y(2) + 44*y(3) + 9*y(4) - 36*y(5))/429., &
               (-5*y(-5) - 4*y(-4) - 3*y(-3) - 2*y(-2) - y(-1) + y(1) + 2*y(2) +  3*y(3) + 4*y(4) + 5*y(5))/(110.*dx), &
               -(-15*y(-5) - 6*y(-4) + y(-3) + 6*y(-2) + 9*y(-1) + 10*y(0) + 9*y(1) + 6*y(2) + y(3) - 6*y(4) - 15*y(5))/(429.*dx**2) /)
       case(3)
          yp(0:3) = (/ (-36*y(-5) + 9*y(-4) + 44*y(-3) + 69*y(-2) + 84*y(-1) + 89*y(0) + 84*y(1) + 69*y(2) + 44*y(3) + 9*y(4) - 36*y(5))/429., &
               (300*y(-5) - 294*y(-4) - 532*y(-3) - 503*y(-2) - 296*y(-1) +  296*y(1) + 503*y(2) + 532*y(3) + 294*y(4) - 300*y(5))/(5148.*dx), &
               -(-15*y(-5) - 6*y(-4) + y(-3) + 6*y(-2) + 9*y(-1) + 10*y(0) +  9*y(1) + 6*y(2) + y(3) - 6*y(4) - 15*y(5))/(429.*dx**2), &
               (-30*y(-5) + 6*y(-4) + 22*y(-3) + 23*y(-2) + 14*y(-1) - 14*y(1) - 23*y(2) - 22*y(3) - 6*y(4) + 30*y(5))/(858.*dx**3) /)
       case(4)
          yp(0:4) = (/ (18*y(-5) - 45*y(-4) - 10*y(-3) + 60*y(-2) + 120*y(-1) + 143*y(0) + 120*y(1) + 60*y(2) - 10*y(3) - 45*y(4) + 18*y(5))/429., &
               (300*y(-5) - 294*y(-4) - 532*y(-3) - 503*y(-2) - 296*y(-1) + 296*y(1) + 503*y(2) + 532*y(3) + 294*y(4) - 300*y(5))/(5148.*dx), &
               (-90*y(-5) + 174*y(-4) + 146*y(-3) + y(-2) - 136*y(-1) - 190*y(0) - 136*y(1) + y(2) + 146*y(3) + 174*y(4) - 90*y(5))/(1716.*dx**2), &
               (-30*y(-5) + 6*y(-4) + 22*y(-3) + 23*y(-2) + 14*y(-1) - 14*y(1) - 23*y(2) - 22*y(3) - 6*y(4) + 30*y(5))/(858.*dx**3), &
               -(-6*y(-5) + 6*y(-4) + 6*y(-3) + y(-2) - 4*y(-1) - 6*y(0) - 4*y(1) + y(2) + 6*y(3) + 6*y(4) - 6*y(5))/(143.*dx**4) /)
       case(5)
          yp(0:5) = (/ (18*y(-5) - 45*y(-4) - 10*y(-3) + 60*y(-2) + 120*y(-1) +  143*y(0) + 120*y(1) + 60*y(2) - 10*y(3) - 45*y(4) + 18*y(5))/429., &
               (-573*y(-5) + 2166*y(-4) - 1249*y(-3) - 3774*y(-2) - 3084*y(-1) + 3084*y(1) + 3774*y(2) + 1249*y(3) - 2166*y(4) + 573*y(5))/(17160.*dx), &
               (-90*y(-5) + 174*y(-4) + 146*y(-3) + y(-2) - 136*y(-1) - 190*y(0) - 136*y(1) + y(2) + 146*y(3) + 174*y(4) - 90*y(5))/(1716.*dx**2), &
               (129*y(-5) - 402*y(-4) - 11*y(-3) + 340*y(-2) + 316*y(-1) - 316*y(1) - 340*y(2) + 11*y(3) + 402*y(4) - 129*y(5))/(2288.*dx**3), &
               -(-6*y(-5) + 6*y(-4) + 6*y(-3) + y(-2) - 4*y(-1) - 6*y(0) - 4*y(1) + y(2) + 6*y(3) + 6*y(4) - 6*y(5))/(143.*dx**4), &
               (-3*y(-5) + 6*y(-4) + y(-3) - 4*y(-2) - 4*y(-1) + 4*y(1) + 4*y(2) - y(3) - 6*y(4) + 3*y(5))/(52.*dx**5) /)
       case(6)
          yp(0:6) = (/ (-28*y(-5) + 161*y(-4) - 308*y(-3) + 28*y(-2) + 784*y(-1) + 1157*y(0) + 784*y(1) + 28*y(2) - 308*y(3) + 161*y(4) - 28*y(5))/ 2431., &
               (-573*y(-5) + 2166*y(-4) - 1249*y(-3) - 3774*y(-2) - 3084*y(-1) + 3084*y(1) + 3774*y(2) + 1249*y(3) - 2166*y(4) + 573*y(5))/(17160.*dx), &
               (68745*y(-5) - 365334*y(-4) + 540907*y(-3) + 441663*y(-2) -  320196*y(-1) - 731570*y(0) - 320196*y(1) + 441663*y(2) +  540907*y(3) - 365334*y(4) + 68745*y(5))/(2.1879e6*dx**2), &
               (129*y(-5) - 402*y(-4) - 11*y(-3) + 340*y(-2) + 316*y(-1) - 316*y(1) - 340*y(2) + 11*y(3) + 402*y(4) - 129*y(5))/(2288.*dx**3), &
               -(915*y(-5) - 4152*y(-4) + 3401*y(-3) + 3624*y(-2) - 1548*y(-1) - 4480*y(0) - 1548*y(1) + 3624*y(2) + 3401*y(3) -  4152*y(4) + 915*y(5))/(13260.*dx**4), &
               (-3*y(-5) + 6*y(-4) + y(-3) - 4*y(-2) - 4*y(-1) + 4*y(1) + 4*y(2) - y(3) - 6*y(4) + 3*y(5))/(52.*dx**5), &
               (15*y(-5) - 48*y(-4) + 29*y(-3) + 36*y(-2) - 12*y(-1) - 40*y(0) -  12*y(1) + 36*y(2) + 29*y(3) - 48*y(4) + 15*y(5))/(170.*dx**6) /)
       end select
    case(6)
       select case(m)
       case(0)
          yp(0) = sum(y(-6:6))/13.d0
       case(1)
          yp(0:1) = (/  sum(y(-6:6))/13.d0, (-6*y(-6) - 5*y(-5) - 4*y(-4) - 3*y(-3) - 2*y(-2) - y(-1) + y(1) +  2*y(2) + 3*y(3) + 4*y(4) + 5*y(5) + 6*y(6))/(182.*dx) /)
       case(2)
          yp(0:2) = (/ (-11*y(-6) + 9*y(-4) + 16*y(-3) + 21*y(-2) + 24*y(-1) + 25*y(0) + 24*y(1) + 21*y(2) + 16*y(3) + 9*y(4) - 11*y(6))/143., &
               (-6*y(-6) - 5*y(-5) - 4*y(-4) - 3*y(-3) - 2*y(-2) - y(-1) + y(1) +  2*y(2) + 3*y(3) + 4*y(4) + 5*y(5) + 6*y(6))/(182.*dx), &
               (22*y(-6) + 11*y(-5) + 2*y(-4) - 5*y(-3) - 10*y(-2) - 13*y(-1) -  14*y(0) - 13*y(1) - 10*y(2) - 5*y(3) + 2*y(4) + 11*y(5) + 22*y(6))/(1001.*dx**2) /)
       case(3)
          yp(0:3) = (/ (-11*y(-6) + 9*y(-4) + 16*y(-3) + 21*y(-2) + 24*y(-1) + 25*y(0) + 24*y(1) + 21*y(2) + 16*y(3) + 9*y(4) - 11*y(6))/143., &
               (1133*y(-6) - 660*y(-5) - 1578*y(-4) - 1796*y(-3) - 1489*y(-2) -  832*y(-1) + 832*y(1) + 1489*y(2) + 1796*y(3) + 1578*y(4) + 660*y(5) - 1133*y(6))/(24024.*dx), &
               (22*y(-6) + 11*y(-5) + 2*y(-4) - 5*y(-3) - 10*y(-2) - 13*y(-1) -  14*y(0) - 13*y(1) - 10*y(2) - 5*y(3) + 2*y(4) + 11*y(5) + 22*y(6) )/(1001.*dx**2), &
               (-11*y(-6) + 6*y(-4) + 8*y(-3) + 7*y(-2) + 4*y(-1) - 4*y(1) - 7*y(2) - 8*y(3) - 6*y(4) + 11*y(6))/(572.*dx**3) /)
       case(4)
          yp(0:4) = (/ (110*y(-6) - 198*y(-5) - 135*y(-4) + 110*y(-3) + 390*y(-2) + 600*y(-1) + 677*y(0) + 600*y(1) + 390*y(2) + 110*y(3) - 135*y(4) - 198*y(5) + 110*y(6))/2431., &
               (1133*y(-6) - 660*y(-5) - 1578*y(-4) - 1796*y(-3) - 1489*y(-2) -  832*y(-1) + 832*y(1) + 1489*y(2) + 1796*y(3) + 1578*y(4) +  660*y(5) - 1133*y(6))/(24024.*dx), &
               -(2211*y(-6) - 2970*y(-5) - 3504*y(-4) - 1614*y(-3) + 971*y(-2) + 3016*y(-1) + 3780*y(0) + 3016*y(1) + 971*y(2) - 1614*y(3) -  3504*y(4) - 2970*y(5) + 2211*y(6))/(58344.*dx**2), &
               (-11*y(-6) + 6*y(-4) + 8*y(-3) + 7*y(-2) + 4*y(-1) - 4*y(1) -  7*y(2) - 8*y(3) - 6*y(4) + 11*y(6))/(572.*dx**3), &
               (99*y(-6) - 66*y(-5) - 96*y(-4) - 54*y(-3) + 11*y(-2) + 64*y(-1) + 84*y(0) + 64*y(1) + 11*y(2) - 54*y(3) - 96*y(4) - 66*y(5) +  99*y(6))/(4862.*dx**4) /)
       case(5)
          yp(0:5) = (/ (110*y(-6) - 198*y(-5) - 135*y(-4) + 110*y(-3) + 390*y(-2) +  600*y(-1) + 677*y(0) + 600*y(1) + 390*y(2) + 110*y(3) -  135*y(4) - 198*y(5) + 110*y(6))/2431., &
               (-9647*y(-6) + 27093*y(-5) - 12*y(-4) - 33511*y(-3) - 45741*y(-2) -  31380*y(-1) + 31380*y(1) + 45741*y(2) + 33511*y(3) + 12*y(4) -  27093*y(5) + 9647*y(6))/(291720.*dx), &
               -(2211*y(-6) - 2970*y(-5) - 3504*y(-4) - 1614*y(-3) + 971*y(-2) + 3016*y(-1) + 3780*y(0) + 3016*y(1) + 971*y(2) - 1614*y(3) - 3504*y(4) - 2970*y(5) + 2211*y(6))/(58344.*dx**2), &
               (1430*y(-6) - 3267*y(-5) - 1374*y(-4) + 1633*y(-3) + 3050*y(-2) + 2252*y(-1) - 2252*y(1) - 3050*y(2) - 1633*y(3) + 1374*y(4) +  3267*y(5) - 1430*y(6))/(38896.*dx**3), &
               (99*y(-6) - 66*y(-5) - 96*y(-4) - 54*y(-3) + 11*y(-2) + 64*y(-1) +  84*y(0) + 64*y(1) + 11*y(2) - 54*y(3) - 96*y(4) - 66*y(5) +  99*y(6))/(4862.*dx**4), &
               (-22*y(-6) + 33*y(-5) + 18*y(-4) - 11*y(-3) - 26*y(-2) - 20*y(-1) + 20*y(1) + 26*y(2) + 11*y(3) - 18*y(4) - 33*y(5) + 22*y(6))/(884.*dx**5) /)
       case(6)
          yp(0:6) = (/ (-770*y(-6) + 3388*y(-5) - 3605*y(-4) - 3500*y(-3) +   4550*y(-2) + 14000*y(-1) + 18063*y(0) + 14000*y(1) + 4550*y(2) -    3500*y(3) - 3605*y(4) + 3388*y(5) - 770*y(6))/46189., &
               (-9647*y(-6) + 27093*y(-5) - 12*y(-4) - 33511*y(-3) - 45741*y(-2) -   31380*y(-1) + 31380*y(1) + 45741*y(2) + 33511*y(3) + 12*y(4) -   27093*y(5) + 9647*y(6))/(291720.*dx), &
               (482977*y(-6) - 1936330*y(-5) + 1403408*y(-4) + 2635618*y(-3) +   836377*y(-2) - 1871480*y(-1) - 3101140*y(0) - 1871480*y(1) +  836377*y(2) + 2635618*y(3) + 1403408*y(4) - 1936330*y(5) +  482977*y(6))/(1.662804e7*dx**2), &
               (1430*y(-6) - 3267*y(-5) - 1374*y(-4) + 1633*y(-3) + 3050*y(-2) +  2252*y(-1) - 2252*y(1) - 3050*y(2) - 1633*y(3) + 1374*y(4) +  3267*y(5) - 1430*y(6))/(38896.*dx**3), &
               -(2068*y(-6) - 7051*y(-5) + 2120*y(-4) + 6607*y(-3) + 2980*y(-2) - 3476*y(-1) - 6496*y(0) - 3476*y(1) + 2980*y(2) + 6607*y(3) + 2120*y(4) - 7051*y(5) + 2068*y(6))/(50388.*dx**4), &
               (-22*y(-6) + 33*y(-5) + 18*y(-4) - 11*y(-3) - 26*y(-2) - 20*y(-1) +  20*y(1) + 26*y(2) + 11*y(3) - 18*y(4) - 33*y(5) + 22*y(6))/ (884.*dx**5), &
               (22*y(-6) - 55*y(-5) + 8*y(-4) + 43*y(-3) + 22*y(-2) - 20*y(-1) - 40*y(0) - 20*y(1) + 22*y(2) + 43*y(3) +  8*y(4) - 55*y(5) + 22*y(6))/(646.*dx**6) /)
       end select
    case(7)
       select case(m)
       case(0)
          yp(0) = sum(y(-7:7))/15.d0
       case(1)
          yp(0:1) = (/  sum(y(-7:7))/15.d0 ,  &
               (-7*y(-7) - 6*y(-6) - 5*y(-5) - 4*y(-4) - 3*y(-3) - 2*y(-2) - y(-1) + y(1) + 2*y(2) + 3*y(3) + 4*y(4) + 5*y(5) + 6*y(6) +  7*y(7))/(280.*dx) /)
       case(2)
          yp(0:2) = (/ (-78*y(-7) - 13*y(-6) + 42*y(-5) + 87*y(-4) + 122*y(-3) + 147*y(-2) + 162*y(-1) + 167*y(0) + 162*y(1) + 147*y(2) + 122*y(3) + 87*y(4) + 42*y(5) - 13*y(6) - 78*y(7))/1105., &
               (-7*y(-7) - 6*y(-6) - 5*y(-5) - 4*y(-4) - 3*y(-3) - 2*y(-2) - y(-1) + y(1) + 2*y(2) + 3*y(3) + 4*y(4) + 5*y(5) + 6*y(6) + 7*y(7))/(280.*dx), &
               -(-91*y(-7) - 52*y(-6) - 19*y(-5) + 8*y(-4) + 29*y(-3) + 44*y(-2) + 53*y(-1) + 56*y(0) + 53*y(1) + 44*y(2) + 29*y(3) + 8*y(4) - 19*y(5) - 52*y(6) - 91*y(7))/(6188.*dx**2) /)
       case(3)
          yp(0:3) = (/ (-78*y(-7) - 13*y(-6) + 42*y(-5) + 87*y(-4) + 122*y(-3) + 147*y(-2) + 162*y(-1) + 167*y(0) + 162*y(1) + 147*y(2) + 122*y(3) + 87*y(4) + 42*y(5) - 13*y(6) - 78*y(7))/1105., &
               (12922*y(-7) - 4121*y(-6) - 14150*y(-5) - 18334*y(-4) -  17842*y(-3) - 13843*y(-2) - 7506*y(-1) + 7506*y(1) + 13843*y(2) + 17842*y(3) + 18334*y(4) + 14150*y(5) + 4121*y(6) -  12922*y(7))/(334152.*dx), &
               -(-91*y(-7) - 52*y(-6) - 19*y(-5) + 8*y(-4) + 29*y(-3) + 44*y(-2) + 53*y(-1) + 56*y(0) + 53*y(1) + 44*y(2) + 29*y(3) + 8*y(4) - 19*y(5) - 52*y(6) - 91*y(7))/(6188.*dx**2), &
               (-91*y(-7) - 13*y(-6) + 35*y(-5) + 58*y(-4) + 61*y(-3) + 49*y(-2) + 27*y(-1) - 27*y(1) - 49*y(2) - 61*y(3) - 58*y(4) - 35*y(5) + 13*y(6) + 91*y(7))/(7956.*dx**3) /) 
       case(4)
          yp(0:4) = (/ (2145*y(-7) - 2860*y(-6) - 2937*y(-5) - 165*y(-4) +  3755*y(-3) + 7500*y(-2) + 10125*y(-1) + 11063*y(0) + 10125*y(1) + 7500*y(2) + 3755*y(3) - 165*y(4) - 2937*y(5) - 2860*y(6) + 2145*y(7))/46189., &
               (12922*y(-7) - 4121*y(-6) - 14150*y(-5) - 18334*y(-4) - 17842*y(-3) - 13843*y(-2) - 7506*y(-1) + 7506*y(1) + 13843*y(2) + 17842*y(3) + 18334*y(4) + 14150*y(5) + 4121*y(6) - 12922*y(7))/(334152.*dx), &
               (-31031*y(-7) + 29601*y(-6) + 44495*y(-5) + 31856*y(-4) + 6579*y(-3) - 19751*y(-2) - 38859*y(-1) - 45780*y(0) - 38859*y(1) - 19751*y(2) + 6579*y(3) + 31856*y(4) + 44495*y(5) + 29601*y(6) - 31031*y(7))/(1.108536e6*dx**2), &
               (-91*y(-7) - 13*y(-6) + 35*y(-5) + 58*y(-4) + 61*y(-3) + 49*y(-2) + 27*y(-1) - 27*y(1) - 49*y(2) - 61*y(3) - 58*y(4) - 35*y(5) +  13*y(6) + 91*y(7))/(7956.*dx**3), &
               (1001*y(-7) - 429*y(-6) - 869*y(-5) - 704*y(-4) - 249*y(-3) + 251*y(-2) + 621*y(-1) + 756*y(0) + 621*y(1) + 251*y(2) - 249*y(3) - 704*y(4) - 869*y(5) - 429*y(6) + 1001*y(7))/(92378.*dx**4) /)
       case(5)
          yp(0:5) = (/ (2145*y(-7) - 2860*y(-6) - 2937*y(-5) - 165*y(-4) +  3755*y(-3) + 7500*y(-2) + 10125*y(-1) + 11063*y(0) +  10125*y(1) + 7500*y(2) + 3755*y(3) - 165*y(4) - 2937*y(5) -  2860*y(6) + 2145*y(7))/46189., &
               (-78351*y(-7) + 169819*y(-6) + 65229*y(-5) - 130506*y(-4) -   266401*y(-3) - 279975*y(-2) - 175125*y(-1) + 175125*y(1) +   279975*y(2) + 266401*y(3) + 130506*y(4) - 65229*y(5) -  169819*y(6) + 78351*y(7))/(2.5194e6*dx), &
               (-31031*y(-7) + 29601*y(-6) + 44495*y(-5) + 31856*y(-4) +  6579*y(-3) - 19751*y(-2) - 38859*y(-1) - 45780*y(0) -  38859*y(1) - 19751*y(2) + 6579*y(3) + 31856*y(4) + 44495*y(5) +  29601*y(6) - 31031*y(7))/(1.108536e6*dx**2), &
               (8281*y(-7) - 14404*y(-6) - 10379*y(-5) + 1916*y(-4) + 11671*y(-3) + 14180*y(-2) + 9315*y(-1) - 9315*y(1) - 14180*y(2) - 11671*y(3) - 1916*y(4) + 10379*y(5) + 14404*y(6) - 8281*y(7))/(335920.*dx**3), &
               (1001*y(-7) - 429*y(-6) - 869*y(-5) - 704*y(-4) - 249*y(-3) +  251*y(-2) + 621*y(-1) + 756*y(0) + 621*y(1) + 251*y(2) -  249*y(3) - 704*y(4) - 869*y(5) - 429*y(6) + 1001*y(7))/ (92378.*dx**4), &
               (-1001*y(-7) + 1144*y(-6) + 979*y(-5) + 44*y(-4) -  751*y(-3) - 1000*y(-2) - 675*y(-1) + 675*y(1) + 1000*y(2) + 751*y(3) - 44*y(4) - 979*y(5) - 1144*y(6) + 1001*y(7))/(83980.*dx**5) /)
       case(6)
          yp(0:6) = (/ (-260*y(-7) + 910*y(-6) - 476*y(-5) - 1085*y(-4) - 140*y(-3) +   1750*y(-2) + 3500*y(-1) + 4199*y(0) + 3500*y(1) + 1750*y(2) -  140*y(3) - 1085*y(4) - 476*y(5) + 910*y(6) - 260*y(7))/12597., &
               (-78351*y(-7) + 169819*y(-6) + 65229*y(-5) - 130506*y(-4) -  266401*y(-3) - 279975*y(-2) - 175125*y(-1) + 175125*y(1) +     279975*y(2) + 266401*y(3) + 130506*y(4) - 65229*y(5) -     169819*y(6) + 78351*y(7))/(2.5194e6*dx), &
               (573118*y(-7) - 1810211*y(-6) + 445570*y(-5) + 2138176*y(-4) +  1798522*y(-3) + 18325*y(-2) - 1850650*y(-1) - 2625700*y(0) -  1850650*y(1) + 18325*y(2) + 1798522*y(3) + 2138176*y(4) +  445570*y(5) - 1810211*y(6) + 573118*y(7))/(2.26746e7*dx**2), &
               (8281*y(-7) - 14404*y(-6) - 10379*y(-5) + 1916*y(-4) +  11671*y(-3) + 14180*y(-2) + 9315*y(-1) - 9315*y(1) -  14180*y(2) - 11671*y(3) - 1916*y(4) + 10379*y(5) + 14404*y(6) -  8281*y(7))/(335920.*dx**3), &
               -(19019*y(-7) - 50908*y(-6) - 3355*y(-5) + 39248*y(-4) +   39521*y(-3) + 7460*y(-2) - 28865*y(-1) - 44240*y(0) -   28865*y(1) + 7460*y(2) + 39521*y(3) + 39248*y(4) - 3355*y(5) -   50908*y(6) + 19019*y(7))/(755820.*dx**4), &
               (-1001*y(-7) + 1144*y(-6) + 979*y(-5) + 44*y(-4) - 751*y(-3) -  1000*y(-2) - 675*y(-1) + 675*y(1) + 1000*y(2) + 751*y(3) -  44*y(4) - 979*y(5) - 1144*y(6) + 1001*y(7))/(83980.*dx**5), &
               (143*y(-7) - 286*y(-6) - 55*y(-5) + 176*y(-4) + 197*y(-3) +  50*y(-2) - 125*y(-1) - 200*y(0) - 125*y(1) + 50*y(2) +  197*y(3) + 176*y(4) - 55*y(5) - 286*y(6) + 143*y(7))/(9690.*dx**6) /)
       end select
    case default
       stop "numderv_polyfit: n too big"
    end select
  end subroutine numderv_polyfit


End Module General_utils


