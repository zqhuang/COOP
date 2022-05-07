!!this module does svd decomposition and least square fit
module coop_svd_mod
  use coop_constants_mod
  implicit none
#include "constants.h"

  private

  public:: coop_svd_least_square_one, coop_svd_decompose, coop_svd_least_square_all, coop_svd_decompose_invw, coop_svd_sol_single, coop_svd_sol_multiple

contains

  Function Coop_outerprod(A,B)
    COOP_REAL, DIMENSION(:), INTENT(IN) :: a,b
    COOP_REAL, DIMENSION(size(a),size(b)) :: coop_outerprod
    Coop_outerprod = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  end function Coop_outerprod
  
  Function coop_pythag(a,b)
    COOP_REAL, INTENT(IN) :: a,b
    COOP_REAL :: coop_pythag
    COOP_REAL :: absa,absb
    absa=abs(a)
    absb=abs(b)
    if (absa > absb) then
       coop_pythag=absa*sqrt(1.+(absb/absa)**2)
    else
       if (absb == 0.0) then
          coop_pythag=0.0
       else
          coop_pythag=absb*sqrt(1.0+(absa/absb)**2)
       end if
    end if
  End function coop_pythag


  !! SVD decomposition a = U diag(w) V^T, the output is w and V, a is replaced with U
  recursive subroutine coop_svd_decompose(m, n, a, w, v)
    COOP_INT,parameter::maxits = 30
    COOP_INT m, n
    COOP_REAL, INTENT(INOUT) :: a(:, :)
    COOP_REAL, INTENT(OUT) :: w(:)
    COOP_REAL, INTENT(OUT) :: v(:, :)
    COOP_INT :: i,its,j,k,l,nm
    COOP_REAL :: anorm,c,f,g,h,s,scale,x,y,z
    COOP_REAL tempm(m), rv1(n),tempn(n)
#if HAS_LAPACK
    COOP_REAL:: u(m, m)
    COOP_INT :: info, lwork
    COOP_INT, dimension(:),allocatable::iwork
    COOP_REAL, dimension(:),allocatable::work
    lwork = (min(M,N)*(6+4*min(M,N))+max(M,N))*2
    allocate(work(lwork), iwork(8*min(m, n)))
    call dgesdd("A", m, n, a,  m, w, u, m, v, n, work, lwork, iwork, info)
    deallocate(work, iwork)
    if(info .gt. 0)then
       write(*,*) "Warning: svd decomposition failed to converge in LAPACK"
       return
    endif
    if(info .lt. 0)then
       print*, "the ", -info, " th argument in svd decomposition is wrong"
       stop
    endif
    a(1:m, 1:min(n,m)) = u(1:m, 1:min(m,n))
    v = transpose(v)
    return
#endif
    if(size(a, 1) .ne. m .or. size(a, 2) .lt. max(m, n) .or. size(w) .lt. max(m, n) .or. size(v,1).ne. n .or. size(v, 2) .lt. n) stop "svd_decompose; find fire"
    if(m .lt. n)then
       v = 0.d0
       v(1:n, 1:m) = transpose(a(1:m, 1:n))
       call coop_svd_decompose(n, m, v, w, a(1:m, 1:m))
       return
    endif
    g=0.0
    scale=0.0
    do i=1,n
       l=i+1
       rv1(i)=scale*g
       g=0.0
       scale=0.0
       if (i <= m) then
          scale=sum(abs(a(i:m,i)))
          if (scale /= 0.0) then
             a(i:m,i)=a(i:m,i)/scale
             s=dot_product(a(i:m,i),a(i:m,i))
             f=a(i,i)
             g=-sign(sqrt(s),f)
             h=f*g-s
             a(i,i)=f-g
             tempn(l:n)=matmul(a(i:m,i),a(i:m,l:n))/h
             a(i:m,l:n)=a(i:m,l:n)+coop_outerprod(a(i:m,i),tempn(l:n))
             a(i:m,i)=scale*a(i:m,i)
          end if
       end if
       w(i)=scale*g
       g=0.0
       scale=0.0
       if ((i <= m) .and. (i /= n)) then
          scale=sum(abs(a(i,l:n)))
          if (scale /= 0.0) then
             a(i,l:n)=a(i,l:n)/scale
             s=dot_product(a(i,l:n),a(i,l:n))
             f=a(i,l)
             g=-sign(sqrt(s),f)
             h=f*g-s
             a(i,l)=f-g
             rv1(l:n)=a(i,l:n)/h
             tempm(l:m)=matmul(a(l:m,l:n),a(i,l:n))
             a(l:m,l:n)=a(l:m,l:n)+coop_outerprod(tempm(l:m),rv1(l:n))
             a(i,l:n)=scale*a(i,l:n)
          end if
       end if
    end do
    anorm=maxval(abs(w)+abs(rv1))
    do i=n,1,-1
       if (i < n) then
          if (g /= 0.0) then
             v(l:n,i)=(a(i,l:n)/a(i,l))/g
             tempn(l:n)=matmul(a(i,l:n),v(l:n,l:n))
             v(l:n,l:n)=v(l:n,l:n)+coop_outerprod(v(l:n,i),tempn(l:n))
          end if
          v(i,l:n)=0.0
          v(l:n,i)=0.0
       end if
       v(i,i)=1.0
       g=rv1(i)
       l=i
    end do
    do i=min(m,n),1,-1
       l=i+1
       g=w(i)
       a(i,l:n)=0.0
       if (g /= 0.0) then
          g=1.0/g
          tempn(l:n)=(matmul(a(l:m,i),a(l:m,l:n))/a(i,i))*g
          a(i:m,l:n)=a(i:m,l:n)+coop_outerprod(a(i:m,i),tempn(l:n))
          a(i:m,i)=a(i:m,i)*g
       else
          a(i:m,i)=0.0
       end if
       a(i,i)=a(i,i)+1.0
    end do
    do k=n,1,-1
       do its=1,maxits
          do l=k,1,-1
             nm=l-1
             if ((abs(rv1(l))+anorm) == anorm) exit
             if ((abs(w(nm))+anorm) == anorm) then
                c=0.0
                s=1.0
                do i=l,k
                   f=s*rv1(i)
                   rv1(i)=c*rv1(i)
                   if ((abs(f)+anorm) == anorm) exit
                   g=w(i)
                   h=coop_pythag(f,g)
                   w(i)=h
                   h=1.0/h
                   c= (g*h)
                   s=-(f*h)
                   tempm(1:m)=a(1:m,nm)
                   a(1:m,nm)=a(1:m,nm)*c+a(1:m,i)*s
                   a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
                end do
                exit
             end if
          end do
          z=w(k)
          if (l == k) then
             if (z < 0.0) then
                w(k)=-z
                v(1:n,k)=-v(1:n,k)
             end if
             exit
          end if
          if (its == maxits)then
             write(*,*) 'error: no convergence in coop_svd_decompose'
             stop
          endif
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
          g=coop_pythag(f,COOP_REAL_OF(1.0))
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0
          s=1.0
          do j=l,nm
             i=j+1
             g=rv1(i)
             y=w(i)
             h=s*g
             g=c*g
             z=coop_pythag(f,h)
             rv1(j)=z
             c=f/z
             s=h/z
             f= (x*c)+(g*s)
             g=-(x*s)+(g*c)
             h=y*s
             y=y*c
             tempn(1:n)=v(1:n,j)
             v(1:n,j)=v(1:n,j)*c+v(1:n,i)*s
             v(1:n,i)=-tempn(1:n)*s+v(1:n,i)*c
             z=coop_pythag(f,h)
             w(j)=z
             if (z /= 0.0) then
                z=1.0/z
                c=f*z
                s=h*z
             end if
             f= (c*g)+(s*y)
             x=-(s*g)+(c*y)
             tempm(1:m)=a(1:m,j)
             a(1:m,j)=a(1:m,j)*c+a(1:m,i)*s
             a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
          end do
          rv1(l)=0.0
          rv1(k)=f
          w(k)=x
       end do
    end do
  end subroutine coop_svd_decompose


  !!SVD decomposition a = U diag(1/w) V^T, the output is w and V, a is replaced with U
  Subroutine coop_svd_decompose_invw(m, n, a, w, v)
    COOP_REAL,parameter::wtol = 1.d-7
    COOP_INT m, n
    COOP_REAL, INTENT(INOUT) :: a(m, n)
    COOP_REAL, INTENT(OUT) :: w(n)
    COOP_REAL, INTENT(OUT) :: v(n, n)
    COOP_REAL anorm
    call coop_svd_decompose(m,n,a,w,v)
    anorm = wtol*maxval(w)
    where(w .lt. anorm)
       w = 0.d0
    elsewhere
       w = 1.d0/w
    end where
  End Subroutine coop_svd_decompose_invw

  Subroutine coop_svd_sol_single(m, n, u, w, v, b, x)
    COOP_INT m, n
    COOP_REAL, INTENT(IN) :: u(m, n), w(n), v(n, n), b(m)
    COOP_REAL, INTENT(OUT) :: x(n)
    x=matmul(v,matmul(b,u)*w)
  End Subroutine coop_svd_sol_single

  subroutine coop_svd_sol_multiple(m, n, nd, u, w, v, b, x)
    COOP_INT m, n, nd
    COOP_REAL, INTENT(IN) :: u(m, n), w(n), v(n, n),  b(m, nd)
    COOP_REAL, INTENT(OUT) :: x(n, nd)
    COOP_INT i
    do i=1,nd
       x(:, i) =matmul(v,matmul(b(:,i), u)*w)
    enddo
  end subroutine coop_svd_sol_multiple

  subroutine coop_svd_least_square_one(m, n, A, b, x)
    COOP_INT m, n
    COOP_REAL,intent(INOUT)::a(m, n)
    COOP_REAL,intent(IN)::b(m)
    COOP_REAL,intent(OUT)::x(n)
    COOP_REAL w(n), v(n, n)
    call coop_svd_decompose_invw(m, n, A, w, v)
    call coop_svd_sol_single(m, n, A, w, v, b, x)
  end subroutine coop_svd_least_square_one

  subroutine coop_svd_least_square_all(m, n, nd, A, b, x)
    COOP_INT m,n, nd
    COOP_REAL,intent(INOUT)::a(m, n)
    COOP_REAL,intent(IN)::b(m, nd)
    COOP_REAL,intent(OUT)::x(n, nd)
    COOP_REAL w(n), v(n, n)
    call coop_svd_decompose_invw(m, n, A, w, v)
    call coop_svd_sol_multiple(m, n, nd, A, w, v, b, x)
  end subroutine coop_svd_least_square_all

end module coop_svd_mod
