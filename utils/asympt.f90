module asymp_utils
  use general_utils
  implicit none

  
  
contains

  subroutine asymp_splint(y,z,n)
!  z = integral from 1 to infinity of s(i)di, where s(i) is the spline fit
!  to y(i).
    integer, intent(in) :: n
    real(dl), intent(in) :: y(n)
    real(dl), intent(out) :: z
    real(dl) :: dy1, dyn
   ! dy1=(-11._dl*y(1)+18._dl*y(2)-9._dl*y(3)+2._dl*y(4))/6._dl
    dy1=0._dl
    dyn=(11._dl*y(n)-18._dl*y(n-1)+9._dl*y(n-2)-2._dl*y(n-3))/6._dl
    z=0.5d0*(y(1)+y(n))+(dy1-dyn)/12._dl
    z= (z + sum(y(2:n-1)))
  end subroutine asymp_splint

  
  function asymp_sum(a) result(s)
    integer n,  ntail
    real(dl) s, a(:), ss, ssmax, err1, err2, rat, alpha
    s = sum(a)
    n = size(a)
    if(any(a.lt.0) .or. n.lt.100)return  !!only for positive arrays we do extrapolation
    ntail = 1
    ssmax = s * 0.03
    ss = a(n)
    do while(ntail .lt. n/5 )
       ss = ss + a(n-ntail)
       if(ss .ge. ssmax)exit
       ntail = ntail + 1
    enddo
    if(ntail .le. 10) return  !!error too big so we do not do any corrections (insead of doing likely wrong ones...)
    if(mod(ntail,2).ne.0) ntail = ntail + 1
    err1 = sum(a(n-ntail+1: n-ntail/2))
    err2 = sum(a(n-ntail/2+1:n))
    if(err2 .lt. ss*1.d-4) return
    rat = err1/err2
    !!if rat ~ 1,  assume a power-law
    alpha = (rat - 1.d0)*(n*2.d0-ntail)/ntail    
    if(alpha .le. 1.1d0) return !!not converging
    if(alpha .gt. 5.d0) then !!large index, maybe exponential is better
       s = s + min(err2/(rat - 1.d0), ssmax)
    else
       s = s + min(err2/(alpha-1)*2.d0*n/ntail*(1.d0-ntail/2.d0/n)**alpha, ssmax)
    endif
  end function asymp_sum
  

end module asymp_utils
