module Gaussianfit_utils
  use random_utils
  use interpolation_utils
  implicit none

  integer,parameter::gaussianfit_maxnum = 4
  integer,parameter::ggfit_maxnum = 3
!!fit a positive function with the sum of n Gaussians
  type gaussianfit_vars
     integer n
     real(dl),dimension(gaussianfit_maxnum)::s,c,h
  end type gaussianfit_vars


  type ggfit_vars
     integer n
     real(dl),dimension(ggfit_maxnum)::s,c,h,a
  end type ggfit_vars

  interface ggfit_tilted
     module procedure ggfit_tilted_s, ggfit_tilted_vx, ggfit_tilted_va, ggfit_tilted_vax
  end interface ggfit_tilted

  interface ggfit_super_tilted
     module procedure ggfit_super_tilted_s, ggfit_super_tilted_vx, ggfit_super_tilted_va, ggfit_super_tilted_vax
  end interface ggfit_super_tilted


  interface gaussianfit_eval
     module procedure gaussianfit_eval_s, gaussianfit_eval_v
  end interface gaussianfit_eval


  interface ggfit_eval
     module procedure ggfit_eval_s, ggfit_eval_v
  end interface ggfit_eval

contains


  subroutine gaussianfit_get_heights(n, x, y, gvar)
    integer,intent(IN)::n
    real(dl),intent(IN)::x(n),y(n)
    type(gaussianfit_vars)::gvar
    real(dl) fx(n, gvar%n)
    integer i
    do i=1,gvar%n
       fx(:, i) = exp(-((x-gvar%c(i))/gvar%s(i))**2/2.d0)
    enddo
    call fit_template(n, gvar%n, y, fx, gvar%h(1:gvar%n))
  end subroutine gaussianfit_get_heights


  subroutine ggfit_get_heights(n, x, y, gvar)
    integer,intent(IN)::n
    real(dl),intent(IN)::x(n),y(n)
    type(ggfit_vars)::gvar
    real(dl) fx(n, gvar%n)
    integer i
    do i = 1, gvar%n
       fx(:, i) = exp(- &
            ( ggfit_tilted(gvar%a(i), x - gvar%c(i) )/gvar%s(i) )**2/2.d0 &
            )
    enddo
    call fit_template(n, gvar%n, y, fx, gvar%h(1:gvar%n))
  end subroutine ggfit_get_heights



  subroutine gaussianfit_eval_s(gvar, x, y)
    real(dl),intent(IN)::x
    real(dl),intent(OUT)::y
    type(gaussianfit_vars)::gvar
    y = sum(gvar%h(1:gvar%n)*exp(-((x-gvar%c(1:gvar%n))/gvar%s(1:gvar%n))**2/2.d0))
  end subroutine gaussianfit_eval_s


  subroutine gaussianfit_eval_v(gvar, x, y)
    real(dl),dimension(:),intent(IN)::x
    real(dl),dimension(:),intent(OUT)::y
    type(gaussianfit_vars)::gvar
    integer i
    do i=1,size(x)
       y(i) = sum(gvar%h(1:gvar%n)*exp(-((x(i)-gvar%c(1:gvar%n))/gvar%s(1:gvar%n))**2/2.d0))
    enddo
  end subroutine gaussianfit_eval_v


  subroutine ggfit_eval_s(gvar, x, y)
    real(dl),intent(IN)::x
    real(dl),intent(OUT)::y
    type(ggfit_vars)::gvar
    y = sum( gvar%h(1:gvar%n) * exp( &
         -  (ggfit_tilted(gvar%a(1:gvar%n), x-gvar%c(1:gvar%n))/gvar%s(1:gvar%n))**2/2.d0 &
         ))
  end subroutine ggfit_eval_s


  subroutine ggfit_eval_v(gvar, x, y)
    real(dl),dimension(:),intent(IN)::x
    real(dl),dimension(:),intent(OUT)::y
    type(ggfit_vars)::gvar
    integer i
    do i=1,size(x)
       y(i) = sum( gvar%h(1:gvar%n) * exp( &
         -  (ggfit_tilted(gvar%a(1:gvar%n), x(i)-gvar%c(1:gvar%n))/gvar%s(1:gvar%n))**2/2.d0 &
         ))
    enddo
  end subroutine ggfit_eval_v

  subroutine gaussianfit_get_diff(n, x, y, gvar, diff, cderv, sderv)
    integer,intent(IN)::n
    real(dl),intent(IN)::x(n),y(n)
    type(gaussianfit_vars),intent(IN):: gvar
    real(dl),intent(OUT)::diff
    real(dl),intent(OUT)::cderv(gvar%n), sderv(gvar%n)
    real(dl) edf(gvar%n), tdf
    integer i
    diff = 0.d0
    cderv = 0.d0
    sderv = 0.d0
    do i=1, n
       edf = gvar%h(1:gvar%n)*exp(-((x(i)-gvar%c(1:gvar%n))/gvar%s(1:gvar%n))**2/2.d0) 
       tdf = sum(edf) - y(i)
       diff = diff + tdf**2
       edf = 2.d0 * tdf * edf * (x(i) - gvar%c(1:gvar%n))/gvar%s(1:gvar%n)**2
       cderv = cderv + edf
       sderv = sderv + edf * (x(i)-gvar%c(1:gvar%n))/gvar%s(1:gvar%n)
    enddo
  end subroutine gaussianfit_get_diff


  subroutine ggfit_get_diff(n, x, y, gvar, diff, cderv, sderv, aderv)
    integer,intent(IN)::n
    real(dl),intent(IN)::x(n),y(n)
    type(ggfit_vars),intent(IN):: gvar
    real(dl),intent(OUT)::diff
    real(dl),intent(OUT)::cderv(gvar%n), sderv(gvar%n),aderv(gvar%n)
    real(dl) edf(gvar%n), tdf
    integer i
    diff = 0.d0
    cderv = 0.d0
    sderv = 0.d0
    aderv = 0.d0
    do i=1, n
       edf = gvar%h(1:gvar%n)*exp(-(ggfit_tilted(gvar%a(1:gvar%n), x(i)-gvar%c(1:gvar%n))/gvar%s(1:gvar%n))**2/2.d0) 
       tdf = sum(edf) - y(i)
       diff = diff + tdf**2
       edf = 2.d0 * tdf * edf * ggfit_tilted(gvar%a(1:gvar%n), x(i) - gvar%c(1:gvar%n))/gvar%s(1:gvar%n)**2
       cderv = cderv + edf * dexp(gvar%a(1:gvar%n)*(x(i) - gvar%c(1:gvar%n)))
       sderv = sderv + edf * ggfit_tilted(gvar%a(1:gvar%n), x(i) - gvar%c(1:gvar%n))/gvar%s(1:gvar%n)
       aderv = aderv + edf  * ggfit_super_tilted(gvar%a(1:gvar%n), x(i) - gvar%c(1:gvar%n))
    enddo
  end subroutine ggfit_get_diff




  subroutine gaussianfit_any(n, x, y, m, gvar)
    integer,parameter::maxfails = 49
    integer,parameter::maxnsteps = 80000
    integer,intent(IN)::n, m
    real(dl),intent(IN)::x(n),y(n)
    type(gaussianfit_vars),intent(OUT)::gvar
    real(dl) sy, csave, hsave(m), ssave
    integer i, nfail, nsteps
    real(dl) diff, initdiff, newdiff, cderv(m), sderv(m), newcderv(m), newsderv(m), step
    call random_init()
    call gaussianfit_alloc(gvar, m)
    step = 0.05d0+0.08/m
100 continue    
    sy = sum(y)
    gvar%c(1:m) = sum(x*y)/sy
    gvar%s(1) = sqrt(sum((x-gvar%c(1))**2*y)/sy)
    do i=2,m
       gvar%s(i) = gvar%s(1) * i
    enddo
    call gaussianfit_get_heights(n, x, y, gvar)
    call gaussianfit_get_diff(n, x, y, gvar, diff, cderv, sderv)    
    initdiff = diff
    nfail = 0
    hsave = gvar%h(1:m)
    nsteps = 0
    do while(diff .gt. initdiff*1.d-6 .and. nfail .lt. maxfails)
       i = random_index(m)
       csave = gvar%c(i)
       ssave = gvar%s(i)          
       nfail = nfail + 1
       if(random_unit().lt.0.5d0)then !!do c
          gvar%c(i) = gvar%c(i) -  sign(min(diff/(abs(cderv(i))+1.d-99), gvar%s(i)), cderv(i))*random_unit()*step/(1+nfail)
       else !!do s
          gvar%s(i) = gvar%s(i) - sign(min(diff/(abs(sderv(i))+1.d-99), gvar%s(i)), sderv(i))*random_unit()*step/(1+nfail)
       endif
       call gaussianfit_get_heights(n, x, y, gvar)
       call gaussianfit_get_diff(n, x, y, gvar, newdiff, newcderv, newsderv)
       if(newdiff .lt. diff)then
          diff = newdiff
          cderv = newcderv
          sderv = newsderv
          hsave = gvar%h(1:m)
          nfail = 0
       else
          gvar%c(i) = csave
          gvar%s(i) = ssave
          gvar%h(1:m) = hsave
          nfail = nfail + 1
       endif
       nsteps = nsteps + 1
       if(nsteps .gt. maxnsteps)then
          nsteps = 0
          step = step/2.d0
          if(step .gt. 0.002d0)then
             write(*,*) "warning: gaussianfit retrying!"
             goto 100
          else
             write(*,*) "warning: gaussianfit may have failed!"
             return
          endif
       endif
    enddo
  end subroutine gaussianfit_any

  subroutine gaussianfit_alloc(gvar, n)
    integer,intent(IN)::n
    type(gaussianfit_vars),intent(INOUT)::gvar
    if(.not.(n.gt.0 .and. n .le. gaussianfit_maxnum))then
       write(*,*) "Overflow in Gaussianfit: number of gaussians = ", gvar%n
       stop 
    endif
    gvar%n = n
  end subroutine gaussianfit_alloc


  subroutine ggfit_alloc(gvar, n)
    integer,intent(IN)::n
    type(ggfit_vars),intent(INOUT)::gvar
    if(.not.(n.gt.0 .and. n .le. ggfit_maxnum))then
       write(*,*) "Overflow in Ggfit: number of ggs = ", gvar%n
       stop 
    endif
    gvar%n = n
  end subroutine ggfit_alloc


  subroutine ggfit_get_tilted(a, x, y)
    real(dl) x, a, y
    y = a * x
    if(abs(y) .lt. 1.d-3)then
       y = x * ( 1.d0 + y*(0.5d0 + y*(1.d0/6.d0+ y/24.d0)))
    else
       y = (dexp(y)-1.d0)/a
    endif
    
  end subroutine ggfit_get_tilted

  function ggfit_tilted_s(a, x) result(y)
    real(dl) x, a, y
    call ggfit_get_tilted(a, x, y)
  end function ggfit_tilted_s

  function ggfit_tilted_vx(a, x) result(y)
    real(dl),intent(IN):: a
    real(dl),dimension(:),intent(IN):: x 
    real(dl) y(size(x))
    integer i
    do i=1,size(x)
       call ggfit_get_tilted(a, x(i), y(i))
    enddo
  end function ggfit_tilted_vx

  function ggfit_tilted_va(a, x) result(y)
    real(dl),intent(IN):: x
    real(dl),dimension(:),intent(IN):: a
    real(dl) y(size(a))
    integer i
    do i=1,size(a)
       call ggfit_get_tilted(a(i), x, y(i))
    enddo
  end function ggfit_tilted_va


  function ggfit_tilted_vax(a, x) result(y)
    real(dl),dimension(:),intent(IN):: x
    real(dl),dimension(:),intent(IN):: a
    real(dl) y(size(x))
    integer i
    do i=1,size(x)
       call ggfit_get_tilted(a(i), x(i), y(i))
    enddo
  end function ggfit_tilted_vax


  subroutine ggfit_get_super_tilted(a, x, y)
    real(dl) x, a, y
    y = a * x
    if(abs(y) .lt. 1.d-2)then
       y = x **2 * ( -0.5d0 + y*( -1.d0/3.d0 + y*( -1.d0/8.d0 - y/30.d0)))
    else
       y = (dexp(y)*(1-y)-1.d0)/a**2
    endif
    
  end subroutine ggfit_get_super_tilted

  function ggfit_super_tilted_s(a, x) result(y)
    real(dl) x, a, y
    call ggfit_get_super_tilted(a, x, y)
  end function ggfit_super_tilted_s

  function ggfit_super_tilted_vx(a, x) result(y)
    real(dl),intent(IN):: a
    real(dl),dimension(:),intent(IN):: x 
    real(dl) y(size(x))
    integer i
    do i=1,size(x)
       call ggfit_get_super_tilted(a, x(i), y(i))
    enddo
  end function ggfit_super_tilted_vx

  function ggfit_super_tilted_va(a, x) result(y)
    real(dl),intent(IN):: x
    real(dl),dimension(:),intent(IN):: a
    real(dl) y(size(a))
    integer i
    do i=1,size(a)
       call ggfit_get_super_tilted(a(i), x, y(i))
    enddo
  end function ggfit_super_tilted_va


  function ggfit_super_tilted_vax(a, x) result(y)
    real(dl),dimension(:),intent(IN):: x
    real(dl),dimension(:),intent(IN):: a
    real(dl) y(size(x))
    integer i
    do i=1,size(x)
       call ggfit_get_super_tilted(a(i), x(i), y(i))
    enddo
  end function ggfit_super_tilted_vax


  subroutine ggfit_any(n, x, y, m, gvar)
    integer,parameter::maxfails = 99
    integer,parameter::maxnsteps = 30000
    integer,intent(IN)::n, m
    real(dl),intent(IN)::x(n),y(n)
    type(ggfit_vars),intent(OUT)::gvar
    real(dl) sy, csave, hsave(m), ssave, asave
    integer i, nfail, nsteps, j
    real(dl) diff, initdiff, newdiff, cderv(m), sderv(m), aderv(m), newcderv(m), newsderv(m), newaderv(m), step
    call random_init()
    call ggfit_alloc(gvar, m)
    step = 0.1d0+0.05/m
100 continue    
    sy = sum(y)
    gvar%c(1:m) = sum(x*y)/sy
    gvar%s(1) = sqrt(sum((x-gvar%c(1))**2*y)/sy)
    do i=2,m
       gvar%s(i) = gvar%s(1) * i
    enddo
    gvar%a(1:m) = 0.d0
    call ggfit_get_heights(n, x, y, gvar)
    call ggfit_get_diff(n, x, y, gvar, diff, cderv, sderv, aderv)    
    initdiff = diff
    nfail = 0
    hsave = gvar%h(1:m)
    nsteps = 0
    do while(diff .gt. initdiff*1.d-6 .and. nfail .lt. maxfails)
       i = random_index(m)
       csave = gvar%c(i)
       ssave = gvar%s(i)          
       asave = gvar%a(i)
       nfail = nfail + 1
       j = random_index(3)
       select case(j)
       case(1)
          gvar%c(i) = gvar%c(i) -  sign(min(diff/(abs(cderv(i))+1.d-99), gvar%s(i)), cderv(i))*random_unit()*step/(1+nfail)
       case(2)
          gvar%s(i) = gvar%s(i) - sign(min(diff/(abs(sderv(i))+1.d-99), gvar%s(i)), sderv(i))*random_unit()*step/(1+nfail)
       case(3)
          gvar%a(i) = gvar%a(i) - sign(min(diff/(abs(aderv(i))+1.d-99), 0.2d0/gvar%s(i)), aderv(i))*random_unit()*step/(1+nfail)
       end select
       call ggfit_get_heights(n, x, y, gvar)
       call ggfit_get_diff(n, x, y, gvar, newdiff, newcderv, newsderv, newaderv)
       if(newdiff .lt. diff)then
          diff = newdiff
          cderv = newcderv
          sderv = newsderv
          aderv = newaderv
          hsave = gvar%h(1:m)
          nfail = 0
       else
          gvar%c(i) = csave
          gvar%s(i) = ssave
          gvar%a(i) = asave
          gvar%h(1:m) = hsave
          nfail = nfail + 1
       endif
       nsteps = nsteps + 1
       if(nsteps .gt. maxnsteps)then
          nsteps = 0
          step = step/2.d0
          if(step .gt. 0.002d0)then
             write(*,*) "warning: ggfit retrying!"
             goto 100
          else
             write(*,*) "warning: ggfit may have failed!"
             return
          endif
       endif
    enddo
  end subroutine ggfit_any

  subroutine ggfit_eval_dervs(gvar, n, x, y, np)
    integer,intent(IN):: n,np
    real(dl),intent(IN)::x(n)
    real(dl),intent(OUT)::y(n,0:np)
    type(ggfit_vars) gvar
    integer i, j
    real(dl) xtbys(n, gvar%n), ss(n, gvar%n), eebys(n, gvar%n), invs(gvar%n)
    do i=1, gvar%n
       do j=1,n
          call ggfit_get_tilted(gvar%a(i), x(j)-gvar%c(i), xtbys(j, i))
          eebys(j, i) = (1.d0 + gvar%a(i)*xtbys(j,i))/gvar%s(i)
          xtbys(j, i) = xtbys(j,i)/gvar%s(i)
          ss(j, i) = gvar%h(i) * dexp(-(xtbys(j,i))**2/2.d0)
       enddo
       invs(i) = 1.d0/gvar%s(i)
    enddo
    do j = 1, n
       y(j, 0) = sum(ss(j,1:gvar%n))
    enddo
    if(np.eq.0)return
    do j = 1 , n
       y(j, 1) = - sum(ss(j, 1:gvar%n)  * eebys(j, 1:gvar%n) * xtbys(j,1:gvar%n) )
    enddo
    if(np.eq.1)return
    do j=1,n
       y(j, 2) = sum(ss(j, 1:gvar%n) * eebys(j, 1:gvar%n)  * (invs(1:gvar%n) + (xtbys(j,1:gvar%n)**2 -2.d0) * eebys(j, 1:gvar%n)))
    enddo
    if(np.eq.2)return
    do j=1,n
       y(j, 3) = sum(ss(j, 1:gvar%n) * eebys(j,1:gvar%n) * &
            eebys(j,1:gvar%n)*xtbys(j, 1:gvar%n)*((4.d0-xtbys(j, 1:gvar%n)**2)*eebys(j, 1:gvar%n)-invs(1:gvar%n)) &
            + gvar%a(1:gvar%n)*(invs(1:gvar%n) + 2.d0*eebys(j, 1:gvar%n)*(xtbys(j,1:gvar%n)**2-2.d0)) &
            )
    enddo
    if(np.eq.3) return
  end subroutine ggfit_eval_dervs


end module Gaussianfit_utils


