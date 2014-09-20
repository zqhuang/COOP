program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT :: l, i
  COOP_INT, parameter::n=2000
  COOP_REAL x(n), jl(n), amp(n), phase(n)
  type(coop_asy)::fig
  call coop_prtsystime(.true.)
  !$omp parallel do
  do l = 1, 5000
     call coop_jl_setup_amp_phase(l)
  enddo
  !$omp end parallel do
  call coop_prtsystime()
  write(*,*) "enter l = "
  read(*,*) l
  call coop_set_uniform(n, x, 1.d-5, max(l*2.5d0, 100.d0))
  jl = coop_jl(l, x)
  do i=1, n
     call coop_jl_get_amp_phase(l, x(i), amp(i), phase(i))
  enddo
  print*, maxval(abs(jl - amp*cos(phase)))/maxval(jl)
  call fig%open("jl.txt")
  call fig%init(xlabel = "$x$", ylabel = "$j_l(x)$")
  call coop_asy_curve(fig, x, jl, color="red")
  call coop_asy_curve(fig, x, amp*cos(phase), color="blue", linetype = "dotted")
  call fig%close()




!!$  subroutine get_max(l, x, jl)
!!$    COOP_INT l
!!$    COOP_REAL x, jl, jltry, dx
!!$    x = l + 0.5d0
!!$    jl = x*coop_jl(l, x)
!!$    dx = 0.1d0
!!$    do while(dx .gt. 1.d-8)
!!$       jltry = (x+dx)*coop_jl(l, x+dx)
!!$       do while(jltry .gt. jl)
!!$          jl = jltry
!!$          x = x + dx
!!$          jltry = (x+dx)*coop_jl(l, x+dx)
!!$       enddo
!!$       jltry = (x-dx)*coop_jl(l, x-dx)
!!$       do while(jltry .gt. jl)
!!$          jl = jltry
!!$          x = x - dx
!!$          jltry = (x-dx)*coop_jl(l, x-dx)
!!$       enddo
!!$       dx = dx/2.d0
!!$    enddo
!!$  end subroutine get_max

end program Test
