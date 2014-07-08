program test
  use coop_wrapper_utils
  use coop_background_mod
  use coop_wrapper
  implicit none
#include "constants.h"
  integer,parameter::n = 30
  integer,parameter::m = 50
  integer i, j, k
  COOP_REAL::w(m), omch2(n), h(n, m)
  COOP_REAL:: plt_h(n), plt_w(n, n)
  type(coop_asy)::fig
  call coop_set_uniform(n, w, -1.25d0, -0.75d0)
  call coop_set_uniform(n, omch2, 0.11d0, 0.125d0)
  do j=1, n
     do i=1,n
        COOP_COSMO_PARAMS = coop_arguments(r= (/ 0.02214d0, omch2(i), 1.041d-2,  0.088d0, 0.d0, 0.06d0, w(j) /), i = (/ COOP_DE_W0, 7, 1 /) )
        call coop_setup_global_cosmology()
        h(i,j) =  COOP_COSMO%h()
     enddo
  enddo
  h = h*100.d0
  call coop_set_uniform(n, plt_h, 66.5d0, 72.5d0)
  do j=1, m
     do i=1, n
        if(plt_h(j) .le. h(i,m) .or. plt_h(j) .ge. h(i, 1))then
           print*, plt_h(j), h(i, 1), h(i, n)
           stop "h over flow"
        else
           k = 1
           do while(h(i, k) .gt. plt_h(j))
              k = k + 1
           enddo
           plt_w(i, j) = (w(k)*(plt_h(j) - h(i, k)) +  w(k-1)*(h(i, k-1)-plt_h(j)))/(h(i, k-1)- h(i, k))
        endif
     enddo
  enddo
  call fig%open("hwomch2.txt")
  call fig%init(xlabel = "$\Omega_c h^2$", ylabel = "$H_0$")
  call coop_asy_density(fig, plt_w, xmin = omch2(1), xmax = omch2(n), ymin = plt_h(1), ymax = plt_h(n), label = "$w$" )
  call coop_asy_label(fig, "$\theta_{\rm CMB} = 0.01041$ (Planck constraint)",  x = omch2(n/2), y = plt_h(n/2), color="black") 
  call fig%close()
end program test
