program fR1d
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  COOP_INT,parameter::nr = 30000, nt = 80000, output_span=400
  COOP_REAL,parameter::eps = 1.d-6, n=1.d0, Lambda = 2.1d0, minR = 1.d-8, varphi_min = -1.d-4, varphi_max = -1.d-12
  COOP_REAL::r(nr), dr, rmax,  dt,  dotvarphi(nr), lapvarphi(nr), rho(nr), varphi(nr), ricci(nr)
  COOP_REAL,dimension(:,:),allocatable::varphilog
  COOP_INT::ir, it, i, j
  type(coop_file)::fp
  allocate(varphilog(nr, 0:nt/output_span))
  rmax = 0.1d0
  dr = rmax/nr
  r = ( (/ (i, i = 1, nr) /) -0.5d0) * dr

!  rho = 50.d0*(1.d0+tanh(( rmax/4 - r)/(rmax/20.d0)))
  rho = 100.d0*exp(-(r/rmax*8.d0)**2)
  call get_ini()
  dotvarphi = 0.d0
  dt = dr/5.d0
  varphilog(:,0) =varphi

  do it = 1, nt
     !$omp parallel do
     do ir = 1, nr
        ricci(ir) = invF(varphi(ir))
     enddo
     !$omp end parallel do
     call get_lap(varphi, lapvarphi)
     dotvarphi = dotvarphi + (lapvarphi + (rho - ricci)/3.d0)*dt
     varphi(1:nr-1) = max(min(varphi(1:nr-1) + dotvarphi(1:nr-1) * dt, varphi_max), varphi_min)
     if(mod(it, output_span).eq.0)then
        varphilog(:,it/output_span) =varphi
!        dotvarphi = 0.d0
     endif
  enddo
  
  call fp%open("phi.txt")
  do i=1, nt/output_span
     write(fp%unit, "("//COOP_STR_OF(nr/output_span+1)//"E14.5)") i*dt*output_span, varphilog(output_span:nr:output_span, i)
  enddo
  call fp%close()


  call fp%open("psi.txt")
  do i=1, nr, 10
     write(fp%unit, "("//COOP_STR_OF(nt/output_span+2)//"E14.5)") r(i), varphilog(i, :)
  enddo
  call fp%close()

contains

  subroutine get_ini()
    call get_QS(rho, varphi)

    dotvarphi = 0.d0
    dt = dr/5.d0
    varphilog(:,0) =varphi

    do it = 1, nt
       !$omp parallel do
       do ir = 1, nr
          ricci(ir) = invF(varphi(ir))
       enddo
       !$omp end parallel do
       call get_lap(varphi, lapvarphi)
       dotvarphi = dotvarphi + (lapvarphi + (rho - ricci)/3.d0)*dt
       varphi(1:nr-1) = max(min(varphi(1:nr-1) + dotvarphi(1:nr-1) * dt, varphi_max), varphi_min)
       if(mod(it, output_span).eq.0)then
          varphilog(:,it/output_span) =varphi
          dotvarphi = 0.d0
       endif
    enddo
  end subroutine get_ini

  function F(R)
    COOP_REAL::R, F
    F = -eps*(R/(4.d0*Lambda) + (-varphi_min/eps)**(-1.d0/(n+1.d0)) )**(-(n+1.d0))
  end function F


  function invF(varphi)
    COOP_REAL::varphi, invF
    invF = 4.d0*Lambda* (max(-varphi/eps, -varphi_max/eps)**(-1.d0/(n+1.d0)) - (-varphi_min/eps)**(-1.d0/(n+1.d0)))
  end function invF

  subroutine get_lap(f, lapf)
    COOP_REAL::f(nr), lapf(nr), rf(0:nr+1)
    COOP_INT::i
    rf(1:nr) = r*f
    rf(0) = 0.d0
    rf(nr+1) = f(nr)*(r(nr)+dr)
    do i = 1, nr
       lapf(i) = (rf(i+1)+rf(i-1)-2.d0*rf(i))/(r(i)*dr**2)
    enddo
  end subroutine get_lap


  subroutine get_QS(rho, varphi)
    COOP_REAL::varphi(nr), rho(nr)
    COOP_INT::i
    do i=1, nr
       varphi(i) = F(rho(i))
    enddo
  end subroutine get_QS
  
end program fR1d
