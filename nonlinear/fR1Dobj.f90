module fR1d_mod
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"

  type fR1D_object
     COOP_INT::nr
     !!number of discretized comoving coordinates
     COOP_INT::ns
     !!number of mass shells
     COOP_REAL::rmax, dr, Omega_m, Omega_Lambda
     COOP_REAL:: a, H
     COOP_REAL,dimension(:),allocatable::q, v, m, r, rho, psum, pdiff
     !!q: comoving positions of mass shells (time dependent)
     !!v = dq/d tau: comoving velocities of mass shells (time dependent)
     !!m: masses of shells (fixed)
     !!r: discretized comoving coordinates (fixed)
     !!rho: discretized density
     !!psum: discretized Psi+Phi
     !!pdiff: discretized Psi-Phi
   contains
     procedure:: free => fR1d_object_free
     procedure:: init => fR1d_object_init
     procedure:: Hofa => fR1D_object_Hofa
  end type fR1D_object

contains

  function fR1D_object_Hofa(this, a) result(H)
    class(fR1d_object)::this
    COOP_REAL::H
    COOP_REAL, optional::a
    if(present(a))then
       H = sqrt(this%Omega_m/a + this%Omega_Lambda*a**2)
    else
       H = sqrt(this%Omega_m/this%a + this%Omega_Lambda*this%a**2)
    endif
  end function fR1D_object_Hofa

  subroutine fR1d_object_free(this)
    class(fR1d_object)::this
    COOP_DEALLOC(this%q)
    COOP_DEALLOC(this%m)
    COOP_DEALLOC(this%r)
    COOP_DEALLOC(this%rho)
    COOP_DEALLOC(this%psum)
    COOP_DEALLOC(this%pdiff)
  end subroutine fR1d_object_free

  subroutine fR1d_object_init(this, Omega_m, nr, rmax, ns, a_ini, delta_ini, r_halo)
    class(fR1d_object)::this
    COOP_REAL::omega_m, rmax, a_ini, delta_ini, r_halo
    COOP_INT::ns, nr
    COOP_REAL::dq, rhobar
    COOP_INT::i, j
    call this%free()
    this%rmax = rmax
    this%nr  = nr
    this%ns = ns
    this%dr = rmax/nr
    this%Omega_m = Omega_m
    this%Omega_Lambda = 1.d0 - Omega_m
    allocate(this%r(nr), this%rho(nr), this%psum(nr), this%pdiff(nr), this%q(ns), this%m(ns))
    this%r = ( (/ ( i, i = 1, nr ) /) - 0.5d0 ) * this%dr
    this%a = a_ini
    this%H = this%Hofa()
    dq = this%rmax/ns
    this%q = ( (/ ( i, i = 1, ns ) /) - 0.5d0 ) * dq
    this%v = 0.d0
    rhobar = coop_pi*4.d0*Omega_m/a_ini**3
    !$omp parallel do
    do i=1, ns
       this%m(i) = rhobar * ((this%q(i)+dq/2)**3 - (this%q(i)-dq/2)**3)
       if(this%q(i) .le. r_halo) this%m(i) = this%m(i)*(1.d0+delta_ini)
    enddo
    !$omp end parallel do
    
  end subroutine fR1d_object_init



end module fR1d_mod
