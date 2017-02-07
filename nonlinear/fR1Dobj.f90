module fR1d_mod
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"

  COOP_INT,parameter::coop_fr1d_rho_maxp = 14

  type coop_fr1d_obj
     COOP_INT::nr
     !!number of discretized comoving coordinates
     COOP_INT::ns
     !!number of mass shells
     COOP_REAL::rmax, dr, Omega_m, Omega_Lambda
     COOP_REAL:: a, H
     COOP_REAL,dimension(:),allocatable::q, v, m, intm, rhobar, r, rho, phi, phip
     !!q: comoving positions of mass shells (time dependent)
     !!v = dq/d tau: comoving velocities of mass shells (time dependent)
     !!m: masses of shells (fixed)
     !!r: discretized comoving coordinates (fixed)
     !!rho: discretized density
     !!intm: integrated mass
     !!phi: discretized scalar field
     !!phip: d phi/ d tau
   contains
     procedure:: free => coop_fr1d_obj_free
     procedure:: init => coop_fr1d_obj_init
     procedure:: Hofa => coop_fr1d_obj_Hofa
     procedure:: get_rho => coop_fr1d_obj_get_rho
     procedure:: sort_shells => coop_fr1d_obj_sort_shells
     procedure:: integrate_mass => coop_fr1d_obj_integrate_mass
  end type coop_fr1d_obj

contains


  subroutine coop_fr1d_obj_get_rho(this)
    class(coop_fr1d_obj)::this
    COOP_INT::i, j1, j2
    COOP_REAL::fr, frp, lambda, low, up
    call this%sort_shells()
    call this%integrate_mass()
    j1 = 0
    
    do i=1, this%nr
       low = this%r(i) - this%dr/4.d0
       up = this%r(i) + this%dr/4.d0
       do while(this%q(j1+1) .lt. low )
          j1 = j1+ 1
       enddo
       j2 = j1 +1 
       do while(this%q(j2).lt. up)
          j2  = j2+1
       enddo
       lambda = (this%r(i) - this%q(j1))/(this%q(j2)-this%q(j1))       
       fr = lambda*this%rhobar(j2) + (1.d0-lambda)*this%rhobar(j1)
       frp = (this%r(i)/3.d0)*(this%rhobar(j2)-this%rhobar(j1))/(this%q(j2)-this%q(j1))
       lambda = (tanh(log((this%q(j2)-this%q(j1))/this%r(i))*2.5d0 + 8.d0)+1.d0)/2.d0
       this%rho(i) = (fr+frp) * lambda &
               +  (this%intm(j2)-this%intm(j1))/(this%q(j2)-this%q(j1))/(coop_pi*(this%q(j1)+this%q(j2))**2)*(1.d0-lambda)

    enddo
  end subroutine coop_fr1d_obj_get_rho


  subroutine coop_fr1d_obj_sort_shells(this)
    class(coop_fr1d_obj)::this
    COOP_INT::i, j1, j2, j3
    COOP_REAL::vtmp, qtmp
    do i = 2, this%ns
       if( this%q(i) .ge. this%q(i-1) ) cycle
       j1 = 1
       j2 = i-1
       do while(j2-j1 .gt. 1)
          j3 = (j1+j2)/2
          if(this%q(i) .ge. this%q(j3))then
             j1 = j3
          else
             j2 = j3
          endif
       enddo
       vtmp = this%v(i)
       qtmp = this%q(i)
       this%v(j2+1:i) = this%v(j2:i-1)
       this%q(j2+1:i) = this%q(j2:i-1)
       this%v(j2) = vtmp
       this%q(j2) = qtmp
    enddo
  end subroutine coop_fr1d_obj_sort_shells


  subroutine coop_fr1d_obj_integrate_mass(this)
    class(coop_fr1d_obj)::this
    COOP_INT::i
    COOP_REAL::lambda1, lambda2, u, v
    lambda1 = (this%q(1)*2.d0/(this%q(1)+this%q(2)))**3
    this%intm(0) = 0.d0
    this%intm(1) = this%m(1)*lambda1
    do i=2, this%ns-1
       u = (this%q(i+1) - this%q(i-1))/(this%q(i)+this%q(i-1))
       v = (this%q(i) - this%q(i-1))/(this%q(i)+this%q(i-1))
       if(u .gt. 0.d0)then
          lambda2 = (v/u)* ((v**2 + 3.d0*(v+1.d0))/(u**2 + 3.d0*(u+1.d0)))
       else
          lambda2 = 0.5d0
       endif
       this%intm(i) = this%intm(i-1) + this%m(i-1)*(1.d0-lambda1) + this%m(i)*lambda2
       lambda1 = lambda2
    enddo
    lambda2 = 0.5d0 * (1.d0 - (this%ns - 0.5d0)/(2.d0*this%ns*(this%ns-2.d0/3.d0)  + 1.d0))  !!this is the initial weight
    this%intm(this%ns) = this%intm(this%ns-1) + this%m(this%ns-1)*(1.d0-lambda1) + this%m(this%ns)*lambda2
    this%intm(this%ns+1) = this%intm(this%ns) + this%m(this%ns)*(1.d0-lambda2)
    this%rhobar(1:this%ns+1) = this%intm(1:this%ns+1)/((coop_pi*4.d0/3.d0)*this%q(1:this%ns+1)**3)
    this%rhobar(0) = this%rhobar(1)
  end subroutine coop_fr1d_obj_integrate_mass


  function coop_fr1d_obj_Hofa(this, a) result(H)
    class(coop_fr1d_obj)::this
    COOP_REAL::H
    COOP_REAL, optional::a
    if(present(a))then
       H = sqrt(this%Omega_m/a + this%Omega_Lambda*a**2)
    else
       H = sqrt(this%Omega_m/this%a + this%Omega_Lambda*this%a**2)
    endif
  end function coop_fr1d_obj_Hofa

  subroutine coop_fr1d_obj_free(this)
    class(coop_fr1d_obj)::this
    COOP_DEALLOC(this%q)
    COOP_DEALLOC(this%m)
    COOP_DEALLOC(this%r)
    COOP_DEALLOC(this%rho)
    COOP_DEALLOC(this%intm)
    COOP_DEALLOC(this%rhobar)
    COOP_DEALLOC(this%phi)
    COOP_DEALLOC(this%phip)
  end subroutine coop_fr1d_obj_free

  subroutine coop_fr1d_obj_init(this, Omega_m, nr, rmax, ns, a_ini, delta_ini, r_halo, bw_halo)
    class(coop_fr1d_obj)::this
    COOP_REAL::omega_m, rmax, a_ini, delta_ini, r_halo, bw_halo
    COOP_INT::ns, nr
    COOP_REAL::dq, rhobar, mfac
    COOP_INT::i, j
    call this%free()
    this%rmax = rmax
    this%nr  = nr
    this%ns = ns
    this%dr = rmax/nr
    this%Omega_m = Omega_m
    this%Omega_Lambda = 1.d0 - Omega_m

    allocate(this%r(nr), this%intm(0:ns+1), this%rhobar(0:ns+1), this%phi(nr), this%phip(nr), this%q(0:ns+1), this%v(ns), this%m(ns), this%rho(nr) )


    this%r = ( (/ ( i, i = 1, nr ) /) - 0.5d0 ) * this%dr

    this%a = a_ini
    this%H = this%Hofa()
    dq = this%rmax/ns
    this%q(1:ns) = ( (/ ( i, i = 1, ns ) /) - 0.5d0 ) * dq
    this%q(0) = 0.d0
    this%q(ns+1) = rmax
    this%v = 0.d0
    rhobar = coop_pi*4.d0*Omega_m


    !set mass
    !$omp parallel do
    do i=1, ns
       this%m(i) = rhobar * ((this%q(i)+dq/2)**3 - (this%q(i)-dq/2)**3)
    enddo
    !$omp end parallel do    

    mfac = (1.d0+delta_ini)
    this%m = this%m * (mfac + (1.d0-mfac)*(0.5d0+tanh((this%q(1:ns)-r_halo)/bw_halo)/2.d0))
    call this%get_rho()
  end subroutine coop_fr1d_obj_init



end module fR1d_mod
