module coop_sphere_mod
  use coop_wrapper_typedef
  use coop_MPI_mod
  use coop_list_mod
  use coop_random_mod
  use coop_special_function_mod
  use coop_cholesky_mod  
  use coop_matrix_mod
  use coop_interpolation_mod
  use coop_integrate_mod
  use coop_ode_mod
  use coop_file_mod
  use coop_asy_mod
  use coop_fft_mod
  use coop_jl_mod
  use coop_gaussian_peak_stat_mod
  use coop_nd_prob_mod
  use coop_evalstr_mod
  implicit none

#include "constants.h"


  private

  public::coop_sphere_disc, coop_sphere_ang2vec, coop_sphere_vec2ang, coop_sphere_disc_constructor


  type coop_sphere_disc
     COOP_REAL, dimension(3)::nx, ny, nz
     COOP_REAL theta, phi
   contains
     procedure::init => coop_sphere_disc_initialize
     procedure::ang2flat => coop_spere_disc_ang2flat
  end type coop_sphere_disc
!!$
!!$  interface coop_sphere_disc
!!$     module procedure coop_sphere_disc_constructor
!!$  end interface coop_sphere_disc

contains

  subroutine coop_sphere_ang2vec(theta, phi, vec)
    COOP_REAL theta, phi, vec(3)
    vec(1:2) = sin(theta)*  (/ cos(phi), sin(phi) /)
    vec(3) = cos(theta)
  end subroutine coop_sphere_ang2vec

  subroutine coop_sphere_vec2ang(vec, theta, phi)
    COOP_REAL vec(3), theta, phi, norm
    if(vec(1) .eq. 0.d0 .and. vec(2).eq.0.d0)then
       if(vec(3).ge.0.d0)then
          theta = 0.d0
       else
          theta = coop_pi
       endif
       phi = 0.d0
    else
       norm = sqrt(sum(vec**2))
       theta = acos(vec(3)/norm)
       phi = atan2(vec(2), vec(1))
    endif
  end subroutine coop_sphere_vec2ang

  subroutine coop_sphere_disc_initialize(this, theta, phi, theta_x, phi_x)
    class(coop_sphere_disc)::this
    COOP_REAL theta, phi, norm
    COOP_REAL,optional:: theta_x, phi_x
    this%theta = theta
    this%phi = phi
    call coop_sphere_ang2vec(this%theta, this%phi, this%nz)
    if(present(theta_x) .and. present(phi_x))then
       call coop_sphere_ang2vec(theta_x, phi_x, this%nx)
       this%nx = this%nx - dot_product(this%nx, this%nz)*this%nz
       norm = sqrt(sum(this%nx**2))
       if(norm .gt. 1.d-10)then
          this%nx  = this%nx / norm
       else
          this%nx = (/ sin(this%phi), - cos(this%phi), COOP_REAL_OF(0.) /)
       endif
    else
       this%nx = (/ sin(this%phi), - cos(this%phi), COOP_REAL_OF(0.) /)
    endif
    call coop_vector_cross_product(this%nz, this%nx, this%ny)
  end subroutine coop_sphere_disc_initialize


  function coop_sphere_disc_constructor(theta, phi, theta_x, phi_x) result(this)
    type(coop_sphere_disc)::this
    COOP_REAL theta, phi, norm
    COOP_REAL,optional:: theta_x, phi_x
    this%theta = theta
    this%phi = phi
    call coop_sphere_ang2vec(this%theta, this%phi, this%nz)
    if(present(theta_x) .and. present(phi_x))then
       call coop_sphere_ang2vec(theta_x, phi_x, this%nx)
       this%nx = this%nx - dot_product(this%nx, this%nz)*this%nz
       norm = sqrt(sum(this%nx**2))
       if(norm .gt. 1.d-10)then
          this%nx  = this%nx / norm
       else
          this%nx = (/ sin(this%phi), - cos(this%phi), COOP_REAL_OF(0.) /)
       endif
    else
       this%nx = (/ sin(this%phi), - cos(this%phi), COOP_REAL_OF(0.) /)
    endif
    call coop_vector_cross_product(this%nz, this%nx, this%ny)
  end function coop_sphere_disc_constructor

  !!input coor = ( theta, phi), output coor = (x, y)
  subroutine coop_spere_disc_ang2flat(this, coor) 
    class(coop_sphere_disc)::this
    COOP_REAL s, coor(2), vec(3)
    call coop_sphere_ang2vec(coor(1), coor(2), vec)
    coor(1) = dot_product(this%nx, vec)
    coor(2) = dot_product(this%ny, vec)
    s = sum(coor**2)
    if(s.lt. coop_tiny)return
    s = sqrt(2.d0*(1.d0 - dot_product(this%nz , vec))/s)
    coor(1) = s*coor(1)
    coor(2) = s*coor(2)
  end subroutine coop_spere_disc_ang2flat
  

end module coop_sphere_mod
