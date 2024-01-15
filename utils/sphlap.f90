module coop_sphlap_mod
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
  use coop_sphere_mod
  implicit none

#include "constants.h"

  type coop_sphlap_obj
     COOP_INT::n
     COOP_REAL::dr, dr2, twodr, c2
     COOP_REAL, dimension(:), allocatable:: r, rdr, c1, yeq
   contains
     procedure::free => coop_sphlap_free
     procedure::init => coop_sphlap_init
     procedure::solve => coop_sphlap_solve
     procedure::solve_exp => coop_sphlap_solve_exp
  end type coop_sphlap_obj

contains

  subroutine coop_sphlap_obj_free()
    class(coop_sphlap_obj)::this
    COOP_DEALLOC(this%r)
    COOP_DEALLOC(this%rdr)
    COOP_DEALLOC(this%c1)
    COOP_DEALLOC(this%yeq)
  end subroutine coop_sphlap_obj_free


  subroutine coop_sphlap_obj_init(n, dr)
    class(coop_sphlap_obj)::this
    COOP_INT::i
    call this%free()
    allocate(this%r(0:n), this%rdr(0:n), this%c1(n), yeq(0:n))
    this%dr = dr
    this%dr2 = this%dr**2
    this%twodr = this%dr*2.d0
    this%c2 =  0.25d0/this%dr2
    this%r = (/ (i = 0, n) /) * dr
    this%rdr = this%r * dr
    this%c1 = 1.d0/this%dr2 + 1.d0/this%rdr(1:n)
  end subroutine coop_sphlap_obj_init

  !!solve \nabla^2 y = f(y) in a sphere, assuming spherical symmetry y = y(r) , y(0) and y(n) are known
  subroutine coop_sphlap_obj_solve(this, y, f)
    class(coop_sphlap_obj)::this
    COOP_REAL::y(0:this%n)
    external f
    COOP_REAL::f
    stop " coop_sphlap_obj_solve to be done."
  end subroutine coop_sphlap_obj_solve


  !!solve \nabla^2 (e^y) = f(y) in a sphere, assuming spherical symmetry y = y(r) and y>0; 
  subroutine coop_sphlap_obj_solve_exp(this, y, f)
    class(coop_sphlap_obj)::this
    COOP_REAL::y(0:this%n)
    external f
    COOP_REAL::f
    COOP_INT::i, ibreak, j
    !$omp parallel do
    do i=0, this%n
       this%yeq(i) = 
    enddo
    !$omp end parallel do
  end subroutine coop_sphlap_obj_solve_exp

  
end module coop_sphlap_mod
