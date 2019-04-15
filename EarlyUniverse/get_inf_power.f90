program getp
  use coop_wrapper_utils
  use coop_lattice_fields_mod
  implicit none
#include "constants.h"
  type(coop_lattice_initial_power)::p
  COOP_INT::ik, i, j
  COOP_REAL::P_R, lnkMpc, eps, lna, deltaphi
  COOP_REAL,dimension(:),allocatable::omega
  call coop_infbg%setup(nflds = 2, epsilon_end = 1.d0, f_ini = (/ 19.d0*coop_lattice_Mp, 0.d0 /) )
  write(*,"(A)") "================================================================================="
  write(*, "(A, E16.7)") "Inflation lasted for :", -coop_infbg%lna(1), " efolds"
  write(*,"(A)") "================================================================================="
  write(*, "(A, E16.7)") "H/M_p at the end of inflation:", exp(coop_infbg%lnH(coop_infbg%nsteps))
  write(*,"(A)") "---------------------------------------------------------------------------------"    
  write(*, "(A)") "Field values (/M_p) at the end of inflation:"
  write(*, "("//COOP_STR_OF(coop_infbg%nflds)//"E16.7)") coop_infbg%f(coop_infbg%nsteps, :)/coop_lattice_Mp
  write(*,"(A)") "---------------------------------------------------------------------------------"      
  write(*, "(A)") "Time derivatives of field values (/M_p^2) at the end of inflation:"
  write(*, "("//COOP_STR_OF(coop_infbg%nflds)//"E16.7)") coop_infbg%fd(coop_infbg%nsteps, :)/coop_lattice_Mpsq
  write(*,"(A)") "---------------------------------------------------------------------------------"
  lnkMpc = coop_infbg%lnkMpc()
  write(*,*) "Mpc => ln a = ", coop_infbg%lna_of_lnk(lnkMpc) 
  write(*,"(A)") "---------------------------------------------------------------------------------"
  !!---------------set up ----------------
  p%nk = 10
  p%nflds = coop_infbg%nflds
  allocate(p%lnk(p%nk))
  allocate(p%fdbyf(p%nflds, p%nk))
  allocate(p%f_cov(p%nflds, p%nflds, p%nk))
  write(*,*) "now calculate perturbatioins"
  !!------------------------------------------
  call p%init(lnkmin = lnkMpc-10.d0, lnkmax = lnkMpc, is_diagonal = .true.)
  write(*,*) "perturbations done", maxval(p%f_cov), minval(p%f_cov)
  !!------------------------------------------  
  allocate(omega(p%nflds))
  omega = coop_infbg%fd(coop_infbg%nsteps,:)/sqrt(sum(coop_infbg%fd(coop_infbg%nsteps,:)**2))
  do ik=1, p%nk
     P_R = dot_product(omega, matmul(p%f_cov(:,:,ik), omega))
     lna = coop_infbg%lna_of_lnk(p%lnk(ik))
     write(*,*) exp(p%lnk(ik)-lnkMpc),  P_R/2./coop_lattice_Mpsq/coop_infbg%eps(coop_infbg%nsteps), coop_infbg%Hubble(lna)**2/(8.*coop_pi2*coop_lattice_Mpsq*coop_infbg%epsilon(lna)) !!check slow roll
  enddo
end program Getp
