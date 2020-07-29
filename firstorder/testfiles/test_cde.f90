program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::cosmology
  type(coop_real_table)::paramtable
  logical success
  COOP_INT,parameter::nz=200, nk = 3
  COOP_REAL,parameter::zmin = 0.d0, zmax = 10.d0
  COOP_REAL:: z(nz), Dbya(nk, nz), Phi(nk, nz), Psi(nk, nz), a(nz), tau
  COOP_REAL::k(nk), kMpc(nk)
  COOP_INT::iz, ik
  type(coop_asy)::figure
  call paramtable%insert("ombh2", 0.022d0)
  call paramtable%insert("omch2", 0.12d0)
  call paramtable%insert("h", 0.67d0)
  call paramtable%insert("tau", 0.07d0)
  call paramtable%insert("As", 2.15d-9)
  call paramtable%insert("ns", 0.969d0)
  call paramtable%insert("de_Q", 0.4082d0)
  call paramtable%insert("de_np_index", 1.d0)
  call paramtable%insert("de_fR_epsilon", 0.05d0)
  call cosmology%set_up(paramtable, success)
  if(success)then
     print*, "cosmology initialized"
  else
     stop "initialization failed"
  endif

  call cosmology%compute_source(0, success)
  if(success) then
     write(*,*) "source computed"
  else
     stop "firstorder solver failed"
  endif

  call coop_set_uniform(nz, z, zmin, zmax)
  a = 1.d0/(1.d0+z)
  kMpc = (/ 0.001d0, 0.01d0, 0.1d0 /) !0.03d0, 0.1d0, 0.2d0 /)
  k = kMpc/cosmology%H0Mpc()
  do iz = 1, nz
     tau = cosmology%tauofa(a(iz))
     call cosmology%source(0)%get_delta_sync_trans(tau, nk, k, Dbya(:, iz))
     Dbya(:, iz) = Dbya(:,iz)/a(iz)
     call cosmology%source(0)%get_Phi_Trans( tau, nk, k, Phi(:, iz))
     call cosmology%source(0)%get_Psi_Trans(tau, nk, k, Psi(:, iz))
  enddo
  do iz=2, nz
     Dbya(:, iz) = Dbya(:, iz)/Dbya(:, 1)
     Phi(:, iz) = Phi(:, iz)/Phi(:, 1)
     Psi(:, iz) = Psi(:, iz)/Psi(:, 1)
  enddo
  Phi(:,1)= 1.d0
  Psi(:,1) = 1.d0
  Dbya(:,1) = 1.d0
  call figure%open("growth.txt")
  call figure%init(xlabel = "$z$", ylabel="$(1+z)D(z)$")
  do ik = 1, nk
     call figure%plot( z, Dbya(ik,:), color=figure%color(ik), linetype =figure%linetype(ik), linewidth = figure%linewidth(ik), legend="$k="//COOP_STR_OF(kMpc(ik))//" \mathrm{Mpc}^{-1}$")
  enddo
  call figure%legend(0.3, 0.5)
  call figure%close()

  call figure%open("NewtonianPotential.txt")
  call figure%init(xlabel = "$z$", ylabel="$(1+z)D(z)$")
  do ik = 1, nk
     call figure%plot( z, Phi(ik,:), color=figure%color(ik), linetype =figure%linetype(ik), linewidth = figure%linewidth(ik), legend="$k="//COOP_STR_OF(kMpc(ik))//" \mathrm{Mpc}^{-1}$")
  enddo
  call figure%legend(0.3, 0.5)
  call figure%close()

  call figure%open("SpatialCurvature.txt")
  call figure%init(xlabel = "$z$", ylabel="$(1+z)D(z)$")
  do ik = 1, nk
     call figure%plot( z, Psi(ik,:), color=figure%color(ik), linetype =figure%linetype(ik), linewidth = figure%linewidth(ik), legend="$k="//COOP_STR_OF(kMpc(ik))//" \mathrm{Mpc}^{-1}$")
  enddo
  call figure%legend(0.3, 0.5)
  call figure%close()


end program test
