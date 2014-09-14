program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::fod
  type(coop_pert_object)::pert
  type(coop_file)fp
  integer i, m, iq, ik
  COOP_REAL ktauc
  COOP_REAL z
  call fod%set_standard_cosmology(h = 0.7d0, Omega_b = 0.045d0, Omega_c = 0.255d0, tau_re = 0.07d0, nu_mass_eV = 0.06d0)
  call coop_prtsystime(.true.)
  call fod%compute_source(m=0)
  call coop_prtsystime()
  ik = 30
  call fp%open('solution.txt', 'w')
  do i=1, fod%source(0)%ntau
     write(fp%unit, "(20E16.7)") fod%source(0)%lna(i), fod%source(0)%s(i, ik, :)
  enddo
  call fp%close()

  


end program test
