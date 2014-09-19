program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::fod
  type(coop_pert_object)::pert
  type(coop_file)fp
  integer i, m, iq, ik,j, l, sterm
  COOP_REAL, dimension(coop_num_Cls, 2:2500)::Cls
  COOP_REAL norm
  COOP_REAL z, a, s, stau

  !!set cosmology
  call fod%Set_Planck_bestfit()
  !if you want extended models
  !  call fod%set_standard_cosmology(Omega_b=0.047d0, Omega_c=0.26d0, h = 0.68d0, tau_re = 0.08d0, nu_mass_eV = 0.06d0, As = 2.15d-9, ns = 0.962d0, nrun = -0.01d0, r = 0.11d0, nt = -0.01d0, YHe = 0.25d0, Nnu = 3.d0)

  !!compute the scalar source
  call fod%compute_source(m=0)

  !!compute the Cl's
  print*, fod%Clzetazeta(100, 1.5d0), fod%Clzetazeta(100, 1.5d0, 1.500001d0)
  call fod%source(0)%get_All_Cls(2, 2500, Cls)

  call fp%open('Cls.txt', 'w')
  norm = 2.726**2*1.d12
  do l=2, 2500
     write(fp%unit, "(I5, 20E16.7)") l, Cls(:, l)*(l*(l+1.d0)/coop_2pi*norm)
  enddo
  call fp%close()

  

end program test
