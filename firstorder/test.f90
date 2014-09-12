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
  call fod%init(h=COOP_REAL_OF(0.68d0))
  call fod%add_species(coop_baryon(COOP_REAL_OF(0.045d0)))
  call fod%add_species(coop_cdm(COOP_REAL_OF(0.255d0)))
  call fod%add_species(coop_radiation(fod%Omega_radiation()))
  call fod%add_species(coop_neutrinos_massless(fod%Omega_massless_neutrinos_per_species()*3))
 ! call fod%add_species(coop_neutrinos_massive(fod%Omega_nu_per_species_from_mnu_eV(1.d-6)*1,fod%Omega_massless_neutrinos_per_species()*1))
  call fod%add_species(coop_de_lambda(fod%Omega_k()))
  call fod%setup_background()
  fod%optre = 0.07
  call fod%set_xe()
  print*, "Omega_b = ", fod%Omega_b
  print*, "Omega_c = ", fod%Omega_c
  print*, "Omega_m = ", fod%Omega_m
  print*, "Omega_g = ", fod%Omega_g
  print*, "Omega_nu = ", fod%Omega_nu
  print*, "Omega_r = ", fod%Omega_r
  print*, "Omega_de = ", O0_DE(fod)%Omega


  call fod%compute_source(0)
  ik = 50
  call fp%open('solution.txt', 'w')
  do i=1, fod%source(0)%ntau
     write(fp%unit, "(20E16.7)") fod%source(0)%lna(i), fod%source(0)%s(i, ik, :)
  enddo
  call fp%close()

  


end program test
