program bgtest
  use coop_wrapper_background
  implicit none
  type(coop_cosmology_background)::bg
  call bg%init(h=0.68d0)
  call bg%add_species(coop_baryon(0.046d0))
  call bg%add_species(coop_radiation(bg%Omega_radiation()))
  call bg%add_species(coop_neutrinos_massless(bg%Omega_massless_neutrinos_per_species()*(bg%Nnu())))
  call bg%add_species(coop_cdm(0.25d0))
  call bg%add_species(coop_de_lambda(bg%Omega_k()))
  call bg%setup_background()
  print*, "at z=1089, aH = ", bg%Hratio(1.d0/1090.d0)/1090.d0 * bg%H0Mpc()
end program bgtest
