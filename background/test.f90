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
  print*, "at z=2, H0 * conformal time = ", bg%conformal_time(1.d0/(1.d0+2.d0))
  print*, "        H(z)/H0 = ",  bg%Hratio(1.d0/(1.d0+2.d0))
  print*, "CDM is the ", bg%index_of("CDM"), "th species"
end program bgtest
