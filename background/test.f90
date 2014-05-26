program test
  use coop_background
  implicit none
  type(coop_cosmology_background) bg
  type(coop_species_basic) baryon, cdm
  
  bg = coop_cosmology_background()
  call baryon%init(name="Baryon", id=1, Omega=0.046d0, w = 0.d0, cs2=0.d0)
  call cdm%init(name="CDM", id=2, Omega=0.26d0, w = 0.d0, cs2=0.d0)
  call bg%add_species(baryon)
  call bg%add_species(cdm)
  call bg%print()
end program test
