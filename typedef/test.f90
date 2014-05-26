program test
  use coop_types
  implicit none

  type(coop_cosmology_null):: cosmo

  cosmo  = coop_cosmology_null("LCDM", 1)
  call cosmo%print
  print*, cosmo%name
end program test
