module coop_wrapper_typedef
  use coop_page_mod
  use coop_constants_mod
  use coop_svd_mod
  use coop_basicutils_mod
  use coop_sort_mod
  use coop_sortrev_mod
  use coop_string_mod
  use coop_arguments_mod
  use coop_function_mod
  use coop_particle_mod
  use coop_species_mod
  use coop_cosmology_mod
  implicit none

contains

  subroutine coop_wrapper_typedef_print()
    write(*,*) "This is COOP VERSION "//trim(coop_version)
    write(*,*) "Wrapper for typedef."
  end subroutine coop_wrapper_typedef_print


end module coop_wrapper_typedef
