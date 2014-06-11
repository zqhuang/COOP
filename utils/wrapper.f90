module coop_wrapper_utils
  use coop_wrapper_typedef
  use coop_MPI_mod
  use coop_sort_mod
  use coop_list_mod
  use coop_random_mod
  use coop_special_function_mod
  use coop_matrix_mod
  use coop_interpolation_mod
  use coop_integrate_mod
  use coop_ode_mod
  use coop_file_mod
  use coop_asy_mod
  use coop_fft_mod
  implicit none

contains

  subroutine coop_wrapper_utils_print()
    write(*,*) "This is COOP VERSION "//trim(coop_version)
    write(*,*) "Wrapper for utils"
  end subroutine coop_wrapper_utils_print


end module coop_wrapper_utils
