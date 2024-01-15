module coop_wrapper_background
  use coop_wrapper_utils
  use coop_background_mod
  implicit none

contains

  subroutine coop_wrapper_background_print()
    write(*,*) "This is COOP VERSION "//trim(coop_version)
    write(*,*) "Wrapper for background"
  end subroutine coop_wrapper_background_print
  
end module coop_wrapper_background
