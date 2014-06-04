module coop_wrapper
  use coop_wrapper_firstorder
  implicit none

contains

  subroutine coop_print_info()
    write(*,*) "This COOP Version "//trim(COOP_VERSION)
    write(*,*) "Author: Zhiqi Huang"
  end subroutine coop_print_info
end module coop_wrapper
