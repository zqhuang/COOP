module coop_wrapper
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
contains

  subroutine coop_print_info()
    write(*,*) "This COOP Version "//trim(coop_version)
    write(*,*) "Author: Zhiqi Huang"
  end subroutine coop_print_info
end module coop_wrapper
