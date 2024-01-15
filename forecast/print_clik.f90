program test
  use coop_wrapper_firstorder
  use coop_forecast_mod
  implicit none
#include "constants.h"
  type(coop_clik_object)::pl
  COOP_STRING::filename
  logical::wp
  call coop_MPI_init()
  call coop_get_command_line_argument(key = 'data', arg = filename)
  call coop_get_command_line_argument(key = 'wp', arg = wp, default = .false.)
  call pl%init(trim(filename))
  call pl%print_names(with_param = wp)
  call coop_MPI_finalize()
end program test
