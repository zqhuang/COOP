module coop_contours_mod
  use coop_wrapper_typedef
  use coop_interpolation_mod
  use coop_file_mod
  use coop_list_mod
  implicit none
  private
#include "constants.h"

  public::coop_contours, coop_image

  type coop_image
     COOP_INT::nx, ny
     COOP_SINGLE,dimension(:,:),allocatable::image
     
  end type coop_image

  type coop_contours
  end type coop_contours

contains

  
end module coop_contours_mod
