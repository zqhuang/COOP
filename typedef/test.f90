program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  print*, COOP_DATAPATH_I("%DATASETDIR%cmb/cmbdata_dataset%u.dat", 2)
end program Test
