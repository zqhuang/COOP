program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::fod
  type(coop_pert_object)::pert
  type(coop_file)fp
  integer i, m, iq, ik,j, l, sterm
  COOP_REAL, dimension(coop_num_Cls, 2:2500)::Cls
  COOP_REAL norm
  COOP_REAL z, a, s, stau
  call fod%Set_Planck_bestfit()
  call coop_prtsystime(.true.)
  call fod%compute_source(m=0)
  call coop_prtsystime()
  call coop_prtsystime(.true.)
  call fod%source(0)%get_All_Cls(2, 2500, Cls)
  call coop_prtsystime()
  call fp%open('Cls.txt', 'w')
  norm = 2.726**2*1.d12
  do l=2, 2500
     write(fp%unit, "(I5, 20E16.7)") l, Cls(:, l)*(l*(l+1.d0)/coop_2pi*norm)
  enddo
  call fp%close()

  

end program test
