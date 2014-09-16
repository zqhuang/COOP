program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::fod
  type(coop_pert_object)::pert
  type(coop_file)fp
  integer i, m, iq, ik,j, l
  COOP_REAL, dimension(coop_num_Cls)::Cls
  COOP_REAL norm
  COOP_REAL z
  call fod%Set_Planck_bestfit()
  call coop_prtsystime(.true.)
  call fod%compute_source(m=0)
  call coop_prtsystime()

  print*, fod%source(0)%k(1), fod%source(0)%k(10)

  norm = 2.726**2*1.d12
  call fp%open('solution.txt', 'w')
  do l=2, 5
     call write_l(l)
  enddo
  do l = 10, 100, 10
     call write_l(l)
  enddo
  do l= 120, 500, 20
     call write_l(l)
  enddo
  do l = 525, 1250, 25
     call write_l(l)
  enddo
  call fp%close()

  
contains

  subroutine write_l(l)
    COOP_INT l
    call fod%source(0)%get_Cls(l, Cls)
    write(*, *) l , l*(l+1.)/coop_2pi*Cls(coop_index_ClTT)*norm
    write(fp%unit, *) l , l*(l+1.)/coop_2pi*Cls(coop_index_ClTT)*norm, l*(l+1.)/coop_2pi*Cls(coop_index_ClTE)*norm, l*(l+1.)/coop_2pi*Cls(coop_index_ClEE)*norm
  end subroutine write_l
   
end program test
