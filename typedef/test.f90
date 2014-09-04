program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  COOP_REAL lnam, lnrho, m
  m = 0.d0
  lnam = log(m)
  call coop_fermion_get_lnrho(lnam, lnrho)
  print*, exp(lnrho)*(7.d0/120.d0*coop_pi**4)
  print*, sum(sqrt( coop_fermion_int_q5 ** 2 + m ** 2 )* coop_fermion_int_q5**3/(exp( coop_fermion_int_q5)+1.d0) * coop_fermion_int_kernel5 )
end program Test
