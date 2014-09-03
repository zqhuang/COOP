program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  COOP_REAL lnam, lnrho, m
  m = 2.d0
  lnam = log(m)
  call coop_fermion_get_lnrho(lnam, lnrho)
  print*, exp(lnrho)
  print*, sum( sqrt( coop_fermion_int_q5 ** 2 + m ** 2 )* coop_fermion_int_q5**2 /(exp(coop_fermion_int_q5)+1)* coop_fermion_int_kernel5 ) /(7.d0/120.d0*coop_pi**4)
end program Test
