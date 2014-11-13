program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  COOP_INT, parameter::n=10005
  COOP_REAL x(n), y(n), a, b
  COOP_INT i, m
  type(coop_function)::f
  open(10, file = 'test_select.txt')
  i = 0
  do
     i = i + 1
     read(10, *, END=100, ERR=100) x(i), y(i)
     y(i) = exp(y(i))
     m = i
  enddo
100 close(10)
  print*, "found ", m, " numbers"
  call f%Init_NonUniform(x=x(1:m), f=y(1:m), ylog = .true., smooth = .true.)
  a = x(1)
  b = x(m)
  call coop_set_uniform(n, x, a, b)
  open(10, file ="f12.txt")
  do i=1, n
     write(10,"(4E16.7)") x(i), f%eval(x(i)), f%derivative(x(i)), f%derivative2(x(i))
  enddo
  close(10)
end program Test
