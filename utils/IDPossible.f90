program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT::month
  COOP_INT::date_start, date_end
  COOP_INT::id(18), i
  COOP_STRING::f10, f4, smonth
10 write(*,*) "Enter the first 10 digits:"
  read(*,*) f10
  f10 = trim(f10)
  if(len_trim(f10) .ne. 10)goto 10
  if(verify(f10(1:10), "0123456789").ne. 0) goto 10
  do i = 1, 10
     read(f10(i:i), *) id(i)
  enddo

20 write(*,*) "Enter the last 4 digits:"
  read(*,*) f4
  f4 = trim(f4)
  if(len_trim(f4) .ne. 4)goto 20
  if(verify(f4(1:3), "0123456789").ne. 0) goto 20
  if(verify(f4(4:4), "0123456789Xx").ne. 0) goto 20
  do i = 1, 3
     read(f4(i:i), *) id(i+14)
  enddo
  if(f4(4:4) .eq. "X" .or. f4(4:4) .eq. "x")then
     id(18) = 10
  else
     read(f4(4:4), *) id(18)
  endif
  
30 write(*,*) "Enter Month (1-12):"
  read(*,*) smonth
  if(verify(smonth(1:len_trim(smonth)), "0123456789") .ne. 0) goto 30
  read(smonth, *) month
  if(month .le. 0 .or. month .gt. 12) goto 30
  
  date_start = 1
  select case(month)  
  case(1,3,5,7,8,10,12)
     date_end = 31
  case(2)
     if( (mod(id(9)*2+id(10), 4) .eq. 0 .and. id(9)*id(10) .ne. 0) .or. (mod(id(7)*2+id(8), 4) .eq. 0) .and. id(9)+id(10) .eq. 0 )then
        date_end = 29
     else
        date_end = 28
     endif
  case default
     date_end = 30
  end select
  id(11) = month/10
  id(12) = mod(month, 10)
  do i=date_start, date_end
     id(13) = i/10
     id(14) = mod(i, 10)
     if(last_digit(id(1:17)) .eq. id(18))then
        write(*,"(18I1)") id
     endif
  enddo

contains

  function last_digit(f17) result(last)
    COOP_INT::f17(17)
    COOP_INT::last, s
    COOP_INT,dimension(17)::check = (/ 7, 9, 10, 5, 8, 4, 2, 1, 6, 3, 7, 9, 10, 5, 8, 4, 2 /)
    COOP_INT,dimension(11)::map = (/ 1, 0, 10, 9, 8, 7, 6, 5, 4, 3, 2 /)
    s = mod(sum(f17*check), 11)
    last=map(s+1)
  end function last_digit

end program Test  
