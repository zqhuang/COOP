module coop_string
  use coop_constants
  use coop_basicutils
  implicit none
#include "constants.h"

  private
  public::coop_num2str,  coop_ndigits, coop_str2int, coop_str2real, coop_str2logical, coop_substr, coop_str_replace

  Interface coop_num2str
     procedure coop_int2str, coop_real2str, coop_logical2str
  end Interface coop_num2str


contains

  function coop_Ndigits(i, ndigits, base) result(str_ndigits)
    COOP_INT, optional::base
    COOP_INT i, ndigits, j 
    COOP_STRING Str_Ndigits, strzero
    if(present(base))then
       str_ndigits = coop_int2str(i, base)
    else
       str_ndigits = coop_int2str(i)
    endif
    strzero = ""
    do j = len_trim(str_ndigits)+1, ndigits
       strzero = trim(strzero)//"0"
    enddo
    str_Ndigits = trim(strzero)//trim(Str_Ndigits)
  end function coop_Ndigits

  function coop_int2str(i, base)
    COOP_INT,INTENT(IN)::i
    COOP_INT, optional::base
    COOP_STRING coop_int2str
    COOP_INT j, k
    character,dimension(36),parameter::dig = (/ "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", &
         "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", &
         "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", &
         "U", "V", "W", "X", "Y", "Z" /)
    coop_int2str = ""
    if(present(base))then
       if(base .gt. 36 .or. base .lt. 2) stop "coop_int2str: base overflow"
       j = i
       do while(j.gt. 0)
          k = j/base
          coop_int2str = dig(j-k*base+1)//trim(coop_int2str)
          j = k
       enddo
    else
       WRITE (coop_int2str,*) i
       coop_int2str= Trim(ADJUSTL(coop_int2str))
    endif
  End function coop_int2str

  Function Coop_real2str(x,fmt) !!This is a smater algorithm that convert number to string
    COOP_STRING Coop_real2str, str
    COOP_REAL x, ax
    COOP_INT ind,  n
    COOP_UNKNOWN_STRING ,optional::fmt
    if(present(fmt))then
       if(fmt .eq. "none")then
          coop_real2str = ""
          return
       endif
       if(trim(fmt).eq. "d")then
          write(str,"(G22.15)") x
          if(scan(str, "e").eq.0 .and. scan(str, "E").eq.0)then
             str=trim(adjustl(str))
             n=scan(str,".")
             if(n>0)then
                n=max(n-1,verify(trim(str),"0",.true.))
             else
                n=Len_Trim(str)
             endif
             Coop_real2str =  Str(1:n)//"d0"
          else
             Coop_real2str = trim(adjustl(str))
             n = scan(Coop_real2str, "e")
             if(n.eq.0) n = scan(Coop_real2str,"E")
             if(n.le.1) stop "error in Coop_real2str: does not look like a valid number"
             ind = n-1
             if(scan(Coop_real2str(1:ind),".") .gt. 0)then
                do while(Coop_real2str(ind:ind).eq."0")
                   ind = ind -1
                enddo
             endif
             if(Coop_real2str(n+1:n+1) == "+")then
                Coop_real2str = Coop_real2str(1:ind)//"d"//trim(Coop_real2str(n+2:))
             else
                Coop_real2str = Coop_real2str(1:ind)//"d"//trim(Coop_real2str(n+1:))
             endif
          endif
          return
       endif
       if(trim(fmt).ne."f" .and. trim(fmt).ne."")then
          write(str,"("//trim(fmt)//")") x
          Coop_real2str = trim(adjustl(str))
          return
       endif
    endif
    ax = abs(x)
    if(ax.eq.0.)then
       Coop_real2str="0"
       return
    endif
    ind=0
    do while(ax.lt.1.d0)
       ind = ind -1
       ax = ax*10.d0
    enddo
    do while(ax.ge.10.)
       ind = ind + 1
       ax = ax/10.d0
    enddo
    if(ind.gt.5 .or. ind.lt.-5)then
       write(Str,*) nint(ax*1.d5)
       Str=trim(adjustl(Str))
       n=verify(trim(str), "0", .true.)
       if(n.gt.1)then
          Coop_real2str=str(1:1)//"."//str(2:n)
       else
          Coop_real2str=Str(1:1)
       endif
       write(str,*) ind
       Coop_real2str=Trim(Coop_real2str)//"e"//Trim(adjustl(Str))
    else
       ax=abs(x)
       write(Str,'(F16.8)') ax
       str=trim(adjustl(str))
       n=scan(str,".")
       if(n>0)then
          n=max(n-1,verify(trim(str),"0",.true.))
       else
          n=Len_Trim(str)
       endif
       if(Str(n:n).eq.'.')then
          Coop_real2str =  Str(1:n-1)
       else
          Coop_real2str = Str(1:n)
       endif
    endif
    if(x.lt.0.d0) Coop_real2str='-'//trim(Coop_real2str)
  End Function coop_real2str


  function Coop_logical2str(x)
    logical x
    COOP_STRING coop_logical2str
    if(x)then
       Coop_logical2str="T"
    else
       Coop_logical2str="F"
    endif
  end function Coop_logical2str


  Function coop_str2int(str)
    COOP_UNKNOWN_STRING  str
    COOP_INT coop_str2int
    read(str,*)coop_str2int
  End Function coop_str2int

  function coop_str2real(str)
    COOP_UNKNOWN_STRING  str
    COOP_REAL coop_str2real
    read(str,*)coop_str2real
  end function coop_str2real

  function coop_str2logical(str)
    COOP_UNKNOWN_STRING str
    logical coop_str2logical
    read(str,*) coop_str2logical
  end function coop_str2logical

!!================= Now the really cool stuff ====================
  function coop_substr(str, pos_start, pos_end, terminator)
    COOP_LONG_STRING coop_substr
    COOP_UNKNOWN_STRING  str
    COOP_INT pos_start, pos_end, i
    character,optional::terminator
    if(present(terminator))then
       do i = pos_start, pos_end
          if(str(i:i) .eq. terminator) exit
       enddo
       coop_substr = str(pos_start:i-1)
    else
       coop_substr = str(pos_start:pos_end)
    endif
  end function coop_substr

  function coop_str_replace(str, needle, replacement) result(repl)
    COOP_UNKNOWN_STRING str, needle, replacement
    COOP_LONG_STRING repl
    COOP_INT i, istart, k, lenr, lenn, lens
    lenr = len(replacement)
    lenn = len(needle)
    lens = len(str)
    repl = ""
    k = 1
    istart = 1
    i = index(str(istart:), needle)
    do while(i.ne.0)
       repl(k:k+i-2+lenr) = str(istart:istart+i-2)//replacement
       k = k + i + lenr-1
       istart = istart + i + lenn-1
       if(istart .gt. lens)exit
       i = index(str(istart:), needle)
    enddo
    repl(k:k+lens - istart) = str(istart:lens)
  end function coop_str_replace

end module coop_string
