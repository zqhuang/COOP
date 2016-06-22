module coop_evalstr_mod
  use coop_wrapper_typedef
  implicit none
!!only support real number expressions
!!support operators: + - * / ^
!!support functions: exp, log, sin, cos, tan, sqrt, asin, acos, atan
!!support constants: pi

#include "constants.h"

  COOP_INT, parameter::coop_math_expression_max_num_vars = 9999
  COOP_INT, parameter::coop_math_expression_error_brackets_not_match = 1
  COOP_INT, parameter::coop_math_expression_error_unknown_function = 2
  COOP_INT, parameter::coop_math_expression_error_operators_together = 3
  COOP_INT, parameter::coop_math_expression_error_numbers_together = 4
  COOP_INT, parameter::coop_math_expression_error_extra_dot = 5
  COOP_INT, parameter::coop_math_expression_error_too_many_numbers = 6
  COOP_INT, parameter::coop_math_expression_error_invalid_characters = 7
  
  type coop_math_expression 
     COOP_LONG_STRING::expr     !!contains ()+-*/^a-z and digits
     !! e is reserved for numbers
     !! x-> exp
     !! l-> ln
     !! r-> sqrt
     !! s-> sin
     !! c-> cos
     !! t-> tan
     !!more mappings to be added
     COOP_INT::n = 0
     COOP_REAL::vars(coop_math_expression_max_num_vars)  
   contains
     procedure::init => coop_math_expression_init
     procedure::simplify => coop_math_expression_simplify
  end type coop_math_expression


contains

  subroutine coop_eval_math(mathexpr, ans, vars)
    type(coop_math_expression)::cme
    COOP_REAL,optional::vars(:)
    COOP_INT::error,length,ind
    COOP_UNKNOWN_STRING::mathexpr
    COOP_REAL::ans
    if(present(vars))then
       call cme%init(mathexpr, error,vars)
    else
       call cme%init(mathexpr, error)
    endif
    if(error.ne.0)then
       write(*,*) trim(mathexpr)
       stop "invalid math expression (init)"
    endif
    length = len_trim(cme%expr)
    call cme%simplify(cme%expr, length, error)
    if(error.ne.0)then
       write(*,*) trim(mathexpr)
       stop "invalid math expression"
    endif
    read(cme%expr(2:5),*)ind
    ans = cme%vars(ind)
  end subroutine coop_eval_math


  subroutine coop_math_expression_init(this, mathexpr, error, vars)
    class(coop_math_expression)::this
    COOP_REAL,optional::vars(:)
    COOP_INT::error
    COOP_UNKNOWN_STRING::mathexpr
    COOP_LONG_STRING::str
    COOP_INT::i, l, j, k, ll, ind, maxind
    if(present(vars))then
       maxind = size(vars)
    else
       maxind = 0
    endif
    this%n = 0
    l = 0
    ll = len_trim(mathexpr)
    error = 0
    i = 1
    
    do while(i.le.ll)
       select case(ichar(mathexpr(i:i)))
       case(coop_ascii_0:coop_ascii_9, coop_ascii_plus, coop_ascii_dash, coop_ascii_slash, coop_ascii_hat, coop_ascii_dot, coop_ascii_left_bracket, coop_ascii_right_bracket, coop_ascii_lower_a:coop_ascii_lower_z)
          l = l + 1
          str(l:l) = mathexpr(i:i)
       case(coop_ascii_upper_a:coop_ascii_upper_z)
          l=l+1
          str(l:l) = char(ichar(mathexpr(i:i)) + coop_ascii_lower_minus_upper)
       case(1:coop_ascii_space)
          !! do nothing
       case(coop_ascii_star)
          if(l.ge.1)then
             if(str(l:l).eq."*")then
                str(l:l) = "^"
             else
                l = l + 1
                str(l:l) = "*"
             endif
          else
             error = coop_math_expression_error_operators_together
             return
          endif
       case(coop_ascii_dollar)
          if(present(vars))then
             l = l + 1
             str(l:l) = "$"
             if(i.ge.ll)then
                error = coop_math_expression_error_invalid_characters
                return
             endif
          else
             error = coop_math_expression_error_invalid_characters
             return
          endif
       case default
          error = coop_math_expression_error_invalid_characters
          return
       end select
       i=i+1
    enddo
    str(l+1:) = ''
    if(index(str, "pi").ne.0)then
       str = coop_str_replace(str(1:l), "pi", "3.14159265359")
       l = len_trim(str)
    endif
    if(index(str, "ln").ne.0)then
       str = coop_str_replace(str(1:l), "ln", "log")
       l = len_trim(str)
    endif
    if(index(str, "log10").ne.0)then
       str = coop_str_replace(str(1:l), "log10", "llo")
       l = len_trim(str)
    endif
    i = 1
    this%expr = ""
    do while(i.le.l)
       j = scan(str(i:l), "0123456789.")
       if(j.eq.0)then
          this%expr = trim(this%expr)//str(i:l)
          return
       endif
       j = i + j - 1
       !!include possible "+/-" sign before the digits
       if(j.gt.1)then
          if(str(j-1:j-1).eq."-" .or. str(j-1:j-1).eq."+")then
             if(j.gt.2)then
                if(scan(str(j-2:j-2), "0123456789.)").eq.0)then !!no number before it
                   j = j -1
                endif
             else
               j = j-1
             endif
          endif
          if(str(j-1:j-1).eq."$")then
             k = j
             do while(coop_is_digit(str(k:k)))
                k = k + 1
                if(k.gt.l)exit
             enddo
             if(k.le.j)then
                error = coop_math_expression_error_invalid_characters
                return
             endif
             read(str(j:k-1),*) ind
             if(ind .gt. maxind)then
                error = coop_math_expression_error_numbers_together
                return
             endif
             this%n = this%n + 1
             if(this%n .gt.  coop_math_expression_max_num_vars)then
                error = coop_math_expression_error_too_many_numbers
                return
             endif
             this%vars(this%n) = vars(ind)
             this%expr = trim(this%expr)//str(i:j-1)//trim(coop_4digits(this%n))
             i = k
             cycle
          endif
       endif
       this%expr = trim(this%expr)//str(i:j-1)
       i = j
       !!now i is the location of the beginning of the number
       j = verify(str(i+1:l), "0123456789.")
       if(j.eq.0 .or. i+j.gt.l)then
          j = l + 1
       else
          j = i+j
          if(str(j:j).eq."e")then
             j = j+ 1
             if(str(j:j).eq."+" .or. str(j:j).eq."-") j = j+1
             do while(coop_is_digit(str(j:j)))
                j = j+1
                if(j.gt.l)exit
             enddo
          endif
       endif
       this%n = this%n+1
       if(this%n .gt. coop_math_expression_max_num_vars)then
          error = coop_math_expression_error_too_many_numbers   
          return
       endif
       read(str(i:j-1), *) this%vars(this%n)
       this%expr = trim(this%expr)//"$"//trim(coop_4digits(this%n))
       i = j
    enddo
    if(scan(this%expr, ".") .ne. 0)then
       error = coop_math_expression_error_extra_dot
       return
    endif
  end subroutine coop_math_expression_init



  !!errors : 
  recursive subroutine coop_math_expression_simplify(this, str, length, error) 
    class(coop_math_expression)::this
    character(len=:),allocatable::tmp
    COOP_UNKNOWN_STRING::str
    COOP_INT::length !!length of the str
    COOP_INT error
    COOP_INT::ipos, i, j, nleft, nright, l, ind
    l = length
    !!simplify all brackets
    ipos = scan(str(1:l), "(")
    do while(ipos .ne. 0)
       nleft = 1
       nright = 0
       i = ipos + 1
       do
          j = scan(str(i:l),"()")
          if(j.eq.0)then
             error = coop_math_expression_error_brackets_not_match
             return
          endif
          j = i + j -1
          if(str(j:j).eq."(")then
             nleft  = nleft + 1
          else
             nright = nright + 1
          endif
          if(nright .eq. nleft)then
            ! allocate(character(len=j-ipos-1)::tmp) !!seems not necessary; fortran does automatic allocation
             tmp = str(ipos+1:j-1)
             call coop_math_expression_simplify(this, tmp, j-ipos-1, error)
             if(error.ne.0)return
             read(tmp(2:5),*) ind
             if(ipos.ge.4)then  !!check functions
                select case(str(ipos-3:ipos-1))
                case("exp")
                   this%vars(ind) = exp(this%vars(ind))
                   ipos = ipos - 3
                case("log")
                   this%vars(ind) = log(this%vars(ind))
                   ipos = ipos - 3
                case("llo")
                   this%vars(ind) = log10(this%vars(ind))
                   ipos = ipos - 3
                case("sin")
                   if(ipos .ge. 5)then
                      if(str(ipos-4:ipos-4).eq."a")then
                         this%vars(ind) = asin(this%vars(ind))                         
                         ipos = ipos - 4
                      else
                         select case(ichar(str(ipos-4:ipos-4)))
                         case(coop_ascii_lower_a : coop_ascii_lower_z)
                            error = coop_math_expression_error_unknown_function 
                            return
                         case(coop_ascii_0:coop_ascii_9)
                            error = coop_math_expression_error_numbers_together
                            return
                         end select
                         this%vars(ind) = sin(this%vars(ind))
                         ipos = ipos - 3
                      endif
                   else
                      this%vars(ind) = sin(this%vars(ind))
                      ipos = ipos - 3
                   endif
                case("cos")
                   if(ipos.ge.5)then
                      if(str(ipos-4:ipos-4).eq."a")then
                         this%vars(ind) = acos(this%vars(ind))                         
                         ipos = ipos - 4
                      else
                         select case(ichar(str(ipos-4:ipos-4)))
                         case(coop_ascii_lower_a : coop_ascii_lower_z)
                            error = coop_math_expression_error_unknown_function 
                            return
                         case(coop_ascii_0:coop_ascii_9)
                            error = coop_math_expression_error_numbers_together
                            return
                         end select
                         this%vars(ind) = cos(this%vars(ind))
                         ipos = ipos - 3
                      endif
                   else
                      this%vars(ind) = cos(this%vars(ind))
                      ipos = ipos - 3
                   endif
                case("tan")
                   if(ipos.ge.5)then
                      if(str(ipos-4:ipos-4).eq."a")then
                         this%vars(ind) = atan(this%vars(ind))                         
                         ipos = ipos - 4
                      else
                         select case(ichar(str(ipos-4:ipos-4)))
                         case(coop_ascii_lower_a : coop_ascii_lower_z)
                            error = coop_math_expression_error_unknown_function 
                            return
                         case(coop_ascii_0:coop_ascii_9)
                            error = coop_math_expression_error_numbers_together
                            return
                         end select
                         this%vars(ind) = tan(this%vars(ind))
                         ipos = ipos - 3
                      endif
                   else
                      this%vars(ind) = tan(this%vars(ind))
                      ipos = ipos - 3
                   endif
                case("qrt")
                   if(ipos .ge. 5)then
                      if(str(ipos-4:ipos-4).eq."s")then
                         this%vars(ind) = sqrt(this%vars(ind))
                         ipos = ipos - 4
                      else
                         error = coop_math_expression_error_unknown_function 
                         return
                      endif
                   else
                      error = coop_math_expression_error_unknown_function 
                      return
                   endif
                case("mma")
                   if(ipos .ge. 6)then
                      if(str(ipos-5:ipos-4).eq."ga")then
                         this%vars(ind) = gamma(this%vars(ind))
                         ipos = ipos - 5
                      else
                         error = coop_math_expression_error_unknown_function 
                         return
                      endif
                   else
                      error = coop_math_expression_error_unknown_function 
                      return
                   endif
                case default
                   select case(ichar(str(ipos-4:ipos-4)))
                   case(coop_ascii_lower_a : coop_ascii_lower_z)
                      error = coop_math_expression_error_unknown_function 
                      return
                   case(coop_ascii_0:coop_ascii_9)
                      error = coop_math_expression_error_numbers_together
                      return
                   end select
                end select
             endif
             str = str(1:ipos-1)//tmp(1:5)//str(j+1:l)
             l = l+(ipos+4-j)
             deallocate(tmp)
             exit
          endif
          i = j + 1
          if(i.gt.l)then
             error = coop_math_expression_error_brackets_not_match
             return
          endif
       enddo
       ipos = scan(str(1:l), "(")
    enddo

    ipos = scan(str(1:l), "^")
    do while(ipos .ne. 0)
       if(ipos .lt. 6 .or. ipos + 5 .gt. l)then
          error = coop_math_expression_error_operators_together
          return
       endif
       if(str(ipos-5:ipos-5) .ne. "$" .or. str(ipos+1:ipos+1).ne."$")then
          error = coop_math_expression_error_operators_together
          return
       endif
       read(str(ipos-4:ipos-1),*) nleft
       read(str(ipos+2:ipos+5),*) nright
       this%vars(nleft) = this%vars(nleft)**this%vars(nright)
       str = str(1:ipos-1)//str(ipos+6:l)
       l = l - 6
       ipos = scan(str(1:l), "^")
    enddo

    ipos = scan(str(1:l), "*/")
    do while(ipos .ne. 0)
       if(ipos .lt. 6 .or. ipos + 5 .gt. l)then
          error = coop_math_expression_error_operators_together
          return
       endif
       if(str(ipos-5:ipos-5) .ne. "$" .or. str(ipos+1:ipos+1).ne."$")then
          error = coop_math_expression_error_operators_together
          return
       endif
       read(str(ipos-4:ipos-1),*) nleft
       read(str(ipos+2:ipos+5),*) nright
       if(str(ipos:ipos).eq."*")then
          this%vars(nleft) = this%vars(nleft)*this%vars(nright)
       else
          this%vars(nleft) = this%vars(nleft)/this%vars(nright)
       endif
       str = str(1:ipos-1)//str(ipos+6:l)
       l = l - 6
       ipos = scan(str(1:l), "*/")
    enddo


    ipos = scan(str(1:l), "+-")
    do while(ipos .ne. 0)
       if(ipos .lt. 6 .or. ipos + 5 .gt. l)then
          if(ipos .eq. 1 .and. 6 .le. l)then
             if(str(ipos+1:ipos+1).ne."$")then
                error = coop_math_expression_error_operators_together
                return
             endif
             read(str(ipos+2:ipos+5),*) nright
             if(str(ipos:ipos).eq."-") this%vars(nright) = - this%vars(nright) 
             str = str(2:l)
             l = l -1
          else
             error = coop_math_expression_error_operators_together
             return
          endif
       else
          if(str(ipos-5:ipos-5) .ne. "$" .or. str(ipos+1:ipos+1).ne."$")then
             error = coop_math_expression_error_operators_together
             return
          endif
          read(str(ipos-4:ipos-1),*) nleft
          read(str(ipos+2:ipos+5),*) nright
          if(str(ipos:ipos).eq."+")then
             this%vars(nleft) = this%vars(nleft)+this%vars(nright)
          else
             this%vars(nleft) = this%vars(nleft)-this%vars(nright)
          endif
          str = str(1:ipos-1)//str(ipos+6:l)
          l = l - 6
       endif
       ipos = scan(str(1:l), "+-")
    enddo

    error = 0
    return
  end subroutine coop_math_expression_simplify

end module coop_evalstr_mod
