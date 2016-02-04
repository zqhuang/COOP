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
  end type coop_math_expression


contains

  !!errors : 
  !! 1 : brackets do not match
  recursive subroutine coop_math_expression_simplify(this, str, length, err) 
    class(coop_math_expression)::this
    character(len=:),allocatable::tmp
    COOP_UNKNOWN_STRING::str
    COOP_INT::length !!length of the str
    COOP_INT err
    COOP_INT::ipos, i, j, nleft, nright, l
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
             err = coop_math_expression_error_brackets_not_match
             return
          endif
          j = i + j -1
          if(str(j:j).eq."(")then
             nleft  = nleft + 1
          else
             nright = nright + 1
          endif
          if(nright .eq. nleft)then
             allocate(character(len=j-ipos-1)::tmp)
             tmp = str(ipos+1:j-1)
             call coop_math_expression_simplify(this, tmp, j-ipos-1, err)
             if(err.ne.0)return
             this%n = this%n + 1
             read(tmp, *) this%vars(this%n)
             if(ipos.ge.4)then  !!check functions
                select case(str(ipos-3:ipos-1))
                case("exp")
                   this%vars(this%n) = exp(this%vars(this%n))
                   ipos = ipos - 3
                case("log")
                   this%vars(this%n) = log(this%vars(this%n))
                   ipos = ipos - 3
                case("sin")
                   if(ipos .ge. 5)then
                      if(str(ipos-4:ipos-4).eq."a")then
                         this%vars(this%n) = asin(this%vars(this%n))                         
                         ipos = ipos - 4
                      else
                         select case(ichar(str(ipos-4:ipos-4))
                         case(coop_ascii_lower_a : coop_ascii_lower_z)
                            error = coop_math_expression_error_unknown_function 
                            return
                         case(coop_ascii_0:coop_ascii_9, coop_ascii_dot)
                            error = coop_math_expression_error_numbers_together
                            return
                         end select
                         this%vars(this%n) = sin(this%vars(this%n))
                         ipos = ipos - 3
                      endif
                   else
                      this%vars(this%n) = sin(this%vars(this%n))
                      ipos = ipos - 3
                   endif
                case("cos")
                   if(ipos.ge.5)then
                      if(str(ipos-4:ipos-4).eq."a")then
                         this%vars(this%n) = acos(this%vars(this%n))                         
                         ipos = ipos - 4
                      else
                         select case(ichar(str(ipos-4:ipos-4))
                         case(coop_ascii_lower_a : coop_ascii_lower_z)
                            error = coop_math_expression_error_unknown_function 
                            return
                         case(coop_ascii_0:coop_ascii_9, coop_ascii_dot)
                            error = coop_math_expression_error_numbers_together
                            return
                         end select
                         this%vars(this%n) = cos(this%vars(this%n))
                         ipos = ipos - 3
                      endif
                   else
                      this%vars(this%n) = cos(this%vars(this%n))
                      ipos = ipos - 3
                   endif
                case("tan")
                   if(ipos.ge.5)then
                      if(str(ipos-4:ipos-4).eq."a")then
                         this%vars(this%n) = atan(this%vars(this%n))                         
                         ipos = ipos - 4
                      else
                         select case(ichar(str(ipos-4:ipos-4))
                         case(coop_ascii_lower_a : coop_ascii_lower_z)
                            error = coop_math_expression_error_unknown_function 
                            return
                         case(coop_ascii_0:coop_ascii_9, coop_ascii_dot)
                            error = coop_math_expression_error_numbers_together
                            return
                         end select
                         this%vars(this%n) = tan(this%vars(this%n))
                         ipos = ipos - 3
                      endif
                   else
                      this%vars(this%n) = tan(this%vars(this%n))
                      ipos = ipos - 3
                   endif
                case("qrt")
                   if(ipos .ge. 5)then
                      if(str(ipos-4:ipos-4).eq."s")then
                         this%vars(this%n) = sqrt(this%vars(this%n))
                         ipos = ipos - 4
                      else
                         error = coop_math_expression_error_unknown_function 
                         return
                      endif
                   else
                      error = coop_math_expression_error_unknown_function 
                      return
                   endif
                case default
                   select case(ichar(str(ipos-4:ipos-4))
                   case(coop_ascii_lower_a : coop_ascii_lower_z)
                      error = coop_math_expression_error_unknown_function 
                      return
                   case(coop_ascii_0:coop_ascii_9, coop_ascii_dot)
                      error = coop_math_expression_error_numbers_together
                      return
                   end select
                end select
             endif
             str = str(1:ipos-1)//"$"//trim(coop_4digits(this%n))//str(j+1:l)
             l = len_trim(str)
             exit
          endif
          i = j + 1
          if(i.gt.l)then
             err = coop_math_expression_error_brackets_not_match
             return
          endif
       enddo
       ipos = scan(str(1:l), "(")
    enddo
    err = 0
    return
  end subroutine coop_math_expression_simplify

  subroutine coop_math_expression_init(this, mathexpr)
    class(coop_math_expression)::this
    COOP_UNKNOWN_STRING::mathexpr
    COOP_LONG_STRING::str
    COOP_INT::i, l, j, k
    l = 0
    do i=1, len_trim(mathexpr)
       select case(ichar(mathexpr(i:i)))
       case(coop_ascii_0:coop_ascii_9, coop_ascii_lower_a:coop_ascii_lower_z, coop_ascii_plus, coop_ascii_dash, coop_ascii_star, coop_ascii_slash, coop_ascii_hat, coop_ascii_dot, coop_ascii_left_bracket, coop_ascii_right_bracket)
          l = l + 1
          str(l:l) = mathexpr(i:i)
       case(coop_ascii_upper_a:coop_ascii_upper_z)
          l=l+1
          str(l:l) = char(ichar(mathexpr(i:i)) + coop_ascii_lower_minus_upper)
       end select
    enddo
    str(l+1:) = ""
    i = 1
    this%n = 0
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
                if(scan(str(j-2:j-2), "0123456789.").eq.0)then !!no number before it
                   j = j -1
                endif
             else
               j = j-1
             endif
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
       read(str(i:j-1), *) this%vars(this%n)
       this%expr = trim(this%expr)//"$"//trim(coop_4digits(this%n))
       i = j
    enddo
  end subroutine coop_math_expression_init


end module coop_evalstr_mod
