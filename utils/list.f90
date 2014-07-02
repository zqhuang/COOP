module coop_list_mod
  use coop_wrapper_typedef
  implicit none

#include "constants.h"
  private

  COOP_INT,parameter::sp = kind(1.)
  COOP_INT,parameter::dl = kind(1.d0)

  public::coop_list_integer, coop_list_real, coop_list_realarr, coop_list_double, coop_list_logical, coop_list_string, coop_list_character, coop_string_to_list, coop_dictionary, coop_dictionary_lookup, coop_get_prime_numbers, coop_list_get_element

  interface coop_list_initialize
     module procedure coop_list_integer_initialize, coop_list_real_initialize, coop_list_double_initialize, coop_list_logical_initialize, coop_list_string_initialize, coop_list_character_initialize, coop_list_realarr_initialize
  end interface coop_list_initialize


  interface coop_list_is_initialized
     module procedure coop_list_integer_is_initialized, coop_list_real_is_initialized, coop_list_double_is_initialized, coop_list_logical_is_initialized, coop_list_string_is_initialized, coop_list_character_is_initialized, coop_list_realarr_is_initialized
  end interface coop_list_is_initialized

  interface coop_list_push
     module procedure coop_list_integer_push, coop_list_real_push, coop_list_double_push, coop_list_logical_push, coop_list_string_push, coop_list_character_push, coop_list_realarr_push
  end interface coop_list_push


  interface coop_list_pop
     module procedure coop_list_integer_pop, coop_list_real_pop, coop_list_double_pop, coop_list_logical_pop, coop_list_string_pop, coop_list_character_pop, coop_list_realarr_pop
  end interface coop_list_pop

  interface coop_list_element
     module procedure coop_list_integer_element, coop_list_real_element, coop_list_double_element, coop_list_logical_element, coop_list_string_element, coop_list_character_element, coop_list_realarr_element
  end interface coop_list_element

  interface coop_list_get_element
     module procedure coop_list_integer_get_element, coop_list_real_get_element, coop_list_double_get_element, coop_list_logical_get_element, coop_list_string_get_element, coop_list_character_get_element, coop_list_realarr_get_element
  end interface coop_list_get_element


  interface coop_string_to_list
     module procedure coop_string_to_list_string, coop_string_to_list_real, coop_string_to_list_double, coop_string_to_list_integer, coop_string_to_list_logical
  end interface coop_string_to_list

  interface coop_dictionary_lookup
     module procedure coop_dictionary_lookup_string, coop_dictionary_lookup_int, coop_dictionary_lookup_real, coop_dictionary_lookup_double, coop_dictionary_lookup_logical, coop_dictionary_lookup_int_array, coop_dictionary_lookup_real_array, coop_dictionary_lookup_double_array, coop_dictionary_lookup_logical_array
  end interface coop_dictionary_lookup


  COOP_INT,parameter:: coop_list_i1_max_length = 2**12
  COOP_INT,parameter:: coop_list_i2_max_length = coop_list_i1_max_length * 4
  COOP_INT,parameter:: coop_list_i3_max_length = coop_list_i2_max_length * 4
  COOP_INT,parameter:: coop_list_i4_max_length = coop_list_i3_max_length * 4


  COOP_INT,parameter:: coop_list_s1_max_length =  coop_list_i1_max_length/coop_string_length
  COOP_INT,parameter:: coop_list_s2_max_length = coop_list_s1_max_length * 4
  COOP_INT,parameter:: coop_list_s3_max_length = coop_list_s2_max_length * 4
  COOP_INT,parameter:: coop_list_s4_max_length = coop_list_s3_max_length * 4

  COOP_INT,parameter::coop_list_unit_len = 8192

  type coop_list_integer
     COOP_INT n, stack, loc
     COOP_INT,dimension(:),allocatable::i1
     COOP_INT,dimension(:),allocatable::i2
     COOP_INT,dimension(:),allocatable::i3
     COOP_INT,dimension(:),allocatable::i4
   contains
     procedure::init => coop_list_integer_initialize
     procedure::isinit => coop_list_integer_is_initialized
     procedure::push => coop_list_integer_push
     procedure::pop => coop_list_integer_pop
     procedure::element => coop_list_integer_element
     procedure::get_element => coop_list_integer_get_element
  end type coop_list_integer

  type coop_list_real
     COOP_INT n, stack, loc
     real(sp),dimension(:),allocatable::i1
     real(sp),dimension(:),allocatable::i2
     real(sp),dimension(:),allocatable::i3
     real(sp),dimension(:),allocatable::i4
   contains
     procedure::init => coop_list_real_initialize
     procedure::isinit => coop_list_real_is_initialized
     procedure::push => coop_list_real_push
     procedure::pop => coop_list_real_pop
     procedure::element => coop_list_real_element
     procedure::get_element => coop_list_real_get_element
  end type coop_list_real

  type coop_list_realarr
     COOP_INT dim
     COOP_INT n, stack, loc
     real(sp),dimension(:,:),allocatable::i1
     real(sp),dimension(:,:),allocatable::i2
     real(sp),dimension(:,:),allocatable::i3
     real(sp),dimension(:,:),allocatable::i4
   contains
     procedure::init => coop_list_realarr_initialize
     procedure::isinit => coop_list_realarr_is_initialized
     procedure::push => coop_list_realarr_push
     procedure::pop => coop_list_realarr_pop
     procedure::element => coop_list_realarr_element
     procedure::get_element => coop_list_realarr_get_element
  end type coop_list_realarr

  type coop_list_double
     COOP_INT n, stack, loc
     real(dl),dimension(:),allocatable::i1
     real(dl),dimension(:),allocatable::i2
     real(dl),dimension(:),allocatable::i3
     real(dl),dimension(:),allocatable::i4
   contains
     procedure::isinit => coop_list_double_is_initialized
     procedure::init => coop_list_double_initialize
     procedure::push => coop_list_double_push
     procedure::pop => coop_list_double_pop
     procedure::element => coop_list_double_element
     procedure::get_element => coop_list_double_get_element
  end type coop_list_double

  type coop_list_logical
     COOP_INT n, stack, loc
     logical,dimension(:),allocatable::i1
     logical,dimension(:),allocatable::i2
     logical,dimension(:),allocatable::i3
     logical,dimension(:),allocatable::i4
   contains
     procedure::isinit => coop_list_logical_is_initialized
     procedure::init => coop_list_logical_initialize
     procedure::push => coop_list_logical_push
     procedure::pop => coop_list_logical_pop
     procedure::element => coop_list_logical_element
     procedure::get_element => coop_list_logical_get_element
  end type coop_list_logical

  type coop_list_character
     COOP_INT n, stack, loc
     character,dimension(:),allocatable::i1
     character,dimension(:),allocatable::i2
     character,dimension(:),allocatable::i3
     character,dimension(:),allocatable::i4
   contains
     procedure::isinit => coop_list_character_is_initialized
     procedure::init => coop_list_character_initialize
     procedure::push => coop_list_character_push
     procedure::pop => coop_list_character_pop
     procedure::element => coop_list_character_element
     procedure::get_element => coop_list_character_get_element
  end type coop_list_character

  type coop_list_string
     COOP_INT n, stack, loc
     COOP_STRING,dimension(:),allocatable::i1
     COOP_STRING,dimension(:),allocatable::i2
     COOP_STRING,dimension(:),allocatable::i3
     COOP_STRING,dimension(:),allocatable::i4
   contains
     procedure::isinit => coop_list_string_is_initialized
     procedure::init => coop_list_string_initialize
     procedure::push => coop_list_string_push
     procedure::pop => coop_list_string_pop
     procedure::element => coop_list_string_element
     procedure::get_element => coop_list_string_get_element
  end type coop_list_string

  type coop_dictionary
     COOP_INT n, capacity
     COOP_SHORT_STRING,dimension(:),allocatable::key
     COOP_STRING,dimension(:),allocatable::val
     COOP_INT, dimension(:),allocatable::id
   contains
     procedure::print => coop_dictionary_print
     procedure::insert => coop_dictionary_insert
     procedure::index => coop_dictionary_key_index
     procedure::value => coop_dictionary_value
     procedure::update => coop_dictionary_update
     procedure::free => coop_dictionary_free
  end type coop_dictionary


contains

!! initialize


  subroutine coop_list_string_initialize(l)
    class(coop_list_string) l
    l%n = 0
    l%stack = 1
    l%loc = 0
    if(allocated(l%i4))deallocate(l%i4)
    if(allocated(l%i3))deallocate(l%i3)
    if(allocated(l%i2))deallocate(l%i2)
    if(allocated(l%i1))deallocate(l%i1)
  end subroutine coop_list_string_initialize



  subroutine coop_list_real_initialize(l)
    class(coop_list_real) l
    l%n = 0
    l%stack = 1
    l%loc = 0
    if(allocated(l%i4))deallocate(l%i4)
    if(allocated(l%i3))deallocate(l%i3)
    if(allocated(l%i2))deallocate(l%i2)
    if(allocated(l%i1))deallocate(l%i1)
  end subroutine coop_list_real_initialize



  subroutine coop_list_realarr_initialize(l)
    class(coop_list_realarr) l
    l%n = 0
    l%stack = 1
    l%loc = 0
    if(allocated(l%i4))deallocate(l%i4)
    if(allocated(l%i3))deallocate(l%i3)
    if(allocated(l%i2))deallocate(l%i2)
    if(allocated(l%i1))deallocate(l%i1)
  end subroutine coop_list_realarr_initialize


  subroutine coop_list_double_initialize(l)
    class(coop_list_double) l
    l%n = 0
    l%stack = 1
    l%loc = 0
    if(allocated(l%i4))deallocate(l%i4)
    if(allocated(l%i3))deallocate(l%i3)
    if(allocated(l%i2))deallocate(l%i2)
    if(allocated(l%i1))deallocate(l%i1)
  end subroutine coop_list_double_initialize


  subroutine coop_list_integer_initialize(l)
    class(coop_list_integer) l
    l%n = 0
    l%stack = 1
    l%loc = 0
    if(allocated(l%i4))deallocate(l%i4)
    if(allocated(l%i3))deallocate(l%i3)
    if(allocated(l%i2))deallocate(l%i2)
    if(allocated(l%i1))deallocate(l%i1)
  end subroutine coop_list_integer_initialize


  subroutine coop_list_logical_initialize(l)
    class(coop_list_logical) l
    l%n = 0
    l%stack = 1
    l%loc = 0
    if(allocated(l%i4))deallocate(l%i4)
    if(allocated(l%i3))deallocate(l%i3)
    if(allocated(l%i2))deallocate(l%i2)
    if(allocated(l%i1))deallocate(l%i1)
  end subroutine coop_list_logical_initialize


  subroutine coop_list_character_initialize(l)
    class(coop_list_character) l
    l%n = 0
    l%stack = 1
    l%loc = 0
    if(allocated(l%i4))deallocate(l%i4)
    if(allocated(l%i3))deallocate(l%i3)
    if(allocated(l%i2))deallocate(l%i2)
    if(allocated(l%i1))deallocate(l%i1)
  end subroutine coop_list_character_initialize


  function coop_list_integer_is_initialized(l) result(ini)
    class(coop_list_integer) l
    logical ini
    ini = allocated(l%i1)
  end function coop_list_integer_is_initialized


  function coop_list_real_is_initialized(l) result(ini)
    class(coop_list_real) l
    logical ini
    ini = allocated(l%i1)
  end function coop_list_real_is_initialized


  function coop_list_double_is_initialized(l) result(ini)
    class(coop_list_double) l
    logical ini
    ini = allocated(l%i1)
  end function coop_list_double_is_initialized


  function coop_list_character_is_initialized(l) result(ini)
    class(coop_list_character) l
    logical ini
    ini = allocated(l%i1)
  end function coop_list_character_is_initialized

  function coop_list_logical_is_initialized(l) result(ini)
    class(coop_list_logical) l
    logical ini
    ini = allocated(l%i1)
  end function coop_list_logical_is_initialized

  function coop_list_string_is_initialized(l) result(ini)
    class(coop_list_string) l
    logical ini
    ini = allocated(l%i1)
  end function coop_list_string_is_initialized

  function coop_list_realarr_is_initialized(l) result(ini)
    class(coop_list_realarr) l
    logical ini
    ini = allocated(l%i1)
  end function coop_list_realarr_is_initialized


!!push and pop
  subroutine coop_list_integer_push(l, i)
    class(coop_list_integer) l
    COOP_INT i
    if(allocated(l%i1))then
       l%n = l%n+1
       select case(l%stack)
       case(1)
          if(l%loc .lt. coop_list_i1_max_length)then
             l%loc = l%loc + 1
             l%i1(l%loc) = i
             return
          else
             l%loc = 1
             l%stack = 2
             allocate(l%i2(coop_list_i2_max_length))
             l%i2(1) = i
             return
          endif
       case(2)
          if(l%loc .lt. coop_list_i2_max_length)then
             l%loc = l%loc + 1
             l%i2(l%loc) = i
             return
          else
             l%loc = 1
             l%stack = 3
             allocate(l%i3(coop_list_i3_max_length))
             l%i3(1) = i
             return
          endif
       case(3)
          if(l%loc .lt. coop_list_i3_max_length)then
             l%loc = l%loc + 1
             l%i3(l%loc) = i
             return
          else
             l%loc = 1
             l%stack = 4
             allocate(l%i4(coop_list_i4_max_length))
             l%i4(1) = i
             return
          endif
       case(4)
          if(l%loc .lt. coop_list_i4_max_length)then
             l%loc = l%loc + 1
             l%i4(l%loc) = i
             return
          else
             write(*,*) "coop_list_integer_push: Coop_list over flow"
             stop
          endif
       end select
    else
       allocate(l%i1(coop_list_i1_max_length))
       l%n = 1
       l%stack = 1
       l%loc = 1
       l%i1(1) = i
    endif
  end subroutine coop_list_integer_push

  subroutine coop_list_integer_pop(l)
    class(coop_list_integer) l
    if(allocated(l%i1))then
       l%n = l%n-1
       if(l%loc.gt.1)then
          l%loc = l%loc-1
          return
       else
          l%stack = l%stack - 1
          select case(l%stack)
          case(0)
             l%loc = 0
             deallocate(l%i1)
          case(1)
             l%loc = coop_list_i1_max_length
             deallocate(l%i2)
          case(2)
             l%loc = coop_list_i2_max_length
             deallocate(l%i3)
          case(3)
             l%loc = coop_list_i3_max_length
             deallocate(l%i4)
          end select
       endif
    endif
  end subroutine coop_list_integer_pop



  subroutine coop_list_real_push(l, i)
    class(coop_list_real) l
    real(sp) i
    if(allocated(l%i1))then
       l%n = l%n+1
       select case(l%stack)
       case(1)
          if(l%loc .lt. coop_list_i1_max_length)then
             l%loc = l%loc + 1
             l%i1(l%loc) = i
             return
          else
             l%loc = 1
             l%stack = 2
             allocate(l%i2(coop_list_i2_max_length))
             l%i2(1) = i
             return
          endif
       case(2)
          if(l%loc .lt. coop_list_i2_max_length)then
             l%loc = l%loc + 1
             l%i2(l%loc) = i
             return
          else
             l%loc = 1
             l%stack = 3
             allocate(l%i3(coop_list_i3_max_length))
             l%i3(1) = i
             return
          endif
       case(3)
          if(l%loc .lt. coop_list_i3_max_length)then
             l%loc = l%loc + 1
             l%i3(l%loc) = i
             return
          else
             l%loc = 1
             l%stack = 4
             allocate(l%i4(coop_list_i4_max_length))
             l%i4(1) = i
             return
          endif
       case(4)
          if(l%loc .lt. coop_list_i4_max_length)then
             l%loc = l%loc + 1
             l%i4(l%loc) = i
             return
          else
             write(*,*) "coop_list_real_push: Coop_list over flow"
             stop
          endif
       end select
    else
       allocate(l%i1(coop_list_i1_max_length))
       l%n = 1
       l%stack = 1
       l%loc = 1
       l%i1(1) = i
    endif
  end subroutine coop_list_real_push

  subroutine coop_list_real_pop(l)
    class(coop_list_real) l
    if(allocated(l%i1))then
       l%n = l%n-1
       if(l%loc.gt.1)then
          l%loc = l%loc-1
          return
       else
          l%stack = l%stack - 1
          select case(l%stack)
          case(0)
             l%loc = 0
             deallocate(l%i1)
          case(1)
             l%loc = coop_list_i1_max_length
             deallocate(l%i2)
          case(2)
             l%loc = coop_list_i2_max_length
             deallocate(l%i3)
          case(3)
             l%loc = coop_list_i3_max_length
             deallocate(l%i4)
          end select
       endif
    endif
  end subroutine coop_list_real_pop



  subroutine coop_list_realarr_push(l, i)
    class(coop_list_realarr) l
    real(sp),dimension(:),intent(IN)::i
    if(allocated(l%i1))then
       if(size(i).ne. l%dim) stop "coop_list_realarr_push: wrong size of input array"
       l%n = l%n+1
       select case(l%stack)
       case(1)
          if(l%loc .lt. coop_list_i1_max_length)then
             l%loc = l%loc + 1
             l%i1(:,l%loc) = i
             return
          else
             l%loc = 1
             l%stack = 2
             allocate(l%i2(l%dim, coop_list_i2_max_length))
             l%i2(:,1) = i
             return
          endif
       case(2)
          if(l%loc .lt. coop_list_i2_max_length)then
             l%loc = l%loc + 1
             l%i2(:,l%loc) = i
             return
          else
             l%loc = 1
             l%stack = 3
             allocate(l%i3(l%dim, coop_list_i3_max_length))
             l%i3(:,1) = i
             return
          endif
       case(3)
          if(l%loc .lt. coop_list_i3_max_length)then
             l%loc = l%loc + 1
             l%i3(:,l%loc) = i
             return
          else
             l%loc = 1
             l%stack = 4
             allocate(l%i4(l%dim, coop_list_i4_max_length))
             l%i4(:, 1) = i
             return
          endif
       case(4)
          if(l%loc .lt. coop_list_i4_max_length)then
             l%loc = l%loc + 1
             l%i4(:, l%loc) = i
             return
          else
             write(*,*) "coop_list_realarr_push: Coop_list over flow"
             stop
          endif
       end select
    else
       l%dim = size(i)
       allocate(l%i1(l%dim, coop_list_i1_max_length))
       l%n = 1
       l%stack = 1
       l%loc = 1
       l%i1(:,1) = i
    endif
  end subroutine coop_list_realarr_push

  subroutine coop_list_realarr_pop(l)
    class(coop_list_realarr) l
    if(allocated(l%i1))then
       l%n = l%n-1
       if(l%loc.gt.1)then
          l%loc = l%loc-1
          return
       else
          l%stack = l%stack - 1
          select case(l%stack)
          case(0)
             l%loc = 0
             deallocate(l%i1)
          case(1)
             l%loc = coop_list_i1_max_length
             deallocate(l%i2)
          case(2)
             l%loc = coop_list_i2_max_length
             deallocate(l%i3)
          case(3)
             l%loc = coop_list_i3_max_length
             deallocate(l%i4)
          end select
       endif
    endif
  end subroutine coop_list_realarr_pop


  subroutine coop_list_double_push(l, i)
    class(coop_list_double) l
    real(dl) i
    if(allocated(l%i1))then
       l%n = l%n+1
       select case(l%stack)
       case(1)
          if(l%loc .lt. coop_list_i1_max_length)then
             l%loc = l%loc + 1
             l%i1(l%loc) = i
             return
          else
             l%loc = 1
             l%stack = 2
             allocate(l%i2(coop_list_i2_max_length))
             l%i2(1) = i
             return
          endif
       case(2)
          if(l%loc .lt. coop_list_i2_max_length)then
             l%loc = l%loc + 1
             l%i2(l%loc) = i
             return
          else
             l%loc = 1
             l%stack = 3
             allocate(l%i3(coop_list_i3_max_length))
             l%i3(1) = i
             return
          endif
       case(3)
          if(l%loc .lt. coop_list_i3_max_length)then
             l%loc = l%loc + 1
             l%i3(l%loc) = i
             return
          else
             l%loc = 1
             l%stack = 4
             allocate(l%i4(coop_list_i4_max_length))
             l%i4(1) = i
             return
          endif
       case(4)
          if(l%loc .lt. coop_list_i4_max_length)then
             l%loc = l%loc + 1
             l%i4(l%loc) = i
             return
          else
             write(*,*) "coop_list_double_push: Coop_list over flow"
             stop
          endif
       end select
    else
       allocate(l%i1(coop_list_i1_max_length))
       l%n = 1
       l%stack = 1
       l%loc = 1
       l%i1(1) = i
    endif
  end subroutine coop_list_double_push

  subroutine coop_list_double_pop(l)
    class(coop_list_double) l
    if(allocated(l%i1))then
       l%n = l%n-1
       if(l%loc.gt.1)then
          l%loc = l%loc-1
          return
       else
          l%stack = l%stack - 1
          select case(l%stack)
          case(0)
             l%loc = 0
             deallocate(l%i1)
          case(1)
             l%loc = coop_list_i1_max_length
             deallocate(l%i2)
          case(2)
             l%loc = coop_list_i2_max_length
             deallocate(l%i3)
          case(3)
             l%loc = coop_list_i3_max_length
             deallocate(l%i4)
          end select
       endif
    endif
  end subroutine coop_list_double_pop

  subroutine coop_list_logical_push(l, i)
    class(coop_list_logical) l
    logical i
    if(allocated(l%i1))then
       l%n = l%n+1
       select case(l%stack)
       case(1)
          if(l%loc .lt. coop_list_i1_max_length)then
             l%loc = l%loc + 1
             l%i1(l%loc) = i
             return
          else
             l%loc = 1
             l%stack = 2
             allocate(l%i2(coop_list_i2_max_length))
             l%i2(1) = i
             return
          endif
       case(2)
          if(l%loc .lt. coop_list_i2_max_length)then
             l%loc = l%loc + 1
             l%i2(l%loc) = i
             return
          else
             l%loc = 1
             l%stack = 3
             allocate(l%i3(coop_list_i3_max_length))
             l%i3(1) = i
             return
          endif
       case(3)
          if(l%loc .lt. coop_list_i3_max_length)then
             l%loc = l%loc + 1
             l%i3(l%loc) = i
             return
          else
             l%loc = 1
             l%stack = 4
             allocate(l%i4(coop_list_i4_max_length))
             l%i4(1) = i
             return
          endif
       case(4)
          if(l%loc .lt. coop_list_i4_max_length)then
             l%loc = l%loc + 1
             l%i4(l%loc) = i
             return
          else
             write(*,*) "coop_list_logical_push: Coop_list over flow"
             stop
          endif
       end select
    else
       allocate(l%i1(coop_list_i1_max_length))
       l%n = 1
       l%stack = 1
       l%loc = 1
       l%i1(1) = i
    endif
  end subroutine coop_list_logical_push

  subroutine coop_list_logical_pop(l)
    class(coop_list_logical) l
    if(allocated(l%i1))then
       l%n = l%n-1
       if(l%loc.gt.1)then
          l%loc = l%loc-1
          return
       else
          l%stack = l%stack - 1
          select case(l%stack)
          case(0)
             l%loc = 0
             deallocate(l%i1)
          case(1)
             l%loc = coop_list_i1_max_length
             deallocate(l%i2)
          case(2)
             l%loc = coop_list_i2_max_length
             deallocate(l%i3)
          case(3)
             l%loc = coop_list_i3_max_length
             deallocate(l%i4)
          end select
       endif
    endif
  end subroutine coop_list_logical_pop


  recursive subroutine coop_list_character_push(l, i)
    class(coop_list_character) l
    COOP_UNKNOWN_STRING i
    COOP_INT slen, cut
    slen = len(i)
    if(slen .eq. 0) return
    if(allocated(l%i1))then
       l%n = l%n+1
       select case(l%stack)
       case(1)
          if(l%loc + slen .le. coop_list_i1_max_length)then
             l%i1(l%loc+1:l%loc+slen) = i
             l%loc = l%loc + slen
             return
          else
             cut = coop_list_i1_max_length - l%loc
             if(cut.ne.0)l%i1(l%loc+1:coop_list_i1_max_length) = i(1 : cut)
             l%loc = 0
             l%stack = 2
             allocate(l%i2(coop_list_i2_max_length))
             call coop_list_character_push(l, i(cut+1:slen))
             return
          endif
       case(2)
          if(l%loc + slen .le. coop_list_i2_max_length)then
             l%i2(l%loc+1:l%loc+slen) = i
             l%loc = l%loc + slen
             return
          else
             cut = coop_list_i2_max_length - l%loc
             if(cut.ne.0)l%i2(l%loc+1:coop_list_i2_max_length) = i(1 : cut)
             l%loc = 0
             l%stack = 3
             allocate(l%i3(coop_list_i3_max_length))
             call coop_list_character_push(l, i(cut+1:slen))
             return
          endif
       case(3)
          if(l%loc + slen .le. coop_list_i3_max_length)then
             l%i3(l%loc+1:l%loc+slen) = i
             l%loc = l%loc + slen
             return
          else
             cut = coop_list_i3_max_length - l%loc
             if(cut.ne.0)l%i3(l%loc+1:coop_list_i3_max_length) = i(1 : cut)
             l%loc = 0
             l%stack = 4
             allocate(l%i4(coop_list_i4_max_length))
             call coop_list_character_push(l, i(cut+1:slen))
             return
          endif
       case(4)
          if(l%loc+slen .le. coop_list_i4_max_length)then
             l%i4(l%loc+1:l%loc+slen) = i
             l%loc = l%loc + slen
             return
          else
             write(*,*) "coop_list_character_push_multiple: Coop_list over flow"
             stop
          endif
       end select
    else
       allocate(l%i1(coop_list_i1_max_length))
       l%n = 0
       l%stack = 1
       l%loc = 0
       call coop_list_character_push(l, i)
    endif
  end subroutine coop_list_character_push


  subroutine coop_list_character_pop(l)
    class(coop_list_character) l
    if(allocated(l%i1))then
       l%n = l%n-1
       if(l%loc.gt.1)then
          l%loc = l%loc-1
          return
       else
          l%stack = l%stack - 1
          select case(l%stack)
          case(0)
             l%loc = 0
             deallocate(l%i1)
          case(1)
             l%loc = coop_list_i1_max_length
             deallocate(l%i2)
          case(2)
             l%loc = coop_list_i2_max_length
             deallocate(l%i3)
          case(3)
             l%loc = coop_list_i3_max_length
             deallocate(l%i4)
          end select
       endif
    endif
  end subroutine coop_list_character_pop

  subroutine coop_list_string_push(l, i)
    class(coop_list_string) l
    COOP_UNKNOWN_STRING i
    if(allocated(l%i1))then
       l%n = l%n+1
       select case(l%stack)
       case(1)
          if(l%loc .lt. coop_list_s1_max_length)then
             l%loc = l%loc + 1
             l%i1(l%loc) = i
             return
          else
             l%loc = 1
             l%stack = 2
             allocate(l%i2(coop_list_s2_max_length))
             l%i2(1) = i
             return
          endif
       case(2)
          if(l%loc .lt. coop_list_s2_max_length)then
             l%loc = l%loc + 1
             l%i2(l%loc) = i
             return
          else
             l%loc = 1
             l%stack = 3
             allocate(l%i3(coop_list_s3_max_length))
             l%i3(1) = i
             return
          endif
       case(3)
          if(l%loc .lt. coop_list_s3_max_length)then
             l%loc = l%loc + 1
             l%i3(l%loc) = i
             return
          else
             l%loc = 1
             l%stack = 4
             allocate(l%i4(coop_list_s4_max_length))
             l%i4(1) = i
             return
          endif
       case(4)
          if(l%loc .lt. coop_list_s4_max_length)then
             l%loc = l%loc + 1
             l%i4(l%loc) = i
             return
          else
             write(*,*) "coop_list_string_push: Coop_list over flow"
             stop
          endif
       end select
    else
       allocate(l%i1(coop_list_s1_max_length))
       l%n = 1
       l%stack = 1
       l%loc = 1
       l%i1(1) = i
    endif
  end subroutine coop_list_string_push

  subroutine coop_list_string_pop(l)
    class(coop_list_string) l
    if(allocated(l%i1))then
       l%n = l%n-1
       if(l%loc.gt.1)then
          l%loc = l%loc-1
          return
       else
          l%stack = l%stack - 1
          select case(l%stack)
          case(0)
             l%loc = 0
             deallocate(l%i1)
          case(1)
             l%loc = coop_list_s1_max_length
             deallocate(l%i2)
          case(2)
             l%loc = coop_list_s2_max_length
             deallocate(l%i3)
          case(3)
             l%loc = coop_list_s3_max_length
             deallocate(l%i4)
          end select
       endif
    endif
  end subroutine coop_list_string_pop

!!dynamic string operations  
  subroutine coop_string_to_list_string(str, l, deliminator, ignore_repeat)
    COOP_UNKNOWN_STRING str
    COOP_UNKNOWN_STRING,optional::deliminator
    logical,optional::ignore_repeat
    class(coop_list_string) l
    COOP_INT istart, idel, slen
    call coop_list_initialize(l)
    slen = len(str)
    if(slen .eq. 0)then
       return
    endif
    istart = 1
100 continue
    if(present(deliminator))then
       idel = scan(str(istart:slen), deliminator)
    else
       idel = scan(str(istart:slen), " "//coop_newline//coop_tab//coop_carriage_return//coop_backspace)
    endif
    select case(idel)
    case(0)
       call coop_list_push(l, str(istart:slen))
       return
    case(1)
       if(present(ignore_repeat))then
          if(.not.(ignore_repeat))then
             call coop_list_push(l, "")
          endif
       endif
       istart = istart + 1
       if(istart .gt. slen) return
       goto 100
    case default
       call coop_list_push(l, str(istart:istart+idel-2))
       istart = istart + idel
       if(istart .gt. slen) return
       goto 100
    end select
  end subroutine coop_string_to_list_string

  subroutine coop_string_to_list_real(str, l)
    class(coop_list_real) l
    COOP_UNKNOWN_STRING str
    type(coop_list_string) lstr
    COOP_STRING tmp
    COOP_INT i
    real(sp) x
    call coop_list_initialize(l)
    call coop_string_to_list_string(str, lstr, " ,;"//coop_newline//coop_tab//coop_backspace//coop_carriage_return)
    do i=1, lstr%n
       call coop_list_get_element(lstr, i, tmp)
       read(tmp, *) x
       call coop_list_push(l, x)
    enddo
    call coop_list_initialize(lstr)
  end subroutine coop_string_to_list_real



  subroutine coop_string_to_list_double(str, l)
    class(coop_list_double) l
    COOP_UNKNOWN_STRING str
    type(coop_list_string) lstr
    COOP_STRING tmp
    COOP_INT i
    real(dl) x
    call coop_list_initialize(l)
    call coop_string_to_list_string(str, lstr, " ,;"//coop_newline//coop_tab//coop_backspace//coop_carriage_return)
    do i=1, lstr%n
       call coop_list_get_element(lstr, i, tmp)
       read(tmp, *) x
       call coop_list_push(l, x)
    enddo
    call coop_list_initialize(lstr)
  end subroutine coop_string_to_list_double



  subroutine coop_string_to_list_logical(str, l)
    class(coop_list_logical) l
    COOP_UNKNOWN_STRING str
    type(coop_list_string) lstr
    COOP_STRING tmp
    COOP_INT i
    logical x
    call coop_list_initialize(l)
    call coop_string_to_list_string(str, lstr, " ,;"//coop_newline//coop_tab//coop_backspace//coop_carriage_return)
    do i=1, lstr%n
       call coop_list_get_element(lstr, i, tmp)
       read(tmp, *) x
       call coop_list_push(l, x)
    enddo
    call coop_list_initialize(lstr)
  end subroutine coop_string_to_list_logical



  subroutine coop_string_to_list_integer(str, l)
    class(coop_list_integer) l
    COOP_UNKNOWN_STRING str
    type(coop_list_string) lstr
    COOP_STRING tmp
    COOP_INT i
    COOP_INT x
    call coop_list_initialize(l)
    call coop_string_to_list_string(str, lstr, " ,;"//coop_newline//coop_tab//coop_backspace//coop_carriage_return)
    do i=1, lstr%n
       call coop_list_get_element(lstr, i, tmp)
       read(tmp, *) x
       call coop_list_push(l, x)
    enddo
    call coop_list_initialize(lstr)
  end subroutine coop_string_to_list_integer

!! element
  function coop_list_integer_element(l, i) result(elem)
    class(coop_list_integer) l
    COOP_INT i, j
    COOP_INT elem
    call coop_list_integer_get_element(l, i, elem)
  end function coop_list_integer_element

  function coop_list_real_element(l, i) result(elem)
    class(coop_list_real) l
    COOP_INT i, j
    real(sp) elem
    call coop_list_real_get_element(l, i, elem)
  end function coop_list_real_element


  function coop_list_double_element(l, i) result(elem)
    class(coop_list_double) l
    COOP_INT i, j
    real(dl) elem
    call coop_list_double_get_element(l, i, elem)
  end function coop_list_double_element


  function coop_list_logical_element(l, i) result(elem)
    class(coop_list_logical) l
    COOP_INT i, j
    logical elem
    call coop_list_logical_get_element(l, i, elem)
  end function coop_list_logical_element


  function coop_list_character_element(l, i) result(elem)
    class(coop_list_character) l
    COOP_INT i, j
    character elem
    call coop_list_character_get_element(l, i, elem)
  end function coop_list_character_element


  function coop_list_string_element(l, i) result(elem)
    class(coop_list_string) l
    COOP_INT i, j
    COOP_STRING elem
    call coop_list_string_get_element(l, i, elem)
  end function coop_list_string_element

  function coop_list_realarr_element(l, i) result(elem)
    class(coop_list_realarr) l
    COOP_INT i, j
    real(sp) elem(l%dim)
    call coop_list_realarr_get_element(l, i, elem)
  end function coop_list_realarr_element



!! get_element
  subroutine coop_list_integer_get_element(l, i, elem)
    class(coop_list_integer) l
    COOP_INT i, j
    COOP_INT elem
    if(i.le. coop_list_i1_max_length)then
       elem = l%i1(i)
       return
    endif
    j = i - coop_list_i1_max_length
    if(j .le. coop_list_i2_max_length)then
       elem = l%i2(j)
       return
    endif
    j = j - coop_list_i2_max_length
    if(j.le. coop_list_i3_max_length)then
       elem = l%i3(j)
       return
    endif
    elem = l%i4(j - coop_list_i3_max_length)
    return
  end subroutine coop_list_integer_get_element

  subroutine coop_list_real_get_element(l, i, elem)
    class(coop_list_real) l
    COOP_INT i, j
    real(sp) elem
    if(i.le. coop_list_i1_max_length)then
       elem = l%i1(i)
       return
    endif
    j = i - coop_list_i1_max_length
    if(j .le. coop_list_i2_max_length)then
       elem = l%i2(j)
       return
    endif
    j = j - coop_list_i2_max_length
    if(j.le. coop_list_i3_max_length)then
       elem = l%i3(j)
       return
    endif
    elem = l%i4(j - coop_list_i3_max_length)
    return
  end subroutine coop_list_real_get_element


  subroutine coop_list_double_get_element(l, i, elem)
    class(coop_list_double) l
    COOP_INT i, j
    real(dl) elem
    if(i.le. coop_list_i1_max_length)then
       elem = l%i1(i)
       return
    endif
    j = i - coop_list_i1_max_length
    if(j .le. coop_list_i2_max_length)then
       elem = l%i2(j)
       return
    endif
    j = j - coop_list_i2_max_length
    if(j.le. coop_list_i3_max_length)then
       elem = l%i3(j)
       return
    endif
    elem = l%i4(j - coop_list_i3_max_length)
    return
  end subroutine coop_list_double_get_element


  subroutine coop_list_logical_get_element(l, i, elem)
    class(coop_list_logical) l
    COOP_INT i, j
    logical elem
    if(i.le. coop_list_i1_max_length)then
       elem = l%i1(i)
       return
    endif
    j = i - coop_list_i1_max_length
    if(j .le. coop_list_i2_max_length)then
       elem = l%i2(j)
       return
    endif
    j = j - coop_list_i2_max_length
    if(j.le. coop_list_i3_max_length)then
       elem = l%i3(j)
       return
    endif
    elem = l%i4(j - coop_list_i3_max_length)
    return
  end subroutine coop_list_logical_get_element


  subroutine coop_list_character_get_element(l, i, elem)
    class(coop_list_character) l
    COOP_INT i, j
    character elem
    if(i.le. coop_list_i1_max_length)then
       elem = l%i1(i)
       return
    endif
    j = i - coop_list_i1_max_length
    if(j .le. coop_list_i2_max_length)then
       elem = l%i2(j)
       return
    endif
    j = j - coop_list_i2_max_length
    if(j.le. coop_list_i3_max_length)then
       elem = l%i3(j)
       return
    endif
    elem = l%i4(j - coop_list_i3_max_length)
    return
  end subroutine coop_list_character_get_element


  subroutine coop_list_string_get_element(l, i, elem)
    class(coop_list_string) l
    COOP_INT i, j
    COOP_STRING elem
    if(i.le. coop_list_s1_max_length)then
       elem = l%i1(i)
       return
    endif
    j = i - coop_list_s1_max_length
    if(j .le. coop_list_i2_max_length)then
       elem = l%i2(j)
       return
    endif
    j = j - coop_list_s2_max_length
    if(j.le. coop_list_i3_max_length)then
       elem = l%i3(j)
       return
    endif
    elem = l%i4(j - coop_list_s3_max_length)
    return
  end subroutine coop_list_string_get_element

  subroutine coop_list_realarr_get_element(l, i, elem)
    class(coop_list_realarr) l
    COOP_INT i, j
    real(sp) elem(:)
    if(i.le. coop_list_i1_max_length)then
       elem = l%i1(:,i)
       return
    endif
    j = i - coop_list_i1_max_length
    if(j .le. coop_list_i2_max_length)then
       elem = l%i2(:,j)
       return
    endif
    j = j - coop_list_i2_max_length
    if(j.le. coop_list_i3_max_length)then
       elem = l%i3(:,j)
       return
    endif
    elem = l%i4(:,j - coop_list_i3_max_length)
    return
  end subroutine coop_list_realarr_get_element

  subroutine coop_dictionary_insert(dict, key, val, overwrite)
    class(coop_dictionary):: dict
    COOP_UNKNOWN_STRING key, val
    logical,optional::overwrite
    COOP_SHORT_STRING,dimension(:),allocatable::tmpkey
    COOP_STRING,dimension(:),allocatable::tmpval
    COOP_INT,dimension(:),allocatable::tmpid
    COOP_INT iup, ilow, imid
    if(trim(adjustl(key)).eq."") return
    if(.not.allocated(dict%key))then
       dict%capacity = coop_list_unit_len
       allocate(dict%key(dict%capacity), dict%val(dict%capacity), dict%id(dict%capacity))
       dict%n = 1
       dict%key(1) = trim(adjustl(key))
       dict%val(1) = trim(adjustl(val))
       dict%id(1) = 1
       return
    endif
    if(dict%n .eq. dict%capacity)then
       allocate(tmpkey(dict%n), tmpval(dict%n), tmpid(dict%n))
       tmpkey = dict%key
       tmpval = dict%val
       tmpid = dict%id
       deallocate(dict%key, dict%val, dict%id)
       dict%capacity = dict%capacity + coop_list_unit_len
       allocate(dict%key(dict%capacity), dict%val(dict%capacity), dict%id(dict%capacity))
       dict%key(1:dict%n) = tmpkey
       dict%val(1:dict%n) = tmpval
       dict%id(1:dict%n) = tmpid
       deallocate(tmpkey, tmpval, tmpid)
    endif
    iup = dict%n
    ilow = 1
    dict%n = dict%n + 1
    dict%key(dict%n) = trim(adjustl(key))
    dict%val(dict%n) = trim(adjustl(val))
    if(lgt(dict%key(dict%n), dict%key(dict%id(iup))))then
       dict%id(dict%n) = dict%n
       return
    endif
    if(llt(dict%key(dict%n), dict%key(dict%id(ilow))))then
       dict%id(2:dict%n) = dict%id(1:dict%n-1)
       dict%id(1) = dict%n
       return
    endif
    do while(iup .gt. ilow+1)
       imid = (iup + ilow)/2
       if(lle(dict%key(dict%n), dict%key(dict%id(imid))))then
          iup = imid
       else
          ilow = imid
       endif
    end do
    if(dict%key(dict%id(ilow)) .eq. dict%key(dict%n))then
       if(present(overwrite))then
          if(.not. overwrite)then
             write(*,*) "key conflict, cannot insert into the dictionary."
             return
          endif
       endif
       dict%val(dict%id(ilow)) = dict%val(dict%n)
       dict%n = dict%n - 1
    elseif(dict%key(dict%id(iup)) .eq. dict%key(dict%n))then
       if(present(overwrite))then
          if(.not. overwrite)then
             write(*,*) "key conflict, cannot insert into the dictionary."
             return
          endif
       endif
       dict%val(dict%id(iup)) = dict%val(dict%n)
       dict%n = dict%n - 1
    else
       dict%id(iup+1:dict%n) = dict%id(iup:dict%n-1)
       dict%id(iup) = dict%n
    endif
  end subroutine coop_dictionary_insert



  function coop_dictionary_key_index(dict, key) result(ind)
    class(coop_dictionary):: dict
    COOP_UNKNOWN_STRING, intent(IN)::key
    COOP_SHORT_STRING thekey
    COOP_INT ind
    COOP_INT iup, ilow, imid
    if(.not. allocated(dict%key))then
       ind = 0
       return
    endif
    iup = dict%n
    ilow = 1
    thekey = trim(adjustl(key))
    if(lgt(thekey, dict%key(dict%id(iup))))then
       ind = 0
       return
    endif
    if(llt(thekey, dict%key(dict%id(ilow))))then
       ind = 0
       return
    endif
    do while(iup .gt. ilow+1)
       imid = (iup + ilow)/2
       if(lle(thekey, dict%key(dict%id(imid))))then
          iup = imid
       else
          ilow = imid
       endif
    end do
    if(thekey .eq. dict%key(dict%id(ilow)))then
       ind = dict%id(ilow)
       return
    endif
    if(thekey .eq. dict%key(dict%id(iup)))then
       ind = dict%id(iup)
       return
    endif
    ind = 0
    return
  end function coop_dictionary_key_index


  subroutine coop_dictionary_update(dict, key, val)
    class(coop_dictionary)::dict
    COOP_UNKNOWN_STRING,intent(IN)::key
    COOP_UNKNOWN_STRING,intent(IN)::val
    call dict%insert(key, val, overwrite = .true.)
  end subroutine coop_dictionary_update


  subroutine coop_dictionary_lookup_string(dict, key, val)
    class(coop_dictionary):: dict
    COOP_UNKNOWN_STRING, intent(IN):: key
    COOP_UNKNOWN_STRING, intent(OUT)::val
    COOP_INT ind
    ind = coop_dictionary_key_index(dict, key)
    if(ind .eq. 0)then
       val = ""
    else
       val = trim(dict%val(ind))
    endif
    return
  end subroutine coop_dictionary_lookup_string

  subroutine coop_dictionary_lookup_int(dict, key, val, default_val)
    class(coop_dictionary):: dict
    COOP_UNKNOWN_STRING key
    COOP_INT ind
    COOP_INT val
    COOP_INT,optional::default_val
    ind = coop_dictionary_key_index(dict, key)
    if(ind .eq. 0)then
       if(present(default_val))then
          val = default_val
       else
          stop "coop_dictionary_lookup_int: key not found; you need to set a default value in this case"
       endif
    else
       read(dict%val(ind), *) val
    endif
  end subroutine coop_dictionary_lookup_int


  subroutine coop_dictionary_lookup_real(dict, key, val, default_val)
    class(coop_dictionary):: dict
    COOP_UNKNOWN_STRING key
    COOP_INT ind
    real(sp) val
    real(sp),optional::default_val
    ind = coop_dictionary_key_index(dict, key)
    if(ind .eq. 0)then
       if(present(default_val))then
          val = default_val
       else
          stop "coop_dictionary_lookup_real: key not found; you need to set a default value in this case"
       endif
    else
       read(dict%val(ind), *) val
    endif
  end subroutine coop_dictionary_lookup_real


  subroutine coop_dictionary_lookup_double(dict, key, val, default_val)
    class(coop_dictionary):: dict
    COOP_UNKNOWN_STRING key
    COOP_INT ind
    real(dl) val
    real(dl),optional::default_val
    ind = coop_dictionary_key_index(dict, key)
    if(ind .eq. 0)then
       if(present(default_val))then
          val = default_val
       else
          stop "coop_dictionary_lookup_double: key not found; you need to set a default value in this case"
       endif
    else
       read(dict%val(ind), *) val
    endif
  end subroutine coop_dictionary_lookup_double


 
 
  subroutine coop_dictionary_lookup_logical(dict, key, val, default_val)
    class(coop_dictionary):: dict
    COOP_UNKNOWN_STRING key
    COOP_INT ind
    logical val
    logical,optional::default_val
    ind = coop_dictionary_key_index(dict, key)
    if(ind .eq. 0)then
       if(present(default_val))then
          val = default_val
       else
          stop "coop_dictionary_lookup_logical: key not found; you need to set a default value in this case"
       endif
    else
       read(dict%val(ind), *) val
    endif
  end subroutine coop_dictionary_lookup_logical



  subroutine coop_dictionary_lookup_int_array(dict, key, val, default_val)
    class(coop_dictionary):: dict
    COOP_UNKNOWN_STRING key
    COOP_INT ind
    COOP_INT,dimension(:),intent(OUT):: val
    COOP_INT,dimension(:),optional::default_val
    ind = coop_dictionary_key_index(dict, key)
    if(ind .eq. 0)then
       if(present(default_val))then
          val = default_val
       else
          stop "coop_dictionary_lookup_int: key not found; you need to set a default value in this case"
       endif
    else
       read(dict%val(ind), *) val
    endif
  end subroutine coop_dictionary_lookup_int_array


  subroutine coop_dictionary_lookup_real_array(dict, key, val, default_val)
    class(coop_dictionary):: dict
    COOP_UNKNOWN_STRING key
    COOP_INT ind
    real(sp),dimension(:),intent(OUT)::  val
    real(sp),dimension(:),optional::default_val
    ind = coop_dictionary_key_index(dict, key)
    if(ind .eq. 0)then
       if(present(default_val))then
          val = default_val
       else
          stop "coop_dictionary_lookup_real: key not found; you need to set a default value in this case"
       endif
    else
       read(dict%val(ind), *) val
    endif
  end subroutine coop_dictionary_lookup_real_array


  subroutine coop_dictionary_lookup_double_array(dict, key, val, default_val)
    class(coop_dictionary):: dict
    COOP_UNKNOWN_STRING key
    COOP_INT ind
    real(dl),dimension(:),intent(OUT)::  val
    real(dl),dimension(:),optional::default_val
    ind = coop_dictionary_key_index(dict, key)
    if(ind .eq. 0)then
       if(present(default_val))then
          val = default_val
       else
          stop "coop_dictionary_lookup_double: key not found; you need to set a default value in this case"
       endif
    else
       read(dict%val(ind), *) val
    endif
  end subroutine coop_dictionary_lookup_double_array
 
  subroutine coop_dictionary_lookup_logical_array(dict, key, val, default_val)
    class(coop_dictionary):: dict
    COOP_UNKNOWN_STRING key
    COOP_INT ind
    logical,dimension(:),intent(OUT)::  val
    logical,dimension(:),optional::default_val
    ind = coop_dictionary_key_index(dict, key)
    if(ind .eq. 0)then
       if(present(default_val))then
          val = default_val
       else
          stop "coop_dictionary_lookup_logical: key not found; you need to set a default value in this case"
       endif
    else
       read(dict%val(ind), *) val
    endif
  end subroutine coop_dictionary_lookup_logical_array


  function coop_dictionary_value(dict, key) result(val)
    class(coop_dictionary):: dict
    COOP_UNKNOWN_STRING key
    COOP_STRING val
    call coop_dictionary_lookup_string(dict, key, val)
  end function coop_dictionary_value


  subroutine coop_dictionary_free(dict)
    class(coop_dictionary):: dict
    if(allocated(dict%key))then
       deallocate(dict%key, dict%val, dict%id)
    endif
  end subroutine coop_dictionary_free

  subroutine coop_dictionary_print(dict)
    class(coop_dictionary):: dict
    COOP_INT i
    if(allocated(dict%key))then
       do i=1, dict%n
          write(*,"(3A)")  trim(dict%key(dict%id(i))), " = ", trim(dict%val(dict%id(i)))
       enddo
    else
       write(*,*) "Empty coop_dictionary"
    endif
  end subroutine coop_dictionary_print


  !!return the at most m largest prime number that is <= n; p(1)<=p(2)<=...
  subroutine coop_get_prime_numbers(n, pl)
    COOP_INT n
    type(coop_list_integer)::pl
    COOP_INT,dimension(:),allocatable::plist
    COOP_INT i, j, np, top, npmax
    logical isp
    npmax = min(n-1, ceiling((n+1.d0)/log(n+1.d0)*1.2) + 10000)
    allocate(plist(npmax))
    np = 0
    do i=2, n
       top = floor(sqrt(i+1.d-10))
       isp = .true.
       do j=1, np
          if(plist(j) .gt. top) exit
          if(mod(i, plist(j)).eq.0)then
             isp = .false.
             exit
          endif
       enddo
       if(isp)then
          np = np+1
          if(np.gt. npmax) stop "Prime numbers overflow"
          plist(np) = i
          call pl%push(i)
       endif
    enddo
    deallocate(plist)
  end subroutine coop_get_prime_numbers



End module coop_list_mod

