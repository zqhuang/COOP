module coop_list_mod
  use coop_wrapper_typedef
  implicit none

#include "constants.h"
  private

  public::coop_list_integer, coop_list_real, coop_list_realarr, coop_list_double, coop_list_logical, coop_list_string, coop_list_character, coop_string_to_list, coop_dictionary, coop_dictionary_lookup, coop_get_prime_numbers, coop_list_get_element,  coop_command_line_has_argument, coop_command_line_argument, coop_get_command_line_argument, coop_dynamic_array_integer, coop_int_table, coop_real_table, coop_int_set

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
     module procedure coop_dictionary_lookup_string, coop_dictionary_lookup_int, coop_dictionary_lookup_real, coop_dictionary_lookup_double, coop_dictionary_lookup_logical, coop_dictionary_lookup_int_array, coop_dictionary_lookup_real_array, coop_dictionary_lookup_double_array, coop_dictionary_lookup_logical_array, coop_dictionary_lookup_list_real, coop_dictionary_lookup_list_integer, coop_dictionary_lookup_list_double, coop_dictionary_lookup_list_logical, coop_dictionary_lookup_list_string
  end interface coop_dictionary_lookup

  interface coop_get_command_line_argument
     module procedure coop_get_command_line_argument_single, coop_get_command_line_argument_real, coop_get_command_line_argument_logical, coop_get_command_line_argument_int, coop_get_command_line_argument_str, coop_get_command_line_argument_intarr, coop_get_command_line_argument_realarr, coop_get_command_line_argument_singlearr
  end interface coop_get_command_line_argument


  COOP_INT,parameter:: coop_list_i1_max_length = 2**14
  COOP_INT,parameter:: coop_list_i2_max_length = coop_list_i1_max_length * 4
  COOP_INT,parameter:: coop_list_i3_max_length = coop_list_i2_max_length * 4
  COOP_INT,parameter:: coop_list_i4_max_length = coop_list_i3_max_length * 4


  COOP_INT,parameter:: coop_list_s1_max_length =  coop_list_i1_max_length/coop_string_length
  COOP_INT,parameter:: coop_list_s2_max_length = coop_list_s1_max_length * 4
  COOP_INT,parameter:: coop_list_s3_max_length = coop_list_s2_max_length * 4
  COOP_INT,parameter:: coop_list_s4_max_length = coop_list_s3_max_length * 4

  COOP_INT,parameter::coop_list_unit_len = 8192


  type coop_dynamic_array_integer
     COOP_INT::n = 0
     COOP_INT,dimension(:),allocatable::i
   contains
     procedure::init => coop_dynamic_array_integer_initialize
     procedure::free => coop_dynamic_array_integer_initialize
     procedure::get_indices => coop_dynamic_array_integer_get_indices
  end type coop_dynamic_array_integer


  type coop_list_integer
     COOP_INT::n =0
     COOP_INT::stack = 1
     COOP_INT::loc = 0
     COOP_INT,dimension(:),allocatable::i1
     COOP_INT,dimension(:),allocatable::i2
     COOP_INT,dimension(:),allocatable::i3
     COOP_INT,dimension(:),allocatable::i4
   contains
     procedure::init => coop_list_integer_initialize
     procedure::free => coop_list_integer_initialize     
     procedure::isinit => coop_list_integer_is_initialized
     procedure::push => coop_list_integer_push
     procedure::pop => coop_list_integer_pop
     procedure::element => coop_list_integer_element
     procedure::get_element => coop_list_integer_get_element
  end type coop_list_integer

  type coop_list_real
     COOP_INT::n =0
     COOP_INT::stack = 1
     COOP_INT::loc = 0
     COOP_SINGLE,dimension(:),allocatable::i1
     COOP_SINGLE,dimension(:),allocatable::i2
     COOP_SINGLE,dimension(:),allocatable::i3
     COOP_SINGLE,dimension(:),allocatable::i4
   contains
     procedure::init => coop_list_real_initialize
     procedure::free => coop_list_real_initialize     
     procedure::isinit => coop_list_real_is_initialized
     procedure::push => coop_list_real_push
     procedure::pop => coop_list_real_pop
     procedure::element => coop_list_real_element
     procedure::get_element => coop_list_real_get_element
  end type coop_list_real

  type coop_list_realarr
     COOP_INT dim
     COOP_INT::n = 0
     COOP_INT::stack = 1
     COOP_INT::loc = 0
     COOP_SINGLE,dimension(:,:),allocatable::i1
     COOP_SINGLE,dimension(:,:),allocatable::i2
     COOP_SINGLE,dimension(:,:),allocatable::i3
     COOP_SINGLE,dimension(:,:),allocatable::i4
   contains
     procedure::init => coop_list_realarr_initialize     
     procedure::free => coop_list_realarr_initialize     
     procedure::isinit => coop_list_realarr_is_initialized
     procedure::push => coop_list_realarr_push
     procedure::pop => coop_list_realarr_pop
     procedure::element => coop_list_realarr_element
     procedure::get_element => coop_list_realarr_get_element
     procedure::sort => coop_list_realarr_sort
  end type coop_list_realarr

  type coop_list_double
     COOP_INT::n = 0
     COOP_INT::stack = 1
     COOP_INT::loc = 0
     COOP_REAL,dimension(:),allocatable::i1
     COOP_REAL,dimension(:),allocatable::i2
     COOP_REAL,dimension(:),allocatable::i3
     COOP_REAL,dimension(:),allocatable::i4
   contains
     procedure::isinit => coop_list_double_is_initialized
     procedure::init => coop_list_double_initialize
     procedure::free => coop_list_double_initialize     
     procedure::push => coop_list_double_push
     procedure::pop => coop_list_double_pop
     procedure::element => coop_list_double_element
     procedure::get_element => coop_list_double_get_element
  end type coop_list_double

  type coop_list_logical
     COOP_INT::n = 0
     COOP_INT::stack=1
     COOP_INT::loc = 0
     logical,dimension(:),allocatable::i1
     logical,dimension(:),allocatable::i2
     logical,dimension(:),allocatable::i3
     logical,dimension(:),allocatable::i4
   contains
     procedure::isinit => coop_list_logical_is_initialized
     procedure::init => coop_list_logical_initialize
     procedure::free => coop_list_logical_initialize     
     procedure::push => coop_list_logical_push
     procedure::pop => coop_list_logical_pop
     procedure::element => coop_list_logical_element
     procedure::get_element => coop_list_logical_get_element
  end type coop_list_logical

  type coop_list_character
     COOP_INT::n = 0
     COOP_INT::stack = 1
     COOP_INT::loc = 0
     character,dimension(:),allocatable::i1
     character,dimension(:),allocatable::i2
     character,dimension(:),allocatable::i3
     character,dimension(:),allocatable::i4
   contains
     procedure::isinit => coop_list_character_is_initialized
     procedure::init => coop_list_character_initialize
     procedure::free => coop_list_character_initialize     
     procedure::push => coop_list_character_push
     procedure::pop => coop_list_character_pop
     procedure::element => coop_list_character_element
     procedure::get_element => coop_list_character_get_element
  end type coop_list_character

  type coop_list_string
     COOP_INT::n = 0
     COOP_INT::stack = 1
     COOP_INT::loc = 0
     COOP_STRING,dimension(:),allocatable::i1
     COOP_STRING,dimension(:),allocatable::i2
     COOP_STRING,dimension(:),allocatable::i3
     COOP_STRING,dimension(:),allocatable::i4
   contains
     procedure::isinit => coop_list_string_is_initialized
     procedure::init => coop_list_string_initialize
     procedure::free => coop_list_string_initialize     
     procedure::push => coop_list_string_push
     procedure::pop => coop_list_string_pop
     procedure::element => coop_list_string_element
     procedure::get_element => coop_list_string_get_element
  end type coop_list_string


  type coop_int_set
     COOP_INT::n = 0
     COOP_INT::capacity = 0
     COOP_INT,dimension(:),allocatable::x
   contains
     procedure::free => coop_int_set_free
     procedure::has => coop_int_set_has
     procedure::insert => coop_int_set_insert
     procedure::delete => coop_int_set_delete
  end type coop_int_set


  type coop_dictionary
     COOP_INT::n = 0
     COOP_INT::capacity = 0
     COOP_SHORT_STRING,dimension(:),allocatable::key
     COOP_STRING,dimension(:),allocatable::val
     COOP_INT, dimension(:),allocatable::id
   contains
     procedure::print => coop_dictionary_print
     procedure::insert => coop_dictionary_insert !!insert a key
     procedure::delete => coop_dictionary_delete  !!delete a key
     procedure::index => coop_dictionary_key_index
     procedure::value => coop_dictionary_value
     procedure::update => coop_dictionary_update
     procedure::free => coop_dictionary_free
     procedure::init => coop_dictionary_free
     procedure::load_input => coop_dictionary_load_input
  end type coop_dictionary


  type coop_int_table
     COOP_INT::n = 0
     COOP_INT::capacity = 0
     COOP_SHORT_STRING,dimension(:),allocatable::key
     COOP_INT,dimension(:),allocatable::val
     COOP_INT, dimension(:),allocatable::id
   contains
     procedure::print => coop_int_table_print
     procedure::load_dictionary => coop_int_table_load_dictionary
     procedure::insert => coop_int_table_insert !!insert a key
     procedure::delete => coop_int_table_delete  !!delete a key
     procedure::index => coop_int_table_key_index
     procedure::value => coop_int_table_value
     procedure::lookup => coop_int_table_lookup
     procedure::update => coop_int_table_update
     procedure::free => coop_int_table_free
     procedure::init => coop_int_table_free
  end type coop_int_table


  type coop_real_table
     COOP_INT::n = 0
     COOP_INT::capacity = 0
     COOP_SHORT_STRING,dimension(:),allocatable::key
     COOP_REAL,dimension(:),allocatable::val
     COOP_INT, dimension(:),allocatable::id
   contains
     procedure::print => coop_real_table_print
     procedure::load_dictionary => coop_real_table_load_dictionary
     procedure::insert => coop_real_table_insert !!insert a key
     procedure::delete => coop_real_table_delete  !!delete a key
     procedure::index => coop_real_table_key_index
     procedure::value => coop_real_table_value
     procedure::update => coop_real_table_update
     procedure::lookup => coop_real_table_lookup
     procedure::free => coop_real_table_free
     procedure::init => coop_real_table_free
  end type coop_real_table



  type(coop_dictionary)::coop_command_line_inputs

contains

!! initialize
  subroutine coop_dynamic_array_integer_initialize(this)
    class(coop_dynamic_array_integer)::this
    COOP_DEALLOC(this%i)
    this%n = 0
  end subroutine coop_dynamic_array_integer_initialize

  subroutine coop_dynamic_array_integer_get_indices(this, mask, start_index)
    class(coop_dynamic_array_integer)::this
    logical,dimension(:),intent(IN)::mask
    COOP_INT::n, i, j
    COOP_INT,optional::start_index
    n = size(mask)
    call this%free()
    this%n = count(mask)
    allocate(this%i(this%n))
    j = 0
    do i=1, n
       if(mask(i))then
          j = j + 1
          this%i(j) = i
       endif
    enddo
    if(present(start_index))then
       this%i = this%i + (start_index - 1)
    endif
  end subroutine coop_dynamic_array_integer_get_indices
  
  subroutine coop_list_string_initialize(l)
    class(coop_list_string) l
    l%n = 0
    l%stack = 1
    l%loc = 0
    COOP_DEALLOC(l%i4)
    COOP_DEALLOC(l%i3)
    COOP_DEALLOC(l%i2)
    COOP_DEALLOC(l%i1)
  end subroutine coop_list_string_initialize



  subroutine coop_list_real_initialize(l)
    class(coop_list_real) l
    l%n = 0
    l%stack = 1
    l%loc = 0
    COOP_DEALLOC(l%i4)
    COOP_DEALLOC(l%i3)
    COOP_DEALLOC(l%i2)
    COOP_DEALLOC(l%i1)
  end subroutine coop_list_real_initialize



  subroutine coop_list_realarr_initialize(l)
    class(coop_list_realarr) l
    l%n = 0
    l%stack = 1
    l%loc = 0
    COOP_DEALLOC(l%i4)
    COOP_DEALLOC(l%i3)
    COOP_DEALLOC(l%i2)
    COOP_DEALLOC(l%i1)
  end subroutine coop_list_realarr_initialize


  subroutine coop_list_double_initialize(l)
    class(coop_list_double) l
    l%n = 0
    l%stack = 1
    l%loc = 0
    COOP_DEALLOC(l%i4)
    COOP_DEALLOC(l%i3)
    COOP_DEALLOC(l%i2)
    COOP_DEALLOC(l%i1)
  end subroutine coop_list_double_initialize


  subroutine coop_list_integer_initialize(l)
    class(coop_list_integer) l
    l%n = 0
    l%stack = 1
    l%loc = 0
    COOP_DEALLOC(l%i4)
    COOP_DEALLOC(l%i3)
    COOP_DEALLOC(l%i2)
    COOP_DEALLOC(l%i1)
  end subroutine coop_list_integer_initialize


  subroutine coop_list_logical_initialize(l)
    class(coop_list_logical) l
    l%n = 0
    l%stack = 1
    l%loc = 0
    COOP_DEALLOC(l%i4)
    COOP_DEALLOC(l%i3)
    COOP_DEALLOC(l%i2)
    COOP_DEALLOC(l%i1)
  end subroutine coop_list_logical_initialize


  subroutine coop_list_character_initialize(l)
    class(coop_list_character) l
    l%n = 0
    l%stack = 1
    l%loc = 0
    COOP_DEALLOC(l%i4)
    COOP_DEALLOC(l%i3)
    COOP_DEALLOC(l%i2)
    COOP_DEALLOC(l%i1)
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
    COOP_SINGLE i
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
    COOP_SINGLE,dimension(:),intent(IN)::i
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
    COOP_REAL i
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
    COOP_SINGLE x
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
    COOP_REAL x
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
    COOP_SINGLE elem
    call coop_list_real_get_element(l, i, elem)
  end function coop_list_real_element


  function coop_list_double_element(l, i) result(elem)
    class(coop_list_double) l
    COOP_INT i, j
    COOP_REAL elem
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
    COOP_SINGLE elem(l%dim)
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
    COOP_SINGLE elem
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
    COOP_REAL elem
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
    COOP_SINGLE elem(:)
    if(i.le. coop_list_i1_max_length)then
       elem(1:l%dim) = l%i1(:,i)
       return
    endif
    j = i - coop_list_i1_max_length
    if(j .le. coop_list_i2_max_length)then
       elem(1:l%dim) = l%i2(:,j)
       return
    endif
    j = j - coop_list_i2_max_length
    if(j.le. coop_list_i3_max_length)then
       elem(1:l%dim) = l%i3(:,j)
       return
    endif
    elem(1:l%dim) = l%i4(:,j - coop_list_i3_max_length)
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
    if(dict%n .eq. 0)then
       dict%n = 1
       dict%key(1) = trim(adjustl(key))
       dict%val(1) = trim(adjustl(val))
       dict%id(1) = 1
       return
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
!             write(*,*) "key "//trim(dict%key(dict%n))//" conflicts with existing keys, cannot insert into the dictionary."
             dict%n = dict%n - 1
             return
          endif
       endif
       dict%val(dict%id(ilow)) = dict%val(dict%n)
       dict%n = dict%n - 1
    elseif(dict%key(dict%id(iup)) .eq. dict%key(dict%n))then
       if(present(overwrite))then
          if(.not. overwrite)then
!             write(*,*) "key "//trim(dict%key(dict%n))//" conflicts with existing keys, cannot insert into the dictionary."
             dict%n = dict%n - 1
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





  subroutine coop_dictionary_delete(dict, key)
    class(coop_dictionary):: dict
    COOP_UNKNOWN_STRING key
    COOP_SHORT_STRING,dimension(:),allocatable::tmpkey
    COOP_INT ind, i, j
    ind = dict%index(key)
    if(ind .eq. 0) return
    dict%key(ind:dict%n-1) = dict%key(ind+1:dict%n)
    dict%val(ind:dict%n-1) = dict%val(ind+1:dict%n)
    do i=1, dict%n
       if(dict%id(i).gt.ind)then
          dict%id(i) = dict%id(i)-1
       elseif(dict%id(i) .eq. ind)then
          j = i
       endif
    enddo
    dict%id(j:dict%n-1) = dict%id(j+1:dict%n)
    dict%n = dict%n-1
    if(dict%n .eq. 0) call dict%free()
  end subroutine coop_dictionary_delete
  


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


  subroutine coop_dictionary_lookup_string(dict, key, val, default_val)
    class(coop_dictionary):: dict
    COOP_UNKNOWN_STRING, intent(IN):: key
    COOP_UNKNOWN_STRING, intent(OUT)::val
    COOP_UNKNOWN_STRING, optional::default_val
    COOP_INT ind
    ind = coop_dictionary_key_index(dict, key)
    if(ind .eq. 0)then
       if(present(default_val))then
          val = trim(default_val)
       else
          val = ""
       endif
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
          write(*,*) "coop_dictionary_lookup: key "//trim(key)//" not found; you need to set a default value in this case"
          stop 
       endif
    else
       read(dict%val(ind), *) val
    endif
  end subroutine coop_dictionary_lookup_int


  subroutine coop_dictionary_lookup_real(dict, key, val, default_val)
    class(coop_dictionary):: dict
    COOP_UNKNOWN_STRING key
    COOP_INT ind
    COOP_SINGLE val
    COOP_SINGLE,optional::default_val
    ind = coop_dictionary_key_index(dict, key)
    if(ind .eq. 0)then
       if(present(default_val))then
          val = default_val
       else
          write(*,*) "coop_dictionary_lookup: key "//trim(key)//" not found; you need to set a default value in this case"
          stop 
       endif
    else
       read(dict%val(ind), *) val
    endif
  end subroutine coop_dictionary_lookup_real


  subroutine coop_dictionary_lookup_double(dict, key, val, default_val)
    class(coop_dictionary):: dict
    COOP_UNKNOWN_STRING key
    COOP_INT ind
    COOP_REAL val
    COOP_REAL,optional::default_val
    ind = coop_dictionary_key_index(dict, key)
    if(ind .eq. 0)then
       if(present(default_val))then
          val = default_val
       else
          write(*,*) "coop_dictionary_lookup: key "//trim(key)//" not found; you need to set a default value in this case"
          stop 
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
          write(*,*) "coop_dictionary_lookup: key "//trim(key)//" not found; you need to set a default value in this case"
          stop 
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
          write(*,*) "coop_dictionary_lookup: key "//trim(key)//" not found; you need to set a default value in this case"
          stop 
       endif
    else
       read(dict%val(ind), *) val
    endif
  end subroutine coop_dictionary_lookup_int_array


  subroutine coop_dictionary_lookup_real_array(dict, key, val, default_val)
    class(coop_dictionary):: dict
    COOP_UNKNOWN_STRING key
    COOP_INT ind
    COOP_SINGLE,dimension(:),intent(OUT)::  val
    COOP_SINGLE,dimension(:),optional::default_val
    ind = coop_dictionary_key_index(dict, key)
    if(ind .eq. 0)then
       if(present(default_val))then
          val = default_val
       else
          write(*,*) "coop_dictionary_lookup: key "//trim(key)//" not found; you need to set a default value in this case"
          stop 
       endif
    else
       read(dict%val(ind), *) val
    endif
  end subroutine coop_dictionary_lookup_real_array


  subroutine coop_dictionary_lookup_double_array(dict, key, val, default_val)
    class(coop_dictionary):: dict
    COOP_UNKNOWN_STRING key
    COOP_INT ind
    COOP_REAL,dimension(:),intent(OUT)::  val
    COOP_REAL,dimension(:),optional::default_val
    ind = coop_dictionary_key_index(dict, key)
    if(ind .eq. 0)then
       if(present(default_val))then
          val = default_val
       else
          write(*,*) "coop_dictionary_lookup: key "//trim(key)//" not found; you need to set a default value in this case"
          stop 
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
          write(*,*) "coop_dictionary_lookup: key "//trim(key)//" not found; you need to set a default value in this case"
          stop 
       endif
    else
       read(dict%val(ind), *) val
    endif
  end subroutine coop_dictionary_lookup_logical_array


  subroutine coop_dictionary_lookup_list_real(dict, key, l)
    class(coop_dictionary)::dict
    COOP_UNKNOWN_STRING::key
    type(coop_list_real)::l
    COOP_STRING::val
    call coop_dictionary_lookup_string(dict, key, val)
    call coop_string_to_list(val, l)
  end subroutine coop_dictionary_lookup_list_real

  subroutine coop_dictionary_lookup_list_double(dict, key, l)
    class(coop_dictionary)::dict
    COOP_UNKNOWN_STRING::key
    type(coop_list_double)::l
    COOP_STRING::val
    call coop_dictionary_lookup_string(dict, key, val)
    call coop_string_to_list(val, l)
  end subroutine coop_dictionary_lookup_list_double

  subroutine coop_dictionary_lookup_list_logical(dict, key, l)
    class(coop_dictionary)::dict
    COOP_UNKNOWN_STRING::key
    type(coop_list_logical)::l
    COOP_STRING::val    
    call coop_dictionary_lookup_string(dict, key, val)
    call coop_string_to_list(val, l)
  end subroutine coop_dictionary_lookup_list_logical

  subroutine coop_dictionary_lookup_list_string(dict, key, l)
    class(coop_dictionary)::dict
    COOP_UNKNOWN_STRING::key
    type(coop_list_string)::l
    COOP_STRING::val    
    call coop_dictionary_lookup_string(dict, key, val)
    call coop_string_to_list(val, l)
  end subroutine coop_dictionary_lookup_list_string
  

  subroutine coop_dictionary_lookup_list_integer(dict, key, l)
    class(coop_dictionary)::dict
    COOP_UNKNOWN_STRING::key
    type(coop_list_integer)::l
    COOP_STRING::val    
    call coop_dictionary_lookup_string(dict, key, val)
    call coop_string_to_list(val, l)    
  end subroutine coop_dictionary_lookup_list_integer
  
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



  subroutine coop_list_realarr_sort(this, key_index)
    COOP_INT key_index
    class(coop_list_realarr) this
    COOP_SINGLE, dimension(:),allocatable::key, tmp
    type(coop_list_realarr) copy
    COOP_INT, dimension(:),allocatable::index
    COOP_INT i, j
    allocate(key(this%n), index(this%n), tmp(this%dim))
    do i=1, this%n
       call this%get_element(i, tmp)
       key(i) = tmp(key_index)
    enddo
    call coop_quicksort_index(key, index)
    select type(this)
    type is(coop_list_realarr)
       copy = this
       call this%init()
       do i = 1, copy%n
          call this%push(copy%element(index(i)))
       enddo
    end select
    deallocate(key, index, tmp)
  end subroutine coop_list_realarr_sort

  subroutine coop_dictionary_load_input(this)
    class(coop_dictionary)::this
    COOP_INT::i, n
    COOP_STRING::key, val
    n = iargc()
    i= 1
    do while(i .le. n)
       call getarg(i, key)
       key = trim(adjustl(key))
       if(key(1:1) .eq. '-')then
          if(i.eq.n)then
             val = ''
          else
             i = i + 1
             call getarg(i, val)
          endif
          call this%insert(key = key(2:), val = val, overwrite = .true.)
       else
          val = key
          key = COOP_STR_OF(i)
          call this%insert(key = key, val = val, overwrite = .true.)          
       endif
       i = i +1
    enddo
  end subroutine coop_dictionary_load_input

  subroutine coop_load_command_line_inputs()
    if(coop_command_line_inputs%n .eq. 0) &
         call coop_command_line_inputs%load_input()
  end subroutine coop_load_command_line_inputs  

  function coop_command_line_has_argument(key) result(has)
    COOP_UNKNOWN_STRING::key
    logical has
    call coop_load_command_line_inputs()
    has = (coop_command_line_inputs%index(key) .ne. 0)
  end function coop_command_line_has_argument

  function coop_command_line_argument(index, key, default) result(str)
    COOP_INT,optional::index
    COOP_UNKNOWN_STRING,optional::key, default
    COOP_STRING::str
    if(present(index))then
       str = coop_InputArgs(index)
    elseif(present(key))then
       call coop_load_command_line_inputs()
       str = coop_command_line_inputs%value(trim(adjustl(key)))
    else
       if(present(default))then
          str = default
       else
          str = coop_InputArgs(1)
       endif
    endif
  end function coop_command_line_argument

  subroutine coop_get_command_line_argument_str(index, key, arg, default)
    COOP_INT,optional::index
    COOP_UNKNOWN_STRING,optional::key, default
    COOP_UNKNOWN_STRING::arg
    COOP_INT::ind
    if(present(index))then
       if(present(default))then
          call coop_get_Input(index, arg, default)
       else
          call coop_get_Input(index, arg)
       endif
    elseif(present(key))then
       call coop_load_command_line_inputs()
       ind = coop_command_line_inputs%index(trim(adjustl(key)))
       if(ind.eq.0)then
          if(present(default))then
             arg = default
             return
          endif
          write(*,*) 'key '//trim(key)//' missing in command line'
          stop
       endif
       arg = trim(coop_command_line_inputs%val(ind))
    else
       if(present(default))then
          arg = default
       else
          stop "cannot find the string argument in command line input"
       endif
    endif
  end subroutine coop_get_command_line_argument_str



  subroutine coop_get_command_line_argument_int(index, key, arg, default)
    COOP_INT,optional::index
    COOP_UNKNOWN_STRING,optional::key
    COOP_INT,optional::default
    COOP_INT::arg
    COOP_STRING::str
    if(present(index))then
       if(present(default))then
          call coop_get_Input(index, arg, default)
       else
          call coop_get_Input(index, arg)
       endif
    elseif(present(key))then
       call coop_load_command_line_inputs()
       str = coop_command_line_inputs%value(trim(adjustl(key)))
       if(trim(str).eq.'')then
          if(present(default))then
             arg = default
             return
          endif
          write(*,*) 'key '//trim(key)//' missing in command line'
          stop
       endif       
       read(str, *) arg
    else
       if(present(default))then
          arg = default
       else
          stop "cannot find the integer argument in command line input"
       endif
    endif
  end subroutine coop_get_command_line_argument_int



  subroutine coop_get_command_line_argument_intarr(index, key, arg, default)
    COOP_INT,optional::index
    COOP_UNKNOWN_STRING,optional::key
    COOP_INT,optional::default(:)
    COOP_INT::arg(:)
    COOP_STRING::str
    if(present(index))then       
       if(present(default))then
          call coop_get_Input(index, str, "")
       else
          call coop_get_Input(index, str)
       endif
    elseif(present(key))then
       call coop_load_command_line_inputs()
       str = coop_command_line_inputs%value(trim(adjustl(key)))
       if(trim(str).eq.'')then
          if(present(default))then
             arg = default
             return
          endif
          write(*,*) 'key '//trim(key)//' missing in command line'
          stop
       endif       
    else
       if(present(default))then
          arg = default
          return
       else
          stop "cannot find the integer argument in command line input"
       endif       
    endif
    read(str, *) arg
  end subroutine coop_get_command_line_argument_intarr

  
  subroutine coop_get_command_line_argument_logical(index, key, arg, default)
    COOP_INT,optional::index
    COOP_UNKNOWN_STRING,optional::key
    logical::arg
    logical,optional::default    
    COOP_STRING::str
    if(present(index))then
       if(present(default))then
          call coop_get_Input(index, arg, default)
       else
          call coop_get_Input(index, arg)
       endif
    elseif(present(key))then
       call coop_load_command_line_inputs()
       str = coop_command_line_inputs%value(trim(adjustl(key)))
       if(trim(str).eq.'')then
          if(present(default))then
             arg = default
             return
          endif          
          write(*,*) 'key '//trim(key)//' missing in command line'
          stop
       endif
       read(str, *) arg
    else
       if(present(default))then
          arg = default
       else
          stop "cannot find the logical argument in command line input"
       endif
    endif
  end subroutine coop_get_command_line_argument_logical



  subroutine coop_get_command_line_argument_real(index, key, arg, default)
    COOP_INT,optional::index
    COOP_UNKNOWN_STRING,optional::key
    COOP_REAL::arg
    COOP_REAL,optional::default    
    COOP_STRING::str
    if(present(index))then
       if(present(default))then
          call coop_get_Input(index, arg, default)
       else
          call coop_get_Input(index, arg)
       endif
    elseif(present(key))then
       call coop_load_command_line_inputs()
       str = coop_command_line_inputs%value(trim(adjustl(key)))
       if(trim(str).eq.'')then
          if(present(default))then
             arg = default
             return
          endif          
          write(*,*) 'key '//trim(key)//' missing in command line'
          stop
       endif       
       read(str, *) arg
    else
       if(present(default))then
          arg = default
       else
          stop "cannot find the REAL argument in command line input"
       endif
    endif
  end subroutine coop_get_command_line_argument_real




  subroutine coop_get_command_line_argument_realarr(index, key, arg, default)
    COOP_INT,optional::index
    COOP_UNKNOWN_STRING,optional::key
    COOP_REAL,optional::default(:)
    COOP_REAL::arg(:)
    COOP_STRING::str
    if(present(index))then       
       if(present(default))then
          call coop_get_Input(index, str, "")
       else
          call coop_get_Input(index, str)
       endif
    elseif(present(key))then
       call coop_load_command_line_inputs()
       str = coop_command_line_inputs%value(trim(adjustl(key)))
       if(trim(str).eq.'')then
          if(present(default))then
             arg = default
             return
          endif
          write(*,*) 'key '//trim(key)//' missing in command line'
          stop
       endif       
    else
       if(present(default))then
          arg = default
          return
       else
          stop "cannot find the integer argument in command line input"
       endif       
    endif
    read(str, *) arg
  end subroutine coop_get_command_line_argument_realarr
  


  subroutine coop_get_command_line_argument_single(index, key, arg, default)
    COOP_INT,optional::index
    COOP_UNKNOWN_STRING,optional::key
    COOP_SINGLE::arg
    COOP_SINGLE,optional::default    
    COOP_STRING::str
    if(present(index))then
       if(present(default))then
          call coop_get_Input(index, arg, default)
       else
          call coop_get_Input(index, arg)
       endif
    elseif(present(key))then
       call coop_load_command_line_inputs()
       str = coop_command_line_inputs%value(trim(adjustl(key)))
       if(trim(str).eq.'')then
          if(present(default))then
             arg = default
             return
          endif          
          write(*,*) 'key '//trim(key)//' missing in command line'
          stop
       endif
       read(str, *) arg
    else
       if(present(default))then
          arg = default
       else
          stop "cannot find the SINGLE argument in command line input"
       endif
    endif
  end subroutine coop_get_command_line_argument_single



  subroutine coop_get_command_line_argument_singlearr(index, key, arg, default)
    COOP_INT,optional::index
    COOP_UNKNOWN_STRING,optional::key
    COOP_SINGLE,optional::default(:)
    COOP_SINGLE::arg(:)
    COOP_STRING::str
    if(present(index))then       
       if(present(default))then
          call coop_get_Input(index, str, "")
       else
          call coop_get_Input(index, str)
       endif
    elseif(present(key))then
       call coop_load_command_line_inputs()
       str = coop_command_line_inputs%value(trim(adjustl(key)))
       if(trim(str).eq.'')then
          if(present(default))then
             arg = default
             return
          endif
          write(*,*) 'key '//trim(key)//' missing in command line'
          stop
       endif       
    else
       if(present(default))then
          arg = default
          return
       else
          stop "cannot find the integer argument in command line input"
       endif       
    endif
    read(str, *) arg
  end subroutine coop_get_command_line_argument_singlearr
  

  subroutine coop_int_table_load_dictionary(this, dict)
    class(coop_int_table)::this
    type(coop_dictionary)::dict
    COOP_INT::i
    do i = 1, dict%n
       if(coop_is_integer(dict%val(i)))then
          call this%insert(dict%key(i), coop_str2int(dict%val(i)))
       endif
    enddo
  end subroutine coop_int_table_load_dictionary
  

  subroutine coop_int_table_insert(this, key, val, overwrite)
    class(coop_int_table):: this
    COOP_UNKNOWN_STRING key
    COOP_INT::val
    logical,optional::overwrite
    COOP_SHORT_STRING,dimension(:),allocatable::tmpkey
    COOP_INT,dimension(:),allocatable::tmpval
    COOP_INT,dimension(:),allocatable::tmpid
    COOP_INT iup, ilow, imid
    if(trim(adjustl(key)).eq."") return
    if(.not.allocated(this%key))then
       this%capacity = coop_list_unit_len
       allocate(this%key(this%capacity), this%val(this%capacity), this%id(this%capacity))
       this%n = 1
       this%key(1) = trim(adjustl(key))
       this%val(1) = val
       this%id(1) = 1
       return
    endif
    if(this%n .eq. this%capacity)then
       allocate(tmpkey(this%n), tmpval(this%n), tmpid(this%n))
       tmpkey = this%key
       tmpval = this%val
       tmpid = this%id
       deallocate(this%key, this%val, this%id)
       this%capacity = this%capacity + coop_list_unit_len
       allocate(this%key(this%capacity), this%val(this%capacity), this%id(this%capacity))
       this%key(1:this%n) = tmpkey
       this%val(1:this%n) = tmpval
       this%id(1:this%n) = tmpid
       deallocate(tmpkey, tmpval, tmpid)
    endif
    if(this%n .eq. 0)then
       this%n = 1
       this%key(1) = trim(adjustl(key))
       this%val(1) = val
       this%id(1) = 1
       return
    endif
    iup = this%n
    ilow = 1
    this%n = this%n + 1
    this%key(this%n) = trim(adjustl(key))
    this%val(this%n) = val
    if(lgt(this%key(this%n), this%key(this%id(iup))))then
       this%id(this%n) = this%n
       return
    endif
    if(llt(this%key(this%n), this%key(this%id(ilow))))then
       this%id(2:this%n) = this%id(1:this%n-1)
       this%id(1) = this%n
       return
    endif
    do while(iup .gt. ilow+1)
       imid = (iup + ilow)/2
       if(lle(this%key(this%n), this%key(this%id(imid))))then
          iup = imid
       else
          ilow = imid
       endif
    end do
    if(this%key(this%id(ilow)) .eq. this%key(this%n))then
       if(present(overwrite))then
          if(.not. overwrite)then
!             write(*,*) "key "//trim(this%key(this%n))//" conflicts with existing keys, cannot insert into the int_table."
             this%n = this%n - 1
             return
          endif
       endif
       this%val(this%id(ilow)) = this%val(this%n)
       this%n = this%n - 1
    elseif(this%key(this%id(iup)) .eq. this%key(this%n))then
       if(present(overwrite))then
          if(.not. overwrite)then
!             write(*,*) "key "//trim(this%key(this%n))//" conflicts with existing keys, cannot insert into the int_table."
             this%n = this%n - 1
             return
          endif
       endif
       this%val(this%id(iup)) = this%val(this%n)
       this%n = this%n - 1
    else
       this%id(iup+1:this%n) = this%id(iup:this%n-1)
       this%id(iup) = this%n
    endif
  end subroutine coop_int_table_insert


  subroutine coop_int_table_delete(this, key)
    class(coop_int_table):: this
    COOP_UNKNOWN_STRING key
    COOP_SHORT_STRING,dimension(:),allocatable::tmpkey
    COOP_INT ind, i, j
    ind = this%index(key)
    if(ind .eq. 0) return
    this%key(ind:this%n-1) = this%key(ind+1:this%n)
    this%val(ind:this%n-1) = this%val(ind+1:this%n)
    do i=1, this%n
       if(this%id(i).gt.ind)then
          this%id(i) = this%id(i)-1
       elseif(this%id(i) .eq. ind)then
          j = i
       endif
    enddo
    this%id(j:this%n-1) = this%id(j+1:this%n)
    this%n = this%n-1
    if(this%n .eq. 0) call this%free()
  end subroutine coop_int_table_delete
  


  function coop_int_table_key_index(this, key) result(ind)
    class(coop_int_table):: this
    COOP_UNKNOWN_STRING, intent(IN)::key
    COOP_SHORT_STRING thekey
    COOP_INT ind
    COOP_INT iup, ilow, imid
    if(.not. allocated(this%key))then
       ind = 0
       return
    endif
    iup = this%n
    ilow = 1
    thekey = trim(adjustl(key))
    if(lgt(thekey, this%key(this%id(iup))))then
       ind = 0
       return
    endif
    if(llt(thekey, this%key(this%id(ilow))))then
       ind = 0
       return
    endif
    do while(iup .gt. ilow+1)
       imid = (iup + ilow)/2
       if(lle(thekey, this%key(this%id(imid))))then
          iup = imid
       else
          ilow = imid
       endif
    end do
    if(thekey .eq. this%key(this%id(ilow)))then
       ind = this%id(ilow)
       return
    endif
    if(thekey .eq. this%key(this%id(iup)))then
       ind = this%id(iup)
       return
    endif
    ind = 0
    return
  end function coop_int_table_key_index


  subroutine coop_int_table_update(this, key, val)
    class(coop_int_table)::this
    COOP_UNKNOWN_STRING,intent(IN)::key
    COOP_INT,intent(IN)::val
    call this%insert(key, val, overwrite = .true.)
  end subroutine coop_int_table_update


  subroutine coop_int_table_lookup(this, key, val, default_val)
    class(coop_int_table):: this
    COOP_UNKNOWN_STRING, intent(IN):: key
    COOP_INT, intent(OUT)::val
    COOP_INT, optional::default_val
    COOP_INT ind
    ind = coop_int_table_key_index(this, key)
    if(ind .eq. 0)then
       if(present(default_val))then
          val = default_val
       else
          write(*,*) "int_table_lookup: key "//trim(key)//" cannot be found."
          call this%print()
          stop
       endif
    else
       val = this%val(ind)
    endif
    return
  end subroutine coop_int_table_lookup

  subroutine coop_int_table_print(this)
    class(coop_int_table):: this
    COOP_INT i
    if(allocated(this%key))then
       do i=1, this%n
          write(*,"(2A,I10)")  trim(this%key(this%id(i))), " = ", this%val(this%id(i))
       enddo
    else
       write(*,*) "Empty coop_int_table"
    endif
  end subroutine coop_int_table_print


  function coop_int_table_value(this, key) result(val)
    class(coop_int_table):: this
    COOP_UNKNOWN_STRING:: key
    COOP_INT val
    call coop_int_table_lookup(this, key, val)
  end function coop_int_table_value


  subroutine coop_int_table_free(this)
    class(coop_int_table):: this
    if(allocated(this%key))then
       deallocate(this%key, this%val, this%id)
    endif
  end subroutine coop_int_table_free



  subroutine coop_real_table_insert(this, key, val, overwrite)
    class(coop_real_table):: this
    COOP_UNKNOWN_STRING key
    COOP_REAL::val
    logical,optional::overwrite
    COOP_SHORT_STRING,dimension(:),allocatable::tmpkey
    COOP_REAL,dimension(:),allocatable::tmpval
    COOP_INT,dimension(:),allocatable::tmpid
    COOP_INT iup, ilow, imid
    if(trim(adjustl(key)).eq."") return
    if(.not.allocated(this%key))then
       this%capacity = coop_list_unit_len
       allocate(this%key(this%capacity), this%val(this%capacity), this%id(this%capacity))
       this%n = 1
       this%key(1) = trim(adjustl(key))
       this%val(1) = val
       this%id(1) = 1
       return
    endif
    if(this%n .eq. this%capacity)then
       allocate(tmpkey(this%n), tmpval(this%n), tmpid(this%n))
       tmpkey = this%key
       tmpval = this%val
       tmpid = this%id
       deallocate(this%key, this%val, this%id)
       this%capacity = this%capacity + coop_list_unit_len
       allocate(this%key(this%capacity), this%val(this%capacity), this%id(this%capacity))
       this%key(1:this%n) = tmpkey
       this%val(1:this%n) = tmpval
       this%id(1:this%n) = tmpid
       deallocate(tmpkey, tmpval, tmpid)
    endif
    if(this%n .eq. 0)then
       this%n = 1
       this%key(1) = trim(adjustl(key))
       this%val(1) = val
       this%id(1) = 1
       return
    endif
    iup = this%n
    ilow = 1
    this%n = this%n + 1
    this%key(this%n) = trim(adjustl(key))
    this%val(this%n) = val
    if(lgt(this%key(this%n), this%key(this%id(iup))))then
       this%id(this%n) = this%n
       return
    endif
    if(llt(this%key(this%n), this%key(this%id(ilow))))then
       this%id(2:this%n) = this%id(1:this%n-1)
       this%id(1) = this%n
       return
    endif
    do while(iup .gt. ilow+1)
       imid = (iup + ilow)/2
       if(lle(this%key(this%n), this%key(this%id(imid))))then
          iup = imid
       else
          ilow = imid
       endif
    end do
    if(this%key(this%id(ilow)) .eq. this%key(this%n))then
       if(present(overwrite))then
          if(.not. overwrite)then
!             write(*,*) "key "//trim(this%key(this%n))//" conflicts with existing keys, cannot insert into the real_table."
             this%n = this%n - 1
             return
          endif
       endif
       this%val(this%id(ilow)) = this%val(this%n)
       this%n = this%n - 1
    elseif(this%key(this%id(iup)) .eq. this%key(this%n))then
       if(present(overwrite))then
          if(.not. overwrite)then
!             write(*,*) "key "//trim(this%key(this%n))//" conflicts with existing keys, cannot insert into the real_table."
             this%n = this%n - 1
             return
          endif
       endif
       this%val(this%id(iup)) = this%val(this%n)
       this%n = this%n - 1
    else
       this%id(iup+1:this%n) = this%id(iup:this%n-1)
       this%id(iup) = this%n
    endif
  end subroutine coop_real_table_insert


  subroutine coop_real_table_delete(this, key)
    class(coop_real_table):: this
    COOP_UNKNOWN_STRING key
    COOP_SHORT_STRING,dimension(:),allocatable::tmpkey
    COOP_INT ind, i, j
    ind = this%index(key)
    if(ind .eq. 0) return
    this%key(ind:this%n-1) = this%key(ind+1:this%n)
    this%val(ind:this%n-1) = this%val(ind+1:this%n)
    do i=1, this%n
       if(this%id(i).gt.ind)then
          this%id(i) = this%id(i)-1
       elseif(this%id(i) .eq. ind)then
          j = i
       endif
    enddo
    this%id(j:this%n-1) = this%id(j+1:this%n)
    this%n = this%n-1
    if(this%n .eq. 0) call this%free()
  end subroutine coop_real_table_delete
  


  function coop_real_table_key_index(this, key) result(ind)
    class(coop_real_table):: this
    COOP_UNKNOWN_STRING, intent(IN)::key
    COOP_SHORT_STRING thekey
    COOP_INT ind
    COOP_INT iup, ilow, imid
    if(.not. allocated(this%key))then
       ind = 0
       return
    endif
    iup = this%n
    ilow = 1
    thekey = trim(adjustl(key))
    if(lgt(thekey, this%key(this%id(iup))))then
       ind = 0
       return
    endif
    if(llt(thekey, this%key(this%id(ilow))))then
       ind = 0
       return
    endif
    do while(iup .gt. ilow+1)
       imid = (iup + ilow)/2
       if(lle(thekey, this%key(this%id(imid))))then
          iup = imid
       else
          ilow = imid
       endif
    end do
    if(thekey .eq. this%key(this%id(ilow)))then
       ind = this%id(ilow)
       return
    endif
    if(thekey .eq. this%key(this%id(iup)))then
       ind = this%id(iup)
       return
    endif
    ind = 0
    return
  end function coop_real_table_key_index


  subroutine coop_real_table_update(this, key, val)
    class(coop_real_table)::this
    COOP_UNKNOWN_STRING,intent(IN)::key
    COOP_REAL,intent(IN)::val
    call this%insert(key, val, overwrite = .true.)
  end subroutine coop_real_table_update


  subroutine coop_real_table_lookup(this, key, val, default_val)
    class(coop_real_table):: this
    COOP_UNKNOWN_STRING, intent(IN):: key
    COOP_REAL, intent(OUT)::val
    COOP_REAL, optional::default_val
    COOP_INT ind
    ind = coop_real_table_key_index(this, key)
    if(ind .eq. 0)then
       if(present(default_val))then
          val = default_val
       else
          write(*,*) "int_table_lookup: key "//trim(key)//" cannot be found."
          call this%print()
          stop
       endif
    else
       val = this%val(ind)
    endif
    return
  end subroutine coop_real_table_lookup

  subroutine coop_real_table_print(this)
    class(coop_real_table):: this
    COOP_INT i
    if(allocated(this%key))then
       do i=1, this%n
          write(*,"(2A,G14.5)")  trim(this%key(this%id(i))), " = ", this%val(this%id(i))
       enddo
    else
       write(*,*) "Empty coop_real_table"
    endif
  end subroutine coop_real_table_print


  function coop_real_table_value(this, key) result(val)
    class(coop_real_table):: this
    COOP_UNKNOWN_STRING:: key
    COOP_REAL val
    call coop_real_table_lookup(this, key, val)
  end function coop_real_table_value


  subroutine coop_real_table_free(this)
    class(coop_real_table):: this
    if(allocated(this%key))then
       deallocate(this%key, this%val, this%id)
    endif
  end subroutine coop_real_table_free

  subroutine coop_real_table_load_dictionary(this, dict)
    class(coop_real_table)::this
    type(coop_dictionary)::dict
    COOP_INT::i
    do i = 1, dict%n
       if(coop_is_number(dict%val(i)))then
          call this%insert(dict%key(i), coop_str2real(dict%val(i)))
       endif
    enddo
  end subroutine coop_real_table_load_dictionary


  subroutine coop_int_set_free(this)
    class(coop_int_set):: this
    COOP_DEALLOC(this%x)
    this%n = 0
    this%capacity = 0
  end subroutine coop_int_set_free

  subroutine coop_int_set_insert(this, x)
    class(coop_int_set):: this
    COOP_INT::x
    COOP_INT, dimension(:),allocatable::tmpx
    COOP_INT::iup, ilow, imid
    if(.not. allocated(this%x))then
       this%capacity = coop_list_unit_len
       allocate(this%x(this%capacity))
       this%n = 1
       this%x(1) = x
       return
    endif
    if(this%n .eq. this%capacity)then
       allocate(tmpx(this%n))
       tmpx = this%x
       deallocate(this%x)
       this%capacity = this%capacity + coop_list_unit_len
       allocate(this%x(this%capacity))
       this%x(1:this%n) = tmpx
       deallocate(tmpx)
    endif
    if(this%n .eq.  0)then
       this%x(1) = x
       this%n = 1
       return
    endif
    if( x .lt. this%x(1))then
       this%x(2:this%n+1) = this%x(1:this%n)
       this%x(1) = x
       this%n = this%n + 1
       return
    endif
    if(x .gt. this%x(this%n))then
       this%x(this%n+1) = x
       this%n = this%n + 1
       return
    endif
    iup = this%n
    ilow = 1
    do while(iup .gt. ilow+1)
       imid = (iup + ilow)/2
       if(this%x(imid) .gt. x)then
          iup = imid
       else
          ilow = imid
       endif
    end do
    if(this%x(ilow).eq. x)return
    if(this%x(iup).eq.x)return
    this%x(iup+1:this%n+1) = this%x(iup:this%n)
    this%x(iup) = x
    this%n  = this%n + 1
  end subroutine coop_int_set_insert


  subroutine coop_int_set_delete(this, x)
    class(coop_int_set):: this
    COOP_INT::x
    COOP_INT::iup, ilow, imid
    if(.not. allocated(this%x) .or. this%n .eq. 0) return
    if( x .lt. this%x(1) .or. x .gt. this%x(this%n))return
    if(x .eq. this%x(this%n))then
       this%n = this%n - 1
       if(this%n .eq. 0)call this%free()
       return
    endif
    iup = this%n  !!this%x(iup) strictly > x
    ilow = 1     !!this%x(ilow) <= x
    do while(iup .gt. ilow+1)
       imid = (iup + ilow)/2
       if(this%x(imid) .gt. x)then
          iup = imid
       else
          ilow = imid
       endif
    end do
    if(this%x(ilow) .eq. x)then
       this%x(ilow:this%n-1) = this%x(ilow+1:this%n)
       this%n = this%n - 1
    endif
  end subroutine coop_int_set_delete

  function coop_int_set_has(this, x) result(has)
    class(coop_int_set):: this
    COOP_INT::x
    COOP_INT::iup, ilow, imid
    logical::has
    has = .false.
    if(.not. allocated(this%x) .or. this%n .eq. 0) return
    if( x .lt. this%x(1) .or. x .gt. this%x(this%n))return
    if(x .eq. this%x(this%n))then
       has = .true.
       return
    endif
    iup = this%n  !!this%x(iup) strictly > x
    ilow = 1     !!this%x(ilow) <= x
    do while(iup .gt. ilow+1)
       imid = (iup + ilow)/2
       if(this%x(imid) .gt. x)then
          iup = imid
       else
          ilow = imid
       endif
    end do
    has = (this%x(ilow) .eq. x)
  end function coop_int_set_has

End module coop_list_mod


