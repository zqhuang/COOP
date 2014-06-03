module coop_list_mod
  use coop_wrapper_typedef
  implicit none
#include "constants.h"

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


  interface string_to_coop_list
     module procedure string_to_coop_list_string, string_to_coop_list_real, string_to_coop_list_double, string_to_coop_list_integer, string_to_coop_list_logical
  end interface string_to_coop_list

  interface coop_dictionary_lookup
     module procedure coop_dictionary_lookup_string, coop_dictionary_lookup_int, coop_dictionary_lookup_real, coop_dictionary_lookup_double, coop_dictionary_lookup_logical, coop_dictionary_lookup_int_array, coop_dictionary_lookup_real_array, coop_dictionary_lookup_double_array, coop_dictionary_lookup_logical_array
  end interface coop_dictionary_lookup


  integer,parameter:: coop_list_i1_max_length = 2**12
  integer,parameter:: coop_list_i2_max_length = coop_list_i1_max_length * 4
  integer,parameter:: coop_list_i3_max_length = coop_list_i2_max_length * 4
  integer,parameter:: coop_list_i4_max_length = coop_list_i3_max_length * 4


  integer,parameter:: coop_list_s1_max_length =  coop_list_i1_max_length/coop_string_length
  integer,parameter:: coop_list_s2_max_length = coop_list_s1_max_length * 4
  integer,parameter:: coop_list_s3_max_length = coop_list_s2_max_length * 4
  integer,parameter:: coop_list_s4_max_length = coop_list_s3_max_length * 4

  integer,parameter::coop_list_unit_len = 8192

  type coop_list_integer
     integer n, stack, loc
     integer,dimension(:),allocatable::i1
     integer,dimension(:),allocatable::i2
     integer,dimension(:),allocatable::i3
     integer,dimension(:),allocatable::i4
  end type coop_list_integer

  type coop_list_real
     integer n, stack, loc
     real,dimension(:),allocatable::i1
     real,dimension(:),allocatable::i2
     real,dimension(:),allocatable::i3
     real,dimension(:),allocatable::i4
  end type coop_list_real

  type coop_list_realarr
     integer dim
     integer n, stack, loc
     real,dimension(:,:),allocatable::i1
     real,dimension(:,:),allocatable::i2
     real,dimension(:,:),allocatable::i3
     real,dimension(:,:),allocatable::i4
  end type coop_list_realarr

  type coop_list_double
     integer n, stack, loc
     COOP_REAL,dimension(:),allocatable::i1
     COOP_REAL,dimension(:),allocatable::i2
     COOP_REAL,dimension(:),allocatable::i3
     COOP_REAL,dimension(:),allocatable::i4
  end type coop_list_double

  type coop_list_logical
     integer n, stack, loc
     logical,dimension(:),allocatable::i1
     logical,dimension(:),allocatable::i2
     logical,dimension(:),allocatable::i3
     logical,dimension(:),allocatable::i4
  end type coop_list_logical

  type coop_list_character
     integer n, stack, loc
     character,dimension(:),allocatable::i1
     character,dimension(:),allocatable::i2
     character,dimension(:),allocatable::i3
     character,dimension(:),allocatable::i4
  end type coop_list_character

  type coop_list_string
     integer n, stack, loc
     COOP_STRING,dimension(:),allocatable::i1
     COOP_STRING,dimension(:),allocatable::i2
     COOP_STRING,dimension(:),allocatable::i3
     COOP_STRING,dimension(:),allocatable::i4
  end type coop_list_string

  type coop_dictionary
     integer n, capacity
     COOP_SHORT_STRING,dimension(:),allocatable::key
     COOP_STRING,dimension(:),allocatable::val
     integer, dimension(:),allocatable::id
   contains
     procedure::print => coop_dictionary_print
     procedure::insert => coop_dictionary_insert
     procedure::index => coop_dictionary_key_index
     procedure::value => coop_dictionary_value
     procedure::free => coop_dictionary_free
  end type coop_dictionary

  type coop_point_3dfield
     real f
     integer x, y, z 
  end type coop_point_3dfield

  type coop_list_3dfield
     integer n
     integer capacity
     type(coop_point_3dfield),dimension(:),allocatable::p
  end type coop_list_3dfield



contains


!!3d filed
  subroutine init_coop_list_3dfield(l3d)
    type(coop_list_3dfield) l3d
    if(allocated(l3d%p))deallocate(l3d%p)
    allocate(l3d%p(coop_list_unit_len))
    l3d%n = 0
    l3d%capacity = coop_list_unit_len
  end subroutine init_coop_list_3dfield

  subroutine print_coop_list_3dfield(l3d)
    type(coop_list_3dfield) l3d
    integer i
    do i=1,l3d%n
       write(*,"(A,I5,A,I5,A,I5,A,E14.5)") "(",l3d%p(i)%x,",",l3d%p(i)%y,",",l3d%p(i)%z," ): ",l3d%p(i)%f
    enddo
  end subroutine print_coop_list_3dfield

  subroutine insert_coop_list_3dfield(l3d, f, x,y,z)
    type(coop_list_3dfield) l3d
    real f
    integer x, y, z
    type(coop_point_3dfield),dimension(:),allocatable:: ptmp
    if(l3d%n .ge. l3d%capacity)then
       allocate(ptmp(l3d%capacity))
       ptmp(1:l3d%capacity)=l3d%p(1:l3d%capacity)
       deallocate(l3d%p)
       allocate(l3d%p(l3d%capacity + coop_list_unit_len))
       l3d%p(1:l3d%capacity) = ptmp(1:l3d%capacity)
       l3d%capacity = l3d%capacity + coop_list_unit_len
       deallocate(ptmp)
    endif
    l3d%n=l3d%n+1
    l3d%p(l3d%n)%f = f
    l3d%p(l3d%n)%x = x
    l3d%p(l3d%n)%y = y
    l3d%p(l3d%n)%z = z
  end subroutine insert_coop_list_3dfield

  subroutine sort_coop_list_3dfield(l3d) !!f descending
    type(coop_list_3dfield)l3d
    integer,dimension(:,:),Allocatable::stack  !!StACK
    type(coop_point_3dfield) x
    integer l1,l2,lpoint,i,j
    Allocate(stack(2,l3d%n))
    l1=1
    l2=l3d%n
    lpoint=0
    do
       do While(l1.lt.l2)
          i=l1
          j=l2
          x=l3d%p(l1)
          do while(i.lt.j)
             do While(i.lt.j .and. x%f.ge. l3d%p(j)%f)
                j=j-1
             enddo
             if(i.lt.j)then
                l3d%p(i)=l3d%p(j)
                i=i+1
                do While(i.lt.j .And. l3d%p(i)%f .gt. X%f)
                   i=i+1
                enddo
                if(i.lt.j)then
                   l3d%p(j)=l3d%p(i)
                   j=j-1
                endif
             endif
          enddo
          l3d%p(i)=x
          lpoint=lpoint+1
          stack(1,lpoint)=i+1
          stack(2,lpoint)=l2
          l2=i-1
       enddo
       if(lpoint.eq.0)eXit
       l1=stack(1,lpoint)
       l2=stack(2,lpoint)
       lpoint=lpoint-1
    enddo
    deallocate(stack)
  end subroutine sort_coop_list_3dfield

!! initialize


  subroutine coop_list_string_initialize(l)
    type(coop_list_string) l
    l%n = 0
    l%stack = 1
    l%loc = 0
    if(allocated(l%i4))deallocate(l%i4)
    if(allocated(l%i3))deallocate(l%i3)
    if(allocated(l%i2))deallocate(l%i2)
    if(allocated(l%i1))deallocate(l%i1)
  end subroutine coop_list_string_initialize



  subroutine coop_list_real_initialize(l)
    type(coop_list_real) l
    l%n = 0
    l%stack = 1
    l%loc = 0
    if(allocated(l%i4))deallocate(l%i4)
    if(allocated(l%i3))deallocate(l%i3)
    if(allocated(l%i2))deallocate(l%i2)
    if(allocated(l%i1))deallocate(l%i1)
  end subroutine coop_list_real_initialize



  subroutine coop_list_realarr_initialize(l)
    type(coop_list_realarr) l
    l%n = 0
    l%stack = 1
    l%loc = 0
    if(allocated(l%i4))deallocate(l%i4)
    if(allocated(l%i3))deallocate(l%i3)
    if(allocated(l%i2))deallocate(l%i2)
    if(allocated(l%i1))deallocate(l%i1)
  end subroutine coop_list_realarr_initialize


  subroutine coop_list_double_initialize(l)
    type(coop_list_double) l
    l%n = 0
    l%stack = 1
    l%loc = 0
    if(allocated(l%i4))deallocate(l%i4)
    if(allocated(l%i3))deallocate(l%i3)
    if(allocated(l%i2))deallocate(l%i2)
    if(allocated(l%i1))deallocate(l%i1)
  end subroutine coop_list_double_initialize


  subroutine coop_list_integer_initialize(l)
    type(coop_list_integer) l
    l%n = 0
    l%stack = 1
    l%loc = 0
    if(allocated(l%i4))deallocate(l%i4)
    if(allocated(l%i3))deallocate(l%i3)
    if(allocated(l%i2))deallocate(l%i2)
    if(allocated(l%i1))deallocate(l%i1)
  end subroutine coop_list_integer_initialize


  subroutine coop_list_logical_initialize(l)
    type(coop_list_logical) l
    l%n = 0
    l%stack = 1
    l%loc = 0
    if(allocated(l%i4))deallocate(l%i4)
    if(allocated(l%i3))deallocate(l%i3)
    if(allocated(l%i2))deallocate(l%i2)
    if(allocated(l%i1))deallocate(l%i1)
  end subroutine coop_list_logical_initialize


  subroutine coop_list_character_initialize(l)
    type(coop_list_character) l
    l%n = 0
    l%stack = 1
    l%loc = 0
    if(allocated(l%i4))deallocate(l%i4)
    if(allocated(l%i3))deallocate(l%i3)
    if(allocated(l%i2))deallocate(l%i2)
    if(allocated(l%i1))deallocate(l%i1)
  end subroutine coop_list_character_initialize


  function coop_list_integer_is_initialized(l) result(ini)
    type(coop_list_integer) l
    logical ini
    ini = allocated(l%i1)
  end function coop_list_integer_is_initialized


  function coop_list_real_is_initialized(l) result(ini)
    type(coop_list_real) l
    logical ini
    ini = allocated(l%i1)
  end function coop_list_real_is_initialized


  function coop_list_double_is_initialized(l) result(ini)
    type(coop_list_double) l
    logical ini
    ini = allocated(l%i1)
  end function coop_list_double_is_initialized


  function coop_list_character_is_initialized(l) result(ini)
    type(coop_list_character) l
    logical ini
    ini = allocated(l%i1)
  end function coop_list_character_is_initialized

  function coop_list_logical_is_initialized(l) result(ini)
    type(coop_list_logical) l
    logical ini
    ini = allocated(l%i1)
  end function coop_list_logical_is_initialized

  function coop_list_string_is_initialized(l) result(ini)
    type(coop_list_string) l
    logical ini
    ini = allocated(l%i1)
  end function coop_list_string_is_initialized

  function coop_list_realarr_is_initialized(l) result(ini)
    type(coop_list_realarr) l
    logical ini
    ini = allocated(l%i1)
  end function coop_list_realarr_is_initialized


!!push and pop
  subroutine coop_list_integer_push(l, i)
    type(coop_list_integer) l
    integer i
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
    type(coop_list_integer) l
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
    type(coop_list_real) l
    real i
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
    type(coop_list_real) l
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
    type(coop_list_realarr) l
    real,dimension(:),intent(IN)::i
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
    type(coop_list_realarr) l
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
    type(coop_list_double) l
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
    type(coop_list_double) l
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
    type(coop_list_logical) l
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
    type(coop_list_logical) l
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
    type(coop_list_character) l
    COOP_UNKNOWN_STRING i
    integer slen, cut
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
    type(coop_list_character) l
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
    type(coop_list_string) l
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
    type(coop_list_string) l
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
  subroutine string_to_coop_list_string(str, l, deliminator, ignore_repeat)
    COOP_UNKNOWN_STRING str
    COOP_UNKNOWN_STRING,optional::deliminator
    logical,optional::ignore_repeat
    type(coop_list_string) l
    integer istart, idel, slen
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
  end subroutine string_to_coop_list_string

  subroutine string_to_coop_list_real(str, l)
    type(coop_list_real) l
    COOP_UNKNOWN_STRING str
    type(coop_list_string) lstr
    COOP_STRING tmp
    integer i
    real x
    call coop_list_initialize(l)
    call string_to_coop_list_string(str, lstr, " ,;"//coop_newline//coop_tab//coop_backspace//coop_carriage_return)
    do i=1, lstr%n
       call coop_list_get_element(lstr, i, tmp)
       read(tmp, *) x
       call coop_list_push(l, x)
    enddo
    call coop_list_initialize(lstr)
  end subroutine string_to_coop_list_real



  subroutine string_to_coop_list_double(str, l)
    type(coop_list_double) l
    COOP_UNKNOWN_STRING str
    type(coop_list_string) lstr
    COOP_STRING tmp
    integer i
    COOP_REAL x
    call coop_list_initialize(l)
    call string_to_coop_list_string(str, lstr, " ,;"//coop_newline//coop_tab//coop_backspace//coop_carriage_return)
    do i=1, lstr%n
       call coop_list_get_element(lstr, i, tmp)
       read(tmp, *) x
       call coop_list_push(l, x)
    enddo
    call coop_list_initialize(lstr)
  end subroutine string_to_coop_list_double



  subroutine string_to_coop_list_logical(str, l)
    type(coop_list_logical) l
    COOP_UNKNOWN_STRING str
    type(coop_list_string) lstr
    COOP_STRING tmp
    integer i
    logical x
    call coop_list_initialize(l)
    call string_to_coop_list_string(str, lstr, " ,;"//coop_newline//coop_tab//coop_backspace//coop_carriage_return)
    do i=1, lstr%n
       call coop_list_get_element(lstr, i, tmp)
       read(tmp, *) x
       call coop_list_push(l, x)
    enddo
    call coop_list_initialize(lstr)
  end subroutine string_to_coop_list_logical



  subroutine string_to_coop_list_integer(str, l)
    type(coop_list_integer) l
    COOP_UNKNOWN_STRING str
    type(coop_list_string) lstr
    COOP_STRING tmp
    integer i
    integer x
    call coop_list_initialize(l)
    call string_to_coop_list_string(str, lstr, " ,;"//coop_newline//coop_tab//coop_backspace//coop_carriage_return)
    do i=1, lstr%n
       call coop_list_get_element(lstr, i, tmp)
       read(tmp, *) x
       call coop_list_push(l, x)
    enddo
    call coop_list_initialize(lstr)
  end subroutine string_to_coop_list_integer

!! element
  function coop_list_integer_element(l, i) result(elem)
    type(coop_list_integer) l
    integer i, j
    integer elem
    call coop_list_integer_get_element(l, i, elem)
  end function coop_list_integer_element

  function coop_list_real_element(l, i) result(elem)
    type(coop_list_real) l
    integer i, j
    real elem
    call coop_list_real_get_element(l, i, elem)
  end function coop_list_real_element


  function coop_list_double_element(l, i) result(elem)
    type(coop_list_double) l
    integer i, j
    COOP_REAL elem
    call coop_list_double_get_element(l, i, elem)
  end function coop_list_double_element


  function coop_list_logical_element(l, i) result(elem)
    type(coop_list_logical) l
    integer i, j
    logical elem
    call coop_list_logical_get_element(l, i, elem)
  end function coop_list_logical_element


  function coop_list_character_element(l, i) result(elem)
    type(coop_list_character) l
    integer i, j
    character elem
    call coop_list_character_get_element(l, i, elem)
  end function coop_list_character_element


  function coop_list_string_element(l, i) result(elem)
    type(coop_list_string) l
    integer i, j
    COOP_STRING elem
    call coop_list_string_get_element(l, i, elem)
  end function coop_list_string_element

  function coop_list_realarr_element(l, i) result(elem)
    type(coop_list_realarr) l
    integer i, j
    real elem(l%dim)
    call coop_list_realarr_get_element(l, i, elem)
  end function coop_list_realarr_element



!! get_element
  subroutine coop_list_integer_get_element(l, i, elem)
    type(coop_list_integer) l
    integer i, j
    integer elem
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
    type(coop_list_real) l
    integer i, j
    real elem
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
    type(coop_list_double) l
    integer i, j
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
    type(coop_list_logical) l
    integer i, j
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
    type(coop_list_character) l
    integer i, j
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
    type(coop_list_string) l
    integer i, j
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
    type(coop_list_realarr) l
    integer i, j
    real elem(:)
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

  subroutine coop_dictionary_insert(dict, key, val)
    class(coop_dictionary):: dict
    COOP_SHORT_STRING,dimension(:),allocatable::tmpkey
    COOP_STRING,dimension(:),allocatable::tmpval
    integer,dimension(:),allocatable::tmpid
    integer iup, ilow, imid
    COOP_UNKNOWN_STRING key, val
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
       dict%val(dict%id(ilow)) = dict%val(dict%n)
       dict%n = dict%n - 1
    elseif(dict%key(dict%id(iup)) .eq. dict%key(dict%n))then
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
    integer ind
    integer iup, ilow, imid
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


  subroutine coop_dictionary_lookup_string(dict, key, val)
    class(coop_dictionary):: dict
    COOP_UNKNOWN_STRING, intent(IN):: key
    COOP_UNKNOWN_STRING, intent(OUT)::val
    integer ind
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
    integer ind
    integer val
    integer,optional::default_val
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
    integer ind
    real val
    real,optional::default_val
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
    integer ind
    COOP_REAL val
    COOP_REAL,optional::default_val
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
    integer ind
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
    integer ind
    integer,dimension(:),intent(OUT):: val
    integer,dimension(:),optional::default_val
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
    integer ind
    real,dimension(:),intent(OUT)::  val
    real,dimension(:),optional::default_val
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
    integer ind
    COOP_REAL,dimension(:),intent(OUT)::  val
    COOP_REAL,dimension(:),optional::default_val
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
    integer ind
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
    integer i
    if(allocated(dict%key))then
       do i=1, dict%n
          write(*,"(3A)")  trim(dict%key(dict%id(i))), " = ", trim(dict%val(dict%id(i)))
       enddo
    else
       write(*,*) "Empty coop_dictionary"
    endif
  end subroutine coop_dictionary_print

End module coop_list_mod

