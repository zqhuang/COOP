module tmp_mod
  type shape
     integer :: color
     logical :: filled
     integer :: x
     integer :: y
   contains
     procedure::init => init_shape
  end type shape

  type, EXTENDS ( shape ) :: rectangle
     integer :: length
     integer :: width
  end type rectangle

  type, EXTENDS ( rectangle ) :: square
  end type square


contains

  subroutine init_shape(sh)
    class(shape) sh
    sh%color = 1
    sh%filled = .false.
    sh%x = 0
    sh%y = 0
    select type(sh)
    type is (shape)
       sh%x = -1
       sh%y = -1
    class is(rectangle)
       sh%length = 1
       sh%width = 1
    class default
       stop "UNKNOWN type"
    end select
  end subroutine init_shape

end module tmp_mod
