module asy_utils
  use file_io_utils
  use list_utils
  use interpolation_utils
  implicit none

#include "utils.h"

  interface asy_error_bar
     module procedure asy_error_bar_s, asy_error_bar_d
  end interface asy_error_bar

  interface asy_interpolate_curve
     module procedure asy_interpolate_curve_d, asy_interpolate_curve_s
  end interface asy_interpolate_curve

  interface asy_gray_color
     module procedure asy_gray_color_s, asy_gray_color_d
  end interface asy_gray_color

  interface asy_rgb_color
     module procedure asy_rgb_color_s, asy_rgb_color_d
  end interface asy_rgb_color

  interface asy_label
     module procedure asy_label_s, asy_label_d
  end interface asy_label

  interface asy_legend
     module procedure asy_legend_s, asy_legend_d, asy_legend_default, asy_legend_location
  end interface asy_legend

  interface asy_dot
     module procedure asy_dot_s, asy_dot_d
  end interface asy_dot


  interface asy_line
     module procedure asy_line_s, asy_line_d
  end interface asy_line

  interface asy_labels
     module procedure asy_labels_s, asy_labels_d
  end interface asy_labels

  interface asy_dots
     module procedure asy_dots_s, asy_dots_d
  end interface asy_dots

  interface asy_lines
     module procedure asy_lines_s, asy_lines_d
  end interface asy_lines

  interface asy_contour
     module procedure asy_contour_s, asy_contour_d, asy_contour_mult_s, asy_contour_mult_d,  asy_contour_arr_d, asy_contour_arr_s, asy_contour_path
  end interface asy_contour

  interface asy_curve
     module procedure asy_curve_s, asy_curve_d
  end interface asy_curve

  interface asy_density
     module procedure asy_density_s, asy_density_d, asy_irregular_density_s, asy_irregular_density_d
  end interface asy_density

  interface asy_path_append
     module procedure asy_path_append_s,  asy_path_append_d
  end interface asy_path_append

  interface asy_topaxis
     module procedure asy_topaxis_s, asy_topaxis_d
  end interface asy_topaxis

  interface asy_rightaxis
     module procedure asy_rightaxis_s, asy_rightaxis_d
  end interface asy_rightaxis




  integer, parameter::asy_path_max_nclosed = 4096
  type asy_path
     integer nclosed
     integer,dimension(:),allocatable:: length
     type(list_realarr) l
  end type asy_path



contains


  subroutine asy_init(fp,  xmin, xmax, ymin, ymax, width, height, caption, xlabel, ylabel, xlog, ylog, zlog, doclip, nblocks)
    type(file_pointer) fp
    real,optional:: width, height
    UNKNOWN_STRING,optional:: caption, xlabel, ylabel
    integer,optional::nblocks
    logical,optional::xlog, ylog, zlog, doclip
    real,optional:: xmin, xmax, ymin, ymax
    real minx, maxx, miny, maxy
    real w, h
    character(len = 5) tmp
    if(present(width))then
       w = width
    else
       w = 6.
    endif
    if(present(height))then
       h = height
    else
       h = 5.
    endif
    write(fp%unit, "(2F12.2)") w, h
    if(present(caption))then
       if(caption.eq."")then
          write(fp%unit, "(A)") "NULL"
       else
          write(fp%unit, "(A)") trim(caption)
       endif
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(xlabel))then
       if(trim(xlabel).eq."")then
          write(fp%unit, "(A)") "NULL"
       else
          write(fp%unit, "(A)") trim(xlabel)
       endif
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(ylabel))then
       if(trim(ylabel).eq."")then
          write(fp%unit, "(A)") "NULL"
       else
          write(fp%unit, "(A)") trim(ylabel)
       endif
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(xlog))then
       if(xlog)then
          tmp(1:2) = "1 "
       else
          tmp(1:2) = "0 "
       endif
    else
       tmp(1:2) = "0 "
    endif
    if(present(ylog))then
       if(ylog)then
          tmp(3:4) = "1 "
       else
          tmp(3:4) = "0 "
       endif
    else
       tmp(3:4) = "0 "
    endif
    if(present(zlog))then
       if(zlog)then
          tmp(5:5) = "1"
       else
          tmp(5:5) = "0"
       endif
    else
       tmp(5:5) = "0"
    endif
    write(fp%unit, "(A)") tmp
    if(present(doclip))then
       if(doclip)then
          write(fp%unit, "(A)") "1"
       else
          write(fp%unit, "(A)") "0"
       endif
    else
       write(fp%unit, "(A)") "0"
    endif
    if(present(xmin))then
       minx = xmin
    else
       minx = 1.1e30
    endif
    if(present(xmax))then
       maxx = xmax
    else
       maxx = -1.1e30
    endif
    if(present(ymin))then
       miny = ymin
    else
       miny = 1.1e30
    endif
    if(present(ymax))then
       maxy = ymax
    else
       maxy = -1.1e30
    endif
    write(fp%unit,"(2G14.5)") minx, maxx
    write(fp%unit,"(2G14.5)") miny, maxy
    if(present(nblocks))then
       write(fp%unit,"(I5)") nblocks    !!plot only n blocks
    else
       write(fp%unit,"(I5)") 0   !!0 means any number of blocks
    endif
  end subroutine asy_init

  subroutine asy_dot_d(fp, x, y, color, symbol)
    type(file_pointer) fp
    integer n
    real(dl) x,y
    UNKNOWN_STRING,optional:: color
    UNKNOWN_STRING,optional::symbol
    write(fp%unit, "(A)") "DOTS"
    write(fp%unit, "(A)") "1"
    if(present(color))then
       write(fp%unit, "(A)") trim(color)
    else
       write(fp%unit, "(A)") "black"
    endif
    if(present(symbol))then
       write(fp%unit, "(A)") trim(symbol)
    else
       write(fp%unit, "(A)") "dot"
    endif
    write(fp%unit, "(2G14.5)") real(x), real(y)
  end subroutine asy_dot_d


  subroutine asy_dot_s(fp, x, y, color, symbol)
    type(file_pointer) fp
    integer n
    real x,y
    UNKNOWN_STRING,optional:: color
    UNKNOWN_STRING,optional::symbol
    write(fp%unit, "(A)") "DOTS"
    write(fp%unit, "(A)") "1"
    if(present(color))then
       write(fp%unit, "(A)") trim(color)
    else
       write(fp%unit, "(A)") "black"
    endif
    if(present(symbol))then
       write(fp%unit, "(A)") trim(symbol)
    else
       write(fp%unit, "(A)") "dot"
    endif
    write(fp%unit, "(2G14.5)") x, y
  end subroutine asy_dot_s

  subroutine asy_dots_d(fp, x, y, color, symbol)
    type(file_pointer) fp
    integer n,i
    real(dl),dimension(:),intent(IN)::x,y
    UNKNOWN_STRING,optional:: color
    UNKNOWN_STRING,optional::symbol
    n = getdim( "asy_dot_block", size(x), size(y))
    write(fp%unit, "(A)") "DOTS"
    write(fp%unit, "(I8)") n
    if(present(color))then
       write(fp%unit, "(A)") trim(color)
    else
       write(fp%unit, "(A)") "black"
    endif
    if(present(symbol))then
       write(fp%unit, "(A)") trim(symbol)
    else
       write(fp%unit, "(A)") "dot"
    endif
    do i=1,n
       write(fp%unit, "(2G14.5)") real(x(i)), real(y(i))
    enddo
  end subroutine asy_dots_d

  subroutine asy_dots_s(fp, x, y, color, symbol)
    type(file_pointer) fp
    integer n,i
    real,dimension(:),intent(IN)::x,y
    UNKNOWN_STRING,optional:: color
    UNKNOWN_STRING,optional::symbol
    n = getdim( "asy_dot_block", size(x), size(y))
    write(fp%unit, "(A)") "DOTS"
    write(fp%unit, "(I8)") n
    if(present(color))then
       write(fp%unit, "(A)") trim(color)
    else
       write(fp%unit, "(A)") "black"
    endif
    if(present(symbol))then
       write(fp%unit, "(A)") trim(symbol)
    else
       write(fp%unit, "(A)") "dot"
    endif
    do i=1,n
       write(fp%unit, "(2G14.5)") x(i), y(i)
    enddo
  end subroutine asy_dots_s


  subroutine asy_line_d(fp, xstart, ystart, xend, yend, color, linetype, linewidth)
    type(file_pointer) fp
    real(dl)::xstart, ystart, xend, yend
    UNKNOWN_STRING, optional::color, linetype
    real,optional::linewidth
    STRING lineproperty
    write(fp%unit, "(A)") "LINES"
    write(fp%unit, "(A)") "1"
    if(present(color))then
       lineproperty=trim(color)
    else
       lineproperty = "black"
    endif
    if(present(linetype))then
       lineproperty = trim(lineproperty)//"_"//trim(linetype)
    else
       lineproperty = trim(lineproperty)//"_solid"
    endif
    if(present(linewidth))then
       lineproperty = trim(lineproperty)//"_"//trim(num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    write(fp%unit, "(4G14.5)") real(xstart), real(ystart), real(xend), real(yend)
  end subroutine asy_line_d


  subroutine asy_line_s(fp, xstart, ystart, xend, yend, color, linetype, linewidth)
    type(file_pointer) fp
    real::xstart, ystart, xend, yend
    UNKNOWN_STRING, optional::color, linetype
    real,optional::linewidth
    STRING lineproperty
    write(fp%unit, "(A)") "LINES"
    write(fp%unit, "(A)") "1"
    if(present(color))then
       lineproperty=trim(color)
    else
       lineproperty = "black"
    endif
    if(present(linetype))then
       lineproperty = trim(lineproperty)//"_"//trim(linetype)
    else
       lineproperty = trim(lineproperty)//"_solid"
    endif
    if(present(linewidth))then
       lineproperty = trim(lineproperty)//"_"//trim(num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    write(fp%unit, "(4G14.5)") xstart, ystart, xend, yend
  end subroutine asy_line_s
  
  subroutine asy_lines_d(fp, xstart, ystart, xend, yend, color, linetype, linewidth)
    type(file_pointer) fp
    real(dl),dimension(:),intent(IN)::xstart, ystart, xend, yend
    UNKNOWN_STRING, optional::color, linetype
    real,optional::linewidth
    STRING lineproperty
    integer n, i
    write(fp%unit, "(A)") "LINES"
    n = getdim("asy_lines", size(xstart), size(ystart), size(xend), size(yend))
    write(fp%unit, "(I8)") n
    if(present(color))then
       lineproperty=trim(color)
    else
       lineproperty = "black"
    endif
    if(present(linetype))then
       lineproperty = trim(lineproperty)//"_"//trim(linetype)
    else
       lineproperty = trim(lineproperty)//"_solid"
    endif
    if(present(linewidth))then
       lineproperty = trim(lineproperty)//"_"//trim(num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    do i = 1, n
       write(fp%unit, "(4G14.5)") real(xstart(i)), real(ystart(i)), real(xend(i)), real(yend(i))
    enddo
  end subroutine asy_lines_d


  subroutine asy_lines_s(fp, xstart, ystart, xend, yend, color, linetype, linewidth)
    type(file_pointer) fp
    real,dimension(:),intent(IN)::xstart, ystart, xend, yend
    UNKNOWN_STRING, optional::color, linetype
    real,optional::linewidth
    STRING lineproperty
    integer n, i
    write(fp%unit, "(A)") "LINES"
    n = getdim("asy_lines", size(xstart), size(ystart), size(xend), size(yend))
    write(fp%unit, "(I8)") n
    if(present(color))then
       lineproperty=trim(color)
    else
       lineproperty = "black"
    endif
    if(present(linetype))then
       lineproperty = trim(lineproperty)//"_"//trim(linetype)
    else
       lineproperty = trim(lineproperty)//"_solid"
    endif
    if(present(linewidth))then
       lineproperty = trim(lineproperty)//"_"//trim(num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    do i = 1, n
       write(fp%unit, "(4G14.5)") xstart(i), ystart(i), xend(i), yend(i)
    enddo
  end subroutine asy_lines_s


  subroutine asy_interpolate_curve_d(fp, xraw, yraw, interpolate, color, linetype, linewidth, legend)
    integer,parameter::n = 256
    real(dl) x(n), y(n),  minx, maxx, dx
    integer w(n)
    real(dl),dimension(:),allocatable::xx, yy, yy2
    type(file_pointer) fp
    real(dl),dimension(:),intent(IN)::xraw, yraw
    UNKNOWN_STRING interpolate
    UNKNOWN_STRING,optional:: color, linetype
    real,optional::linewidth
    UNKNOWN_STRING, optional::legend
    integer i, j, m, nn, npt
    STRING lineproperty
    logical do_raw
    m = GetDim("asy_interpolate_curve", size(xraw), size(yraw))
    minx = minval(xraw)
    maxx = maxval(xraw)
    y = 0.
    w = 0
    if(minx .ge. maxx)then
       npt = 1
       x(1) = minx
       select case(trim(interpolate))
       case("LinearLog", "LogLog")
          y(1) = exp(sum(log(yraw))/m)
       case default
          y(1) = sum(yraw)/m
       end select
       goto 100
    endif
    npt = n
    do_raw = .false.
    select case(trim(interpolate))
    case("LinearLinear")
       call set_uniform(n, x, minx, maxx)
       dx = (maxx - minx)/(n-1.)
       do i=1, m
          j = nint((xraw(i) - minx)/dx)+1
          y(j) = y(j) + yraw(i)
          w(j) = w(j) + 1
       enddo
       where(w.ne.0)
          y = y/w
       end where
       if(any(w.eq.0))then
          m = count(w.ne.0)
          allocate(xx(m), yy(m), yy2(m))
          j = 1
          do i=1, n
             if(w(i).ne.0)then
                xx(j) = x(i)
                yy(j) = y(i)
                j = j + 1
             endif
          enddo
          call splines(xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call splints(xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
    case("LinearLog")
       call set_uniform(n, x, minx, maxx)
       dx = (maxx - minx)/(n-1.)
       do i=1, m
          j = nint((xraw(i) - minx)/dx)+1
          y(j) = y(j) + log(yraw(i))
          w(j) = w(j) + 1
       enddo
       where(w.ne.0)
          y = y/w
       end where
       if(any(w.eq.0))then
          m = count(w.ne.0)
          allocate(xx(m), yy(m), yy2(m))
          j = 1
          do i=1, n
             if(w(i).ne.0)then
                xx(j) = x(i)
                yy(j) = y(i)
                j = j + 1
             endif
          enddo
          call splines(xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call splints(xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
       y = exp(y)
    case("LogLinear")
       minx = log(minx)
       maxx = log(maxx)
       call set_uniform(n, x, minx, maxx)
       dx = (maxx - minx)/(n-1.)
       do i=1, m
          j = nint((log(xraw(i)) - minx)/dx)+1
          y(j) = y(j) + yraw(i)
          w(j) = w(j) + 1
       enddo
       where(w.ne.0)
          y = y/w
       end where
       if(any(w.eq.0))then
          m = count(w.ne.0)
          allocate(xx(m), yy(m), yy2(m))
          j = 1
          do i=1, n
             if(w(i).ne.0)then
                xx(j) = x(i)
                yy(j) = y(i)
                j = j + 1
             endif
          enddo
          call splines(xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call splints(xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
       x = exp(x)
    case("LogLog")
       minx = log(minx)
       maxx = log(maxx)
       call set_uniform(n, x, minx, maxx)
       dx = (maxx - minx)/(n-1.)
       do i=1, m
          j = nint((log(xraw(i)) - minx)/dx)+1
          y(j) = y(j) + log(yraw(i))
          w(j) = w(j) + 1
       enddo
       where(w.ne.0)
          y = y/w
       end where
       if(any(w.eq.0))then
          m = count(w.ne.0)
          allocate(xx(m), yy(m), yy2(m))
          j = 1
          do i=1, n
             if(w(i).ne.0)then
                xx(j) = x(i)
                yy(j) = y(i)
                j = j + 1
             endif
          enddo
          call splines(xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call splints(xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
       x = exp(x)
       y = exp(y)
    case default
       do_raw = .true.
       npt = m
    end select

100 write(fp%unit, "(A)") "CURVE"
    write(fp%unit, "(I8)") npt
    if(present(legend))then
       if(trim(legend).ne."")then
          write(fp%unit, "(A)") trim(legend)
       else
          write(fp%unit, "(A)") "NULL"
       endif
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(color))then
       lineproperty=trim(color)
    else
       lineproperty = "black"
    endif
    if(present(linetype))then
       lineproperty = trim(lineproperty)//"_"//trim(linetype)
    else
       lineproperty = trim(lineproperty)//"_solid"
    endif
    if(present(linewidth))then
       lineproperty = trim(lineproperty)//"_"//trim(num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    write(fp%unit, "(A)") "0"
    if(do_raw)then
       do i=1,npt
          write(fp%unit,"(2G14.5)") real(xraw(i)), real(yraw(i))
       enddo
    else
       do i=1,npt
          write(fp%unit,"(2G14.5)") real(x(i)), real(y(i))
       enddo
    endif
  end subroutine asy_interpolate_curve_d

  subroutine asy_interpolate_curve_s(fp, xraw, yraw, interpolate, color, linetype, linewidth, legend)
    integer,parameter::n = 256
    real(dl) x(n), y(n),  minx, maxx, dx
    integer w(n)
    real(dl),dimension(:),allocatable::xx, yy, yy2
    type(file_pointer) fp
    real,dimension(:),intent(IN)::xraw, yraw
    UNKNOWN_STRING interpolate
    UNKNOWN_STRING,optional:: color, linetype
    real,optional::linewidth
    UNKNOWN_STRING, optional::legend
    integer i, j, m, nn, npt
    STRING lineproperty
    logical do_raw
    m = GetDim("asy_interpolate_curve", size(xraw), size(yraw))
    minx = minval(xraw)
    maxx = maxval(xraw)
    y = 0.
    w = 0
    if(minx .ge. maxx)then
       npt = 1
       x(1) = minx
       select case(trim(interpolate))
       case("LinearLog", "LogLog")
          y(1) = exp(sum(log(yraw))/m)
       case default
          y(1) = sum(yraw)/m
       end select
       goto 100
    endif
    npt = n
    do_raw = .false.
    select case(trim(interpolate))
    case("LinearLinear")
       call set_uniform(n, x, minx, maxx)
       dx = (maxx - minx)/(n-1.)
       do i=1, m
          j = nint((xraw(i) - minx)/dx)+1
          y(j) = y(j) + yraw(i)
          w(j) = w(j) + 1
       enddo
       where(w.ne.0)
          y = y/w
       end where
       if(any(w.eq.0))then
          m = count(w.ne.0)
          allocate(xx(m), yy(m), yy2(m))
          j = 1
          do i=1, n
             if(w(i).ne.0)then
                xx(j) = x(i)
                yy(j) = y(i)
                j = j + 1
             endif
          enddo
          call splines(xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call splints(xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
    case("LinearLog")
       call set_uniform(n, x, minx, maxx)
       dx = (maxx - minx)/(n-1.)
       do i=1, m
          j = nint((xraw(i) - minx)/dx)+1
          y(j) = y(j) + log(yraw(i))
          w(j) = w(j) + 1
       enddo
       where(w.ne.0)
          y = y/w
       end where
       if(any(w.eq.0))then
          m = count(w.ne.0)
          allocate(xx(m), yy(m), yy2(m))
          j = 1
          do i=1, n
             if(w(i).ne.0)then
                xx(j) = x(i)
                yy(j) = y(i)
                j = j + 1
             endif
          enddo
          call splines(xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call splints(xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
       y = exp(y)
    case("LogLinear")
       minx = log(minx)
       maxx = log(maxx)
       call set_uniform(n, x, minx, maxx)
       dx = (maxx - minx)/(n-1.)
       do i=1, m
          j = nint((log(xraw(i)) - minx)/dx)+1
          y(j) = y(j) + yraw(i)
          w(j) = w(j) + 1
       enddo
       where(w.ne.0)
          y = y/w
       end where
       if(any(w.eq.0))then
          m = count(w.ne.0)
          allocate(xx(m), yy(m), yy2(m))
          j = 1
          do i=1, n
             if(w(i).ne.0)then
                xx(j) = x(i)
                yy(j) = y(i)
                j = j + 1
             endif
          enddo
          call splines(xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call splints(xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
       x = exp(x)
    case("LogLog")
       minx = log(minx)
       maxx = log(maxx)
       call set_uniform(n, x, minx, maxx)
       dx = (maxx - minx)/(n-1.)
       do i=1, m
          j = nint((log(xraw(i)) - minx)/dx)+1
          y(j) = y(j) + log(yraw(i))
          w(j) = w(j) + 1
       enddo
       where(w.ne.0)
          y = y/w
       end where
       if(any(w.eq.0))then
          m = count(w.ne.0)
          allocate(xx(m), yy(m), yy2(m))
          j = 1
          do i=1, n
             if(w(i).ne.0)then
                xx(j) = x(i)
                yy(j) = y(i)
                j = j + 1
             endif
          enddo
          call splines(xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call splints(xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
       x = exp(x)
       y = exp(y)
    case default
       do_raw = .true.
       npt = m
    end select

100 write(fp%unit, "(A)") "CURVE"
    write(fp%unit, "(I8)") npt
    if(present(legend))then
       if(trim(legend).ne."")then
          write(fp%unit, "(A)") trim(legend)
       else
          write(fp%unit, "(A)") "NULL"
       endif
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(color))then
       lineproperty=trim(color)
    else
       lineproperty = "black"
    endif
    if(present(linetype))then
       lineproperty = trim(lineproperty)//"_"//trim(linetype)
    else
       lineproperty = trim(lineproperty)//"_solid"
    endif
    if(present(linewidth))then
       lineproperty = trim(lineproperty)//"_"//trim(num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    write(fp%unit, "(A)") "0"
    if(do_raw)then
       do i=1,npt
          write(fp%unit,"(2G14.5)") xraw(i), yraw(i)
       enddo
    else
       do i=1,npt
          write(fp%unit,"(2G14.5)") real(x(i)), real(y(i))
       enddo
    endif
  end subroutine asy_interpolate_curve_s


  subroutine asy_curve_from_file(fp, fname, interpolate, xcol, ycol, color, linetype, linewidth, legend)
    integer,parameter::n = 256
    integer, optional::xcol, ycol
    real(dl) x(n), y(n),  minx, maxx, dx
    integer w(n)
    real(dl),dimension(:),allocatable::xx, yy, yy2
    type(file_pointer) fp, fp2
    real(dl),dimension(:),allocatable::xraw, yraw, line
    UNKNOWN_STRING interpolate, fname
    UNKNOWN_STRING,optional:: color, linetype
    real,optional::linewidth
    UNKNOWN_STRING, optional::legend
    integer i, j, m, nn, npt, ncols, xl, yl
    STRING lineproperty
    logical do_raw
    if(.not. file_exists(fname))then
       write(*,*) "asy_curve_from_file: the data file "//trim(fname)//" does not exist"
       stop 
    endif
    m = file_numlines(fname)
    ncols = file_numColumns(fname)
    if(ncols .eq. 0)then
       write(*,*) "asy_curve_from_file: the data file "//trim(fname)//" is empty"
       stop 
    endif
    if(present(xcol))then
       if(xcol .gt. ncols) stop "asy_curve_from_file: xcol > ncol"
       xl = xcol
    else
       x = min(ncols-1, 1)
    endif
    if(present(ycol))then
       if(ycol .gt. ncols) stop "asy_curve_from_file: ycol > ncol"
       yl = ycol
    else
       yl = min(ncols, 2)
    endif
    allocate(line(0:ncols))
    allocate(xraw(m), yraw(m))
    fp2 = open_file(trim(fname), "r")
    do i = 1, m
       if(file_readoneline_double_arr(fp2, line(1:ncols)))then
          line(0) = i
          xraw(i) = line(xl)
          yraw(i) = line(yl)
       else
          write(*,*) trim(fname)//" has a wrong format"
          stop
       endif
    enddo
    call close_file(fp2)
    minx = minval(xraw)
    maxx = maxval(xraw)
    y = 0.
    w = 0
    if(minx .ge. maxx)then
       npt = 1
       x(1) = minx
       select case(trim(interpolate))
       case("LinearLog", "LogLog")
          y(1) = exp(sum(log(yraw))/m)
       case default
          y(1) = sum(yraw)/m
       end select
       goto 100
    endif
    npt = n
    do_raw = .false.
    select case(trim(interpolate))
    case("LinearLinear")
       call set_uniform(n, x, minx, maxx)
       dx = (maxx - minx)/(n-1.)
       do i=1, m
          j = nint((xraw(i) - minx)/dx)+1
          y(j) = y(j) + yraw(i)
          w(j) = w(j) + 1
       enddo
       where(w.ne.0)
          y = y/w
       end where
       if(any(w.eq.0))then
          m = count(w.ne.0)
          allocate(xx(m), yy(m), yy2(m))
          j = 1
          do i=1, n
             if(w(i).ne.0)then
                xx(j) = x(i)
                yy(j) = y(i)
                j = j + 1
             endif
          enddo
          call splines(xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call splints(xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
    case("LinearLog")
       call set_uniform(n, x, minx, maxx)
       dx = (maxx - minx)/(n-1.)
       do i=1, m
          j = nint((xraw(i) - minx)/dx)+1
          y(j) = y(j) + log(yraw(i))
          w(j) = w(j) + 1
       enddo
       where(w.ne.0)
          y = y/w
       end where
       if(any(w.eq.0))then
          m = count(w.ne.0)
          allocate(xx(m), yy(m), yy2(m))
          j = 1
          do i=1, n
             if(w(i).ne.0)then
                xx(j) = x(i)
                yy(j) = y(i)
                j = j + 1
             endif
          enddo
          call splines(xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call splints(xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
       y = exp(y)
    case("LogLinear")
       minx = log(minx)
       maxx = log(maxx)
       call set_uniform(n, x, minx, maxx)
       dx = (maxx - minx)/(n-1.)
       do i=1, m
          j = nint((log(xraw(i)) - minx)/dx)+1
          y(j) = y(j) + yraw(i)
          w(j) = w(j) + 1
       enddo
       where(w.ne.0)
          y = y/w
       end where
       if(any(w.eq.0))then
          m = count(w.ne.0)
          allocate(xx(m), yy(m), yy2(m))
          j = 1
          do i=1, n
             if(w(i).ne.0)then
                xx(j) = x(i)
                yy(j) = y(i)
                j = j + 1
             endif
          enddo
          call splines(xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call splints(xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
       x = exp(x)
    case("LogLog")
       minx = log(minx)
       maxx = log(maxx)
       call set_uniform(n, x, minx, maxx)
       dx = (maxx - minx)/(n-1.)
       do i=1, m
          j = nint((log(xraw(i)) - minx)/dx)+1
          y(j) = y(j) + log(yraw(i))
          w(j) = w(j) + 1
       enddo
       where(w.ne.0)
          y = y/w
       end where
       if(any(w.eq.0))then
          m = count(w.ne.0)
          allocate(xx(m), yy(m), yy2(m))
          j = 1
          do i=1, n
             if(w(i).ne.0)then
                xx(j) = x(i)
                yy(j) = y(i)
                j = j + 1
             endif
          enddo
          call splines(xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call splints(xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
       x = exp(x)
       y = exp(y)
    case default
       do_raw = .true.
       npt = m
    end select

100 write(fp%unit, "(A)") "CURVE"
    write(fp%unit, "(I8)") npt
    if(present(legend))then
       if(trim(legend).ne."")then
          write(fp%unit, "(A)") trim(legend)
       else
          write(fp%unit, "(A)") "NULL"
       endif
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(color))then
       lineproperty=trim(color)
    else
       lineproperty = "black"
    endif
    if(present(linetype))then
       lineproperty = trim(lineproperty)//"_"//trim(linetype)
    else
       lineproperty = trim(lineproperty)//"_solid"
    endif
    if(present(linewidth))then
       lineproperty = trim(lineproperty)//"_"//trim(num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    write(fp%unit, "(A)") "0"
    if(do_raw)then
       do i=1,npt
          write(fp%unit,"(2G14.5)") real(xraw(i)), real(yraw(i))
       enddo
    else
       do i=1,npt
          write(fp%unit,"(2G14.5)") real(x(i)), real(y(i))
       enddo
    endif
    deallocate(xraw, yraw, line)
  end subroutine asy_curve_from_file


  subroutine asy_curve_d(fp, x, y, smooth, color, linetype, linewidth, legend)
    type(file_pointer) fp
    logical,optional::smooth
    real(dl),dimension(:),intent(IN)::x,y
    UNKNOWN_STRING,optional:: color, linetype
    real,optional::linewidth
    UNKNOWN_STRING, optional::legend
    integer i,n
    STRING lineproperty
    n = getdim("asy_curve", size(x), size(y))
    write(fp%unit, "(A)") "CURVE"
    write(fp%unit, "(I8)") n
    if(present(legend))then
       if(trim(legend).ne."")then
          write(fp%unit, "(A)") trim(legend)
       else
          write(fp%unit, "(A)") "NULL"
       endif
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(color))then
       lineproperty=trim(color)
    else
       lineproperty = "black"
    endif
    if(present(linetype))then
       lineproperty = trim(lineproperty)//"_"//trim(linetype)
    else
       lineproperty = trim(lineproperty)//"_solid"
    endif
    if(present(linewidth))then
       lineproperty = trim(lineproperty)//"_"//trim(num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    if(present(smooth))then
       if(smooth)then
          write(fp%unit, "(A)") "1"
       else
          write(fp%unit, "(A)") "0"
       endif
    else
       write(fp%unit, "(A)") "0"  !!no smoothing by default
    endif
    do i=1,n
       write(fp%unit,"(2G14.5)") real(x(i)), real(y(i))
    enddo
  end subroutine asy_curve_d

  subroutine asy_curve_s(fp, x, y, smooth, color, linetype, linewidth, legend)
    type(file_pointer) fp
    logical,optional::smooth
    real,dimension(:),intent(IN)::x,y
    UNKNOWN_STRING,optional:: color, linetype
    real,optional::linewidth
    integer i,n
    UNKNOWN_STRING,optional::legend
    STRING lineproperty
    n = getdim("asy_curve", size(x), size(y))
    write(fp%unit, "(A)") "CURVE"
    write(fp%unit, "(I8)") n
    if(present(legend))then
       if(trim(legend).ne."")then
          write(fp%unit, "(A)") trim(legend)
       else
          write(fp%unit, "(A)") "NULL"
       endif
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(color))then
       lineproperty=trim(color)
    else
       lineproperty = "black"
    endif
    if(present(linetype))then
       lineproperty = trim(lineproperty)//"_"//trim(linetype)
    else
       lineproperty = trim(lineproperty)//"_solid"
    endif
    if(present(linewidth))then
       lineproperty = trim(lineproperty)//"_"//trim(num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    if(present(smooth))then
       if(smooth)then
          write(fp%unit, "(A)") "1"
       else
          write(fp%unit, "(A)") "0"
       endif
    else
       write(fp%unit, "(A)") "0"  !!no smoothing by default
    endif
    do i=1,n
       write(fp%unit,"(2G14.5)") x(i), y(i)
    enddo
  end subroutine asy_curve_s

  subroutine asy_labels_d(fp, labels, x, y, color)
    STRING, dimension(:),intent(IN)::labels
    real(dl),dimension(:),intent(IN)::x, y
    type(file_pointer) fp
    UNKNOWN_STRING,optional::color
    integer n, i
    write(fp%unit, "(A)") "LABELS"
    n = getdim("asy_labels", size(x), size(y), size(labels))
    write(fp%unit,"(I8)") n
    if(present(color))then
       write(fp%unit, "(A)") trim(color)
    else
       write(fp%unit, "(A)") "black"
    endif
    do i=1,n
       write(fp%unit, "(2G14.5)") real(x(i)), real(y(i))
       if(trim(labels(i)).eq."")then
          write(fp%unit, "(A)") "NULLL"
       else
          write(fp%unit, "(A)") trim(labels(i))
       endif
    enddo
  end subroutine asy_labels_d

  subroutine asy_labels_s(fp, labels, x, y, color)
    STRING, dimension(:),intent(IN)::labels
    real,dimension(:),intent(IN)::x, y
    type(file_pointer) fp
    UNKNOWN_STRING,optional::color
    integer n, i
    write(fp%unit, "(A)") "LABELS"
    n = getdim("asy_labels", size(x), size(y), size(labels))
    write(fp%unit,"(I8)") n
    if(present(color))then
       write(fp%unit, "(A)") trim(color)
    else
       write(fp%unit, "(A)") "black"
    endif
    do i=1,n
       write(fp%unit, "(2G14.5)") x(i), y(i)
       if(trim(labels(i)).eq."")then
          write(fp%unit, "(A)") "NULLL"
       else
          write(fp%unit, "(A)") trim(labels(i))
       endif
    enddo
  end subroutine asy_labels_s


  subroutine asy_label_d(fp, label, x, y, color)
    type(file_pointer) fp
    UNKNOWN_STRING label
    UNKNOWN_STRING,optional::color
    real(dl) x, y
    write(fp%unit, "(A)") "LABELS"
    write(fp%unit, "(A)") "1"
    if(present(color))then
       write(fp%unit, "(A)") trim(color)
    else
       write(fp%unit, "(A)") "black"
    endif
    write(fp%unit, "(2G14.5)") real(x), real(y)
    if(trim(label).eq."")then
       write(fp%unit, "(A)") "NULL"
    else
       write(fp%unit, "(A)") trim(label)
    endif
  end subroutine asy_label_d


  subroutine asy_label_s(fp, label, x, y, color)
    type(file_pointer) fp
    UNKNOWN_STRING label
    UNKNOWN_STRING,optional::color
    real x, y
    write(fp%unit, "(A)") "LABELS"
    write(fp%unit, "(A)") "1"
    if(present(color))then
       write(fp%unit, "(A)") trim(color)
    else
       write(fp%unit, "(A)") "black"
    endif
    write(fp%unit, "(2G14.5)") x, y
    if(trim(label).eq."")then
       write(fp%unit, "(A)") "NULL"
    else
       write(fp%unit, "(A)") trim(label)
    endif
  end subroutine asy_label_s



  subroutine asy_contour_d(fp, x, y, colorfill, smooth, linecolor, linetype, linewidth)
    type(file_pointer) fp
    logical,optional::smooth
    real(dl),dimension(:),intent(IN)::x,y
    UNKNOWN_STRING,optional::colorfill
    UNKNOWN_STRING,optional:: linecolor, linetype
    real,optional::linewidth
    integer i,n
    STRING lineproperty
    n = getdim("asy_contour", size(x), size(y))
    write(fp%unit, "(A)") "CONTOUR"
    write(fp%unit, "(A)") "1" !!type 1 contour
    if(present(colorfill))then
       write(fp%unit, "(A)") trim(colorfill)
    else
       write(fp%unit, "(A)") "black"
    endif
    if(present(linecolor).or. present(linetype).or.present(linewidth))then
       if(present(linecolor))then
          lineproperty=trim(linecolor)
       else
          lineproperty = "black"
       endif
       if(present(linetype))then
          lineproperty = trim(lineproperty)//"_"//trim(linetype)
       else
          lineproperty = trim(lineproperty)//"_solid"
       endif
       if(present(linewidth))then
          lineproperty = trim(lineproperty)//"_"//trim(num2str(linewidth))
       endif
       write(fp%unit, "(A)") trim(lineproperty)
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(smooth))then
       if(smooth)then
          write(fp%unit, "(A)") "1"
       else
          write(fp%unit, "(A)") "0"
       endif
    else
       write(fp%unit, "(A)") "0"  !!no smoothing by default
    endif
    write(fp%unit, "(A)") "1"  !!1 path 
    write(fp%unit, "(I8)") n
    do i=1,n
       write(fp%unit,"(2G14.5)") real(x(i)), real(y(i))
    enddo
  end subroutine asy_contour_d

  subroutine asy_contour_s(fp, x, y, colorfill, smooth, linecolor, linetype, linewidth)
    type(file_pointer) fp
    logical,optional::smooth
    real,dimension(:),intent(IN)::x,y
    UNKNOWN_STRING,optional::colorfill
    UNKNOWN_STRING,optional:: linecolor, linetype
    real,optional::linewidth
    integer i,n
    STRING lineproperty
    n = getdim("asy_contour", size(x), size(y))
    write(fp%unit, "(A)") "CONTOUR"
    write(fp%unit, "(A)") "1"   !!type I contour
    if(present(colorfill))then
       write(fp%unit, "(A)") trim(colorfill)
    else
       write(fp%unit, "(A)") "black"
    endif
    if(present(linecolor).or. present(linetype).or.present(linewidth))then
       if(present(linecolor))then
          lineproperty=trim(linecolor)
       else
          lineproperty = "black"
       endif
       if(present(linetype))then
          lineproperty = trim(lineproperty)//"_"//trim(linetype)
       else
          lineproperty = trim(lineproperty)//"_solid"
       endif
       if(present(linewidth))then
          lineproperty = trim(lineproperty)//"_"//trim(num2str(linewidth))
       endif
       write(fp%unit, "(A)") trim(lineproperty)
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(smooth))then
       if(smooth)then
          write(fp%unit, "(A)") "1"
       else
          write(fp%unit, "(A)") "0"
       endif
    else
       write(fp%unit, "(A)") "0"  !!no smoothing by default
    endif
    write(fp%unit, "(A)") "1"  !!1 path 
    write(fp%unit, "(I8)") n
    do i=1,n
       write(fp%unit,"(2G14.5)") x(i), y(i)
    enddo
  end subroutine asy_contour_s


  subroutine asy_contour_path(fp, path, colorfill, smooth, linecolor, linetype, linewidth)
    type(file_pointer) fp
    logical,optional::smooth
    type(asy_path) path
    UNKNOWN_STRING,optional::colorfill
    UNKNOWN_STRING,optional:: linecolor, linetype
    real,optional::linewidth
    integer i, ipath, pl
    STRING lineproperty
    if(path%nclosed.eq.0 .or. (.not.list_is_initialized(path%l)) .or. path%l%n .eq. 0)return
    if(path%l%dim .ne.2) stop "asy_contour_path: wrong dimension of list in path"
    write(fp%unit, "(A)") "CONTOUR"
    write(fp%unit, "(A)") "1"   !!type I contour
    if(present(colorfill))then
       write(fp%unit, "(A)") trim(colorfill)
    else
       write(fp%unit, "(A)") "black"
    endif
    if(present(linecolor).or. present(linetype).or.present(linewidth))then
       if(present(linecolor))then
          lineproperty=trim(linecolor)
       else
          lineproperty = "black"
       endif
       if(present(linetype))then
          lineproperty = trim(lineproperty)//"_"//trim(linetype)
       else
          lineproperty = trim(lineproperty)//"_solid"
       endif
       if(present(linewidth))then
          lineproperty = trim(lineproperty)//"_"//trim(num2str(linewidth))
       endif
       write(fp%unit, "(A)") trim(lineproperty)
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(smooth))then
       if(smooth)then
          write(fp%unit, "(A)") "1"
       else
          write(fp%unit, "(A)") "0"
       endif
    else
       write(fp%unit, "(A)") "0"  !!no smoothing by default
    endif
    write(fp%unit, "(I8)") path%nclosed
    i = 1
    do ipath = 1, path%nclosed
       write(fp%unit, "(I8)") path%length(ipath)
       do pl = 1, path%length(ipath)
          write(fp%unit,"(2G14.5)") list_element(path%l, i)
          i = i + 1
       enddo
    enddo
  end subroutine asy_contour_path

  subroutine asy_contour_mult_d(fp, x, y, colorfill, smooth, linecolor, linetype, linewidth)
    type(file_pointer) fp
    logical,optional::smooth
    real(dl),dimension(:,:),intent(IN)::x,y
    UNKNOWN_STRING,optional::colorfill
    UNKNOWN_STRING,optional:: linecolor, linetype
    real,optional::linewidth
    integer i,n, m,ipath
    STRING lineproperty
    n = getdim("asy_contour", size(x,1), size(y,1))
    m = getdim("asy_contour", size(x,2), size(y,2))
    write(fp%unit, "(A)") "CONTOUR"
    write(fp%unit, "(A)") "1"   !!type I contour
    if(present(colorfill))then
       write(fp%unit, "(A)") trim(colorfill)
    else
       write(fp%unit, "(A)") "black"
    endif
    if(present(linecolor).or. present(linetype).or.present(linewidth))then
       if(present(linecolor))then
          lineproperty=trim(linecolor)
       else
          lineproperty = "black"
       endif
       if(present(linetype))then
          lineproperty = trim(lineproperty)//"_"//trim(linetype)
       else
          lineproperty = trim(lineproperty)//"_solid"
       endif
       if(present(linewidth))then
          lineproperty = trim(lineproperty)//"_"//trim(num2str(linewidth))
       endif
       write(fp%unit, "(A)") trim(lineproperty)
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(smooth))then
       if(smooth)then
          write(fp%unit, "(A)") "1"
       else
          write(fp%unit, "(A)") "0"
       endif
    else
       write(fp%unit, "(A)") "0"  !!no smoothing by default
    endif
    write(fp%unit, "(I8)") m
    do ipath = 1, m
       write(fp%unit, "(I8)") n
       do i=1,n
          write(fp%unit,"(2G14.5)") real(x(i, ipath)), real(y(i, ipath))
       enddo
    enddo
  end subroutine asy_contour_mult_d

  subroutine asy_contour_mult_s(fp, x, y, colorfill, smooth, linecolor, linetype, linewidth)
    type(file_pointer) fp
    logical,optional::smooth
    real,dimension(:,:),intent(IN)::x,y
    UNKNOWN_STRING,optional::colorfill
    UNKNOWN_STRING,optional:: linecolor, linetype
    real,optional::linewidth
    integer i,n, m,ipath
    STRING lineproperty
    n = getdim("asy_contour", size(x,1), size(y,1))
    m = getdim("asy_contour", size(x,2), size(y,2))
    write(fp%unit, "(A)") "CONTOUR"
    write(fp%unit, "(A)") "1"   !!type I contour
    if(present(colorfill))then
       write(fp%unit, "(A)") trim(colorfill)
    else
       write(fp%unit, "(A)") "black"
    endif
    if(present(linecolor).or. present(linetype).or.present(linewidth))then
       if(present(linecolor))then
          lineproperty=trim(linecolor)
       else
          lineproperty = "black"
       endif
       if(present(linetype))then
          lineproperty = trim(lineproperty)//"_"//trim(linetype)
       else
          lineproperty = trim(lineproperty)//"_solid"
       endif
       if(present(linewidth))then
          lineproperty = trim(lineproperty)//"_"//trim(num2str(linewidth))
       endif
       write(fp%unit, "(A)") trim(lineproperty)
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(smooth))then
       if(smooth)then
          write(fp%unit, "(A)") "1"
       else
          write(fp%unit, "(A)") "0"
       endif
    else
       write(fp%unit, "(A)") "0"  !!no smoothing by default
    endif
    write(fp%unit, "(I8)") m
    do ipath = 1, m
       write(fp%unit, "(I8)") n
       do i=1,n
          write(fp%unit,"(2G14.5)") x(i, ipath), y(i, ipath)
       enddo
    enddo
  end subroutine asy_contour_mult_s


  subroutine asy_contour_arr_d(fp, z, xmin, xmax, ymin, ymax, cvals, colorfill, smooth, linecolor, linetype, linewidth)
    type(file_pointer) fp
    real(dl),dimension(:,:)::z
    real(dl) xmin, xmax, ymin, ymax
    real(dl) cvals(:)
    integer nx, ny, i, nc
    STRING,dimension(:)::colorfill
    STRING,dimension(:)::linecolor, linetype
    real,dimension(:),optional::linewidth
    logical smooth
    nx = size(z, 1)
    ny = size(z, 2)
    write(fp%unit, "(A)") "CONTOUR"
    write(fp%unit, "(A)") "2" !!type 2 contour
    write(fp%unit, "(2G14.5)") real(xmin), real(xmax)
    write(fp%unit, "(2G14.5)") real(ymin), real(ymax)
    if(present(linewidth))then
       nc = getdim("asy_contour_arr_d", size(cvals), size(colorfill), size(linecolor), size(linetype), size(linewidth))
    else
       nc = getdim("asy_contour_arr_d", size(cvals), size(colorfill), size(linecolor), size(linetype))
    endif
    write(fp%unit, "(I8)") nc
    write(fp%unit, "("//trim(num2str(nc))//"G14.5)") cvals
    do i = 1, nc
       write(fp%unit, "(A)") trim(colorfill(i))
       if(present(linewidth))then
          write(fp%unit, "(A)") trim(asy_linestyle(linecolor(i), linetype(i), linewidth(i)))
       else
          write(fp%unit, "(A)") trim(asy_linestyle(linecolor(i), linetype(i)))
       endif
    enddo
    if(smooth)then
       write(fp%unit, "(A)") "1"
    else
       write(fp%unit, "(A)") "0"
    endif
    write(fp%unit, "(2I8)") nx, ny
    do i=1,nx
       write(fp%unit, "("//trim(num2str(ny))//"G14.5)") real(z(i,:))
    enddo
  end subroutine asy_contour_arr_d
  


  subroutine asy_contour_arr_s(fp, z, xmin, xmax, ymin, ymax, cvals, colorfill, smooth, linecolor, linetype, linewidth)
    type(file_pointer) fp
    real,dimension(:,:)::z
    real xmin, xmax, ymin, ymax
    real cvals(:)
    integer nx, ny, i, nc
    STRING,dimension(:)::colorfill
    STRING,dimension(:)::linecolor, linetype
    real,dimension(:),optional::linewidth
    logical smooth
    nx = size(z, 1)
    ny = size(z, 2)
    write(fp%unit, "(A)") "CONTOUR"
    write(fp%unit, "(A)") "2" !!type 2 contour
    write(fp%unit, "(2G14.5)") real(xmin), real(xmax)
    write(fp%unit, "(2G14.5)") real(ymin), real(ymax)
    if(present(linewidth))then
       nc = getdim("asy_contour_arr_s", size(cvals), size(colorfill), size(linecolor), size(linetype), size(linewidth))
    else
       nc = getdim("asy_contour_arr_s", size(cvals), size(colorfill), size(linecolor), size(linetype))
    endif
    write(fp%unit, "(I8)") nc
    write(fp%unit, "("//trim(num2str(nc))//"G14.5)") cvals
    do i = 1, nc
       write(fp%unit, "(A)") trim(colorfill(i))
       if(present(linewidth))then
          write(fp%unit, "(A)") trim(asy_linestyle(linecolor(i), linetype(i), linewidth(i)))
       else
          write(fp%unit, "(A)") trim(asy_linestyle(linecolor(i), linetype(i)))
       endif
    enddo
    if(smooth)then
       write(fp%unit, "(A)") "1"
    else
       write(fp%unit, "(A)") "0"
    endif
    write(fp%unit, "(2I8)") nx, ny
    do i=1,nx
       write(fp%unit, "("//trim(num2str(ny))//"G14.5)") z(i,:)
    enddo
  end subroutine asy_contour_arr_s




  subroutine asy_density_d(fp, z, xmin, xmax, ymin, ymax, label, zmin, zmax)
    type(file_pointer) fp
    real(dl),dimension(:,:)::z
    real(dl) xmin, xmax, ymin, ymax
    real(dl),optional::zmin, zmax
    UNKNOWN_STRING, optional::label
    real(dl) minz, maxz
    integer nx, ny, i
    nx = size(z, 1)
    ny = size(z, 2)
    write(fp%unit, "(A)") "DENSITY"
    if(present(label))then
       if(trim(label).ne."")then
          write(fp%unit,"(A)") trim(label)
       else
          write(fp%unit,"(A)") "NULL"
       endif
    else
       write(fp%unit,"(A)") "NULL"
    endif
    write(fp%unit, "(2G14.5)") real(xmin), real(xmax)
    write(fp%unit, "(2G14.5)") real(ymin), real(ymax)
    if(present(zmin))then
       minz = zmin
    else
       call array_get_threshold(z, nx*ny, 0.99d0, minz)
    endif
    if(present(zmax))then
       maxz = zmax
    else
       call array_get_threshold(z, nx*ny, 0.01d0, maxz)
    endif

    write(fp%unit, "(2G14.5)") minz, maxz
    write(fp%unit, "(A)") "0"
    write(fp%unit, "(2I8)") nx, ny
    do i=1,nx
       write(fp%unit, "("//trim(num2str(ny))//"G14.5)") real(z(i,:))
    enddo
  end subroutine asy_density_d

  subroutine asy_density_s(fp, z, xmin, xmax, ymin, ymax, label, zmin, zmax)
    type(file_pointer) fp
    real,dimension(:,:)::z
    real xmin, xmax, ymin, ymax
    real,optional::zmin, zmax
    real minz, maxz
    integer nx, ny, i
    UNKNOWN_STRING, optional::label
    nx = size(z, 1)
    ny = size(z, 2)
    write(fp%unit, "(A)") "DENSITY"
    if(present(label))then
       if(trim(label).ne."")then
          write(fp%unit,"(A)") trim(label)
       else
          write(fp%unit,"(A)") "NULL"
       endif
    else
       write(fp%unit,"(A)") "NULL"
    endif

    write(fp%unit, "(2G14.5)") xmin, xmax
    write(fp%unit, "(2G14.5)") ymin, ymax
    if(present(zmin))then
       minz = zmin
    else
       call float_array_get_threshold(z, nx*ny, 0.99, minz)
    endif
    if(present(zmax))then
       maxz = zmax
    else
       call float_array_get_threshold(z, nx*ny, 0.01, maxz)
    endif
    write(fp%unit, "(2G14.5)") minz, maxz
    write(fp%unit, "(A)") "0"
    write(fp%unit, "(2I8)") nx, ny
    do i=1,nx
       write(fp%unit, "("//trim(num2str(ny))//"G14.5)") z(i,:)
    enddo
  end subroutine asy_density_s


 subroutine asy_irregular_density_d(fp, x, y, z, label, xmin, xmax, ymin, ymax, zmin, zmax)
    type(file_pointer) fp
    real(dl),dimension(:)::x, y, z
    UNKNOWN_STRING, optional::label
    real(dl),optional::xmin, xmax, ymin, ymax, zmin, zmax
    integer i, n
    real(dl) minx, miny, minz, maxx, maxy, maxz
    n = getdim("asy_irregular_density", size(x), size(y), size(z))
    if(present(xmin))then
       minx = xmin
    else
       call array_get_threshold(x, n, 0.99d0, minx)
    endif
    if(present(xmax))then
       maxx = xmax
    else
       call array_get_threshold(x, n, 0.01d0, maxx)
    endif
    if(present(ymin))then
       miny = ymin
    else
       call array_get_threshold(y, n, 0.99d0, miny)
    endif
    if(present(ymax))then
       maxy = ymax
    else
       call array_get_threshold(y, n, 0.01d0, maxy)
    endif
    if(present(zmin))then
       minz = zmin
    else
       call array_get_threshold(z, n, 0.99d0, minz)
    endif
    if(present(zmax))then
       maxz = zmax
    else
       call array_get_threshold(z, n, 0.01d0, maxz)
    endif
    write(fp%unit,"(A)") "DENSITY"
    if(present(label))then
       if(trim(label).ne."")then
          write(fp%unit,"(A)") trim(label)
       else
          write(fp%unit,"(A)") "NULL"
       endif
    else
       write(fp%unit,"(A)") "NULL"
    endif
    write(fp%unit, "(2G14.5)") minx, maxx
    write(fp%unit, "(2G14.5)") miny, maxy
    write(fp%unit, "(2G14.5)") minz, maxz
    write(fp%unit,"(A)") "1"
    write(fp%unit,"(I10)") n
    do i=1,n
       write(fp%unit, "(3G14.5)") real(x(i)), real(y(i)), real(z(i))
    enddo
  end subroutine asy_irregular_density_d

  subroutine asy_irregular_density_s(fp, x, y, z, label, xmin, xmax, ymin, ymax, zmin, zmax)
    type(file_pointer) fp
    real,dimension(:)::x, y, z
    UNKNOWN_STRING, optional::label
    real,optional::xmin, xmax, ymin, ymax, zmin, zmax
    integer i, n
    real minx, miny, minz, maxx, maxy, maxz
    n = getdim("asy_irregular_density", size(x), size(y), size(z))
    if(present(xmin))then
       minx = xmin
    else
       call float_array_get_threshold(x, n, 0.99, minx)
    endif
    if(present(xmax))then
       maxx = xmax
    else
       call float_array_get_threshold(x, n, 0.01, maxx)
    endif
    if(present(ymin))then
       miny = ymin
    else
       call float_array_get_threshold(y, n, 0.99, miny)
    endif
    if(present(ymax))then
       maxy = ymax
    else
       call float_array_get_threshold(y, n, 0.01, maxy)
    endif
    if(present(zmin))then
       minz = zmin
    else
       call float_array_get_threshold(z, n, 0.99, minz)
    endif
    if(present(zmax))then
       maxz = zmax
    else
       call float_array_get_threshold(z, n, 0.01, maxz)
    endif
    write(fp%unit,"(A)") "DENSITY"
    if(present(label))then
       if(trim(label).ne."")then
          write(fp%unit,"(A)") trim(label)
       else
          write(fp%unit,"(A)") "NULL"
       endif
    else
       write(fp%unit,"(A)") "NULL"
    endif
    write(fp%unit, "(2G14.5)") minx, maxx
    write(fp%unit, "(2G14.5)") miny, maxy
    write(fp%unit, "(2G14.5)") minz, maxz
    write(fp%unit,"(A)") "1"
    write(fp%unit,"(I10)") n
    do i=1,n
       write(fp%unit, "(3G14.5)") x(i), y(i), z(i)
    enddo
  end subroutine asy_irregular_density_s

  subroutine asy_clip_d(fp, x, y, smooth)
    type(file_pointer) fp
    logical,optional::smooth
    real(dl),dimension(:),intent(IN)::x,y
    integer i,n
    n = getdim("asy_clip", size(x), size(y))
    write(fp%unit, "(A)") "CLIP"
    write(fp%unit, "(I8)") n
    if(present(smooth))then
       if(smooth)then
          write(fp%unit, "(A)") "1"
       else
          write(fp%unit, "(A)") "0"
       endif
    else
       write(fp%unit, "(A)") "0"  !!no smoothing by default
    endif
    do i=1,n
       write(fp%unit,"(2G14.5)") real(x(i)), real(y(i))
    enddo

  end subroutine asy_clip_d

  subroutine asy_clip_s(fp, x, y, smooth)
    type(file_pointer) fp
    logical,optional::smooth
    real,dimension(:),intent(IN)::x,y
    integer i,n
    n = getdim("asy_clip", size(x), size(y))
    write(fp%unit, "(A)") "CLIP"
    write(fp%unit, "(I8)") n
    if(present(smooth))then
       if(smooth)then
          write(fp%unit, "(A)") "1"
       else
          write(fp%unit, "(A)") "0"
       endif
    else
       write(fp%unit, "(A)") "0"  !!no smoothing by default
    endif
    do i=1,n
       write(fp%unit,"(2G14.5)") x(i), y(i)
    enddo
  end subroutine asy_clip_s

  function asy_linestyle(color, linetype, linewidth)
    UNKNOWN_STRING color
    UNKNOWN_STRING,optional::linetype
    real,optional::linewidth
    STRING asy_linestyle
    if(trim(color).eq."" .or. trim(color).eq."NULL")then
       asy_linestyle = "NULL"
       return
    endif
    if(present(linetype))then
       asy_linestyle = trim(color)//"_"//trim(linetype)
    else
       asy_linestyle = trim(color)//"_solid"
    endif
    if(present(linewidth))then
       asy_linestyle = trim(asy_linestyle)//"_"//trim(num2str(linewidth))
    endif
  end function asy_linestyle

  subroutine asy_path_append_s(path, x, y)
    real x, y
    type(asy_path) path
    call list_push(path%l, (/ x, y /) )
  end subroutine asy_path_append_s

  subroutine asy_path_append_d(path, x, y)
    real(dl) x, y
    type(asy_path) path
    call asy_path_append_s(path, real(x), real(y))
  end subroutine asy_path_append_d

  subroutine asy_path_close(path)
    type(asy_path) path
    if(allocated(path%length))then
       path%nclosed = path%nclosed + 1
       if(path%nclosed .gt. asy_path_max_nclosed) stop "asy_path_close: too many closed paths"
       path%length(path%nclosed) = path%l%n - sum(path%length(1:path%nclosed-1))
    else
       allocate(path%length(asy_path_max_nclosed))
       path%nclosed = 1
       path%length(1) = path%l%n
    endif
  end subroutine asy_path_close

  subroutine asy_path_clear(path)
    type(asy_path) path
    if(allocated(path%length))deallocate(path%length)
    call list_initialize(path%l)
  end subroutine asy_path_clear

  !!generate a path that encloses the region of f(x, y) > threshold
  subroutine asy_path_2dcl(path, f, xmin, xmax, ymin, ymax, threshold, resolution)
    integer,parameter::default_resolution = 64
    integer,optional:: resolution
    type(asy_path) path
    external f
    real f
    real xmin, xmax, ymin, ymax, threshold, dx, dy, dxby2, dyby2
    integer n
    logical,dimension(:,:),allocatable::above
    integer,dimension(:,:,:),allocatable::lines
    integer i, j, imin, jmin, last
    call asy_path_clear(path) !!clear the path
    if(present(resolution))then
       n = resolution
    else
       n = default_resolution
    endif
    dx = (xmax-xmin)/n
    dxby2 = dx/2.
    dy = (ymax-ymin)/n
    dyby2 = dy/2.
    allocate(above(0:n+1, 0:n+1),lines(2, -1:n, -1:n))
    above(0,:) = .false.
    above(:,0) = .false.
    above(n+1,:) = .false.
    above(:, n+1) = .false.
    do j = 1, n
       do i=1,n
          above(i, j) = (f(xmin+dx*(i-0.5), ymin+dy*(j-0.5)) .gt. threshold)
       enddo
    enddo
    lines = 0
    do j=0, n
       do i=0, n
          if(above(i+1, j+1) .and. .not. above(i+1, j))then
             lines(1, i, j) = 1
          elseif(.not. above(i+1, j+1) .and. above(i+1, j))then
             lines(1, i, j)= -1
          endif
          if(above(i, j+1) .and. .not. above(i+1, j+1))then
             lines(2, i, j) =  1
          elseif(.not. above(i, j+1) .and. above(i+1, j+1))then
             lines(2, i, j) = -1
          endif
       enddo
    enddo
    jmin = 0
    imin = 0
100 continue
    do 
       if(lines(1, imin, jmin).eq.0)then
          imin = imin + 1
          if(imin .gt. n)then
             jmin = jmin+1
             if(jmin .gt. n)then
                deallocate(above, lines)
                return
             endif
             imin = 0
          endif
       else
          exit
       endif
    enddo
    j = jmin
    if(lines(1, imin, jmin) .eq. 1)then
       lines(1, imin, jmin) = 0
       i = imin
       call insert_point()
       i =  i + 1
       call insert_point()
       last = -1
    else
       lines(1, imin, jmin) = 0
       imin = imin + 1
       i = imin
       call insert_point()
       i = i - 1
       call insert_point()
       last = 1
    endif

    do 
       select case(last)
       case(1)
          if(lines(2, i, j-1) .eq. -1)then
             lines(2, i, j-1) = 0
             j = j - 1
             last = 2
          elseif(lines(1, i-1, j) .eq. -1)then
             lines(1, i-1, j) = 0
             i = i - 1             
          else
             lines(2, i, j) = 0
             j = j + 1
             last = -2
          endif
       case(-1)
          if(lines(2, i, j) .eq. 1)then
             lines(2, i, j) = 0
             j = j + 1
             last = -2
          elseif(lines(1, i, j) .eq. 1)then
             lines(1, i, j) = 0
             i = i + 1
          else
             lines(2, i, j-1) = 0
             j = j -1
             last =2
          endif
       case(2)
          if(lines(1, i, j) .eq. 1)then
             lines(1, i, j) = 0
             i = i + 1
             last  = -1
          elseif(lines(2, i, j-1).eq.-1)then
             lines(2, i, j-1) = 0
             j = j - 1
          else
             lines(1, i-1, j) = 0
             i = i - 1
             last = 1
          endif
       case default !!actually must be -2
          if(lines(1, i-1, j).eq. -1)then
             lines(1, i-1, j) = 0
             i = i-1
             last = 1
          elseif(lines(2, i, j) .eq. 1)then
             lines(2, i, j) = 0
             j = j + 1
          else 
             lines(1, i, j) = 0
             i = i + 1
             last = -1
          endif
       end select
       if(i.eq.imin .and. j.eq.jmin)then !!this happens in finite steps
          call asy_path_close(path)
          goto 100 
       else
          call insert_point()
       endif
    enddo

  contains

    subroutine insert_point()
      real x, y
      real fp(2), s
      x =  xmin+i*dx
      y = ymin+j*dy
      if(i.eq.0 .or. i.eq.n)then
         fp(1) = 0 
      else
         fp(1) = (f(x+dxby2, y) - f(x-dxby2, y))
      endif
      if(j.eq.0 .or. j.eq.n)then
         fp(2) = 0
      else
         fp(2) =  (f(x, y+dyby2) - f(x, y-dyby2))
      endif
      s= sum(fp**2)
      if(s.gt.0)then
         s = (threshold - f(x, y))/s
         x = min(xmax, max(x + max(min(fp(1)*s, 0.5), -0.5)*dx, xmin))
         y = min(ymax, max(y + max(min(fp(2)*s, 0.5), -0.5)*dy, ymin))
      endif
      call asy_path_append(path, x, y)
    end subroutine insert_point

  end subroutine asy_path_2dcl


  subroutine asy_path_from_array_center(path, f, xmin, xmax, ymin, ymax, threshold)
    type(asy_path) path
    real f(:,:)
    real xmin, xmax, ymin, ymax, threshold, dx, dy
    integer nx, ny, n
    nx = size(f, 1)
    ny = size(f, 2)
    n = min(max(nx, ny,16)*2, 256) !!higher resolution is useless since we do bilinear interpolation
    dx = (xmax-xmin)/nx
    dy = (ymax- ymin)/ny
    call asy_path_2dcl(path, finterp, xmin, xmax, ymin, ymax, threshold, n)
  contains
    
    function  finterp(x, y)
      real x, y, finterp
      integer i, j, ip1, jp1
      real ri, rj
      ri = (x-xmin)/dx + 0.5
      rj = (y-ymin)/dy + 0.5
      i = floor(ri)
      j = floor(rj)
      ri = ri - i
      rj = rj - j
      ip1 = i + 1
      jp1 = j + 1
      i = min(max(i,1), nx)
      ip1 = min(max(ip1,1), nx)
      j = min(max(j, 1), ny)
      jp1 = min(max(jp1, 1), ny)
      finterp = (f(i, j) * (1.-ri) + f(ip1, j)*ri)*(1.-rj)+(f(i, jp1)*(1.-ri) + f(ip1, jp1)*ri)*rj
    end function finterp
    
  end subroutine asy_path_from_array_center



  subroutine asy_path_from_array(path, f, xmin, xmax, ymin, ymax, threshold)
    type(asy_path) path
    real f(:,:)
    real xmin, xmax, ymin, ymax, threshold, dx, dy
    integer nx, ny, n
    nx = size(f, 1)
    ny = size(f, 2)
    n = min(max(nx, ny,16)*2, 256) !!higher resolution is useless since we do bilinear interpolation
    dx = (xmax-xmin)/(nx-1)
    dy = (ymax- ymin)/(ny-1)
    call asy_path_2dcl(path, finterp, xmin, xmax, ymin, ymax, threshold, n)
  contains
    
    function  finterp(x, y)
      real x, y, finterp
      integer i, j
      real ri, rj
      ri = (x-xmin)/dx + 1.
      rj = (y-ymin)/dy + 1.
      i = min(max(floor(ri), 1), nx-1)
      j = min(max(floor(rj), 1), ny-1)
      ri = ri - i
      rj = rj - j
      finterp = (f(i, j)*(1.-ri) + f(i+1, j)*ri)*(1.-rj)+(f(i, j+1)*(1.-ri) + f(i+1, j+1)*ri)*rj
    end function finterp
    
  end subroutine asy_path_from_array


  subroutine asy_legend_default(fp)
    type(file_pointer) fp
    call asy_legend_location(fp, "N")
  end subroutine asy_legend_default

  subroutine asy_legend_location(fp, loc, cols)
    type(file_pointer) fp
    UNKNOWN_STRING :: loc
    integer, optional::cols
    write(fp%unit, "(A)") "LEGEND"
    if(trim(loc).ne."")then
       write(fp%unit, "(A)") trim(loc)
    else
       write(fp%unit, "(A)") "E"
    end if
    if(present(cols))then
       write(fp%unit, "(I8)")  cols
    else
       write(fp%unit, "(I8)")  1
    endif
  end subroutine asy_legend_location

  subroutine asy_legend_d(fp, x, y, cols)
    type(file_pointer) fp
    real(dl)x, y
    integer,optional::cols
    write(fp%unit, "(A)") "LEGEND"
    write(fp%unit, "(A)") "NULL"
    write(fp%unit, "(2G15.4)") x, y
    if(present(cols))then
       write(fp%unit, "(I8)")  cols
    else
       write(fp%unit, "(I8)")  1
    endif
  end subroutine asy_legend_d

  subroutine asy_legend_s(fp, x, y, cols)
    type(file_pointer) fp
    real x, y
    integer,optional::cols
    write(fp%unit, "(A)") "LEGEND"
    write(fp%unit, "(A)") "NULL"
    write(fp%unit, "(2G15.4)") x, y
    if(present(cols))then
       write(fp%unit, "(I8)")  cols
    else
       write(fp%unit, "(I8)")  1
    endif
  end subroutine asy_legend_s


  subroutine asy_topaxis_s(fp, xmin, xmax, islog,label)
    real xmin, xmax
    logical islog
    type(file_pointer) fp
    UNKNOWN_STRING label
    write(fp%unit, "(A)")"EXTRA_AXIS"
    write(fp%unit, "(A)") "top"
    if(trim(label).ne."")then
       write(fp%unit, "(A)") trim(label)
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(islog)then
       write(fp%unit, "(A)") "1"
    else
       write(fp%unit, "(A)") "0"
    endif
    write(fp%unit, "(2G14.5)") xmin, xmax
  end subroutine asy_topaxis_s

  subroutine asy_topaxis_d(fp, xmin, xmax, islog,label)
    real(dl) xmin, xmax
    logical islog
    type(file_pointer) fp
    UNKNOWN_STRING label
    call asy_topaxis_s(fp, real(xmin), real(xmax), islog,label)
  end subroutine asy_topaxis_d


  subroutine asy_rightaxis_s(fp, ymin, ymax, islog, label)
    real ymin, ymax
    logical islog
    type(file_pointer) fp
    UNKNOWN_STRING label
    write(fp%unit, "(A)")"EXTRA_AXIS"
    write(fp%unit, "(A)") "right"
    if(trim(label).ne."")then
       write(fp%unit, "(A)") trim(label)
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(islog)then
       write(fp%unit, "(A)") "1"
    else
       write(fp%unit, "(A)") "0"
    endif
    write(fp%unit, "(2G14.5)") ymin, ymax
  end subroutine asy_rightaxis_s

  subroutine asy_rightaxis_d(fp, ymin, ymax, islog, label)
    real(dl) ymin, ymax
    logical islog
    UNKNOWN_STRING label
    type(file_pointer) fp
    call asy_rightaxis_s(fp, real(ymin), real(ymax), islog,label)
  end subroutine asy_rightaxis_d

  subroutine asy_error_bar_d(fp, x, y, dy_minus, dy_plus, dx_minus, dx_plus, color)
    type(file_pointer) fp
    real(dl) x, y
    real(dl),optional::dy_minus, dy_plus, dx_minus, dx_plus
    UNKNOWN_STRING, optional::color
    if(present(color))then
       call asy_label(fp, "$\bullet$", x, y, color)
       if(present(dy_minus))then
          call asy_line(fp, x, y, x, y-dy_minus, linewidth = 1., color=color)
          call asy_label(fp, "-", x, y-dy_minus, color=color)
       endif
       if(present(dy_plus))then
          call asy_line(fp, x, y, x, y+dy_plus, linewidth = 1., color= color)
          call asy_label(fp, "-", x, y+dy_plus, color=color)
       endif
       if(present(dx_minus))then
          call asy_line(fp, x, y, x-dx_minus, y, linewidth = 1., color=color)
          call asy_label(fp, "{\tiny $|$}", x-dx_minus, y, color=color)
       endif
       if(present(dx_plus))then
          call asy_line(fp, x, y, x+dx_plus, y, linewidth = 1., color=color)
          call asy_label(fp, "{\tiny $|$}", x+dx_plus, y, color=color)
       endif
    else
       call asy_label(fp, "$\bullet$", x, y)
       if(present(dy_minus))then
          call asy_line(fp, x, y, x, y-dy_minus, linewidth = 1.)
          call asy_label(fp, "-", x, y-dy_minus)
       endif
       if(present(dy_plus))then
          call asy_line(fp, x, y, x, y+dy_plus, linewidth = 1.)
          call asy_label(fp, "-", x, y+dy_plus)
       endif
       if(present(dx_minus))then
          call asy_line(fp, x, y, x-dx_minus, y, linewidth = 1.)
          call asy_label(fp, "{\tiny $|$}", x-dx_minus, y)
       endif
       if(present(dx_plus))then
          call asy_line(fp, x, y, x+dx_plus, y, linewidth = 1.)
          call asy_label(fp, "{\tiny $|$}", x+dx_plus, y)
       endif
    endif
  end subroutine asy_error_bar_d

  subroutine asy_error_bar_s(fp, x, y, dy_minus, dy_plus, dx_minus, dx_plus, color)
    type(file_pointer) fp
    real x, y
    real,optional::dy_minus, dy_plus, dx_minus, dx_plus
    UNKNOWN_STRING, optional::color
    if(present(color))then
       call asy_label(fp, "$\bullet$", x, y, color)
       if(present(dy_minus))then
          call asy_line(fp, x, y, x, y-dy_minus, linewidth = 1., color=color)
          call asy_label(fp, "-", x, y-dy_minus, color=color)
       endif
       if(present(dy_plus))then
          call asy_line(fp, x, y, x, y+dy_plus, linewidth = 1., color= color)
          call asy_label(fp, "-", x, y+dy_plus, color=color)
       endif
       if(present(dx_minus))then
          call asy_line(fp, x, y, x-dx_minus, y, linewidth = 1., color=color)
          call asy_label(fp, "{\tiny $|$}", x-dx_minus, y, color=color)
       endif
       if(present(dx_plus))then
          call asy_line(fp, x, y, x+dx_plus, y, linewidth = 1., color=color)
          call asy_label(fp, "{\tiny $|$}", x+dx_plus, y, color=color)
       endif
    else
       call asy_label(fp, "$\bullet$", x, y)
       if(present(dy_minus))then
          call asy_line(fp, x, y, x, y-dy_minus, linewidth = 1.)
          call asy_label(fp, "-", x, y-dy_minus)
       endif
       if(present(dy_plus))then
          call asy_line(fp, x, y, x, y+dy_plus, linewidth = 1.)
          call asy_label(fp, "-", x, y+dy_plus)
       endif
       if(present(dx_minus))then
          call asy_line(fp, x, y, x-dx_minus, y, linewidth = 1.)
          call asy_label(fp, "{\tiny $|$}", x-dx_minus, y)
       endif
       if(present(dx_plus))then
          call asy_line(fp, x, y, x+dx_plus, y, linewidth = 1.)
          call asy_label(fp, "{\tiny $|$}", x+dx_plus, y)
       endif
    endif
  end subroutine asy_error_bar_s


  function asy_rgb_color_s(r, g, b) result(color)
    real r, g, b
    SHORT_STRING color
    color = "RGB:"//trim(num2str(nint(min(1., max(0., r))*255.*8.)/8.))//":"//trim(num2str(nint(min(1., max(0., g))*255.*8.)/8.))//":"//trim(num2str(nint(min(1., max(0., b))*255.*8.)/8.))
  end function asy_rgb_color_s

  function asy_rgb_color_d(r,g,b) result(color)
    real(dl) r, g, b
    SHORT_STRING color
    color = asy_rgb_color_s(real(r), real(g), real(b))
  end function asy_rgb_color_d


  function asy_gray_color_s(gray) result(color)
    real gray
    SHORT_STRING color
    color = "GRAY:"//trim(num2str(nint(min(1., max(0., gray))*255.*8.)/8.))
  end function asy_gray_color_s


  function asy_gray_color_d(gray) result(color)
    real(dl) gray
    SHORT_STRING color
    color = asy_gray_color_s(real(gray))
  end function asy_gray_color_d

end module asy_utils
