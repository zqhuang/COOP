module coop_asy_mod
  use coop_wrapper_typedef
  use coop_interpolation_mod
  use coop_file_mod
  use coop_list_mod
  implicit none
  private

  public::coop_asy, coop_asy_error_bar, coop_asy_interpolate_curve, coop_asy_gray_color, coop_asy_rgb_color, coop_asy_label, coop_asy_legend, coop_asy_dot, coop_asy_line, coop_asy_labels, coop_asy_dots, coop_asy_lines, coop_asy_contour, coop_asy_curve, coop_asy_density,  coop_asy_topaxis, coop_asy_rightaxis, coop_asy_clip


#include "constants.h"

  integer,parameter::dl = kind(1.d0)
  integer,parameter::sp = kind(1.)
  real(sp),parameter::coop_asy_default_width = 6.6
  real(sp),parameter::coop_asy_default_height = 5.5

  type, extends(coop_file) :: coop_asy
     real(sp) xmin, xmax, ymin, ymax, width, height
   contains
     procedure::init => coop_asy_init
     procedure::write_coor => coop_asy_write_coor
     procedure::write_limits => coop_asy_write_limits
  end type coop_asy



  COOP_INT , parameter::coop_asy_path_max_nclosed = 4096

  interface coop_asy_error_bar
     module procedure coop_asy_error_bar_s, coop_asy_error_bar_d
  end interface coop_asy_error_bar

  interface coop_asy_interpolate_curve
     module procedure coop_asy_interpolate_curve_d, coop_asy_interpolate_curve_s
  end interface coop_asy_interpolate_curve

  interface coop_asy_gray_color
     module procedure coop_asy_gray_color_s, coop_asy_gray_color_d
  end interface coop_asy_gray_color

  interface coop_asy_rgb_color
     module procedure coop_asy_rgb_color_s, coop_asy_rgb_color_d
  end interface coop_asy_rgb_color

  interface coop_asy_label
     module procedure coop_asy_label_s, coop_asy_label_d
  end interface coop_asy_label

  interface coop_asy_legend
     module procedure coop_asy_legend_s, coop_asy_legend_d, coop_asy_legend_default, coop_asy_legend_location
  end interface coop_asy_legend

  interface coop_asy_dot
     module procedure coop_asy_dot_s, coop_asy_dot_d
  end interface coop_asy_dot

  interface coop_asy_clip
     module procedure coop_asy_clip_s, coop_asy_clip_d
  end interface coop_asy_clip


  interface coop_asy_line
     module procedure coop_asy_line_s, coop_asy_line_d
  end interface coop_asy_line

  interface coop_asy_labels
     module procedure coop_asy_labels_s, coop_asy_labels_d
  end interface coop_asy_labels

  interface coop_asy_dots
     module procedure coop_asy_dots_s, coop_asy_dots_d
  end interface coop_asy_dots

  interface coop_asy_lines
     module procedure coop_asy_lines_s, coop_asy_lines_d
  end interface coop_asy_lines

  interface coop_asy_contour
     module procedure coop_asy_contour_s, coop_asy_contour_d, coop_asy_contour_mult_s, coop_asy_contour_mult_d,  coop_asy_contour_arr_d, coop_asy_contour_arr_s, coop_asy_contour_path
  end interface coop_asy_contour

  interface coop_asy_curve
     module procedure coop_asy_curve_s, coop_asy_curve_d
  end interface coop_asy_curve

  interface coop_asy_density
     module procedure coop_asy_density_s, coop_asy_density_d, coop_asy_irregular_density_s, coop_asy_irregular_density_d
  end interface coop_asy_density

  interface coop_asy_path_append
     module procedure coop_asy_path_append_s,  coop_asy_path_append_d
  end interface coop_asy_path_append

  interface coop_asy_topaxis
     module procedure coop_asy_topaxis_s, coop_asy_topaxis_d
  end interface coop_asy_topaxis

  interface coop_asy_rightaxis
     module procedure coop_asy_rightaxis_s, coop_asy_rightaxis_d
  end interface coop_asy_rightaxis

  type coop_asy_path
     COOP_INT  nclosed
     COOP_INT ,dimension(:),allocatable:: length
     type(coop_list_realarr) l
  end type coop_asy_path



contains


  subroutine coop_asy_init(fp,  xmin, xmax, ymin, ymax, width, height, caption, xlabel, ylabel, xlog, ylog, zlog, doclip, nblocks)
    class(coop_asy) fp
    real(sp),optional:: width, height
    COOP_UNKNOWN_STRING,optional:: caption, xlabel, ylabel
    COOP_INT ,optional::nblocks
    logical,optional::xlog, ylog, zlog, doclip
    real(sp),optional:: xmin, xmax, ymin, ymax
    character(len = 5) tmp
    if(present(width))then
       fp%width = width
    else
       fp%width = coop_asy_default_width
    endif
    if(present(height))then
       fp%height = height
    else
       fp%height = coop_asy_default_height
    endif
    write(fp%unit, "(2F12.2)") fp%width, fp%height
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
       fp%xmin = xmin
    else
       fp%xmin = 1.1e30
    endif
    if(present(xmax))then
       fp%xmax = xmax
    else
       fp%xmax = -1.1e30
    endif
    if(present(ymin))then
       fp%ymin = ymin
    else
       fp%ymin = 1.1e30
    endif
    if(present(ymax))then
       fp%ymax = ymax
    else
       fp%ymax = -1.1e30
    endif
    write(fp%unit,"(2G14.5)") fp%xmin, fp%xmax
    write(fp%unit,"(2G14.5)") fp%ymin, fp%ymax
    if(present(nblocks))then
       write(fp%unit,"(I5)") nblocks    !!plot only n blocks
    else
       write(fp%unit,"(I5)") 0   !!0 means any number of blocks
    endif
  end subroutine coop_asy_init

  subroutine coop_asy_write_limits(fp, xmin, xmax, ymin, ymax)
    class(coop_asy) fp
    real(sp) xmin, xmax, ymin, ymax
    write(fp%unit, "(2G14.5)") xmin, xmax
    write(fp%unit, "(2G14.5)") ymin, ymax
    fp%xmin = min(xmin, fp%xmin)
    fp%xmax = min(xmax, fp%xmax)
    fp%ymin = min(ymin, fp%ymin)
    fp%ymax = min(ymax, fp%ymax)
  end subroutine coop_asy_write_limits

  subroutine coop_asy_write_coor(fp, x, y, x2, y2)
    class(coop_asy) fp
    real(sp) x, y
    real(sp),optional:: x2, y2
    if(present(x2) .and. present(y2))then
       write(fp%unit, "(4G14.5)") x, y, x2, y2
       fp%xmin = min(x2, fp%xmin)
       fp%xmax = max(x2, fp%xmax)
       fp%ymin = min(y2, fp%ymin)
       fp%ymax = max(y2, fp%ymax)
    else
       write(fp%unit, "(2G14.5)") x, y
    endif
    fp%xmin = min(x, fp%xmin)
    fp%xmax = max(x, fp%xmax)
    fp%ymin = min(y, fp%ymin)
    fp%ymax = max(y, fp%ymax)
  end subroutine coop_asy_write_coor

  subroutine coop_asy_dot_d(fp, x, y, color, symbol)
    class(coop_asy) fp
    COOP_INT  n
    real(dl) x,y
    COOP_UNKNOWN_STRING,optional:: color
    COOP_UNKNOWN_STRING,optional::symbol
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
    call fp%write_coor(real(x, sp), real(y, sp))
  end subroutine coop_asy_dot_d


  subroutine coop_asy_dot_s(fp, x, y, color, symbol)
    class(coop_asy) fp
    COOP_INT  n
    real(sp) x,y
    COOP_UNKNOWN_STRING,optional:: color
    COOP_UNKNOWN_STRING,optional::symbol
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
    call fp%write_coor(x, y)
  end subroutine coop_asy_dot_s



  subroutine coop_asy_dots_d(fp, x, y, color, symbol)
    class(coop_asy) fp
    COOP_INT  n,i
    real(dl),dimension(:),intent(IN)::x,y
    COOP_UNKNOWN_STRING,optional:: color
    COOP_UNKNOWN_STRING,optional::symbol
    n = coop_getdim( "coop_asy_dot_block", size(x), size(y))
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
       call fp%write_coor(real(x(i),sp), real(y(i),sp))
    enddo
  end subroutine coop_asy_dots_d

  subroutine coop_asy_dots_s(fp, x, y, color, symbol)
    class(coop_asy) fp
    COOP_INT  n,i
    real(sp),dimension(:),intent(IN)::x,y
    COOP_UNKNOWN_STRING,optional:: color
    COOP_UNKNOWN_STRING,optional::symbol
    n = coop_getdim( "coop_asy_dot_block", size(x), size(y))
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
       call fp%write_coor(x(i),y(i))
    enddo
  end subroutine coop_asy_dots_s


  subroutine coop_asy_line_d(fp, xstart, ystart, xend, yend, color, linetype, linewidth)
    class(coop_asy) fp
    real(dl)::xstart, ystart, xend, yend
    COOP_UNKNOWN_STRING, optional::color, linetype
    real(sp),optional::linewidth
    COOP_STRING lineproperty
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
       lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    call fp%write_coor(real(xstart, sp), real(ystart, sp), real(xend, sp), real(yend, sp))
  end subroutine coop_asy_line_d


  subroutine coop_asy_line_s(fp, xstart, ystart, xend, yend, color, linetype, linewidth)
    class(coop_asy) fp
    real::xstart, ystart, xend, yend
    COOP_UNKNOWN_STRING, optional::color, linetype
    real(sp),optional::linewidth
    COOP_STRING lineproperty
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
       lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    call fp%write_coor(xstart, ystart, xend, yend)
  end subroutine coop_asy_line_s
  
  subroutine coop_asy_lines_d(fp, xstart, ystart, xend, yend, color, linetype, linewidth)
    class(coop_asy) fp
    real(dl),dimension(:),intent(IN)::xstart, ystart, xend, yend
    COOP_UNKNOWN_STRING, optional::color, linetype
    real(sp),optional::linewidth
    COOP_STRING lineproperty
    COOP_INT  n, i
    write(fp%unit, "(A)") "LINES"
    n = coop_getdim("coop_asy_lines", size(xstart), size(ystart), size(xend), size(yend))
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
       lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    do i = 1, n
       call fp%write_coor(real(xstart(i),sp), real(ystart(i),sp), real(xend(i), sp), real(yend(i), sp))
    enddo
  end subroutine coop_asy_lines_d


  subroutine coop_asy_lines_s(fp, xstart, ystart, xend, yend, color, linetype, linewidth)
    class(coop_asy) fp
    real(sp),dimension(:),intent(IN)::xstart, ystart, xend, yend
    COOP_UNKNOWN_STRING, optional::color, linetype
    real(sp),optional::linewidth
    COOP_STRING lineproperty
    COOP_INT  n, i
    write(fp%unit, "(A)") "LINES"
    n = coop_getdim("coop_asy_lines", size(xstart), size(ystart), size(xend), size(yend))
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
       lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    do i = 1, n
       call fp%write_coor( xstart(i), ystart(i), xend(i), yend(i))
    enddo
  end subroutine coop_asy_lines_s


  subroutine coop_asy_interpolate_curve_d(fp, xraw, yraw, interpolate, color, linetype, linewidth, legend)
    COOP_INT ,parameter::n = 256
    real(dl) x(n), y(n),  minx, maxx, dx
    COOP_INT  w(n)
    real(dl),dimension(:),allocatable::xx, yy, yy2
    class(coop_asy) fp
    real(dl),dimension(:),intent(IN)::xraw, yraw
    COOP_UNKNOWN_STRING interpolate
    COOP_UNKNOWN_STRING,optional:: color, linetype
    real(sp),optional::linewidth
    COOP_UNKNOWN_STRING, optional::legend
    COOP_INT  i, j, m, nn, npt
    COOP_STRING lineproperty
    logical do_raw
    m = Coop_getdim("coop_asy_interpolate_curve", size(xraw), size(yraw))
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
       call coop_set_uniform(n, x, minx, maxx)
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
          call coop_spline(m, xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call coop_splint(m, xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
    case("LinearLog")
       call coop_set_uniform(n, x, minx, maxx)
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
          call coop_spline(m, xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call coop_splint(m, xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
       y = exp(y)
    case("LogLinear")
       minx = log(minx)
       maxx = log(maxx)
       call coop_set_uniform(n, x, minx, maxx)
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
          call coop_spline(m, xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call coop_splint(m, xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
       x = exp(x)
    case("LogLog")
       minx = log(minx)
       maxx = log(maxx)
       call coop_set_uniform(n, x, minx, maxx)
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
          call coop_spline(m, xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call coop_splint(m, xx, yy, yy2, x(i), y(i))
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
       lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    write(fp%unit, "(A)") "0"
    if(do_raw)then
       do i=1,npt
          call fp%write_coor( real(xraw(i), sp), real(yraw(i), sp) )
       enddo
    else
       do i=1,npt
          call fp%write_coor(real(x(i), sp), real(y(i), sp))
       enddo
    endif
  end subroutine coop_asy_interpolate_curve_d

  subroutine coop_asy_interpolate_curve_s(fp, xraw, yraw, interpolate, color, linetype, linewidth, legend)
    COOP_INT ,parameter::n = 256
    real(dl) x(n), y(n),  minx, maxx, dx
    COOP_INT  w(n)
    real(dl),dimension(:),allocatable::xx, yy, yy2
    class(coop_asy) fp
    real(sp),dimension(:),intent(IN)::xraw, yraw
    COOP_UNKNOWN_STRING interpolate
    COOP_UNKNOWN_STRING,optional:: color, linetype
    real(sp),optional::linewidth
    COOP_UNKNOWN_STRING, optional::legend
    COOP_INT  i, j, m, nn, npt
    COOP_STRING lineproperty
    logical do_raw
    m = Coop_getdim("coop_asy_interpolate_curve", size(xraw), size(yraw))
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
       call coop_set_uniform(n, x, minx, maxx)
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
          call coop_spline(m, xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call coop_splint(m, xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
    case("LinearLog")
       call coop_set_uniform(n, x, minx, maxx)
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
          call coop_spline(m, xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call coop_splint(m,xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
       y = exp(y)
    case("LogLinear")
       minx = log(minx)
       maxx = log(maxx)
       call coop_set_uniform(n, x, minx, maxx)
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
          call coop_spline(m, xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call coop_splint(m, xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
       x = exp(x)
    case("LogLog")
       minx = log(minx)
       maxx = log(maxx)
       call coop_set_uniform(n, x, minx, maxx)
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
          call coop_spline(m, xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call coop_splint(m, xx, yy, yy2, x(i), y(i))
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
       lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    write(fp%unit, "(A)") "0"
    if(do_raw)then
       do i=1,npt
          call fp%write_coor( xraw(i), yraw(i) )
       enddo
    else
       do i=1,npt
          call fp%write_coor( real(x(i), sp), real(y(i), sp) )
       enddo
    endif
  end subroutine coop_asy_interpolate_curve_s


  subroutine coop_asy_curve_from_file(fp, fname, interpolate, xcol, ycol, color, linetype, linewidth, legend)
    COOP_INT ,parameter::n = 256
    COOP_INT , optional::xcol, ycol
    real(dl) x(n), y(n),  minx, maxx, dx
    COOP_INT  w(n)
    real(dl),dimension(:),allocatable::xx, yy, yy2
    class(coop_asy) fp
    type(coop_file) fp2
    real(dl),dimension(:),allocatable::xraw, yraw, line
    COOP_UNKNOWN_STRING interpolate, fname
    COOP_UNKNOWN_STRING,optional:: color, linetype
    real(sp),optional::linewidth
    COOP_UNKNOWN_STRING, optional::legend
    COOP_INT  i, j, m, nn, npt, ncols, xl, yl
    COOP_STRING lineproperty
    logical do_raw
    if(.not. coop_file_exists(fname))then
       write(*,*) "coop_asy_curve_from_file: the data file "//trim(fname)//" does not exist"
       stop 
    endif
    m = coop_file_numlines(fname)
    ncols = coop_file_numColumns(fname)
    if(ncols .eq. 0)then
       write(*,*) "coop_asy_curve_from_file: the data file "//trim(fname)//" is empty"
       stop 
    endif
    if(present(xcol))then
       if(xcol .gt. ncols) stop "coop_asy_curve_from_file: xcol > ncol"
       xl = xcol
    else
       x = min(ncols-1, 1)
    endif
    if(present(ycol))then
       if(ycol .gt. ncols) stop "coop_asy_curve_from_file: ycol > ncol"
       yl = ycol
    else
       yl = min(ncols, 2)
    endif
    allocate(line(0:ncols))
    allocate(xraw(m), yraw(m))
    call fp2%open(trim(fname), "r")
    do i = 1, m
       if(fp2%read_real_array(line(1:ncols)))then
          line(0) = i
          xraw(i) = line(xl)
          yraw(i) = line(yl)
       else
          write(*,*) trim(fname)//" has a wrong format"
          stop
       endif
    enddo
    call fp2%close()
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
       call coop_set_uniform(n, x, minx, maxx)
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
          call coop_spline(m, xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call coop_splint(m, xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
    case("LinearLog")
       call coop_set_uniform(n, x, minx, maxx)
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
          call coop_spline(m, xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call coop_splint(m, xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
       y = exp(y)
    case("LogLinear")
       minx = log(minx)
       maxx = log(maxx)
       call coop_set_uniform(n, x, minx, maxx)
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
          call coop_spline(m, xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call coop_splint(m, xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
       x = exp(x)
    case("LogLog")
       minx = log(minx)
       maxx = log(maxx)
       call coop_set_uniform(n, x, minx, maxx)
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
          call coop_spline(m, xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call coop_splint(m, xx, yy, yy2, x(i), y(i))
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
       lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    write(fp%unit, "(A)") "0"
    if(do_raw)then
       do i=1,npt
          call fp%write_coor( real(xraw(i), sp), real(yraw(i), sp))
       enddo
    else
       do i=1,npt
          call fp%write_coor(real(x(i), sp), real(y(i), sp))
       enddo
    endif
    deallocate(xraw, yraw, line)
  end subroutine coop_asy_curve_from_file


  subroutine coop_asy_curve_d(fp, x, y, smooth, color, linetype, linewidth, legend)
    class(coop_asy) fp
    logical,optional::smooth
    real(dl),dimension(:),intent(IN)::x,y
    COOP_UNKNOWN_STRING,optional:: color, linetype
    real(sp),optional::linewidth
    COOP_UNKNOWN_STRING, optional::legend
    COOP_INT  i,n
    COOP_STRING lineproperty
    n = coop_getdim("coop_asy_curve", size(x), size(y))
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
       lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
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
       call fp%write_coor(real(x(i), sp), real(y(i), sp))
    enddo
  end subroutine coop_asy_curve_d

  subroutine coop_asy_curve_s(fp, x, y, smooth, color, linetype, linewidth, legend)
    class(coop_asy) fp
    logical,optional::smooth
    real(sp),dimension(:),intent(IN)::x,y
    COOP_UNKNOWN_STRING,optional:: color, linetype
    real(sp),optional::linewidth
    COOP_INT  i,n
    COOP_UNKNOWN_STRING,optional::legend
    COOP_STRING lineproperty
    n = coop_getdim("coop_asy_curve", size(x), size(y))
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
       lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
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
       call fp%write_coor( x(i), y(i) )
    enddo
  end subroutine coop_asy_curve_s

  subroutine coop_asy_labels_d(fp, labels, x, y, color)
    COOP_STRING, dimension(:),intent(IN)::labels
    real(dl),dimension(:),intent(IN)::x, y
    class(coop_asy) fp
    COOP_UNKNOWN_STRING,optional::color
    COOP_INT  n, i
    write(fp%unit, "(A)") "LABELS"
    n = coop_getdim("coop_asy_labels", size(x), size(y), size(labels))
    write(fp%unit,"(I8)") n
    if(present(color))then
       write(fp%unit, "(A)") trim(color)
    else
       write(fp%unit, "(A)") "black"
    endif
    do i=1,n
       call fp%write_coor(real(x(i), sp), real(y(i), sp))
       if(trim(labels(i)).eq."")then
          write(fp%unit, "(A)") "NULLL"
       else
          write(fp%unit, "(A)") trim(labels(i))
       endif
    enddo
  end subroutine coop_asy_labels_d

  subroutine coop_asy_labels_s(fp, labels, x, y, color)
    COOP_STRING, dimension(:),intent(IN)::labels
    real(sp),dimension(:),intent(IN)::x, y
    class(coop_asy) fp
    COOP_UNKNOWN_STRING,optional::color
    COOP_INT  n, i
    write(fp%unit, "(A)") "LABELS"
    n = coop_getdim("coop_asy_labels", size(x), size(y), size(labels))
    write(fp%unit,"(I8)") n
    if(present(color))then
       write(fp%unit, "(A)") trim(color)
    else
       write(fp%unit, "(A)") "black"
    endif
    do i=1,n
       call fp%write_coor( x(i), y(i))
       if(trim(labels(i)).eq."")then
          write(fp%unit, "(A)") "NULLL"
       else
          write(fp%unit, "(A)") trim(labels(i))
       endif
    enddo
  end subroutine coop_asy_labels_s


  subroutine coop_asy_label_d(fp, label, x, y, color)
    class(coop_asy) fp
    COOP_UNKNOWN_STRING label
    COOP_UNKNOWN_STRING,optional::color
    real(dl) x, y
    write(fp%unit, "(A)") "LABELS"
    write(fp%unit, "(A)") "1"
    if(present(color))then
       write(fp%unit, "(A)") trim(color)
    else
       write(fp%unit, "(A)") "black"
    endif
    call fp%write_coor( real(x, sp), real(y, sp) )
    if(trim(label).eq."")then
       write(fp%unit, "(A)") "NULL"
    else
       write(fp%unit, "(A)") trim(label)
    endif
  end subroutine coop_asy_label_d


  subroutine coop_asy_label_s(fp, label, x, y, color)
    class(coop_asy) fp
    COOP_UNKNOWN_STRING label
    COOP_UNKNOWN_STRING,optional::color
    real(sp) x, y
    write(fp%unit, "(A)") "LABELS"
    write(fp%unit, "(A)") "1"
    if(present(color))then
       write(fp%unit, "(A)") trim(color)
    else
       write(fp%unit, "(A)") "black"
    endif
    call fp%write_coor( x, y)
    if(trim(label).eq."")then
       write(fp%unit, "(A)") "NULL"
    else
       write(fp%unit, "(A)") trim(label)
    endif
  end subroutine coop_asy_label_s



  subroutine coop_asy_contour_d(fp, x, y, colorfill, smooth, linecolor, linetype, linewidth)
    class(coop_asy) fp
    logical,optional::smooth
    real(dl),dimension(:),intent(IN)::x,y
    COOP_UNKNOWN_STRING,optional::colorfill
    COOP_UNKNOWN_STRING,optional:: linecolor, linetype
    real(sp),optional::linewidth
    COOP_INT  i,n
    COOP_STRING lineproperty
    n = coop_getdim("coop_asy_contour", size(x), size(y))
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
          lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
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
       call fp%write_coor(real(x(i), sp), real(y(i), sp))
    enddo
  end subroutine coop_asy_contour_d

  subroutine coop_asy_contour_s(fp, x, y, colorfill, smooth, linecolor, linetype, linewidth)
    class(coop_asy) fp
    logical,optional::smooth
    real(sp),dimension(:),intent(IN)::x,y
    COOP_UNKNOWN_STRING,optional::colorfill
    COOP_UNKNOWN_STRING,optional:: linecolor, linetype
    real(sp),optional::linewidth
    COOP_INT  i,n
    COOP_STRING lineproperty
    n = coop_getdim("coop_asy_contour", size(x), size(y))
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
          lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
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
       call fp%write_coor( x(i), y(i) )
    enddo
  end subroutine coop_asy_contour_s


  subroutine coop_asy_contour_path(fp, path, colorfill, smooth, linecolor, linetype, linewidth)
    class(coop_asy) fp
    logical,optional::smooth
    type(coop_asy_path) path
    real(sp) xy(2)
    COOP_UNKNOWN_STRING,optional::colorfill
    COOP_UNKNOWN_STRING,optional:: linecolor, linetype
    real(sp),optional::linewidth
    COOP_INT  i, ipath, pl
    COOP_STRING lineproperty
    if(path%nclosed.eq.0 .or. (.not. path%l%isinit()) .or. path%l%n .eq. 0)return
    if(path%l%dim .ne.2) stop "coop_asy_contour_path: wrong dimension of list in path"
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
          lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
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
          call path%l%get_element(i, xy)
          call fp%write_coor(xy(1), xy(2))
          i = i + 1
       enddo
    enddo
  end subroutine coop_asy_contour_path

  subroutine coop_asy_contour_mult_d(fp, x, y, colorfill, smooth, linecolor, linetype, linewidth)
    class(coop_asy) fp
    logical,optional::smooth
    real(dl),dimension(:,:),intent(IN)::x,y
    COOP_UNKNOWN_STRING,optional::colorfill
    COOP_UNKNOWN_STRING,optional:: linecolor, linetype
    real(sp),optional::linewidth
    COOP_INT  i,n, m,ipath
    COOP_STRING lineproperty
    n = coop_getdim("coop_asy_contour", size(x,1), size(y,1))
    m = coop_getdim("coop_asy_contour", size(x,2), size(y,2))
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
          lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
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
          call fp%write_coor(real(x(i, ipath), sp), real(y(i, ipath), sp))
       enddo
    enddo
  end subroutine coop_asy_contour_mult_d

  subroutine coop_asy_contour_mult_s(fp, x, y, colorfill, smooth, linecolor, linetype, linewidth)
    class(coop_asy) fp
    logical,optional::smooth
    real(sp),dimension(:,:),intent(IN)::x,y
    COOP_UNKNOWN_STRING,optional::colorfill
    COOP_UNKNOWN_STRING,optional:: linecolor, linetype
    real(sp),optional::linewidth
    COOP_INT  i,n, m,ipath
    COOP_STRING lineproperty
    n = coop_getdim("coop_asy_contour", size(x,1), size(y,1))
    m = coop_getdim("coop_asy_contour", size(x,2), size(y,2))
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
          lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
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
          call fp%write_coor( x(i, ipath), y(i, ipath) )
       enddo
    enddo
  end subroutine coop_asy_contour_mult_s


  subroutine coop_asy_contour_arr_d(fp, z, xmin, xmax, ymin, ymax, cvals, colorfill, smooth, linecolor, linetype, linewidth)
    class(coop_asy) fp
    real(dl),dimension(:,:)::z
    real(dl) xmin, xmax, ymin, ymax
    real(dl) cvals(:)
    COOP_INT  nx, ny, i, nc
    COOP_STRING,dimension(:)::colorfill
    COOP_STRING,dimension(:)::linecolor, linetype
    real(sp),dimension(:),optional::linewidth
    logical smooth
    nx = size(z, 1)
    ny = size(z, 2)
    write(fp%unit, "(A)") "CONTOUR"
    write(fp%unit, "(A)") "2" !!type 2 contour
    call fp%write_limits( real(xmin, sp), real(xmax, sp), real(ymin, sp), real(ymax, sp))
    if(present(linewidth))then
       nc = coop_getdim("coop_asy_contour_arr_d", size(cvals), size(colorfill), size(linecolor), size(linetype), size(linewidth))
    else
       nc = coop_getdim("coop_asy_contour_arr_d", size(cvals), size(colorfill), size(linecolor), size(linetype))
    endif
    write(fp%unit, "(I8)") nc
    write(fp%unit, "("//trim(coop_num2str(nc))//"G14.5)") cvals
    do i = 1, nc
       write(fp%unit, "(A)") trim(colorfill(i))
       if(present(linewidth))then
          write(fp%unit, "(A)") trim(coop_asy_linestyle(linecolor(i), linetype(i), linewidth(i)))
       else
          write(fp%unit, "(A)") trim(coop_asy_linestyle(linecolor(i), linetype(i)))
       endif
    enddo
    if(smooth)then
       write(fp%unit, "(A)") "1"
    else
       write(fp%unit, "(A)") "0"
    endif
    write(fp%unit, "(2I8)") nx, ny
    do i=1,nx
       write(fp%unit, "("//trim(coop_num2str(ny))//"G14.5)") real(z(i,:))
    enddo
  end subroutine coop_asy_contour_arr_d
  


  subroutine coop_asy_contour_arr_s(fp, z, xmin, xmax, ymin, ymax, cvals, colorfill, smooth, linecolor, linetype, linewidth)
    class(coop_asy) fp
    real(sp),dimension(:,:)::z
    real(sp) xmin, xmax, ymin, ymax
    real(sp) cvals(:)
    COOP_INT  nx, ny, i, nc
    COOP_STRING,dimension(:)::colorfill
    COOP_STRING,dimension(:)::linecolor, linetype
    real(sp),dimension(:),optional::linewidth
    logical smooth
    nx = size(z, 1)
    ny = size(z, 2)
    write(fp%unit, "(A)") "CONTOUR"
    write(fp%unit, "(A)") "2" !!type 2 contour
    call fp%write_limits(xmin, xmax, ymin, ymax)
    if(present(linewidth))then
       nc = coop_getdim("coop_asy_contour_arr_s", size(cvals), size(colorfill), size(linecolor), size(linetype), size(linewidth))
    else
       nc = coop_getdim("coop_asy_contour_arr_s", size(cvals), size(colorfill), size(linecolor), size(linetype))
    endif
    write(fp%unit, "(I8)") nc
    write(fp%unit, "("//trim(coop_num2str(nc))//"G14.5)") cvals
    do i = 1, nc
       write(fp%unit, "(A)") trim(colorfill(i))
       if(present(linewidth))then
          write(fp%unit, "(A)") trim(coop_asy_linestyle(linecolor(i), linetype(i), linewidth(i)))
       else
          write(fp%unit, "(A)") trim(coop_asy_linestyle(linecolor(i), linetype(i)))
       endif
    enddo
    if(smooth)then
       write(fp%unit, "(A)") "1"
    else
       write(fp%unit, "(A)") "0"
    endif
    write(fp%unit, "(2I8)") nx, ny
    do i=1,nx
       write(fp%unit, "("//trim(coop_num2str(ny))//"G14.5)") z(i,:)
    enddo
  end subroutine coop_asy_contour_arr_s




  subroutine coop_asy_density_d(fp, z, xmin, xmax, ymin, ymax, label, zmin, zmax)
    class(coop_asy) fp
    real(dl),dimension(:,:)::z
    real(dl) xmin, xmax, ymin, ymax
    real(dl),optional::zmin, zmax
    COOP_UNKNOWN_STRING, optional::label
    real(dl) minz, maxz
    COOP_INT  nx, ny, i
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
    call fp%write_limits(real(xmin, sp), real(xmax, sp), real(ymin, sp), real(ymax, sp))
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
       write(fp%unit, "("//trim(coop_num2str(ny))//"G14.5)") real(z(i,:))
    enddo
  end subroutine coop_asy_density_d

  subroutine coop_asy_density_s(fp, z, xmin, xmax, ymin, ymax, label, zmin, zmax)
    class(coop_asy) fp
    real(sp),dimension(:,:)::z
    real(sp) xmin, xmax, ymin, ymax
    real(sp),optional::zmin, zmax
    real(sp) minz, maxz
    COOP_INT  nx, ny, i
    COOP_UNKNOWN_STRING, optional::label
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
    call fp%write_limits(xmin, xmax, ymin, ymax)
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
       write(fp%unit, "("//trim(coop_num2str(ny))//"G14.5)") z(i,:)
    enddo
  end subroutine coop_asy_density_s


 subroutine coop_asy_irregular_density_d(fp, x, y, z, label, xmin, xmax, ymin, ymax, zmin, zmax)
    class(coop_asy) fp
    real(dl),dimension(:)::x, y, z
    COOP_UNKNOWN_STRING, optional::label
    real(dl),optional::xmin, xmax, ymin, ymax, zmin, zmax
    COOP_INT  i, n
    real(dl) minx, miny, minz, maxx, maxy, maxz
    n = coop_getdim("coop_asy_irregular_density", size(x), size(y), size(z))
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
    call fp%write_limits(real(minx, sp), real(maxx, sp), real(miny, sp), real(maxy, sp))
    write(fp%unit, "(2G14.5)") minz, maxz
    write(fp%unit,"(A)") "1"
    write(fp%unit,"(I10)") n
    do i=1,n
       write(fp%unit, "(3G14.5)") real(x(i)), real(y(i)), real(z(i))
    enddo
  end subroutine coop_asy_irregular_density_d

  subroutine coop_asy_irregular_density_s(fp, x, y, z, label, xmin, xmax, ymin, ymax, zmin, zmax)
    class(coop_asy) fp
    real(sp),dimension(:)::x, y, z
    COOP_UNKNOWN_STRING, optional::label
    real(sp),optional::xmin, xmax, ymin, ymax, zmin, zmax
    COOP_INT  i, n
    real(sp) minx, miny, minz, maxx, maxy, maxz
    n = coop_getdim("coop_asy_irregular_density", size(x), size(y), size(z))
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
    call fp%write_limits(minx, maxx, miny, maxy) !real(minx, sp), real(maxx, sp), real(miny, sp), real(maxy, sp))
    write(fp%unit, "(2G14.5)") minz, maxz
    write(fp%unit,"(A)") "1"
    write(fp%unit,"(I10)") n
    do i=1,n
       write(fp%unit, "(3G14.5)") x(i), y(i), z(i)
    enddo
  end subroutine coop_asy_irregular_density_s

  subroutine coop_asy_clip_d(fp, x, y, smooth)
    class(coop_asy) fp
    logical,optional::smooth
    real(dl),dimension(:),intent(IN)::x,y
    COOP_INT  i,n
    n = coop_getdim("coop_asy_clip", size(x), size(y))
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
       call fp%write_coor(real(x(i), sp), real(y(i), sp))
    enddo

  end subroutine coop_asy_clip_d

  subroutine coop_asy_clip_s(fp, x, y, smooth)
    class(coop_asy) fp
    logical,optional::smooth
    real(sp),dimension(:),intent(IN)::x,y
    COOP_INT  i,n
    n = coop_getdim("coop_asy_clip", size(x), size(y))
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
       call fp%write_coor( x(i), y(i))
    enddo
  end subroutine coop_asy_clip_s

  function coop_asy_linestyle(color, linetype, linewidth)
    COOP_UNKNOWN_STRING color
    COOP_UNKNOWN_STRING,optional::linetype
    real(sp),optional::linewidth
    COOP_STRING coop_asy_linestyle
    if(trim(color).eq."" .or. trim(color).eq."NULL")then
       coop_asy_linestyle = "NULL"
       return
    endif
    if(present(linetype))then
       coop_asy_linestyle = trim(color)//"_"//trim(linetype)
    else
       coop_asy_linestyle = trim(color)//"_solid"
    endif
    if(present(linewidth))then
       coop_asy_linestyle = trim(coop_asy_linestyle)//"_"//trim(coop_num2str(linewidth))
    endif
  end function coop_asy_linestyle

  subroutine coop_asy_path_append_s(path, x, y)
    real(sp) x, y
    type(coop_asy_path) path
    call path%l%push( (/ x, y /) )
  end subroutine coop_asy_path_append_s

  subroutine coop_asy_path_append_d(path, x, y)
    real(dl) x, y
    type(coop_asy_path) path
    call coop_asy_path_append_s(path, real(x), real(y))
  end subroutine coop_asy_path_append_d

  subroutine coop_asy_path_close(path)
    type(coop_asy_path) path
    if(allocated(path%length))then
       path%nclosed = path%nclosed + 1
       if(path%nclosed .gt. coop_asy_path_max_nclosed) stop "coop_asy_path_close: too many closed paths"
       path%length(path%nclosed) = path%l%n - sum(path%length(1:path%nclosed-1))
    else
       allocate(path%length(coop_asy_path_max_nclosed))
       path%nclosed = 1
       path%length(1) = path%l%n
    endif
  end subroutine coop_asy_path_close

  subroutine coop_asy_path_clear(path)
    type(coop_asy_path) path
    if(allocated(path%length))deallocate(path%length)
    call path%l%init()
  end subroutine coop_asy_path_clear

  !!generate a path that encloses the region of f(x, y) > threshold
  subroutine coop_asy_path_2dcl(path, f, xmin, xmax, ymin, ymax, threshold, resolution)
    COOP_INT ,parameter::default_resolution = 64
    COOP_INT ,optional:: resolution
    type(coop_asy_path) path
    external f
    real(sp) f
    real(sp) xmin, xmax, ymin, ymax, threshold, dx, dy, dxby2, dyby2
    COOP_INT  n
    logical,dimension(:,:),allocatable::above
    COOP_INT ,dimension(:,:,:),allocatable::lines
    COOP_INT  i, j, imin, jmin, last
    call coop_asy_path_clear(path) !!clear the path
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
          call coop_asy_path_close(path)
          goto 100 
       else
          call insert_point()
       endif
    enddo

  contains

    subroutine insert_point()
      real(sp) x, y
      real(sp) fp(2), s
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
      call coop_asy_path_append(path, x, y)
    end subroutine insert_point

  end subroutine coop_asy_path_2dcl


  subroutine coop_asy_path_from_array_center(path, f, xmin, xmax, ymin, ymax, threshold)
    type(coop_asy_path) path
    real(sp) f(:,:)
    real(sp) xmin, xmax, ymin, ymax, threshold, dx, dy
    COOP_INT  nx, ny, n
    nx = size(f, 1)
    ny = size(f, 2)
    n = min(max(nx, ny,16)*2, 256) !!higher resolution is useless since we do bilinear interpolation
    dx = (xmax-xmin)/nx
    dy = (ymax- ymin)/ny
    call coop_asy_path_2dcl(path, finterp, xmin, xmax, ymin, ymax, threshold, n)
  contains
    
    function  finterp(x, y)
      real(sp) x, y, finterp
      COOP_INT  i, j, ip1, jp1
      real(sp) ri, rj
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
    
  end subroutine coop_asy_path_from_array_center



  subroutine coop_asy_path_from_array(path, f, xmin, xmax, ymin, ymax, threshold)
    type(coop_asy_path) path
    real(sp) f(:,:)
    real(sp) xmin, xmax, ymin, ymax, threshold, dx, dy
    COOP_INT  nx, ny, n
    nx = size(f, 1)
    ny = size(f, 2)
    n = min(max(nx, ny,16)*2, 256) !!higher resolution is useless since we do bilinear interpolation
    dx = (xmax-xmin)/(nx-1)
    dy = (ymax- ymin)/(ny-1)
    call coop_asy_path_2dcl(path, finterp, xmin, xmax, ymin, ymax, threshold, n)
  contains
    
    function  finterp(x, y)
      real(sp) x, y, finterp
      COOP_INT  i, j
      real(sp) ri, rj
      ri = (x-xmin)/dx + 1.
      rj = (y-ymin)/dy + 1.
      i = min(max(floor(ri), 1), nx-1)
      j = min(max(floor(rj), 1), ny-1)
      ri = ri - i
      rj = rj - j
      finterp = (f(i, j)*(1.-ri) + f(i+1, j)*ri)*(1.-rj)+(f(i, j+1)*(1.-ri) + f(i+1, j+1)*ri)*rj
    end function finterp
    
  end subroutine coop_asy_path_from_array


  subroutine coop_asy_legend_default(fp)
    class(coop_asy) fp
    call coop_asy_legend_location(fp, "N")
  end subroutine coop_asy_legend_default

  subroutine coop_asy_legend_location(fp, loc, cols)
    class(coop_asy) fp
    COOP_UNKNOWN_STRING :: loc
    COOP_INT , optional::cols
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
  end subroutine coop_asy_legend_location

  subroutine coop_asy_legend_d(fp, x, y, cols)
    class(coop_asy) fp
    real(dl)x, y
    COOP_INT ,optional::cols
    write(fp%unit, "(A)") "LEGEND"
    write(fp%unit, "(A)") "NULL"
    write(fp%unit, "(2G15.4)") x, y
    if(present(cols))then
       write(fp%unit, "(I8)")  cols
    else
       write(fp%unit, "(I8)")  1
    endif
  end subroutine coop_asy_legend_d

  subroutine coop_asy_legend_s(fp, x, y, cols)
    class(coop_asy) fp
    real(sp) x, y
    COOP_INT ,optional::cols
    write(fp%unit, "(A)") "LEGEND"
    write(fp%unit, "(A)") "NULL"
    write(fp%unit, "(2G15.4)") x, y
    if(present(cols))then
       write(fp%unit, "(I8)")  cols
    else
       write(fp%unit, "(I8)")  1
    endif
  end subroutine coop_asy_legend_s


  subroutine coop_asy_topaxis_s(fp, xmin, xmax, islog,label)
    real(sp) xmin, xmax
    logical islog
    class(coop_asy) fp
    COOP_UNKNOWN_STRING label
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
  end subroutine coop_asy_topaxis_s

  subroutine coop_asy_topaxis_d(fp, xmin, xmax, islog,label)
    real(dl) xmin, xmax
    logical islog
    class(coop_asy) fp
    COOP_UNKNOWN_STRING label
    call coop_asy_topaxis_s(fp, real(xmin), real(xmax), islog,label)
  end subroutine coop_asy_topaxis_d


  subroutine coop_asy_rightaxis_s(fp, ymin, ymax, islog, label)
    real(sp) ymin, ymax
    logical islog
    class(coop_asy) fp
    COOP_UNKNOWN_STRING label
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
  end subroutine coop_asy_rightaxis_s

  subroutine coop_asy_rightaxis_d(fp, ymin, ymax, islog, label)
    real(dl) ymin, ymax
    logical islog
    COOP_UNKNOWN_STRING label
    class(coop_asy) fp
    call coop_asy_rightaxis_s(fp, real(ymin), real(ymax), islog,label)
  end subroutine coop_asy_rightaxis_d

  subroutine coop_asy_error_bar_d(fp, x, y, dy_minus, dy_plus, dx_minus, dx_plus, color)
    class(coop_asy) fp
    real(dl) x, y
    real(dl),optional::dy_minus, dy_plus, dx_minus, dx_plus
    COOP_UNKNOWN_STRING, optional::color
    if(present(color))then
       call coop_asy_label(fp, "$\bullet$", x, y, color)
       if(present(dy_minus))then
          call coop_asy_line(fp, x, y, x, y-dy_minus, linewidth = 1., color=color)
          call coop_asy_label(fp, "-", x, y-dy_minus, color=color)
       endif
       if(present(dy_plus))then
          call coop_asy_line(fp, x, y, x, y+dy_plus, linewidth = 1., color= color)
          call coop_asy_label(fp, "-", x, y+dy_plus, color=color)
       endif
       if(present(dx_minus))then
          call coop_asy_line(fp, x, y, x-dx_minus, y, linewidth = 1., color=color)
          call coop_asy_label(fp, "{\tiny $|$}", x-dx_minus, y, color=color)
       endif
       if(present(dx_plus))then
          call coop_asy_line(fp, x, y, x+dx_plus, y, linewidth = 1., color=color)
          call coop_asy_label(fp, "{\tiny $|$}", x+dx_plus, y, color=color)
       endif
    else
       call coop_asy_label(fp, "$\bullet$", x, y)
       if(present(dy_minus))then
          call coop_asy_line(fp, x, y, x, y-dy_minus, linewidth = 1.)
          call coop_asy_label(fp, "-", x, y-dy_minus)
       endif
       if(present(dy_plus))then
          call coop_asy_line(fp, x, y, x, y+dy_plus, linewidth = 1.)
          call coop_asy_label(fp, "-", x, y+dy_plus)
       endif
       if(present(dx_minus))then
          call coop_asy_line(fp, x, y, x-dx_minus, y, linewidth = 1.)
          call coop_asy_label(fp, "{\tiny $|$}", x-dx_minus, y)
       endif
       if(present(dx_plus))then
          call coop_asy_line(fp, x, y, x+dx_plus, y, linewidth = 1.)
          call coop_asy_label(fp, "{\tiny $|$}", x+dx_plus, y)
       endif
    endif
  end subroutine coop_asy_error_bar_d

  subroutine coop_asy_error_bar_s(fp, x, y, dy_minus, dy_plus, dx_minus, dx_plus, color)
    class(coop_asy) fp
    real(sp) x, y
    real(sp),optional::dy_minus, dy_plus, dx_minus, dx_plus
    COOP_UNKNOWN_STRING, optional::color
    if(present(color))then
       call coop_asy_label(fp, "$\bullet$", x, y, color)
       if(present(dy_minus))then
          call coop_asy_line(fp, x, y, x, y-dy_minus, linewidth = 1., color=color)
          call coop_asy_label(fp, "-", x, y-dy_minus, color=color)
       endif
       if(present(dy_plus))then
          call coop_asy_line(fp, x, y, x, y+dy_plus, linewidth = 1., color= color)
          call coop_asy_label(fp, "-", x, y+dy_plus, color=color)
       endif
       if(present(dx_minus))then
          call coop_asy_line(fp, x, y, x-dx_minus, y, linewidth = 1., color=color)
          call coop_asy_label(fp, "{\tiny $|$}", x-dx_minus, y, color=color)
       endif
       if(present(dx_plus))then
          call coop_asy_line(fp, x, y, x+dx_plus, y, linewidth = 1., color=color)
          call coop_asy_label(fp, "{\tiny $|$}", x+dx_plus, y, color=color)
       endif
    else
       call coop_asy_label(fp, "$\bullet$", x, y)
       if(present(dy_minus))then
          call coop_asy_line(fp, x, y, x, y-dy_minus, linewidth = 1.)
          call coop_asy_label(fp, "-", x, y-dy_minus)
       endif
       if(present(dy_plus))then
          call coop_asy_line(fp, x, y, x, y+dy_plus, linewidth = 1.)
          call coop_asy_label(fp, "-", x, y+dy_plus)
       endif
       if(present(dx_minus))then
          call coop_asy_line(fp, x, y, x-dx_minus, y, linewidth = 1.)
          call coop_asy_label(fp, "{\tiny $|$}", x-dx_minus, y)
       endif
       if(present(dx_plus))then
          call coop_asy_line(fp, x, y, x+dx_plus, y, linewidth = 1.)
          call coop_asy_label(fp, "{\tiny $|$}", x+dx_plus, y)
       endif
    endif
  end subroutine coop_asy_error_bar_s


  function coop_asy_rgb_color_s(r, g, b) result(color)
    real(sp) r, g, b
    COOP_SHORT_STRING color
    color = "RGB:"//trim(coop_num2str(nint(min(1., max(0., r))*255.*8.)/8.))//":"//trim(coop_num2str(nint(min(1., max(0., g))*255.*8.)/8.))//":"//trim(coop_num2str(nint(min(1., max(0., b))*255.*8.)/8.))
  end function coop_asy_rgb_color_s

  function coop_asy_rgb_color_d(r,g,b) result(color)
    real(dl) r, g, b
    COOP_SHORT_STRING color
    color = coop_asy_rgb_color_s(real(r), real(g), real(b))
  end function coop_asy_rgb_color_d


  function coop_asy_gray_color_s(gray) result(color)
    real(sp) gray
    COOP_SHORT_STRING color
    color = "GRAY:"//trim(coop_num2str(nint(min(1., max(0., gray))*255.*8.)/8.))
  end function coop_asy_gray_color_s


  function coop_asy_gray_color_d(gray) result(color)
    real(dl) gray
    COOP_SHORT_STRING color
    color = coop_asy_gray_color_s(real(gray))
  end function coop_asy_gray_color_d

end module coop_asy_mod