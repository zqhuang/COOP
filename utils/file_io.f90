module coop_file_mod
  use coop_wrapper_typedef
  use coop_list_mod
  use coop_matrix_mod
  implicit none
#include "constants.h"


private

public::coop_file, coop_copy_file, coop_delete_file, coop_create_file, coop_create_directory, coop_delete_directory, coop_file_numcolumns, coop_file_numlines, coop_load_dictionary, coop_free_file_unit, coop_file_exists, coop_file_encrypt, coop_file_decrypt, coop_string_encrypt, coop_string_decrypt, coop_file_load_function, coop_dir_exists, coop_export_dictionary, coop_import_matrix, coop_export_matrix, coop_file_load_realarr, coop_open_file, coop_dynamic_array_real

  character,parameter::coop_text_comment_symbol = "#"


  Type coop_file
     COOP_INT::unit = 0
     integer(8) pt  !!for binary mode
     COOP_SHORT_STRING mode
     COOP_LONG_STRING path
   contains
     procedure::open => coop_file_open
     procedure::close => coop_file_close
     procedure::skip_lines => coop_file_skiplines
     procedure::read_int => coop_file_readline_int
     procedure::read_int_array => coop_file_readline_int_arr
     procedure::read_real => coop_file_readline_real
     procedure::read_real_array => coop_file_readline_real_arr
     procedure::read_string => coop_file_readline_string
  End type coop_file

  COOP_INT,parameter::coop_file_unit_min = coop_tmp_file_unit + 10
  COOP_INT,parameter::coop_file_unit_max = coop_file_unit_min + 100

!!$  interface coop_file
!!$     procedure coop_open_file
!!$  end interface coop_file

  interface coop_import_matrix
     module procedure coop_import_matrix_s, coop_import_matrix_d
  end interface coop_import_matrix

  interface coop_export_matrix
     module procedure coop_export_matrix_s, coop_export_matrix_d
  end interface coop_export_matrix
  
    type coop_dynamic_array_real
       COOP_INT::nrows = 0
       COOP_INT::ncols = 0
       COOP_REAL, dimension(:, :),allocatable::f
     contains
       procedure::init => coop_dynamic_array_real_init
       procedure::load_txt => coop_dynamic_array_real_load_txt
       procedure::dump_txt => coop_dynamic_array_real_dump_txt
       procedure::free => coop_dynamic_array_real_free
    end type coop_dynamic_array_real

contains


  
  subroutine coop_dynamic_array_real_init(this, nrows, ncols)
    class(coop_dynamic_array_real)::this
    COOP_INT::nrows, ncols
    if(this%nrows .eq.  nrows  .and. this%ncols .eq. ncols)return
    if(allocated(this%f))deallocate(this%f)
    this%nrows = nrows
    this%ncols = ncols
    allocate(this%f(nrows, ncols))
  end subroutine coop_dynamic_array_real_init

  subroutine coop_dynamic_array_real_free(this)
    class(coop_dynamic_array_real)::this
    if(allocated(this%f))deallocate(this%f)
    this%nrows = 0
    this%ncols = 0
  end subroutine coop_dynamic_array_real_free

  subroutine coop_dynamic_array_real_load_txt(this, filename)
    class(coop_dynamic_array_real)::this
    COOP_UNKNOWN_STRING::filename
    type(coop_file)::fp
    COOP_INT::i
    call this%free()
    this%nrows = coop_file_NumLines(filename)
    this%ncols = coop_file_NumColumns(filename)
    if(this%nrows .gt. 0 .and. this%ncols .gt. 0)then
       allocate(this%f(this%nrows, this%ncols))
       call fp%open(filename, "r")
       do i=1, this%nrows
          if(.not. fp%read_real_array(this%f(i, :)))then
             write(*,*) trim(filename)
             stop "data file is broken"
          endif
       enddo
       call fp%close()
    endif
  end subroutine coop_dynamic_array_real_load_txt


    subroutine coop_dynamic_array_real_dump_txt(this, filename)
    class(coop_dynamic_array_real)::this
    COOP_UNKNOWN_STRING::filename
    type(coop_file)::fp
    COOP_STRING::fmt
    COOP_INT::i
    fmt = "("//COOP_STR_OF(this%ncols)//"E16.7)"
    call fp%open(filename, "w")
    do i=1, this%nrows
       write(fp%unit, trim(fmt)) this%f(i, :)
    enddo
    call fp%close()    
  end subroutine coop_dynamic_array_real_dump_txt

  

  subroutine coop_file_load_function(fname, col1, col2, func, check_boundary)
    COOP_UNKNOWN_STRING fname
    type(coop_file)fp
    logical check_boundary, xlog, ylog
    COOP_INT col1, col2, nrows, ncols, i
    type(coop_function)::func
    COOP_REAL,dimension(:),allocatable::x, f, line
    COOP_REAL::meanx, meanf
    if(.not.coop_file_exists(fname))call coop_return_error("file_load_function","file "//trim(fname)//" does not exists", "stop")
    ncols = coop_file_NumColumns(fname)
    if(ncols .lt. col1 .or. ncols .lt. col2) call coop_return_error("file_load_function","Missing column", "stop")
    nrows = coop_file_NumLines(fname)
    if(nrows.le.1) call coop_return_error("file_load_function","Missing column", "stop")
    allocate(x(nrows), f(nrows), line(ncols))
    call fp%open(fname, "r")
    do i=1, nrows
       read(fp%unit, *, ERR=100) line
       x(i) = line(col1)
       f(i) = line(col2)
    enddo
    call fp%close()
    meanx= sum(x)/nrows
    meanf = sum(f)/nrows
    xlog = all(x.gt.0.d0) .and. (count(x.lt.meanx) .gt. (nrows*4/5))
    ylog = all(f.gt.0.d0) .and. (count(f.lt.meanf).gt. (nrows*4/5))
    if(xlog)then
       if(is_uniform(nrows, log(x)))then
          call func%init(n=nrows, xmin=x(1), xmax=x(nrows),f=f, check_boundary = check_boundary, xlog=xlog, ylog=ylog)
          return
       endif
    else
       if(is_uniform(nrows, x))then
          call func%init(n=nrows, xmin=x(1), xmax=x(nrows),f=f, check_boundary = check_boundary, xlog=xlog, ylog=ylog)
          return
       endif
    endif
    call func%init_NonUniform(x = x, f=f,  xlog=xlog, ylog=ylog, check_boundary = check_boundary)
    return
100 call coop_return_error("file_load_function","Missing row", "stop")

  contains
    function is_uniform(n, a)
      COOP_INT n, i
      logical is_uniform
      COOP_REAL a(n), da
      da = a(2)-a(1)
      do i=3, n
         if(abs((a(i)-a(i-1))/da-1.d0).gt. 1.d-4)then
            is_uniform=.false.
            return
         endif
      enddo
      is_uniform  = .true.
      return
    end function is_uniform
  end subroutine coop_file_load_function
  
  function coop_free_file_unit() result(lun)
    COOP_INT lun
    logical used
    do lun = coop_file_unit_min, coop_file_unit_max
       Inquire(Unit=lun, Opened=used)
       If (.not. used) return
    enddo
    lun = 0
  End function coop_free_file_unit

  subroutine coop_delete_file(filename)
    COOP_UNKNOWN_STRING filename
    if(coop_file_exists(trim(filename))) call system("rm -f "//trim(filename))
  end subroutine coop_delete_file

  subroutine coop_create_file(filename)
    COOP_UNKNOWN_STRING filename
    if(.not.coop_file_exists(trim(filename))) call system("touch "//trim(filename))
  end subroutine coop_create_file

  subroutine coop_create_directory(directory)
    COOP_UNKNOWN_STRING directory
    if(.not. coop_dir_exists(directory)) call system("mkdir "//trim(directory))
  end subroutine coop_create_directory

  subroutine coop_delete_directory(directory)
    COOP_UNKNOWN_STRING directory
    call system("rm -rf "//trim(directory))
  end subroutine coop_delete_directory


  subroutine coop_copy_file(file1, file2)
    COOP_UNKNOWN_STRING file1, file2
    call system("cp "//trim(file1)//" "//trim(file2))
  end subroutine coop_copy_file


  subroutine coop_file_open(fp, filename, mode, recl)
    class(coop_file)::fp
    COOP_UNKNOWN_STRING filename
    COOP_UNKNOWN_STRING,optional::mode
    COOP_INT,optional:: recl
    !!supported modes:
    !!    'w' (text mode for read/write),
    !!    'b' (binary mode for read/write), 
    !!    'r' (text mode, read only),
    !!    'rb' (binary, read only)
    !!    'a' (text mode for append)
    !!    'u' (unformated)
    !!    'ua' (unformated, append)    
    !!    'ur' (unformated, read only)
#include "file_init.h"    
  end subroutine coop_file_open

  function coop_open_file(filename,  mode, recl) result(fp)
    type(coop_file) fp
    COOP_UNKNOWN_STRING filename
    COOP_UNKNOWN_STRING,optional::mode
    COOP_INT,optional:: recl
#include "file_init.h"    
  end function coop_open_file

  function coop_open_file_skip_comments(fname) result(fp)
    COOP_UNKNOWN_STRING fname
    type(coop_file) fp
    COOP_INT i, nc
    COOP_LONG_STRING inline
    call fp%open(fname, "r")
    nc = 0
50  read(fp%unit, "(A)", ERR=100, END=100) inline
    inline = adjustl(inline)
    if(trim(inline) .eq. "" .or. inline(1:1).eq.coop_text_comment_symbol)then
       nc = nc + 1
       goto 50
    endif
100 call fp%close()
    call fp%open(fname, "r")
    do i=1, nc
       read(fp%unit, "(A)") inline
    enddo
  end function coop_open_file_skip_comments
  
  subroutine coop_file_close(fp)
    class(coop_file) fp
    close(fp%unit)
    fp%pt = 1
    fp%unit = 0
    fp%path = ""
  end subroutine coop_file_close


  function coop_file_exists(FileName) result(file_exists)
    COOP_UNKNOWN_STRING, INTENT(IN) :: FileName
    logical File_Exists
    if(trim(adjustl(filename)).eq."")then
       file_exists = .false.
       return
    endif
    inquire(FILE=trim(adjustl(fileName)), EXIST= File_Exists)
  end function coop_file_exists

  function coop_dir_exists(dirName) result(dir_exists)
    COOP_UNKNOWN_STRING, INTENT(IN) :: dirName
    logical dir_exists
    if(trim(adjustl(dirname)).eq."")then
       dir_exists = .false.
       return
    endif
#ifdef __GFORTRAN__
    inquire(FILE=trim(adjustl(dirName)), EXIST = dir_exists)    
#else
    inquire(DIRECTORY=trim(adjustl(dirName)), EXIST = dir_exists)
#endif
  end function coop_dir_exists
  

  subroutine coop_File_OverWrite_Query(FileName)
    COOP_UNKNOWN_STRING, INTENT(IN) :: FileName
    character(LEN=3) ans
    if(coop_File_Exists(trim(Filename)))then
       write(*,*) "The file "//trim(FileName)//" already exists"
       write(*,*) "Do you want to overwrite it (yes/no)?"
       read(*,*) ans
       if(ans.ne."yes" .and. ans.ne."YES") stop       
    endif
  end subroutine coop_File_OverWrite_Query

  Function coop_File_NumColumns(FileName) result(NumColumns)
    COOP_UNKNOWN_STRING FileName
    COOP_LONG_STRING FirstLine
    COOP_INT NumColumns,I,Length,J
    type(coop_file) fp
    NumColumns=0
    if(.not. coop_File_Exists(filename)) return
    call fp%open(trim(filename), "r")
    if( coop_file_readline_string(fp, firstline) )then
       numcolumns = coop_string_contain_numbers(firstline)
    endif
100 call fp%close()
  End Function Coop_File_NumColumns

  subroutine coop_file_load_realarr(filename, rl, num)
    COOP_UNKNOWN_STRING, INTENT(IN)::filename
    Type(Coop_file) fp
    type(coop_list_realarr)::rl
    COOP_INT, optional::num
    COOP_INT::ncols
    COOP_LONG_STRING inline    
    COOP_SINGLE, dimension(:),allocatable::s
    if(.not. coop_file_exists(filename))return
    if(present(num))then !!unformatted
       ncols = num
    else
       ncols = coop_file_numColumns(filename)
    endif
    if(ncols .le. 0)return
    allocate(s(ncols))
    if(present(num))then  !!unformatted file, with num entries per line
       call fp%open(filename, "ru")
       do
          read(fp%unit, ERR=100, END=100) s
          call rl%push(s)          
       enddo
    else
       call fp%open(filename, "r")
       do while(coop_file_readline_string(fp, inline))
          read(inline, *, ERR=100) s
          call rl%push(s)
       enddo
    end if
100 call fp%close()
    deallocate(s)
  end subroutine coop_file_load_realarr

  Function coop_File_NumLines(filename, num) result(Lines) !! # comment, empty lines are skipped
    COOP_UNKNOWN_STRING, INTENT(IN)::filename
    COOP_LONG_STRING inline
    COOP_INT Lines
    COOP_INT, optional::num
    COOP_REAL, dimension(:), allocatable::s
    Type(Coop_file) fp
    Lines = 0
    if(.not. coop_File_Exists(filename)) return
    if(present(num))then  !!unformatted file, with num entries per line
       allocate(s(num))
       call fp%open(filename, "ru")
       do
          read(fp%unit, ERR=100, END=100) s
          lines = lines + 1
       enddo
    else
       call fp%open(filename, "r")
       do while(coop_file_readline_string(fp, inline))
          lines = lines + 1
       enddo
    end if
100 call fp%close()
    if(present(num))deallocate(s)
  End Function Coop_File_NumLines


  Function Coop_file_readline_Int(fp, param) result(success)
    COOP_INT param
    logical success
    class(coop_file) fp
    COOP_LONG_STRING line
    Line = ""
    success = .false.
    do while(Line(1:1) .eq. coop_text_comment_symbol .or. Trim(Line) .eq. "")
       read(fp%unit, '(A)', End=200, Err=200) line
       line = adjustl(trim(line))
    enddo
    read(Line, *, Err = 200, End=200) param
    success = .true.
200 return
  End Function Coop_file_readline_Int

  Function Coop_file_readline_Real(fp, param) result(success)
    COOP_REAL param
    logical success
    class(coop_file) fp
    COOP_LONG_STRING line
    Line = ""
    success = .false.
    do while(Line(1:1) .eq. coop_text_comment_symbol .or. Trim(Line) .eq. "")
       read(fp%unit, '(A)', End=200, Err=200) line
       line = adjustl(trim(line))
    enddo
    read(Line, *, Err = 200, End=200) param
    success = .true.
200 return
  End Function Coop_file_readline_Real


  Function Coop_file_readline_Int_arr(fp, params) result(success)
    COOP_INT,dimension(:)::params
    logical success
    class(coop_file) fp
    COOP_LONG_STRING line
    Line = ""
    success = .false.
    do while(Line(1:1) .eq. coop_text_comment_symbol .or. Trim(Line) .eq. "")
       read(fp%unit, '(A)', End=200, Err=200) line
       line = adjustl(trim(line))
    enddo
    read(Line, *, Err = 200, End=200) params
    success = .true.
200 return
  End Function Coop_file_readline_Int_arr


  Function Coop_file_readline_Real_arr(fp, params) result(success)
    COOP_REAL,dimension(:)::params
    logical success
    class(coop_file) fp
    COOP_LONG_STRING line
    Line = ""
    success = .false.
    do while(Line(1:1) .eq. coop_text_comment_symbol .or. Trim(Line) .eq. "")
       read(fp%unit, '(A)', End=200, Err=200) line
       line = adjustl(trim(line))
    enddo
    read(Line, *, Err = 200, End=200) params
    success = .true.
200 return
  End Function Coop_file_readline_Real_arr

  Function Coop_file_readline_String(fp, line) result(success)
    logical success
    class(coop_file) fp
    COOP_UNKNOWN_STRING line
    Line = ""
    success = .false.
    do while(Line(1:1) .eq. coop_text_comment_symbol .or. Trim(Line) .eq. "")
       read(fp%unit, '(A)', End=200, Err=200) line
       line = adjustl(trim(line))
    enddo
    success = .true.
    return
200 line = ""
    return
  End Function Coop_file_readline_String

  Subroutine coop_File_SkipLines(fp, numlines)
    class(coop_file) fp
    COOP_INT numlines, lines
    COOP_LONG_STRING inline
    if(numlines.le.0)return
    Lines = 0
    do while(coop_file_readline_string(fp, inline))
       lines = lines + 1
       if(lines .ge. numlines) exit
    enddo
200 return
  end Subroutine Coop_File_SkipLines


  recursive subroutine coop_load_dictionary(fname, dict, delimitor, col_key, col_value)
    COOP_UNKNOWN_STRING fname
    COOP_STRING::path
    type(coop_dictionary) dict
    COOP_STRING line
    integer i
    type(coop_file) fp
    type(coop_list_string) l
    COOP_INT eqloc
    COOP_UNKNOWN_STRING, optional::delimitor
    COOP_INT, optional:: col_key, col_value
    COOP_INT colm, linelen
    COOP_SHORT_STRING:: delim
    COOP_STRING::val, subval
    if(present(col_key) .and. present(col_value))then  !!key and value from column col_key and column col_value
       if(present(delimitor))then
          delim = delimitor
       else
          delim = " ,;"//coop_tab
       endif
       colm = max(col_key, col_value)
       call fp%open(trim(fname), "r")
       do while(fp%read_string(line))
          call coop_string_to_list(line, l, delim)
          if(l%n .ge. colm)then
             call dict%insert(l%element(col_key), l%element(col_value))
          endif
       enddo
       call fp%close()
       call l%init()
    elseif(present(col_value))then !!key = # of lines
       if(present(delimitor))then
          delim = delimitor
       else
          delim = " ,;"//coop_tab
       endif
       call fp%open(trim(fname), "r")
       i = 0
       do while(fp%read_string(line))
          call coop_string_to_list(line, l, delim)
          if(l%n .ge. col_value)then
             i = i + 1
             call dict%insert(COOP_STR_OF(i), l%element(col_value))
          endif
       enddo
       call fp%close()
       call l%init()
    elseif(present(col_key))then  !!any thing after col_key is treated as value
       if(present(delimitor))then
          delim = delimitor
       else
          delim = " ,;"//coop_tab
       endif
       call fp%open(trim(fname), "r")
       do while(fp%read_string(line))
          call coop_string_to_list(line, l, delim)
          if(l%n .gt. col_key)then
             val = ""
             do i = col_key+1, l%n
                call l%get_element(i, subval)
                if(subval(1:1).eq. coop_text_comment_symbol)exit
                val = trim(val)//" "//trim(subval)
             enddo
             if(trim(val).ne."")call dict%insert(l%element(col_key), val)
          endif
       enddo
       call fp%close()
       call l%init()       
    else  !!use "=" to specify (key, value) pairs; free stype; allow DEFAULT(subfile) entry
       path = trim(coop_file_path_of(fname))
       if(present(delimitor))then
          delim = delimitor
       else
          delim = "="
       endif
       call fp%open(trim(fname), "r")
       do while(fp%read_string(line))
          linelen = len_trim(line)
          if(line(1:8) .eq. "DEFAULT(" .and. line(linelen:linelen) .eq. ")")then
             call coop_load_dictionary(trim(path)//adjustl(trim(line(9:linelen-1))), dict, trim(delim))
          end if
          eqloc = scan(line, trim(delim))
          if(eqloc .gt. 1 .and. eqloc .lt. len_trim(line))then
             if(trim(line(eqloc+1:)).ne."")call dict%insert(line(1:eqloc-1), line(eqloc+1:))
          endif
       enddo
       call fp%close()
    endif
  end subroutine coop_load_dictionary

  subroutine coop_export_dictionary(fname, dict)
    COOP_UNKNOWN_STRING ::fname
    type(coop_dictionary)::dict
    COOP_INT i
    type(coop_file)::fp
    call fp%open(fname,"w")
    do i = 1, dict%n
       write(fp%unit, "(A)") trim(dict%key(i))//" = "//trim(dict%val(i))
    enddo
    call fp%close()
  end subroutine coop_export_dictionary

  subroutine coop_file_encrypt(input, output)
    COOP_UNKNOWN_STRING input, output
    COOP_LONG_STRING line
    COOP_INT il
    type(coop_file)::fin, fout
    call fin%open(input)
    call fout%open(output)
    il = 0
    do
       il = il + 1
       read(fin%unit, "(A)", END=100, ERR=100) line
       line = trim(adjustl(line))
       if(trim(line).ne."") &
            call coop_string_encrypt(line, il**2+2017)
       write(fout%unit, "(A)") trim(line)
    enddo
100 call fin%close()
    call fout%close()
  end subroutine coop_file_encrypt


  subroutine coop_file_decrypt(input, output)
    COOP_UNKNOWN_STRING input, output
    COOP_LONG_STRING line
    COOP_INT il
    type(coop_file)::fin, fout
    call fin%open(input)
    call fout%open(output)
    il = 0
    do
       il = il + 1
       read(fin%unit, "(A)", END=100, ERR=100) line
       if(trim(line).ne."") call coop_string_decrypt(line, il**2+2017)
       write(fout%unit, "(A)") trim(line)
    enddo
100 call fin%close()
    call fout%close()
  end subroutine coop_file_decrypt



  subroutine coop_string_encrypt(str, seed)
    COOP_UNKNOWN_STRING str
    COOP_INT seed, p1, p2, ic, i
    type(coop_list_integer)::pl
    call coop_get_prime_numbers( mod(seed, 999983) + 128, pl)
    p1 = pl%element(pl%n)
    p2 = pl%element(pl%n-1)
    call pl%init()
    do i=1, len(str)
       ic = ichar(str(i:i))
       if(ic .ge. 39 .and. ic .le. 126)then
          ic = mod((ic-38)*p1, 89) + 38
          str(i:i) = char(ic)
       elseif(ic .ge. 33 .and. ic .le. 38)then
          ic = mod((ic-32)*p2, 7) + 32
          str(i:i) = char(ic)
       endif
    enddo
  end subroutine coop_string_encrypt

  subroutine coop_string_decrypt(str, seed)
    COOP_UNKNOWN_STRING str
    COOP_INT seed, p1, p2, q1, q2, ic, i
    type(coop_list_integer)::pl
    call coop_get_prime_numbers( mod(seed, 999983) + 128, pl)
    p1 = pl%element(pl%n)
    p2 = pl%element(pl%n-1)
    q1 = 1
    q2 = 1
    do while(q1.le.88)
       if(mod(q1*p1, 89).eq.1)exit
       q1 = q1 + 1
    enddo
    do while(q2.le.6)
       if(mod(q2*p2, 7) .eq. 1) exit
       q2 = q2 + 1
    enddo
    call pl%init()
    do i=1, len(str)
       ic = ichar(str(i:i))
       if(ic .ge. 39 .and. ic .le. 126)then
          ic = mod((ic-38)*q1, 89) + 38
          str(i:i) = char(ic)
       elseif(ic .ge. 33 .and. ic .le. 38)then
          ic = mod((ic-32)*q2, 7) + 32
          str(i:i) = char(ic)
       endif
    enddo
  end subroutine coop_string_decrypt


  subroutine coop_import_matrix_d(filename, mat, nx, ny, success)
    COOP_UNKNOWN_STRING::filename
    COOP_INT nx, ny
    COOP_REAL::mat(nx, ny)
    logical,optional::success
    type(coop_file)::fp
    call fp%open(filename, "r")
    if(present(success))then
       call coop_read_matrix(fp%unit,  mat, nx, ny, success)
    else
       call coop_read_matrix(fp%unit, mat, nx, ny)       
    endif
    call fp%close()
  end subroutine coop_import_matrix_d

  subroutine coop_import_matrix_s(filename, mat, nx, ny, success)
    COOP_UNKNOWN_STRING::filename
    COOP_INT nx, ny
    COOP_SINGLE::mat(nx, ny)
    logical,optional::success
    type(coop_file)::fp
    call fp%open(filename, "r")
    if(present(success))then
       call coop_read_matrix(fp%unit, mat, nx, ny, success)
    else
       call coop_read_matrix(fp%unit, mat, nx, ny)       
    endif
    call fp%close()
  end subroutine coop_import_matrix_s




  subroutine coop_export_matrix_d(filename, mat, nx, ny)
    COOP_UNKNOWN_STRING::filename
    COOP_INT nx, ny
    COOP_REAL::mat(nx, ny)
    type(coop_file)::fp
    call fp%open(filename, "w")
    call coop_write_matrix(fp%unit, mat, nx, ny)       
    call fp%close()
  end subroutine coop_export_matrix_d

  subroutine coop_export_matrix_s(filename, mat, nx, ny)
    COOP_UNKNOWN_STRING::filename
    COOP_INT nx, ny
    COOP_SINGLE::mat(nx, ny)
    type(coop_file)::fp
    call fp%open(filename, "w")
    call coop_write_matrix(fp%unit, mat, nx, ny)       
    call fp%close()
  end subroutine coop_export_matrix_s  


end module coop_file_mod
