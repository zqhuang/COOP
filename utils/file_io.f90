module file_io_utils
  use basic_utils
  use string_utils
  use matrix_utils
  use list_utils
  implicit none
#include "utils.h"

  character,parameter::text_comment_symbol = "#"

  Type file_pointer
     integer unit
     integer(8) pt  !!for binary mode
     SHORT_STRING mode
     LONG_STRING path
  End type file_pointer

  Interface Write_File
     Module Procedure Int_to_File, float_To_File,float_arr_to_file, int_arr_to_file
  End Interface

  Interface Read_File
     Module Procedure Int_From_File, float_From_File, Int_arr_From_File, Float_Arr_From_File
  ENd Interface

  InterFace File_ReadOneLine
     Module Procedure File_ReadOneLine_Int, File_ReadOneLine_Double, File_ReadOneLine_Int_arr, File_ReadOneLine_Double_arr, File_ReadOneLine_String
  End InterFace File_ReadOneLine

  integer,parameter::file_min_unit = 11
  integer,parameter::file_max_unit = 99

  logical::unit_opened(file_min_unit:file_max_unit) = .false.  !!1-6 may be used by the system; 7,8,9,10 are reserved for temporary operations
  private::unit_opened


contains

  function new_file_unit()
    integer new_file_unit
    new_file_unit = file_min_unit
    do while(unit_opened(new_file_unit))
       new_file_unit = new_file_unit + 1
       if(new_file_unit .gt. file_max_unit)stop "All file units have been used."
    enddo
    unit_opened(new_file_unit)=.true.
    return
  end function new_file_unit

  subroutine close_file_unit(funit)
    integer funit
    close(funit)
    unit_opened(funit) = .false.
  end subroutine close_file_unit

  subroutine delete_file(filename)
    UNKNOWN_STRING filename
    if(file_exists(trim(filename))) call system("rm -f "//trim(filename))
  end subroutine delete_file

  subroutine create_file(filename)
    UNKNOWN_STRING filename
    if(.not.file_exists(trim(filename))) call system("touch "//trim(filename))
  end subroutine create_file

  subroutine create_directory(directory)
    UNKNOWN_STRING directory
    call system("mkdir "//trim(directory))
  end subroutine create_directory

  subroutine copy_file(file1, file2)
    UNKNOWN_STRING file1, file2
    call system("cp "//trim(file1)//" "//trim(file2))
  end subroutine copy_file


  function open_file(filename,mode, recl) 
    !!supported modes:
    !!    'w' (text mode for read/write),
    !!    'b' (binary mode for read/write), 
    !!    'r' (text mode, read only),
    !!    'rb' (binary, read only)
    !!    'a' (text mode for append)
    !!    'u' (unformated)
    !!    'ur' (unformated, read only)
    type(file_pointer) open_file
    UNKNOWN_STRING filename
    UNKNOWN_STRING,optional::mode
    character(LEN=2)::open_mode
    integer funit
    integer,optional:: recl
    funit = file_min_unit
    !$omp critical
    do while (unit_opened(funit))
       funit = funit+1
       if (funit.gt.file_max_unit)then
          stop "error in open_txt_file: all file units are occupied"
       endif
    enddo
    !$omp end critical
    open_file%unit = funit
    open_file%pt = 1
    if(present(mode))then
       open_mode = adjustl(trim(mode))
    else
       open_mode = 'w'
    endif
    select case(trim(open_mode))
    case("w","W")
       open(funit,FILE=trim(filename),FORM="FORMATTED",STATUS="UNKNOWN",ACCESS='SEQUENTIAL',ACTION="READWRITE", Err=200)
       open_file%mode = 'txt'
    case("wq","WQ", "qw", "QW")
       call file_overwrite_query(trim(filename))
       open(funit,FILE=trim(filename),FORM="FORMATTED",STATUS="UNKNOWN",ACCESS='SEQUENTIAL',ACTION="READWRITE", Err=200)
       open_file%mode = 'txt'
    case("r", "R")
       open(funit,FILE=trim(filename),FORM="FORMATTED",STATUS="UNKNOWN", ACCESS='SEQUENTIAL', ACTION="READ", Err=200)
       open_file%mode = 'txt'
    case("a", "A")
       open(funit,FILE=trim(filename),FORM="FORMATTED",STATUS="UNKNOWN", ACCESS='APPEND',ACTION="READWRITE", Err=200)
       open_file%mode = 'txt'
    case("b", "B")
       if(present(recl))then
          open(UNIT=funit,FILE=trim(filename),FORM="UNFORMATTED",STATUS="UNKNOWN", ACCESS="DIRECT",ACTION="READWRITE",RECL=recl,ERR=200) 
          open_file%mode = 'bin'
       else
          open(UNIT=funit,FILE=trim(filename),FORM="UNFORMATTED",STATUS="UNKNOWN", ACCESS="DIRECT",ACTION="READWRITE",RECL=IB,ERR=200) 
          open_file%mode = 'bin'
       endif
    case("rb", "RB", "br", "BR")
       if(present(recl))then
          open(UNIT=funit,FILE=trim(filename),FORM="UNFORMATTED",STATUS="UNKNOWN", ACCESS="DIRECT",ACTION="READ",RECL=recl, ERR=200) 
          open_file%mode = 'bin'
       else
          open(UNIT=funit,FILE=trim(filename),FORM="UNFORMATTED",STATUS="UNKNOWN", ACCESS="DIRECT",ACTION="READ",RECL=IB, ERR=200) 
          open_file%mode = 'bin'
       endif
    case('u', "U")
       open(UNIT=funit,FILE=trim(filename),FORM="UNFORMATTED",STATUS="UNKNOWN", ACCESS="SEQUENTIAL",ACTION="READWRITE",ERR=200)        
       open_file%mode = "unf"
    case('uq', "UQ", "qu", "QU")
       call file_overwrite_query(trim(filename))
       open(UNIT=funit,FILE=trim(filename),FORM="UNFORMATTED",STATUS="UNKNOWN", ACCESS="SEQUENTIAL",ACTION="READWRITE",ERR=200)        
       open_file%mode = "unf"
    case('ru', "RU", "ur", "UR")
       open(UNIT=funit,FILE=trim(filename),FORM="UNFORMATTED",STATUS="UNKNOWN", ACCESS="SEQUENTIAL",ACTION="READ",ERR=200)        
       open_file%mode = "unf"
    case('AU', "au", "ua", "UA")
       open(UNIT=funit,FILE=trim(filename),FORM="UNFORMATTED",STATUS="UNKNOWN", ACCESS="APPEND",ACTION="READWRITE",ERR=200)        
       open_file%mode = "unf"
    case default
       write(*,*) "Error in open_file: unknow mode "//trim(mode)
       goto 200
    end select
    open_file%path = trim(filename)
    unit_opened(funit) = .true.
    return
200 write(*,*) "Error: can not open file "//trim(filename)
    stop
  end function open_file

  function open_file_skip_comments(fname) result(fp)
    UNKNOWN_STRING fname
    type(file_pointer) fp
    integer i, nc
    LONG_STRING inline
    fp = open_file(fname, "r")
    nc = 0
50  read(fp%unit, "(A)", ERR=100, END=100) inline
    inline = adjustl(inline)
    if(trim(inline) .eq. "" .or. inline(1:1).eq.text_comment_symbol)then
       nc = nc + 1
       goto 50
    endif
100 call close_file(fp)
    fp = open_file(fname, "r")
    do i=1, nc
       read(fp%unit, "(A)") inline
    enddo
  end function open_file_skip_comments
  
  subroutine dump_xy_to(fname, x, y)
    real(dl),dimension(:)::x,y
    integer i
    UNKNOWN_STRING fname
    type(file_pointer) fp
    fp = open_file(fname,"w")
    do i=1, min(size(x),size(y))
       write(fp%unit,'(2E16.7)') x(i),y(i)
    enddo
    call close_file(fp)
  end subroutine dump_xy_to

  subroutine dump_arrs(funit, x1, x2, x3)
    integer,optional::funit
    real(dl),dimension(:)::x1
    real(dl),dimension(:),optional::x2, x3
    integer i, n
    if(present(funit))then
       if(present(x2))then
          if(present(x3))then
             n= GetDim("dump_arrs", size(x1), size(x2), size(x3))
             do i=1,n
                write(funit,'(I6, 3G16.7)') i, x1(i), x2(i), x3(i)
             enddo
          else
             n= GetDim("dump_arrs", size(x1), size(x2))
             do i=1,n
                write(funit,'(I6, 2G16.7)') i, x1(i), x2(i)
             enddo
          endif
       else
          n= size(x1)
          do i=1,n
             write(funit,'(I6, G16.7)') i, x1(i)
          enddo
       endif
    else
       if(present(x2))then
          if(present(x3))then
             n= GetDim("dump_arrs", size(x1), size(x2), size(x3))
             do i=1,n
                write(*,'(I6, 3G16.7)') i, x1(i), x2(i), x3(i)
             enddo
          else
             n= GetDim("dump_arrs", size(x1), size(x2))
             do i=1,n
                write(*,'(I6, 2G16.7)') i, x1(i), x2(i)
             enddo
          endif
       else
          n= size(x1)
          do i=1,n
             write(*,'(I6, G16.7)') i, x1(i)
          enddo
       endif
    endif
  end subroutine dump_arrs

  subroutine dump_to_file(filename, arr, trans)
    type(file_pointer)::fp
    UNKNOWN_STRING filename
    real(dl),dimension(:,:),intent(in)::arr
    integer i, m, n
    logical,optional::trans
    logical do_trans
    m = size(arr, 1)
    n = size(arr, 2)
    fp = open_file(filename, "w")
    if(present(trans))then
       do_trans = trans
    else
       do_trans = .false.
    endif
    if(do_trans)then
       do i=1, m
          write(fp%unit, '(I8,'//trim(num2str(n))//'E16.7)') i, arr(i, :)
       enddo
    else
       do i=1, n
          write(fp%unit, '(I8,'//trim(num2str(m))//'E16.7)') i, arr(:, i)
       enddo
    endif
    call close_file(fp)
  end subroutine dump_to_file

  subroutine int_to_txt_file(fp,i,fmt)
    type(file_pointer) fp
    integer i
    UNKNOWN_STRING,optional::fmt
    if(present(fmt))then
       write(fp%unit,trim(fmt)) i
    else
       write(fp%unit,*) i
    endif
  end subroutine int_to_txt_file

  subroutine float_to_txt_file(fp,x,fmt)
    type(file_pointer) fp
    real(dl) x
    UNKNOWN_STRING,optional::fmt
    if(present(fmt))then
       write(fp%unit,trim(fmt)) x
    else
       write(fp%unit,*) x
    endif
  end subroutine float_to_txt_file

  subroutine int_arr_to_txt_file(fp,i,fmt)
    type(file_pointer)fp
    integer,dimension(:),intent(in)::i
    UNKNOWN_STRING,optional::fmt
    if(present(fmt))then
       write(fp%unit,trim(fmt)) i
    else
       write(fp%unit,"("//trim(Num2Str(size(i)))//"I12)") i
    endif
  end subroutine int_arr_to_txt_file

  subroutine float_arr_to_txt_file(fp,x,fmt)
    type(file_pointer) fp
    real(dl),dimension(:),intent(in)::x
    UNKNOWN_STRING,optional::fmt
    if(present(fmt))then
       write(fp%unit,trim(fmt)) x
    else
       write(fp%unit,"("//trim(Num2str(size(x)))//"G20.9)") x
    endif
  end subroutine float_arr_to_txt_file

  subroutine int_to_binary_file(fp, i)
    type(file_pointer) fp
    integer i
    write(fp%unit,REC=fp%pt) i
    fp%pt=fp%pt+1
  end subroutine int_to_binary_file

  subroutine int_arr_to_binary_file(fp,x)
    type(file_pointer)fp
    integer,dimension(:)::x
    integer i
    do i=1,size(x)
       call int_to_binary_file(fp,x(i))
    enddo
  end subroutine int_arr_to_binary_file


!!if you use "bin" format, I assume you want to save disk space; so convert to single precision 
  subroutine float_to_binary_file(fp, x)
    type(file_pointer) fp
    real(dl) x
    write(fp%unit,REC=fp%pt)  real(x, IB)
    fp%pt=fp%pt +1
  end subroutine float_to_binary_file

  subroutine float_arr_to_binary_file(fp,x)
    type(file_pointer)fp
    real(dl),dimension(:)::x
    integer i
    do i=1,size(x)
       call float_to_binary_file(fp,x(i))
    enddo
  end subroutine float_arr_to_binary_file

  subroutine int_to_file(fp,i,fmt)
    type(file_pointer) fp
    integer i
    UNKNOWN_STRING,optional::fmt
    select case(fp%mode)
    case("txt") 
       if(present(fmt))then
          call int_to_txt_file(fp,i,fmt)
       else
          call int_to_txt_file(fp,i)
       endif
    case("bin")
       call int_to_binary_file(fp,i)
    case("unf")
       write(fp%unit) i
    case default
       write(*,*) "File IO mode "//fp%mode//" is not supported."
       stop
    end select
  end subroutine int_to_file

  subroutine int_arr_to_file(fp,i,fmt)
    type(file_pointer) fp
    integer,dimension(:):: i
    UNKNOWN_STRING,optional::fmt
    select case(fp%mode)
    case("txt") 
       if(present(fmt))then
          call int_arr_to_txt_file(fp,i,fmt)
       else
          call int_arr_to_txt_file(fp,i)
       endif
    case("bin")
       call int_arr_to_binary_file(fp,i)
    case("unf")
       write(fp%unit) i
    case default
       write(*,*) "File IO mode "//fp%mode//" is not supported."
    end select
  end subroutine int_arr_to_file


  subroutine float_to_file(fp,x,fmt)
    type(file_pointer) fp
    real(dl) x
    UNKNOWN_STRING,optional:: fmt
    select case(fp%mode)
    case("txt") 
       if(present(fmt))then
          call float_to_txt_file(fp,x,fmt)
       else
          call float_to_txt_file(fp,x)
       endif
    case("bin")
       call float_to_binary_file(fp,x)
    case("unf")
       write(fp%unit) x
    case default
       write(*,*) "File IO mode "//fp%mode//" is not supported."
    end select
  end subroutine float_to_file

  subroutine float_arr_to_file(fp,x,fmt)
    type(file_pointer) fp
    real(dl),dimension(:):: x
    UNKNOWN_STRING,optional::fmt
    select case(fp%mode)
    case("txt") 
       if(present(fmt))then
          call float_arr_to_txt_file(fp,x,fmt)
       else
          call float_arr_to_txt_file(fp,x)
       endif
    case("bin")
       call float_arr_to_binary_file(fp,x)
    case("unf")
       write(fp%unit) x
    case default
       write(*,*) "File IO mode "//fp%mode//" is not supported."
    end select
  end subroutine float_arr_to_file

  subroutine int_from_binary_file(fp, i)
    type(file_pointer) fp
    integer i
    read(fp%unit, REC=fp%pt) i
    fp%pt=fp%pt+1
  end subroutine int_from_binary_file

  subroutine int_from_txt_file(fp,i)
    type(file_pointer) fp
    integer i
    read(fp%unit,*) i
  end subroutine int_from_txt_file
 
  subroutine int_arr_from_binary_file(fp,x)
    type(file_pointer) fp
    integer,dimension(:):: x
    integer i
    do i=1,size(x)
       call int_from_binary_file(fp,x(i))
    enddo
  end subroutine int_arr_from_binary_file

  subroutine int_arr_from_txt_file(fp,x)
    type(file_pointer) fp
    integer,dimension(:):: x
    read(fp%unit,*) x
  end subroutine int_arr_from_txt_file

  subroutine float_from_binary_file(fp, x)
    type(file_pointer) fp
    real(dl) x
    real(IB) xapp
    read(fp%unit, REC=fp%pt) xapp
    x = xapp
    fp%pt=fp%pt+1
  end subroutine float_from_binary_file

  subroutine float_from_txt_file(fp,x)
    type(file_pointer) fp
    real(dl) x
    read(fp%unit,*) x
  end subroutine float_from_txt_file

  subroutine float_arr_from_txt_file(fp,x)
    type(file_pointer) fp
    real(dl),dimension(:):: x
    read(fp%unit,*) x    
  end subroutine float_arr_from_txt_file

  subroutine float_arr_from_binary_file(fp,x)
    type(file_pointer) fp
    real(dl),dimension(:) :: x
    integer i
    do i=1,size(x)
       call float_from_binary_file(fp,x(i))
    enddo
  end subroutine float_arr_from_binary_file

  subroutine int_from_file(fp,i)
    type(file_pointer) fp
    integer i
    select case(fp%mode)
    case("txt")
       call int_from_txt_file(fp,i)
    case("bin")
       call int_from_binary_file(fp,i)
    case("unf")
       read(fp%unit) i
    case default
       write(*,*) "File IO mode "//fp%mode//" is not supported."
       stop
    end select
  end subroutine int_from_file

  subroutine int_arr_from_file(fp,i)
    type(file_pointer) fp
    integer,dimension(:):: i
    select case(fp%mode)
    case("txt")
       call int_arr_from_txt_file(fp,i)
    case("bin")
       call int_arr_from_binary_file(fp,i)
    case("unf")
       read(fp%unit) i
    case default
       write(*,*) "File IO mode "//fp%mode//" is not supported."
       stop
    end select
  end subroutine int_arr_from_file

  subroutine float_from_file(fp,x)
    type(file_pointer) fp
    real(dl) x
    select case(fp%mode)
    case("txt")
       call float_from_txt_file(fp,x)
    case("bin")
       call float_from_binary_file(fp,x)
    case("unf")
       read(fp%unit) x
    case default
       write(*,*) "File IO mode "//fp%mode//" is not supported."
       stop
    end select   
  end subroutine float_from_file

  subroutine float_arr_from_file(fp,x)
    type(file_pointer) fp
    real(dl),dimension(:):: x
    select case(fp%mode)
    case("txt")
       call float_arr_from_txt_file(fp,x)
    case("bin")
       call float_arr_from_binary_file(fp,x)
    case("unf")
       read(fp%unit) x
    case default
       write(*,*) "File IO mode "//fp%mode//" is not supported."
       stop
    end select
  end subroutine float_arr_from_file

  subroutine close_file(fp)
    type(file_pointer) fp
    if( fp%unit.ge.file_min_unit .and. fp%unit.le.file_max_unit)then
       close(fp%unit)
       unit_opened(fp%unit)=.false.
       fp%pt = 1
       fp%unit = 0
       fp%path = ""
    else
       write(*,*) "can not close unknow file unit: "//trim(int2str(fp%unit))
       stop
    endif    
  end subroutine close_file


  FUNCTION File_Exists(FileName)
    UNKNOWN_STRING, INTENT(IN) :: FileName
    LOGICAL File_Exists
    INQUIRE(FILE=trim(FileName), EXIST = File_Exists)
  END FUNCTION File_Exists




  subroutine File_OverWrite_Query(FileName)
    UNKNOWN_STRING, INTENT(IN) :: FileName
    character(LEN=3) ans
    if(File_Exists(trim(Filename)))then
       write(*,*) "The file "//trim(FileName)//" already exists"
       write(*,*) "Do you want to overwrite it (yes/no)?"
       read(*,*) ans
       if(ans.ne."yes" .and. ans.ne."YES") stop       
    endif
  end subroutine File_OverWrite_Query

  Function File_NumColumns(FileName) result(NumColumns)
    UNKNOWN_STRING FileName
    LONG_STRING FirstLine
    Integer NumColumns,I,Length,J
    type(file_pointer) fp
    NumColumns=0
    if(.not. File_Exists(filename)) return
    fp = open_file(trim(filename), "r")
    if( file_readOneline_string(fp, firstline) )then
       Length=scan(FirstLine,"1234567890+-.eE",BACK=.true.)+1
       FirstLine(Length:Length) = " "
       i = 1
       do while(i.lt. length)
          NumColumns=NumColumns+1
          J=verify(FirstLine(i:Length),"1234567890-+.eE")
          I = I + J
          J=scan(FirstLine(I:Length),"1234567890-+.eE")
          if(J.eq.0) goto 100
          I = I + J - 1
       enddo
    endif
100 call close_file(fp)
  End Function File_NumColumns

  Function File_NumLines(filename) result(Lines) !! # comment, empty lines are skipped
    UNKNOWN_STRING, INTENT(IN)::filename
    LONG_STRING inline
    integer Lines
    Type(File_pointer) fp
    Lines = 0
    if(.not. File_Exists(filename)) return
    fp = open_file(filename, "r")
    do while(file_readoneline_string(fp, inline))
       lines = lines + 1
    enddo
100 call close_file(fp)
  End Function File_NumLines


  Function File_ReadOneLine_Int(fp, param) result(success)
    integer param
    logical success
    type(file_pointer) fp
    LONG_STRING line
    Line = ""
    success = .false.
    do while(Line(1:1) .eq. text_comment_symbol .or. Trim(Line) .eq. "")
       read(fp%unit, '(A)', End=200, Err=200) line
       line = adjustl(trim(line))
    enddo
    read(Line, *, Err = 200, End=200) param
    success = .true.
200 return
  End Function File_ReadOneLine_Int

  Function File_ReadOneLine_Double(fp, param) result(success)
    real(dl) param
    logical success
    type(file_pointer) fp
    LONG_STRING line
    Line = ""
    success = .false.
    do while(Line(1:1) .eq. text_comment_symbol .or. Trim(Line) .eq. "")
       read(fp%unit, '(A)', End=200, Err=200) line
       line = adjustl(trim(line))
    enddo
    read(Line, *, Err = 200, End=200) param
    success = .true.
200 return
  End Function File_ReadOneLine_Double


  Function File_ReadOneLine_Int_arr(fp, params) result(success)
    integer,dimension(:)::params
    logical success
    type(file_pointer) fp
    LONG_STRING line
    Line = ""
    success = .false.
    do while(Line(1:1) .eq. text_comment_symbol .or. Trim(Line) .eq. "")
       read(fp%unit, '(A)', End=200, Err=200) line
       line = adjustl(trim(line))
    enddo
    read(Line, *, Err = 200, End=200) params
    success = .true.
200 return
  End Function File_ReadOneLine_Int_arr


  Function File_ReadOneLine_Double_arr(fp, params) result(success)
    real(dl),dimension(:)::params
    logical success
    type(file_pointer) fp
    LONG_STRING line
    Line = ""
    success = .false.
    do while(Line(1:1) .eq. text_comment_symbol .or. Trim(Line) .eq. "")
       read(fp%unit, '(A)', End=200, Err=200) line
       line = adjustl(trim(line))
    enddo
    read(Line, *, Err = 200, End=200) params
    success = .true.
200 return
  End Function File_ReadOneLine_Double_arr

  Function File_ReadOneLine_String(fp, line) result(success)
    logical success
    type(file_pointer) fp
    UNKNOWN_STRING line
    Line = ""
    success = .false.
    do while(Line(1:1) .eq. text_comment_symbol .or. Trim(Line) .eq. "")
       read(fp%unit, '(A)', End=200, Err=200) line
       line = adjustl(trim(line))
    enddo
    success = .true.
    return
200 line = ""
    return
  End Function File_ReadOneLine_String


  subroutine file_readdata(fname, nrows, ncols, f)
    UNKNOWN_STRING fname
    type(file_pointer) fp
    integer nrows, ncols
    real(dl) f(nrows, ncols), params(ncols)
    integer i
    logical success
    fp = open_file(trim(fname), 'r')
    do i=1, nrows
       success =  file_readoneline_double_arr(fp, params)
       if(.not. success) stop "file does not contain the expected data..."
       f(i,:) = params
    enddo
    call close_file(fp)
  end subroutine file_readdata

  Function File_ReadMatrix(fp, matrix) result(success)
    integer nx, ny
    logical success
    real(dl),dimension(:,:),intent(inout)::matrix
    type(file_pointer) fp
    nx = size(matrix, 1)
    ny = size(matrix, 2)
    call Read_Matrix(fp%unit, nx, ny, matrix, success)
  end Function File_ReadMatrix

  subroutine File_WriteMatrix(fp, matrix)
    real(dl) matrix(:,:)
    type(file_pointer) fp
    call Write_matrix(fp%unit, matrix)
  end subroutine File_WriteMatrix

  Subroutine File_SkipLines(fp, numlines)
    type(file_pointer) fp
    integer numlines, i
    LONG_STRING line
    do i = 1, numlines
       read(fp%unit, '(A)', End=200, Err=200) line
    enddo
200 return
  end Subroutine File_SkipLines

  subroutine file_read_csv_line()
  end subroutine file_read_csv_line

  subroutine list_character_dump(l, fname)
    DYNAMIC_STRING l
    UNKNOWN_STRING fname
    type(file_pointer) fp
    fp = open_file(trim(fname), "w")
    if(allocated(l%i1))then
       select case(l%stack)
       case(1)
          write(fp%unit,"(A)")l%i1(1:l%loc)
       case(2)
          write(fp%unit,"(A)")l%i1
          write(fp%unit,"(A)")l%i2(1:l%loc)
       case(3)
          write(fp%unit,"(A)")l%i1
          write(fp%unit,"(A)")l%i2
          write(fp%unit,"(A)")l%i3(1:l%loc)
       case(4)
          write(fp%unit,"(A)")l%i1
          write(fp%unit,"(A)")l%i2
          write(fp%unit,"(A)")l%i3
          write(fp%unit,"(A)")l%i4(1:l%loc)
       end select
    endif
    call close_file(fp)
  end subroutine list_character_dump


  subroutine load_dictionary(fname, dict)
    UNKNOWN_STRING fname
    type(dictionary) dict
    STRING line
    type(file_pointer) fp
    integer eqloc
    fp = open_file(trim(fname), "r")
    do while(file_readoneline_string(fp, line))
       eqloc = scan(line, "=")
       if(eqloc .gt. 1 .and. eqloc .lt. len_trim(line))then
          call dictionary_insert(dict, line(1:eqloc-1), line(eqloc+1:))
       endif
    enddo
    call close_file(fp)
  end subroutine load_dictionary


end module file_io_utils
