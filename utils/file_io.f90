module coop_file_mod
  use coop_wrapper_typedef
  use coop_list_mod
  implicit none
#include "constants.h"


private

public::coop_file, coop_copy_file, coop_delete_file, coop_create_file, coop_create_directory, coop_delete_directory, coop_file_numcolumns, coop_file_numlines, coop_load_dictionary, coop_free_file_unit, coop_file_exists

  character,parameter::text_comment_symbol = "#"


  Type coop_file
     COOP_INT unit
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

  COOP_INT,parameter::coop_file_unit_min = 11
  COOP_INT,parameter::coop_file_unit_max = 99

  interface coop_file
     procedure coop_open_file
  end interface coop_file

contains
  
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
    call system("mkdir "//trim(directory))
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
    if(trim(inline) .eq. "" .or. inline(1:1).eq.text_comment_symbol)then
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
    LOGICAL File_Exists
    INQUIRE(FILE=trim(FileName), EXIST = File_Exists)
  end function coop_file_exists


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
100 call fp%close()
  End Function Coop_File_NumColumns

  Function coop_File_NumLines(filename) result(Lines) !! # comment, empty lines are skipped
    COOP_UNKNOWN_STRING, INTENT(IN)::filename
    COOP_LONG_STRING inline
    COOP_INT Lines
    Type(Coop_file) fp
    Lines = 0
    if(.not. coop_File_Exists(filename)) return
    call fp%open(filename, "r")
    do while(coop_file_readline_string(fp, inline))
       lines = lines + 1
    enddo
100 call fp%close()
  End Function Coop_File_NumLines


  Function Coop_file_readline_Int(fp, param) result(success)
    COOP_INT param
    logical success
    class(coop_file) fp
    COOP_LONG_STRING line
    Line = ""
    success = .false.
    do while(Line(1:1) .eq. text_comment_symbol .or. Trim(Line) .eq. "")
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
    do while(Line(1:1) .eq. text_comment_symbol .or. Trim(Line) .eq. "")
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
    do while(Line(1:1) .eq. text_comment_symbol .or. Trim(Line) .eq. "")
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
    do while(Line(1:1) .eq. text_comment_symbol .or. Trim(Line) .eq. "")
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
    do while(Line(1:1) .eq. text_comment_symbol .or. Trim(Line) .eq. "")
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
    COOP_INT numlines, i
    COOP_LONG_STRING line
    do i = 1, numlines
       read(fp%unit, '(A)', End=200, Err=200) line
    enddo
200 return
  end Subroutine Coop_File_SkipLines


  subroutine coop_load_dictionary(fname, dict)
    COOP_UNKNOWN_STRING fname
    type(coop_dictionary) dict
    COOP_STRING line
    type(coop_file) fp
    COOP_INT eqloc
    call fp%open(trim(fname), "r")
    do while(fp%read_string(line))
       eqloc = scan(line, "=")
       if(eqloc .gt. 1 .and. eqloc .lt. len_trim(line))then
          call dict%insert(line(1:eqloc-1), line(eqloc+1:))
       endif
    enddo
    call fp%close()
  end subroutine coop_load_dictionary


end module coop_file_mod
