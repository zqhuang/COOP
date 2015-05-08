    character(LEN=2)::open_mode
    if(trim(adjustl(filename)).eq."") stop "coop_file_open: file name cannot be a null string"    
    fp%pt = 1
    if(present(mode))then
       open_mode = adjustl(trim(mode))
    else
       open_mode = 'w'
    endif
    fp%unit = coop_free_file_unit()
    select case(trim(open_mode))
    case("w","W")
       open(fp%unit,FILE=trim(adjustl(filename)),FORM="FORMATTED",STATUS="UNKNOWN",ACCESS='SEQUENTIAL',ACTION="READWRITE", Err=200)
       fp%mode = 'txt'
    case("wq","WQ", "qw", "QW")
       call coop_file_overwrite_query(trim(adjustl(filename)))
       open(fp%unit,FILE=trim(adjustl(filename)),FORM="FORMATTED",STATUS="UNKNOWN",ACCESS='SEQUENTIAL',ACTION="READWRITE", Err=200)
       fp%mode = 'txt'
    case("r", "R")
       open(fp%unit,FILE=trim(adjustl(filename)),FORM="FORMATTED",STATUS="UNKNOWN", ACCESS='SEQUENTIAL', ACTION="READ", Err=200)
       fp%mode = 'txt'
    case("a", "A")
       open(fp%unit,FILE=trim(adjustl(filename)),FORM="FORMATTED",STATUS="UNKNOWN", ACCESS='APPEND',ACTION="READWRITE", Err=200)
       fp%mode = 'txt'
    case("b", "B")
       if(present(recl))then
          open(UNIT=fp%unit,FILE=trim(adjustl(filename)),FORM="UNFORMATTED",STATUS="UNKNOWN", ACCESS="DIRECT",ACTION="READWRITE",RECL=recl,ERR=200) 
          fp%mode = 'bin'
       else
          open(UNIT=fp%unit,FILE=trim(adjustl(filename)),FORM="UNFORMATTED",STATUS="UNKNOWN", ACCESS="DIRECT",ACTION="READWRITE",RECL=coop_integer_length,ERR=200) 
          fp%mode = 'bin'
       endif
    case("rb", "RB", "br", "BR")
       if(present(recl))then
          open(UNIT=fp%unit,FILE=trim(adjustl(filename)),FORM="UNFORMATTED",STATUS="UNKNOWN", ACCESS="DIRECT",ACTION="READ",RECL=recl, ERR=200) 
          fp%mode = 'bin'
       else
          open(UNIT=fp%unit,FILE=trim(adjustl(filename)),FORM="UNFORMATTED",STATUS="UNKNOWN", ACCESS="DIRECT",ACTION="READ",RECL=coop_integer_length, ERR=200) 
          fp%mode = 'bin'
       endif
    case('u', "U")
       open(UNIT=fp%unit,FILE=trim(adjustl(filename)),FORM="UNFORMATTED",STATUS="UNKNOWN", ACCESS="SEQUENTIAL",ACTION="READWRITE",ERR=200)        
       fp%mode = "unf"
       fp%mode = "unf"       
    case('uq', "UQ", "qu", "QU")
       call coop_file_overwrite_query(trim(adjustl(filename)))
       open(UNIT=fp%unit,FILE=trim(adjustl(filename)),FORM="UNFORMATTED",STATUS="UNKNOWN", ACCESS="SEQUENTIAL",ACTION="READWRITE",ERR=200)        
       fp%mode = "unf"
    case('ru', "RU", "ur", "UR")
       open(UNIT=fp%unit,FILE=trim(adjustl(filename)),FORM="UNFORMATTED",STATUS="UNKNOWN", ACCESS="SEQUENTIAL",ACTION="READ",ERR=200)        
       fp%mode = "unf"
    case('AU', "au", "ua", "UA")
       open(UNIT=fp%unit,FILE=trim(adjustl(filename)),FORM="UNFORMATTED",STATUS="UNKNOWN", ACCESS="APPEND",ACTION="READWRITE",ERR=200)        
       fp%mode = "unf"
    case default
       write(*,*) "Error in fp: unknow mode "//trim(mode)
       goto 200
    end select
    fp%path = trim(adjustl(filename))
    return
200 write(*,*) "Error: can not open file "//trim(adjustl(filename))
    stop
