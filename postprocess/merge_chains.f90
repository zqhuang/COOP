program getdist
  use coop_wrapper_utils
  use coop_inifile_mod
  implicit none
#include "constants.h"
  COOP_INT,parameter::min_lines = 50
  COOP_STRING root1, root2, root_out, tmpstr
  COOP_LONG_STRING:: line
  type(coop_file)::fp1, fp2, fp
  COOP_INT::row1, col1, row2, col2, ind, i,perc, j
  logical::copy_config
  root1 = coop_InputArgs(1)
  if(trim(coop_file_postfix_of(root1)).eq.'ini')then     
     call fp%open_skip_comments(root1)
     if(.not. fp%read_string(root_out)) call coop_return_error("MergeChains", "You need to set root_out in the 1st line of ini file.", "stop")
     if(.not. fp%read_int(perc)) call coop_return_error("MergeChains", "You need to set ignore_percent (20) in the 2nd line of ini file.", "stop")
     if(.not. fp%read_int(ind)) call coop_return_error("MergeChains", "You need to set maximum index (8) in the 3rd line of ini file.", "stop")
     copy_config = .true.     
     do while( fp%read_string(root2))
        if(trim(root2).eq.trim(root_out)) call coop_return_error("MergeChains", "root to be merged is the same as root_output", "stop")        
        if(copy_config)then
           call system("cp "//trim(root2)//".inputparams "//trim(root_out)//".inputparams")
           call system("cp "//trim(root2)//".ranges "//trim(root_out)//".ranges")
           call system("cp "//trim(root2)//".paramnames "//trim(root_out)//".paramnames")
           copy_config = .false.
        endif
        do i=1, ind
           row2 = coop_file_numlines(trim(root2)//"_"//COOP_STR_OF(i)//".txt")
           if(row2.le. min_lines)cycle
           call fp1%open(trim(root_out)//"_"//COOP_STR_OF(i)//".txt", "a")
           call fp2%open(trim(root2)//"_"//COOP_STR_OF(i)//".txt")
           do j=1, max(row2*perc/100, min_lines)
              read(fp2%unit, "(a)")line
           enddo
           do j= max(row2*perc/100, min_lines)+1, row2
              read(fp2%unit, "(a)")line
              write(fp1%unit, "(a)") trim(line)
           enddo
           call fp1%close()
           call fp2%close()
        enddo
     enddo
     call fp%close()
  else
     if(trim(root1).eq.'')stop "./MergeChains chainname1 chainname2 output_name [num_chains] [discard_percent]"
     root2 = coop_InputArgs(2)
     if(trim(root2).eq.'' .or. trim(root2).eq.trim(root1))stop "./MergeChains chainname1 chainname2 [numchains]"
     root_out = coop_InputArgs(3)
     if(trim(root_out).eq.'' .or. trim(root_out).eq.trim(root1)  .or. trim(root_out) .eq. trim(root2) )stop "./MergeChains chainname1 chainname2 [numchains]"  
     tmpstr = coop_InputArgs(4)
     if(trim(tmpstr).eq.'')then
        ind = 8
     else
        read(tmpstr, *) ind
     endif
     tmpstr = coop_InputArgs(5)
     if(trim(tmpstr).eq.'')then
        perc = 20
     else
        read(tmpstr, *) perc
     endif

     if(trim(root_out) .eq. "overwrite")then
        do i=1, ind
           row2 = coop_file_numlines(trim(root2)//"_"//COOP_STR_OF(i)//".txt")
           if(row2 .le. min_lines)cycle  !!ignore small files           
           call fp1%open(trim(root1)//"_"//COOP_STR_OF(i)//".txt", "a")
           call fp2%open(trim(root2)//"_"//COOP_STR_OF(i)//".txt")
           do j=1, max(row2*perc/100, min_lines)
              read(fp2%unit, "(a)")line
           enddo
           do j= max(row2*perc/100, min_lines)+1, row2
              read(fp2%unit, "(a)")line
              write(fp1%unit, "(a)") trim(line)
           enddo
           call fp1%close()
           call fp2%close()
        enddo

     else
        call system("cp "//trim(root1)//".inputparams "//trim(root_out)//".inputparams")
        call system("cp "//trim(root1)//".ranges "//trim(root_out)//".ranges")
        call system("cp "//trim(root1)//".paramnames "//trim(root_out)//".paramnames")

        do i=1, ind
           row2 = coop_file_numlines(trim(root2)//"_"//COOP_STR_OF(i)//".txt")
           if(row2 .le. min_lines)cycle  !!ignore small files
           call system("cp "//trim(root1)//"_"//COOP_STR_OF(i)//".txt "//trim(root_out)//"_"//COOP_STR_OF(i)//".txt")
           call fp1%open(trim(root_out)//"_"//COOP_STR_OF(i)//".txt", "a")
           call fp2%open(trim(root2)//"_"//COOP_STR_OF(i)//".txt")
           do j=1, max(row2*perc/100, min_lines)
              read(fp2%unit, "(a)")line
           enddo
           do j= max(row2*perc/100, min_lines)+1, row2
              read(fp2%unit, "(a)")line
              write(fp1%unit, "(a)") trim(line)
           enddo
           call fp1%close()
           call fp2%close()
        enddo
     endif
  endif
end program getdist


