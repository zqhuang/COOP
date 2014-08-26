module coop_MPI_mod
  use coop_wrapper_typedef
  implicit none

  private

#ifdef MPI
#include "mpif.h"
#endif

#define MPI2 YES  

  integer,parameter::dl = coop_real_length

  Interface coop_MPI_Sum
     module procedure coop_MPI_Sum_single, coop_MPI_Sum_Arr, coop_MPI_Sum_matrix, coop_MPI_Sum_Lattice
  End Interface coop_MPI_Sum

  public::coop_MPI_init, coop_MPI_finalize, coop_MPI_abort, coop_MPI_rank, coop_MPI_NumProc, coop_MPI_sum, coop_MPI_Barrier, coop_MPI_Send, coop_MPI_Recv, coop_MPI_sendrecv, coop_MPI_LeftRank, coop_MPI_RightRank, coop_MPI_print, coop_MPI_echo

   
contains

  

  subroutine coop_MPI_init()
#ifdef MPI
    integer ierror
    call mpi_init(ierror)
    if (ierror .ne. MPI_SUCCESS) stop 'MPI fail: cannot initialize'
#else
    call coop_feedback("MPI not used...ignoring the MPI related operations...", 2)
#endif
  end subroutine coop_MPI_init



  subroutine coop_MPI_Abort(msg)
    character(LEN=*),optional::msg
#ifdef MPI
    integer errcode, ierror
#endif
    if(present(msg))then
       call coop_MPI_Echo(msg)
    else
       call coop_MPI_Echo("aborting")
    endif
#ifdef MPI
    call MPI_Abort(MPI_COMM_WORLD, errcode, ierror)
    if (ierror .ne. MPI_SUCCESS)then
       stop 'MPI_abort failed...'
    else
       call coop_MPI_Finalize()
    endif
#else
    stop
#endif
  End subroutine coop_MPI_Abort

  subroutine coop_MPI_Finalize()
#ifdef MPI
    integer ierror
    call MPI_Finalize(ierror)
    if (ierror .ne. MPI_SUCCESS) stop 'MPI fail: cannot finalize'
#else
    stop 
#endif    
  end subroutine coop_MPI_Finalize

  function coop_MPI_rank() result(thisrank)
    integer,save::rank = 0
    integer thisrank
#ifdef MPI
    logical,save::init=.true.
    integer ierror
    if(init)then
       call mpi_comm_rank(MPI_COMM_WORLD, rank,ierror)
       if (ierror .ne. MPI_SUCCESS) stop 'MPI fail: cannot get the rank'
       init = .false.
    endif
#endif
    thisrank = rank
    return
  end function coop_MPI_rank

  function coop_MPI_NumProc() result(num)
    integer num
    integer,save::numproc = 1
#ifdef MPI
    logical,save::init=.true.
    integer ierror
    if(init)then
       call mpi_comm_size(MPI_COMM_WORLD, numproc, ierror)
       if (ierror .ne. MPI_SUCCESS) stop 'MPI fail: cannot get the # of processors'
       init = .false.
    endif
#endif
    num = numproc
    return
  end function coop_MPI_NumProc

  function coop_MPI_RightRank()
    integer coop_MPI_RightRank
    coop_MPI_RightRank = mod(coop_MPI_Rank() + 1, coop_MPI_NumProc())
  end function coop_MPI_RightRank

  function coop_MPI_LeftRank()
    integer coop_MPI_LeftRank
    coop_MPI_LeftRank = mod(coop_MPI_Rank() + coop_MPI_NumProc() - 1, coop_MPI_NumProc())
  end function coop_MPI_LeftRank

  subroutine coop_MPI_Print(str)
    character(LEN=*) str
    if(coop_MPI_Rank().eq.0) write(*,*) str
  end subroutine coop_MPI_Print

  subroutine coop_MPI_Echo(str)
    character(LEN=*) str
    write(*,'(A5, I3, A)') "Rank#", coop_MPI_Rank(), ": "//trim(str)
  end subroutine coop_MPI_Echo

  subroutine coop_MPI_Send(buf, count, dest, tag)
    real(dl) buf
    integer count, dest
    integer,optional::tag
#ifdef MPI
    integer  ierror
    if(present(tag))then
       call MPI_Send(buf, count, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierror)
    else
       call MPI_Send(buf, count, MPI_DOUBLE_PRECISION, dest, coop_MPI_Rank()*coop_MPI_NumProc() + dest, MPI_COMM_WORLD, ierror)
    endif
    if(ierror .ne. MPI_SUCCESS) call coop_MPI_Abort("MPI_Send failed")
#else
    if(dest.ne.0 .or. count .lt.0) call coop_MPI_Abort("MPI_Send wrong input")
#endif
  end subroutine coop_MPI_Send

  subroutine coop_MPI_Recv(buf, count, source, tag)
    real(dl) buf
    integer count, source
    integer, optional::tag
#ifdef MPI
    integer status(MPI_STATUS_SIZE), ierror
    if(present(tag))then
       call MPI_Recv(buf, count, MPI_DOUBLE_PRECISION, source, tag, MPI_COMM_WORLD, status, ierror)
    else
       call MPI_Recv(buf, count, MPI_DOUBLE_PRECISION, source, source*coop_MPI_NumProc() + coop_MPI_Rank(), MPI_COMM_WORLD, status, ierror)
    endif
    if(ierror .ne. MPI_SUCCESS) call coop_MPI_Abort("MPI_Recv failed")
#else
    if(source .ne. 0) call coop_MPI_Abort("MPI wrong input")
#endif    
  end subroutine coop_MPI_Recv

  subroutine coop_MPI_SendRecv(sendbuf, sendcount, dest, sendtag, recvbuf, recvcount, source, recvtag)
    real(dl) sendbuf, recvbuf
    integer sendcount, dest, sendtag, recvcount, source, recvtag
#ifdef MPI
    integer status(MPI_STATUS_SIZE), ierror
    call MPI_SendRecv(sendbuf, sendcount, MPI_DOUBLE_PRECISION, dest, sendtag, recvbuf, recvcount , MPI_DOUBLE_PRECISION, source, recvtag, MPI_COMM_WORLD, status, ierror)   
    if(ierror .ne. MPI_SUCCESS) call coop_MPI_Abort("MPI_SendRecv failed")
#else
    if( dest .ne. 0) call coop_MPI_Abort("MPI_SendRecv wrong input!")
#endif
  end subroutine coop_MPI_SendRecv


  subroutine coop_MPI_Sum_Single(buf)
    real(dl) buf
#ifdef MPI
    integer ierror
#ifdef MPI2
    call MPI_Allreduce(MPI_IN_PLACE, buf, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
#else
    real(dl) sumbuf
    call MPI_Allreduce(buf, sumbuf, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
    buf = sumbuf
#endif
    if(ierror .ne. MPI_SUCCESS) call coop_MPI_Abort("coop_MPI_Sum_Single failed")
#endif
  end subroutine coop_MPI_Sum_Single

  subroutine coop_MPI_Sum_Arr(buf)
    real(dl),dimension(:),intent(inout)::buf
#ifdef MPI
    integer ierror
#ifdef MPI2
    call MPI_Allreduce(MPI_IN_PLACE, buf,  size(buf), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
#else
    real(dl) sumbuf(size(buf))
    call MPI_Allreduce(buf, sumbuf,  size(buf), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
    buf = sumbuf
#endif
    if(ierror .ne. MPI_SUCCESS) call coop_MPI_Abort("coop_MPI_Sum_Arr failed")    
#endif
  end subroutine coop_MPI_Sum_Arr

  subroutine coop_MPI_Sum_Matrix(buf)
    real(dl),dimension(:,:),intent(inout)::buf
#ifdef MPI
    integer ierror
#ifdef MPI2
    call MPI_Allreduce(MPI_IN_PLACE, buf, size(buf,1)*size(buf,2), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
#else
    real(dl) sumbuf(size(buf,1), size(buf,2))
    call MPI_Allreduce(buf,sumbuf, size(buf,1)*size(buf,2), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
    buf = sumbuf
#endif
    if(ierror .ne. MPI_SUCCESS) call coop_MPI_Abort("coop_MPI_Sum_Matrix failed")    
#endif
  end subroutine coop_MPI_Sum_Matrix

  subroutine coop_MPI_Sum_Lattice(buf)
    real(dl),dimension(:,:,:),intent(inout)::buf
#ifdef MPI
    integer ierror
#ifdef MPI2    
    call MPI_Allreduce(MPI_IN_PLACE, buf, size(buf,1)*size(buf,2)*size(buf,3), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
#else
    real(dl),dimension(size(buf,1),size(buf,2), size(buf,3))::sumbuf
    call MPI_Allreduce(buf, sumbuf, size(buf,1)*size(buf,2)*size(buf,3), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
    buf = sumbuf
#endif
    if(ierror .ne. MPI_SUCCESS) call coop_MPI_Abort("coop_MPI_Sum_Lattice failed")    
#endif
  end subroutine coop_MPI_Sum_Lattice

  subroutine coop_MPI_Barrier()
#ifdef MPI
    integer ierror
    call MPI_Barrier(MPI_COMM_WORLD, ierror)
    if(ierror .ne. MPI_SUCCESS) call coop_MPI_Abort("coop_MPI_Barrier failed")
#endif
  end subroutine coop_MPI_Barrier

end module Coop_MPI_mod
