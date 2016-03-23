module coop_fft_mod
  use coop_basicutils_mod
  implicit none
#include "constants.h"
  private
  public::coop_fft_forward, coop_fft_backward

  interface coop_fft_forward
     module procedure::coop_fft_forward_1d, coop_fft_forward_2d, coop_fft_forward_3d, coop_fft_forward_2d_ss, coop_fft_forward_3d_ss
  end interface coop_fft_forward

  interface coop_fft_backward
     module procedure::coop_fft_backward_1d, coop_fft_backward_2d, coop_fft_backward_3d, coop_fft_backward_2d_ss, coop_fft_backward_3d_ss
  end interface coop_fft_backward

contains

  subroutine coop_fft_forward_1d(n, fx, fk)
    COOP_INT n
    COOP_REAL fx(n)
    COOP_COMPLEX fk(n/2+1)
#if HAS_FFTW
    call fft_1d_forward(n, fx, fk)
#else
    call coop_return_error("coop_fft_forward_1d", "coop_fft: Cannot find FFTW library. Please change the configura.in file.", "stop")
#endif    
  end subroutine coop_fft_forward_1d

  subroutine coop_fft_backward_1d(n, fk, fx, normalize)
    COOP_INT n
    COOP_REAL fx(n)
    COOP_COMPLEX fk(n/2+1)
    logical, optional::normalize
#if HAS_FFTW
    call fft_1d_backward(n, fk, fx)
    if(present(normalize))then
       if(normalize)then
          fx = fx/n
       endif
    else
       fx = fx/n
    endif
#else
    call coop_return_error("coop_fft", "Cannot find FFTW library. Please change the configura.in file.", "stop")
#endif    
  end subroutine coop_fft_backward_1d

  subroutine coop_fft_forward_2d(ny, nx, fx, fk)
    COOP_INT nx, ny
    COOP_REAL fx(ny,nx)
    COOP_COMPLEX fk(ny/2+1,nx)
#if HAS_FFTW
    call fft_2d_forward(nx, ny, fx, fk)
#else
    call coop_return_error("coop_fft", "Cannot find FFTW library. Please change the configura.in file.", "stop")
#endif    
  end subroutine coop_fft_forward_2d

  subroutine coop_fft_backward_2d(ny, nx, fk, fx, normalize)
    COOP_INT nx, ny
    COOP_REAL fx(ny, nx)
    COOP_COMPLEX fk(ny/2+1, nx)
    logical, optional::normalize
#if HAS_FFTW
    call fft_2d_backward(nx, ny, fk, fx)
    if(present(normalize))then
       if(normalize)then
          fx = fx/(dble(nx)*dble(ny))
       endif
    else
       fx = fx/(dble(nx)*dble(ny))
    endif
#else
    call coop_return_error("coop_fft", "Cannot find FFTW library. Please change the configura.in file.", "stop")
#endif    
  end subroutine coop_fft_backward_2d


  subroutine coop_fft_forward_2d_ss(ny, nx, fx, fk)
    COOP_INT nx, ny
    COOP_REAL fx(ny*nx)
    COOP_COMPLEX fk(ny/2+1, nx)
#if HAS_FFTW
    call fft_2d_forward(nx, ny, fx, fk)
#else
    call coop_return_error("coop_fft", "Cannot find FFTW library. Please change the configura.in file.", "stop")
#endif    
  end subroutine coop_fft_forward_2d_ss

  subroutine coop_fft_backward_2d_ss(ny, nx, fk, fx, normalize)
    COOP_INT nx, ny
    COOP_REAL fx(ny*nx)
    COOP_COMPLEX fk(ny/2+1, nx)
    logical, optional::normalize
#if HAS_FFTW
    call fft_2d_backward(nx, ny, fk, fx)
    if(present(normalize))then
       if(normalize)fx = fx/(dble(nx)*dble(ny))
    else
       fx = fx/(dble(nx)*dble(ny))
    endif
#else
    call coop_return_error("coop_fft", "Cannot find FFTW library. Please change the configura.in file.", "stop")
#endif    
  end subroutine coop_fft_backward_2d_ss


  subroutine coop_fft_forward_3d(nz, ny, nx, fx, fk)
    COOP_INT nz, ny, nx
    COOP_REAL fx(nz, ny, nx)
    COOP_COMPLEX fk(nz/2+1, ny, nx)
#if HAS_FFTW
    call fft_3d_forward(nx, ny, nz, fx, fk)
#else
    call coop_return_error("coop_fft", "Cannot find FFTW library. Please change the configura.in file.", "stop")
#endif    
  end subroutine coop_fft_forward_3d

  subroutine coop_fft_backward_3d(nz, ny, nx, fk, fx, normalize)
    COOP_INT nz, ny, nx
    COOP_REAL fx(nz, ny, nx)
    COOP_COMPLEX fk(nz/2+1, ny, nx)
    logical,optional::normalize
#if HAS_FFTW
    call fft_3d_backward(nx, ny, nz, fk, fx)
    if(present(normalize))then
       if(normalize)fx = fx/(dble(nx)*dble(ny)*dble(nz))
    else
       fx = fx/(dble(nx)*dble(ny)*dble(nz))
    endif
#else
    call coop_return_error("coop_fft", "Cannot find FFTW library. Please change the configura.in file.", "stop")
#endif    
  end subroutine coop_fft_backward_3d



  subroutine coop_fft_forward_3d_ss(nz, ny, nx, fx, fk)
    COOP_INT nz, ny, nx
    COOP_REAL fx(nz* ny*nx)
    COOP_COMPLEX fk(nz/2+1, ny, nx)
#if HAS_FFTW
    call fft_3d_forward(nx, ny, nz, fx, fk)
#else
    call coop_return_error("coop_fft", "Cannot find FFTW library. Please change the configura.in file.", "stop")
#endif    
  end subroutine coop_fft_forward_3d_ss

  subroutine coop_fft_backward_3d_ss(nz, ny, nx, fk, fx, normalize)
    COOP_INT nz, ny, nx
    COOP_REAL fx(nz*ny*nx)
    COOP_COMPLEX fk(nz/2+1, ny, nx)
    logical,optional::normalize
#if HAS_FFTW
    call fft_3d_backward(nx, ny, nz, fk, fx)
    if(present(normalize))then
       if(normalize)fx = fx/(dble(nx)*dble(ny)*dble(nz))
    else
       fx = fx/(dble(nx)*dble(ny)*dble(nz))
    endif
#else
    call coop_return_error("coop_fft", "Cannot find FFTW library. Please change the configura.in file.", "stop")
#endif    
  end subroutine coop_fft_backward_3d_ss

end module coop_fft_mod
