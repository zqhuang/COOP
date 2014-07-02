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

  integer,parameter::dlc = kind( (1.d0, 1.d0) )
  integer,parameter::dl = kind(1.d0)

contains

  subroutine coop_fft_forward_1d(n, fx, fk)
    COOP_INT n
    real(dl) fx(n)
    complex(dlc) fk(n/2+1)
#ifdef HAS_FFTW
    call fft_1d_forward(n, fx, fk)
#else
    write(*,*) "coop_fft: Cannot find FFTW library. Please change the configura.in file."
#endif    
  end subroutine coop_fft_forward_1d

  subroutine coop_fft_backward_1d(n, fk, fx)
    COOP_INT n
    real(dl) fx(n)
    complex(dlc) fk(n/2+1)
#ifdef HAS_FFTW
    call fft_1d_backward(n, fk, fx)
    fx = fx/n
#else
    write(*,*) "coop_fft: Cannot find FFTW library. Please change the configura.in file."
#endif    
  end subroutine coop_fft_backward_1d


  subroutine coop_fft_forward_2d(ny, nx, fx, fk)
    COOP_INT nx, ny
    real(dl) fx(ny,nx)
    complex(dlc) fk(ny/2+1,nx)
#ifdef HAS_FFTW
    call fft_2d_forward(nx, ny, fx, fk)
#else
    write(*,*) "coop_fft: Cannot find FFTW library. Please change the configura.in file."
#endif    
  end subroutine coop_fft_forward_2d

  subroutine coop_fft_backward_2d(ny, nx, fk, fx)
    COOP_INT nx, ny
    real(dl) fx(ny, nx)
    complex(dlc) fk(ny/2+1, nx)
#ifdef HAS_FFTW
    call fft_2d_backward(nx, ny, fk, fx)
    fx = fx/(dble(nx)*dble(ny))
#else
    write(*,*) "coop_fft: Cannot find FFTW library. Please change the configura.in file."
#endif    
  end subroutine coop_fft_backward_2d


  subroutine coop_fft_forward_2d_ss(ny, nx, fx, fk)
    COOP_INT nx, ny
    real(dl) fx(ny*nx)
    complex(dlc) fk(ny/2+1, nx)
#ifdef HAS_FFTW
    call fft_2d_forward(nx, ny, fx, fk)
#else
    write(*,*) "coop_fft: Cannot find FFTW library. Please change the configura.in file."
#endif    
  end subroutine coop_fft_forward_2d_ss

  subroutine coop_fft_backward_2d_ss(ny, nx, fk, fx)
    COOP_INT nx, ny
    real(dl) fx(ny*nx)
    complex(dlc) fk(ny/2+1, nx)
#ifdef HAS_FFTW
    call fft_2d_backward(nx, ny, fk, fx)
    fx = fx/(dble(nx)*dble(ny))
#else
    write(*,*) "coop_fft: Cannot find FFTW library. Please change the configura.in file."
#endif    
  end subroutine coop_fft_backward_2d_ss


  subroutine coop_fft_forward_3d(nz, ny, nx, fx, fk)
    COOP_INT nz, ny, nx
    real(dl) fx(nz, ny, nx)
    complex(dlc) fk(nz/2+1, ny, nx)
#ifdef HAS_FFTW
    call fft_3d_forward(nx, ny, nz, fx, fk)
#else
    write(*,*) "coop_fft: Cannot find FFTW library. Please change the configura.in file."
#endif    
  end subroutine coop_fft_forward_3d

  subroutine coop_fft_backward_3d(nz, ny, nx, fk, fx)
    COOP_INT nz, ny, nx
    real(dl) fx(nz, ny, nx)
    complex(dlc) fk(nz/2+1, ny, nx)
#ifdef HAS_FFTW
    call fft_3d_backward(nx, ny, nz, fk, fx)
    fx = fx/(dble(nx)*dble(ny)*dble(nz))
#else
    write(*,*) "coop_fft: Cannot find FFTW library. Please change the configura.in file."
#endif    
  end subroutine coop_fft_backward_3d



  subroutine coop_fft_forward_3d_ss(nz, ny, nx, fx, fk)
    COOP_INT nz, ny, nx
    real(dl) fx(nz* ny*nx)
    complex(dlc) fk(nz/2+1, ny, nx)
#ifdef HAS_FFTW
    call fft_3d_forward(nx, ny, nz, fx, fk)
#else
    write(*,*) "coop_fft: Cannot find FFTW library. Please change the configura.in file."
#endif    
  end subroutine coop_fft_forward_3d_ss

  subroutine coop_fft_backward_3d_ss(nz, ny, nx, fk, fx)
    COOP_INT nz, ny, nx
    real(dl) fx(nz*ny*nx)
    complex(dlc) fk(nz/2+1, ny, nx)
#ifdef HAS_FFTW
    call fft_3d_backward(nx, ny, nz, fk, fx)
    fx = fx/(dble(nx)*dble(ny)*dble(nz))
#else
    write(*,*) "coop_fft: Cannot find FFTW library. Please change the configura.in file."
#endif    
  end subroutine coop_fft_backward_3d_ss

end module coop_fft_mod
