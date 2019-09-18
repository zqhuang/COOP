module coop_fft_mod
  use coop_basicutils_mod
  use coop_random_mod
  implicit none
#include "constants.h"
  private
  public::coop_fft_forward, coop_fft_backward, coop_generate_white_noise

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
    call coop_return_error("coop_fft_forward_1d", "coop_fft: Cannot find FFTW library. Please change the configure.in file.", "stop")
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
    call coop_return_error("coop_fft", "Cannot find FFTW library. Please change the configure.in file.", "stop")
#endif    
  end subroutine coop_fft_backward_1d

  subroutine coop_fft_forward_2d(ny, nx, fx, fk)
    COOP_INT nx, ny
    COOP_REAL fx(ny,nx)
    COOP_COMPLEX fk(ny/2+1,nx)
#if HAS_FFTW
    call fft_2d_forward(nx, ny, fx, fk)
#else
    call coop_return_error("coop_fft", "Cannot find FFTW library. Please change the configure.in file.", "stop")
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
    call coop_return_error("coop_fft", "Cannot find FFTW library. Please change the configure.in file.", "stop")
#endif    
  end subroutine coop_fft_backward_2d


  subroutine coop_fft_forward_2d_ss(ny, nx, fx, fk)
    COOP_INT nx, ny
    COOP_REAL fx(ny*nx)
    COOP_COMPLEX fk(ny/2+1, nx)
#if HAS_FFTW
    call fft_2d_forward(nx, ny, fx, fk)
#else
    call coop_return_error("coop_fft", "Cannot find FFTW library. Please change the configure.in file.", "stop")
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
    call coop_return_error("coop_fft", "Cannot find FFTW library. Please change the configure.in file.", "stop")
#endif    
  end subroutine coop_fft_backward_2d_ss


  subroutine coop_fft_forward_3d(nz, ny, nx, fx, fk)
    COOP_INT nz, ny, nx
    COOP_REAL fx(nz, ny, nx)
    COOP_COMPLEX fk(nz/2+1, ny, nx)
#if HAS_FFTW
    call fft_3d_forward(nx, ny, nz, fx, fk)
#else
    call coop_return_error("coop_fft", "Cannot find FFTW library. Please change the configure.in file.", "stop")
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
    call coop_return_error("coop_fft", "Cannot find FFTW library. Please change the configure.in file.", "stop")
#endif    
  end subroutine coop_fft_backward_3d



  subroutine coop_fft_forward_3d_ss(nz, ny, nx, fx, fk)
    COOP_INT nz, ny, nx
    COOP_REAL fx(nz* ny*nx)
    COOP_COMPLEX fk(nz/2+1, ny, nx)
#if HAS_FFTW
    call fft_3d_forward(nx, ny, nz, fx, fk)
#else
    call coop_return_error("coop_fft", "Cannot find FFTW library. Please change the configure.in file.", "stop")
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
    call coop_return_error("coop_fft", "Cannot find FFTW library. Please change the configure.in file.", "stop")
#endif    
  end subroutine coop_fft_backward_3d_ss


  subroutine coop_generate_white_noise(n, output, norm, mean)
    COOP_REAL,optional::norm, mean
    COOP_INT::n
    COOP_REAL::output(n, n, n)
#if HAS_FFTW    
    COOP_INT::i,j,k      
    COOP_COMPLEX::fc(0:n/2, 0:n-1, 0:n-1)
    COOP_REAL::amp, ksq1d, ksq2d
    if(present(norm))then
       amp = norm/sqrt(coop_4pi)
    else
       amp = 1.d0/sqrt(coop_4pi)
    endif
    do k=1, n-1
       ksq1d =  dble(min(k, n-k))**2
       do j=1, n-1
          ksq2d = dble(min(j, n-j))**2 + ksq1d
          do i=1, (n-1)/2
             fc(i,j,k) = coop_random_complex_Gaussian()*amp/(ksq2d + dble(i)**2)**0.75
          enddo
       enddo
    enddo
    do i=0, n/2, (n-1)/2+1
       ksq1d = dble(i)**2
       do j=1, (n-1)/2
          ksq2d = ksq1d + dble(j)**2
          do k=1, n-1
             fc(i,j,k) =  coop_random_complex_Gaussian()*amp/dble(ksq2d + dble(min(k, n-k))**2)**0.75
             fc(i,n-j,n-k) = conjg(fc(i,j,k))
          enddo
          fc(i,j,0) =  coop_random_complex_Gaussian()*amp/ksq2d**0.75
          fc(i,n-j,0) = conjg(fc(i,j,0))
       enddo
       do j=0, n/2, (n-1)/2+1
          ksq2d = ksq1d + dble(j)**2          
          do k=1, (n-1)/2
             fc(i,j,k) =  coop_random_complex_Gaussian()*amp/(ksq2d+dble(k)**2)**0.75
             fc(i,j,n-k) = conjg(fc(i,j,k))
          enddo
          do k=0, n/2, (n-1)/2+1
             fc(i,j,k) =  coop_random_complex_Gaussian(.true.)/(ksq2d+dble(k)**2+1.d-50)**0.75 !!must be real, add 1.d-50 to avoid divide by zero error, the 0,0,0 component will be fixed later
          enddo
       enddo
    enddo
    if(present(mean))then
       fc(0,0,0) = mean
    else
       fc(0,0,0) = 0.d0
    endif
    call coop_fft_backward(n, n, n, fc, output, .false.)
#else
    call coop_return_error("coop_generate_white_noise", "Cannot find FFTW library. Please change the configure.in file.", "stop")    
#endif    
  end subroutine coop_generate_white_noise
  

end module coop_fft_mod
