program test
  use coop_wrapper_firstorder
  use coop_page_mod
  implicit none
#include "constants.h"
  integer, parameter::nw = 200, nomm = 1000
  real*8,parameter::wmin = -3.d0, wmax = 0.d0, ommmin= 0.01d0, ommmax = 1.d0
  real*8::w = -9.d0
  real*8::omegam = 0.3d0
  real*8,parameter::h = 0.7d0  
  real*8,parameter::omegak = 0.d0
  real*8::t0 
  real*8::eta, q0, dt, muerr
  integer:: i, n, im, iw, nused
  integer::imin(0:nw), imax(0:nw)
  real*8::t, a, z
  real*8,parameter::xlower = -1.3, xupper = 1., ylower=0.8, yupper = 1.2
  logical::firstleft = .true., firstright = .true.

  omegam = 0.3
  w=-1.
  t0 = wcdm_tofa(omegam, omegak, w, 1.d0)
  eta =1.d0 -  1.5 * (omegam/2 + (1.+3.*w)/2. * (1.-omegam-omegak) + 1.) * t0**2
  print*, eta, t0
  print*, page_chiofa(t0, eta, 0.001d0)/(0.674/3000.d0)/1.e3, wcdm_chiofa(omegam, omegak, w, 0.001d0)/(0.674/3000.d0)/1.e3
  stop
  
  imin = -1
  imax = -1
  do iw = 0, nw
     w = wmin + (wmax-wmin)*iw/nw
     do im= 0, nomm
        omegam = ommmin + (ommmax-ommmin)*im/nomm
        t0 = wcdm_tofa(omegam, omegak, w, 1.d0)
        eta =1.d0 -  1.5 * (omegam/2 + (1.+3.*w)/2. * (1.-omegam-omegak) + 1.) * t0**2
        if(eta .ge. xlower .and. eta .le. xupper .and. t0 .ge. ylower .and. t0 .le. yupper)then
           if(imin(iw) .lt. 0) imin(iw) = im
           imax(iw) = im
        endif
     enddo
  enddo

  do iw = 0, nw
     if(imin(iw) .ge. 0)then
        w = wmin + (wmax-wmin)*iw/nw                
        if(firstleft)then
           do im=  imax(iw), imin(iw),-1
              omegam = ommmin + (ommmax-ommmin)*im/nomm
              t0 = wcdm_tofa(omegam, omegak, w, 1.d0)
              eta =1.d0 -  1.5 * (omegam/2 + (1.+3.*w)/2. * (1.-omegam-omegak) + 1.) * t0**2
              write(*,"(2G14.5)") eta, t0              
           enddo
           firstleft = .false.
        else
           omegam = ommmin + (ommmax-ommmin)*imin(iw)/nomm
           t0 = wcdm_tofa(omegam, omegak, w, 1.d0)
           eta =1.d0 -  1.5 * (omegam/2 + (1.+3.*w)/2. * (1.-omegam-omegak) + 1.) * t0**2
           write(*,"(2G14.5)") eta, t0
        endif
     endif
  enddo
  do iw = nw,0,-1
     if(imax(iw) .ge. 0)then
        w = wmin + (wmax-wmin)*iw/nw
        if(firstright)then
           do im=  imin(iw), imax(iw)
              omegam = ommmin + (ommmax-ommmin)*im/nomm
              t0 = wcdm_tofa(omegam, omegak, w, 1.d0)
              eta =1.d0 -  1.5 * (omegam/2 + (1.+3.*w)/2. * (1.-omegam-omegak) + 1.) * t0**2
              write(*,"(2G14.5)") eta, t0              
           enddo
           firstright = .false.
        else
           omegam = ommmin + (ommmax-ommmin)*imax(iw)/nomm
           t0 = wcdm_tofa(omegam, omegak, w, 1.d0)
           eta =1.d0 -  1.5 * (omegam/2 + (1.+3.*w)/2. * (1.-omegam-omegak) + 1.) * t0**2
           write(*,"(2G14.5)") eta, t0
        endif
     endif
  enddo
  


!!$  n  = 200  
!!$     q0 = omegam/2 + (1.+3.*w)/2. * (1.-omegam-omegak) 
!!$     muerr = 0.
!!$     do i=1, n
!!$        z = dble(i)/n * 1.5d0
!!$        !     print*, z,
!!$        if(muerr .lt. abs(page_distance_moduli(t0, eta, omegak, h, z, z ) - wcdm_distance_moduli(omegam, omegak, w, h, z, z))) muerr =  abs(page_distance_moduli(t0, eta, omegak, h, z, z ) - wcdm_distance_moduli(omegam, omegak, w, h, z, z))
!!$     enddo
!!$     print*, omegam, t0, eta, q0, muerr
  
end program test
