program test
  use coop_wrapper_firstorder
  use coop_ellipse_collapse_mod
  implicit none
#include "constants.h"
  !!table format:
  !!  (F_pk, e_nu, p_nu, zvir1)
  COOP_INT, parameter::nthreads = 8
  type(coop_ellipse_collapse_params),dimension(nthreads)::params
  COOP_REAL::Omega_m,w, Omega_k, h
  COOP_REAL,dimension(:,:,:),allocatable::zvir1
  COOP_REAL,dimension(:),allocatable::f, p, e
  COOP_INT::nf, np, ne, if, ie, ip, ithread
  COOP_REAL::fmin, fmax, pmin, pmax, emin, emax
  logical::logf, binary
  type(coop_file)::fp
  COOP_STRING::output
  call coop_get_command_line_argument(key = "out", arg = output)
  call coop_get_command_line_argument(key = "numf", arg = nf, default = 50)
  call coop_get_command_line_argument(key = "nume", arg = ne, default = 50)
  call coop_get_command_line_argument(key = "nump", arg = np, default = 50)
  call coop_get_command_line_argument(key = "omm", arg = Omega_m, default = 0.3d0)
  call coop_get_command_line_argument(key = "omk", arg = Omega_k, default = 0.d0)
  call coop_get_command_line_argument(key = "h", arg = h, default = 0.7d0)
  call coop_get_command_line_argument(key = "w", arg = w, default = -1.d0)
  call coop_get_command_line_argument(key = "fr1", arg = params(1)%collapse_a_ratio(1), default = 0.18d0)
  call coop_get_command_line_argument(key = "fr2", arg = params(1)%collapse_a_ratio(2), default = 0.18d0)
  call coop_get_command_line_argument(key = "fr3", arg = params(1)%collapse_a_ratio(3), default = 0.18d0)
  params(1)%collapse_a_ratio = max(params(1)%collapse_a_ratio, 5.d-4)
  call coop_get_command_line_argument(key = "binary", arg = binary, default = .true.)
  call coop_get_command_line_argument(key = "logf", arg = logf, default = .true.)
  call coop_get_command_line_argument(key = "fmin", arg = fmin, default = 1.5d0)
  call coop_get_command_line_argument(key = "fmax", arg = fmax, default = 15.d0)
  call coop_get_command_line_argument(key = "emin", arg = emin, default = 0.d0)
  call coop_get_command_line_argument(key = "emax", arg = emax, default = 1.d0)
  call coop_get_command_line_argument(key = "pmin", arg = pmin, default = -emax)
  call coop_get_command_line_argument(key = "pmax", arg = pmax, default = emax)
  write(*,*) "========================================================"
  write(*,*) "./TabZ1 -out ... [-binary ...(T)] [-numf ...(100)] [-nume ...(50)] [-nump ...(100)] [-logf ...(T)] [-fmax ...(15)] [-fmin ...(1.5)] [-emax ...(1)] [-emin ...(0)] [-pmax ...(emax)] [-pmin ...(-emax)] ... [-omm ...(0.3)] [-omk ...(0.)] [-h ...(0.7)] [-w ...(-1)] [-fr1 ...(0.18)] [-fr2 ... (0.18)] [-fr3 ... (0.18)]"
  write(*,*) "Examples:"
  write(*,*) "./TabZ1 -out mytable.dat"
  write(*,*) "./TabZ1 -out mytable.dat -numf 200 -nume 200 -nump 200"
  write(*,*) "./TabZ1 -out sphcol.txt -binary F -numf 20 -nume 1 -nump 1 -fmin 1.6 -fmax 10. -emax 0 -omm 0.29 -w -0.9"
  write(*,*) "========================================================"
  if(fmin .gt. fmax .or. emin .gt. emax .or. pmin.gt.pmax .or. fmin .lt. 0.d0 .or. emin .lt. 0.d0 .or. abs(pmin).gt. emax .or. abs(pmax) .gt. emax)then
     write(*,*) "check error: the input must satisfy "
     write(*,*) "numf>=1 and nume>=1 and nump>=1 and fmin<=fmax and emin<=emax and pmin<=pmax and fmin>0 and emin>=0 and abs(pmin)<=emax and abs(pmax)<=emax"
     stop
  endif
  if(nf*ne*np .gt. 10000 .and. .not. binary)then
     write(*,*) "For large numf x nume x nump force binary format."
     binary = .true.
  endif
  if(nf*ne*np .gt. 100000000)then
     write(*,*) "check your numf x nume x nump; table size is too big"
     stop
  endif
  call params(1)%init(Omega_m = Omega_m, Omega_k = Omega_k, h = h, w = w)
  allocate(f(nf), p(np), e(ne), zvir1(np, ne, nf))
  call coop_set_uniform(nf, f, fmin, fmax, logscale = logf)
  call coop_set_uniform(ne, e, emin, emax)
  call coop_set_uniform(np, p, pmin, pmax)
  call coop_prtsystime(.true.)
  if( nf .gt. nthreads)then
     do ithread = 2, nthreads
        params(ithread) = params(1)
     enddo
     !$omp parallel do private(ithread, if, ie, ip)
     do ithread = 1, nthreads
        do if = ithread, nf, nthreads
           do ie = 1, ne
              do ip = 1, np
                 if(abs(p(ip)).gt. e(ie))then
                    zvir1(ip, ie, if) = coop_ellipse_collapse_bad_zvir
                 else
                    call params(ithread)%init(F_pk = f(if), e_nu = e(ie), p_nu=p(ip))
                    zvir1(ip, ie, if) = params(ithread)%zvir1()
                 endif
              enddo
           enddo
        enddo
     enddo
     !$omp end parallel do
  else
     do if = 1, nf
        do ie = 1, ne
           do ip = 1, np
              if(abs(p(ip)).gt. e(ie))then
                 zvir1(ip, ie, if) = coop_ellipse_collapse_bad_zvir
              else
                 call params(1)%init(F_pk = f(if), e_nu = e(ie), p_nu=p(ip))
                 zvir1(ip, ie, if) = params(1)%zvir1()
              endif
           enddo
        enddo
     enddo
  endif
  call coop_prtsystime()
  if(binary)then
     call fp%open(output, "u")
     write(fp%unit) zvir1
     call fp%close()
  else
     call fp%open(output)
     write(fp%unit, "(4A16)") "#  F_pk         ",  "    e_nu  ", "   p_nu   ", "  z_1  "
     do if = 1, nf
        do ie = 1, ne
           do ip = 1, np
              if(zvir1(ip, ie, if) .ge. 0.d0)then
                 write(fp%unit, "(4E16.7)") f(if), e(ie), p(ip), zvir1(ip, ie, if)
              endif
           enddo
        enddo
     enddo
     call fp%close()
  endif
  write(*,*) "The table is dumped to "//trim(output)

end program test
