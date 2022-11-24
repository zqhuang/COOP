!!in linux environment compile the code with
!! gfortran -O3 particle_trajectory.f90 -o calcp
!!and run it with
!! ./calcp 
!! output: t, x, y, z, vx, vy vz (units: microsecond, mm, mm, mm, m/s, m/s, m/s)
module pt_subs
  implicit none
  integer,parameter::IB = selected_int_kind(8)
  integer,parameter::DL = selected_real_kind(12)
  real(kind=DL),parameter::pi  = 3.14159265358979323846264338327950288_dl
  
  !!---------------------USER DEFINED VARIABLES IN SI UNITS-----------------
  real(kind=DL),parameter::mass = 171.0_dl * 1.6605e-27_dl
  real(kind=DL),parameter::charge = 1.6e-19_dl
  !!potential normalization  
  real(kind=DL),parameter::norm = 43503.0_dl 
  !!grid samples in each dimension  
  integer(kind=IB),parameter::nx(3) = (/ 101, 101, 126 /) 
  !!domain boundaries   
  real(kind=DL),parameter::xmin(3) = (/ -1.e-3_dl, -1.e-3_dl, 1.e-3_dl /) 
  real(kind=DL),parameter::xmax(3) = (/ 1.e-3_dl, 1.e-3_dl, 3.5e-3_dl /)
  !!initial position
  real(kind=DL),dimension(3),parameter::xini = (/ 0.0_dl, 0.0_dl, 2.248e-3_dl /)
  !!initial velocity
  real(kind=DL),dimension(3),parameter::vini = (/ 1.0_dl, 1.0_dl, 1.0_dl /)
  !!sine wave frequency (=1/period)
  real(kind=DL),parameter::freq = 2.e7_dl
  !!output file 
  character(LEN=*),parameter::output_file = "ptraj_old.txt"
  !!end time
  real(kind=DL),parameter::tend = 10 * 1.e-6_dl  
  !!--------------------------------------------------------------------  

  real(kind=DL),parameter::charge_by_mass = charge/mass  
  real(kind=DL),parameter::omega = 2.0_dl*pi*freq
  real(kind=DL)::E_field(3, nx(1), nx(2), nx(3))
  real(kind=DL),parameter::dx(3) = (xmax-xmin)/(nx-1)
  
contains

  subroutine load_data(filename)
    character(LEN=*)::filename
    integer(kind=IB)::i,j,k
    character(LEN=1024)::str
    real(kind=DL) line(7)
    open(10, FILE=filename)
    read(10,"(a)") str
    do while(str(1:5) .ne. "x,y,z")
       read(10,"(A)") str
       !       print*, trim(str)
    enddo
    do i=1, nx(1)
       do j=1, nx(2)
          do k=1, nx(3)
             read(10,*, END=100, ERR=100) line
             line(1:3) = line(1:3)*1.e-3_dl
             if(nint((line(1)-xmin(1))/dx(1))+1 .ne. i &
                  .or. nint((line(2)-xmin(2))/dx(2))+1 .ne. j &
                  .or. nint((line(3)-xmin(3))/dx(3))+1 .ne. k )then
                goto 100
             endif
             E_field(:, i, j, k) = line(4:6)
          enddo
       enddo
    enddo
    close(10)
    E_field = E_field * (norm * charge_by_mass*1.e3_dl) 
    return
100 stop "Error in the data file"
  end subroutine load_data


  !!use cubic interpolation in z direction only (motion in x, y directions is negligible anyway)
  subroutine get_acceleration_interp3(x, a)
    real(kind=DL)::x(3), a(3), u(3), d(3)
    integer(kind=IB)::i(3)
    u = (x - xmin)/dx+1.0_dl
    i = floor(u)
    if(any(i.le.1) .or. any(i.ge. nx-1))then
       write(*, *) x
       stop "Error: particle hits the domain boundary."
    endif
    u = u - i
    d = 1.0_dl - u
    !!use E, biliner interpolation    
    a = (1.0_dl + u(3))*(1._dl - u(3)/2._dl)* d(3) *( &
         d(2)*( &
             E_field(:, i(1), i(2), i(3)) * d(1) &
           + E_field(:, i(1)+1, i(2), i(3)) * u(1) &
         ) &
         + u(2)*(  &
         E_field(:, i(1), i(2)+1, i(3)) * d(1) &
         + E_field(:, i(1)+1, i(2)+1, i(3)) * u(1) &
         ) &
         ) &
         + u(3)*(1.0_dl + u(3) * d(3) / 2._dl ) * ( &
         d(2)*( &
             E_field(:, i(1), i(2), i(3)+1) * d(1) &
           + E_field(:, i(1)+1, i(2), i(3)+1) * u(1) &
         ) &
         + u(2)*(  &
         E_field(:, i(1), i(2)+1, i(3)+1) * d(1) &
         + E_field(:, i(1)+1, i(2)+1, i(3)+1) * u(1) &
         ) &
         ) &
         + u(3)*(-1.0_dl/3.0_dl + u(3) *(0.5_dl -  u(3)/ 6._dl)) * ( &
         d(2)*( &
             E_field(:, i(1), i(2), i(3)-1) * d(1) &
           + E_field(:, i(1)+1, i(2), i(3)-1) * u(1) &
         ) &
         + u(2)*(  &
         E_field(:, i(1), i(2)+1, i(3)-1) * d(1) &
         + E_field(:, i(1)+1, i(2)+1, i(3)-1) * u(1) &
         ) &
         ) &
         + u(3)/6._dl * (u(3)**2 - 1.0_dl)  * ( &
         d(2)*( &
             E_field(:, i(1), i(2), i(3)+2) * d(1) &
           + E_field(:, i(1)+1, i(2), i(3)+2) * u(1) &
         ) &
         + u(2)*(  &
         E_field(:, i(1), i(2)+1, i(3)+2) * d(1) &
         + E_field(:, i(1)+1, i(2)+1, i(3)+2) * u(1) &
         ) &
         ) 
  end subroutine get_acceleration_interp3
  
  

  !!y(1:3) particle position
  !!y(4:6) particle velocity
  subroutine eom(n, t, y, yp)
    integer(kind=IB)::n
    real(kind=DL)::t, y(n), yp(n)
    yp(1:3) = y(4:6)
    call get_acceleration_interp3(y(1:3), yp(4:6))
    yp(4:6) = yp(4:6) *sin(omega*t)
  end subroutine eom


  !! given fcn(n, t, y(1:n), yp(1:n)) (input n,t, y, return yp = dy/dt);
  !! evolve t to t+h and update y
  subroutine RungeKutta4th(n, fcn, t, y, h)
    external fcn
    integer(kind=IB)  n
    real(kind=DL) y(n)
    real(kind=DL) t, h
    real(kind=DL) k1(n), k2(n), k3(n), k4(n)
    call fcn(n, t, y, k1)
    call fcn(n, t+h/2._dl, y+k1*(h/2._dl), k2)
    call fcn(n, t+h/2._dl, y+k2*(h/2._dl), k3)
    t = t + h
    call fcn(n, t, y+k3*h, k4)
    y = y+(k1+2._dl*(k2+k3)+k4)*(h/6._dl)
  end subroutine RungeKutta4th

end module pt_subs


program pt
  use pt_subs
  integer(kind=IB),parameter::n=6
  real(kind=DL)::y(n), t, yp(n), Eini
  integer(kind=IB)::i
  real(kind=DL),parameter::stepsize = 2.e-4_dl/freq
  call load_data("19_05_18_01.csv")
  y(1:3) = xini
  y(4:6) = vini
  t = 0.0_dl

  open(11, FILE=output_file)
  write(11, "(7G14.5)") t*1.e6_dl, y(1:3)*1.e3_dl, y(4:6) !!output unit: microsecond, mm, m/s
  do i=1, nint(tend/stepsize)
     call eom(n, t, y, yp)
     call RungeKutta4th(n, eom, t, y, stepsize)
     if(mod(i, 100).eq.0)   write(11, "(7G14.5)") t*1.e6_dl, y(1:3)*1.e3_dl, y(4:6) !!output unit: microsecond, mm, m/s
  enddo
  close(11)
end program pt


