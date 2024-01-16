!!QuN: QUadruple-precision delta N calculator (if your machine support quadruple-precision)
!!Please compile with gfortran
!!Solve the background evoluation and obtain delta N as a function of initial conditions
!!Use delta N formula to compute power spectrum and NonGaussianity of the primordial curvature perturbations.
!!by Zhiqi Huang
!!May 2019

module QuN_subs
  implicit none
  integer,parameter::IB = selected_int_kind(8)
  integer,parameter::DL = selected_real_kind(20)
  integer(kind=IB),parameter::verbose = 1
  real(kind=DL),parameter::Mpl = 1024.0_dl  !!this can be an arbitrary constant  
  real(kind=DL),parameter::pi  = 3.14159265358979323846264338327950288_dl
  
  !!======== User defined variable ===================
  integer(kind=IB),parameter::User_nflds = 2 !!number of fields
  !!parameters used in the potential
  real(kind=DL),parameter::lambda = 1.e-13_dl
  real(kind=DL),parameter::gsq = lambda * 200.0_dl 
  !!--------------------------------------------------
contains

  !!============== user defined functions =============
  !!calculate potential and its derivatives w.r.t. the fields
  !!I do not need second derivatives because this code does not do perturbations
  subroutine User_calc_potential(f, V, Vp)
    real(kind=DL),INTENT(IN)::f(User_nflds)
    !!by default f(1) is the inflaton;
    !!define your macros below
    real(kind=DL),optional::V, Vp(User_nflds)
    if(present(V))then
       V = (lambda/4.0_dl) * f(1) ** 4 + (gsq/2.0_dl) * f(1)**2 * f(2) ** 2
    endif
    if(present(Vp))then
       Vp(1) = lambda * f(1) ** 3 + gsq * f(1) * f(2) ** 2
       Vp(2) = gsq * f(1)**2 * f(2) 
    endif
  end subroutine User_calc_potential
  
  !!------------------end of User section --------------

  function int2str(i)
    integer(kind=IB)::i
    character(LEN=64)::int2str
    write(int2str,*) i
    int2str= trim(adjustl(int2str))
  End function int2str
  

  subroutine QuN_canonical_inflation(y, yp)
    integer(kind=IB),parameter::n = 2*User_nflds+2
    real(kind=DL)::y(n), yp(n)
    !! t cosmological time
    !!y(1:User_nflds) the field values
    !!y(User_nflds+1:2*User_nflds) dot (field values)
    !!y(2*User_nflds+1)  ln a
    !!y(2*User_nflds+2)  H
    !!yp = dy/dt
    integer(kind=IB)::i
    yp(1:User_nflds) = y(User_nflds+1:2*User_nflds)  !! d f/ dt  = \dot f 
    call User_calc_potential(f = y(1:User_nflds), Vp = yp(User_nflds+1:2*User_nflds)) !!get V and Vp (saved in yp(User_nflds+1:2*User_nflds))
    yp(User_nflds+1:2*User_nflds) = - 3.0_dl * y(2*User_nflds+2) *  yp(1:User_nflds) - yp(User_nflds+1:2*User_nflds) !!\ddot f + 3H \dot f + dV/df = 0
    yp(2*User_nflds+1) = y(2*User_nflds+2)  !! d(ln a)/dt = H
    yp(2*User_nflds+2) = -sum(yp(1:User_nflds)**2)/(2.0_dl * Mpl**2) !!\dot H = - (\rho + p)/(2M_p^2) = -(\dot\phi^2)/(2M_p^2)
  end subroutine QuN_canonical_inflation

  subroutine QuN_canonical_RK4th(y, h)
    integer(kind=IB),parameter::n = 2*User_nflds+2    
    real(kind=DL)::y(n)   !!variables
    real(kind=DL) h(3) !!time, h(1) = step size, h(2) = h(1)/2, h(3) = h(1)/6
    real(kind=DL),save::workspace(n, 4)
    call QuN_canonical_inflation(y, workspace(1:n, 1))
    call QuN_canonical_inflation(y + workspace(1:n,1) * h(2), workspace(1:n,2))
    call QuN_canonical_inflation(y + workspace(1:n,2) * h(2), workspace(1:n,3))
    call QuN_canonical_inflation(y + workspace(1:n,3)*h(1), workspace(1:n,4))
    y = y + (workspace(1:n, 1) + 2._dl * (workspace(1:n,2) + workspace(1:n,3)) + workspace(1:n,4)) * h(3)
  end subroutine QuN_Canonical_RK4th

  function QuN_canonical_Hubble(f, dot_f) result(H)
    real(kind=DL),dimension(User_nflds)::f
    real(kind=DL),dimension(User_nflds),optional::dot_f
    real(kind=DL)::H, lastH
    real(kind=DL)::V, Vp(User_nflds), a, b, c
    integer(kind=IB)::i
    if(present(dot_f))then
       call User_calc_potential(f = f, V=V)           
       H =  sqrt((sum(dot_f**2)/2.0_dl + V)/(3.0_dl*Mpl**2))
    else
       call User_calc_potential(f = f, V=V, Vp=Vp)
       a = 54.0_dl*Mpl**2
       b = -18.0_dl*V
       c = sum(Vp**2)
       H = sqrt((-b + sqrt(b**2+4.0_dl*a*c))/(2.0_dl*a))
    endif
  end function QuN_Canonical_Hubble


  subroutine QuN_calc_canonical_efolds(Hend, efolds, f_ini, dot_f_ini)
    real(kind=DL),parameter::accuracy = 4.e-4_dl    
    integer(kind=IB),parameter::n = 2*User_nflds+2        
    real(kind=DL),dimension(User_nflds)::f_ini
    real(kind=DL),dimension(User_nflds),optional::dot_f_ini    
    real(kind=DL),dimension(n)::y, yp
    real(kind=DL)::Hend, efolds, V, Vp(User_nflds), h(3)
    integer(kind=IB)::count, update_freq
    y(1:User_nflds) = f_ini
    y(n-1) = 0.0_dl        
    if(present(dot_f_ini))then
       y(User_nflds+1:2*User_nflds) = dot_f_ini
       y(n) = QuN_canonical_Hubble(f_ini, dot_f_ini)
    else
       y(n) = QuN_canonical_Hubble(f_ini)
       call User_calc_potential(f = f_ini, V=V, Vp=Vp)
       y(User_nflds+1:2*User_nflds) = - Vp/(3.0_dl*y(n))       
    endif
    call Qun_canonical_inflation(y, yp)
    if(verbose>1)then
       write(*,*) "Initial field values: ", y(1:User_nflds)/Mpl
       write(*,*) "Initial field derivatives: ", y(User_Nflds+1:2*User_nflds)/Mpl**2
       write(*,*) "Initial Hubble: ", y(n)/Mpl
    endif
    if(y(n) .lt. Hend .or. yp(n) .gt. 0.0_dl)then
       if(verbose>1)write(*,*) "Error: H_ini < H_end or dH/dt>0"
       efolds = -1.0_dl
       return
    endif
    if(yp(n) .lt. 0.0_dl)then
       call update_step( max(min( accuracy / y(n) , (Hend-y(n))/yp(n)/(update_freq*1.1_dl) ), 1.e-13_dl/y(n)) )
    else
       call update_step( accuracy / y(n) )
    endif
    update_freq = 20
    count = 0
    do while(y(n) .gt. Hend)
       call QuN_canonical_RK4th(y, h)
       call check_step()
    enddo
    !! roll back the overshoot
    call Qun_canonical_inflation(y, yp)        
    if(yp(n) .lt. 0.0_dl )then
       call update_step( (Hend-y(n))/yp(n) )
       if( h(1) * y(n) .lt. accuracy)then
          call QuN_canonical_RK4th(y, h)
       else
          if(verbose>1)write(*,*) "Warning: unstable solution around H_end!"
       endif
    endif
    if(verbose>0)call check_energy_conservation()
    efolds = y(n-1)
  contains
    subroutine update_step(hwant)
      real(kind=DL)::hwant
      h(1) = hwant
      h(2) = hwant/2.0_dl
      h(3) = hwant/6.0_dl
    end subroutine update_step
    
    subroutine check_energy_conservation()
      call User_calc_potential(f = y(1:User_nflds), V=V)
      write(*,"(A, E12.3)") "Energy conservation accuracy:", &
           (V+sum(y(User_nflds+1:User_nflds*2)**2)/2.0_dl)/(3.0_dl*Mpl**2)/y(n)**2-1.0_dl
    end subroutine check_energy_conservation

    subroutine check_step()
      real(kind=DL)::eps
      count = count + 1
      if(mod(count, update_freq) .eq. 0)then
         call Qun_canonical_inflation(y, yp)
         eps = sum(yp(1:User_nflds)**2) / (2.0_dl*Mpl**2) / y(n)**2
         if(yp(n) .lt. 0.0_dl)then
            call update_step( max(min( accuracy / y(n) / (1.0_dl + 50.0_dl*eps) , &
                 (Hend-y(n))/yp(n)/(update_freq*1.1_dl) ), 1.e-13_dl/y(n)) )
         else
            call update_step( accuracy / y(n) / (1.0_dl + 50.0_dl*eps) )
         endif
      endif

    end subroutine check_step
    
  end subroutine QuN_calc_canonical_efolds
  
end module QuN_subs


program QuN_main
  use QuN_subs
  implicit none
  !!central initial values  
  real(kind=DL),dimension(User_nflds),parameter::f_ini = (/ 25.0_dl * Mpl, 1.e-5_dl*Mpl /)
  !!====scan a 2D lattice of (f(1), f(2))============
  !!lattice boundary in unit of H/(2pi)  
  real(kind=DL),parameter::scan_boundary = 2._dl
  integer(kind=IB),parameter::n = 5 !!lattice size (-n:n, -n:n) totally (2n+1)^2 points
  real(kind=DL),dimension(-n:n, -n:n)::efolds  !!save the result
  !!----------------------------------------------------
  real(kind=DL)::H_end, V, H_ini, Vp(1:User_nflds), step
  integer(kind=IB)::i, j
  if(verbose>1) write(*,"(A, I3, A)") "Working with ", DL*8, "bits precision."
  H_ini = QuN_Canonical_Hubble( f_ini )
  call User_calc_potential( f = f_ini, V = V, Vp = Vp)
  H_end = QuN_Canonical_Hubble( f_ini/10.0_dl )
  step = H_ini/(2.0_dl*pi)*scan_boundary/n
  !$omp parallel do private(i, j)
  do i=-n, n
     do j=-n, n
        call QuN_calc_canonical_efolds(H_end, efolds(i, j), f_ini + (/ i, j /)*step)
     enddo
  enddo
  !$omp end parallel do
  efolds = efolds - efolds(0, 0)
  open(10, file="deltaN.txt")
  do i=-n, n
     write(10, "("//trim(int2str(2*n+1))//"E32.20)") efolds(i, :)
  enddo
  close(10)
end program QuN_main
