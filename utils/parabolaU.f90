module parabolaU_mod
#include "constants.h"      
  use coop_wrapper_utils

  private

  public::z_equilibrium
  
  !!--------input parameters ---------------------------------
  COOP_REAL::z1 = 0.1d0, z2 = 5.d0, V1=0.d0, V2=-1.d0, V3 = 0.d0
  COOP_REAL::f = 2.1d0, z3 = 29.195d0, gap12=0.1d0, gap23=0.1d0
  !!----------------------------------------------------------
  COOP_REAL::c, rmax, rcut12, rcut23, pivot2, pivot3
  COOP_INT::is1, is2, is3, n_samples
  COOP_INT,parameter:: n = 200, nbuff = 1
  COOP_INT,parameter:: pow1 = 20, pow2 = 20, pow3 = 20
  COOP_INT,parameter:: n_params = pow1+pow2+pow3+3 + 4*nbuff
  COOP_REAL,parameter:: eps = 1.d-3   !!buffer size
  COOP_REAL,dimension(:),allocatable::params, V_samples
  COOP_REAL,dimension(:,:),allocatable::p2V
  COOP_INT::n1, n2, n3, i1min, i1max, i2min, i2max, i3min, i3max
  COOP_REAL::sigma(n), r(n), dr(n), dsdr(n)

contains

  subroutine set_up_params()
    c=1.d0/(4.d0*f)
    rmax=sqrt(z3/c)
    rcut12 = sqrt(z1/c)
    rcut23 = sqrt(z2/c)
    is1 = max(ceiling(n * (rcut12/rmax)*1.2), 5)
    is2=  max(ceiling(n * ((rcut23-rcut12)/rmax)*1.2), 5)
    is3 = max(ceiling(n * ((rmax-rcut23)/rmax)*1.2), 5)
    n_samples = is1 + is2 + is3
    COOP_DEALLOC(params)
    allocate(params(n_params))
    COOP_DEALLOC(V_samples)
    allocate(V_samples(n_samples))
    COOP_DEALLOC(p2V)
    allocate(p2V(n_samples, n_params))
    pivot2 = (rcut12+rcut23)/2.d0
    pivot3 = rmax    
  end subroutine set_up_params

  function dlbydr(r)
    COOP_REAL::r, dlbydr
    dlbydr = sqrt(1.d0 + (4.d0*c**2) * r**2)
  end function dlbydr


  subroutine set_up_rs()
    COOP_REAL::totl, dl, l1, l2, l3, dl1, dl2, dl3, delta_l, r1max, r2min, r2max, r3min
    COOP_INT::i
    l1 = coop_integrate(dlbydr, 0.d0, rcut12, 1.d-7) - gap12/2.d0
    l2 = coop_integrate(dlbydr, rcut12, rcut23, 1.d-7) - (gap12+gap23)/2.d0
    l3 = coop_integrate(dlbydr, rcut23, rmax, 1.d-7) - gap23/2.d0
    
    totl = l1 + l2 + l3
    delta_l = eps*totl !!width of the gap
    dl = totl / (n-4*nbuff)
    n1 = nint(l1/dl) + nbuff
    n3 = nint(l3/dl) + nbuff
    n2 = n - n1 - n3
    dl1 = (l1 - delta_l*nbuff)/(n1-nbuff)
    dl2 = (l2 - delta_l*2*nbuff)/(n2-2*nbuff)
    dl3 = (l3 - delta_l*nbuff)/(n3-nbuff)    

    dr(1) = dl1
    dr(1) = dr(1) + (dl1 - coop_integrate(dlbydr, 0.d0, dr(1), 1.d-6))/dlbydr(dr(1))
    dr(1) = dr(1) + (dl1 - coop_integrate(dlbydr, 0.d0, dr(1), 1.d-6))/dlbydr(dr(1))    
    r(1) = dr(1)/2.d0
    do i=2, n1-nbuff
       dr(i) = dl1/dlbydr(r(i-1)+dr(i-1))
       dr(i) = dr(i) + (dl1 - coop_integrate(dlbydr, r(i-1)+dr(i-1)/2.d0, r(i-1)+dr(i-1)/2.d0+dr(i), 1.d-6))/dlbydr(r(i-1)+dr(i-1)/2.d0+dr(i))
       dr(i) = dr(i) + (dl1 - coop_integrate(dlbydr, r(i-1)+dr(i-1)/2.d0, r(i-1)+dr(i-1)/2.d0+dr(i), 1.d-6))/dlbydr(r(i-1)+dr(i-1)/2.d0+dr(i))
       r(i) = r(i-1) + (dr(i-1)+dr(i))/2.d0 
    enddo
    do i=n1-nbuff+1, n1
       dr(i) = delta_l/dlbydr(r(i-1)+dr(i-1)/2.d0)
       dr(i) = dr(i) + (delta_l - coop_integrate(dlbydr, r(i-1)+dr(i-1)/2.d0, r(i-1)+dr(i-1)/2.d0+dr(i), 1.d-5))/dlbydr(r(i-1)+dr(i-1)/2.d0+dr(i))
       dr(i) = dr(i) + (delta_l - coop_integrate(dlbydr, r(i-1)+dr(i-1)/2.d0, r(i-1)+dr(i-1)/2.d0+dr(i), 1.d-5))/dlbydr(r(i-1)+dr(i-1)/2.d0+dr(i))
       r(i) = r(i-1) + (dr(i-1)+dr(i))/2.d0 
    enddo
    r1max = r(n1)+ dr(n1)/2.d0

    r2min = r1max + gap12/dlbydr(r1max)
    r2min = r1max + gap12*2.d0/(dlbydr(r1max)+dlbydr(r2min))
    r2min = r2min + (gap12 - coop_integrate(dlbydr, r1max, r2min, 1.d-6))/dlbydr(r2min)
    r2min = r2min + (gap12 - coop_integrate(dlbydr, r1max, r2min, 1.d-6))/dlbydr(r2min)

    dr(n1+1) = delta_l/dlbydr(r2min)
    dr(n1+1) = dr(n1+1) + (delta_l - coop_integrate(dlbydr, r2min, r2min + dr(n1+1), 1.d-5))/dlbydr(r2min+dr(n1+1))
    dr(n1+1) = dr(n1+1) + (delta_l - coop_integrate(dlbydr, r2min, r2min + dr(n1+1), 1.d-5))/dlbydr(r2min+dr(n1+1))    
    r(n1+1) = r2min + dr(n1+1)/2.d0
    do i = n1+2, n1+nbuff
       dr(i) = delta_l/dlbydr(r(i-1)+dr(i-1))
       dr(i) = dr(i) + (delta_l - coop_integrate(dlbydr, r(i-1)+dr(i-1)/2.d0, r(i-1)+dr(i-1)/2.d0+dr(i), 1.d-6))/dlbydr(r(i-1)+dr(i-1)/2.d0+dr(i))
       dr(i) = dr(i) + (delta_l - coop_integrate(dlbydr, r(i-1)+dr(i-1)/2.d0, r(i-1)+dr(i-1)/2.d0+dr(i), 1.d-6))/dlbydr(r(i-1)+dr(i-1)/2.d0+dr(i))       
       r(i) = r(i-1) + (dr(i-1)+dr(i))/2.d0 
    enddo
    do i = n1+nbuff+1, n1+n2-nbuff
       dr(i) = dl2/dlbydr(r(i-1)+dr(i-1))
       dr(i) = dr(i) + (dl2 - coop_integrate(dlbydr, r(i-1)+dr(i-1)/2.d0, r(i-1)+dr(i-1)/2.d0+dr(i), 1.d-6))/dlbydr(r(i-1)+dr(i-1)/2.d0+dr(i))
       dr(i) = dr(i) + (dl2 - coop_integrate(dlbydr, r(i-1)+dr(i-1)/2.d0, r(i-1)+dr(i-1)/2.d0+dr(i), 1.d-6))/dlbydr(r(i-1)+dr(i-1)/2.d0+dr(i))       
       r(i) = r(i-1) + (dr(i-1)+dr(i))/2.d0 
    enddo
    do i=n1+n2-nbuff+1, n1+n2
       dr(i) = delta_l/dlbydr(r(i-1)+dr(i-1)/2.d0)
       dr(i) = dr(i) + (delta_l - coop_integrate(dlbydr, r(i-1)+dr(i-1)/2.d0, r(i-1)+dr(i-1)/2.d0+dr(i), 1.d-5))/dlbydr(r(i-1)+dr(i-1)/2.d0+dr(i))
       dr(i) = dr(i) + (delta_l - coop_integrate(dlbydr, r(i-1)+dr(i-1)/2.d0, r(i-1)+dr(i-1)/2.d0+dr(i), 1.d-5))/dlbydr(r(i-1)+dr(i-1)/2.d0+dr(i))
       r(i) = r(i-1) + (dr(i-1)+dr(i))/2.d0 
    enddo
    r2max = r(n1+n2)+dr(n1+n2)/2.d0

    r3min = r2max + gap23/dlbydr(r2max)
    r3min = r2max + gap23*2.d0/(dlbydr(r2max) + dlbydr(r3min))
    r3min = r3min + (gap23 - coop_integrate(dlbydr, r2max, r3min, 1.d-6))/dlbydr(r3min)
    r3min = r3min + (gap23 - coop_integrate(dlbydr, r2max, r3min, 1.d-6))/dlbydr(r3min)

    dr(n1+n2+1) = delta_l/dlbydr(r3min)
    dr(n1+n2+1) = dr(n1+n2+1) + (delta_l - coop_integrate(dlbydr, r3min, r3min+dr(n1+n2+1), 1.d-5))/dlbydr(r3min + dr(n1+n2+1))
    dr(n1+n2+1) = dr(n1+n2+1) + (delta_l - coop_integrate(dlbydr, r3min, r3min+dr(n1+n2+1), 1.d-5))/dlbydr(r3min + dr(n1+n2+1))
    r(n1+n2+1) = r3min + dr(n1+n2+1)/2.d0
    do i=n1+n2+2, n1+n2+nbuff
       dr(i) = delta_l/dlbydr(r(i-1)+dr(i-1))
       dr(i) = dr(i) + (delta_l - coop_integrate(dlbydr, r(i-1)+dr(i-1)/2.d0, r(i-1)+dr(i-1)/2.d0+dr(i), 1.d-6))/dlbydr(r(i-1)+dr(i-1)/2.d0+dr(i))
       dr(i) = dr(i) + (delta_l - coop_integrate(dlbydr, r(i-1)+dr(i-1)/2.d0, r(i-1)+dr(i-1)/2.d0+dr(i), 1.d-6))/dlbydr(r(i-1)+dr(i-1)/2.d0+dr(i))       
       r(i) = r(i-1) + (dr(i-1)+dr(i))/2.d0 
    enddo
    do i = n1+n2+nbuff+1, n
       dr(i) = dl3/dlbydr(r(i-1)+dr(i-1))
       dr(i) = dr(i) + (dl3 - coop_integrate(dlbydr, r(i-1)+dr(i-1)/2.d0, r(i-1)+dr(i-1)/2.d0+dr(i), 1.d-6))/dlbydr(r(i-1)+dr(i-1)/2.d0+dr(i))
       dr(i) = dr(i) + (dl3 - coop_integrate(dlbydr, r(i-1)+dr(i-1)/2.d0, r(i-1)+dr(i-1)/2.d0+dr(i), 1.d-6))/dlbydr(r(i-1)+dr(i-1)/2.d0+dr(i))       
       r(i) = r(i-1) + (dr(i-1)+dr(i))/2.d0 
    enddo
    i1min = 1
    i1max = n1
    i2min = n1+1
    i2max = n1+n2
    i3min = n1+n2+1
    i3max = n
  end subroutine set_up_rs
  

  function pot_per_theta(theta, args)
    type(coop_arguments)::args
    COOP_REAL::theta, pot_per_theta
    pot_per_theta = 1.d0/sqrt( args%r(1) - args%r(2) * cos(theta) )
  end function pot_per_theta


  function pot_per_dr(r, a, sigma, theta_min, theta_max)
    COOP_REAL::pot_per_dr, sigma, r, a, tmin, tmax
    COOP_REAL,optional::theta_min, theta_max
    type(coop_arguments)::args
    if(present(theta_min))then
       tmin = theta_min
    else
       tmin  = 0.d0
    endif
    if(present(theta_max))then
       tmax = theta_max
    else
       tmax = coop_2pi
    endif
    if(r .eq. 0.d0 .and. a .eq. 0.d0)then
       pot_per_dr = (tmax-tmin) * sigma
    else
       call args%init(r = (/ r ** 2 + a **2 + (c*(r-a)*(r+a))**2, &
             2 * r * a /) )
       pot_per_dr = sigma * r * dlbydr(r) * coop_integrate(pot_per_theta, tmin, tmax, args, 1.d-6)
       call args%free()
    endif
  end function pot_per_dr

  
  function pot_at_a(a, n, r, dr, sigma, dsdr)
    COOP_REAL::a, pot_at_a
    COOP_INT::n, i
    COOP_REAL::r(n), dr(n)
    COOP_REAL::sigma(n), dsdr(n)
    pot_at_a = 0.d0
    !$omp parallel do reduction(+:pot_at_a) 
    do i = 1, n
       if( r(i) - dr(i)/2.d0 .gt. a .or. r(i)+dr(i)/2.d0 .lt. a)then
          pot_at_a =  pot_at_a + ( &
               pot_per_dr(r(i), a, sigma(i)) &
               + pot_per_dr(r(i) + (2.d0/3.d0)*dr(i), a, sigma(i) + (2.d0/3.d0)*dr(i)*dsdr(i))  &
               + pot_per_dr(r(i) - (2.d0/3.d0)*dr(i), a, sigma(i) - (2.d0/3.d0)*dr(i)*dsdr(i))  &               
               )* dr(i)/3.d0
       else
          pot_at_a = pot_at_a +  ( &
               pot_per_dr(r(i), a, sigma(i), dr(i)/r(i)/2.d0, coop_2pi - dr(i)/r(i)/2.d0) &
               + pot_per_dr(r(i)+dr(i)*(2.d0/3.d0), a, sigma(i)+dsdr(i)*dr(i)*(2.d0/3.d0), dr(i)/r(i)/2.d0, coop_2pi - dr(i)/r(i)/2.d0) &
               + pot_per_dr(r(i)-dr(i)*(2.d0/3.d0), a, sigma(i)-dsdr(i)*dr(i)*(2.d0/3.d0), dr(i)/r(i)/2.d0, coop_2pi - dr(i)/r(i)/2.d0) &               
               ) * dr(i)/3.d0 + coop_2pi*sigma(i)*dr(i)
       endif
    enddo
    !$omp end parallel do
  end function pot_at_a


  

  subroutine get_all_pot(n, r, dr, sigma, dsdr) 
    COOP_INT::n, m, msum
    COOP_REAL::r(n), dr(n), sigma(n), dsdr(n)
    COOP_REAL::chi2, a
    COOP_INT::i
    logical,save::init
    chi2 = 0.d0
    m = is1
    msum = m
    do i=1, m
       a = rcut12*(i-1)/m
       V_samples(i) = pot_at_a(a, n, r, dr, sigma, dsdr)
    enddo
    m = is2
    do i=1, m
       a = rcut12 + (rcut23 - rcut12) * i/(m+1)
       V_samples(msum+i) = pot_at_a(a, n, r, dr, sigma, dsdr)
    enddo
    msum = msum + m       
    m = is3
    do i=1, m
       a = rcut23 + (rmax - rcut23) * i / m
       V_samples(msum+i) =  pot_at_a(a, n, r, dr, sigma, dsdr)
    enddo
  end subroutine get_all_pot


  function axis_potential(z) result(V)
    COOP_REAL::z, V
    V = coop_2pi*sum(sigma*sqrt((1.d0+4*c**2*r**2)/(r**2+(c*r**2-z)**2))*r*dr)
  end function axis_potential

  function axis_potential_derv(z) result(dVdz)
    COOP_REAL::z, dVdz
    dVdz = coop_2pi*sum(sigma*sqrt((1.d0+4*c**2*r**2)/((r**2+(c*r**2-z)**2))**3)*(c*r**2-z)*r*dr)
  end function axis_potential_derv
  

  function full_potential(z, a) result(V)
    COOP_INT,parameter::nsteps = 500
    COOP_REAL::z, V, theta, dtheta, a
    COOP_INT::i
    dtheta = coop_2pi/nsteps
    V=0.d0
    do i=1, nsteps
       theta = (i-1)*dtheta
       V = V + sum(sigma*sqrt((1.d0+4*c**2*r**2)/(r**2+a**2-2*a*cos(theta)*r+(c*r**2-z)**2))*r*dr)
    enddo
    V = V*dtheta
  end function full_potential
  

  subroutine get_sigma()
    COOP_INT:: j

    COOP_REAL::a1(0:pow1), a2(0:pow2), a3(0:pow3)    
    a1 = params(1:pow1+1)
    a2 = params(pow1+2:pow1+pow2+2)
    a3 = params(pow1+pow2+3:pow1+pow2+pow3+3)
    sigma(i1min:i1max) = a1(0)
    dsdr = 0.d0
    do j = 1, pow1
       if(a1(j).ne. 0.d0)then
          sigma(i1min:i1max) =  sigma(i1min:i1max) + a1(j)*(r(i1min:i1max)/rcut12)**(2*j)
          dsdr(i1min:i1max) = dsdr(i1min:i1max) + a1(j)/rcut12*(2*j)*(r(i1min:i1max)/rcut12)**(2*j-1)
       endif
    enddo
    sigma(i2min:i2max) = a2(0)
    do j = 1, pow2
       if(a2(j).ne. 0.d0)then       
          sigma(i2min:i2max) =  sigma(i2min:i2max) + a2(j)*((r(i2min:i2max)-pivot2)/(rcut23-rcut12))**j
          dsdr(i2min:i2max)  = dsdr(i2min:i2max) + a2(j)/(rcut23-rcut12)*j*((r(i2min:i2max)-pivot2)/(rcut23-rcut12))**(j-1)
       endif
    enddo
    sigma(i3min:i3max) = a3(0)
    do j = 1, pow3
       if(a3(j).ne.0.d0)then
          sigma(i3min:i3max) =  sigma(i3min:i3max) + a3(j)*((r(i3min:i3max)-pivot3)/(rmax-rcut23))**j
          dsdr(i3min:i3max) = dsdr(i3min:i3max)  + a3(j)/(rmax-rcut23)*j*((r(i3min:i3max)-pivot3)/(rmax-rcut23))**(j-1)
       endif
    enddo
    do j=1, nbuff
       sigma(i1max+1-j) = sigma(i1max+1-j)+params(n_params-4*j+1)
       sigma(i2min+j-1) = sigma(i2min+j-1)+params(n_params-4*j+2)
       sigma(i2max+1-j) = sigma(i2max+1-j)+params(n_params-4*j+3)
       sigma(i3min+j-1) = sigma(i3min+j-1)+params(n_params-4*j+4)
    enddo
  end subroutine get_sigma

  function z_equilibrium(focal_length, z_gap12, z_gap23, z_max, size_gap12, size_gap23, potential_1, potential_2, potential_3) result(zeq)
    COOP_REAL::focal_length, z_gap12, z_gap23, z_max, size_gap12, size_gap23, potential_1, potential_2, potential_3
    COOP_REAL::zmin, zmax, zeq, Vpmin, Vpmax, Vpmid
    COOP_INT::i
    f = focal_length
    z1 = z_gap12
    z2 = z_gap23
    z3 = z_max
    gap12 = size_gap12
    gap23 = size_gap23
    V1 = potential_1
    V2 = potential_2
    V3 = potential_3
    call set_up_params()
    call set_up_rs()
    do i=1, n_params
       params = 0.d0
       params(i) = 1.d0
       call get_sigma()
       call get_all_pot(n, r, dr, sigma, dsdr) 
       p2V(:, i) = V_samples
    enddo
    V_samples(1:is1) = V1
    V_samples(is1+1: is1+is2) = V2
    V_samples(is1+is2+1:is1+is2+is3) = V3
    call coop_svd_least_square_one(n_samples, n_params, p2V, V_samples, params)
    call get_sigma()
    zmin = 0.5d0*z1
    zmax = z2*0.2+z3*0.8
    Vpmin = axis_potential_derv(zmin)
    Vpmax = axis_potential_derv(zmax)
    if(Vpmax * Vpmin .ge. 0.d0)then
       zmax = z2*0.1+z3*0.9
       Vpmax = axis_potential_derv(zmax)
    endif
    if(Vpmin * Vpmax .ge. 0.d0)then
       write(*,*) "Weird configuration; Failed to search z_eq; You may want to have a look at the V(z) function."
       do i=1, 101
          zmin = z3*(i-1)/100.d0
          write(*,*) zmin, axis_potential(zmin), axis_potential_derv(zmin)
       enddo
       zeq = -1.d0
       return
    endif
    do while(zmax - zmin .gt. 1.d-5)
       zeq = (zmin+zmax)/2.d0
       Vpmid = axis_potential_derv(zeq)
       if(Vpmid * Vpmin .lt. 0.d0)then
          zmax = zeq
          Vpmax = Vpmid
       else
          zmin = zeq
          Vpmin = Vpmid
       endif
    enddo
  end function z_equilibrium

end module parabolaU_mod
