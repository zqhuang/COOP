module coop_coupledDE_collapse_mod
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"

#if DO_COUPLED_DE
  !!by default everything is a function of scale factor a;  a = 1 today.
  !!physical time t is used in the definition of some variables, its unit = 1/H_0 

  !!this is a global accuracy parameter
  !!for high accuracy test you can use something like 1.e-4
  !!for normal runs you can use something ~ 1.e-3
  COOP_REAL, parameter::coop_coupledDE_collapse_accuracy = 1.d-3

  !!set z_vir = this value if not collapsed
  COOP_REAL, parameter::coop_coupledDE_collapse_bad_zvir = -1.d0

  type coop_coupledDE_collapse_params
     COOP_REAL,dimension(3)::lambda = (/ 0.d0, 0.d0, 0.d0 /) !!lambda's
     logical::is_spherical = .false.
     COOP_INT::num_ode_vars = 8
     COOP_REAL,dimension(3)::collapse_a_ratio =  (/ 0.178, 0.178, 0.178 /) ! (18pi^2)^{-1/3}
     type(coop_2dfunction)::bprime  
     type(coop_cosmology_firstorder)::cosmology
     logical::normalize_to_early = .false.  !!if true, set the initial delta = F_pk * a
   contains
     procedure::free => coop_coupledDE_collapse_params_free  !!subroutine; no argument; release the memory allocated for this object
     procedure::init => coop_coupledDE_collapse_params_init  !!subroutine; this%init(Omega_m, w, Omega_k, h, F_pk, e_nu, p_nu);  initialize the object with these parameters
     procedure::get_bprime => coop_coupledDE_collapse_params_get_bprime !!subroutine; this%get_bprime(x, bprime)  where x(1:3) is the input array x_1, x_2, x_3, bprime(1:3) is the output b'_1, b'_2, b'_3; 
     procedure::Growth_D => coop_coupledDE_collapse_params_Growth_D !! this%Growth_D(a) is a function that returns the growth factor D; input a is the scale factor
     procedure::Growth_H_D => coop_coupledDE_collapse_params_Growth_H_D  !!this%Growth_H_D(a) is a function that returns d ln D/d t, where a is the scale factor
     procedure::dadt => coop_coupledDE_collapse_params_aH  !!this%dadt(a) is a function that  returns da/dt
     procedure::ddotabya => coop_coupledDE_collapse_params_ddotabya !!this%ddotabya is a function that returns ( \ddot a / a )
     procedure::ddotphi => coop_coupledDE_collapse_params_ddotphi !!this%ddotphi(a, \phi, \dot\phi, \rho_c, intQofphi) returns \ddot \phi
     procedure::set_initial_conditions => coop_coupledDE_collapse_params_set_initial_conditions  !!this%set_initial_conditions(y) is a subroutine set the initial conditions for y = (x_1, x_2, x_3, d x_1/dt, d x_2/dt, d x_3/dt) 
     procedure::evolve=> coop_coupledDE_collapse_params_evolve  !!this%evolve(a, y, a_end) is a subroutine that evolves the current  y = (x_1, x_2, x_3, d x_1/dt, d x_2/dt, d x_3/dt) from current scale factor a to the scale factor a_end; the input is y at a; after the call a becomes a_end.
     procedure::get_solution => coop_coupledDE_collapse_params_get_solution  !!this%get_solution(a_arr, x_arr) is a subroutine; input a_arr(1:n) is an array of scale factors in ascending order; return the result in x_arr(1:3, 1:n), where x_arr(1:3, i) is the solution of (x_1, x_2, x_3) at scale factor a_arr(i).
     procedure::zvir1 => coop_coupledDE_collapse_params_zvir1  !!this%zvir1() is a function that returns zvir1

  end type coop_coupledDE_collapse_params


contains

  subroutine coop_coupledDE_collapse_params_set_initial_conditions(this, a_ini, y)
    !!set inital vector y = ( x_1, x_2,  x_3,  d x_1/dt, d x_2/dt, d x_3/dt, phi, d phi/dt ) at a = a_ini
    class(coop_coupledDE_collapse_params)::this
    COOP_REAL,intent(IN)::a_ini
    COOP_REAL::y(:), D_ini,  dadt_ini, corr(3), suml, suml2, H_D
    if(this%normalize_to_early)then
       D_ini = coop_growth_fitting(0.3d0, -1.d0, 1.d0/a_ini - 1.d0) 
       H_D = sqrt(0.3/a_ini**3 + 0.7)
    else
       D_ini = this%Growth_D(a_ini)
       H_D = this%Growth_H_D(a_ini)
    endif
    dadt_ini = this%dadt(a_ini)
    suml = sum(this%lambda)
    suml2 = sum(this%lambda**2)
    corr = (1.2d0*this%lambda*suml + 0.6d0*suml2 + 11.d0/35.d0*(suml**2-suml2)-3.d0*this%lambda**2)/10.d0  !!this is for 2nd-order correction of initial conditions; in spherical case corr = 3/7 lambda^2
    y(1:3) = a_ini *(1.d0-this%lambda*D_ini - corr*D_ini**2)
    y(4:6) = y(1:3)*dadt_ini/a_ini - a_ini*(D_ini*H_D*(this%lambda + 2.d0*corr*D_ini))  
    y(7) = O0_DE(this%cosmology)%cplDE_phi_lna%eval(log(a_ini))  !phi
    y(8) = O0_DE(this%cosmology)%cplDE_phi_prime_lna%eval(log(a_ini))*this%cosmology%Hratio(a_ini)   !!d phi/dt
  end subroutine coop_coupledDE_collapse_params_set_initial_conditions

  subroutine coop_coupledDE_collapse_odes(n, a, y, dyda, params)
  !!the ODE that evolves  y = ( x_1, x_2,  x_3,  d x_1/dt, d x_2/dt, d x_3/dt, phi, d phi/dt )
  !!other inputs: 
  !!n = 8 is the dimension of y
  !!a is the scale factor
  !!params is the object containing all the parameters and methods
  !! return dyda = d y/d a
    COOP_INT::n
    COOP_REAL::a, y(n), dyda(n), dadt, bprime(3), growthD, delta, dark_Energy_term, radiation_term,  rho_m, rhombar, rhombarby3, dlnmdt, q, dvdphi
    type(coop_coupledDE_collapse_params)::params
    COOP_REAL,parameter::eps = coop_coupledDE_collapse_accuracy, max_decay_rate = 1.d5
    COOP_REAL::suppression_factor, arat(3)
    dadt = params%dadt(a)
    radiation_term = 0.d0
    dark_Energy_term =  -(y(8)**2-O0_DE(params%cosmology)%cplde_Vofphi%eval(y(7)))*2.d0/3.d0
    if(.not. params%cosmology%baryon_is_coupled)stop "The current version cannot solve halo collapse for models where baryon is not coupled."
    rho_m = 3.d0*(params%cosmology%Omega_c_bare+params%cosmology%Omega_b_bare)*exp(O0_DE(params%cosmology)%cplde_intQofphi%eval(y(7)))/(y(1)*y(2)*y(3))
    dlnmdt = O0_DE(params%cosmology)%cplde_intQofphi%derivative(y(7))*y(8)
    if(params%is_spherical)then
       arat(1) = (y(1)/a/params%collapse_a_ratio(1) - 1.d0)/eps
       if(arat(1) .lt. -1.d0 .and. y(4) .lt. 0.d0)then  !!collapsed; freeze it
          dyda(1) = 0.d0
          dyda(4) = 0.d0
       else  !!still collapsing
          !!----------- equation for x_1 -------------------------
          dyda(1) = y(4) / dadt    !!d x_1/da = d x_1/dt /(da/dt)
          dyda(4) = y(1)/dadt/2.d0 * ( &   !! d( dx_1/dt)/da =(d^2 x_1/dt^2)/(da/dt)
               dark_Energy_term  &
               + radiation_term &
               -  rho_m/3.d0 &  !!matter contribution
               ) - dlnmdt*y(4)/dadt
          !!-----------end of equation for x_1 -------------------------
          if(arat(1) .lt. 0.d0  .and. y(4) .lt. 0.d0)then  !!do suppression around x_i/a = fr_i so that the derivative is continuous; no need to change
             suppression_factor =  sin(coop_pio2*(1.d0+arat(1)))**4
             dyda(1) = dyda(1)*suppression_factor
             dyda(4) = dyda(4)*suppression_factor
             
          endif
       endif
       dyda(2:3) = dyda(1)
       dyda(5:6) = dyda(4)
    else
       rhombar = O0_CDM(params%cosmology)%density(a) + O0_BARYON(params%cosmology)%density(a)
       rhombarby3 = rhombar/3.d0
       delta = rho_m/rhombar-1.d0  !a**3/(y(1)*y(2)*y(3))-1.d0
       call params%get_bprime(y(1:3), bprime)       
       growthD = params%growth_D(a) !!this isn't accurate because in general D is scale dependent in the coupled DE model
       arat = (y(1:3)/a/params%collapse_a_ratio - 1.d0)/eps
       if(arat(1) .lt. -1.d0  .and. y(4) .lt. 0.d0)then  !!collapsed; freeze it
          dyda(1) = 0.d0
          dyda(4) = 0.d0
       else  !!still collapsing
          !!----------- equation for x_1 -------------------------
          dyda(1) = y(4) / dadt    !!d x_1/da = d x_1/dt /(da/dt)
          dyda(4) = y(1)/dadt/2.d0 * ( &   !! d( dx_1/dt)/da =(d^2 x_1/dt^2)/(da/dt)
               dark_Energy_term  &
               + radiation_term &
               -  rhombarby3*(1.d0 + delta *(1.d0+bprime(1)*1.5d0) + (3.d0*params%lambda(1)-sum(params%lambda))*growthD ) &  !!matter contribution
               )- dlnmdt*y(4)/dadt
          !!----------- end of equation for x_1 -------------------------
          if(arat(1) .lt. 0.d0)then !!do suppression around x_i/a = fr_i so that the derivative is continuous; this is for stability of the ode solver; 
             suppression_factor =  sin(coop_pio2*(1.d0+arat(1)))**4 
             dyda(1) = dyda(1)*suppression_factor
             dyda(4) = dyda(4)*suppression_factor
          endif
       endif
       if(arat(2).lt.-1.d0 .and. y(5) .lt. 0.d0 )then !!collapsed; freeze it
          dyda(2) = 0.d0
          dyda(5) = 0.d0
       else !!still collapsing
          !!----------- equation for x_2 -------------------------
          dyda(2) = y(5) / dadt    !!d x_1/da = d x_1/dt /(da/dt)
          dyda(5) = y(2)/dadt/2.d0 * ( &   !! d( dx_1/dt)/da =(d^2 x_1/dt^2)/(da/dt)
               dark_Energy_term  &
               + radiation_term &  
               - rhombarby3*(1.d0 + delta *(1.d0+bprime(2)*1.5d0) + (3.d0*params%lambda(2)-sum(params%lambda))*growthD ) &  !!matter contribution
               )- dlnmdt*y(5)/dadt
          !!----------- end of equation for x_2 -------------------------
          if(arat(2) .lt. 0.d0)then !!do suppression around x_i/a = fr_i so that the derivative is continuous; this is for stability of the ode solver;
             suppression_factor =  sin(coop_pio2*(1.d0+arat(2)))**4 
             dyda(2) = dyda(2)*suppression_factor
             dyda(5) = dyda(5)*suppression_factor
          endif
       endif
       if(arat(3).lt.-1.d0 .and. y(6) .lt. 0.d0)then !!collapsed; freeze it
          dyda(3) = 0.d0
          dyda(6) = 0.d0
       else !!still collapsing
          !!----------- equation for x_3 -------------------------
          dyda(3) = y(6) / dadt    !!d x_1/da = d x_1/dt /(da/dt)
          dyda(6) = y(3)/dadt/2.d0 * ( &   !! d( dx_1/dt)/da =(d^2 x_1/dt^2)/(da/dt)
               dark_Energy_term  &
               + radiation_term &  
               -  rhombarby3*(1.d0 + delta *(1.d0+bprime(3)*1.5d0) + (3.d0*params%lambda(3)-sum(params%lambda))*growthD ) &  !!matter contribution
               ) - dlnmdt*y(6)/dadt      
          !!----------- end of equation for x_3 -------------------------
          if(arat(3) .lt. 0.d0)then  !!do suppression around x_i/a = fr_i so that the derivative is continuous; this is for stability of the ode solver;
             suppression_factor =  sin(coop_pio2*(1.d0+arat(3)))**4
             dyda(3) = dyda(3)*suppression_factor
             dyda(6) = dyda(6)*suppression_factor
          endif
       endif
    endif
    !!the scalar field
    if(y(7) .gt. 1.d-30)then
       dyda(7) = y(8)/dadt
       q =  O0_DE(params%cosmology)%cplde_intQofphi%derivative(y(7))
       dvdphi = O0_DE(params%cosmology)%cplde_Vofphi%derivative(y(7))
       dyda(8) = params%ddotphi(a, y(7), y(8), rho_m, q, dvdphi, dadt/a)/dadt
       dyda(7) = max(-y(7)/a*max_decay_rate, dyda(7))
    else
       dyda(7) = 1.d-29
       dyda(8) = 1.d-29
    endif
  end subroutine coop_coupledDE_collapse_odes

!!=====================You don't need to read anything below ==================

!!==================== methods of class coop_coupledDE_collapse_params ============
  subroutine coop_coupledDE_collapse_params_evolve(this, a, y, a_end)
  !!evolve y from a to a_end
  !!vector y = ( x_1, x_2,  x_3,  d x_1/dt, d x_2/dt, d x_3/dt )
  !!the object this contains all the parameters for the model
    class(coop_coupledDE_collapse_params)::this
    COOP_REAL::a, y(this%num_ode_vars), a_end
    COOP_REAL,parameter::tol = max(1.d-10, coop_coupledDE_collapse_accuracy*1.d-5)
    COOP_INT::ind
    COOP_REAL::c(24), w(this%num_ode_vars, 9)
    select type(this)
    type is(coop_coupledDE_collapse_params)
       ind = 1
       call coop_dverk_with_coupledDE_collapse_params(this%num_ode_vars, coop_coupledDE_collapse_odes, this, a, y, a_end, tol, ind, c, this%num_ode_vars, w)
    class default
       stop "Evolve: Extended class this has not been implemented"
    end select
  end subroutine coop_coupledDE_collapse_params_evolve

  subroutine coop_coupledDE_collapse_params_get_solution(this, a_arr, x_arr)
    !!get multiple snapshots of of x 
    !!input a_arr(1:n) is an array of scale factors in ascending order;
    !! return the result in x_arr(1:3, 1:n), where x_arr(1:3, i) is the solution of (x_1, x_2, x_3) at scale factor a_arr(i).    
    class(coop_coupledDE_collapse_params)::this
    COOP_REAL::a_arr(:), x_arr(:,:), y(8), a
    COOP_REAL,parameter::tol = max(1.d-10, coop_coupledDE_collapse_accuracy*1.d-5)
    COOP_INT::ind, i, n, m
    COOP_REAL::c(24), w(this%num_ode_vars, 9)
    n = size(a_arr)
    if(size(x_arr, 2) .ne. n) stop "get_solution: input a_arr and x_arr have different sizes"
    m = min(size(x_arr, 1), 8)  !!normally m = 3 but not always
    select type(this)
    type is(coop_coupledDE_collapse_params)
       ind = 1
       if(this%normalize_to_early)then
          a = 0.02d0
       else
          a = min(max(0.005d0, min(0.05d0, 100.d0*coop_coupledDE_collapse_accuracy))/maxval(abs(this%lambda)), a_arr(1)*0.99d0, 0.03d0)
       endif
       call this%set_initial_conditions(a, y)
       do i=1, n
          call coop_dverk_with_coupledDE_collapse_params(this%num_ode_vars, coop_coupledDE_collapse_odes, this, a, y, a_arr(i), tol, ind, c, this%num_ode_vars, w)
          x_arr(1:m, i) = y(1:m)
          if(size(x_arr,1).ge.9)then
             x_arr(9, i) = this%cosmology%Hratio(a_arr(i))*a_arr(i)**1.5
          endif
          if(size(x_arr,1).ge.10)then
             x_arr(10, i) = this%Growth_D(a_arr(i))/a_arr(i)
          endif
       enddo
    class default
       stop "Evolve: Extended class this has not been implemented"
    end select
  end subroutine coop_coupledDE_collapse_params_get_solution

!!return the redshift where the last virialized axis collapses
  function coop_coupledDE_collapse_params_zvir1(this) result(zvir1)
    class(coop_coupledDE_collapse_params)::this
    COOP_REAL::a, a_last, a_next, y(this%num_ode_vars), Frho, ycopy(this%num_ode_vars), incr, zvir1
    COOP_REAL,parameter::tol = max(1.d-10, coop_coupledDE_collapse_accuracy*1.d-5)
    COOP_INT::ind
    COOP_REAL::c(24), w(this%num_ode_vars, 9)
    COOP_INT::indcopy
    COOP_REAL::ccopy(24), wcopy(this%num_ode_vars, 9)
    select type(this)
    type is(coop_coupledDE_collapse_params)
       ind = 1
       Frho = sum(this%lambda)
       if(Frho .lt. 1.d0 .or. this%lambda(1) .gt. this%lambda(2) .or. this%lambda(2) .gt. this%lambda(3) .or. this%lambda(1) .lt. -50.d0 )then !!bad input or obviously won't collapse
          zvir1 = -1.d0
          return
       endif
       if(this%normalize_to_early)then
          a = 0.02d0
       else
          a = min(max(0.005d0, min(0.05d0, 50.d0*coop_coupledDE_collapse_accuracy)/maxval(abs(this%lambda))), 0.03d0)
       endif
       call this%set_initial_conditions(a, y)
       a_next = 0.1d0/Frho
       call coop_dverk_with_coupledDE_collapse_params(this%num_ode_vars, coop_coupledDE_collapse_odes, this, a, y, a_next, tol/10.d0, ind, c, this%num_ode_vars, w)
       incr = 0.05d0
       do 
          a_last = a
          indcopy = ind
          ccopy = c
          wcopy = w
          ycopy = y
          a_next = min(a*(1.d0+incr), 1.d0)
          call coop_dverk_with_coupledDE_collapse_params(this%num_ode_vars, coop_coupledDE_collapse_odes, this, a, y, a_next, tol, ind, c, this%num_ode_vars, w)
          if(y(1)/a .lt.  this%collapse_a_ratio(1)  .and. y(4) .le. 0.d0)then
             if(incr .lt. coop_coupledDE_collapse_accuracy)then
                zvir1 = 2.d0/(a+a_last) - 1.d0
                return
             endif
             a = a_last
             ind = indcopy
             c = ccopy
             w = wcopy
             y = ycopy
             incr = incr/2.d0
             cycle
          endif
          if(a_next .gt. 1.d0-coop_coupledDE_collapse_accuracy/2.d0)then
             zvir1 =  coop_coupledDE_collapse_bad_zvir
             return
          endif
       enddo
    class default
       stop "Evolve: Extended class this has not been implemented"
    end select
  end function coop_coupledDE_collapse_params_zvir1



  subroutine coop_coupledDE_collapse_params_init(this, params, update_cosmology)
    class(coop_coupledDE_collapse_params)::this
    type(coop_dictionary)::params
    logical, optional::update_cosmology
    logical::success
    COOP_INT::i
    COOP_REAL::F_pk, e_nu, p_nu
    call coop_dictionary_lookup(params, "normalize_to_early", this%normalize_to_early, .false.)
    !set up cosmology
    if(present(update_cosmology))then
       if(update_cosmology)then
          call this%cosmology%init_from_dictionary(params, level = coop_init_level_set_pert, success=success)
          if(.not. success)then
             stop "cannot initialize cosmology"
          else
             write(*,*) "cosmology is initialized"
          endif
       endif
    endif
    
    !set up lambda
    do i = 1, 3
       call coop_dictionary_lookup(params, "lambda"//COOP_STR_OF(i), this%lambda(i), -1.1d30)
    enddo
    if(any(this%lambda .lt. -1.d30))then
       call coop_dictionary_lookup(params, "collapse_fpk", F_pk)
       call coop_dictionary_lookup(params, "collapse_e", e_nu)
       call coop_dictionary_lookup(params, "collapse_p", p_nu)
       this%lambda(1) = (F_pk/3.d0)*(1.d0 - 3.d0*e_nu + p_nu)
       this%lambda(2) = (F_pk/3.d0)*(1.d0 - 2.d0*p_nu)
       this%lambda(3) = (F_pk/3.d0)*(1.d0 + 3.d0*e_nu + p_nu)
    endif
    if(.not. this%lambda(3) .ge. 0.d0) stop "error: lambda_3 must be positive"
    call coop_dictionary_lookup(params, "collapse_ratio_1",  this%collapse_a_ratio(1))
    call coop_dictionary_lookup(params, "collapse_ratio_2",  this%collapse_a_ratio(1))
    call coop_dictionary_lookup(params, "collapse_ratio_3",  this%collapse_a_ratio(1))
    this%collapse_a_ratio = max(this%collapse_a_ratio, 0.01d0)
    this%is_spherical = abs(this%lambda(1) - this%lambda(2)) .lt. 1.d-8 .and. abs(this%lambda(1) - this%lambda(3)) .lt. 1.d-8 .and. abs(this%collapse_a_ratio(1) - this%collapse_a_ratio(2)) .lt. 1.d-3 .and. abs(this%collapse_a_ratio(1)-this%collapse_a_ratio(3)) .lt. 1.d-3

    !!set b' function
    if(.not. this%bprime%initialized)call this%bprime%init_symmetric(f = coop_coupledDE_collapse_bprime_reduced, nx = 501, xmin = 1.d-6, xmax = 1.d6, xlog = .true., name = "BPRIME")
  end subroutine coop_coupledDE_collapse_params_init


  subroutine coop_coupledDE_collapse_params_get_bprime(this, x, bprime)
    class(coop_coupledDE_collapse_params)::this
    COOP_REAL,intent(IN)::x(3)
    COOP_REAL,intent(OUT)::bprime(3)
    bprime(1) = this%bprime%eval( (x(1)/x(2))**2,  (x(1)/x(3))**2 )
    bprime(2) = this%bprime%eval( (x(2)/x(1))**2,  (x(2)/x(3))**2 )
    bprime(3) = this%bprime%eval( (x(3)/x(2))**2,  (x(3)/x(1))**2 )
  end subroutine coop_coupledDE_collapse_params_get_bprime

  !!D(a)
  function coop_coupledDE_collapse_params_Growth_D(this, a) result(D)
    class(coop_coupledDE_collapse_params)::this
    COOP_REAL::a, D
    D = this%cosmology%growth_of_z(1.d0/a-1.d0)/this%cosmology%growth_of_z(0.d0)
  end function coop_coupledDE_collapse_params_Growth_D

  !!d ln D/dt 
  function coop_coupledDE_collapse_params_Growth_H_D(this, a) result(H_D)
    class(coop_coupledDE_collapse_params)::this
    COOP_REAL::a, H_D
    H_D = this%cosmology%fgrowth_of_z(1.d0/a-1.d0)*this%cosmology%Hratio(a)
  end function coop_coupledDE_collapse_params_Growth_H_D
  
  !! H a / (H_0 a_0)
  function coop_coupledDE_collapse_params_aH(this, a) result(aH)
    class(coop_coupledDE_collapse_params)::this
    COOP_REAL::a, aH
    aH = this%cosmology%aHratio(a)
  end function coop_coupledDE_collapse_params_aH

!!return H_0^2  \ddot a / a
  function coop_coupledDE_collapse_params_ddotabya(this, a) result(ddotabya)
    class(coop_coupledDE_collapse_params)::this
    COOP_REAL::a, ddotabya
    ddotabya = (1.d0+this%cosmology%HdotbyHsq(a))*this%cosmology%H2a4(a)/a**4
  end function coop_coupledDE_collapse_params_ddotabya

  subroutine coop_coupledDE_collapse_params_free(this)
    class(coop_coupledDE_collapse_params)::this
    call this%bprime%free()
    call this%cosmology%free()
  end subroutine coop_coupledDE_collapse_params_free


!!=========================================================================
!!utilities; no need to understand or change anything below
  subroutine coop_dverk_with_coupledDE_collapse_params(n, fcn, params, x, y, xend, tol, ind, c, nw, w)
    type(coop_coupledDE_collapse_params) params
#define DVERK_ARGUMENTS ,params
#include "dverk.h"    
#undef DVERK_ARGUMENTS
  end subroutine coop_dverk_with_coupledDE_collapse_params


  function coop_coupledDE_collapse_bprime_reduced(lambda1, lambda2) result(bprime)
    COOP_REAL::lambda1, lambda2, bprime
    if(abs(lambda1-1.d0) .lt. 1.d-4 .and. abs(lambda2-1.d0) .lt. 1.d-4)then !!both close to 0
       bprime = (2.d0/15.d0)*(2.d0-lambda1 - lambda2) + (3.d0/28.d0-0.05d0)*((1.d0-lambda1)**2+(1.d0-lambda2)**2+(2.d0/3.d0)*(1.d0-lambda1)*(1.d0-lambda2))
       return
    endif
    bprime = (coop_elliptic_Rd(1.d0/lambda1, 1.d0/lambda2, 1.d0)/sqrt(lambda1*lambda2)-1.d0)*2.d0/3.d0
  end function coop_coupledDE_collapse_bprime_reduced

  function coop_coupledDE_collapse_params_ddotphi(this, a, phi, phidot, rho_m, Q, dVdphi, H) result(ddotphi)
    class(coop_coupledDE_collapse_params)::this
    COOP_REAL::a, phi, phidot, rho_m, Q, dVdphi, H
    COOP_REAL::ddotphi
    !!put some upper bound to avoid going to netative phi
    ddotphi =-3.d0*H*phidot - dVdphi - Q*rho_m
  end function coop_coupledDE_collapse_params_ddotphi

#else
#endif
end module coop_coupledDE_collapse_mod

