module solvep
  use ode_utils
  implicit none

  private

  real(dl),parameter::Mpl = 1024.d0
  real(dl):: spw_mphi = 5.8d-6 * Mpl
  real(dl):: spw_phi0 = 16. * Mpl
  real(dl):: spw_deltaphi = 0.15 * Mpl
  real(dl):: spw_amp = 0.1d0
  integer,parameter::ntbl = 2048+1
  integer,parameter::ipivot = (ntbl+1)/2
  real(dl),dimension(ntbl)::lna_tbl, fric_tbl, fric2_tbl, h_tbl, h2_tbl, intfric_tbl, intfric2_tbl, tenfric_tbl, tenfric2_tbl, tenintfric_tbl, tenintfric2_tbl
  real(dl),parameter::lna_start = 0.d0
  real(dl),parameter::lna_end = 30.d0
  real(dl),parameter::lna_pivot = (lna_start+lna_end)/2.d0
  real(dl)::bgy(3, ntbl), bgyp(3, ntbl)


  public spw_obtain_power
  

contains


  

  subroutine spw_obtain_power(m, ps, pt, lnkmin, lnkmax, mphi, nefolds, bump_dn, bump_amp, bump_width )
    integer m
    real(dl) mphi, nefolds, bump_dn, bump_amp, bump_width
    real(dl)   ps(m), pt(m), lnk(m), phi_start, lnkmin, lnkmax
    integer i
    call set_uniform(ntbl, lna_tbl, lna_start, lna_end)
    !!get the background  
    phi_start  = sqrt(4.d0*(nefolds + lna_pivot))*Mpl
    spw_mphi = mphi*Mpl
    spw_amp = bump_amp
    spw_phi0 = sqrt(4.d0*(nefolds+bump_dn))*Mpl
    spw_deltaphi = bump_width*Mpl
    bgy(2, 1) = phi_start 
    bgy(1, 1) = sqrt(potential(phi_start)/3.d0/Mpl**2)
    do i=1, 3
       bgy(3,1) = -dVdphi(phi_start)/3.d0/bgy(1,1)
       bgy(1, 1) = sqrt((potential(phi_start)+bgy(3,1)**2/2.d0 )/3.d0/Mpl**2)    
    enddo
    call genericdverk(bg_eq, lna_tbl, bgy, 1.d-10)
    do i=1, ntbl
       call bg_eq(3, lna_tbl(i), bgy(:, i), bgyp(:, i))
       h_tbl(i) = bgy(1, i)
       bgyp(:, i) = bgyp(:, i)*h_tbl(i)  !!convert to dot y
       fric_tbl(i) = 3.d0*h_tbl(i) - 2.d0*bgyp(1, i)/h_tbl(i) + 2.d0*bgyp(3,i)/bgy(3, i)
       intfric_tbl(i) = 6.d0*lna_tbl(i) - 4.d0*log(abs(h_tbl(i)/bgy(3,i)))
       tenfric_tbl(i) = 3.d0*h_tbl(i)
       tenintfric_tbl(i) = 6.d0*lna_tbl(i)
    enddo
    call splines(lna_tbl, intfric_tbl, intfric2_tbl)
    call splines(lna_tbl, fric_tbl, fric2_tbl)
    call splines(lna_tbl, tenintfric_tbl, tenintfric2_tbl)
    call splines(lna_tbl, tenfric_tbl, tenfric2_tbl)
    call splines(lna_tbl, h_tbl, h2_tbl)

  !!solve the perturbations
    call set_uniform(m, lnk, lnkmin+log(h_tbl(ipivot))+lna_tbl(ipivot), lnkmax+log(h_tbl(ipivot))+lna_tbl(ipivot))
    !$omp parallel do
    do i=1, m
       call get_scal_power(exp(lnk(i)), ps(i))
       call get_tens_power(exp(lnk(i)), pt(i))
    enddo
    !$omp end parallel do
  end subroutine spw_obtain_power


  subroutine get_scal_power(k, pk)
    type(params_pass) params
    integer istart, iend
    real(dl) k, pk,  y(2), lna, lna_end, c(24), w(2, 9)
    integer ind
    istart = 1
    do while(k/exp(lna_tbl(istart))/h_tbl(istart) .gt. 80.d0)
       istart = istart  + 1
       if(istart .ge. ntbl) stop "k too big"
    enddo
    iend = istart
    do while(k/exp(lna_tbl(iend))/h_tbl(iend) .gt. 0.02d0)
       iend = iend + 1
       if(iend .ge. ntbl) stop "k too big"
    enddo
    y(1) = (k/const_2pi) * h_tbl(istart)/bgy(3, istart)/exp(lna_tbl(istart))
    y(2) = y(1) * ( bgyp(1,istart)/bgy(1, istart) - bgyp(3, istart)/bgy(3, istart) - h_tbl(istart) )
    params%fp(1) = k**2
    params%fp(2) = 2.d0* (log((k/exp(lna_tbl(istart)))*y(1)**2)) + intfric_tbl(istart)
    lna = lna_tbl(istart)
    ind = 1
    c = 0
    w = 0
    call dverk_with_params(params, scal_eq, lna, y, lna_tbl(iend), ind, c, w, 1.d-9)
    pk = y(1)**2 
  end subroutine get_scal_power


  subroutine get_tens_power(k, pk)
    type(params_pass) params
    integer istart, iend
    real(dl) k, pk,  y(2), lna, lna_end, c(24), w(2, 9)
    integer ind
    istart = 1
    do while(k/exp(lna_tbl(istart))/h_tbl(istart) .gt. 60.d0)
       istart = istart  + 1
       if(istart .ge. ntbl) stop "k too big"
    enddo
    iend = istart
    do while(k/exp(lna_tbl(iend))/h_tbl(iend) .gt. 0.02d0)
       iend = iend + 1
       if(iend .ge. ntbl) stop "k too big"
    enddo
    y(1) = (k*const_sqrt2/const_pi)/exp(lna_tbl(istart))/Mpl
    y(2) = -y(1) * h_tbl(istart)
    params%fp(1) = k**2
    params%fp(2) = 2.d0* (log((k/exp(lna_tbl(istart)))*y(1)**2)) + tenintfric_tbl(istart)
    lna = lna_tbl(istart)
    ind = 1
    c = 0
    w = 0
    call dverk_with_params(params, tens_eq, lna, y, lna_tbl(iend), ind, c, w, 1.d-7)
    pk = y(1)**2 
  end subroutine get_tens_power


  subroutine bg_eq(n, lna, y, yprime)
!!y(1)  H
!!y(2)  phi
!!y(3)  dot phi
!!y(4) \int f_1(t) dt, where f_1(t) =  3 H - 2 \dot H/H + 2 \ddot \phi /\dot\phi
    type(params_pass) params
    integer n
    real(dl) lna, y(n), yprime(n), V, K
    V = potential(y(2))
    K = y(3)**2/2.d0
    yprime(1) = (-(4.d0*K - 2.d0*V)/6.d0/Mpl**2 - y(1)**2)/y(1)
    yprime(2) = y(3)/y(1)
    yprime(3) = -dVdphi(y(2))/y(1) - 3.d0*y(3)
  end subroutine bg_eq


  subroutine scal_eq(params, n, lna, y, yprime)
    type(params_pass) params
    integer n
    real(dl) lna, y(n), yprime(n), h
!!params%fp(1) k^2
!!params%fp(2) initial dot theta ^2 A^4 
!!params%fp(3) initial intfric
!!y(1) A
!!y(2) dot A
    h = hubble(lna)
    yprime(1) = y(2)/h
    yprime(2) = (-fric(lna)*y(2) - (params%fp(1)/exp(2.d0*lna))*y(1) + exp(params%fp(2) - intfric(lna))/y(1)**3)/h
  end subroutine scal_eq


  subroutine tens_eq(params, n, lna, y, yprime)
    type(params_pass) params
    integer n
    real(dl) lna, y(n), yprime(n), h
!!params%fp(1) k^2
!!params%fp(2) initial dot theta ^2 A^4 
!!params%fp(3) initial intfric
!!y(1) A
!!y(2) dot A
    h = hubble(lna)
    yprime(1) = y(2)/h
    yprime(2) = (-tenfric(lna)*y(2) - (params%fp(1)/exp(2.d0*lna))*y(1) + exp(params%fp(2) - tenintfric(lna))/y(1)**3)/h
  end subroutine tens_eq





  function potential(phi) result(V)
    real(dl) phi, V
    V = (spw_mphi**2/2.d0)*phi**2 *( 1.d0 + spw_amp*(1.d0+tanh((phi-spw_phi0)/spw_deltaphi))/2.d0)
  end function potential
  
  function dVdphi(phi) 
    real(dl) phi, dVdphi
    dVdphi = (spw_mphi**2*phi) *( 1.d0 + spw_amp*(1.d0+tanh((phi-spw_phi0)/spw_deltaphi))/2.d0 + phi/4.d0/spw_deltaphi * spw_amp/cosh((phi-spw_phi0)/spw_deltaphi)**2)
  end function dVdphi

  function fric(lna)
    real(dl) lna, fric
    call splints(lna_tbl, fric_tbl, fric2_tbl, lna, fric)
  end function fric

  function intfric(lna)
    real(dl) lna, intfric
    call splints(lna_tbl, intfric_tbl, intfric2_tbl, lna, intfric)
  end function intfric

  function tenfric(lna)
    real(dl) lna, tenfric
    call splints(lna_tbl, tenfric_tbl, tenfric2_tbl, lna, tenfric)
  end function tenfric

  function tenintfric(lna)
    real(dl) lna, tenintfric
    call splints(lna_tbl, tenintfric_tbl, tenintfric2_tbl, lna, tenintfric)
  end function tenintfric


  function hubble(lna) 
    real(dl) lna, hubble
    call splints(lna_tbl, h_tbl, h2_tbl, lna, hubble)
  end function hubble

end module solvep
