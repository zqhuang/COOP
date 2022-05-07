module coop_wl_mod
  use coop_wrapper_firstorder
  use coop_halofit_mod
  implicit none
#include "constants.h"

  type coop_wl_object
     COOP_STRING::name = ""
     COOP_INT::num_z = 0
     COOP_INT::num_z_bins = 0
     COOP_INT::num_theta_bins = 0
     COOP_INT:: num_mask = 0
     COOP_INT::num_z_p = 0
     COOP_INT::num_obs = 0
     COOP_REAL,dimension(:),allocatable::theta_bins, z_bins, xi_obs, xi_theory, z_p  
     COOP_INT,dimension(:),allocatable::mask_indices
     COOP_REAL,dimension(:,:),allocatable::p, wl_invcov, wl_cov !!p is the source galaxy distribution p(z, bin)
     COOP_REAL::ah_factor = 1.d0
     COOP_REAL::max_z, kmax
     logical::cut_theta  !!cut nonlinear scales
   contains
     procedure::Loglike => coop_wl_object_LogLike
     procedure::Init => Coop_WL_Object_Init
     procedure::free => coop_wl_object_free
     procedure::get_theory => coop_wl_object_get_theory
  end type coop_wl_object
  

contains

  subroutine coop_wl_object_free(this)
    class(coop_wl_object)::this
    if(allocated(this%theta_bins))deallocate(this%theta_bins)
    if(allocated(this%z_bins))deallocate(this%z_bins)
    if(allocated(this%xi_obs))deallocate(this%xi_obs)
    if(allocated(this%xi_theory))deallocate(this%xi_theory)    
    if(allocated(this%z_p))deallocate(this%z_p)
    if(allocated(this%p))deallocate(this%p)
    if(allocated(this%mask_indices))deallocate(this%mask_indices)
    this%num_z_bins = 0
    this%num_theta_bins = 0
    this%num_mask = 0
    this%num_z_p = 0
    this%num_obs = 0
  end subroutine coop_wl_object_free

  subroutine coop_wl_object_init(this, filename)
    class(coop_wl_object)::this
    COOP_UNKNOWN_STRING::filename
    type(coop_dictionary)::dict
    COOP_STRING :: measurements_file, window_file, measurements_format, cov_file, cut_file
    type(coop_file)::f
    COOP_REAL :: dummy1,dummy2,pnorm
    COOP_REAL, allocatable, dimension(:,:) :: temp
    COOP_REAL, allocatable, dimension(:,:) :: cut_values
    COOP_INT, allocatable, dimension(:) :: mask
    COOP_INT i,iz,it,ib,j,k,nt,izl,izh
    COOP_REAL :: xi_plus_cut, xi_minus_cut
    logical::success
    call this%free()
    call coop_load_dictionary(COOP_DATAPATH(filename), dict)
    call coop_dictionary_lookup(dict, "name", this%name, "Weak Lensing Data")
    call coop_dictionary_lookup(dict, "nz_wl", this%num_z, 100)
    call coop_dictionary_lookup(dict, "max_z", this%max_z)
    call coop_dictionary_lookup(dict, "cut_theta", this%cut_theta, .false.)
    call coop_dictionary_lookup(dict, "num_z_bins", this%num_z_bins, 1)
    call coop_dictionary_lookup(dict, "num_theta_bins", this%num_theta_bins)
    call coop_dictionary_lookup(dict, "num_z_p", this%num_z_p)
    nt = this%num_z_bins*(1+this%num_z_bins)/2
    allocate(this%theta_bins(this%num_theta_bins))
    allocate(this%z_bins(this%num_z_bins))
    this%num_obs = this%num_theta_bins*nt*2
    allocate(this%xi_obs(this%num_obs), this%xi_theory(this%num_obs))
    allocate(this%z_p(this%num_z_p))
    allocate(this%p(this%num_z_p,this%num_z_bins), this%wl_cov(this%num_obs, this%num_obs))
    call coop_dictionary_lookup(dict, "kmax", this%kmax)
    call coop_dictionary_lookup(dict, "ah_factor", this%ah_factor, 1.d0)
    call coop_dictionary_lookup(dict, "measurements_file", measurements_file)
    call coop_dictionary_lookup(dict, "measurements_format", measurements_format)
    if(trim(measurements_format).eq."")stop "WL_init: key measurement_format is not found"
    call coop_dictionary_lookup(dict, "window_file", window_file)
    call coop_dictionary_lookup(dict, "cov_file", cov_file)
    if(trim(cov_file) .eq. "") stop "WL_init: key cov_file is not found"
    call coop_import_matrix(COOP_DATAPATH(cov_file), this%wl_cov, this%num_obs, this%num_obs, success)
    if(.not. success) stop "WL_init: cov_file format incorrect"

    if (this%cut_theta) then
       call coop_dictionary_lookup(dict, 'cut_file', cut_file)
       if(trim(cut_file).eq."") stop "WL_init: key cut_file is not found"
       allocate(cut_values(this%num_z_bins, 2))
       call coop_import_matrix(COOP_DATAPATH(cut_file), cut_values, this%num_z_bins, 2, success)
       if(.not. success) stop "WL_init: cut_file format incorrect"
    end if

    do ib=1,this%num_z_bins
       call f%open(COOP_DATAPATH_I(window_file, ib), "r")
       do iz=1,this%num_z_p
          read (f%unit,*) this%z_p(iz),this%p(iz,ib)
       enddo
       call f%Close()
    enddo

    call f%open(COOP_DATAPATH(measurements_file))
    if(measurements_format == '1bin') then
       do it = 1,this%num_theta_bins
          read (f%unit,*) this%theta_bins(it),this%xi_obs(it),dummy1,this%xi_obs(it+this%num_theta_bins),dummy2
       end do
    elseif (measurements_format == 'nbin') then
       k = 1
       allocate(temp(2*this%num_theta_bins,nt))
       do i=1,2*this%num_theta_bins
          read (F%unit,*) dummy1,temp(i,:)
          if (i.le.this%num_theta_bins) this%theta_bins(i)=dummy1
       end do
       do j=1,nt
          do i=1,2*this%num_theta_bins
             this%xi_obs(k) = temp(i,j)
             k = k + 1
          end do
       end do
       deallocate(temp)
    else
       write(*,*) 'ERROR: Not yet implemented WL measurements format: '//trim(measurements_format)
       stop
    end if
    call f%Close()

    !Normalize window functions p so \int p(z) dz = 1
    do ib=1,this%num_z_bins
       pnorm = 0
       do iz=2,this%num_z_p
          pnorm = pnorm + 0.5d0*(this%p(iz-1,ib)+this%p(iz,ib))*(this%z_p(iz)-this%z_p(iz-1))
       end do
       this%p(:,ib) = this%p(:,ib)/pnorm
    end do
    
    ! Apply Anderson-Hartap correction
    this%wl_cov = this%wl_cov/this%ah_factor
    
    ! Compute theta mask
    allocate(mask(this%num_theta_bins*nt*2))
    if (allocated(cut_values)) then
       mask = 0
       iz = 0
       do izl = 1,this%num_z_bins
          do izh = izl,this%num_z_bins
             iz = iz + 1 ! this counts the bin combinations iz=1 =>(1,1), iz=1 =>(1,2) etc
             do i = 1,this%num_theta_bins
                j = (iz-1)*2*this%num_theta_bins
                xi_plus_cut = max(cut_values(izl,1),cut_values(izh,1))
                xi_minus_cut = max(cut_values(izl,2),cut_values(izh,2))
                if (this%theta_bins(i)>xi_plus_cut) mask(j+i) = 1
                if (this%theta_bins(i)>xi_minus_cut) mask(this%num_theta_bins + j+i) = 1
                ! Testing
                !write(*,'(5i4,3E15.3,2i4)') izl,izh,i,i+j,this%num_theta_bins + j+i,xi_plus_cut,&
                !     xi_minus_cut,this%theta_bins(i),mask(j+i),mask(this%num_theta_bins + j+i)
             end do
          end do
       end do
    else
       mask = 1
    end if
    this%num_mask = sum(mask)
    allocate(this%mask_indices(this%num_mask))
    j = 1
    do i=1,this%num_theta_bins*nt*2
       if (mask(i) == 1) then
          this%mask_indices(j) = i
          j = j+1
       end if
    end do
    
    ! Precompute masked inverse
    allocate(this%wl_invcov(this%num_mask,this%num_mask))
    this%wl_invcov = this%wl_cov(this%mask_indices,this%mask_indices)
    call coop_matrix_Inverse(this%wl_invcov)
    
  end subroutine coop_wl_object_init

  function coop_wl_object_LogLike(this, cosmology) result(loglike)
    class(coop_wl_object)::this
    type(coop_cosmology_firstorder)::cosmology
    COOP_REAL::loglike
    COOP_REAL::vec(this%num_mask)
    loglike = 0.d0
    if(this%num_obs .eq. 0) return
    call this%get_theory(cosmology)
    vec = this%xi_theory(this%mask_indices) - this%xi_obs(this%mask_indices)
    loglike = dot_product(vec, matmul(this%wl_invcov, vec))/2.d0
  end function coop_wl_object_LogLike


  subroutine coop_wl_object_get_theory(this, cosmology)
    class(coop_wl_object)::this    
    type(coop_cosmology_firstorder)::cosmology
    COOP_INT::num_l, num_theta, i_l, i_theta, ib, jb, ncross, icross, iz, izl, izh, i, j
    COOP_REAL, parameter::dlntheta = 0.2d0, dlnl = 0.15d0, xstop = 200.d0, dx = 0.02d0
    COOP_REAL,dimension(:),allocatable::l, lnl, theta, lntheta
    COOP_REAL,dimension(:,:),allocatable::Cl, lnl2Cl, lnl2Cl2
    COOP_REAL::theta_min, theta_max, l_min, l_max, xmin, xmax, x, bessel0, bessel4, Cval
    COOP_REAL::lp, i1, i2, khmin, khmax, lll
    
    COOP_REAL, dimension(:),allocatable::i1p, i2p
    COOP_REAL, dimension(:,:),allocatable::xi1, xi2, xi1sp, xi2sp
    if(this%num_mask.le.0)return !!nothing to calculate
    ncross = this%num_z_bins*(this%num_z_bins+1)/2
    !!compute the convergence power spectrum
    theta_min = minval(this%theta_bins)*0.99
    theta_max = maxval(this%theta_bins)*1.01
    
    l_min = max(0.05d0/(theta_max*coop_SI_arcmin), 1.d0)
    l_max = max(xstop/(theta_min*coop_SI_arcmin), 10.d0, l_min*1.5d0)
    
    num_l = ceiling(log(l_max/l_min)/dlnl)
    num_theta = ceiling(log(theta_max/theta_min)/ dlntheta)
    
    allocate(l(num_l), lnl(num_l), Cl(num_l, ncross), theta(num_theta), lntheta(num_theta), lnl2Cl(num_l, ncross), lnl2Cl2(num_l, ncross))
    
    allocate(xi1(num_theta,ncross),xi2(num_theta,ncross), xi1sp(num_theta, ncross), xi2sp(num_theta, ncross))
    
    allocate(i1p(ncross), i2p(ncross))
    
    call coop_set_uniform(num_l, lnl, log(l_min), log(l_max))
    l = exp(lnl)
    call coop_set_uniform(num_theta, lntheta, log(theta_min), log(theta_max))
    theta = exp(lntheta)
    call coop_halofit_get_weaklensing_power(cosmology, this%num_z_bins, this%num_z_p, this%z_p, this%p, num_l, l, Cl)
    do icross = 1, ncross
       lnl2Cl(:,icross) = log(Cl(:, icross)*l**2 + 1.d-15)
    enddo
    do icross = 1, ncross
       call coop_spline(num_l, lnl, lnl2Cl(:, icross), lnl2Cl2(:, icross))
    enddo
    xi1 = 0.
    xi2= 0.
    do i_theta = 1, num_theta
       xmin = l_min*theta(i_theta)*coop_SI_arcmin
       xmax = l_max*theta(i_theta)*coop_SI_arcmin
       x = xmin
       lp = 0
       i1p = 0
       i2p = 0
       do while(x .lt. xstop .and. x .lt. xmax)
          lll=x/(theta(i_theta)*coop_SI_arcmin)
          Bessel0 = Bessel_J0(x)
          Bessel4 = Bessel_JN(4,x)
          do ib=1,this%num_z_bins
             do jb=1, ib
                icross = COOP_MATSYM_INDEX(this%num_z_bins, ib, jb)
                !interpolate
                call coop_splint(num_l, lnl, lnl2Cl, lnl2Cl2, log(lll), Cval)
                Cval = exp(Cval)/lll
                i1 = Cval*Bessel0
                i2 = Cval*Bessel4
                xi1(i_theta,icross) = xi1(i_theta,icross)+0.5*(i1p(icross)+i1)*(lll-lp)
                xi2(i_theta,icross) = xi2(i_theta,icross)+0.5*(i2p(icross)+i2)*(lll-lp)
                i1p(icross) = i1
                i2p(icross) = i2
             end do
          end do
          x = x+dx
          lp = lll
       end do
    enddo
    xi1 = xi1/coop_2pi
    xi2 = xi2/coop_2pi
    do icross = 1, ncross
       call coop_spline(num_theta, lntheta, xi1(:, icross), xi1sp(:, icross))
       call coop_spline(num_theta, lntheta, xi2(:, icross), xi2sp(:, icross))
    enddo

    iz = 0
    do izl = 1,this%num_z_bins
       do izh = izl,this%num_z_bins
          iz = iz + 1 ! this counts the bin combinations iz=1 =>(1,1), iz=1 =>(1,2) etc
          icross = COOP_MATSYM_INDEX(this%num_z_bins, izl, izh)
          do i = 1,this%num_theta_bins
             j = (iz-1)*2*this%num_theta_bins
             Call coop_splint(num_theta, lntheta, xi1(:, icross), xi1sp(:, icross), log(this%theta_bins(i)),  this%xi_theory(j+i))
             Call coop_splint(num_theta, lntheta, xi2(:, icross), xi2sp(:, icross), log(this%theta_bins(i)),  this%xi_theory(this%num_theta_bins+j+i))
          end do
       end do
    end do
    
  end subroutine coop_wl_object_get_theory
  
end module coop_wl_mod
