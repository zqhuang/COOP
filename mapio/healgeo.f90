module coop_healpix_mod
!!I always assume ring order
  use coop_wrapper_firstorder
  use coop_sphere_mod
#ifdef HAS_HEALPIX
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
  use udgrade_nr
#endif
  implicit none

  
#include "constants.h"

  private

  public::coop_healpix_maps, coop_healpix_disc, coop_healpix_patch, coop_healpix_split,  coop_healpix_plot_spots,  coop_healpix_inpainting, coop_healpix_smooth_maskfile, coop_healpix_output_map, coop_healpix_get_disc, coop_healpix_export_spots, coop_healpix_smooth_mapfile, coop_healpix_patch_get_fr0, coop_healpix_lb2ang, coop_healpix_ang2lb, coop_healpix_fetch_patch, coop_healpix_mask_tol,  coop_healpix_mask_hemisphere, coop_healpix_index_TT,  coop_healpix_index_EE,  coop_healpix_index_BB,  coop_healpix_index_TE,  coop_healpix_index_TB,  coop_healpix_index_EB, coop_healpix_flip_mask, coop_healpix_diffuse_into_mask, coop_healpix_alm_check_done, coop_healpix_want_cls, coop_healpix_default_lmax, coop_healpix_max_threshold
  
  logical::coop_healpix_alm_check_done = .false.
  logical::coop_healpix_want_cls = .true.
  COOP_REAL,parameter::coop_healpix_max_threshold = 10.d0

  COOP_INT,parameter::dlc = kind( (1.d0,1.d0) )
  COOP_INT,parameter::coop_inpainting_lowl_max = 20
  COOP_INT,parameter::coop_inpainting_lowl_min = 5

  COOP_INT, parameter::coop_healpix_default_lmax=2500
  COOP_REAL,parameter::coop_healpix_mask_tol = 0.95  !!default mask tolerance
  COOP_INT::coop_healpix_inpainting_lowl=5
  COOP_REAL,parameter::coop_healpix_diffuse_scale = 10.d0*coop_SI_arcmin
  COOP_INT,parameter::coop_healpix_index_TT = 1
  COOP_INT,parameter::coop_healpix_index_EE = 2
  COOP_INT,parameter::coop_healpix_index_BB = 3
  COOP_INT,parameter::coop_healpix_index_TE = 4
  COOP_INT,parameter::coop_healpix_index_EB = 5
  COOP_INT,parameter::coop_healpix_index_TB = 6

  type, extends(coop_sphere_disc):: coop_healpix_disc
     COOP_INT nside
     COOP_INT center
   contains
     procedure :: pix2ang => coop_healpix_disc_pix2ang
     procedure :: ang2pix => coop_healpix_disc_ang2pix
     procedure :: pix2xy => coop_healpix_disc_pix2xy
     procedure :: xy2pix => coop_healpix_disc_xy2pix
  end type coop_healpix_disc

  type coop_healpix_maps
     COOP_INT npix, nside, nmaps, ordering, lmax, iq, iu, mask_npix, maskpol_npix
     character(LEN=80),dimension(64)::header
     COOP_INT,dimension(:),allocatable::spin
     COOP_SINGLE, dimension(:,:),allocatable::map
     complex, dimension(:,:,:),allocatable::alm
     COOP_SINGLE, dimension(:,:),allocatable::Cl
     COOP_INT,dimension(:),allocatable::mask_listpix, maskpol_listpix
     COOP_REAL chisq, mcmc_temperature
     logical,dimension(:),allocatable::alm_done
     COOP_SINGLE,dimension(:,:),allocatable::checksum
   contains     
     procedure :: init => coop_healpix_maps_init
     procedure :: free => coop_healpix_maps_free
     procedure :: write => coop_healpix_maps_write
     procedure :: read => coop_healpix_maps_read
     procedure :: extend => coop_healpix_maps_extend
     procedure :: import => coop_healpix_maps_import
     procedure :: open => coop_healpix_maps_read
     procedure :: mask => coop_healpix_maps_mask
     procedure :: allocate_alms => coop_healpix_maps_allocate_alms
     procedure :: get_cls =>   coop_healpix_maps_get_cls
     procedure :: get_fullcls => coop_healpix_maps_get_fullcls
     procedure :: map2alm => coop_healpix_maps_map2alm
     procedure :: alm2map => coop_healpix_maps_alm2map
     procedure :: udgrade => coop_healpix_maps_udgrade
     procedure :: simulate => coop_healpix_maps_simulate
     procedure :: simulate_Tmaps => coop_healpix_maps_simulate_Tmaps
     procedure :: simulate_TQUmaps => coop_healpix_maps_simulate_TQUmaps
     procedure :: iqu2TEB => coop_healpix_maps_iqu2TEB
     procedure :: qu2EB => coop_healpix_maps_qu2EB
     procedure :: iqu2LapTEB => coop_healpix_maps_iqu2LapTEB
     procedure :: teb2iqu => coop_healpix_maps_teb2iqu
     procedure :: iqu2TQTUT => coop_healpix_maps_iqu2TQTUT
     procedure :: smooth => coop_healpix_maps_smooth
     procedure :: smooth_with_window => coop_healpix_maps_smooth_with_window
     procedure :: smooth_mask => coop_healpix_smooth_mask
     procedure :: t2zeta => coop_healpix_maps_t2zeta
     procedure :: trim_mask => coop_healpix_trim_mask
     procedure :: convert2nested => coop_healpix_convert_to_nested
     procedure :: convert2ring => coop_healpix_convert_to_ring
     procedure :: filter_alm =>  coop_healpix_filter_alm
     procedure :: get_spots => coop_healpix_maps_get_spots
     procedure :: get_listpix => coop_healpix_maps_get_listpix
     procedure :: stack =>     coop_healpix_maps_stack
     procedure :: multstack =>     coop_healpix_maps_multstack
     procedure :: stack_with_covariance => coop_healpix_maps_stack_with_covariance
     procedure :: stack_with_listpix => coop_healpix_maps_stack_with_listpix
     procedure :: stack_north_south => coop_healpix_maps_stack_north_south
  end type coop_healpix_maps

  type coop_healpix_patch
     COOP_STRING::caption
     COOP_SHORT_STRING,dimension(:),allocatable::label
     COOP_SHORT_STRING::color_table
     COOP_SHORT_STRING:: genre
     COOP_INT::n, mmax, nmaps, npix, nstack_raw
     COOP_REAL::dr, zmin, zmax
     logical::headless_vectors
     COOP_REAL,dimension(:,:,:),allocatable::image
     COOP_REAL,dimension(:),allocatable::r
     COOP_REAL,dimension(:,:,:),allocatable::fr
     COOP_REAL,dimension(:,:,:),allocatable::wcm
     COOP_REAL,dimension(:,:,:),allocatable::wsm
     COOP_INT, dimension(:,:,:),allocatable::icm
     COOP_REAL,dimension(:,:),allocatable::nstack, indisk
     COOP_REAL::num_indisk_tol
   contains
     procedure::free => coop_healpix_patch_free
     procedure::init => coop_healpix_patch_init
     procedure::get_radial_profile => coop_healpix_patch_get_radial_profile
     procedure::get_all_radial_profiles => coop_healpix_patch_get_all_radial_profiles
     procedure::plot => coop_healpix_patch_plot
     procedure::plot_fft => coop_healpix_patch_plot_fft
  end type coop_healpix_patch
  


#define COS2RADIUS(cosx) (sqrt(2.d0*(1.d0 - (cosx))))
#define RADIUS2COS(r)  (1.d0-(r)**2/2.d0)

contains

  subroutine coop_healpix_diffuse_into_mask(this, mask, smoothscale, pol)  !!lambda<1
    real,parameter::nefolds = 2.
    type(coop_healpix_maps) mask, this
    type(coop_healpix_maps) masknew, hgs, maskcopy
    COOP_INT i, j, nsteps, nmaps, istart, iend
    COOP_REAL smoothscale, decay
    logical, optional::pol
    nsteps = ceiling(smoothscale/sqrt(coop_pi*4./this%npix))
    if(nsteps .le. 0 .or. nsteps .gt. 200)stop "coop_healpix_smooth_mask: invalid input of smoothscale"
    decay = 1.d0 !exp(-nefolds/nsteps/2.)
    maskcopy = mask
    call maskcopy%convert2nested()
    call this%convert2nested()
    maskcopy%mask_npix = count(maskcopy%map(:,1) .lt. 0.5)
    allocate(maskcopy%mask_listpix(maskcopy%mask_npix))
    j = 0
    do i = 0, this%npix - 1
       if(maskcopy%map(i,1).lt. 0.5)then
          j = j + 1
          maskcopy%mask_listpix(j) = i
       endif
    enddo
    hgs = this
    masknew = maskcopy
    istart = 1
    iend = 1
    if(present(pol))then
       if(pol .and. this%iq .gt. 0)then
          istart = this%iq
          iend = this%iu
       endif
    endif
    do i=1, nsteps
       call coop_healpix_iterate_mask(this, hgs, maskcopy, masknew, decay)
       call coop_healpix_iterate_mask(hgs, this, maskcopy, maskcopy, decay)
    enddo
    call hgs%free
    call masknew%free
    call maskcopy%free
    call this%convert2ring()
    
  contains 
    
    subroutine coop_healpix_iterate_mask(this_from, this_to, mask_from, mask_to, decay)  
      type(coop_healpix_maps) this_from, this_to
      type(coop_healpix_maps)  mask_from , mask_to
      COOP_INT list(8), nneigh
      COOP_INT i, j
      COOP_REAL decay, summask
#ifdef HAS_HEALPIX
      !$omp parallel do private(list, nneigh, i, summask, j)
      do i = 1, mask_from%mask_npix
         call neighbours_nest(mask_from%nside, mask_from%mask_listpix(i), list, nneigh)
         summask = sum(mask_from%map(list(1:nneigh), 1))
         if(mask_from%map(mask_from%mask_listpix(i),1) .lt. 0.5 .and. summask .gt. 0.)then
            do j = istart, iend
               this_to%map(mask_from%mask_listpix(i), j) =  sum(mask_from%map(list(1:nneigh), 1)*this_from%map(list(1:nneigh), j))/summask*decay
            enddo
            mask_to%map(mask_from%mask_listpix(i),1) = 1.
         endif
      enddo
      !$omp end parallel do
#else
      stop "CANNOT FIND HEALPIX"
#endif
    end subroutine coop_healpix_iterate_mask

  end subroutine coop_healpix_diffuse_into_mask

  subroutine coop_healpix_patch_plot(this, imap, output, use_degree, label)
    COOP_INT::bgrids
    class(coop_healpix_patch)::this
    COOP_INT imap
    COOP_UNKNOWN_STRING::output
    type(coop_asy)::fig
    COOP_INT nb, i, j, k, ns
    COOP_REAL  xc, yc,  norm, r, theta, minz, maxz
    COOP_REAL,dimension(:),allocatable::xstart, xend, ystart, yend
    logical,optional::use_degree
    logical use_rad
    COOP_SHORT_STRING::xlabel, ylabel
    COOP_UNKNOWN_STRING,optional::label    
    COOP_STRING::zlabel
    if(present(label))then
       zlabel = trim(adjustl(label))
    else
       zlabel = trim(adjustl(this%label(imap)))
    endif
    call fig%open(output)
    if(present(use_degree))then
       if(use_degree)then
          this%r = this%r/coop_SI_degree
       endif
       use_rad = .not. use_degree
    else
       use_rad = .true.
    endif
    if(use_rad)then
       xlabel = "$\varpi\cos\phi$"
       ylabel = "$\varpi\sin\phi$"       
    else
       xlabel = "$\varpi\cos\phi (\mathrm{deg})$"
       ylabel = "$\varpi\sin\phi (\mathrm{deg})$"
    endif
    call fig%init(caption = trim(this%caption), xlabel =trim(xlabel), ylabel =trim(ylabel), width = 5., height = 3.9, xmin = -real(this%r(this%n)), xmax = real(this%r(this%n)), ymin = -real(this%r(this%n)), ymax = real(this%r(this%n)))              
    if(imap .le. 0 .or. imap .gt. this%nmaps) stop "coop_healpix_patch_plot: imap overflow"
    if(this%zmin .lt.0.99e30)then
       minz = this%zmin
    else
       call coop_array_get_threshold(this%image(:,:,imap), COOP_REAL_OF(0.99), minz)
    endif
    if(this%zmax .gt. -0.99e30)then
       maxz = this%zmax
    else
       call coop_array_get_threshold(this%image(:,:,imap), COOP_REAL_OF(0.01), maxz)
    endif
    call coop_asy_density(fig, this%image(:,:,imap), -this%r(this%n), this%r(this%n), -this%r(this%n), this%r(this%n), label = trim(zlabel), zmax = maxz, zmin = minz, color_table = trim(this%color_table))
    if(use_rad)then
       theta = nint(2.d0*asin(this%r(this%n)/2.d0)/coop_SI_degree*10.d0)/10.d0
       call coop_asy_label(fig, "$\mathbf{-"//COOP_STR_OF(theta)//"}^\circ$", -this%r(this%n), -this%r(this%n)*1.15, color="blue")
       call coop_asy_label(fig, "$\mathbf{"//COOP_STR_OF(theta)//"}^\circ$", this%r(this%n), -this%r(this%n)*1.15, color="blue")
       call fig%arrow(this%r(this%n),  -this%r(this%n)*1.08, this%r(this%n),  -this%r(this%n)*1.01)
       call fig%arrow(-this%r(this%n),  -this%r(this%n)*1.08, -this%r(this%n),  -this%r(this%n)*1.01)
    endif
    if(this%headless_vectors .and. this%nmaps .eq. 2)then
       bgrids = max(this%n/7, 2)
       norm = maxval(this%image(:,:,1)**2+this%image(:,:,2)**2)
       if(norm .gt. 0.d0)then
          norm = bgrids*this%dr/2./sqrt(norm)*0.96
       else
          goto 100
       endif
       ns = floor((this%n-0.5d0*bgrids)/bgrids)
       nb = (2*ns+1)**2
       allocate(xstart(nb),  ystart(nb), xend(nb), yend(nb))
       k = 0
       ns = ns*bgrids
       select case(this%genre)
       case("QU")
          do j = -ns, ns, bgrids
             do i = -ns, ns, bgrids
                xc = i*this%dr
                yc = j*this%dr
                r = sqrt(this%image(i,j,1)**2+this%image(i,j,2)**2)*norm
                theta = 0.5d0*COOP_POLAR_ANGLE(this%image(i,j,1), this%image(i,j,2))
                k = k + 1
                xstart(k) = xc - r*cos(theta)
                ystart(k) = yc - r*sin(theta)
                xend(k) = 2*xc - xstart(k)
                yend(k) = 2*yc - ystart(k)
             enddo
          enddo
       case("QrUr")
          do j = -ns, ns, bgrids
             do i = -ns, ns, bgrids
                xc = i*this%dr
                yc = j*this%dr
                r = sqrt(this%image(i,j,1)**2+this%image(i,j,2)**2)*norm
                theta = 0.5d0*COOP_POLAR_ANGLE(this%image(i,j,1), this%image(i,j,2)) + COOP_POLAR_ANGLE(xc, yc)
                k = k + 1
                xstart(k) = xc - r*cos(theta)
                ystart(k) = yc - r*sin(theta)
                xend(k) = 2*xc - xstart(k)
                yend(k) = 2*yc - ystart(k)
             enddo
          enddo
       case default
          write(*,"(A)") trim(this%genre)
          stop "Unknown genre for headless vectors"
       end select
       call coop_asy_lines(fig, xstart, ystart, xend, yend, "black", "solid", 2.)
       deallocate(xstart, xend, ystart, yend)
    endif
100 call fig%close()
    if(present(use_degree))then
       if(use_degree)then
          this%r = this%r*coop_SI_degree
       endif
    endif    
  end subroutine coop_healpix_patch_plot



  subroutine coop_healpix_patch_plot_fft(this, imap, output,label)
    COOP_INT,parameter::lmax=100
    COOP_INT::bgrids
    class(coop_healpix_patch)::this
    COOP_INT imap
    COOP_UNKNOWN_STRING::output
    type(coop_asy)::fig
    COOP_INT i,j,ns
    COOP_REAL   minz, maxz, dk
    COOP_COMPLEX, dimension(:,:),allocatable:: fftmap
    COOP_REAL,dimension(:,:),allocatable::remap,immap
    COOP_SINGLE,parameter::width=6.5, height=5.
    COOP_SHORT_STRING::xlabel, ylabel
    COOP_UNKNOWN_STRING,optional::label    
    
    if(imap .le. 0 .or. imap .gt. this%nmaps) stop "coop_healpix_patch_plot: imap overflow"

    allocate(fftmap(0:this%n, 0:this%n*2))
    
    call coop_fft_forward(this%n*2+1, this%n*2+1, this%image(:,:,imap), fftmap)
    fftmap = fftmap/(2*this%n+1)**2

    dk = coop_2pi/(this%dr*(2*this%n+1))
    ns = min(nint(lmax/dk), nint(this%n/coop_sqrt2))

    allocate(remap(-ns:ns,-ns:ns), immap(-ns:ns,-ns:ns))    
    remap(0:ns, 0:ns) = real(fftmap(0:ns, 0:ns) )
    remap(0:ns, -ns:-1) = real(fftmap(0:ns, this%n*2+1-ns:this%n*2))
    
    immap(0:ns, 0:ns) = aimag(fftmap(0:ns, 0:ns) )
    immap(0:ns, -ns:-1) = aimag(fftmap(0:ns, this%n*2+1-ns:this%n*2))

    remap(-ns:-1,-ns:ns)=remap(ns:1:-1,ns:-ns:-1)
    immap(-ns:-1,-ns:ns)=-immap(ns:1:-1,ns:-ns:-1)

    
    call fig%open(coop_file_add_postfix(output, "_FFTRe"))
    
    xlabel = "$\ell_x$"
    ylabel = "$\ell_y$"
    call fig%init(caption = trim(this%caption)//", Re($f_\ell$)", xlabel =trim(xlabel), ylabel =trim(ylabel), width = width, height = height, xmin = -real(ns*dk), xmax = real(ns*dk), ymin = -real(ns*dk), ymax = real(ns*dk))
    call coop_array_get_threshold(remap, COOP_REAL_OF(0.995), minz)
    call coop_array_get_threshold(remap, COOP_REAL_OF(0.005), maxz)
    if(present(label))then
       call coop_asy_density(fig, remap, -ns*dk, ns*dk,-ns*dk, ns*dk, label = "$\mathrm{Re} "//trim(adjustl(label))//"_\ell$", zmax = maxz, zmin = minz, color_table = trim(this%color_table))
    else
       call coop_asy_density(fig, remap, -ns*dk, ns*dk,-ns*dk, ns*dk, label = trim(adjustl(this%label(imap))), zmax = maxz, zmin = minz, color_table = trim(this%color_table))       
    endif
    call fig%close()

    call fig%open(coop_file_add_postfix(output, "_FFTIm"))    
    xlabel = "$\ell_x$"
    ylabel = "$\ell_y$"
    call fig%init(caption = trim(this%caption)//", Im($f_\ell$)", xlabel =trim(xlabel), ylabel =trim(ylabel), width = width, height = height, xmin = -real(ns*dk), xmax = real(ns*dk), ymin = -real(ns*dk), ymax = real(ns*dk))
    call coop_array_get_threshold(immap, COOP_REAL_OF(0.995), minz)
    call coop_array_get_threshold(immap, COOP_REAL_OF(0.005), maxz)
    if(present(label))then
       call coop_asy_density(fig, immap, -ns*dk, ns*dk,-ns*dk, ns*dk, zmax = maxz, zmin = minz, label = "$\mathrm{Im}"//trim(adjustl(label))//"_\ell$", color_table = trim(this%color_table))
    else
       call coop_asy_density(fig, immap, -ns*dk, ns*dk,-ns*dk, ns*dk, zmax = maxz, zmin = minz, label =trim(adjustl(this%label(imap))), color_table = trim(this%color_table))       
    endif
    call fig%close()

    call fig%open(coop_file_add_postfix(output, "_FFTl2Cl"))
    immap = (remap**2+immap**2)*dk**2
    do i=-ns,ns
       do j=-ns,ns
          immap(i,j)=immap(i,j)*dble(i**2+j**2)
       enddo
    enddo
    xlabel = "$\ell_x$"
    ylabel = "$\ell_y$"
    call fig%init(caption = trim(this%caption)//", $\ell^2|f_\ell|^2$", xlabel =trim(xlabel), ylabel =trim(ylabel), width = width, height = height, xmin = -real(ns*dk), xmax = real(ns*dk), ymin = -real(ns*dk), ymax = real(ns*dk))
    call coop_array_get_threshold(immap, COOP_REAL_OF(0.995), minz)
    call coop_array_get_threshold(immap, COOP_REAL_OF(0.005), maxz)
    if(present(label))then
       call coop_asy_density(fig, immap, -ns*dk, ns*dk, -ns*dk, ns*dk, zmax = maxz, zmin = minz, label = "$\ell^2 |"//trim(adjustl(label))//"_\ell|^2 $", color_table = trim(this%color_table))       
    else
       call coop_asy_density(fig, immap, -ns*dk, ns*dk, -ns*dk, ns*dk, zmax = maxz, zmin = minz, label = trim(this%label(imap)), color_table = trim(this%color_table))
    endif
    call fig%close()    
    
    deallocate(fftmap,remap,immap)
  end subroutine coop_healpix_patch_plot_fft
 

  subroutine coop_healpix_patch_free(this)
    class(coop_healpix_patch) this
    if(allocated(this%image))deallocate(this%image)
    if(allocated(this%label))deallocate(this%label)
    if(allocated(this%r))deallocate(this%r)
    if(allocated(this%fr))deallocate(this%fr)
    if(allocated(this%wcm))deallocate(this%wcm)
    if(allocated(this%wsm))deallocate(this%wsm)
    if(allocated(this%icm))deallocate(this%icm)
    if(allocated(this%nstack))deallocate(this%nstack)
    if(allocated(this%indisk))deallocate(this%indisk)
    this%n = -1
    this%mmax = -1
  end subroutine coop_healpix_patch_free

  subroutine coop_healpix_patch_init(this, genre, n, dr, mmax)
    class(coop_healpix_patch) this
    COOP_UNKNOWN_STRING::genre
    COOP_INT n
    COOP_REAL dr, cosmt, sinmt, theta
    COOP_INT i,j,m
    COOP_INT, optional::mmax
    COOP_REAL sumr(0:n+1)
    call this%free()
    this%caption = ""
    this%headless_vectors = .true.
    this%color_table = "Rainbow"
    this%zmin = 1.e30
    this%zmax = -1.e30
    this%genre = trim(adjustl(genre))
    this%n = n
    this%npix = (2*this%n+1)**2
    this%dr = dr
    if(present(mmax))then
       this%mmax = mmax
    else
       this%mmax = 4
    endif
    if(this%n .lt. 0) return
    select case(trim(this%genre))
    case("QU", "QrUr")
       this%nmaps = 2
    case("T","E","B", "I", "zeta")
       this%nmaps = 1
    case default
       write(*,*) "Unknown stacking genre: "//trim(this%genre)
       write(*,*) "Only supports: QU, QrUr, T, E, B"
       stop
    end select
    allocate(this%label(this%nmaps))
    do i=1, this%nmaps
       this%label(i) = ""
    enddo
    select case(trim(this%genre))
    case("QU")
       this%label(1) = "$Q(\mu K)$"
       this%label(2) = "$U(\mu K)$"
    case("QrUr")
       this%label(1) = "$Q_r(\mu K)$"
       this%label(2) = "$U_r(\mu K)$"
    case("T","E","B", "I")
       this%label(1) = "$"//trim(this%genre)//"(\mu K)$"
    case("zeta")
       this%label(1) = "$\zeta (10^{-5})$"
    case default
       write(*,*) "Unknown stacking genre: "//trim(this%genre)
       write(*,*) "Only supports: QU, QrUr, T, E, B"
       stop
    end select

    allocate(this%image(-this%n:this%n, -this%n:this%n, this%nmaps))
    allocate(this%nstack(-this%n:this%n, -this%n:this%n))
    allocate(this%indisk(-this%n:this%n, -this%n:this%n))
    allocate(this%r(0:this%n))
    allocate(this%fr(0:this%n, 0:this%mmax/2, this%nmaps))
    allocate(this%wcm(-this%n:this%n, -this%n:this%n, 0:this%mmax+1))
    allocate(this%wsm(-this%n:this%n, -this%n:this%n, 2:this%mmax+1))
    allocate(this%icm(-this%n:this%n, -this%n:this%n, 0:1))
    this%image = 0.
    this%wcm = 0.
    this%wsm = 0.
    this%fr = 0.

    this%indisk = 1.d0
    do j=1, this%n
       i = ceiling(sqrt((this%n-j)*dble(this%n+j)+1.d-20))
       this%indisk(i:this%n, j) = 0.d0
       this%indisk(-this%n:-i, j) = 0.d0
       this%indisk(i:this%n, -j) = 0.d0
       this%indisk(-this%n:-i, -j) = 0.d0
    enddo
    !$omp parallel do
    do i=0, this%n
       this%r(i) = this%dr * i
    enddo
    !$omp end parallel do


    this%num_indisk_tol = count(this%indisk .ne. 0.d0)*coop_healpix_mask_tol

    !$omp parallel do private(i, j)
    do j=-this%n, this%n
       do i=-this%n, this%n
          this%wcm(i, j, 0) = sqrt(dble(i)**2+dble(j)**2)
          this%icm(i, j, 0) = floor(this%wcm(i, j, 0))
          this%icm(i, j, 1) = this%icm(i, j, 0) + 1
          this%wcm(i, j, 1) = this%wcm(i, j, 0) - this%icm(i, j, 0)
          this%wcm(i, j, 0) = 1.d0 - this%wcm(i, j, 1)
       enddo
    enddo
    !$omp end parallel do
    sumr = 0.
    do j=-this%n, this%n
       do i=-this%n, this%n
          if(this%icm(i,j,0).le. this%n)then
             sumr(this%icm(i,j,0)) = sumr(this%icm(i,j,0)) + this%wcm(i,j,0)
             sumr(this%icm(i,j,1)) = sumr(this%icm(i,j,1)) + this%wcm(i,j,1)
          endif
       enddo
    enddo
    !$omp parallel do private(i, j)
    do j=-this%n, this%n
       do i=-this%n, this%n
          if(this%icm(i, j, 0) .le. this%n)then
             this%wcm(i, j, 0) = this%wcm(i, j, 0)/sumr(this%icm(i, j, 0))
          else
             this%icm(i, j, 0) = 0
             this%wcm(i, j, 0) = 0.d0
          endif
          if(this%icm(i, j, 1) .le. this%n)then
             this%wcm(i, j, 1) = this%wcm(i, j, 1)/sumr(this%icm(i, j, 1))
          else
             this%icm(i, j, 1) = 0
             this%wcm(i, j, 1) = 0.d0
          endif
       enddo
    enddo
    !$omp end parallel do

    !$omp parallel do private(m, i, j, cosmt, sinmt, theta)
    do m = 2, this%mmax, 2
       do j=-this%n, this%n
          do i = -this%n, this%n
             if(this%icm(i,j,0).ne.0)then
                theta = atan2(dble(j), dble(i))
                cosmt = cos(m*theta)*2.d0
                sinmt = sin(m*theta)*2.d0
                this%wcm(i, j, m) = this%wcm(i, j, 0)*cosmt
                this%wcm(i, j, m+1) = this%wcm(i, j, 1)*cosmt
                this%wsm(i, j, m) = this%wcm(i, j, 0)*sinmt
                this%wsm(i, j, m+1) = this%wcm(i, j, 1)*sinmt
             else
                this%wcm(i, j, m) = 0.d0
                this%wcm(i, j, m+1) = 0.d0
                this%wsm(i, j, m) = 0.d0
                this%wsm(i, j, m+1) = 0.d0
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do
  end subroutine coop_healpix_patch_init


  subroutine coop_healpix_patch_get_all_radial_profiles(this)
    class(coop_healpix_patch)::this
    COOP_INT i, j, imap, m
    if(this%mmax .lt. 0) return
    do imap = 1, this%nmaps
       do m = 0, this%mmax, 2
          call this%get_radial_profile(imap, m)
       enddo
    enddo
  end subroutine coop_healpix_patch_get_all_radial_profiles

  subroutine coop_healpix_patch_get_radial_profile(this, imap, m)
    class(coop_healpix_patch)::this
    COOP_INT m, imap, halfm
    COOP_INT i,j
    if(m.gt. this%mmax .or. mod(m,2).ne.0 .or. imap.gt.this%nmaps .or. imap.le.0) stop "coop_healpix_patch_get_radial_profile: wrong input arguments"
    halfm = m/2
    this%fr(:, halfm, imap) = 0.d0
    select case(imap)
    case(1)
       do i=-this%n, this%n
          do j=-this%n, this%n
             this%fr(this%icm(i, j, 0), halfm, imap) =  this%fr(this%icm(i, j, 0), halfm, imap) + this%image(i, j, imap) * this%wcm(i, j, m)
             this%fr(this%icm(i, j, 1), halfm, imap) =  this%fr(this%icm(i, j, 1), halfm, imap) + this%image(i, j, imap)* this%wcm(i, j, m+1) 
          enddo
       enddo
    case(2)
       do i=-this%n, this%n
          do j=-this%n, this%n
             this%fr(this%icm(i, j, 0), halfm, imap) =  this%fr(this%icm(i, j, 0), halfm, imap) + this%image(i, j, imap) * this%wsm(i, j, m)
             this%fr(this%icm(i, j, 1), halfm, imap) =  this%fr(this%icm(i, j, 1), halfm, imap) + this%image(i, j, imap)* this%wsm(i, j, m+1) 
          enddo
       enddo
    case default
       stop "Cannot get radial profile for more than 2 maps."
    end select
  end subroutine coop_healpix_patch_get_radial_profile


  subroutine coop_healpix_maps_simulate(this)
    class(coop_healpix_maps) this
    real,dimension(:),allocatable::sqrtCls
    real,dimension(:, :),allocatable::Cls_sqrteig
    real,dimension(:,:,:),allocatable::Cls_rot
    COOP_INT l
    if(this%nmaps.eq.1 .and. this%spin(1).eq.0)then
       allocate(sqrtCls(0:this%lmax))
       !$omp parallel do
       do l = 0, this%lmax
          sqrtCls(l) = sqrt(this%Cl(l,1))
       enddo
       !$omp end parallel do
       call coop_healpix_maps_simulate_Tmaps(this, this%nside, this%lmax, sqrtCls)
       deallocate(sqrtCls)
    elseif(this%nmaps.eq.3 .and. this%iq .eq.2)then
       allocate(Cls_sqrteig(3, 0:this%lmax), Cls_rot(3,3,0:this%lmax))
       call coop_healpix_Cls2Rot(this%lmax, this%Cl, Cls_sqrteig, Cls_rot)
       call coop_healpix_maps_simulate_TQUmaps(this, this%nside, this%lmax, Cls_sqrteig, Cls_rot)
       deallocate(Cls_sqrteig, Cls_rot)
    else
       stop "unknown coop_healpix_maps_simulate mode"
    endif
  end subroutine coop_healpix_maps_simulate


  subroutine coop_healpix_maps_simulate_Tmaps(this, nside, lmax, sqrtCls)
    class(coop_healpix_maps) this
    COOP_INT nside
    COOP_INT lmax
    COOP_SINGLE sqrtCls(0:lmax)
    COOP_INT l,m
    call this%init( nside = nside, nmaps = 1, spin = (/ 0 /), lmax = lmax)
    !$omp parallel do private(l, m)
    do l=0, lmax
       this%alm(l, 0, 1) = coop_random_complex_Gaussian(.true.)*SqrtCls(l)     
       do m = 1, l
          this%alm(l, m, 1) = coop_random_complex_Gaussian()*SqrtCls(l)
       enddo
    enddo
    !$omp end parallel do
    call this%alm2map
  end subroutine coop_healpix_maps_simulate_Tmaps


  subroutine coop_healpix_maps_get_cls(this) !!I assume you have already called    this_map2alm(this)
    class(coop_healpix_maps)this
    COOP_INT l, m, i, j, k
    if(.not.allocated(this%alm)) stop "coop_healpix_maps_get_cls: you have to call coop_healpix_maps_map2alm before calling this subroutine"
    !$omp parallel do private(i,j,k,l)
    do i=1, this%nmaps
       do j=1, i
          k = coop_matsym_index(this%nmaps, i, j)
          do l = 0, this%lmax
             this%Cl(l, k) = (sum(COOP_MULT_REAL(this%alm(l, 1:l, i), this%alm(l, 1:l, j))) + 0.5d0 * COOP_MULT_REAL(this%alm(l,0,i), this%alm(l,0,j)) )/(l+0.5d0)
          enddo
       enddo
    enddo
    !$omp end parallel do
  end subroutine coop_healpix_maps_get_cls

  subroutine coop_healpix_Cls2Rot(lmax, Cls, Cls_sqrteig, Cls_rot)
    COOP_INT lmax
    real,dimension(0:lmax, 6),intent(IN)::Cls !!ordering is TT, EE, BB, TE, EB, TB
    real, dimension(3, 0:lmax),intent(OUT)::Cls_sqrteig
    real, dimension(3, 3, 0:lmax),intent(OUT)::Cls_rot
    COOP_INT l
    COOP_REAL a2(2,2), a3(3,3)
    COOP_REAL psi2(2,2), psi3(3,3)
    if(all(Cls(:,coop_healpix_index_EB).eq.0.d0) .and. all(Cls(:, coop_healpix_index_TB).eq.0.d0))then
       Cls_sqrteig(coop_healpix_index_BB,:) = sqrt(Cls(:, coop_healpix_index_BB))
       Cls_rot(coop_healpix_index_BB, coop_healpix_index_TT, :) = 0
       Cls_rot(coop_healpix_index_BB, coop_healpix_index_EE, :) = 0
       Cls_rot(coop_healpix_index_TT, coop_healpix_index_BB, :) = 0
       Cls_rot(coop_healpix_index_EE, coop_healpix_index_BB, :) = 0
       Cls_rot(coop_healpix_index_BB, coop_healpix_index_BB, :) = 1.
       do l=0, lmax
          a2(1,1) = Cls(l, coop_healpix_index_TT)
          a2(2,2) = Cls(l, coop_healpix_index_EE)
          a2(1,2) = Cls(l, coop_healpix_index_TE)
          a2(2,1) = a2(1,2)
          call coop_matsymdiag_small(2, a2, psi2)
          Cls_sqrteig(coop_healpix_index_TT, l) = sqrt(a2(1,1))
          Cls_sqrteig(coop_healpix_index_EE, l) = sqrt(a2(2,2))
          Cls_rot( coop_healpix_index_TT:coop_healpix_index_EE, coop_healpix_index_TT:coop_healpix_index_EE, l) = psi2
       end do
    else
       do l=0, lmax
          a3(1,1) = Cls(l,coop_healpix_index_TT)
          a3(2,2) = Cls(l,coop_healpix_index_EE)
          a3(3,3) = Cls(l,coop_healpix_index_BB)
          a3(1,2) = Cls(l,coop_healpix_index_TE)
          a3(2,3) = Cls(l,coop_healpix_index_EB)
          a3(3,1) = Cls(l,coop_healpix_index_TB)
          a3(2,1) = a3(1,2)
          a3(3,2) = a3(2,3)
          a3(1,3) = a3(3,1)
          call coop_matsymdiag_small(3, a3, psi3)
          Cls_sqrteig(coop_healpix_index_TT, l) = sqrt(a3(1,1))
          Cls_sqrteig(coop_healpix_index_EE, l) = sqrt(a3(2,2))
          Cls_sqrteig(coop_healpix_index_BB, l) = sqrt(a3(3,3))
          Cls_rot(:, :, l) = psi3
       end do
    endif
  end subroutine coop_healpix_Cls2Rot

  subroutine coop_healpix_maps_iqu2TQTUT(this)
    class(coop_healpix_maps) this
    call this%map2alm(index_list = (/ 1 /) )
    this%alm(:,:,2) = this%alm(:,:,1)
    this%alm(:,:,3) = 0
    this%spin(2:3) = 2
    this%spin(1) = 0
    call this%alm2map()
  end subroutine coop_healpix_maps_iqu2TQTUT

  subroutine coop_healpix_maps_iqu2TEB(this)
    class(coop_healpix_maps) this
    call this%map2alm( index_list = (/ 2 ,3 /) )
    this%spin = 0
    this%iq = 0
    this%iu = 0
    call this%alm2map(index_list = (/ 2, 3 /) )
  end subroutine coop_healpix_maps_iqu2TEB

  subroutine coop_healpix_maps_qu2EB(this)
    class(coop_healpix_maps) this
    this%spin = 2
    this%iq = 1
    this%iu = 2
    call this%map2alm(index_list = (/ 1, 2 /) )
    this%spin = 0
    this%iq = 0
    this%iu = 0
    call this%alm2map(index_list = (/ 1, 2 /) )
  end subroutine coop_healpix_maps_qu2EB



  subroutine coop_healpix_maps_iqu2LapTEB(this)
    class(coop_healpix_maps) this
    COOP_INT l
    call this%map2alm( index_list = (/ 1, 2, 3 /) )
    this%spin = 0
    this%iq = 0
    this%iu = 0
    do l=0, this%lmax
       this%alm(l, :, :) = this%alm(l, :, :)*l*(l+1.d0)
    enddo
    call this%alm2map( index_list = (/ 1, 2, 3 /) )
  end subroutine coop_healpix_maps_iqu2LapTEB


  subroutine coop_healpix_maps_teb2iqu(this)
    class(coop_healpix_maps) this
    call this%map2alm( index_list = (/ 2, 3 /) )
    this%spin = (/ 0, 2, 2 /)
    this%iq = 2
    this%iu = 3
    call this%alm2map( index_list = (/ 2, 3 /) )
  end subroutine coop_healpix_maps_teb2iqu


  subroutine coop_healpix_maps_simulate_TQUmaps(this, nside, lmax, Cls_sqrteig, Cls_rot)
    class(coop_healpix_maps) this
    COOP_INT lmax, nside
    real,dimension(3, 0:lmax)::Cls_sqrteig
    real,dimension(3, 3, 0:lmax)::Cls_rot
    COOP_INT l, m
    call this%init(nside = nside, nmaps = 3, spin = (/ 0, 2, 2 /), lmax = lmax)    
    !$omp parallel do private(l, m)
    do l=0, lmax
       this%alm(l, 0, :) = matmul(Cls_rot(:, :, l), Cls_sqrteig(:,l) * (/ coop_random_complex_Gaussian(.true.), coop_random_complex_Gaussian(.true.), coop_random_complex_Gaussian(.true.) /) )
       do m = 1, l
          this%alm(l, m, :) = matmul(Cls_rot(:, :, l), Cls_sqrteig(:,l) * (/ coop_random_complex_Gaussian(), coop_random_complex_Gaussian(), coop_random_complex_Gaussian() /) )
       enddo
    enddo
    !$omp end parallel do
    call this%alm2map()
  end subroutine coop_healpix_maps_simulate_TQUmaps


  subroutine coop_healpix_maps_init(this, nside, nmaps, spin, lmax)
    class(coop_healpix_maps) this
    COOP_INT:: nside, nmaps
    COOP_INT:: spin(nmaps)
    COOP_INT, optional::lmax
#ifdef HAS_HEALPIX
    if(allocated(this%map))then
       if(this%nside .eq. nside .and. this%nmaps.eq.nmaps)then
          goto 100
       endif
       deallocate(this%map)
    endif
    if(allocated(this%spin))deallocate(this%spin)
    this%nside = nside
    this%nmaps = nmaps
    this%npix = nside2npix(nside)
    allocate(this%map(0:this%npix - 1, nmaps))
    allocate(this%spin(this%nmaps))
100 this%spin = spin
    if(all(this%spin(1:this%nmaps-1) .ne. 2))then
       this%iq = 0
       this%iu = 0
    else
       this%iq = 1
       do while(this%spin(this%iq) .ne. 2)
          this%iq = this%iq + 1
       enddo
       this%iu = this%iq + 1
       if(this%spin(this%iu).ne.2)then
          this%iq = 0
          this%iu = 0
       endif
    endif
    if(present(lmax)) call this%allocate_alms(lmax)
200 this%ordering = COOP_RING !!default ordering
    call write_minimal_header(this%header,dtype = 'MAP', nside=this%nside, order = this%ordering, creator='Zhiqi Huang', version = 'COOP', units='muK', polar=any(this%spin.eq.2) )
    this%maskpol_npix = 0
#else
    stop "DID not find healpix"
#endif
  end subroutine coop_healpix_maps_init


  subroutine coop_healpix_maps_extend(this, nmaps)
    class(coop_healpix_maps) this
    COOP_INT nmaps
    COOP_SINGLE,dimension(:,:),allocatable::map
    COOP_INT, dimension(:), allocatable::spin
#ifdef HAS_HEALPIX
    if(nmaps .le. this%nmaps) return
    allocate(map(0:this%npix-1, this%nmaps))
    map = this%map
    allocate(spin(this%nmaps))
    spin  = this%spin
    deallocate(this%map, this%spin)
    allocate(this%map(0:this%npix-1, nmaps), this%spin(nmaps))
    this%map(:, 1:this%nmaps) = map
    this%spin(1:this%nmaps) = spin
    if(allocated(this%alm))then
       deallocate(this%alm)
       call this%allocate_alms(this%lmax)
    endif
    deallocate(map, spin)
#else
    stop "DID not find healpix"
#endif
  end subroutine coop_healpix_maps_extend


  subroutine coop_healpix_maps_allocate_alms(this, lmax)
    class(coop_healpix_maps) this
    COOP_INT lmax
    if(allocated(this%alm))then
       if(this%lmax  .eq. lmax )then
          return
       endif
       deallocate(this%alm, this%alm_done, this%checksum)
    endif
    if(allocated(this%cl))deallocate(this%cl)
    this%lmax = lmax
    allocate(this%alm(0:this%lmax, 0:this%lmax, this%nmaps))
    allocate(this%cl(0:this%lmax, this%nmaps*(this%nmaps+1)/2))
    allocate(this%alm_done(this%nmaps), this%checksum(4, this%nmaps))
    this%alm_done = .false.
    this%checksum = 1.e30
    this%alm = 0.
  end subroutine coop_healpix_maps_allocate_alms

  subroutine coop_healpix_maps_free(this)
    class(coop_healpix_maps) this
    if(allocated(this%map))deallocate(this%map)
    if(allocated(this%alm))deallocate(this%alm)
    if(allocated(this%cl))deallocate(this%cl)
    if(allocated(this%spin))deallocate(this%spin)
    if(allocated(this%mask_listpix))deallocate(this%mask_listpix)
    if(allocated(this%maskpol_listpix))deallocate(this%maskpol_listpix)
    if(allocated(this%alm_done))deallocate(this%alm_done)
    if(allocated(this%checksum))deallocate(this%checksum)
  end subroutine coop_healpix_maps_free

  subroutine coop_healpix_maps_read(this, filename, nmaps_wanted, spin, nmaps_to_read, known_size)
    class(coop_healpix_maps) this
    COOP_UNKNOWN_STRING filename
    COOP_INT,optional::nmaps_wanted, nmaps_to_read
    COOP_INT,dimension(:),optional::spin
    integer(8) npixtot
    COOP_INT nmaps_actual
    logical,optional::known_size
#ifdef HAS_HEALPIX
    if(.not. coop_file_exists(filename))then
       write(*,*) trim(filename)
       stop "cannot find the file"
    endif
    if(present(known_size))then
       if(known_size)then
          nmaps_actual = this%nmaps
          goto 200
       endif
    endif
    npixtot = getsize_fits(trim(filename), nmaps = nmaps_actual, nside = this%nside, ordering = this%ordering)
    this%npix =nside2npix(this%nside)
    if(present(nmaps_wanted))then       
       this%nmaps = nmaps_wanted
       if(nmaps_wanted .lt. nmaps_actual)then          
          nmaps_actual = nmaps_wanted
       endif
    else
       this%nmaps = nmaps_actual
    endif
    if(allocated(this%spin))then
       if(size(this%spin) .ne. this%nmaps)then
          deallocate(this%spin)
          allocate(this%spin(this%nmaps))
       endif
    else
       allocate(this%spin(this%nmaps))
    endif
    if(present(spin))then
       if(size(spin).ne. this%nmaps)then
          stop "coop_healpix_maps_read: the list of spins should have the same size as nmaps"
       else
          this%spin = spin
          this%iq = 1
          do while(this%spin(this%iq) .ne.2 .and. this%iq .lt. this%nmaps)
             this%iq = this%iq + 1
          enddo
          if(this%iq .ge. this%nmaps)then
             this%iq = 0
             this%iu = 0
          else
             this%iu = this%iq + 1
          endif
       endif
    else
       select case(this%nmaps)
       case(3)
          if(index(filename, "TEB") .eq. 0 .and. index(filename, "teb") .eq. 0)then
             this%spin(1) = 0
             this%spin(2:3) = 2
             this%iq = 2
             this%iu = 3
             write(*,*) "I assume it is an IQU map, specify spins otherwise"
          else
             this%spin = 0
             this%iq = 0
             this%iu = 0
             write(*,*) "I assume all maps are scalar, specify spins otherwise"
          endif
       case(2)  
          if(index(filename, "TE") .eq. 0 .and. index(filename, "EB").eq.0)then
             this%spin(1:2) = 2
             this%iq = 1
             this%iu = 2
             write(*,*) "I assume it is an QU map, specify spins otherwise"
          else
             this%iq = 0
             this%iu = 0
             this%spin = 0
             write(*,*) "I assume all maps are scalar, specify spins otherwise"
          endif
       case default
          this%iq = 0
          this%iu = 0
          this%spin = 0
          write(*,*) "I assume all maps are scalar, specify spins otherwise"
       end select
    endif
    if(allocated(this%map))then
       if(size(this%map, 1).ne. this%npix .or. size(this%map, 2).ne.this%nmaps)then
          deallocate(this%map)
          allocate(this%map(0:this%npix-1, this%nmaps))
       endif
    else
       allocate(this%map(0:this%npix-1, this%nmaps))
    endif
200 if(present(nmaps_to_read))then
       call input_map(trim(filename), this%map, this%npix, min(nmaps_actual, nmaps_to_read), fmissval = 0.)
    else
       call input_map(trim(filename), this%map, this%npix, nmaps_actual, fmissval = 0.)
    endif
    call this%convert2ring
    call write_minimal_header(this%header,dtype = 'MAP', nside=this%nside, order = this%ordering, creator='Zhiqi Huang', version = 'CosmoLib', units='muK', polar=any(this%spin.eq.2) )
#else
    stop "DID NOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_maps_read



    subroutine coop_healpix_maps_import(this, filename, index_start, index_end, spin)
    class(coop_healpix_maps) this
    COOP_UNKNOWN_STRING filename
    type(coop_healpix_maps)::tmp
    COOP_INT index_start, index_end, spin(index_end-index_start + 1)
#ifdef HAS_HEALPIX
    call tmp%read(filename = trim(filename), nmaps_wanted = index_end-index_start + 1, spin = spin )
    call this%convert2ring()    
    if(tmp%nside .ne. this%nside) stop "import: nside must be the same"
    if(index_end .gt. this%nmaps)then
       call this%extend(index_end)
    endif
    this%map(:, index_start:index_end) = tmp%map
    this%spin(index_start:index_end) = spin
    call tmp%free()
#else
    stop "DID NOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_maps_import


  subroutine coop_healpix_maps_write(this, filename, index_list)
    class(coop_healpix_maps)this
    COOP_UNKNOWN_STRING filename
    COOP_INT,dimension(:),optional::index_list
    logical pol
    if(present(index_list))then
       if(any(index_list .lt. 1 .or. index_list .gt. this%nmaps)) stop "coop_healpix_write_map: index out of range"
       pol = any(this%spin(index_list).eq.2)
    else
       pol =any(this%spin.eq.2)
    endif
    call coop_delete_file(trim(filename))
    if(allocated(this%alm))then
#ifdef HAS_HEALPIX
       call write_minimal_header(this%header,dtype = 'MAP', nside=this%nside, order = this%ordering, creator='Zhiqi Huang', version = 'CosmoLib', units='muK', nlmax = this%lmax, nmmax = this%lmax, polar= pol)
#else
       stop "DID NOT FIND HEALPIX"
#endif
    else
#ifdef HAS_HEALPIX
       call write_minimal_header(this%header,dtype = 'MAP', nside=this%nside, order = this%ordering, creator='Zhiqi Huang', version = 'CosmoLib', units='muK', polar= pol )
#else
       stop "DID NOT FIND HEALPIX"
#endif
    endif
#ifdef HAS_HEALPIX
    if(present(index_list))then
       call output_map(this%map(:, index_list), this%header, trim(filename))
    else
       call output_map(this%map, this%header, trim(filename))
    endif
#endif
  end subroutine coop_healpix_maps_write

  subroutine coop_healpix_convert_to_nested(this)
    class(coop_healpix_maps) this
#ifdef HAS_HEALPIX
    if(.not. allocated(this%map)) stop "coop_healpix_convert_to_nested: map is not allocated yet"
    if(this%ordering .eq. COOP_RING) then
       call convert_ring2nest(this%nside, this%map)
       this%ordering = COOP_NESTED
    elseif(this%ordering .ne. COOP_NESTED)then
       write(*,*) "ordering = ", this%ordering
       stop "coop_healpix_convert_to_nested: unknown ordering"
    endif
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_convert_to_nested

  subroutine coop_healpix_convert_to_ring(this)
    class(coop_healpix_maps) this
#ifdef HAS_HEALPIX
    if(.not. allocated(this%map)) stop "coop_healpix_convert_to_ring: map is not allocated yet"
    if(this%ordering .eq. COOP_NESTED)then
       call convert_nest2ring(this%nside, this%map)
       this%ordering = COOP_RING
    elseif(this%ordering .ne. COOP_RING)then
       write(*,*) "ordering = ", this%ordering
       stop "coop_healpix_convert_to_ring: UNKNOWN ordering"
    endif
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_convert_to_ring

  subroutine coop_healpix_maps_map2alm(this, lmax, index_list)
    class(coop_healpix_maps) this
    COOP_INT,optional::lmax
    COOP_INT i, l, j, lm, n
    COOP_INT,dimension(:),optional::index_list
    complex, dimension(:,:,:),allocatable::alm
#ifdef HAS_HEALPIX
    call this%convert2ring()
    if(present(lmax))then
       if(lmax .gt. this%nside*3)then
          write(*,*) "lmax > nside x 3 is not recommended"
          stop
       else
          lm = lmax          
       endif
    else
       lm =  min(coop_healpix_default_lmax, this%nside*2)
    endif
    call this%allocate_alms(lm)
    if(present(index_list))then
       j = 1
       n = size(index_list)
       do while(j.le.n)
          i = index_list(j)
          if(i.gt.this%nmaps)then
             j  = j + 1
             cycle
          endif
          if(coop_healpix_alm_check_done)then
             if(this%alm_done(i))then
                if( this%checksum(1, i).eq. this%map(0, i) .and. this%checksum(2, i).eq. this%map(this%npix/2, i) .and. this%checksum(3, i) .eq. this%map(this%npix-1, i) .and. nint(this%checksum(4, i)) .eq. this%spin(i)  )then
                   j  = j + 1
                   cycle
                endif
             endif
             this%alm_done(i) = .true.
             this%checksum(1, i) = this%map(0, i)
             this%checksum(2, i) = this%map(this%npix/2, i)
             this%checksum(3, i) = this%map(this%npix-1, i)
             this%checksum(4, i) = this%spin(i)
          endif
          if(this%spin(i).eq.0)then
             call map2alm(this%nside, this%lmax, this%lmax, this%map(:,i), this%alm(:,:,i:i))
             j = j + 1
          else
             if(i.lt. this%nmaps .and. j.lt. n)then
                if(index_list(j+1).eq.i+1 .and. this%spin(i+1) .eq. this%spin(i))then
                   if(.not. allocated(alm))allocate(alm(2, 0:this%lmax, 0:this%lmax))
                   call map2alm_spin(this%nside, this%lmax, this%lmax, this%spin(i), this%map(:,i:i+1), alm)
                   this%alm(:,:,i) = alm(1, :, :)
                   this%alm(:,:,i+1) = alm(2, :, :)
                   j = j + 2
                   cycle
                endif
             endif
             write(*,*) this%spin
             stop "coop_healpix_maps_map2alm: nonzero spin maps must appear in pairs"
          endif
       enddo
    else
       i = 1
       do while(i.le. this%nmaps)
          if(this%spin(i).eq.0)then
             call map2alm(this%nside, this%lmax, this%lmax, this%map(:,i), this%alm(:,:,i:i))
             i = i + 1
          else
             if(.not. allocated(alm))allocate(alm(2, 0:this%lmax, 0:this%lmax))
             if(i.lt. this%nmaps)then
                if(this%spin(i+1) .eq. this%spin(i))then
                   call map2alm_spin(this%nside, this%lmax, this%lmax, this%spin(i), this%map(:,i:i+1), alm)
                   this%alm(:,:,i) = alm(1, :, :)
                   this%alm(:,:,i+1) = alm(2, :, :)
                   i = i + 2
                   cycle
                endif
             endif
             write(*,*) this%spin
             stop "coop_healpix_maps_map2alm: nonzero spin maps must appear in pairs"
          endif
       enddo
    endif
    if(allocated(alm))deallocate(alm)
    if(coop_healpix_want_cls)call this%get_cls
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_maps_map2alm


  subroutine coop_healpix_maps_alm2map(this, index_list)
    class(coop_healpix_maps) this
    COOP_INT i, j, n
    COOP_INT,dimension(:),optional::index_list
    complex,dimension(:,:,:),allocatable::alm
#ifdef HAS_HEALPIX
    if(present(index_list))then
       j = 1
       n = size(index_list)
       do while(j.le.n)
          i = index_list(j)
          if(i.gt. this%nmaps)then
             j = j + 1
             cycle
          endif
          if(this%spin(i).eq.0)then
             call alm2map(this%nside, this%lmax, this%lmax, this%alm(:,:,i:i), this%map(:,i))
             j = j + 1
             cycle
          endif
          if(j .lt. n .and. i .lt. this%nmaps)then
             if(index_list(j+1).eq.i+1 .and. this%spin(i+1) .eq. this%spin(i))then
                if(.not.allocated(alm))allocate(alm(2,0:this%lmax, 0:this%lmax))
                alm(1,:,:) = this%alm(:,:,i)
                alm(2,:,:) = this%alm(:,:,i+1)
                call alm2map_spin(this%nside, this%lmax, this%lmax, this%spin(i), alm, this%map(:,i:i+1))
                j = j + 2
                cycle
             endif
          endif
          stop "coop_healpix_maps_alm2map: nonzero spin maps must appear in pairs"
       enddo
    else
       i = 1
       do while(i.le. this%nmaps)
          if(this%spin(i).eq.0)then
             call alm2map(this%nside, this%lmax, this%lmax, this%alm(:,:,i:i), this%map(:,i))
             i = i + 1
          else
             if(i.lt. this%nmaps)then
                if(this%spin(i+1) .eq. this%spin(i))then
                   if(.not.allocated(alm))allocate(alm(2,0:this%lmax, 0:this%lmax))
                   alm(1,:,:) = this%alm(:,:,i)
                   alm(2,:,:) = this%alm(:,:,i+1)
                
                   call alm2map_spin(this%nside, this%lmax, this%lmax, this%spin(i), alm, this%map(:,i:i+1))
                   i = i + 2
                   cycle
                endif
             endif
             stop "coop_healpix_maps_alm2map: nonzero spin maps must appear in pairs"
          endif
       enddo
    endif
    if(allocated(alm))deallocate(alm)
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_maps_alm2map

  subroutine coop_healpix_filter_alm(this, fwhm, lpower, window, index_list)
    class(coop_healpix_maps) this
    real,optional::window(0:this%lmax)
    real,optional::fwhm
    real,optional::lpower
    COOP_INT,dimension(:), optional::index_list
    COOP_INT l
    COOP_SINGLE c, w(0:this%lmax)
    w = 1.
    if(present(fwhm))then
       c = sign((coop_sigma_by_fwhm * fwhm)**2/2., dble(fwhm))
       !$omp parallel do
       do l = 0,  this%lmax
          w(l) = w(l)*exp(-l*(l+1.)*c)
       enddo
       !$omp end parallel do
    endif
    if(present(lpower))then
       !$omp parallel do
       do l = 0,  this%lmax
          w(l) = w(l)*(l*(l+1.))**(lpower/2.)
       enddo
       !$omp end parallel do
    endif
    if(present(window))then
       !$omp parallel do
       do l = 0,  this%lmax
          w(l) = w(l)*window(l)
       enddo
       !$omp end parallel do       
    endif
    if(present(index_list))then
       !$omp parallel do
       do l = 0, this%lmax
          this%alm(l,:,index_list) = this%alm(l,:,index_list)*w(l)
       enddo
       !$omp end parallel do
    else
       !$omp parallel do
       do l = 0, this%lmax
          this%alm(l,:,:) = this%alm(l,:,:)*w(l)
       enddo
       !$omp end parallel do
    endif
  end subroutine coop_healpix_filter_alm


  subroutine split_angular_mode(n, qmap, umap, m, nr, fr)
    COOP_INT n, m, nr
    COOP_SINGLE qmap(-n:n,-n:n), umap(-n:n,-n:n)
    COOP_SINGLE fr(0:nr), w(0:nr), q, u, r, phi, fpoint
    COOP_INT i, j, ir
    fr = 0
    w = 0
    do i = -n, n
       do j = -n, n
          r = sqrt(real(i)**2 + real(j)**2)
          ir = floor(r)          
          if(ir .le. nr)then
             phi = m*COOP_POLAR_ANGLE(real(i), real(j))
             r = r - ir
             fpoint = qmap(i,j)*cos(phi) + umap(i,j)*sin(phi)
             fr(ir) = fr(ir) + fpoint*(1.d0-r)
             w(ir) = w(ir) + 1.d0-r
             if(ir.ne.nr)then
                fr(ir+1) = fr(ir+1)+fpoint*r
                w(ir+1) = w(ir+1)+r
             endif
          endif
       enddo
    enddo
    do ir = 0, nr
       if(w(ir).gt.0.)then
          fr(ir) = fr(ir)/w(ir)
       endif
    enddo
    do ir = n, nr !!damp the amplitude in the corners
       fr(ir) = fr(ir)*exp(-((ir-n+1.)/n*(1.414/0.414))**2)
    enddo
    if(m.ne.0) fr(0) = 0
  end subroutine split_angular_mode


  subroutine map_filter_modes(n, qmap, umap, ms)
    COOP_INT n, nm
    COOP_INT ms(:)
    COOP_SINGLE qmap(-n:n, -n:n), umap(-n:n, -n:n)
    real,dimension(:,:),allocatable::fr
    COOP_SINGLE r, phi, s1, s2
    COOP_INT i, j, ir, nr, im
    nm = size(ms)
    nr = ceiling(coop_sqrt2*n)+1
    allocate(fr(0:nr, nm))
    do im = 1, nm
       call split_angular_mode(n, qmap, umap, ms(im), nr, fr(0:nr, im))
    enddo
    qmap = 0
    umap = 0
    do i=-n, n
       do j=-n,n
          r = sqrt(real(i)**2 + real(j)**2)
          ir = floor(r)
          r = r - ir
          phi = COOP_POLAR_ANGLE(real(i), real(j))
          do  im =1, nm
             qmap(i,j) = qmap(i,j) + (fr(ir,im)*(1.-r)+fr(ir+1,im)*r)*cos(ms(im)*phi)
             umap(i,j) = umap(i,j) + (fr(ir,im)*(1.-r)+fr(ir+1,im)*r)*sin(ms(im)*phi)
             
          enddo
       enddo
    enddo
    deallocate(fr)
  end subroutine map_filter_modes


  subroutine coop_healpix_get_disc(nside, pix, disc)
    COOP_INT pix, nside
    type(coop_healpix_disc) disc
    COOP_REAL r
    disc%nside  = nside
    disc%center = pix
#ifdef HAS_HEALPIX
    call pix2ang_ring(nside, pix, disc%theta, disc%phi)
    call ang2vec(disc%theta, disc%phi, disc%nz)
#else
    stop "CANNOT FIND HEALPIX"
#endif
    disc%nx = (/  sin(disc%phi) , - cos(disc%phi) , 0.d0 /)
    call coop_vector_cross_product(disc%nz, disc%nx, disc%ny)

  end subroutine coop_healpix_get_disc

  subroutine coop_healpix_disc_pix2ang(disc, pix, r, phi)
    class(coop_healpix_disc) disc
    COOP_INT pix
    COOP_REAL r, phi, vec(3), x, y
    if(pix .eq. disc%center)then
       r = 0
       phi = 0
       return
    endif
#ifdef HAS_HEALPIX
    call pix2vec_ring(disc%nside, pix, vec)
#else
    stop "CANNOT FIND HEALPIX"
#endif
    r = COS2RADIUS( dot_product(vec, disc%nz) )
    x = dot_product(vec, disc%nx)
    y = dot_product(vec, disc%ny)
    phi = COOP_POLAR_ANGLE(x, y)
  end subroutine coop_healpix_disc_pix2ang

  subroutine coop_healpix_disc_ang2pix(disc, r, phi, pix)
    class(coop_healpix_disc) disc
    COOP_REAL r !!in unit of radian
    COOP_REAL phi, vec(3), cost, sint
    COOP_INT pix
    cost = RADIUS2COS(r)
    sint = sqrt(1.d0 - cost**2)
    vec = sint*cos(phi)* disc%nx + sint*sin(phi)*disc%ny + cost*disc%nz
#ifdef HAS_HEALPIX
    call vec2pix_ring(disc%nside, vec, pix)
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_disc_ang2pix


  subroutine coop_healpix_disc_pix2xy(disc, pix, x, y)
    class(coop_healpix_disc) disc
    COOP_INT pix
    COOP_REAL r, phi, vec(3), x, y
    if(pix .eq. disc%center)then
       x = 0
       y = 0
       return
    endif
#ifdef HAS_HEALPIX
    call pix2vec_ring(disc%nside, pix, vec)
#else
    stop "CANNOT FIND HEALPIX"
#endif
    r = COS2RADIUS( dot_product(vec, disc%nz) )
    x = dot_product(vec, disc%nx)
    y = dot_product(vec, disc%ny)
    r = r/sqrt(x**2+y**2)
    x = r * x
    y = r * y
  end subroutine coop_healpix_disc_pix2xy


  subroutine coop_healpix_disc_xy2pix(disc, x, y, pix)
    class(coop_healpix_disc) disc
    COOP_REAL x, y !!in unit of radian
    COOP_REAL vec(3), cost, sint, r
    COOP_INT pix
    r = sqrt(x**2+y**2)
    if(r.lt.1.d-8)then
       pix = disc%center
       return
    endif
    cost = RADIUS2COS(r)
    sint = sqrt(1.d0 - cost**2)
    vec = sint*(x/r)* disc%nx + sint*(y/r)*disc%ny + cost*disc%nz
#ifdef HAS_HEALPIX
    call vec2pix_ring(disc%nside, vec, pix)
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_disc_xy2pix

  subroutine coop_healpix_rotate_qu(qu, phi)
    COOP_SINGLE qu(2)
    COOP_REAL phi, cosp, sinp
    cosp = cos(2.d0*phi)
    sinp = sin(2.d0*phi)
    qu = (/ qu(1)*cosp + qu(2)*sinp,  -qu(1)*sinp + qu(2)*cosp /)
  end subroutine coop_healpix_rotate_qu


  subroutine coop_healpix_patch_get_fr0(patch, nvar, var)
    COOP_INT nvar
    type(coop_healpix_patch)::patch
    COOP_REAL var(nvar)
    call patch%get_radial_profile(1, 0)
    var = patch%fr(0:patch%n, 0, 1)
  end subroutine coop_healpix_patch_get_fr0





  subroutine coop_healpix_maps_stack(this, patch, spots_file, mask, hemisphere_direction, do_weight)
    COOP_UNKNOWN_STRING::spots_file
    COOP_INT,parameter::n_threads = 4
    class(coop_healpix_maps)::this
    type(coop_healpix_disc),dimension(n_threads)::disc
    type(coop_healpix_patch)::patch
    type(coop_healpix_patch),dimension(n_threads)::p, tmp
    type(coop_healpix_maps),optional::mask
    COOP_INT::ns
    COOP_REAL,dimension(:),allocatable::theta, phi, angle
    type(coop_file)::fp
    COOP_INT imap, ithread, i, pix, iaccept, ireject
    COOP_REAL,optional:: hemisphere_direction(2)
    logical,optional::do_weight
    COOP_REAL hcos, hsin
#ifdef HAS_HEALPIX
    if(.not. coop_file_exists(spots_file))then
       write(*,*) "Spots file not found: "//trim(spots_file)
       stop
    endif
    ns = coop_file_numlines(spots_file)
    if(ns .eq. 0)then
       write(*,*) "Spots file empty"
       stop
    endif
    allocate(theta(ns), phi(ns), angle(ns))
    call fp%open(spots_file)    
    do i=1, ns
       read(fp%unit, *) theta(i), phi(i), angle(i)
    enddo
    call fp%close()
    write(*, *) "stacking "//trim(coop_num2str(ns))//" patchs"
    patch%image = 0.d0
    patch%nstack = 0.d0
    patch%nstack_raw = 0
    do ithread=1, n_threads
       p(ithread) = patch
       tmp(ithread) = patch
    enddo
    if(present(hemisphere_direction))then
       hcos = cos(hemisphere_direction(1))
       hsin = sin(hemisphere_direction(1))
       !$omp parallel do private(i, pix) 
       do ithread = 1, n_threads
          do i=ithread, ns, n_threads
             if(hcos * cos(theta(i)) + hsin * sin(theta(i)) * cos(phi(i) - hemisphere_direction(2)) .gt. 0.d0)then
                call ang2pix_ring(this%nside, theta(i), phi(i), pix)
                call coop_healpix_get_disc(this%nside, pix, disc(ithread))
                if(present(mask))then
                   call coop_healpix_stack_on_patch(this, disc(ithread), angle(i), p(ithread), tmp(ithread), mask)    
                else
                   call coop_healpix_stack_on_patch(this, disc(ithread), angle(i), p(ithread), tmp(ithread) )
                endif
             endif
          enddo
       enddo
       !$omp end parallel do   
    else
       !$omp parallel do private(i, pix)
       do ithread = 1, n_threads
          do i=ithread, ns, n_threads
             call ang2pix_ring(this%nside, theta(i), phi(i), pix)
             call coop_healpix_get_disc(this%nside, pix, disc(ithread))
             if(present(mask))then
                call coop_healpix_stack_on_patch(this, disc(ithread), angle(i), p(ithread), tmp(ithread), mask)    
             else
                call coop_healpix_stack_on_patch(this, disc(ithread), angle(i), p(ithread), tmp(ithread) )
             endif
          enddo
       enddo
       !$omp end parallel do
    endif
    do ithread = 1, n_threads
       patch%image = patch%image + p(ithread)%image
       patch%nstack = patch%nstack + p(ithread)%nstack
       patch%nstack_raw = patch%nstack_raw + p(ithread)%nstack_raw
       call p(ithread)%free()
       call tmp(ithread)%free()
    enddo
    if(patch%nstack_raw .ne. 0)then
       if(present(do_weight))then
          if(do_weight)then
             do i=1, patch%nmaps
                patch%image(:, :, i) = patch%image(:, :, i)/max(patch%nstack, 1.d0)
             enddo
          else
             patch%image = patch%image/patch%nstack_raw
          endif
       else
          do i=1, patch%nmaps
             patch%image(:, :, i) = patch%image(:, :, i)/max(patch%nstack, 1.d0)
          enddo
       endif
    else
       write(*,*) "warning: no patches has been found"
       patch%image = 0.d0
    endif
    deallocate(theta, phi, angle)
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_maps_stack



  subroutine coop_healpix_maps_multstack(this, npatches, patch, lowercut, uppercut, spots_file, mask)
    COOP_UNKNOWN_STRING::spots_file
    COOP_INT,parameter::n_threads = 4
    class(coop_healpix_maps)::this
    type(coop_healpix_disc),dimension(n_threads)::disc
    COOP_INT npatches
    type(coop_healpix_patch)::patch(npatches)
    type(coop_healpix_patch)::p(n_threads, npatches), tmp(n_threads)
    type(coop_healpix_maps),optional::mask
    COOP_INT::ns
    COOP_REAL,dimension(:),allocatable::theta, phi, angle, col4
    type(coop_file)::fp
    COOP_INT imap, ithread, i, pix, iaccept, ireject, ip
    COOP_REAL hcos, hsin, maxcol4, mincol4, dcol4
    COOP_SINGLE uppercut, lowercut
#ifdef HAS_HEALPIX
    if(.not. coop_file_exists(spots_file))then
       write(*,*) "Spots file not found: "//trim(spots_file)
       stop
    endif
    ns = coop_file_numlines(spots_file)
    if(ns .eq. 0)then
       write(*,*) "Spots file empty"
       stop
    endif
    allocate(theta(ns), phi(ns), angle(ns), col4(ns))
    call fp%open(spots_file)    
    do i=1, ns
       read(fp%unit, *) theta(i), phi(i), angle(i), col4(i)
    enddo
    call fp%close()
    write(*, *) "stacking "//trim(coop_num2str(ns))//" patchs"
    maxcol4 = min(maxval(col4), uppercut)
    mincol4 = max(minval(col4), lowercut)
    dcol4 = (maxcol4 - mincol4)/npatches
    do ip = 1, npatches
       patch(ip)%image = 0.d0
       patch(ip)%nstack = 0.d0
       patch(ip)%nstack_raw = 0
       do ithread=1, n_threads
          p(ithread, ip) = patch(ip)
       enddo
    enddo
    do ithread = 1, n_threads
       tmp(ithread) = patch(1)
    enddo
    if(dcol4 .gt. 0.d0)then
       !$omp parallel do private(i, pix, ip)
       do ithread = 1, n_threads
          do i=ithread, ns, n_threads
             ip = floor((col4(i) - mincol4)/dcol4 + 1)
             if(ip .ge. 1 .and. ip .le. npatches)then
                call ang2pix_ring(this%nside, theta(i), phi(i), pix)
                call coop_healpix_get_disc(this%nside, pix, disc(ithread))
                if(present(mask))then
                   call coop_healpix_stack_on_patch(this, disc(ithread), angle(i), p(ithread, ip), tmp(ithread), mask)    
                else
                   call coop_healpix_stack_on_patch(this, disc(ithread), angle(i), p(ithread, ip), tmp(ithread) )
                endif
             endif
          enddo
       enddo
       !$omp end parallel do
       do ip=1, npatches
          patch(ip)%caption =  "$I \le "//trim(coop_num2str(nint(mincol4 + ip*dcol4)))//" \mu K$"
          write(*,*) trim(coop_num2str(ip))//": I < "//trim(coop_num2str(mincol4 + ip*dcol4))
       enddo
    else
       dcol4 = (ns + 0.1d0)/npatches
       !$omp parallel do private(i, pix, ip)
       do ithread = 1, n_threads
          do i=ithread, ns, n_threads
             ip = ceiling(i/dcol4)
             if(ip .ge. 1 .and. ip .le. npatches)then
                call ang2pix_ring(this%nside, theta(i), phi(i), pix)
                call coop_healpix_get_disc(this%nside, pix, disc(ithread))
                if(present(mask))then
                   call coop_healpix_stack_on_patch(this, disc(ithread), angle(i), p(ithread, ip), tmp(ithread), mask)    
                else
                   call coop_healpix_stack_on_patch(this, disc(ithread), angle(i), p(ithread, ip), tmp(ithread) )
                endif
             endif
          enddo
       enddo
       !$omp end parallel do
       do ip=1, npatches
          i  = floor(dcol4*ip)
          patch(ip)%caption =  "$I \le "//trim(coop_num2str(nint(col4(i))))//" \mu K$"
          write(*,*) trim(coop_num2str(ip))//": I < "//trim(coop_num2str(nint(col4(i))))//" \mu K$"
       enddo
    endif
    do ithread=1, n_threads
       call tmp(ithread)%free()
    enddo
    do ip=1, npatches
       do ithread = 1, n_threads
          patch(ip)%image = patch(ip)%image + p(ithread, ip)%image
          patch(ip)%nstack = patch(ip)%nstack + p(ithread, ip)%nstack
          patch(ip)%nstack_raw = patch(ip)%nstack_raw + p(ithread, ip)%nstack_raw
          call p(ithread, ip)%free()
       enddo
       if(patch(ip)%nstack_raw .gt. 0)then
          do imap = 1, patch(ip)%nmaps
             patch(ip)%image(:,:, imap) = patch(ip)%image(:,:, imap)/max(patch(ip)%nstack, 1.d0)
          enddo
       endif
    enddo
    deallocate(theta, phi, angle, col4)
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_maps_multstack


  subroutine coop_healpix_maps_stack_with_listpix(this, patch, listpix, listangle, mask, hemisphere_direction)
    type(coop_list_integer)::listpix
    type(coop_list_real)::listangle
    COOP_INT,parameter::n_threads = 4
    class(coop_healpix_maps)::this
    type(coop_healpix_disc),dimension(n_threads)::disc
    type(coop_healpix_patch)::patch
    type(coop_healpix_patch),dimension(n_threads)::p, tmp
    type(coop_healpix_maps),optional::mask
    type(coop_file)::fp
    COOP_REAL,optional::hemisphere_direction(2)
    COOP_INT imap, ithread, i, pix, iaccept, ireject
    COOP_REAL angle, theta, phi, hsint, hcost
    COOP_REAL hcos, hsin
#ifdef HAS_HEALPIX
    patch%image = 0.d0
    patch%nstack = 0.d0
    patch%nstack_raw = 0
    do ithread=1, n_threads
       p(ithread) = patch
       tmp(ithread) = patch
    enddo
    if(present(hemisphere_direction))then
       hsint = sin(hemisphere_direction(1))
       hcost = cos(hemisphere_direction(1))
       !$omp parallel do private(i, pix, angle)
       do ithread = 1, n_threads
          do i=ithread, listpix%n, n_threads
             pix = listpix%element(i)
             call pix2ang_ring(this%nside, pix, theta, phi)
             if(cos(theta)*hcost+sin(theta)*hsint*cos(phi-hemisphere_direction(2)) .ge. 0.d0)then
                angle = listangle%element(i)
                call coop_healpix_get_disc(this%nside, pix, disc(ithread))
                if(present(mask))then
                   call coop_healpix_stack_on_patch(this, disc(ithread), angle, p(ithread), tmp(ithread), mask)    
                else
                   call coop_healpix_stack_on_patch(this, disc(ithread), angle, p(ithread), tmp(ithread) )
                endif
             endif
          enddo
       enddo
       !$omp end parallel do
    else
       !$omp parallel do private(i, pix, angle)
       do ithread = 1, n_threads
          do i=ithread, listpix%n, n_threads
             pix = listpix%element(i)
             angle = listangle%element(i)
             call coop_healpix_get_disc(this%nside, pix, disc(ithread))
             if(present(mask))then
                call coop_healpix_stack_on_patch(this, disc(ithread), angle, p(ithread), tmp(ithread), mask)    
             else
                call coop_healpix_stack_on_patch(this, disc(ithread), angle, p(ithread), tmp(ithread) )
             endif
          enddo
       enddo
       !$omp end parallel do
    endif
    do ithread = 1, n_threads
       patch%image = patch%image + p(ithread)%image
       patch%nstack = patch%nstack + p(ithread)%nstack
       patch%nstack_raw = patch%nstack_raw + p(ithread)%nstack_raw
       call p(ithread)%free()
       call tmp(ithread)%free()
    enddo
    do imap = 1, patch%nmaps
       patch%image(:, :, imap) = patch%image(:,:,imap)/max(patch%nstack, 1.d0)
    enddo
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_maps_stack_with_listpix



  subroutine coop_healpix_maps_stack_north_south(this, patch_north, patch_south, listpix, listangle, hemisphere_direction, mask)
    type(coop_list_integer)::listpix
    type(coop_list_real)::listangle
    COOP_INT,parameter::n_threads = 4
    class(coop_healpix_maps)::this
    type(coop_healpix_disc),dimension(n_threads)::disc
    type(coop_healpix_patch)::patch_north, patch_south
    type(coop_healpix_patch),dimension(n_threads)::p_south, tmp, p_north
    type(coop_healpix_maps),optional::mask
    type(coop_file)::fp
    COOP_REAL::hemisphere_direction(2)
    COOP_INT imap, ithread, i, pix, iaccept, ireject
    COOP_REAL angle, theta, phi, hsint, hcost
    COOP_REAL hcos, hsin
#ifdef HAS_HEALPIX
    patch_north%image = 0.d0
    patch_north%nstack = 0.d0
    patch_north%nstack_raw = 0
    patch_south%image = 0.d0
    patch_south%nstack = 0.d0
    patch_south%nstack_raw = 0
    do ithread=1, n_threads
       p_north(ithread) = patch_north
       p_south(ithread) = patch_south
       tmp(ithread) = patch_north
    enddo
    hsint = sin(hemisphere_direction(1))
    hcost = cos(hemisphere_direction(1))
    !$omp parallel do private(i, pix, angle)
    do ithread = 1, n_threads
       do i=ithread, listpix%n, n_threads
          pix = listpix%element(i)
          call pix2ang_ring(this%nside, pix, theta, phi)
          angle = listangle%element(i)
          call coop_healpix_get_disc(this%nside, pix, disc(ithread))
          if(cos(theta)*hcost+sin(theta)*hsint*cos(phi-hemisphere_direction(2)) .lt. 0.d0)then
             if(present(mask))then
                call coop_healpix_stack_on_patch(this, disc(ithread), angle, p_north(ithread), tmp(ithread), mask)    
             else
                call coop_healpix_stack_on_patch(this, disc(ithread), angle, p_north(ithread), tmp(ithread) )
             endif
          else
             if(present(mask))then
                call coop_healpix_stack_on_patch(this, disc(ithread), angle, p_south(ithread), tmp(ithread), mask)    
             else
                call coop_healpix_stack_on_patch(this, disc(ithread), angle, p_south(ithread), tmp(ithread) )
             endif
          endif
       enddo
    enddo
    !$omp end parallel do
    do ithread = 1, n_threads
       patch_north%image = patch_north%image + p_north(ithread)%image
       patch_north%nstack = patch_north%nstack + p_north(ithread)%nstack
       patch_north%nstack_raw = patch_north%nstack_raw + p_north(ithread)%nstack_raw
       call p_north(ithread)%free()
       patch_south%image = patch_south%image + p_south(ithread)%image
       patch_south%nstack = patch_south%nstack + p_south(ithread)%nstack
       patch_south%nstack_raw = patch_south%nstack_raw + p_south(ithread)%nstack_raw
       call p_south(ithread)%free()
       call tmp(ithread)%free()
    enddo
    do imap = 1, patch_north%nmaps
       patch_north%image(:,:,imap) = patch_north%image(:,:,imap)/max(patch_north%nstack, 1.d0)
       patch_south%image(:,:,imap) = patch_south%image(:,:,imap)/max(patch_south%nstack, 1.d0)
    enddo
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_maps_stack_north_south



  subroutine coop_healpix_maps_stack_with_covariance(this, patch, spots_file, getvar, nvar, mean, cov, mask, hemisphere_direction)
    COOP_UNKNOWN_STRING::spots_file
    COOP_INT,parameter::n_threads = 4
    class(coop_healpix_maps)::this
    type(coop_healpix_disc),dimension(n_threads)::disc
    type(coop_healpix_patch)::patch
    type(coop_healpix_patch),dimension(n_threads)::p, tmp
    COOP_INT nvar
    external getvar
    COOP_REAL cov(nvar, nvar), mean(nvar),  covtmp(nvar, nvar, n_threads), meantmp(nvar, n_threads)
    type(coop_healpix_maps), optional::mask
    COOP_INT::ns
    COOP_REAL,dimension(:),allocatable::theta, phi, angle
    type(coop_file)::fp
    COOP_INT imap, ithread, i, pix, j
    COOP_REAL, optional::hemisphere_direction(2)
    COOP_REAL hcos, hsin
#ifdef HAS_HEALPIX
    ns = coop_file_numlines(spots_file)
    allocate(theta(ns), phi(ns), angle(ns))
    call fp%open(spots_file)
    do i=1, ns
       read(fp%unit, *) theta(i), phi(i), angle(i)
    enddo
    call fp%close()
    patch%image = 0.d0
    patch%nstack = 0.d0
    patch%nstack_raw = 0
    do ithread=1, n_threads
       p(ithread) = patch
       tmp(ithread) = patch
    enddo
    covtmp = 0.d0
    meantmp = 0.d0
    cov = 0.d0
    mean = 0.d0
    if(present(hemisphere_direction))then
       hcos = cos(hemisphere_direction(1))
       hsin = sin(hemisphere_direction(1))
       !$omp parallel do private(i, pix)
       do ithread = 1, n_threads
          do i=ithread, ns, n_threads
             if(hcos*cos(theta(i)) + hsin*sin(theta(i))*cos(phi(i)-hemisphere_direction(2)) .gt. 0.d0)then
                call ang2pix_ring(this%nside, theta(i), phi(i), pix)
                call coop_healpix_get_disc(this%nside, pix, disc(ithread))
                if(present(mask))then
                   call coop_healpix_stack_on_patch_with_covariance(this, disc(ithread), angle(i), p(ithread), tmp(ithread), getvar, nvar, meantmp(:, ithread), covtmp(:,:,ithread), mask)  
                else
                   call coop_healpix_stack_on_patch_with_covariance(this, disc(ithread), angle(i), p(ithread), tmp(ithread), getvar, nvar, meantmp(:, ithread), covtmp(:,:,ithread))
                endif
             endif
          enddo
       enddo
       !$omp end parallel do
    else
       !$omp parallel do private(i, pix)
       do ithread = 1, n_threads
          do i=ithread, ns, n_threads
             call ang2pix_ring(this%nside, theta(i), phi(i), pix)
             call coop_healpix_get_disc(this%nside, pix, disc(ithread))
             if(present(mask))then
                call coop_healpix_stack_on_patch_with_covariance(this, disc(ithread), angle(i), p(ithread), tmp(ithread), getvar, nvar, meantmp(:, ithread), covtmp(:,:,ithread), mask)  
             else
                call coop_healpix_stack_on_patch_with_covariance(this, disc(ithread), angle(i), p(ithread), tmp(ithread), getvar, nvar, meantmp(:, ithread), covtmp(:,:,ithread))
             endif
          enddo
       enddo
       !$omp end parallel do
    endif
    do ithread = 1, n_threads
       patch%image = patch%image + p(ithread)%image
       patch%nstack = patch%nstack + p(ithread)%nstack
       patch%nstack_raw = patch%nstack_raw + p(ithread)%nstack_raw
       cov = cov+covtmp(:,:,ithread)
       mean = mean + meantmp(:, ithread)
       call p(ithread)%free()
       call tmp(ithread)%free()
    enddo
    if(patch%nstack_raw .eq. 0) stop "Nothing stacked"
    mean = mean/patch%nstack_raw
    do imap = 1, patch%nmaps
       patch%image(:,:,imap) = patch%image(:,:,imap)/max(patch%nstack, 1.d0)
    enddo
!!$    patch%image = patch%image/patch%nstack_raw
    do j=1, nvar
       do i=1, j
          cov(i,j) = cov(i,j)/patch%nstack_raw - mean(i)*mean(j)
       enddo
    enddo
    do j=1, nvar
       do i=j+1, nvar
          cov(i,j) =cov(j, i)
       enddo
    enddo
    deallocate(theta, phi, angle)
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_maps_stack_with_covariance
  

  subroutine coop_healpix_stack_on_patch_with_covariance(this, disc, angle, patch, tmp_patch, getvar, nvar, mean, cov, mask)
    class(coop_healpix_maps)::this
    type(coop_healpix_disc)::disc
    COOP_REAL angle
    type(coop_healpix_patch)::patch, tmp_patch
    type(coop_healpix_maps),optional::mask
    COOP_INT i, j
    COOP_INT::nvar
    COOP_REAL::cov(nvar, nvar), mean(nvar)
    external::getvar
    COOP_REAL var(nvar)
    if(present(mask))then
       call coop_healpix_fetch_patch(this, disc, angle, tmp_patch, mask)
       if(present(mask) .and. sum(tmp_patch%nstack*tmp_patch%indisk) .lt. patch%num_indisk_tol)then
          return
       endif
    else
       call coop_healpix_fetch_patch(this, disc, angle, tmp_patch)
    endif
    patch%image = patch%image + tmp_patch%image
    patch%nstack = patch%nstack + tmp_patch%nstack
    patch%nstack_raw = patch%nstack_raw + tmp_patch%nstack_raw
    call getvar(tmp_patch, nvar, var)
    mean = mean + var
    do j=1, nvar
       do i=1, j
          cov(i, j) = cov(i, j) + var(i)*var(j)
       enddo
    enddo
  end subroutine coop_healpix_stack_on_patch_with_covariance


  subroutine coop_healpix_stack_on_patch(this, disc, angle, patch, tmp_patch, mask)
    class(coop_healpix_maps) this
    type(coop_healpix_disc) disc
    type(coop_healpix_maps),optional::mask
    COOP_REAL angle
    type(coop_healpix_patch) patch, tmp_patch
    if(present(mask))then
       call coop_healpix_fetch_patch(this, disc, angle, tmp_patch, mask)
       if(present(mask) .and. sum(tmp_patch%nstack*tmp_patch%indisk) .lt. patch%num_indisk_tol)return
    else
       call coop_healpix_fetch_patch(this, disc, angle, tmp_patch)
    endif
    patch%image = patch%image + tmp_patch%image
    patch%nstack = patch%nstack + tmp_patch%nstack
    patch%nstack_raw = patch%nstack_raw + tmp_patch%nstack_raw
  end subroutine coop_healpix_stack_on_patch


  subroutine coop_healpix_fetch_patch(this, disc, angle, patch, mask)
    class(coop_healpix_maps)::this
    type(coop_healpix_disc) disc
    type(coop_healpix_maps),optional::mask
    COOP_REAL angle
    type(coop_healpix_patch) patch
    COOP_INT i, j, pix
    COOP_REAL x, y, r, phi
    COOP_SINGLE qu(2)
    if(.not. present(mask))patch%nstack = 1.d0
    patch%nstack_raw  = 1
    select case(trim(patch%genre))
    case("T", "E", "B", "I", "zeta")
       do j = -patch%n, patch%n
          do i = -patch%n, patch%n
             x = patch%dr * i
             y = patch%dr * j
             r = sqrt(x**2+y**2)
             phi = COOP_POLAR_ANGLE(x, y) + angle
             call coop_healpix_disc_ang2pix(disc, r, phi, pix)
             if(present(mask))then
                patch%nstack(i, j) = mask%map(pix, 1)
                patch%image(i, j, 1) = this%map(pix,1)*mask%map(pix,1)
             else 
                patch%image(i, j, 1) = this%map(pix,1)
             endif
          enddo
       enddo
    case("QU")
       do j = -patch%n, patch%n
          do i = -patch%n, patch%n
             x = patch%dr * i
             y = patch%dr * j
             r = sqrt(x**2+y**2)
             phi = COOP_POLAR_ANGLE(x, y) + angle
             call coop_healpix_disc_ang2pix(disc, r, phi, pix)
             qu = this%map(pix, this%iq:this%iu)
             call coop_healpix_rotate_qu(qu, angle)
             if(present(mask))then
                patch%nstack(i, j) = mask%map(pix, 1)
                patch%image(i, j, 1:2) = qu * mask%map(pix,1)
             else 
                patch%image(i, j, 1:2) = qu
             endif
          enddo
       enddo
    case("QrUr")
       do j = -patch%n, patch%n
          do i = -patch%n, patch%n
             x = patch%dr * i
             y = patch%dr * j
             r = sqrt(x**2+y**2)
             phi = COOP_POLAR_ANGLE(x, y) + angle
             call coop_healpix_disc_ang2pix(disc, r, phi, pix)
             qu = this%map(pix, this%iq:this%iu)
             call coop_healpix_rotate_qu(qu, phi)
             qu = -qu  !!the nonsense - sign from WMAP convention!
             if(present(mask))then
                patch%nstack(i, j) = mask%map(pix, 1)
                patch%image(i, j, 1:2) = qu * mask%map(pix,1)
             else 
                patch%image(i, j, 1:2) = qu
             endif
          enddo
       enddo
    case default
       write(*,*) "coop_healpix_fetch_patch: Unknown stack_option"//trim(patch%genre)
       stop
    end select
  end subroutine coop_healpix_fetch_patch

  subroutine coop_healpix_smooth_mapfile(mapfile, fwhm)
    COOP_UNKNOWN_STRING mapfile
    type(coop_healpix_maps) map
    COOP_REAL fwhm
    call map%read(mapfile)
    call coop_healpix_maps_smooth(map, fwhm)
    call map%write(trim(coop_file_add_postfix(trim(mapfile),"_smoothed_fwhm"//trim(coop_num2str(nint(fwhm/coop_SI_arcmin)))//"arcmin")))
    write(*,*) "output: "//trim(coop_file_add_postfix(trim(mapfile),"_smoothed_fwhm"//trim(coop_num2str(nint(fwhm/coop_SI_arcmin)))//"arcmin"))
    call map%free()
  end subroutine coop_healpix_smooth_mapfile

  subroutine coop_healpix_maps_smooth(map, fwhm, index_list, l_lower, l_upper)
    class(coop_healpix_maps) map
    COOP_REAL fwhm
    COOP_INT, optional::l_lower, l_upper
    COOP_INT,dimension(:),optional::index_list
    COOP_INT lmax
    if(fwhm .gt. 0.d0)then
       lmax = min(ceiling(3./max(abs(fwhm)*coop_sigma_by_fwhm, 1.d-6)), map%nside*2, coop_healpix_default_lmax)
    else
       lmax = min(map%nside*2, coop_healpix_default_lmax)
    endif
    if(present(l_upper))then
       lmax = min(lmax, l_upper)
    endif
    if(lmax .lt. 2) stop "Huge smoothing scale cannot be done!"    
    if(lmax*abs(fwhm).lt.0.01)return
    if(present(index_list))then
       if(any(index_list .gt. map%nmaps)) stop "smooth: index_list overflow"
       call map%map2alm(lmax, index_list)
       if(present(l_lower)) map%alm(0:l_lower-1, :, :) = 0.
       call map%filter_alm(fwhm = real(fwhm), index_list = index_list)
       call map%alm2map(index_list)
    else
       call map%map2alm(lmax)
       if(present(l_lower)) map%alm(0:l_lower-1, :, :) = 0.       
       call map%filter_alm(fwhm = real(fwhm))
       call map%alm2map()
    endif
  end subroutine coop_healpix_maps_smooth

    subroutine coop_healpix_maps_smooth_with_window(map, fwhm, lmax, window, index_list)
    class(coop_healpix_maps) map
    COOP_REAL fwhm
    COOP_INT lmax
    COOP_INT,dimension(:),optional::index_list
    COOP_SINGLE window(0:lmax)
    if(present(index_list))then
       if(any(index_list .gt. map%nmaps)) stop "smooth: index_list overflow"
       call map%map2alm(lmax, index_list)
       call map%filter_alm(fwhm = real(fwhm), window = window, index_list = index_list)
       call map%alm2map(index_list)
    else
       call map%map2alm(lmax)
       call map%filter_alm(fwhm = real(fwhm), window = window)
       call map%alm2map()
    endif
  end subroutine coop_healpix_maps_smooth_with_window


  subroutine coop_healpix_getQU(Emap_file, QUmap_file)
    COOP_UNKNOWN_STRING Emap_file, QUmap_file
    type(coop_healpix_maps) hge, hgqu
    call hge%read(Emap_file,  nmaps_wanted = 1)
    call hge%map2alm()
    call hgqu%init(nside = hge%nside, nmaps = 2, spin =  (/ 2, 2 /), lmax=hge%lmax)
    hgqu%alm(:, :, 1) = hge%alm(:, :, 1)
    hgqu%alm(:, :, 2) = 0
    call hgqu%alm2map()
    call hgqu%write(QUmap_file)
    call hgqu%free()
    call hge%free()
  end subroutine coop_healpix_getQU


  !!note that map and mask are converted to NESTED order after calling this
  subroutine coop_healpix_maps_get_spots(map, spots,  spot_type, threshold, threshold_pol, mask)
    COOP_UNKNOWN_STRING spot_type
    type(coop_list_realarr)::spots
    COOP_REAL,optional::threshold, threshold_pol
    COOP_REAL theta, phi, rotate_angle, fcut, fcut2
    class(coop_healpix_maps)map
    type(coop_healpix_maps),optional::mask
    COOP_REAL total_weight
    COOP_INT i, iq, iu, j
    COOP_INT nneigh, list(8)
    logical do_mask
#ifdef HAS_HEALPIX
    call spots%init()
    do_mask = .false.
    if(present(mask))then
       if(mask%nside .ne. map%nside)then
          write(*,*) "map nside = ", map%nside
          write(*,*) "mask nside = ", mask%nside
          stop "coop_healpix_export_spots: mask and map must have the same nside"
       endif
       total_weight = count(mask%map(:,1).gt.0.5)
       do_mask = .true.
    else
       total_weight = map%npix
    endif
    call map%convert2nested
    if(do_mask) call mask%convert2nested

    select case(trim(spot_type))
    case("Tmax", "Emax", "Bmax", "zetamax")  !!random orientation

       if(present(threshold))then
          if(threshold.lt.coop_healpix_max_threshold)then          
             if(do_mask)then
                fcut = threshold*sqrt(sum(dble(map%map(:,1))**2, mask%map(:,1).gt.0.5)/total_weight)
             else
                fcut = threshold*sqrt(sum(dble(map%map(:,1))**2)/total_weight)
             endif
          else
             fcut = -1.e30
          endif
       else
          fcut = -1.e30
       endif

       do i=0, map%npix-1
          if(do_mask)then
             if( map%map(i,1) .lt. fcut .or. mask%map(i, 1) .le. 0.5 ) cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if ( all(map%map(list(1:nneigh),1).lt.map%map(i,1)) .and. all(mask%map(list(1:nneigh),1) .gt. 0.5 ) )then
                call random_number(rotate_angle)
                rotate_angle = rotate_angle*coop_2pi
                call pix2ang_nest(map%nside, i, theta, phi)
                call spots%push( real( (/ theta, phi, rotate_angle /) ) )
             endif
          else
             if(map%map(i,1).lt.fcut )cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if(all(map%map(list(1:nneigh),1).lt. map%map(i,1)))then
                call random_number(rotate_angle)
                rotate_angle = rotate_angle*coop_2pi
                call pix2ang_nest(map%nside, i, theta, phi)
                call spots%push( real( (/ theta, phi, rotate_angle /) ) )
             endif
          endif
       enddo
    case("Tmin", "Emin", "Bmin", "zetamin")  !!random orientation
       if(present(threshold))then
          if(threshold .lt. coop_healpix_max_threshold)then
             if(do_mask)then
                fcut = - threshold*sqrt(sum(dble(map%map(:,1))**2, mask%map(:,1).gt.0.5)/total_weight)
             else
                fcut = - threshold*sqrt(sum(dble(map%map(:,1))**2)/total_weight)
             endif
          else
             fcut = 1.e30
          endif
       else
          fcut = 1.e30
       endif
       do i=0, map%npix-1
          if(do_mask)then
             if(map%map(i,1).gt.fcut .or. mask%map(i, 1) .le. 0.5)cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if( all(map%map(list(1:nneigh),1) .gt. map%map(i,1)) .and. all(mask%map(list(1:nneigh), 1) .gt. 0.5 ) ) then
                call random_number(rotate_angle)
                rotate_angle = rotate_angle*coop_2pi
                call pix2ang_nest(map%nside, i, theta, phi)
                call spots%push( real( (/ theta, phi, rotate_angle /) ) )
             endif
          else
             if(map%map(i,1).gt.fcut )cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if(all(map%map(list(1:nneigh),1).gt.map%map(i,1)))then
                call random_number(rotate_angle)
                rotate_angle = rotate_angle*coop_2pi
                call pix2ang_nest(map%nside, i, theta, phi)
                call spots%push( real( (/ theta, phi, rotate_angle /) ) )
             endif
          endif
       enddo
    case("Tmax_QTUTOrient", "zetamax_qzuzOrient") !!oriented with QU, maxima of T
       if(map%nmaps .lt. 3) stop "For Tmax_QTUTOrient mode you need 3 maps"
       if(present(threshold))then
          if(threshold.lt.coop_healpix_max_threshold)then
             if(do_mask)then
                fcut = threshold*sqrt(sum(dble(map%map(:,1))**2, mask%map(:,1).gt.0.5)/total_weight)
             else
                fcut = threshold*sqrt(sum(dble(map%map(:,1))**2)/total_weight)
             endif
          else
             fcut=-1.e30
          endif
       else
          fcut = -1.e30
       endif
       if(present(threshold_pol))then
          if(threshold_pol.lt.coop_healpix_max_threshold)then
             if(do_mask)then
                fcut2 = threshold_pol**2*(sum(dble(map%map(:,2))**2+dble(map%map(:,3))**2, mask%map(:,1).gt.0.5)/total_weight)
             else
                fcut2 = threshold_pol**2*(sum(dble(map%map(:,2))**2+dble(map%map(:,3))**2)/total_weight)
             endif
          else
             fcut2 = -1.e30
          endif
       else
          fcut2 = -1.e30
       endif
       
       do i=0, map%npix-1
          if(do_mask)then
             if( map%map(i,1) .lt. fcut .or. map%map(i,2)**2+map%map(i,3)**2 .lt.fcut2 .or. mask%map(i,1) .le. 0.5)cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if(all(map%map(list(1:nneigh),1) .lt. map%map(i,1)) .and. all(mask%map(list(1:nneigh),1) .gt. 0.5)) then
                call pix2ang_nest(map%nside, i, theta, phi)
                rotate_angle = COOP_POLAR_ANGLE(map%map(i, 2), map%map(i, 3))/2.
                if(coop_random_unit().gt.0.5d0)rotate_angle=rotate_angle+coop_pi                                
                call spots%push( real( (/ theta, phi, rotate_angle /) ) )
             endif
          else
             if(map%map(i,1) .lt. fcut  .or. map%map(i,2)**2+map%map(i,3)**2 .lt.fcut2  )cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if(all(map%map(list(1:nneigh),1) .lt. map%map(i,1)))then
                call pix2ang_nest(map%nside, i, theta, phi)
                rotate_angle = COOP_POLAR_ANGLE(map%map(i, 2), map%map(i, 3))/2.
                if(coop_random_unit().gt.0.5d0)rotate_angle=rotate_angle+coop_pi                                
                call spots%push( real( (/ theta, phi, rotate_angle /) ) )
             endif
          endif
       enddo
    case("Tmin_QTUTOrient", "zetamin_ztutOrient") !!oriented with QU, minima of T
       if(map%nmaps .lt. 3) stop "For Tmin_QTUTOrient mode you need 3 maps"
       if(present(threshold))then
          if(threshold.lt.coop_healpix_max_threshold)then
             if(do_mask)then
                fcut = -threshold*sqrt(sum(dble(map%map(:,1))**2, mask%map(:,1).gt.0.5)/total_weight)
             else
                fcut = -threshold*sqrt(sum(dble(map%map(:,1))**2)/total_weight)
             endif
          else
             fcut = 1.e30
          endif
       else
          fcut = 1.e30
       endif
       do i=0, map%npix-1
          if(do_mask)then
             if( map%map(i,1) .gt. fcut .or. mask%map(i,1) .le. 0.5) cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if( all( map%map(list(1:nneigh), 1) .gt. map%map(i,1) ) .and. all(mask%map(list(1:nneigh),1) .gt. 0.5) ) then
                call pix2ang_nest( map%nside, i, theta, phi )
                rotate_angle = COOP_POLAR_ANGLE(map%map(i, 2), map%map(i, 3))/2.
                if(coop_random_unit().gt.0.5d0)rotate_angle=rotate_angle+coop_pi                                
                call spots%push( real( (/ theta, phi, rotate_angle /) ) )
             endif
          else
             if(map%map(i,1) .gt. fcut)cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if(all(map%map(list(1:nneigh),1).gt. map%map(i,1)))then
                call pix2ang_nest(map%nside, i, theta, phi)
                rotate_angle = COOP_POLAR_ANGLE(map%map(i, 2), map%map(i, 3))/2.
                if(coop_random_unit().gt.0.5d0)rotate_angle=rotate_angle+coop_pi                                
                call spots%push( real( (/ theta, phi, rotate_angle /) ) )
             endif
          endif
       enddo
    case("Pmax", "PTmax", "PZmax")
       select case(map%nmaps)
       case(2:3)
          iq = map%nmaps - 1
          iu = iq + 1
       case default
          stop "For polarization you need to specify Q, U maps in the input file."
       end select
       if(present(threshold))then
          if(threshold.lt.coop_healpix_max_threshold)then
             if(do_mask)then
                fcut = threshold**2*(sum(dble(map%map(:,iq))**2 + dble(map%map(:,iu))**2, mask%map(:,1).gt.0.5)/total_weight)
             else
                fcut = threshold**2*(sum(dble(map%map(:,iq))**2 + dble(map%map(:,iu))**2)/total_weight)
             endif
          else
             threshold = -1.e30
          endif
       else
          fcut = -1.e30
       endif
       do i=0, map%npix-1
          if(do_mask)then
             if( map%map(i,iq)**2 + map%map(i,iu)**2 .lt. fcut  .or.  mask%map(i, 1) .le. 0.5 ) cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if( all( map%map(list(1:nneigh),iq)**2 + map%map(list(1:nneigh),iu)**2 .lt. map%map(i,iq)**2 + map%map(i, iu)**2 ) .and. all(mask%map(list(1:nneigh), 1) .gt. 0.5) ) then
                rotate_angle = COOP_POLAR_ANGLE(map%map(i, iq), map%map(i, iu))/2.
                if(coop_random_unit().gt.0.5d0)rotate_angle=rotate_angle+coop_pi                                
                call pix2ang_nest(map%nside, i, theta, phi)
                call spots%push( real( (/ theta, phi, rotate_angle /) ) )
             endif
          else
             if( map%map(i, iq)**2 + map%map(i, iu)**2 .lt. fcut) cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if( all( map%map(list(1:nneigh), iq)**2 + map%map(list(1:nneigh), iu)**2 .lt. map%map(i, iq)**2 + map%map(i, iu)**2 ) )then
                rotate_angle = COOP_POLAR_ANGLE(map%map(i, iq), map%map(i, iu))/2.
                if(coop_random_unit().gt.0.5d0)rotate_angle=rotate_angle+coop_pi                                
                call pix2ang_nest(map%nside, i, theta, phi)
                call spots%push( real( (/ theta, phi, rotate_angle /) ) )
             endif
          endif
       enddo

    case("PmaxSortT", "PTmaxSortT")
       select case(map%nmaps)
       case(2:3)
          iq = map%nmaps - 1
          iu = iq + 1
       case default
          stop "For polarization you need to specify Q, U maps in the input file."
       end select
       if(present(threshold))then
          if(threshold.lt.coop_healpix_max_threshold)then
             if(do_mask)then
                fcut = threshold**2*(sum(dble(map%map(:,iq))**2 + dble(map%map(:,iu))**2, mask%map(:,1).gt.0.5)/total_weight)
             else
                fcut = threshold**2*(sum(dble(map%map(:,iq))**2 + dble(map%map(:,iu))**2)/total_weight)
             endif
          else
             fcut = -1.e30
          endif
       else
          fcut = -1.e30
       endif
       do i=0, map%npix-1
          if(do_mask)then
             if( map%map(i,iq)**2 + map%map(i,iu)**2 .lt. fcut  .or.  mask%map(i, 1) .le. 0.5 ) cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if( all( map%map(list(1:nneigh),iq)**2 + map%map(list(1:nneigh),iu)**2 .lt. map%map(i,iq)**2 + map%map(i, iu)**2 ) .and. all(mask%map(list(1:nneigh), 1) .gt. 0.5) ) then
                rotate_angle = COOP_POLAR_ANGLE(map%map(i, iq), map%map(i, iu))/2.
                if(coop_random_unit().gt.0.5d0)rotate_angle=rotate_angle+coop_pi
                call pix2ang_nest(map%nside, i, theta, phi)
                call spots%push(  (/ real(theta), real(phi), real(rotate_angle), map%map(i, 1) /) ) 
             endif
          else
             if( map%map(i, iq)**2 + map%map(i, iu)**2 .lt. fcut) cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if( all( map%map(list(1:nneigh), iq)**2 + map%map(list(1:nneigh), iu)**2 .lt. map%map(i, iq)**2 + map%map(i, iu)**2 ) )then
                rotate_angle = COOP_POLAR_ANGLE(map%map(i, iq), map%map(i, iu))/2.
                if(coop_random_unit().gt.0.5d0)rotate_angle=rotate_angle+coop_pi
                call pix2ang_nest(map%nside, i, theta, phi)
                call spots%push( (/ real(theta), real(phi), real(rotate_angle), map%map(i, 1) /) )
             endif
          endif
       enddo
       call spots%sort(4)
    case("Pmin", "PTmin", "PZmin")
       select case(map%nmaps)
       case(2:3)
          iq = map%nmaps - 1
          iu = iq + 1
       case default
          stop "For polarization you need to specify Q, U maps in the input file."
       end select
       if(present(threshold))then
          if(threshold .lt.coop_healpix_max_threshold)then
             if(do_mask)then
                fcut = threshold**2*(sum(dble(map%map(:, iq))**2+dble(map%map(:, iu))**2, mask%map(:,1).gt.0.5)/total_weight)
             else
                fcut = threshold**2*(sum(dble(map%map(:, iq))**2+dble(map%map(:, iu))**2)/total_weight)
             endif
          else
             fcut = 1.d30
          endif
       else
          fcut = 1.d30
       endif
       do i=0, map%npix-1
          if(do_mask)then
             if(map%map(i,iq)**2 + map%map(i,iu)**2 .gt. fcut .or. mask%map(i, 1) .le. 0.5 ) cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if(all(map%map(list(1:nneigh), iq)**2 + map%map(list(1:nneigh), iu)**2 .gt. map%map(i, iq)**2 + map%map(i, iu)**2) .and. all(mask%map(list(1:nneigh), 1) .gt. 0.5) ) then
                rotate_angle = COOP_POLAR_ANGLE(map%map(i,iq), map%map(i,iu))/2.
                if(coop_random_unit().gt.0.5d0)rotate_angle=rotate_angle+coop_pi                
                call pix2ang_nest(map%nside, i, theta, phi)
                call spots%push( real( (/ theta, phi, rotate_angle /) ) )
             endif
          else
             if(map%map(i, iq)**2 + map%map(i, iu)**2 .gt. fcut) cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if(all(map%map(list(1:nneigh), iq)**2 + map%map(list(1:nneigh), iu)**2 .gt. map%map(i, iq)**2 + map%map(i, iu)**2))then
                rotate_angle = COOP_POLAR_ANGLE(map%map(i, iq), map%map(i, iu))/2.
                if(coop_random_unit().gt.0.5d0)rotate_angle=rotate_angle+coop_pi
                call pix2ang_nest(map%nside, i, theta, phi)
                call spots%push( real( (/ theta, phi, rotate_angle /) ) )
             endif
          endif
       enddo
    case default
       write(*,*) trim(spot_type)
       stop "unknown spot type"
    end select
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_maps_get_spots


  subroutine coop_healpix_maps_get_listpix(map, listpix, listangle,  spot_type, threshold, mask)
    COOP_UNKNOWN_STRING spot_type
    type(coop_list_integer)::listpix
    type(coop_list_real)::listangle
    COOP_REAL,optional::threshold
    COOP_REAL theta, phi, rotate_angle, fcut
    class(coop_healpix_maps)map
    type(coop_healpix_maps),optional::mask
    COOP_REAL total_weight
    COOP_INT i, iq, iu, j, ipix
    COOP_INT nneigh, list(8)
    logical do_mask
#ifdef HAS_HEALPIX
    call listpix%init()
    call listangle%init()
    do_mask = .false.
    if(present(mask))then
       if(mask%nside .ne. map%nside) stop "coop_healpix_export_spots: mask and map must have the same nside"
       total_weight = count(mask%map(:,1).gt. 0.5)
       do_mask = .true.
    else
       total_weight = map%npix
    endif
    call map%convert2nested
    if(do_mask)call mask%convert2nested

    select case(trim(spot_type))
    case("Tmax", "Emax", "Bmax", "zetamax")  !!random orientation
       if(present(threshold))then
          if(threshold .lt. coop_healpix_max_threshold)then
             if(do_mask)then
                fcut = threshold*sqrt(sum(dble(map%map(:,1))**2, mask%map(:,1).gt.0.5)/total_weight)
             else
                fcut = threshold*sqrt(sum(dble(map%map(:,1))**2)/total_weight)
             endif
          else
             fcut = -1.e30
          endif
       else
          fcut = -1.e30
       endif
       do i=0, map%npix-1
          if(do_mask)then
             if( map%map(i,1) .lt. fcut .or. mask%map(i, 1) .le. 0.5 ) cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if ( all(map%map(list(1:nneigh),1).lt.map%map(i,1)) .and. all(mask%map(list(1:nneigh),1) .gt. 0.5 ) )then
                call random_number(rotate_angle)
                rotate_angle = rotate_angle*coop_2pi
                call nest2ring(map%nside, i, ipix)
                call listpix%push(ipix)
                call listangle%push(real(rotate_angle))
             endif
          else
             if(map%map(i,1).lt.fcut )cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if(all(map%map(list(1:nneigh),1).lt. map%map(i,1)))then
                call random_number(rotate_angle)
                rotate_angle = rotate_angle*coop_2pi
                call nest2ring(map%nside, i, ipix)
                call listpix%push(ipix)
                call listangle%push(real(rotate_angle))
             endif
          endif
       enddo
    case("Tmin", "Emin", "Bmin")  !!random orientation
       if(present(threshold))then
          if(threshold.lt.coop_healpix_max_threshold)then
             if(do_mask)then
                fcut = - threshold*sqrt(sum(dble(map%map(:,1))**2, mask%map(:,1).gt.0.5)/total_weight)
             else
                fcut = - threshold*sqrt(sum(dble(map%map(:,1))**2)/total_weight)
             endif
          else
             fcut = 1.e30
          endif
       else
          fcut = 1.e30
       endif
       do i=0, map%npix-1
          if(do_mask)then
             if(map%map(i,1).gt.fcut .or. mask%map(i, 1) .le. 0.5)cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if( all(map%map(list(1:nneigh),1) .gt. map%map(i,1)) .and. all(mask%map(list(1:nneigh), 1) .gt. 0.5 ) ) then
                call random_number(rotate_angle)
                rotate_angle = rotate_angle*coop_2pi
                call nest2ring(map%nside, i, ipix)
                call listpix%push(ipix)
                call listangle%push(real(rotate_angle))
             endif
          else
             if(map%map(i,1).gt.fcut )cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if(all(map%map(list(1:nneigh),1).gt.map%map(i,1)))then
                call random_number(rotate_angle)
                rotate_angle = rotate_angle*coop_2pi
                call nest2ring(map%nside, i, ipix)
                call listpix%push(ipix)
                call listangle%push(real(rotate_angle))
             endif
          endif
       enddo
    case("Tmax_QTUTOrient", "zetamax_qzuzOrient") !!oriented with QU, maxima of T
       if(map%nmaps .lt. 3) stop "For Tmax_QTUTOrient mode you need 3 maps"
       if(present(threshold))then
          if(threshold.lt.coop_healpix_max_threshold)then
             if(do_mask)then
                fcut = threshold*sqrt(sum(dble(map%map(:,1))**2,mask%map(:,1).gt.0.5)/total_weight)
             else
                fcut = threshold*sqrt(sum(dble(map%map(:,1))**2)/total_weight)
             endif
          else
             fcut = -1.e30
          endif
       else
          fcut = -1.e30
       endif
       do i=0, map%npix-1
          if(do_mask)then
             if( map%map(i,1) .lt. fcut .or. mask%map(i,1) .le. 0.5)cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if(all(map%map(list(1:nneigh),1) .lt. map%map(i,1)) .and. all(mask%map(list(1:nneigh),1) .gt. 0.5)) then
                rotate_angle = COOP_POLAR_ANGLE(map%map(i, 2), map%map(i, 3))/2.
                if(coop_random_unit().gt.0.5d0)rotate_angle=rotate_angle+coop_pi                
                call nest2ring(map%nside, i, ipix)
                call listpix%push(ipix)
                call listangle%push(real(rotate_angle))
             endif
          else
             if(map%map(i,1) .lt. fcut)cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if(all(map%map(list(1:nneigh),1) .lt. map%map(i,1)))then
                rotate_angle = COOP_POLAR_ANGLE(map%map(i, 2), map%map(i, 3))/2.
                if(coop_random_unit().gt.0.5d0)rotate_angle=rotate_angle+coop_pi                
                call nest2ring(map%nside, i, ipix)
                call listpix%push(ipix)
                call listangle%push(real(rotate_angle))
             endif
          endif
       enddo
    case("Tmin_QTUTOrient", "zetamin_qzuzOrient") !!oriented with QU, minima of T
       if(map%nmaps .lt. 3) stop "For Tmin_QTUTOrient mode you need 3 maps"
       if(present(threshold))then
          if(threshold.lt.coop_healpix_max_threshold)then
             if(do_mask)then
                fcut = -threshold*sqrt(sum(dble(map%map(:,1))**2, mask%map(:,1).gt.0.5)/total_weight)
             else
                fcut = -threshold*sqrt(sum(dble(map%map(:,1))**2)/total_weight)
             endif
          else
             fcut = 1.e30
          endif
       else
          fcut = 1.e30
       endif
       do i=0, map%npix-1
          if(do_mask)then
             if( map%map(i,1) .gt. fcut .or. mask%map(i,1) .le. 0.5) cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if( all( map%map(list(1:nneigh), 1) .gt. map%map(i,1) ) .and. all(mask%map(list(1:nneigh),1) .gt. 0.5) ) then
                rotate_angle = COOP_POLAR_ANGLE(map%map(i, 2), map%map(i, 3))/2.
                if(coop_random_unit().gt.0.5d0)rotate_angle=rotate_angle+coop_pi                                
                call nest2ring(map%nside, i, ipix)
                call listpix%push(ipix)
                call listangle%push(real(rotate_angle))
             endif
          else
             if(map%map(i,1) .gt. fcut)cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if(all(map%map(list(1:nneigh),1).gt. map%map(i,1)))then
                rotate_angle = COOP_POLAR_ANGLE(map%map(i, 2), map%map(i, 3))/2.
                if(coop_random_unit().gt.0.5d0)rotate_angle=rotate_angle+coop_pi                                
                call nest2ring(map%nside, i, ipix)
                call listpix%push(ipix)
                call listangle%push(real(rotate_angle))
             endif
          endif
       enddo
    case("Pmax", "PTmax", "PZmax")
       select case(map%nmaps)
       case(2:3)
          iq = map%nmaps - 1
          iu = iq + 1
       case default
          stop "For polarization you need to specify Q, U maps in the input file."
       end select
       if(present(threshold))then
          if(threshold.lt. coop_healpix_max_threshold)then
             if(do_mask)then
                fcut = threshold**2*(sum(dble(map%map(:,iq))**2 + dble(map%map(:,iu))**2, mask%map(:,1).gt.0.5)/total_weight)
             else
                fcut = threshold**2*(sum(dble(map%map(:,iq))**2 + dble(map%map(:,iu))**2)/total_weight)
             endif
          else
             fcut = -1.e30
          endif
       else
          fcut = -1.e30
       endif
       do i=0, map%npix-1
          if(do_mask)then
             if( map%map(i,iq)**2 + map%map(i,iu)**2 .lt. fcut  .or.  mask%map(i, 1) .le. 0.5 ) cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if( all( map%map(list(1:nneigh),iq)**2 + map%map(list(1:nneigh),iu)**2 .lt. map%map(i,iq)**2 + map%map(i, iu)**2 ) .and. all(mask%map(list(1:nneigh), 1) .gt. 0.5) ) then
                rotate_angle = COOP_POLAR_ANGLE(map%map(i, iq), map%map(i, iu))/2.
                if(coop_random_unit().gt.0.5d0)rotate_angle=rotate_angle+coop_pi                                
                call nest2ring(map%nside, i, ipix)
                call listpix%push(ipix)
                call listangle%push(real(rotate_angle))
             endif
          else
             if( map%map(i, iq)**2 + map%map(i, iu)**2 .lt. fcut) cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if( all( map%map(list(1:nneigh), iq)**2 + map%map(list(1:nneigh), iu)**2 .lt. map%map(i, iq)**2 + map%map(i, iu)**2 ) )then
                rotate_angle = COOP_POLAR_ANGLE(map%map(i, iq), map%map(i, iu))/2.
                if(coop_random_unit().gt.0.5d0)rotate_angle=rotate_angle+coop_pi                                
                call nest2ring(map%nside, i, ipix)
                call listpix%push(ipix)
                call listangle%push(real(rotate_angle))
             endif
          endif
       enddo
    case("Pmin", "PTmin", "PZmin")
       select case(map%nmaps)
       case(2:3)
          iq = map%nmaps - 1
          iu = iq + 1
       case default
          stop "For polarization you need to specify Q, U maps in the input file."
       end select
       if(present(threshold))then
          if(threshold.lt. coop_healpix_max_threshold)then          
             if(do_mask)then
                fcut = threshold**2*(sum(dble(map%map(:, iq))**2+dble(map%map(:, iu))**2, mask%map(:,1).gt.0.5)/total_weight)
             else
                fcut = threshold**2*(sum(dble(map%map(:, iq))**2+dble(map%map(:, iu))**2)/total_weight)
             endif
          else
             fcut = 1.d30
          endif
       else
          fcut = 1.d30
       endif
       do i=0, map%npix-1
          if(do_mask)then
             if(map%map(i,iq)**2 + map%map(i,iu)**2 .gt. fcut .or. mask%map(i, 1) .le. 0.5 ) cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if(all(map%map(list(1:nneigh), iq)**2 + map%map(list(1:nneigh), iu)**2 .gt. map%map(i, iq)**2 + map%map(i, iu)**2) .and. all(mask%map(list(1:nneigh), 1) .gt. 0.5) ) then
                rotate_angle = COOP_POLAR_ANGLE(map%map(i,iq), map%map(i,iu))/2.
                if(coop_random_unit().gt.0.5d0)rotate_angle=rotate_angle+coop_pi                                
                call nest2ring(map%nside, i, ipix)
                call listpix%push(ipix)
                call listangle%push(real(rotate_angle))
             endif
          else
             if(map%map(i, iq)**2 + map%map(i, iu)**2 .gt. fcut) cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if(all(map%map(list(1:nneigh), iq)**2 + map%map(list(1:nneigh), iu)**2 .gt. map%map(i, iq)**2 + map%map(i, iu)**2))then
                rotate_angle = COOP_POLAR_ANGLE(map%map(i, iq), map%map(i, iu))/2.
                if(coop_random_unit().gt.0.5d0)rotate_angle=rotate_angle+coop_pi                                
                call nest2ring(map%nside, i, ipix)
                call listpix%push(ipix)
                call listangle%push(real(rotate_angle))
             endif
          endif
       enddo
    case default
       write(*,*) trim(spot_type)
       stop "unknown spot type"
    end select
    call map%convert2ring()
    if(do_mask)    call mask%convert2ring()
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_maps_get_listpix

  subroutine coop_healpix_export_spots(map_file, spots_file, spot_type, threshold, threshold_pol, mask_file, fwhm)
    COOP_UNKNOWN_STRING map_file, spots_file, spot_type
    COOP_UNKNOWN_STRING, optional::mask_file
    COOP_REAL,optional::threshold, threshold_pol
    type(coop_file) fp
    type(coop_healpix_maps) mask, map
    COOP_REAL,optional::fwhm
    COOP_INT i
    type(coop_list_realarr) spots
    COOP_SINGLE arr(10)
    logical domask
    select case(trim(spot_type))
    case("Tmax_QTUTOrient", "Tmin_QTUTOrient", "PTmax", "PTmin", "PZmin", "PZmax", "zetamax_qzuzOrient", "zetamin_qzuzOrient")
       call map%read(trim(map_file), nmaps_wanted = 3)
    case default
       call map%read(trim(map_file))
    end select
    if(present(fwhm))then
       if(abs(fwhm) .gt. 1.d-3) &
            call map%smooth(fwhm)
    endif
    if(present(mask_file))then
       domask = (trim(mask_file).ne."")
    else
       domask = .false.
    endif
    if(domask) call mask%read(mask_file, nmaps_wanted = 1)
    if(present(threshold_pol))then
       if(present(threshold))then
          if(domask)then
             call map%get_spots(spots = spots,  spot_type = spot_type, threshold=threshold, threshold_pol=threshold_pol, mask = mask)
          else
             call map%get_spots(spots = spots,  spot_type = spot_type, threshold=threshold, threshold_pol=threshold_pol)
          endif
       else
          if(domask)then
             call map%get_spots(spots = spots,  spot_type = spot_type, threshold_pol=threshold_pol, mask = mask)
          else
             call map%get_spots(spots = spots,  spot_type = spot_type, threshold_pol=threshold_pol)
          endif
       endif       
    else
       if(present(threshold))then
          if(domask)then
             call map%get_spots(spots = spots,  spot_type = spot_type, threshold=threshold, mask = mask)
          else
             call map%get_spots(spots = spots,  spot_type = spot_type, threshold=threshold)
          endif
       else
          if(domask)then
             call map%get_spots(spots = spots,  spot_type = spot_type, mask = mask)
          else
             call map%get_spots(spots = spots,  spot_type = spot_type)
          endif
       endif
    endif
    call fp%open(trim(spots_file),"w")
    do i=1, spots%n
       call spots%get_element(i, arr)
       write(fp%unit, "(10G16.7)") arr(1:spots%dim)
    enddo
    call fp%close()
    call map%free
    call mask%free
    call spots%init
  end subroutine coop_healpix_export_spots

  subroutine coop_healpix_output_map(map, header, fname)
    real, dimension(:,:):: map
    character(LEN=80),dimension(:):: header
    COOP_UNKNOWN_STRING fname
#ifdef HAS_HEALPIX
    call coop_delete_file(trim(fname))
    call output_map(map, header, fname)
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_output_map
  
  subroutine coop_healpix_mask_map(mapfile, maskfile, output, index_list)
    type(coop_healpix_maps) map, mask
    COOP_INT i
    COOP_UNKNOWN_STRING mapfile, maskfile, output
    COOP_INT,dimension(:),optional::index_list
    call map%read(mapfile)
    call mask%read(maskfile,  nmaps_wanted = 1)
    if(mask%ordering .eq. COOP_RING)then
       call map%convert2ring()
    elseif(mask%ordering .eq. COOP_NESTED)then
       call map%convert2nested()
    else
       write(*,*) "coop_healpix_mask_map: mask file has unknown ordering."
       stop
    endif
    call coop_delete_file(output)
    if(present(index_list))then
       if(any(index_list .lt. 1 .or. index_list .gt. map%nmaps)) stop "coop_healpix_write_map: index out of range"
       do i=1, size(index_list)
          map%map(:,index_list(i)) = map%map(:, index_list(i))*mask%map(:,1)
       enddo
       call map%write(trim(coop_file_add_postfix(output,"_masked")), index_list)
       call map%write(output)
    else
       do i=1, map%nmaps
          map%map(:,i) = map%map(:, i)*mask%map(:,1)
       enddo
       call map%write(output)
    endif
    call map%free
    call mask%free
  end subroutine coop_healpix_mask_map


  subroutine coop_healpix_trim_maskfile(mask_file, smoothscale, output)
    COOP_UNKNOWN_STRING mask_file
    COOP_SINGLE smoothscale
    COOP_UNKNOWN_STRING, optional::output
    type(coop_healpix_maps) this
    call this%read(mask_file, nmaps_wanted = 1)
    call this%trim_mask( smoothscale)
    if(present(output))then
       call this%write(trim(output))
    else
       call this%write(trim(coop_file_add_postfix(trim(mask_file), "_smoothed")))
    endif
    call this%free
  end subroutine coop_healpix_trim_maskfile

  subroutine coop_healpix_trim_mask(this, smoothscale)
    real,parameter::nefolds = 4
    class(coop_healpix_maps) this
    type(coop_healpix_maps) hgs
    COOP_INT,dimension(:),allocatable::listpix
    COOP_INT i, j, nsteps
    COOP_SINGLE smoothscale
    nsteps = ceiling(smoothscale/sqrt(4.*coop_pi/this%npix))/2
    if(nsteps .le. 0 .or. nsteps .gt. 200)stop "coop_healpix_smooth_mask: invalid input of smoothscale"
    call this%convert2nested()
    this%mask_npix = count(this%map(:,1) .gt. 0.)
    allocate(this%mask_listpix(this%mask_npix))
    j = 0
    do i = 0, this%npix - 1
       if(this%map(i,1).gt. 0.)then
          j = j + 1
          this%mask_listpix(j) = i
       endif
    enddo
    select type(this)
    type is (coop_healpix_maps)
       hgs = this
    class default
       stop "the mask must be basic coop_healpix_maps type"
    end select
    do i=1, nsteps
       call coop_healpix_iterate_mask(this, hgs)
       call coop_healpix_iterate_mask(hgs, this)
    enddo
    call hgs%free
    call this%convert2ring()
  contains 
    subroutine coop_healpix_iterate_mask(this_from, this_to)  
      type(coop_healpix_maps) this_from, this_to
      COOP_INT list(8), nneigh
      COOP_INT i
#ifdef HAS_HEALPIX
      !$omp parallel do private(list, nneigh, i)
      do i = 1, this%mask_npix
         call neighbours_nest(this_from%nside, this%mask_listpix(i), list, nneigh)
         this_to%map(this_from%mask_listpix(i), 1) = min(minval(this_from%map(list(1:nneigh), 1)) , this_from%map(this_from%mask_listpix(i), 1))
      enddo
      !$omp end parallel do
#else
      stop "CANNOT FIND HEALPIX"
#endif
    end subroutine coop_healpix_iterate_mask
  end subroutine coop_healpix_trim_mask


  subroutine coop_healpix_smooth_maskfile(mask_file, smoothscale, output)
    COOP_UNKNOWN_STRING mask_file
    COOP_SINGLE smoothscale
    COOP_UNKNOWN_STRING, optional::output
    type(coop_healpix_maps) this
    call this%read(mask_file, nmaps_wanted = 1)
    call this%smooth_mask( smoothscale)
    if(present(output))then
       call this%write(trim(output))
    else
       call this%write(trim(coop_file_add_postfix(trim(mask_file), "_smoothed")))
    endif
    call this%free
  end subroutine coop_healpix_smooth_maskfile

  subroutine coop_healpix_smooth_mask(this, smoothscale)
    real,parameter::nefolds = 4
    class(coop_healpix_maps) this
    type(coop_healpix_maps) hgs
    COOP_INT i, j, nsteps
    COOP_SINGLE smoothscale, decay
    nsteps = ceiling(smoothscale*this%nside*nefolds/2.)
    if(nsteps .le. 0 .or. nsteps .gt. 200)stop "coop_healpix_smooth_mask: invalid input of smoothscale"
    decay = exp(-nefolds/nsteps/2.)
    call this%convert2nested()
    this%mask_npix = count(this%map(:,1) .lt. 1.)
    allocate(this%mask_listpix(this%mask_npix))
    j = 0
    do i = 0, this%npix - 1
       if(this%map(i,1).lt. 1.)then
          j = j + 1
          this%mask_listpix(j) = i
       endif
    enddo
    select type(this)
    type is (coop_healpix_maps)
       hgs = this
    class default
       stop "the mask must be basic coop_healpix_maps type"
    end select
    do i=1, nsteps
       call coop_healpix_iterate_mask(this, hgs, decay)
       call coop_healpix_iterate_mask(hgs, this, decay)
    enddo
    call hgs%free
    call this%convert2ring()
    deallocate(this%mask_listpix)
  contains 
    subroutine coop_healpix_iterate_mask(this_from, this_to, decay)  
      type(coop_healpix_maps) this_from, this_to
      COOP_INT list(8), nneigh
      COOP_INT i
      COOP_SINGLE decay
#ifdef HAS_HEALPIX
      !$omp parallel do private(list, nneigh, i)
      do i = 1, this%mask_npix
         call neighbours_nest(this_from%nside, this%mask_listpix(i), list, nneigh)
         this_to%map(this_from%mask_listpix(i), 1) = max(maxval(this_from%map(list(1:nneigh), 1)) * decay , this_from%map(this_from%mask_listpix(i), 1))
      enddo
      !$omp end parallel do
#else
      stop "CANNOT FIND HEALPIX"
#endif
    end subroutine coop_healpix_iterate_mask
  end subroutine coop_healpix_smooth_mask

  subroutine coop_healpix_mask_hemisphere(mask, theta, phi)
    COOP_REAL vec(3), theta, phi
    COOP_INT nlist
    type(coop_healpix_maps)::mask
    COOP_INT,dimension(:),allocatable::listpix
#ifdef HAS_HEALPIX
    allocate(listpix(0:mask%npix-1))
    call ang2vec(theta, phi, vec)
    call query_disc(mask%nside, vec, coop_pio2, listpix, nlist, nest = 0, inclusive = 0)
    mask%map(listpix(0:nlist-1),:) = 0.
    deallocate(listpix)
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_mask_hemisphere



  subroutine coop_healpix_inpainting(mode, map_file, mask_file, maskpol_file, output_freq, output_types, mask_smooth_scale)
    COOP_UNKNOWN_STRING map_file, mask_file, mode
    COOP_INT, optional::output_freq 
    COOP_UNKNOWN_STRING,optional:: maskpol_file, output_types
    COOP_INT,parameter:: total_steps = 1000, burnin = 30
    COOP_INT output_steps 
    COOP_INT step, naccept, weight
    logical accept
    type(coop_healpix_maps) map, simumap, mask, maskpol, mapmean, map_diffused
    COOP_SHORT_STRING::ot
    COOP_REAL prev_chisq
    COOP_REAL, optional:: mask_smooth_scale
    COOP_SINGLE::mss
    if(present(mask_smooth_scale))then
       mss = mask_smooth_scale
    else
       mss = 2.*coop_SI_degree
    endif
    if(present(output_freq))then
       output_steps = output_freq
    else
       output_steps = 50
    endif
    if(present(output_types))then
       ot = trim(output_types)
    else
       ot = "IQU"
    end if
    if(trim(mode) .eq. "I")then
       call map%read(map_file, nmaps_wanted = 1)
    else
       call map%read(map_file, nmaps_wanted = 3)
    endif
    mapmean = map   
    mapmean%map = 0
    call mask%read(mask_file, nmaps_wanted = 1)
    where(mask%map(:,1) .lt. 0.5)
       map%map(:,1) = 0.
    end where

    
    if(present(maskpol_file)) then
       call maskpol%read(maskpol_file, nmaps_wanted = 1)
       where(maskpol%map(:,1) .lt. 0.5)
          map%map(:,2) = 0.
          map%map(:,3) = 0.
       end where
       map%maskpol_npix = count(maskpol%map(:,1).lt. 0.5)
    else
       map%maskpol_npix = 0
    endif
    map_diffused = map
    call coop_healpix_diffuse_into_mask(map_diffused, mask, coop_healpix_diffuse_scale)
    call coop_healpix_smooth_mask(mask, mss)
    if(present(maskpol_file))then
       call coop_healpix_diffuse_into_mask(map_diffused, maskpol, coop_healpix_diffuse_scale, pol=.true.)
       call coop_healpix_smooth_mask(maskpol, mss)
    endif
    map%map = map_diffused%map
    write(*,*) "mask smoothed"
    !call mask%write(trim(coop_file_add_postfix(mask_file, "_smoothed")))
    map%mask_npix = count(mask%map(:,1).lt. 0.5)
    call coop_healpix_inpainting_init(map, simumap)

    write(*,*) "Inpaining workspace initialized."
    prev_chisq = map%chisq
    naccept = 0
    weight = 0
    step = 0
    do while(step .lt. total_steps)
       step = step + 1
       if(present(maskpol_file))then
          call coop_healpix_inpainting_step(accept, map, map_diffused, simumap, mask, maskpol)
       else
          call coop_healpix_inpainting_step(accept, map, map_diffused, simumap, mask)
       endif
       if(accept) naccept = naccept + 1
       if(naccept .gt. burnin/2 .and. naccept .le. burnin .and. accept)then
          map%Cl = simumap%Cl
       endif
       if(naccept .ge. burnin)then
          mapmean%map = mapmean%map + map%map
          if(weight.eq.0)then
             write(*,*) "initial trials done: now start sampling"
             step = 0
          endif
          weight = weight + 1
       endif
       if(naccept .lt. burnin)then
          write(*,*) "Trial #"//trim(coop_num2str(step))//" before burn-in. Estimated progress: "//trim(coop_num2str(100.*naccept/real(burnin), "(F10.1)"))//"%"
       endif
       write(*, "(I6, A)") step, " accept = "//trim(coop_num2str(accept))//" temperature = "//trim(coop_num2str(map%mcmc_temperature,"(F10.2)"))//" chisq = "//trim(coop_num2str(prev_chisq, "(G14.3)"))//" --> "//trim(coop_num2str(map%chisq, "(G14.3)"))

       prev_chisq = map%chisq

       
       if(mod(step, output_steps) .eq.0 .and. naccept.ge.burnin)then
          call map%write(trim(coop_file_add_postfix(trim(map_file), "_inp"//trim(coop_ndigits(step, 4)))))         
       endif
       if(mod(step, output_steps*10) .eq. 0 .and. naccept.gt.burnin)then !!less mean maps to save disk space
          simumap%map =  mapmean%map/weight
          call simumap%write(trim(coop_file_add_postfix(trim(map_file), "_inp_mean"//trim(coop_ndigits(step, 4)))))
       endif
    enddo
    call map%free()
    call mapmean%free()
    call simumap%free()
    call mask%free()
    call maskpol%free()
  end subroutine coop_healpix_inpainting

  subroutine coop_healpix_inpainting_init(map, simumap)
    type(coop_healpix_maps) map, simumap,tmpmap
    call map%map2alm()
    map%cl(0:1, :) = 0.
    map%cl(:, coop_healpix_index_TT) = max(map%cl(:, coop_healpix_index_TT), 1.e-20)*(map%npix/real(map%npix - map%mask_npix))
    if(map%nmaps.ge.3)then
       map%cl(:, coop_healpix_index_EE) = max(map%cl(:, coop_healpix_index_EE), 1.e-20)*(map%npix/real(map%npix-map%maskpol_npix))
       map%cl(:, coop_healpix_index_TE) = map%cl(:, coop_healpix_index_TE)* sqrt( (map%npix/real(map%npix - map%mask_npix)) * (map%npix/real(map%npix-map%maskpol_npix)) ) * (1.d0 - 1.d-6)  !!make sure TE^2 < TT * EE
       map%cl(:, coop_healpix_index_BB) = max(map%cl(:, coop_healpix_index_BB), 1.e-12)*(map%npix/real(map%npix-map%maskpol_npix))
       map%cl(:, coop_healpix_index_EB) = 0. 
       map%Cl(:, coop_healpix_index_TB) = 0.
    endif
    map%Cl(0:1,:) = 0.
    map%mcmc_temperature = 20.     !!start with a high temperature
    coop_healpix_inpainting_lowl = 5
    simumap = map
    call coop_healpix_inpainting_get_chisq(map, simumap)

    map%chisq = simumap%chisq
  end subroutine coop_healpix_inpainting_init
  
  subroutine coop_healpix_inpainting_step(accept, map, map_diffused, simumap, mask, maskpol)
    type(coop_healpix_maps) map, simumap, mask, map_diffused
    type(coop_healpix_maps),optional::maskpol
    logical accept
    COOP_INT i
    simumap%Cl = map%Cl
    call simumap%simulate()
    simumap%map(:,1) =  map%map(:,1) * mask%map(:,1) + simumap%map(:,1) * sqrt(1.-mask%map(:,1)**2) + map_diffused%map(:,1)*(1. - mask%map(:,1))
    if(map%nmaps.eq.3)then
       if(present(maskpol))then
          simumap%map(:,2) =  map%map(:,2) * maskpol%map(:,1) + simumap%map(:,2) * sqrt(1.-maskpol%map(:,1)**2) + map_diffused%map(:,2)*(1. - maskpol%map(:,1))
          simumap%map(:,3) =  map%map(:,3) * maskpol%map(:,1) + simumap%map(:,3) * sqrt(1.-maskpol%map(:,1)**2) + map_diffused%map(:,3)*(1. - maskpol%map(:,1))
       else
          stop "coop_healpix_inpainting_step: need polarization mask"
       endif
    endif
    call coop_healpix_inpainting_accept_reject(accept, map, simumap)
  end subroutine coop_healpix_inpainting_step


  subroutine coop_healpix_inpainting_get_chisq(map, simumap)
    type(coop_healpix_maps) map, simumap
    COOP_INT l
    COOP_SINGLE chisq
    call simumap%map2alm()
    if(map%nmaps .eq. 1)then
       chisq = 0.
       !$omp parallel do reduction(+:chisq)
       do l = 2, coop_healpix_inpainting_lowl 
          chisq = chisq + (2*l+1)*(simumap%Cl(l,1)/map%Cl(l,1) - log(simumap%Cl(l,1)/map%Cl(l,1))-1.d0)
       enddo
       !$omp end parallel do
       simumap%chisq = chisq 

    else
       chisq = 0.
       !$omp parallel do reduction(+:chisq)
       do l = 2, coop_healpix_inpainting_lowl 
          chisq = chisq + (2*l+1)*(simumap%Cl(l,coop_healpix_index_BB)/map%Cl(l,coop_healpix_index_BB) - log(simumap%Cl(l,coop_healpix_index_BB)/map%Cl(l,coop_healpix_index_BB))-3.d0 - log((simumap%Cl(l, coop_healpix_index_TT)*simumap%Cl(l, coop_healpix_index_EE) - simumap%Cl(l, coop_healpix_index_TE)**2)/(map%Cl(l, coop_healpix_index_TT)*map%Cl(l, coop_healpix_index_EE) - map%Cl(l, coop_healpix_index_TE)**2)) + (map%Cl(l, coop_healpix_index_TT)*simumap%Cl(l, coop_healpix_index_EE) + simumap%Cl(l, coop_healpix_index_TT)*map%Cl(l, coop_healpix_index_EE) - 2.*map%Cl(l, coop_healpix_index_TE)*simumap%Cl(l, coop_healpix_index_TE) )/(map%Cl(l, coop_healpix_index_TT)*map%Cl(l, coop_healpix_index_EE) - map%Cl(l, coop_healpix_index_TE)**2))
       enddo
       !$omp end parallel do
       simumap%chisq = chisq 
    endif
  end subroutine coop_healpix_inpainting_get_chisq

  subroutine coop_healpix_inpainting_accept_reject(accept, map, simumap)
    type(coop_healpix_maps) map, simumap
    logical accept
    call coop_healpix_inpainting_get_chisq(map, simumap)
    if(simumap%chisq - map%chisq .lt. coop_random_exp()*2*map%mcmc_temperature)then
       accept = .true.
       map%chisq = simumap%chisq
       map%map = simumap%map

       map%alm = simumap%alm
       map%mcmc_temperature = max(1., map%mcmc_temperature*0.95)
       coop_healpix_inpainting_lowl = min(coop_inpainting_lowl_max, coop_healpix_inpainting_lowl + 1)
    else
       accept = .false.
       map%mcmc_temperature = map%mcmc_temperature*1.01
       coop_healpix_inpainting_lowl = max(coop_inpainting_lowl_min, coop_healpix_inpainting_lowl - 1)
    endif
  end subroutine coop_healpix_inpainting_accept_reject


  subroutine coop_healpix_split(filename)
    COOP_UNKNOWN_STRING:: filename
    type(coop_healpix_maps) this
    COOP_INT i
    call this%read(filename)
    do i=1, this%nmaps
       call this%write(trim(coop_file_add_postfix(filename, "_submap"//trim(coop_ndigits(i, 3)))), (/ i /) )
    enddo
    call this%free()
  end subroutine coop_healpix_split

  subroutine coop_healpix_plot_spots(spotsfile, mapfile)
    COOP_UNKNOWN_STRING spotsfile, mapfile
    type(coop_file) fp
    COOP_REAL theta, phi, angle_rotate
    COOP_INT pix
    type(coop_healpix_maps) this
#ifdef HAS_HEALPIX
    call this%init(64, 1, (/ 0 /) )
    this%map = 0
    call fp%open(trim(spotsfile), "r")
    do
       read(fp%unit, *, END=100, ERR=100) theta, phi, angle_rotate
       call ang2pix_ring(this%nside, theta, phi, pix)
       this%map(pix, 1) = 1.d0
    enddo
100 call fp%close()
    call this%write(trim(mapfile))
    call this%free
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_plot_spots

  subroutine coop_healpix_lb2ang(l_deg, b_deg, theta, phi)
    COOP_REAL l_deg, b_deg
    COOP_REAL theta, phi
    theta = coop_pio2 - b_deg*coop_SI_degree
    phi = l_deg * coop_SI_degree
  end subroutine coop_healpix_lb2ang


  subroutine coop_healpix_ang2lb(theta, phi, l_deg, b_deg)
    COOP_REAL l_deg, b_deg
    COOP_REAL theta, phi
    b_deg = (coop_pio2 - theta)/coop_SI_degree
    l_deg = phi/coop_SI_degree
    do while(b_deg .gt. 90.d0)
       b_deg = b_deg - 90.d0
    enddo
    do while(b_deg .lt. -90.d0)
       b_deg = b_deg + 90.d0
    enddo
    do while(l_deg .gt. 360.d0)
       l_deg = l_deg - 360.d0
    enddo
    do while(l_deg .lt. 0.d0)
       l_deg = l_deg + 360.d0
    enddo
  end subroutine coop_healpix_ang2lb


  subroutine coop_healpix_flip_mask(mask, flip)
    type(coop_healpix_maps)mask, flip
    COOP_INT pix, conjpix
    COOP_REAL theta, phi
#ifdef HAS_HEALPIX
    if(allocated(flip%map))then
       if(flip%nside .ne. mask%nside)then
          call flip%free()
          call flip%init(nside = mask%nside, nmaps = 1, spin = (/ 0 /) )
       endif
    else
       call flip%free()
       call flip%init(nside = mask%nside, nmaps = 1, spin = (/ 0 /) )       
    endif
    call mask%convert2ring()
    !$omp parallel do private(pix, conjpix, theta, phi)
    do pix=0, mask%npix-1
       call pix2ang_ring(mask%nside, pix, theta, phi)
       theta = coop_pi-theta
       phi = coop_pi +  phi
       call ang2pix_ring(mask%nside, theta, phi, conjpix)
       flip%map(pix, 1) = mask%map(conjpix, 1)
    enddo
    !$omp end parallel do
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_flip_mask

  subroutine coop_healpix_maps_mask(map, mask, polmask)
    class(coop_healpix_maps) map
    type(coop_healpix_maps) mask
    type(coop_healpix_maps),optional::polmask
    if(map%nside  .ne. mask%nside) stop "mask: nside must be the same as the map"
    map%map(:,1) = map%map(:,1)*mask%map(:,1)
    if(map%nmaps.ge.3 .and. present(polmask))then
       if(map%nside  .ne. polmask%nside) stop "pol mask: nside must be the same as the map"
       map%map(:,2) = map%map(:,2)*polmask%map(:,1)
       map%map(:,3) = map%map(:,3)*polmask%map(:,1)
    endif
  end subroutine coop_healpix_maps_mask

  subroutine coop_healpix_maps_get_fullCls(map, mask, polmask)
    COOP_INT,parameter::lrange = 20
    class(coop_healpix_maps) map
    type(coop_healpix_maps) mask, mapcopy
    type(coop_healpix_maps),optional::polmask
    COOP_REAL, dimension(:),allocatable::ifsky, polfsky
    COOP_REAL:: w(-lrange:lrange)
    COOP_INT i, j, ncls, l, isim
    if(mask%nside .ne. map%nside) stop "nside must agree for get_fullCl"
    if(present(polmask) .and. map%nmaps .eq. 3 .and. map%iq .eq. 2)then
       if(polmask%nside .ne. map%nside) stop "nside must agree for get_fullCl"
       call map%mask(mask, polmask)
       ncls = 6
    else
       call map%mask(mask)
       ncls = 1
    endif
    call map%map2alm()

    allocate(ifsky(0:map%lmax))
    ifsky = count(mask%map(:,1).gt.0.)/dble(mask%npix)
    if(ncls .eq. 6)then
       allocate(polfsky(0:map%lmax))
       polfsky = count(polmask%map(:, 1).gt.0.)/dble(polmask%npix)
    endif

    select type(map)
    type is (coop_healpix_maps)
       mapcopy = map
    class default
       stop "get_fullCls: can only process healpix_maps structure"
    end select
    do l=-lrange, lrange
       w(l) = exp(-(l*2.d0/lrange)**2)
    enddo
    w = w/sum(w)
    call do_cl_map()
    do isim = 1, 3
       call mapcopy%simulate()
       if(present(polmask) .and. map%nmaps .eq. 3 .and. map%iq .eq. 2)then
          call mapcopy%mask(mask, polmask)
       else
          call mapcopy%mask(mask)
       endif
       call mapcopy%map2alm()
       do l= lrange, map%lmax - lrange
          ifsky(l) = ifsky(l)/sum(map%cl(l-lrange:l+lrange, coop_healpix_index_TT)*w)*sum(mapcopy%cl(l-lrange:l+lrange, coop_healpix_index_TT)*w)
       enddo
       do l=0, lrange-1
          ifsky(l) = ifsky(l)/sum(map%cl(0:l+lrange, coop_healpix_index_TT)*w(-l:lrange))*sum(mapcopy%cl(0:l+lrange, coop_healpix_index_TT)*w(-l:lrange))
       enddo
       do l = map%lmax - lrange + 1, map%lmax
          ifsky(l) = ifsky(l)/sum(map%cl(l-lrange:map%lmax, coop_healpix_index_TT)*w(-lrange:map%lmax-l))*sum(mapcopy%cl(l-lrange:map%lmax, coop_healpix_index_TT)*w(-lrange:map%lmax-l))
       enddo
       if(ncls.eq.6)then
          do l= lrange, map%lmax - lrange
             polfsky(l) = polfsky(l)/sum(map%cl(l-lrange:l+lrange, coop_healpix_index_EE)*w)*sum(mapcopy%cl(l-lrange:l+lrange, coop_healpix_index_EE)*w)
          enddo
          do l=0, lrange-1
             polfsky(l) = polfsky(l)/sum(map%cl(0:l+lrange, coop_healpix_index_EE)*w(-l:lrange))*sum(mapcopy%cl(0:l+lrange, coop_healpix_index_EE)*w(-l:lrange))
          enddo
          do l = map%lmax - lrange + 1, map%lmax
             polfsky(l) = polfsky(l)/sum(map%cl(l-lrange:map%lmax, coop_healpix_index_EE)*w(-lrange:map%lmax-l))*sum(mapcopy%cl(l-lrange:map%lmax, coop_healpix_index_EE)*w(-lrange:map%lmax-l))
          enddo
       endif
       call do_cl_map()
    enddo
    map%cl = mapcopy%cl 
    call mapcopy%free()
    deallocate(ifsky)
    if(ncls.eq.6)deallocate(polfsky)

    contains

      subroutine do_cl_map()
       mapcopy%cl(:,coop_healpix_index_TT) = map%cl(:,coop_healpix_index_TT)/ifsky
       select case(mapcopy%nmaps)
       case(1)  !!I
       case(3)  !!I, Q, U
          if(.not. present(polmask))stop "for IQU mapcopys you need to specify polarization mask"
          if( map%iq .ne. 2) stop "get_fullcls only supports IQU map"
          mapcopy%cl(:,coop_healpix_index_EE) = map%cl(:,coop_healpix_index_EE)/polfsky
          mapcopy%cl(:,coop_healpix_index_BB) = map%cl(:,coop_healpix_index_BB)/polfsky
          mapcopy%cl(:,coop_healpix_index_TE) = map%cl(:,coop_healpix_index_TE)/sqrt(ifsky*polfsky)*(1.d0-1.d-6)
          mapcopy%cl(:,coop_healpix_index_TB) = map%cl(:,coop_healpix_index_TB)/sqrt(ifsky*polfsky) * (1.d0 - 1.d-6)
          mapcopy%cl(:,coop_healpix_index_EB) = map%cl(:,coop_healpix_index_EB)/polfsky
       case default
          stop "get_fullCl only supports I or IQU maps"
       end select
     end subroutine do_cl_map
  end subroutine coop_healpix_maps_get_fullCls

  subroutine coop_healpix_maps_udgrade(map, nside)
    class(coop_healpix_maps)::map
    real, dimension(:,:),allocatable::newmap
    COOP_INT nside, npix
#ifdef HAS_HEALPIX
    npix = nside2npix(nside)
    allocate(newmap(0:npix-1, map%nmaps))
    if(map%ordering .eq. COOP_NESTED)then
       call udgrade_nest(map%map, map%nside, newmap, nside, 0., .true.)
    else
       call udgrade_ring(map%map, map%nside, newmap, nside, 0., .true.)
    endif
    deallocate(map%map)
    map%nside = nside
    map%npix = npix
    allocate(map%map(0:npix-1, map%nmaps))
    map%map = newmap
    map%header  = ""
    call write_minimal_header(map%header,dtype = 'MAP', nside=nside, order = map%ordering, creator='Zhiqi Huang', version = 'COOP', units='muK', polar=any(map%spin.eq.2) )
    
    deallocate(newmap)
#else
    stop "HEALPIX library is missing."
#endif    
  end subroutine coop_healpix_maps_udgrade



  subroutine coop_healpix_maps_t2zeta(this, fwhm_arcmin)
    class(coop_healpix_maps)::this
    COOP_INT::l
    type(coop_cosmology_firstorder)::fod
    COOP_REAL,dimension(:,:),allocatable::Cls, Cls_lensed
    COOP_REAL::fwhm_arcmin, norm
#if DO_ZETA_TRANS
    call this%map2alm(index_list = (/ 1 /) )
    
    call fod%Set_Planck_bestfit()
    call fod%compute_source(0)
    allocate(Cls(coop_num_cls, 2:this%lmax), Cls_lensed(coop_num_cls, 2:this%lmax))
    call fod%source(0)%get_All_Cls(2, this%lmax, Cls)
    call coop_get_lensing_Cls(2, this%lmax, Cls, Cls_lensed)
    
    norm = 1.d5 / 2.726d6
    this%alm(0:1, :, 1) = 0.
    !$omp parallel do
    do l = 2, this%lmax
       this%alm(l, :, 1) = this%alm(l,:,1)*(Cls(coop_index_ClTzeta,l)/(Cls(coop_index_ClTT, l) + abs(Cls_lensed(coop_index_ClTT,l))+ Noise(l))/sqrt(beam(l)) * norm )
    enddo
    !$omp end parallel do

    call this%alm2map( index_list = (/ 1 /) )
    this%alm_done(1) = .false.
    deallocate(Cls, Cls_lensed)
  contains

    function beam(l) result(bl)
      COOP_INT l
      COOP_REAL bl
      bl = 1./(1.+l*(l+1)*(fwhm_arcmin*coop_SI_arcmin*coop_sigma_by_fwhm)**2)
    end function beam

    function Noise(l) result(Nl)
      COOP_REAL,parameter::C(6) = (/ 0.000209589 , &
           -2.60632e-6, &
           0.000102449 , & 
           -3.45531e-5, &
           -1.90758e-5, &
           8.00965e-6 /)
      COOP_REAL, parameter:: alpha = 0.75, &
           lmax = 4000
      COOP_INT l
      COOP_REAL Nl, x, fx, t
      x =( (l*(l+1.0))/(lmax*(lmax+1.0)) )**0.25
      t = 2.d0*x - 1.d0
      call coop_get_cheb_value(6, c, t, fx)
      Nl = fx / x ** alpha / (2.726e6)**2 / beam(l)
    end function Noise
#else
    stop "To use t2zeta you have to turn on DO_ZETA_TRANS in include/constants.h"
#endif
  end subroutine coop_healpix_maps_t2zeta



end module coop_healpix_mod


