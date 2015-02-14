module coop_healpix_mod
!!I always assume ring order
  use coop_wrapper_firstorder
  use coop_sphere_mod
  use coop_stacking_mod
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

  !!some stupid conventions...
  !!the direction of headless vector rotate by pi/2 if you use IAU convention
  logical::coop_healpix_IAU_headless_vector = .true.
  COOP_REAL::coop_healpix_QrUrSign = -1.d0  !!WMAP7 invented this...
  
  COOP_SINGLE::coop_healpix_patch_default_figure_width = 4.
  COOP_SINGLE::coop_healpix_patch_default_figure_height = 3.
  logical::coop_healpix_patch_default_want_caption = .false.
  logical::coop_healpix_patch_default_want_label = .false.
  logical::coop_healpix_patch_default_want_arrow = .false.
  logical::coop_healpix_warning = .true.
  
  public::coop_healpix_maps, coop_healpix_disc, coop_healpix_patch, coop_healpix_split,  coop_healpix_inpainting, coop_healpix_smooth_maskfile, coop_healpix_output_map, coop_healpix_smooth_mapfile, coop_healpix_patch_get_fr0, coop_healpix_fetch_patch, coop_healpix_mask_tol,  coop_healpix_mask_hemisphere, coop_healpix_index_TT,  coop_healpix_index_EE,  coop_healpix_index_BB,  coop_healpix_index_TE,  coop_healpix_index_TB,  coop_healpix_index_EB, coop_healpix_flip_mask, coop_healpix_diffuse_into_mask, coop_healpix_alm_check_done, coop_healpix_want_cls, coop_healpix_default_lmax, coop_planck_TNoise, coop_planck_ENoise, coop_Planck_BNoise, coop_highpass_filter, coop_lowpass_filter, coop_gaussian_filter,coop_healpix_latitude_cut_mask, coop_healpix_trim_maskfile, coop_healpix_IAU_headless_vector,  coop_healpix_latitude_cut_smoothmask, coop_healpix_spot_select_mask, coop_healpix_spot_cut_mask, coop_healpix_merge_masks, coop_healpix_patch_default_figure_width, coop_healpix_patch_default_figure_height, coop_healpix_patch_default_want_caption, coop_healpix_patch_default_want_label, coop_healpix_patch_default_want_arrow, coop_healpix_warning 
  

  logical::coop_healpix_alm_check_done = .false.
  logical::coop_healpix_want_cls = .true.
  COOP_REAL,parameter::coop_healpix_zeta_normalization = 1.d-5

  COOP_INT,parameter::dlc = kind( (1.d0,1.d0) )
  COOP_INT,parameter::coop_inpainting_lowl_max = 20
  COOP_INT,parameter::coop_inpainting_lowl_min = 5

  COOP_INT, parameter::coop_healpix_default_lmax=2500
  COOP_REAL,parameter::coop_healpix_mask_tol = 0.5  !!default mask tolerance
  COOP_INT::coop_healpix_inpainting_lowl=5
  COOP_REAL,parameter::coop_healpix_diffuse_scale = 10.d0*coop_SI_arcmin
  COOP_INT,parameter::coop_healpix_index_TT = 1
  COOP_INT,parameter::coop_healpix_index_EE = 2
  COOP_INT,parameter::coop_healpix_index_BB = 3
  COOP_INT,parameter::coop_healpix_index_TE = 4
  COOP_INT,parameter::coop_healpix_index_EB = 5
  COOP_INT,parameter::coop_healpix_index_TB = 6

  type, extends(coop_sphere_disc):: coop_healpix_disc
     COOP_INT::nside = 0
     COOP_INT::center
     COOP_INT::ordering = COOP_RING
   contains
     procedure :: pix2ang => coop_healpix_disc_pix2ang
     procedure :: ang2pix => coop_healpix_disc_ang2pix
     procedure :: pix2xy => coop_healpix_disc_pix2xy
     procedure :: xy2pix => coop_healpix_disc_xy2pix
  end type coop_healpix_disc

  type coop_healpix_maps
     COOP_INT::nside = 0
     COOP_INT::nmaps = 0
     COOP_INT:: npix = 0
     COOP_INT::ordering = COOP_RING
     COOP_INT::lmax = -1
     COOP_INT::iq = 0
     COOP_INT::iu = 0
     COOP_INT::mask_npix= 0
     COOP_INT::maskpol_npix = 0
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
     procedure :: ang2pix => coop_healpix_maps_ang2pix
     procedure :: pix2ang => coop_healpix_maps_pix2ang
     procedure :: vec2pix => coop_healpix_maps_vec2pix
     procedure :: pix2vec => coop_healpix_maps_pix2vec
     procedure :: query_disc => coop_healpix_maps_query_disc
     procedure :: query_strip => coop_healpix_maps_query_strip
     procedure :: get_disc => coop_healpix_maps_get_disc
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
     procedure :: iqu2TQTUT1 => coop_healpix_maps_iqu2TQTUT1     
     procedure :: smooth => coop_healpix_maps_smooth
     procedure :: smooth_with_window => coop_healpix_maps_smooth_with_window
     procedure :: smooth_mask => coop_healpix_smooth_mask
     procedure :: T2zeta => coop_healpix_maps_t2zeta
     procedure :: TE2zeta => coop_healpix_maps_te2zeta
     procedure :: E2zeta => coop_healpix_maps_E2zeta
     procedure :: trim_mask => coop_healpix_trim_mask
     procedure :: convert2nested => coop_healpix_convert_to_nested
     procedure :: convert2ring => coop_healpix_convert_to_ring
     procedure :: filter_alm =>  coop_healpix_filter_alm
     !!stacking stuff
     procedure :: get_peaks => coop_healpix_maps_get_peaks
     procedure :: mark_peaks => coop_healpix_maps_mark_peaks
     procedure :: stack_on_peaks  =>     coop_healpix_maps_stack_on_peaks     
  end type coop_healpix_maps

  type coop_healpix_patch
     COOP_STRING::caption=""
     type(coop_to_be_stacked):: tbs
     COOP_SHORT_STRING::color_table="Rainbow"
     COOP_INT::n = 0
     COOP_INT::mmax = 0
     COOP_INT::nmaps = 0
     COOP_INT::npix = 0
     COOP_INT::nstack_raw = 0
     COOP_REAL::dr
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


  subroutine coop_healpix_maps_ang2pix(this, theta, phi, pix)
    class(coop_healpix_maps)::this
    COOP_REAL theta, phi
    COOP_INT pix
#ifdef HAS_HEALPIX    
    if(this%ordering .eq. COOP_RING)then
       call ang2pix_ring(this%nside, theta, phi, pix)
    else
       call ang2pix_nest(this%nside, theta, phi, pix)       
    endif
#else
    stop "You need to compile with Healpix to use ang2pix" 
#endif       
  end subroutine coop_healpix_maps_ang2pix

  subroutine coop_healpix_maps_pix2ang(this, pix, theta, phi)
    class(coop_healpix_maps)::this
    COOP_REAL theta, phi
    COOP_INT pix
#ifdef HAS_HEALPIX    
    if(this%ordering .eq. COOP_RING)then
       call pix2ang_ring(this%nside, pix, theta, phi)
    else
       call pix2ang_nest(this%nside, pix, theta, phi)
    endif
#else
    stop "You need to compile with Healpix to use pix2ang" 
#endif       
  end subroutine coop_healpix_maps_pix2ang
  

    subroutine coop_healpix_maps_vec2pix(this, vec, pix)
    class(coop_healpix_maps)::this
    COOP_REAL vec(3)
    COOP_INT pix
#ifdef HAS_HEALPIX    
    if(this%ordering .eq. COOP_RING)then
       call vec2pix_ring(this%nside, vec, pix)
    else
       call vec2pix_nest(this%nside, vec, pix)       
    endif
#else
    stop "You need to compile with Healpix to use ang2pix" 
#endif       
  end subroutine coop_healpix_maps_vec2pix

  subroutine coop_healpix_maps_pix2vec(this, pix, vec)
    class(coop_healpix_maps)::this
    COOP_REAL vec(3)
    COOP_INT pix
#ifdef HAS_HEALPIX    
    if(this%ordering .eq. COOP_RING)then
       call pix2vec_ring(this%nside, pix, vec)
    else
       call pix2vec_nest(this%nside, pix, vec)
    endif
#else
    stop "You need to compile with Healpix to use pix2ang" 
#endif       
  end subroutine coop_healpix_maps_pix2vec
  
  subroutine coop_healpix_spot_cut_mask(nside, l_deg, b_deg, r_deg, filename)
    COOP_INT::nside
    COOP_REAL::l_deg, b_deg, r_deg
    COOP_UNKNOWN_STRING::filename
    type(coop_healpix_maps)::mask
    COOP_INT::listpix(0:nside**2*12-1), nlist, pix
    COOP_REAL::theta, phi
#ifdef HAS_HEALPIX
    call coop_healpix_lb2ang(l_deg, b_deg, theta, phi)
    call mask%init(nside = nside, nmaps = 1, spin = (/ 0 /))
    mask%map(:,1) = 1.
    call mask%ang2pix(theta, phi, pix)
    call mask%query_disc(pix, r_deg*coop_SI_degree, listpix, nlist)
    mask%map(listpix(0:nlist-1), 1) = 0.
    call mask%write(trim(filename))
    call mask%free()
#else
    stop "cannot find HEALPIX"
#endif    
  end subroutine coop_healpix_spot_cut_mask


  subroutine coop_healpix_spot_select_mask(nside, l_deg, b_deg, r_deg, filename)
    COOP_INT::nside
    COOP_REAL::l_deg, b_deg, r_deg
    COOP_UNKNOWN_STRING::filename
    type(coop_healpix_maps)::mask
    COOP_INT::listpix(0:nside**2*12-1), nlist, pix
    COOP_REAL::theta, phi
#ifdef HAS_HEALPIX
    call coop_healpix_lb2ang(l_deg, b_deg, theta, phi)
    call mask%init(nside = nside, nmaps = 1, spin = (/ 0 /))
    mask%map(:,1) = 0.
    call mask%ang2pix(theta, phi, pix)
    call mask%query_disc(pix, r_deg*coop_SI_degree, listpix, nlist)
    mask%map(listpix(0:nlist-1), 1) = 1.
    call mask%write(trim(filename))
    call mask%free()
#else
    stop "cannot find HEALPIX"
#endif    
  end subroutine coop_healpix_spot_select_mask
  

  subroutine coop_healpix_merge_masks(mask1, mask2, maskout)
    COOP_UNKNOWN_STRING::mask1,mask2,maskout
    type(coop_healpix_maps)::m1,m2
    call m1%read(mask1, nmaps_wanted=1, spin = (/ 0 /))
    call m2%read(mask2, nmaps_wanted=1, spin = (/ 0 /))
    m1%map = m1%map * m2%map
    call m1%write(maskout)
  end subroutine coop_healpix_merge_masks


  subroutine coop_healpix_latitude_cut_mask(nside, latitude_degree, filename)
    COOP_REAL::latitude_degree
    COOP_UNKNOWN_STRING::filename
    COOP_INT:: nside
    type(coop_healpix_maps)::mask
    COOP_INT::listpix(0:nside**2*12-1), nlist
#ifdef HAS_HEALPIX    
    call mask%init(nside = nside, nmaps = 1, spin = (/ 0 /))
    mask%map(:,1) = 1.
    call query_strip(nside, coop_pio2 - latitude_degree*coop_SI_degree, coop_pio2 + latitude_degree*coop_SI_degree, listpix, nlist, nest = 0, inclusive = 1)
    mask%map(listpix(0:nlist-1), 1) = 0.
    call mask%write(trim(filename))
    call mask%free()
#else
    stop "cannot find HEALPIX"
#endif    
  end subroutine coop_healpix_latitude_cut_mask

  subroutine coop_healpix_latitude_cut_smoothmask(nside, latitude_degree, depth, filename)
    COOP_REAL::latitude_degree, depth
    COOP_UNKNOWN_STRING::filename
    type(coop_healpix_maps)::mask
    COOP_INT nside
    COOP_INT i
    COOP_REAL theta, phi, lat, dep
    COOP_INT::listpix(0:nside**2*12-1), nlist
#ifdef HAS_HEALPIX    
    lat = latitude_degree*coop_SI_degree
    dep = depth*coop_SI_degree
    call mask%init(nside = nside, nmaps = 1, spin = (/ 0 /))
    mask%map(:,1) = 1.
    call query_strip(nside, coop_pio2 - lat, coop_pio2 + lat, listpix, nlist, nest = 0, inclusive = 1)
    do i = 0, nlist-1
       call pix2ang_ring(nside, listpix(i), theta, phi)
       mask%map(listpix(i),1) = exp(- ((lat - abs(theta - coop_pio2))/dep)**2)
    enddo
    call mask%write(trim(filename))
    call mask%free()
#else
    stop "you need to install healpix"
#endif    
  end subroutine coop_healpix_latitude_cut_smoothmask
  

  subroutine coop_healpix_diffuse_into_mask(this, mask, smoothscale, pol)  !!lambda<1
   ! real,parameter::nefolds = 2.
    type(coop_healpix_maps) mask, this
    type(coop_healpix_maps) masknew, hgs, maskcopy
    COOP_INT i, j, nsteps, nmaps, istart, iend
    COOP_REAL smoothscale, decay
    logical, optional::pol
#ifdef HAS_HEALPIX    
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
#else
    stop "You need to install healpix"
#endif    
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

  subroutine coop_healpix_patch_plot(this, imap, output, use_degree)
    COOP_INT::bgrids
    class(coop_healpix_patch)::this
    COOP_INT imap
    COOP_UNKNOWN_STRING::output
    type(coop_asy)::fig
    COOP_INT nb, i, j, k, ns, iq, iu
    COOP_REAL  xc, yc,  norm, r, theta, minz, maxz
    COOP_REAL,dimension(:),allocatable::xstart, xend, ystart, yend
    logical,optional::use_degree
    logical use_rad
    COOP_SHORT_STRING::xlabel, ylabel
    if(imap .le. 0 .or. imap .gt. this%nmaps) stop "coop_healpix_patch_plot: imap overflow"    
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
    if(coop_healpix_patch_default_want_caption)then
       call fig%init(caption = trim(this%caption), xlabel =trim(xlabel), ylabel =trim(ylabel), width = coop_healpix_patch_default_figure_width, height = coop_healpix_patch_default_figure_height, xmin = -real(this%r(this%n)), xmax = real(this%r(this%n)), ymin = -real(this%r(this%n)), ymax = real(this%r(this%n)))
    else
       call fig%init(xlabel =trim(xlabel), ylabel =trim(ylabel), width = coop_healpix_patch_default_figure_width, height = coop_healpix_patch_default_figure_height, xmin = -real(this%r(this%n)), xmax = real(this%r(this%n)), ymin = -real(this%r(this%n)), ymax = real(this%r(this%n)))          
    endif
    if(this%tbs%zmin(imap) .lt.0.99e30)then
       minz = this%tbs%zmin(imap)
    else
       call coop_array_get_threshold(this%image(:,:,imap), COOP_REAL_OF(0.99), minz)
    endif
    if(this%tbs%zmax(imap) .gt. -0.99e30)then
       maxz = this%tbs%zmax(imap)
    else
       call coop_array_get_threshold(this%image(:,:,imap), COOP_REAL_OF(0.01), maxz)
    endif
    if(coop_healpix_patch_default_want_label)then
       call coop_asy_density(fig, this%image(:,:,imap), -this%r(this%n), this%r(this%n), -this%r(this%n), this%r(this%n), label = trim(this%tbs%label(imap)), zmax = maxz, zmin = minz, color_table = trim(this%color_table))
    else
       call coop_asy_density(fig, this%image(:,:,imap), -this%r(this%n), this%r(this%n), -this%r(this%n), this%r(this%n), zmax = maxz, zmin = minz, color_table = trim(this%color_table))       
    endif
    if(use_rad .and. coop_healpix_patch_default_want_arrow)then
       theta = nint(2.d0*asin(this%r(this%n)/2.d0)/coop_SI_degree*10.d0)/10.d0
       call coop_asy_label(fig, "$\mathbf{-"//COOP_STR_OF(theta)//"}^\circ$", -this%r(this%n)*1.02, -this%r(this%n)*1.15, color="blue")
       call coop_asy_label(fig, "$\mathbf{"//COOP_STR_OF(theta)//"}^\circ$", this%r(this%n)*1.02, -this%r(this%n)*1.15, color="blue")
       call fig%arrow(this%r(this%n),  -this%r(this%n)*1.08, this%r(this%n),  -this%r(this%n)*1.01)
       call fig%arrow(-this%r(this%n),  -this%r(this%n)*1.08, -this%r(this%n),  -this%r(this%n)*1.01)
    endif
    if(this%tbs%spin(imap).eq.2 .and. this%tbs%headless_vector(imap) .and. this%nmaps .ge. 2)then
       if(imap .eq. this%nmaps)then
          iq = imap-1
          iu = imap
       elseif(this%tbs%spin(imap+1) .eq. 2)then
          iq = imap
          iu = imap+1
       elseif(imap .gt. 1)then
          if(this%tbs%spin(imap-1).eq.2)then
             iq = imap-1
             iu = imap
          else
             stop "patch_plot: spin 2 must go in pairs"
          endif
       else
          stop "patch_plot: spin 2 must go in pairs"
       endif
       bgrids = max(this%n/7, 2)
       norm = maxval(this%image(:,:,iq)**2+this%image(:,:,iu)**2)
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
       do j = -ns, ns, bgrids
          do i = -ns, ns, bgrids
             xc = i*this%dr
             yc = j*this%dr
             r = sqrt(this%image(i,j,iq)**2+this%image(i,j,iu)**2)*norm
             if(coop_healpix_IAU_headless_vector)then
                theta = 0.5d0*COOP_POLAR_ANGLE(this%image(i,j,iq), this%image(i,j,iu)) + coop_pio2                   
             else
                theta = 0.5d0*COOP_POLAR_ANGLE(this%image(i,j,iq), this%image(i,j,iu))
             endif
             if(this%tbs%local_rotation(imap))then
                theta  = theta + COOP_POLAR_ANGLE(xc, yc)
             endif
             k = k + 1
             xstart(k) = xc - r*cos(theta)
             ystart(k) = yc - r*sin(theta)
             xend(k) = 2*xc - xstart(k)
             yend(k) = 2*yc - ystart(k)
          enddo
       enddo
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
       call coop_asy_density(fig, remap, -ns*dk, ns*dk,-ns*dk, ns*dk, label = trim(adjustl(this%tbs%label(imap))), zmax = maxz, zmin = minz, color_table = trim(this%color_table))       
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
       call coop_asy_density(fig, immap, -ns*dk, ns*dk,-ns*dk, ns*dk, zmax = maxz, zmin = minz, label =trim(adjustl(this%tbs%label(imap))), color_table = trim(this%color_table))       
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
       call coop_asy_density(fig, immap, -ns*dk, ns*dk, -ns*dk, ns*dk, zmax = maxz, zmin = minz, label = trim(this%tbs%label(imap)), color_table = trim(this%color_table))
    endif
    call fig%close()    
    
    deallocate(fftmap,remap,immap)
  end subroutine coop_healpix_patch_plot_fft
 

  subroutine coop_healpix_patch_free(this)
    class(coop_healpix_patch) this
    if(allocated(this%image))deallocate(this%image)
    if(allocated(this%r))deallocate(this%r)
    if(allocated(this%fr))deallocate(this%fr)
    if(allocated(this%wcm))deallocate(this%wcm)
    if(allocated(this%wsm))deallocate(this%wsm)
    if(allocated(this%icm))deallocate(this%icm)
    if(allocated(this%nstack))deallocate(this%nstack)
    if(allocated(this%indisk))deallocate(this%indisk)
    call this%tbs%free()
    this%nmaps = 0
  end subroutine coop_healpix_patch_free

  subroutine coop_healpix_patch_init(this, genre, n, dr, mmax)
    class(coop_healpix_patch) this
    logical,parameter::do_norm = .true.
    logical,parameter::remove_m0 = .true.
    COOP_UNKNOWN_STRING::genre
    COOP_INT n
    COOP_REAL dr, cosmt, sinmt, theta, rij
    COOP_INT i,j,m
    COOP_INT, optional::mmax
    COOP_REAL sumrc(0:n+1), sumrs(0:n+1), weight(0:n+1)
    call this%free()
    call this%tbs%init(genre)
    this%nmaps = this%tbs%nmaps    
    this%n = n
    this%npix = (2*this%n+1)**2
    this%dr = dr
    if(present(mmax))then
       this%mmax = mmax
    else
       this%mmax = 4
    endif
    if(this%n .lt. 0) return
    allocate(this%image(-this%n:this%n, -this%n:this%n, this%nmaps))
    allocate(this%nstack(-this%n:this%n, -this%n:this%n))
    allocate(this%indisk(-this%n:this%n, -this%n:this%n))
    allocate(this%r(0:this%n))
    allocate(this%fr(0:this%n, 0:this%mmax/2, this%nmaps))
    allocate(this%wcm(-this%n:this%n, -this%n:this%n, 0:this%mmax+1))
    allocate(this%wsm(-this%n:this%n, -this%n:this%n, 0:this%mmax+1))
    allocate(this%icm(-this%n:this%n, -this%n:this%n, 0:1))
    this%image = 0.d0
    this%wcm = 0.d0
    this%wsm = 0.d0
    this%wcm(0,0,0) = 1.d0
    this%fr = 0.d0
    
    this%indisk = 1.d0
    do j=1, this%n
       i = ceiling(sqrt((this%n-j)*dble(this%n+j)+1.d-8))
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


    this%num_indisk_tol = floor(count(this%indisk .gt. 0.d0)*coop_healpix_mask_tol)+0.5d0
    !!------------------------------------------------------------------
    !!start doing m=0 
    !$omp parallel do private(i, j, rij)
    do j=-this%n, this%n
       do i=-this%n, this%n
          rij =  sqrt(dble(i**2+j**2))
          this%icm(i, j, 0) = floor(rij)
          rij = rij - this%icm(i,j,0)
          this%icm(i, j, 1) = this%icm(i, j, 0) + 1
          if(this%icm(i,j,0) .gt. n) &
               this%icm(i,j,0) = 0
          if(this%icm(i,j,1) .gt. n) &
               this%icm(i,j,1) = 0
          if(this%icm(i,j,0) .ne. 0) &
               this%wcm(i, j, 0) = 1.d0  - rij
          if(this%icm(i, j, 1) .ne. 0) &
               this%wcm(i, j, 1) = rij
       enddo
    enddo
    !$omp end parallel do

    sumrc = 1.d-30
    
    do j=-this%n, this%n
       do i=-this%n, this%n
          sumrc(this%icm(i,j,0)) = sumrc(this%icm(i,j,0)) + this%wcm(i,j,0)
          sumrc(this%icm(i,j,1)) = sumrc(this%icm(i,j,1)) + this%wcm(i,j,1)
       enddo
    enddo
    
    !$omp parallel do private(i, j)
    do j=-this%n, this%n
       do i=-this%n, this%n
          this%wcm(i, j, 0) = this%wcm(i, j, 0)/sumrc(this%icm(i, j, 0))
          this%wcm(i, j, 1) = this%wcm(i, j, 1)/sumrc(this%icm(i, j, 1))
       enddo
    enddo
    !$omp end parallel do

    !!m = 0 is done
    !!------------------------------------------------------------------
    
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
             endif
          enddo
       enddo

       if(do_norm)then
          !!now normalize
          sumrc = 1.d-6
          sumrs = 1.d-6
          do j=-this%n, this%n
             do i=-this%n, this%n
                if(this%icm(i, j, 0) .ne. 0)then
                   theta = atan2(dble(j), dble(i))
                   cosmt = cos(m*theta)
                   sumrc(this%icm(i,j,0)) = sumrc(this%icm(i,j,0)) + this%wcm(i,j,m)*cosmt
                   sumrc(this%icm(i,j,1)) = sumrc(this%icm(i,j,1)) + this%wcm(i,j,m+1)*cosmt
                   sinmt = sin(m*theta)                
                   sumrs(this%icm(i,j,0)) = sumrs(this%icm(i,j,0)) + this%wsm(i,j,m)*sinmt                

                   sumrs(this%icm(i,j,1)) = sumrs(this%icm(i,j,1)) + this%wsm(i,j,m+1)*sinmt
                endif
             enddo
          enddo
          !$omp parallel do private(i, j)
          do j=-this%n, this%n
             do i=-this%n, this%n
                if(this%icm(i, j, 0).ne. 0)then
                   if(sumrc(this%icm(i, j, 0)) .gt. 0.01d0) this%wcm(i, j, m) = this%wcm(i, j, m)/sumrc(this%icm(i, j, 0))
                   if(sumrs(this%icm(i, j, 0)) .gt. 0.01d0) this%wsm(i, j, m) = this%wsm(i, j, m)/sumrs(this%icm(i, j, 0))
                   if(sumrc(this%icm(i, j, 1)) .gt. 0.01d0) this%wcm(i, j, m+1) = this%wcm(i, j, m+1)/sumrc(this%icm(i, j, 1))
                   if(sumrs(this%icm(i, j, 1)) .gt. 0.01d0) this%wsm(i, j, m+1) = this%wsm(i, j, m+1)/sumrs(this%icm(i, j, 1))
                endif
             enddo
          enddo
          !$omp end parallel do
       end if
       if(remove_m0)then             !!remove m = 0 degeneracy
          sumrc = 0.d0
          sumrs = 0.d0
          do j=-this%n, this%n
             do i = -this%n, this%n
                if(this%icm(i,j,0).ne.0)then
                   theta = atan2(dble(j), dble(i))
                   sumrc(this%icm(i, j, 0)) = sumrc(this%icm(i, j, 0)) + this%wcm(i, j, m)
                   sumrc(this%icm(i, j, 1)) = sumrc(this%icm(i, j, 1)) + this%wcm(i, j, m+1)                
                   sumrs(this%icm(i, j, 0)) = sumrs(this%icm(i, j, 0)) + this%wsm(i, j, m)
                   sumrs(this%icm(i, j, 1)) = sumrs(this%icm(i, j, 1)) + this%wsm(i, j, m+1)
                endif
             enddo
          enddo
          !$omp parallel do private(i, j)
          do j=-this%n, this%n
             do i=-this%n, this%n
                if(this%icm(i, j, 0) .ne. 0)then
                   this%wcm(i, j, m) = this%wcm(i, j, m) - sumrc(this%icm(i, j, 0))*this%wcm(i, j, 0)
                   this%wcm(i, j, m+1) = this%wcm(i, j, m+1) - sumrc(this%icm(i, j, 1))*this%wcm(i, j, 1)
                   this%wsm(i, j, m) = this%wsm(i, j, m) - sumrs(this%icm(i, j, 0))*this%wcm(i, j, 0)
                   this%wsm(i, j, m+1) = this%wsm(i, j, m+1) - sumrs(this%icm(i, j, 1))*this%wcm(i, j, 1)                

                endif
             enddo
          enddo
          !$omp end parallel do
       endif
       
    enddo
    
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

  subroutine coop_healpix_maps_iqu2TQTUT(this, idone)
    class(coop_healpix_maps) this
    logical, optional::idone
    if(.not. present(idone))then
       call this%map2alm(index_list = (/ 1 /) )
    else
       if(.not. idone)then
          call this%map2alm(index_list = (/ 1 /) )
       endif
    endif
    this%alm(:,:,2) = this%alm(:,:,1)
    this%alm(:,:,3) = 0
    this%spin(2:3) = 2
    this%spin(1) = 0
    call this%alm2map( index_list = (/ 2, 3 /) )
  end subroutine coop_healpix_maps_iqu2TQTUT

  subroutine coop_healpix_maps_iqu2TQTUT1(this, idone)
    class(coop_healpix_maps) this
    logical,optional::idone
    COOP_INT l
    if(.not.present(idone))then
       call this%map2alm(index_list = (/ 1 /) )
    else
       if(.not. idone)then
          call this%map2alm(index_list = (/ 1 /) )
       endif
    endif
    this%alm(:,:, 3) = 0.
    do l = 2, this%lmax
       this%alm(l,:,2) = this%alm(l,:,1)*(l*(l+1.))
    enddo
    this%spin(2:3) = 2
    this%spin(1) = 0
    call this%alm2map( index_list = (/ 2, 3 /) )
  end subroutine coop_healpix_maps_iqu2TQTUT1


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

  subroutine coop_healpix_maps_read(this, filename, nmaps_wanted, spin, nmaps_to_read, known_size, nested)
    class(coop_healpix_maps) this
    COOP_UNKNOWN_STRING filename
    COOP_INT,optional::nmaps_wanted, nmaps_to_read
    COOP_INT,dimension(:),optional::spin
    integer(8) npixtot
    COOP_INT nmaps_actual
    logical,optional::known_size, nested
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
    if(this%nmaps.eq.0) stop "cannot read 0 maps"
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
       case(4)
          if(index(filename, "TQUL").ne.0)then
             this%spin(1) = 0
             this%spin(2:3) = 2
             this%spin(4) = 0
             this%iq = 2
             this%iu = 3
             if(coop_healpix_warning)write(*,*) "TQUL map default spin: 0 2 2 0"
          else
             this%spin = 0
             this%iq = 0
             this%iu = 0
             if(coop_healpix_warning)write(*,*) "I assume all maps are scalar, specify spins otherwise"
          endif
       case(3)
          if(index(filename, "TEB") .eq. 0 .and. index(filename, "teb") .eq. 0 )then
             this%spin(1) = 0
             this%spin(2:3) = 2
             this%iq = 2
             this%iu = 3
             if(coop_healpix_warning)write(*,*) "I assume it is an IQU map, specify spins otherwise"
          else
             this%spin = 0
             this%iq = 0
             this%iu = 0
             if(coop_healpix_warning)write(*,*) "I assume all maps are scalar, specify spins otherwise"
          endif
       case(2)  
          if(index(filename, "TE") .eq. 0 .and. index(filename, "EB").eq.0 .and. index(filename, "te") .eq. 0 .and. index(filename, "eb").eq.0 )then
             this%spin(1:2) = 2
             this%iq = 1
             this%iu = 2
             if(coop_healpix_warning)write(*,*) "I assume it is an QU map, specify spins otherwise"
          else
             this%iq = 0
             this%iu = 0
             this%spin = 0
             if(coop_healpix_warning)write(*,*) "I assume all maps are scalar, specify spins otherwise"
          endif
       case(1)
          this%iq = 0
          this%iu = 0
          this%spin = 0          
       case default
          this%iq = 0
          this%iu = 0
          this%spin = 0
          if(coop_healpix_warning)write(*,*) "I assume all maps are scalar, specify spins otherwise"
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
    if(present(nested))then
       if(nested)then
          call this%convert2nested()
       else
          call this%convert2ring()
       endif
    else
       call this%convert2ring()
    endif
    call write_minimal_header(this%header,dtype = 'MAP', nside=this%nside, order = this%ordering, creator='Zhiqi Huang', version = 'CosmoLib', units='muK', polar=any(this%spin.eq.2) )
    this%mask_npix = 0
    this%maskpol_npix = 0
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
             write(*,*) this%nmaps
             write(*,*) j
             write(*,*) index_list
             print*, this%spin
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



  subroutine coop_healpix_maps_get_disc(this, pix, disc)
    class(coop_healpix_maps)::this
    COOP_INT pix
    type(coop_healpix_disc) disc
    COOP_REAL r
#ifdef HAS_HEALPIX
    disc%nside  = this%nside
    disc%center = pix
    disc%ordering = this%ordering
    call this%pix2ang(pix, disc%theta, disc%phi)
    call ang2vec(disc%theta, disc%phi, disc%nz)
    disc%nx = (/  sin(disc%phi) , - cos(disc%phi) , 0.d0 /)
    call coop_vector_cross_product(disc%nz, disc%nx, disc%ny)
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_maps_get_disc

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
    if(disc%ordering .eq. COOP_RING)then
       call pix2vec_ring(disc%nside, pix, vec)
    else
       call pix2vec_nest(disc%nside, pix, vec)       
    endif
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
    if(disc%ordering .eq. COOP_RING)then
       call vec2pix_ring(disc%nside, vec, pix)
    else
       call vec2pix_nest(disc%nside, vec, pix)
    endif
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
    if(disc%ordering .eq. COOP_RING)then
       call pix2vec_ring(disc%nside, pix, vec)
    else
       call pix2vec_nest(disc%nside, pix, vec)
    endif
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
    if(disc%ordering .eq. COOP_RING)then
       call vec2pix_ring(disc%nside, vec, pix)
    else
       call vec2pix_nest(disc%nside, vec, pix)       
    endif
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_disc_xy2pix

  subroutine coop_healpix_rotate_qu(qu, phi, spin)
    COOP_SINGLE qu(2)
    COOP_INT spin
    COOP_REAL phi, cosp, sinp
    cosp = cos(spin*phi)
    sinp = sin(spin*phi)
    qu = (/ qu(1)*cosp + qu(2)*sinp,  -qu(1)*sinp + qu(2)*cosp /)
  end subroutine coop_healpix_rotate_qu


  subroutine coop_healpix_patch_get_fr0(patch, nvar, var)
    COOP_INT nvar
    type(coop_healpix_patch)::patch
    COOP_REAL var(nvar)
    call patch%get_radial_profile(1, 0)
    var = patch%fr(0:patch%n, 0, 1)
  end subroutine coop_healpix_patch_get_fr0


  subroutine coop_healpix_stack_on_patch(this, disc, angle, patch, tmp_patch, mask)
    class(coop_healpix_maps) this
    type(coop_healpix_disc) disc
    type(coop_healpix_maps),optional::mask
    COOP_REAL angle
    type(coop_healpix_patch) patch, tmp_patch
    if(present(mask))then
       call coop_healpix_fetch_patch(this, disc, angle, tmp_patch, mask)
       if(sum(tmp_patch%nstack*tmp_patch%indisk) .lt. patch%num_indisk_tol)return
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
    COOP_INT i, j, pix,  k
    COOP_REAL x, y, r, phi
    COOP_SINGLE qu(2)
    if(.not. present(mask))patch%nstack = 1.d0
    patch%nstack_raw  = 1
    if(all(patch%tbs%spin .eq. 0))then
       do j = -patch%n, patch%n
          do i = -patch%n, patch%n
             x = patch%dr * i
             y = patch%dr * j
             r = sqrt(x**2+y**2)
             phi = COOP_POLAR_ANGLE(x, y) + angle
             call disc%ang2pix( r, phi, pix)
             if(present(mask))then
                patch%nstack(i, j) = mask%map(pix, 1)
                patch%image(i, j, :) = this%map(pix, patch%tbs%ind)*mask%map(pix,1)
             else 
                patch%image(i, j, :) = this%map(pix,patch%tbs%ind)
             endif
          enddo
       enddo
    elseif(patch%nmaps .eq. 2 .and. all(patch%tbs%spin .ne. 0))then
       if(all(patch%tbs%local_rotation))then
          do j = -patch%n, patch%n
             do i = -patch%n, patch%n
                x = patch%dr * i
                y = patch%dr * j
                r = sqrt(x**2+y**2)
                phi = COOP_POLAR_ANGLE(x, y) + angle
                call disc%ang2pix( r, phi, pix)
                qu = this%map(pix, patch%tbs%ind)
                call coop_healpix_rotate_qu(qu, phi, patch%tbs%spin(1))
                qu = qu * coop_healpix_QrUrSign
                if(present(mask))then
                   patch%nstack(i, j) = mask%map(pix, 1)
                   patch%image(i, j, :) = qu*mask%map(pix,1)
                else 
                   patch%image(i, j, :) = qu
                endif
             enddo
          enddo
       else
          do j = -patch%n, patch%n
             do i = -patch%n, patch%n
                x = patch%dr * i
                y = patch%dr * j
                r = sqrt(x**2+y**2)
                phi = COOP_POLAR_ANGLE(x, y) + angle
                call disc%ang2pix( r, phi, pix)
                qu = this%map(pix, patch%tbs%ind)
                call coop_healpix_rotate_qu(qu, angle, patch%tbs%spin(1))
                if(present(mask))then
                   patch%nstack(i, j) = mask%map(pix, 1)
                   patch%image(i, j, :) = qu*mask%map(pix,1)
                else 
                   patch%image(i, j, :) = qu
                endif
             enddo
          enddo
       endif
    else
       do j = -patch%n, patch%n
          do i = -patch%n, patch%n
             x = patch%dr * i
             y = patch%dr * j
             r = sqrt(x**2+y**2)
             phi = COOP_POLAR_ANGLE(x, y) + angle
             call disc%ang2pix( r, phi, pix)
             k = 1
             do while(k.le. patch%nmaps)
                select case(patch%tbs%spin(k))
                case(0)
                   if(present(mask))then
                      patch%nstack(i, j) = mask%map(pix, 1)
                      patch%image(i, j, k) = this%map(pix, patch%tbs%ind(k))*mask%map(pix,1)
                   else 
                      patch%image(i, j, k) = this%map(pix,patch%tbs%ind(k))
                   endif
                case default
                   qu = this%map(pix, patch%tbs%ind(k:k+1))
                   if(patch%tbs%local_rotation(k))then
                      call coop_healpix_rotate_qu(qu, phi, patch%tbs%spin(k))
                      qu = qu * coop_healpix_QrUrSign
                   else
                      call coop_healpix_rotate_qu(qu, angle,patch%tbs%spin(k))
                   endif
                   if(present(mask))then
                      patch%nstack(i, j) = mask%map(pix, 1)
                      patch%image(i, j, k:k+1) = qu*mask%map(pix,1)
                   else 
                      patch%image(i, j, k:k+1) = qu
                   endif
                end select
                if(patch%tbs%spin(k).ne.0)then
                   k = k + 2
                else
                   k = k + 1
                endif
             enddo
          enddo
       enddo
    endif
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
       if(mask%nside .ne. map%nside) stop "coop_healpix_get_listpix: mask and map must have the same nside"
       if(mask%mask_npix .gt. 0)then
          total_weight = mask%mask_npix
       elseif(mask%maskpol_npix .gt. 0)then
          total_weight = mask%maskpol_npix
       else
          mask%mask_npix = count(mask%map(:,1).gt. 0.5)
          total_weight = mask%mask_npix          
       endif
       do_mask = .true.
    else
       total_weight = map%npix
    endif
    call map%convert2nested()
    if(do_mask)call mask%convert2nested()

    select case(trim(spot_type))
    case("Tmax", "Emax", "Bmax", "zetamax")  !!random orientation
       if(present(threshold))then
          if(threshold .lt. coop_stacking_max_threshold)then
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
          if(threshold.lt.coop_stacking_max_threshold)then
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
          if(threshold.lt.coop_stacking_max_threshold)then
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
          if(threshold.lt.coop_stacking_max_threshold)then
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
          if(threshold.lt. coop_stacking_max_threshold)then
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
          if(threshold.lt. coop_stacking_max_threshold)then          
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
    if(do_mask)call mask%convert2ring()
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_maps_get_listpix


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
    COOP_REAL smoothscale
    COOP_UNKNOWN_STRING, optional::output
    type(coop_healpix_maps) this
    call this%read(mask_file, nmaps_wanted = 1)
    call this%trim_mask( smoothscale)
    if(present(output))then
       call this%write(trim(output))
    else
       call this%write(trim(coop_file_add_postfix(trim(mask_file), "_trimmed")))
    endif
    call this%free
  end subroutine coop_healpix_trim_maskfile

  subroutine coop_healpix_trim_mask(this, smoothscale)
    real,parameter::nefolds = 4
    class(coop_healpix_maps) this
    type(coop_healpix_maps) hgs
    COOP_INT,dimension(:),allocatable::listpix
    COOP_INT i, j, nsteps
    COOP_REAL smoothscale
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

  subroutine coop_healpix_maps_query_disc(this, pix, radius, listpix, nlist)
    class(coop_healpix_maps)::this
    COOP_INT pix, listpix(0:), nlist
    COOP_REAL radius, vec(3)
#ifdef HAS_HEALPIX
    call this%pix2vec(pix, vec)
    if(this%ordering .eq. COOP_RING)then
       call query_disc(this%nside, vec, radius, listpix, nlist, nest = 0, inclusive = 0)
    else
       call query_disc(this%nside, vec, radius, listpix, nlist, nest = 1, inclusive = 0)
    endif
#else
      stop "CANNOT FIND HEALPIX"    
#endif    
    end subroutine coop_healpix_maps_query_disc

  subroutine coop_healpix_maps_query_strip(this, pix, theta1, theta2, listpix, nlist)
    class(coop_healpix_maps)::this
    COOP_INT pix, listpix(0:), nlist
    COOP_REAL theta1, theta2, vec(3)
#ifdef HAS_HEALPIX
    call this%pix2vec(pix, vec)
    if(this%ordering .eq. COOP_RING)then
       call query_strip(this%nside, theta1, theta2, listpix, nlist, nest = 0, inclusive = 0)
    else
       call query_strip(this%nside, theta1, theta2, listpix, nlist, nest = 1, inclusive = 0)
    endif
#else
      stop "CANNOT FIND HEALPIX"    
#endif    
    end subroutine coop_healpix_maps_query_strip
    

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
       call mask%pix2ang(pix, theta, phi)
       theta = coop_pi-theta
       phi = coop_pi +  phi
       call mask%ang2pix(theta, phi, conjpix)
       flip%map(pix, 1) = mask%map(conjpix, 1)
    enddo
    !$omp end parallel do
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_flip_mask

  subroutine coop_healpix_maps_mask(map, mask, polmask, remove_l_upper)
    class(coop_healpix_maps) map
    type(coop_healpix_maps) mask
    type(coop_healpix_maps),optional::polmask
    COOP_INT, optional::remove_l_upper
    COOP_SINGLE, dimension(:),allocatable::m01
    complex,dimension(:,:,:),allocatable::alms
    
    if(map%nside  .ne. mask%nside) stop "mask: nside must be the same as the map"
    call map%convert2ring()
    call mask%convert2ring()
    map%map(:,1) = map%map(:,1)*mask%map(:,1)
    if(present(remove_l_upper))then
       if(remove_l_upper .ge. 0)then
          allocate(alms(0:remove_l_upper, 0:remove_l_upper,1))
          if(mask%mask_npix .eq. 0)then
             mask%mask_npix = count(mask%map(:,1).gt. 0.5)
          endif
          allocate(m01(0:map%npix-1))       
          call map2alm(map%nside, remove_l_upper, remove_l_upper, map%map(:,1), alms)
          alms = (alms * mask%npix)/mask%mask_npix
          call alm2map(map%nside, remove_l_upper, remove_l_upper, alms, m01)
          map%map(:,1) = map%map(:,1) - m01*mask%map(:,1)
          deallocate(m01)
       endif
    endif
    if(map%nmaps.ge.3 .and. present(polmask))then
       call polmask%convert2ring()
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



  subroutine coop_healpix_maps_t2zeta(this, fwhm_arcmin, want_unconstrained)
    class(coop_healpix_maps)::this
    COOP_INT::l, i
    logical,optional::want_unconstrained 
    type(coop_cosmology_firstorder)::fod
    COOP_REAL,dimension(:,:),allocatable::Cls, Cls_lensed
    COOP_REAL::fwhm_arcmin, norm, ucnorm, CTT, Knorm
    logical douc
#if DO_ZETA_TRANS
    if(present(want_unconstrained))then
       douc = want_unconstrained
       if(douc .and. this%nmaps .lt. 3) &
            stop "maps_t2zeta for unconstrained map you need nmaps>=3"
    else
       douc = .false.
    endif
    if(maxval(this%map(:,1)) .gt. 1.)then !!unit is muK
       Knorm  = 1.d6
    else
       Knorm = 1.d0
    endif
    call this%map2alm(index_list = (/ 1 /) )
    
    call fod%Set_Planck_bestfit()
    call fod%compute_source(0)
    allocate(Cls(coop_num_cls, 2:this%lmax), Cls_lensed(coop_num_cls, 2:this%lmax))
    call fod%source(0)%get_All_Cls(2, this%lmax, Cls)
    call coop_get_lensing_Cls(2, this%lmax, Cls, Cls_lensed)
    
    norm = 1.d0/coop_healpix_zeta_normalization/ COOP_DEFAULT_TCMB / Knorm

    !$omp parallel do private(i, ucnorm, CTT)
    do l = 2, this%lmax
       CTT = ( Cls(coop_index_ClTT, l) + abs(Cls_lensed(coop_index_ClTT,l)) +  coop_Planck_TNoise(l)/coop_healpix_beam(l, fwhm_arcmin) )
       
       this%alm(l, :, 1) = this%alm(l,:,1) *  (Cls(coop_index_ClTzeta,l)/CTT * norm)!!the lensing Cl is added as a "noise"
       if(douc)then

          ucnorm = Cls(coop_index_Clzetazeta, l) -  Cls(coop_index_ClTzeta,l)**2/CTT
          if(ucnorm .lt. 0.d0) stop "clzetazeta - cltzeta**2/cltt negative?"
          ucnorm = sqrt(ucnorm)*norm
          this%alm(l, 0, 2) = coop_random_complex_Gaussian(.true.)*ucnorm
          do i=1, l
             this%alm(l, i, 2) = coop_random_complex_Gaussian()*ucnorm
          enddo
       endif
    enddo
    !$omp end parallel do
    if(douc)then
       this%alm(:,:,3) = this%alm(:,:,2) + this%alm(:,:,1)
       this%spin(1:3) = 0
       this%alm(0:1,:,3) = 0.
       call this%alm2map( index_list = (/ 1, 2, 3 /) )
       this%alm_done(1:3) = .false.
    else
       this%spin(1) = 0
       this%alm(0:1,:,1) = 0.
       call this%alm2map( index_list = (/ 1 /) )
       this%alm_done(1) = .false.
    endif
    deallocate(Cls, Cls_lensed)

#else
    stop "To use t2zeta you have to turn on DO_ZETA_TRANS in include/constants.h"
#endif
  end subroutine coop_healpix_maps_t2zeta


   subroutine coop_healpix_maps_E2zeta(this, fwhm_arcmin, want_unconstrained)
    class(coop_healpix_maps)::this
    COOP_INT::l, i
    logical,optional::want_unconstrained 
    type(coop_cosmology_firstorder)::fod
    COOP_REAL,dimension(:,:),allocatable::Cls, Cls_lensed
    COOP_REAL::fwhm_arcmin, norm, ucnorm, CEE, Knorm
    logical douc
#if DO_ZETA_TRANS
    if(present(want_unconstrained))then
       douc = want_unconstrained
       if(douc .and. this%nmaps .lt. 3) &
            stop "maps_t2zeta for unconstrained map you need nmaps>=3"
    else
       douc = .false.
    endif
    if(maxval(this%map(:,1)) .gt. 1.)then !!unit is muK
       Knorm  = 1.d6
    else
       Knorm = 1.d0
    endif
    if(this%spin(1).eq.0)then
       call this%map2alm(index_list = (/ 1 /) )
    else
       call this%map2alm(index_list = (/ 1, 2 /) )
    endif
    call fod%Set_Planck_bestfit()
    call fod%compute_source(0)
    allocate(Cls(coop_num_cls, 2:this%lmax), Cls_lensed(coop_num_cls, 2:this%lmax))
    call fod%source(0)%get_All_Cls(2, this%lmax, Cls)
    call coop_get_lensing_Cls(2, this%lmax, Cls, Cls_lensed)
    
    norm = 1.d0/coop_healpix_zeta_normalization/ COOP_DEFAULT_TCMB / Knorm
    !$omp parallel do private(i, ucnorm, CEE, l)
    do l = 2, this%lmax
       CEE = ( Cls(coop_index_ClEE, l) + abs(Cls_lensed(coop_index_ClEE,l)) +  coop_Planck_ENoise(l)/coop_healpix_beam(l, fwhm_arcmin) )
       
       this%alm(l, :, 1) = this%alm(l,:,1) *  (Cls(coop_index_ClEzeta,l)/CEE * norm)!!the lensing Cl is added as a "noise"
       if(douc)then

          ucnorm = Cls(coop_index_Clzetazeta, l) -  Cls(coop_index_ClEzeta,l)**2/CEE
          if(ucnorm .lt. 0.d0) stop "clzetazeta - cltzeta**2/cltt negative?"
          ucnorm = sqrt(ucnorm)*norm
          this%alm(l, 0, 2) = coop_random_complex_Gaussian(.true.)*ucnorm
          do i=1, l
             this%alm(l, i, 2) = coop_random_complex_Gaussian()*ucnorm
          enddo
       endif
    enddo
    !$omp end parallel do
    if(douc)then
       this%alm(:,:,3) = this%alm(:,:,2) + this%alm(:,:,1)
       this%spin(1:3) = 0
       this%alm(0:1,:,3) = 0.
       call this%alm2map( index_list = (/ 1, 2, 3 /) )
       this%alm_done(1:3) = .false.
    else
       this%spin(1) = 0
       this%alm(0:1,:,1) = 0.
       call this%alm2map( index_list = (/ 1 /) )
       this%alm_done(1) = .false.
    endif
    deallocate(Cls, Cls_lensed)

#else
    stop "To use t2zeta you have to turn on DO_ZETA_TRANS in include/constants.h"
#endif
  end subroutine coop_healpix_maps_E2zeta



  subroutine coop_healpix_maps_te2zeta(this, fwhm_arcmin, want_unconstrained)
    class(coop_healpix_maps)::this
    COOP_INT::l, i
    logical, optional::want_unconstrained
    type(coop_cosmology_firstorder)::fod
    COOP_REAL,dimension(:,:),allocatable::Cls, Cls_lensed
    COOP_REAL CTT, CTE, CEE, delta, bl
    COOP_REAL::fwhm_arcmin, norm, ucnorm, coef_T, coef_E, Knorm
    logical douc
#if DO_ZETA_TRANS    
    if(this%nmaps .lt. 3) stop "TE2Zeta not enough maps"
    if(present(want_unconstrained))then
       douc = want_unconstrained
    else
       douc = .false.
    endif
    if(this%spin(1).eq.0 .and. this%spin(2).eq.0)then  !!TEB maps
       call this%map2alm(index_list = (/ 1, 2 /) )
    else  !!IQU maps
       if(this%nmaps .lt. 3) stop "TE2Zeta requires IQU maps or TE maps"
       if(this%spin(1).ne. 0 .or. this%spin(2).ne.2 .or. this%spin(3) .ne. 2) stop "TE2Zeta spin does not match"
       call this%map2alm(index_list = (/ 1, 2, 3 /) )
    endif
    call fod%Set_Planck_bestfit()
    call fod%compute_source(0)
    allocate(Cls(coop_num_cls, 2:this%lmax), Cls_lensed(coop_num_cls, 2:this%lmax))
    call fod%source(0)%get_All_Cls(2, this%lmax, Cls)
    call coop_get_lensing_Cls(2, this%lmax, Cls, Cls_lensed)
    if(maxval(this%map(:,1)) .gt. 1.)then !!unit is muK
       Knorm  = 1.d6
    else
       Knorm = 1.d0
    endif
    
    norm = 1.d0/coop_healpix_zeta_normalization/COOP_DEFAULT_TCMB/Knorm

    !$omp parallel do private(CTT, CTE, CEE, delta, bl, ucnorm, coef_T, coef_E)
    do l = 2, this%lmax
       bl = coop_healpix_beam(l, fwhm_arcmin)
       CTT = Cls(coop_index_ClTT, l) + abs(Cls_lensed(coop_index_ClTT,l))+ coop_Planck_TNoise(l)/bl
       CEE = Cls(coop_index_ClEE, l) + abs(Cls_lensed(coop_index_ClEE, l)) + Coop_Planck_Enoise(l)/bl
       CTE = Cls(coop_index_ClTE, l)
       delta =  CTT*CEE-CTE**2
       coef_T = (Cls(coop_index_ClTzeta,l)*CEE - Cls(coop_index_ClEzeta,l) * CTE)/delta
       coef_E = (-Cls(coop_index_ClTzeta,l)*CTE + Cls(coop_index_ClEzeta,l) * CTT)/delta
       this%alm(l, :, 1) = (coef_T * norm) * this%alm(l,:,1) + (coef_E * norm) * this%alm(l,:,2)
       if(douc)then
          ucnorm = Cls(coop_index_Clzetazeta, l) - (CTT*coef_T**2 +CEE*coef_E**2 +  2.d0*coef_E*coef_T*CTE)
          if(ucnorm .lt. 0.d0) stop "clzetazeta - cltzeta**2/cltt negative?"
          ucnorm = sqrt(ucnorm)*norm
          this%alm(l, 0, 2) = coop_random_complex_Gaussian(.true.)*ucnorm
          do i=1, l
             this%alm(l, i, 2) = coop_random_complex_Gaussian()*ucnorm
          enddo
       endif
    enddo
    !$omp end parallel do
    if(douc)then
       this%alm(:,:,3) = this%alm(:,:,2) + this%alm(:,:,1)
       this%spin(1:3) = 0
       this%alm(0:1,:,1:3) = 0
       call this%alm2map( index_list = (/ 1, 2, 3 /) )
       this%alm_done(1:3) = .false.
    else
       this%spin(1) = 0
       this%alm(0:1,:,1) = 0
       call this%alm2map( index_list = (/ 1 /) )
       this%alm_done(1) = .false.
    endif
    deallocate(Cls, Cls_lensed)
#else
    stop "To use t2zeta you have to turn on DO_ZETA_TRANS in include/constants.h"
#endif
  end subroutine coop_healpix_maps_te2zeta


  function coop_healpix_beam(l, fwhm_arcmin) result(bl)
    COOP_INT l
    COOP_REAL fwhm_arcmin, bl
    bl = exp(-min(l*(l+1)*(fwhm_arcmin*coop_SI_arcmin*coop_sigma_by_fwhm)**2, 10.d0))
  end function coop_healpix_beam

  function coop_Planck_TNoise(l) result(Nl)
    COOP_INT l
    COOP_REAL Nl
    COOP_INT, parameter::fit_n = 10
    COOP_REAL, parameter,dimension(fit_n)::coef =  (/   -34.541512859596317      , &
  -2.3327044684438025      , &
  -5.2434038408357253E-2 , &
  0.10730432003605284      , &
  0.50599237614744652      , &
   4.4555688282625905E-2 , &
 -0.11595894402202822      , &
  -7.9077770071418474E-3 , &
 -0.26911077968221031      , &
  0.16444457464651557  /)
    COOP_REAL, parameter::rmin = coop_ln2, rmax = log(3000.d0)
    call coop_chebeval(fit_n, rmin, rmax, coef, log(dble(l)), Nl)
    Nl = exp(Nl)/COOP_DEFAULT_TCMB **2
  end function Coop_Planck_TNoise

  function coop_Planck_ENoise(l) result(Nl)
    COOP_INT l
    COOP_REAL Nl
    COOP_INT, parameter::fit_n = 10
    COOP_REAL, parameter,dimension(fit_n)::coef = (/  &
         -34.960562790571522      , &
         -0.87917841773973660      , &
         0.80533716121595145      , &
         0.19774868605364659      , &
         6.2456251840707022E-002 , &
         -6.4067689277299111E-002 , &
         6.9826114409022866E-003 , &
         -5.2937498702857466E-002 , &
         -1.3724094895074757E-002 , &
         -3.1217087209044592E-002 /)
    COOP_REAL, parameter::rmin = coop_ln2, rmax = log(3000.d0)
    call coop_chebeval(fit_n, rmin, rmax, coef, log(dble(l)), Nl)
    Nl = exp(Nl)/2.726**2
  end function Coop_Planck_ENoise
  
  function coop_Planck_BNoise(l) result(Nl)
    COOP_INT l
    COOP_REAL Nl
    COOP_INT, parameter::fit_n = 10
    COOP_REAL, parameter,dimension(fit_n)::coef =  (/ &
           -34.981682219713093      , &
           -0.88852205874105650      , &
           0.87526333776073173      , &
           0.13522507399381523      , &
           0.12803135081705186      , &
           -8.4600554185378374E-2 , &
           -1.3976431844054060E-2 , &
           -2.2053627296518163E-2 , &
           -3.6967260171160810E-2 , &
           -1.4665139293283336E-2 /)
    COOP_REAL, parameter::rmin = coop_ln2, rmax = log(3000.d0)
    call coop_chebeval(fit_n, rmin, rmax, coef, log(dble(l)), Nl)
    Nl = exp(Nl)/2.726**2
  end function Coop_Planck_BNoise

  function Coop_highpass_filter(l1, l2, l) result(w)
    COOP_INT l1, l2, l
    COOP_REAL w
    if(l.le. l1)then
       w = 0.d0
       return
    endif
    if(l.ge.l2)then
       w = 1.d0
       return
    endif
    w = sin(dble(l-l1)/dble(l2-l1)*coop_pio2)
  end function Coop_highpass_filter

  function coop_gaussian_filter(fwhm_arcmin, l) result(w)
    COOP_INT l
    COOP_REAL fwhm_arcmin, w
    w = exp(-((coop_sigma_by_fwhm*coop_SI_arcmin/coop_sqrt2)*fwhm_arcmin)**2*l*(l+1.d0))
  end function coop_gaussian_filter

  function coop_lowpass_filter(l1, l2, l) result(w)
    COOP_INT l1, l2, l
    COOP_REAL w
    if(l.ge. l2)then
       w = 0.d0
       return
    endif
    if(l.le.l1)then
       w = 1.d0
       return
    endif
    w = sin(dble(l2-l)/dble(l2-l1)*coop_pio2)
  end function coop_lowpass_filter


  subroutine coop_healpix_maps_get_peaks(this, sto, mask, restore)
    class(coop_healpix_maps)::this  
    type(coop_stacking_options)::sto
    type(coop_healpix_maps),optional::mask
    logical::domask
    COOP_INT i, nneigh, list(8), index_peak, ip
    COOP_INT::total_weight
    COOP_REAL::thetaphi(2)
    logical,optional::restore
#ifdef HAS_HEALPIX
    if(sto%nmaps .ne. this%nmaps)stop "get_peaks: nmaps mismatch"
    call sto%free()
    sto%nside = this%nside
    call this%convert2nested()
    if(present(mask))then
       if(mask%nside .ne. this%nside)then
          write(*,*) "map nside = ", this%nside
          write(*,*) "mask nside = ", mask%nside
          stop "coop_healpix_get_peaks: mask and map must have the same nside"
       endif
       call mask%convert2nested()
       total_weight = count(mask%map(:,1) .gt. 0.5)
       domask = .true.
    else
       domask = .false.
       total_weight = this%npix
    endif
    if(sto%index_I .ne. 0 .and. (abs(sto%I_lower_nu).lt.coop_stacking_max_threshold .or. abs(sto%I_upper_nu).lt. coop_stacking_max_threshold))then  !!rescale I
       if(domask)then
          sto%sigma_I  = sqrt(sum(dble(this%map(:,sto%index_I))**2, mask%map(:,1).gt.0.5)/total_weight)
       else
          sto%sigma_I = sqrt(sum(dble(this%map(:,sto%index_I))**2)/total_weight)
       endif
       sto%I_lower = sto%I_lower_nu* sto%sigma_I
       sto%I_upper = sto%I_upper_nu* sto%sigma_I
    endif
    if(sto%index_L .ne. 0 .and. (abs(sto%L_lower_nu).lt.coop_stacking_max_threshold .or. abs(sto%L_upper_nu).lt. coop_stacking_max_threshold) )then  !!rescale L
       if(domask)then
          sto%sigma_L  = sqrt(sum(dble(this%map(:,sto%index_L))**2, mask%map(:,1).gt.0.5)/total_weight)
       else
          sto%sigma_L = sqrt(sum(dble(this%map(:, sto%index_L))**2)/total_weight)
       endif
       sto%L_lower = sto%L_lower_nu* sto%sigma_L
       sto%L_upper = sto%L_upper_nu* sto%sigma_L
    endif
    if(sto%index_Q .ne. 0 .and. sto%index_U .ne. 0 .and. (abs(sto%P_lower_nu).lt.coop_stacking_max_threshold .or. abs(sto%P_upper_nu).lt. coop_stacking_max_threshold) )then  !!rescale P2
       if(domask)then
          sto%sigma_P  = sqrt(sum(dble(this%map(:,sto%index_Q)**2+this%map(:,sto%index_U)**2), mask%map(:,1).gt.0.5)/total_weight)
       else
          sto%sigma_P = sqrt(sum(dble(this%map(:,sto%index_Q)**2+this%map(:,sto%index_U)**2))/total_weight)
       endif
       sto%P2_lower = (sto%P_lower_nu * sto%sigma_P)**2
       sto%P2_upper = (sto%P_upper_nu * sto%sigma_P)**2
    endif

    select case(sto%genre)
    case(coop_stacking_genre_Imax, coop_stacking_genre_Imin, coop_stacking_genre_Imax_Oriented, coop_stacking_genre_Imin_Oriented)
       index_peak = sto%index_I
       if(index_peak .le. 0 .or. index_peak .gt. this%nmaps) stop "map index of peak overflow"
    case(coop_stacking_genre_Lmax, coop_stacking_genre_Lmin, coop_stacking_genre_Lmax_Oriented, coop_stacking_genre_Lmin_Oriented)
       index_peak = sto%index_L
       if(index_peak .le. 0 .or. index_peak .gt. this%nmaps) stop "map index of peak overflow"       
    case default
       index_peak = sto%index_Q
    end select
    if(sto%nested)then
       if(domask)then    
          select case(sto%genre)
          case(coop_stacking_genre_Imax, coop_stacking_genre_Imax_Oriented, coop_stacking_genre_Lmax, coop_stacking_genre_Lmax_Oriented)
             do i=0, this%npix-1
                if( mask%map(i, 1) .le. 0.5  .or. sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), index_peak).lt.this%map(i, index_peak)) .and. all(mask%map(list(1:nneigh),1) .gt. 0.5 ) )then
                   call sto%peak_pix%push(i)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_Imin, coop_stacking_genre_Imin_Oriented, coop_stacking_genre_Lmin, coop_stacking_genre_Lmin_Oriented)
             do i=0, this%npix-1
                if( mask%map(i, 1) .le. 0.5  .or. sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), index_peak) .gt. this%map(i, index_peak)) .and. all(mask%map(list(1:nneigh),1) .gt. 0.5 ) )then
                   call sto%peak_pix%push(i)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_Pmax_Oriented)
             do i=0, this%npix-1
                if( mask%map(i, 1) .le. 0.5  .or. sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), sto%index_Q)**2 + this%map(list(1:nneigh), sto%index_U)**2 .lt. this%map(i, sto%index_Q)**2 + this%map(i, sto%index_U)**2) .and. all(mask%map(list(1:nneigh),1) .gt. 0.5 ) )then
                   call sto%peak_pix%push(i)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_Pmin_Oriented)
             do i=0, this%npix-1
                if( mask%map(i, 1) .le. 0.5  .or. sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), sto%index_Q)**2 + this%map(list(1:nneigh), sto%index_U)**2 .gt. this%map(i, sto%index_Q)**2 + this%map(i, sto%index_U)**2) .and. all(mask%map(list(1:nneigh),1) .gt. 0.5 ) )then
                   call sto%peak_pix%push(i)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          end select
       else
          select case(sto%genre)
          case(coop_stacking_genre_Imax, coop_stacking_genre_Imax_Oriented, coop_stacking_genre_Lmax, coop_stacking_genre_Lmax_Oriented)
             do i=0, this%npix-1
                if( sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), index_peak).lt.this%map(i, index_peak)) )then
                   call sto%peak_pix%push(i)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_Imin, coop_stacking_genre_Imin_Oriented, coop_stacking_genre_Lmin, coop_stacking_genre_Lmin_Oriented)
             do i=0, this%npix-1
                if( sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), index_peak) .gt. this%map(i, index_peak))  )then
                   call sto%peak_pix%push(i)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_Pmax_Oriented)
             do i=0, this%npix-1
                if( sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), sto%index_Q)**2 + this%map(list(1:nneigh), sto%index_U)**2 .lt. this%map(i, sto%index_Q)**2 + this%map(i, sto%index_U)**2)  )then
                   call sto%peak_pix%push(i)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_Pmin_Oriented)
             do i=0, this%npix-1
                if(sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), sto%index_Q)**2 + this%map(list(1:nneigh), sto%index_U)**2 .gt. this%map(i, sto%index_Q)**2 + this%map(i, sto%index_U)**2) )then
                   call sto%peak_pix%push(i)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          end select          
       end if
    else
       if(domask)then    
          select case(sto%genre)
          case(coop_stacking_genre_Imax, coop_stacking_genre_Imax_Oriented, coop_stacking_genre_Lmax, coop_stacking_genre_Lmax_Oriented)
             do i=0, this%npix-1
                if( mask%map(i, 1) .le. 0.5  .or. sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), index_peak).lt.this%map(i, index_peak)) .and. all(mask%map(list(1:nneigh),1) .gt. 0.5 ) )then
                   call nest2ring(this%nside, i, ip)
                   call sto%peak_pix%push(ip)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_Imin, coop_stacking_genre_Imin_Oriented, coop_stacking_genre_Lmin, coop_stacking_genre_Lmin_Oriented)
             do i=0, this%npix-1
                if( mask%map(i, 1) .le. 0.5  .or. sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), index_peak) .gt. this%map(i, index_peak)) .and. all(mask%map(list(1:nneigh),1) .gt. 0.5 ) )then
                   call nest2ring(this%nside, i, ip)
                   call sto%peak_pix%push(ip)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_Pmax_Oriented)
             do i=0, this%npix-1
                if( mask%map(i, 1) .le. 0.5  .or. sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), sto%index_Q)**2 + this%map(list(1:nneigh), sto%index_U)**2 .lt. this%map(i, sto%index_Q)**2 + this%map(i, sto%index_U)**2) .and. all(mask%map(list(1:nneigh),1) .gt. 0.5 ) )then
                   call nest2ring(this%nside, i, ip)
                   call sto%peak_pix%push(ip)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_Pmin_Oriented)
             do i=0, this%npix-1
                if( mask%map(i, 1) .le. 0.5  .or. sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), sto%index_Q)**2 + this%map(list(1:nneigh), sto%index_U)**2 .gt. this%map(i, sto%index_Q)**2 + this%map(i, sto%index_U)**2) .and. all(mask%map(list(1:nneigh),1) .gt. 0.5 ) )then
                   call nest2ring(this%nside, i, ip)
                   call sto%peak_pix%push(ip)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          end select
       else
          select case(sto%genre)
          case(coop_stacking_genre_Imax, coop_stacking_genre_Imax_Oriented, coop_stacking_genre_Lmax, coop_stacking_genre_Lmax_Oriented)
             do i=0, this%npix-1
                if( sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), index_peak).lt.this%map(i, index_peak)) )then
                   call nest2ring(this%nside, i, ip)
                   call sto%peak_pix%push(ip)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_Imin, coop_stacking_genre_Imin_Oriented, coop_stacking_genre_Lmin, coop_stacking_genre_Lmin_Oriented)
             do i=0, this%npix-1
                if( sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), index_peak) .gt. this%map(i, index_peak))  )then
                   call nest2ring(this%nside, i, ip)
                   call sto%peak_pix%push(ip)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_Pmax_Oriented)
             do i=0, this%npix-1
                if( sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), sto%index_Q)**2 + this%map(list(1:nneigh), sto%index_U)**2 .lt. this%map(i, sto%index_Q)**2 + this%map(i, sto%index_U)**2) )then
                   call nest2ring(this%nside, i, ip)
                   call sto%peak_pix%push(ip)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_Pmin_Oriented)
             do i=0, this%npix-1
                if(  sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), sto%index_Q)**2 + this%map(list(1:nneigh), sto%index_U)**2 .gt. this%map(i, sto%index_Q)**2 + this%map(i, sto%index_U)**2) )then
                   call nest2ring(this%nside, i, ip)
                   call sto%peak_pix%push(ip)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          end select          
       end if
    endif
    if(present(restore))then
       if(restore)then
          call this%convert2ring()
          if(domask) call mask%convert2ring()
       endif
    else  !!by default I do restore
       call this%convert2ring()
       if(domask) call mask%convert2ring()       
    endif
#else
    stop "Cannot find HEALPIX"
#endif    
  end subroutine coop_healpix_maps_get_peaks



  



  
  subroutine coop_healpix_maps_mark_peaks(this, sto, imap)
    COOP_INT,parameter::max_n = 30
    class(coop_healpix_maps) this
    COOP_INT imap
    type(coop_stacking_options)::sto
    type(coop_healpix_disc) disc    
    COOP_INT i, pix, nx, ny, ix, iy, listpix(0:this%npix-1), nlist
    COOP_SINGLE::thetaphi(2)
    COOP_REAL  angle, e, r, rsq, dr, x, y, cost, sint, eplus, eminus
    dr = sqrt(coop_4pi/this%npix)/3.d0  !!just to be conservative
    this%map(:, imap) = 0.
    do i = 1, sto%peak_ang%n
       call sto%peak_ang%get_element(i, thetaphi)
       call this%ang2pix(dble(thetaphi(1)), dble(thetaphi(2)), pix)
       call this%get_disc(pix, disc)
       call sto%peak_get_angle_r_e(i, angle, r, e)
       cost = cos(angle)
       sint = sin(angle)
       eplus = 1.d0+e
       eminus = 1.d0-e
       nx = min(nint(r/sqrt(eminus)/dr), max_n)
       ny = min(nint(r/sqrt(eplus)/dr), max_n)
       rsq = r**2
       do ix = -nx, nx
          do iy = -ny, ny
             x = ix*dr
             y = iy*dr
             if(x**2*eminus+y**2*eplus .gt. rsq)cycle
             call disc%xy2pix(x*cost - y*sint, x*sint + y*cost, pix)
             this%map(pix, imap) = 1.
          enddo
       enddo
    enddo
  end subroutine coop_healpix_maps_mark_peaks

  subroutine coop_healpix_maps_stack_on_peaks(this, sto, patch, mask)
    COOP_INT,parameter::n_threads = 8
    class(coop_healpix_maps)::this
    type(coop_stacking_options)::sto
    type(coop_healpix_disc),dimension(n_threads)::disc
    type(coop_healpix_patch)::patch
    type(coop_healpix_patch),dimension(n_threads)::p, tmp
    type(coop_healpix_maps),optional::mask
    COOP_INT ithread, i
#ifdef HAS_HEALPIX
    patch%image = 0.d0
    patch%nstack = 0.d0
    patch%nstack_raw = 0
    if(sto%nested)then
       call this%convert2nested()
       if(present(mask))call mask%convert2nested()
    else
       call this%convert2ring()
       if(present(mask))call mask%convert2ring()
    endif
    !!adjust the index of maps
    i = 1
    do while(i.le.patch%nmaps)
       if(patch%tbs%ind(i) .gt. this%nmaps) stop "stack_on_peaks: index overflow"       
       do while(patch%tbs%spin(i) .ne. this%spin(patch%tbs%ind(i)))
          patch%tbs%ind(i)  = patch%tbs%ind(i) + 1
          if(patch%tbs%ind(i) .gt. this%nmaps) stop "stack_on_peaks: cannot find the map with corresponding spin"
       enddo
       if(patch%tbs%spin(i) .ne. 0)then
          if(i.ge.patch%nmaps .or. patch%tbs%ind(i) .ge. this%nmaps) stop "stack_on_peaks: nonzero spin must go in pairs"
          patch%tbs%ind(i+1) = patch%tbs%ind(i) + 1
          i = i+2
       else
          i = i+1
       endif
    enddo
    do ithread=1, n_threads
       p(ithread) = patch
       tmp(ithread) = patch
    enddo
    !$omp parallel do private(i, ithread)
    do ithread = 1, n_threads
       do i=ithread, sto%peak_pix%n, n_threads
          call this%get_disc(sto%pix(this%nside, i), disc(ithread))
          if(present(mask))then
             call coop_healpix_stack_on_patch(this, disc(ithread), sto%rotate_angle(i), p(ithread), tmp(ithread), mask)    
          else
             call coop_healpix_stack_on_patch(this, disc(ithread), sto%rotate_angle(i), p(ithread), tmp(ithread))
          endif
       enddo
    enddo
    !$omp end parallel do
    do ithread = 1, n_threads
       patch%image = patch%image + p(ithread)%image
       patch%nstack = patch%nstack + p(ithread)%nstack
       patch%nstack_raw = patch%nstack_raw + p(ithread)%nstack_raw
       call p(ithread)%free()
       call tmp(ithread)%free()
    enddo
    if(patch%nstack_raw .ne. 0)then
       do i=1, patch%nmaps
          patch%image(:, :, i) = patch%image(:, :, i)/max(patch%nstack, 1.d0)
       enddo
    else
       write(*,*) "warning: no patches has been found"
       patch%image = 0.d0
    endif
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_maps_stack_on_peaks
  
  
end module coop_healpix_mod









