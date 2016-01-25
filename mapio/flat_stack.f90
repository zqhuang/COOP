program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  logical,parameter::do_convert = .true.
  COOP_INT,parameter::lmin = 250
  COOP_INT,parameter::lmax = 2500
  COOP_INT,parameter::irepeat = 1
  COOP_INT,parameter::fwhm_arcmin = 5
  COOP_REAL, parameter::reg_limit = 0.02
  COOP_UNKNOWN_STRING,parameter::mapdir = "act15/"
  COOP_UNKNOWN_STRING, parameter::postfix="4"
  COOP_UNKNOWN_STRING,parameter::Ifile = mapdir//"dataCoadd_I_"//postfix//".fits"
  COOP_UNKNOWN_STRING,parameter::Qfile = mapdir//"dataCoadd_Q_"//postfix//".fits"
  COOP_UNKNOWN_STRING,parameter::Ufile = mapdir//"dataCoadd_U_"//postfix//".fits"
  COOP_UNKNOWN_STRING,parameter::I_Hitsfile = mapdir//"weightMap_"//postfix//".fits"
  COOP_UNKNOWN_STRING,parameter::Q_Hitsfile = mapdir//"weightMap_"//postfix//"q.fits"
  COOP_UNKNOWN_STRING,parameter::U_Hitsfile = mapdir//"weightMap_"//postfix//"u.fits"
  !!  COOP_UNKNOWN_STRING,parameter::Hitsfile = mapdir//"NULLFILE.fits"
  !!  COOP_UNKNOWN_STRING,parameter::PSfile = mapdir//"NULLFILE.fits"
  COOP_UNKNOWN_STRING,parameter::PSfile = mapdir//"joinedClusterMasks_"//postfix//".fits"
  !!  COOP_UNKNOWN_STRING,parameter::PSfile = mapdir//"rahulMasksS2_"//postfix//"_cutlevel20.fits"    
  type(coop_fits_image_cea)::imap, umap, qmap, I_hits,Q_hits, U_hits, psmask
  type(coop_asy)::asy
  COOP_INT i, l
  type(coop_file) fp
  COOP_REAL, parameter::patchsize = 90.d0*coop_SI_arcmin
  COOP_UNKNOWN_STRING,parameter::output_dir = "ACTstacking/"
  COOP_REAL::mask_threshold
  COOP_REAL, parameter::smooth_scale = coop_SI_arcmin * fwhm_arcmin
  type(coop_healpix_maps)::hp, mask
  logical:: has_mask = .false.
  logical::has_hits = .false.
  call coop_MPI_Init()
  call hp%init(nside=2048, nmaps=3, genre="IQU", lmax=lmax)
  call mask%init(nside=2048, nmaps=1, genre="MASK", lmax=lmax)  

  if(coop_file_exists(PSfile))then
     call psmask%open(PSFile)
     has_mask = .true.
     imap%image = imap%image
  else
     write(*,*) "Cluster mask file "//trim(psfile)//" is not found; skipping..."
     has_mask = .false.
  endif
  !!imap 
  call imap%open(Ifile)
  if(has_mask)imap%image = imap%image*psmask%image
  if(coop_file_exists(I_hitsFile))then
     call I_hits%open(I_Hitsfile)
     where(I_hits%image .lt. 10. )
        imap%image = 0.
        I_hits%image = 0.
     end where
     if(has_mask) I_hits%image = I_hits%image*psmask%image
     has_hits = .true.
  else
     write(*,*) "Hits file "//trim(I_hitsfile)//" is not found; skipping..." 
     has_hits = .false.
  endif
  write(*,*) "Before regularization, I map min, max:", minval(imap%image), maxval(imap%image)
  call imap%regularize(reg_limit)
  write(*,*) "After regularization, I map min, max:", minval(imap%image), maxval(imap%image)
  if(has_hits)then
     call imap%convert2healpix(hp, 1, mask, hits=I_hits)
     call I_hits%free()
  else
     if(has_mask)then
        call imap%convert2healpix(hp, 1, mask, hits=psmask)
     else
        call imap%convert2healpix(hp, 1, mask)
     endif
  endif
  call imap%free()
  print*, "I map done"
  print*, "==== I fsky = "//trim(coop_num2str(count(mask%map(:,1).gt.0.5)/dble(mask%npix)*100., "(F10.2)"))//"%======="
  call mask%write(mapdir//"act_imask.fits")

  !!qmap
  call qmap%open(Qfile)
  if(has_mask)qmap%image = qmap%image*psmask%image
  if(coop_file_exists(Q_hitsFile))then
     call Q_hits%open(Q_Hitsfile)
     where(Q_hits%image .lt. 10.)
        qmap%image = 0.
        Q_hits%image = 0.
     end where
     if(has_mask)Q_hits%image = Q_hits%image * psmask%image
     has_hits = .true.
  else
     write(*,*) "Hits file "//trim(Q_hitsfile)//" is not found; skipping..."    
     has_hits = .false.
  endif
  write(*,*) "Before regularization, Q map min, max:", minval(qmap%image), maxval(qmap%image)
  call qmap%regularize(reg_limit)
  write(*,*) "After regularization, Q map min, max:", minval(qmap%image), maxval(qmap%image)
  if(has_hits)then
     call qmap%convert2healpix(hp, 2, mask, hits=Q_hits)
     call Q_hits%free()
  else
     if(has_mask)then
        call qmap%convert2healpix(hp, 2, mask, hits=psmask)
     else
        call qmap%convert2healpix(hp, 2, mask)
     endif
  endif
  call qmap%free()
  print*, "==== Q fsky = "//trim(coop_num2str(count(mask%map(:,1).gt.0.5)/dble(mask%npix)*100., "(F10.2)"))//"%======="

  call mask%write(mapdir//"act_qmask.fits")

  !! u map
  call umap%open(Ufile)
  if(coop_file_exists(U_hitsFile))then
     call U_hits%open(U_Hitsfile)
     where(U_hits%image .lt. 10.)
        umap%image = 0.
        U_hits%image = 0.
     end where
     if(has_mask)U_hits%image = U_hits%image * psmask%image
     has_hits = .true.
  else
     write(*,*) "Hits file "//trim(U_hitsfile)//" is not found; skipping..."
     has_hits = .false.
  endif
  write(*,*) "Before regularization, U map min, max:", minval(umap%image), maxval(umap%image)
  call umap%regularize(reg_limit)
  write(*,*) "After regularization, U map min, max:", minval(umap%image), maxval(umap%image)

  if(has_hits)then
     call umap%convert2healpix(hp, 3, mask, hits=U_hits)
     call U_hits%free()
  else
     if(has_mask)then
        call umap%convert2healpix(hp, 3, mask, hits=psmask)
        call psmask%free()
     else
        call umap%convert2healpix(hp, 3, mask)
     endif
  endif
  call umap%free()
  print*, "==== U fsky = "//trim(coop_num2str(count(mask%map(:,1).gt.0.5)/dble(mask%npix)*100., "(F10.2)"))//"%======="
  call mask%write(mapdir//"act_umask.fits")
  call mask%free()
  call hp%smooth(fwhm = smooth_scale, l_lower =lmin, l_upper = lmax)
  print*,"===== smoothed map max min ====="  
  print*, maxval(hp%map(:,1)), minval(hp%map(:,1))
  print*, maxval(hp%map(:,2)), minval(hp%map(:,2))
  print*, maxval(hp%map(:,3)), minval(hp%map(:,3))
  print*,"=================================="  

  call hp%write(mapdir//"act_iqu_"//COOP_STR_OF(fwhm_arcmin)//"a_l"//COOP_STR_OF(lmin)//"-"//COOP_STR_OF(lmax)//".fits")
  call hp%write(mapdir//"act_qu_"//COOP_STR_OF(fwhm_arcmin)//"a_l"//COOP_STR_OF(lmin)//"-"//COOP_STR_OF(lmax)//".fits", index_list=(/ 2, 3/) )     
  call hp%get_QU()
  call hp%write(mapdir//"act_TQTUT_"//COOP_STR_OF(fwhm_arcmin)//"a_l"//COOP_STR_OF(lmin)//"-"//COOP_STR_OF(lmax)//".fits")
  call hp%free()
!!$  else
!!$     call I_hits%get_flatmap(smooth_scale)
!!$     call Q_hits%get_flatmap(smooth_scale)
!!$     call U_hits%get_flatmap(smooth_scale)
!!$     imap%image = imap%image*I_hits%image
!!$     qmap%image = qmap%image*Q_hits%image
!!$     umap%image = umap%image*U_hits%image    
!!$     call imap%get_flatmap(smooth_scale)
!!$     call qmap%get_flatmap(smooth_scale)
!!$     call umap%get_flatmap(smooth_scale)
!!$     call imap%get_QU(qmap, umap)
!!$     call imap%smooth_flat(lmin = lmin, lmax = lmax)
!!$     call qmap%smooth_flat(lmin = lmin, lmax = lmax)
!!$     call umap%smooth_flat(lmin = lmin, lmax = lmax)
!!$     call coop_random_init()
!!$     call imap%find_extrema(I_hits, "spots/act_Tmax.txt", "Tmax", patchsize, irepeat)
!!$     call imap%stack2fig("spots/act_Tmax.txt", "T", patchsize, output_dir//"act_T_onTmax.txt", caption="$T$ on $T_{\max}$", label = "$T (\mu K)$", color_table = "Rainbow")
!!$     call imap%stack2fig("spots/act_Tmax.txt", "Qr", patchsize, output_dir//"act_Qr_onTmax.txt", caption="$Q_r$ on $T_{\max}$", label = "$Q_r (\mu K)$", color_table = "Rainbow")
!!$     call imap%stack2fig("spots/act_Tmax.txt", "Q", patchsize, output_dir//"act_Q_onTmax.txt", caption="$Q$ on $T_{\max}$", label = "$Q (\mu K)$", color_table = "Rainbow")
!!$  endif

   end program test
