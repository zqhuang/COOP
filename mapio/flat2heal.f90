program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  COOP_INT::lmin
  COOP_INT,parameter::lmax = 2500
  COOP_INT,parameter::fwhm_arcmin = 3
  COOP_REAL, parameter::reg_limit = 0.0005
  COOP_UNKNOWN_STRING,parameter::mapdir = "act16/"
  COOP_UNKNOWN_STRING,parameter::map_prefix = "deep56_coadd"
  COOP_UNKNOWN_STRING,parameter::weight_prefix = "deep56_weight"
  COOP_UNKNOWN_STRING,parameter::output_prefix = "act"

  COOP_UNKNOWN_STRING,parameter::Ifile = mapdir//map_prefix//"_I.fits"
  COOP_UNKNOWN_STRING,parameter::Qfile = mapdir//map_prefix//"_Q.fits"
  COOP_UNKNOWN_STRING,parameter::Ufile = mapdir//map_prefix//"_U.fits"
  COOP_UNKNOWN_STRING,parameter::I_Hitsfile = mapdir//weight_prefix//"_I.fits"
  COOP_UNKNOWN_STRING,parameter::Q_Hitsfile = mapdir//weight_prefix//"_Q.fits"
  COOP_UNKNOWN_STRING,parameter::U_Hitsfile = mapdir//weight_prefix//"_U.fits"
  COOP_UNKNOWN_STRING,parameter::PSfile = "NULL.fits"
  COOP_UNKNOWN_STRING,parameter::beam_file = mapdir//"beam_7ar2.txt"
  type(coop_fits_image_cea)::imap, umap, qmap, I_hits,Q_hits, U_hits, psmask
  type(coop_asy)::asy
  COOP_INT i, l
  type(coop_file) fp
  COOP_REAL::beam(0:lmax)
  COOP_REAL, parameter::smooth_scale = coop_SI_arcmin * fwhm_arcmin
  type(coop_healpix_maps)::hp, mask, polmask
  logical:: has_mask = .false.
  logical::has_hits = .false.
  call coop_get_Input(1, lmin)
  call hp%init(nside=2048, nmaps=3, genre="IQU", lmax=lmax)
  call mask%init(nside=2048, nmaps=1, genre="MASK")
  call polmask%init(nside = 2048, nmaps=1, genre = "MASK")
  call fp%open_skip_comments(beam_file)
  do l = 0, lmax
     read(fp%unit, *) i, beam(l)
     if(i.ne.l) stop "beam file error"
     beam(l) = coop_highpass_filter(lmin-10, lmin+10, l)/beam(l)
  enddo
  call fp%close()
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
     where(I_hits%image .lt. 3000.d0)
        imap%image = 0.d0
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
  call mask%write(mapdir//output_prefix//"_imask_"//COOP_STR_OF(fwhm_arcmin)//"a_l"//COOP_STR_OF(lmin)//"-"//COOP_STR_OF(lmax)//".fits")

  !!qmap
  call qmap%open(Qfile)
  if(has_mask)qmap%image = qmap%image*psmask%image
  if(coop_file_exists(Q_hitsFile))then
     call Q_hits%open(Q_Hitsfile)
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



  !! u map
  call umap%open(Ufile)
  if(coop_file_exists(U_hitsFile))then
     call U_hits%open(U_Hitsfile)
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
     call umap%convert2healpix(hp, 3, polmask, hits=U_hits)
     call U_hits%free()
  else
     if(has_mask)then
        call umap%convert2healpix(hp, 3, polmask, hits=psmask)
        call psmask%free()
     else
        call umap%convert2healpix(hp, 3, polmask)
     endif
  endif
  call umap%free()
  print*, "==== U fsky = "//trim(coop_num2str(count(mask%map(:,1).gt.0.5)/dble(mask%npix)*100., "(F10.2)"))//"%======="
  polmask%map(:,1) = polmask%map(:,1)*mask%map(:,1)
  call mask%free()
  call polmask%write(mapdir//output_prefix//"_polmask_"//COOP_STR_OF(fwhm_arcmin)//"a_l"//COOP_STR_OF(lmin)//"-"//COOP_STR_OF(lmax)//".fits")
  print*, "==== pol fsky = "//trim(coop_num2str(count(polmask%map(:,1).gt.0.5)/dble(polmask%npix)*100., "(F10.2)"))//"%======="
  call polmask%free()
  call hp%smooth_with_window(fwhm = smooth_scale, window = beam, lmax = lmax)
  print*,"===== smoothed map max min ====="  
  print*, maxval(hp%map(:,1)), minval(hp%map(:,1))
  print*, maxval(hp%map(:,2)), minval(hp%map(:,2))
  print*, maxval(hp%map(:,3)), minval(hp%map(:,3))
  print*,"=================================="  

  call hp%write(mapdir//output_prefix//"_qu_"//COOP_STR_OF(fwhm_arcmin)//"a_l"//COOP_STR_OF(lmin)//"-"//COOP_STR_OF(lmax)//".fits", index_list=(/ 2, 3/) )     
  call hp%get_QU()
  call hp%write(mapdir//output_prefix//"_TQTUT_"//COOP_STR_OF(fwhm_arcmin)//"a_l"//COOP_STR_OF(lmin)//"-"//COOP_STR_OF(lmax)//".fits")
  call hp%free()

   end program test
