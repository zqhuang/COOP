program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  COOP_INT::lmin
  COOP_INT,parameter::lmax = 2500
  COOP_INT::fwhm_arcmin
  COOP_INT,parameter::lx_cut =10
  COOP_INT,parameter::ly_cut =10
  COOP_UNKNOWN_STRING,parameter::mapdir = "act16/"
  COOP_UNKNOWN_STRING,parameter::weight_prefix = "deep56_weight"
  COOP_UNKNOWN_STRING,parameter::I_Weightsfile = mapdir//weight_prefix//"_I.fits"
  COOP_UNKNOWN_STRING,parameter::Q_Weightsfile = mapdir//weight_prefix//"_Q.fits"
  COOP_UNKNOWN_STRING,parameter::U_Weightsfile = mapdir//weight_prefix//"_U.fits"
  COOP_UNKNOWN_STRING,parameter::PSfile = "NULL.fits"
  COOP_UNKNOWN_STRING,parameter::beam_file = mapdir//"beam_7ar2.txt"
  type(coop_fits_image_cea)::imap, umap, qmap, I_weights,Q_weights, U_weights, psmask
  COOP_INT i, l
  type(coop_file) fp
  COOP_REAL::beam(0:lmax)
  COOP_REAL::fwhm
  type(coop_healpix_maps)::hp, mask, polmask
  logical:: has_mask = .false.
  logical::has_weights = .false.
  COOP_STRING::Ifile, Qfile, Ufile, map_prefix, output_prefix

  call coop_get_Input(1, map_prefix)
  call coop_get_Input(2, output_prefix)
  call coop_get_Input(3, lmin)
  call coop_get_Input(4, fwhm_arcmin)
  Ifile = mapdir//trim(map_prefix)//"_I.fits"
  Qfile = mapdir//trim(map_prefix)//"_Q.fits"
  Ufile = mapdir//trim(map_prefix)//"_U.fits"
  fwhm = coop_SI_arcmin * fwhm_arcmin

  call hp%init(nside=2048, nmaps=3, genre="IQU", lmax=lmax)
  call mask%init(nside=2048, nmaps=1, genre="MASK")
  call polmask%init(nside = 2048, nmaps=1, genre = "MASK")
  call fp%open_skip_comments(beam_file)
  do l = 0, lmax
     read(fp%unit, *) i, beam(l)
     if(i.ne.l) stop "beam file error"
  enddo
  call fp%close()
  if(coop_file_exists(PSfile))then
     call psmask%open(PSFile)
     has_mask = .true.
  else
     write(*,*) "Cluster mask file "//trim(psfile)//" is not found; skipping..."
     has_mask = .false.
  endif
  !!imap 
  call imap%open(Ifile)
  if(has_mask)imap%image = imap%image*psmask%image
  if(coop_file_exists(I_weightsFile))then
     call I_weights%open(I_Weightsfile)
     if(has_mask) I_weights%image = I_weights%image*psmask%image
     has_weights = .true.
  else
     write(*,*) "Weights file "//trim(I_weightsfile)//" is not found; skipping..." 
     has_weights = .false.
  endif
  call imap%smooth(fwhm = fwhm, highpass_l1 = lmin-20, highpass_l2 = lmin + 20, lmax = lmax, lx_cut = lx_cut, ly_cut = ly_cut, beam = beam)
  if(has_weights)then
     call imap%convert2healpix(hp, 1, mask, weights=I_weights)
     call I_weights%free()
  else
     if(has_mask)then
        call imap%convert2healpix(hp, 1, mask, weights=psmask)
     else
        call imap%convert2healpix(hp, 1, mask)
     endif
  endif
  call imap%free()
  print*, "I map done"
  print*, "==== I fsky = "//trim(coop_num2str(count(mask%map(:,1).gt.0.5)/dble(mask%npix)*100., "(F10.2)"))//"%======="
  call mask%write(mapdir//trim(output_prefix)//"_imask_"//COOP_STR_OF(fwhm_arcmin)//"a_l"//COOP_STR_OF(lmin)//"-"//COOP_STR_OF(lmax)//".fits")

  !!qmap
  call qmap%open(Qfile)
  if(has_mask)qmap%image = qmap%image*psmask%image
  if(coop_file_exists(Q_weightsFile))then
     call Q_weights%open(Q_Weightsfile)
     if(has_mask)Q_weights%image = Q_weights%image * psmask%image
     has_weights = .true.
  else
     write(*,*) "Weights file "//trim(Q_weightsfile)//" is not found; skipping..."    
     has_weights = .false.
  endif
  call qmap%smooth(fwhm = fwhm, highpass_l1 = lmin-20, highpass_l2 = lmin + 20, lmax = lmax, lx_cut = lx_cut, ly_cut = ly_cut, beam = beam)
  if(has_weights)then
     call qmap%convert2healpix(hp, 2, mask, weights=Q_weights)
     call Q_weights%free()
  else
     if(has_mask)then
        call qmap%convert2healpix(hp, 2, mask, weights=psmask)
     else
        call qmap%convert2healpix(hp, 2, mask)
     endif
  endif
  call qmap%free()
  print*, "==== Q fsky = "//trim(coop_num2str(count(mask%map(:,1).gt.0.5)/dble(mask%npix)*100., "(F10.2)"))//"%======="



  !! u map
  call umap%open(Ufile)
  if(coop_file_exists(U_weightsFile))then
     call U_weights%open(U_Weightsfile)
     if(has_mask)U_weights%image = U_weights%image * psmask%image
     has_weights = .true.
  else
     write(*,*) "Weights file "//trim(U_weightsfile)//" is not found; skipping..."
     has_weights = .false.
  endif
  call umap%smooth(fwhm = fwhm, highpass_l1 = lmin-20, highpass_l2 = lmin + 20, lmax = lmax, lx_cut = lx_cut, ly_cut = ly_cut, beam = beam)
  if(has_weights)then
     call umap%convert2healpix(hp, 3, polmask, weights=U_weights)
     call U_weights%free()
  else
     if(has_mask)then
        call umap%convert2healpix(hp, 3, polmask, weights=psmask)
        call psmask%free()
     else
        call umap%convert2healpix(hp, 3, polmask)
     endif
  endif
  call umap%free()
  print*, "==== U fsky = "//trim(coop_num2str(count(mask%map(:,1).gt.0.5)/dble(mask%npix)*100., "(F10.2)"))//"%======="
  polmask%map(:,1) = polmask%map(:,1)*mask%map(:,1)
  call mask%free()
  call polmask%write(mapdir//trim(output_prefix)//"_polmask_"//COOP_STR_OF(fwhm_arcmin)//"a_l"//COOP_STR_OF(lmin)//"-"//COOP_STR_OF(lmax)//".fits")
  print*, "==== pol fsky = "//trim(coop_num2str(count(polmask%map(:,1).gt.0.5)/dble(polmask%npix)*100., "(F10.2)"))//"%======="
  call polmask%free()
  print*,"===== smoothed map max min ====="  
  print*, maxval(hp%map(:,1)), minval(hp%map(:,1))
  print*, maxval(hp%map(:,2)), minval(hp%map(:,2))
  print*, maxval(hp%map(:,3)), minval(hp%map(:,3))
  print*,"=================================="  

  call hp%write(mapdir//trim(output_prefix)//"_qu_"//COOP_STR_OF(fwhm_arcmin)//"a_l"//COOP_STR_OF(lmin)//"-"//COOP_STR_OF(lmax)//".fits", index_list=(/ 2, 3/) )     
  call hp%get_QU()
  call hp%write(mapdir//trim(output_prefix)//"_TQTUT_"//COOP_STR_OF(fwhm_arcmin)//"a_l"//COOP_STR_OF(lmin)//"-"//COOP_STR_OF(lmax)//".fits")
  call hp%free()

   end program test
