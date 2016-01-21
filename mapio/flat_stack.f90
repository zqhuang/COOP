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
  COOP_REAL, parameter::reg_limit = 0.01
  COOP_UNKNOWN_STRING,parameter::mapdir = "act16/"
  COOP_UNKNOWN_STRING, parameter::postfix="7ar2"
  COOP_UNKNOWN_STRING,parameter::Ifile = mapdir//"dataCoadd_I_"//postfix//".fits"
  COOP_UNKNOWN_STRING,parameter::Qfile = mapdir//"dataCoadd_Q_"//postfix//".fits"
  COOP_UNKNOWN_STRING,parameter::Ufile = mapdir//"dataCoadd_U_"//postfix//".fits"
  COOP_UNKNOWN_STRING,parameter::Hitsfile = mapdir//"mask_"//postfix//".fits"
  COOP_UNKNOWN_STRING,parameter::PSfile = mapdir//"joinedClusterMasks_"//postfix//".fits"  
  type(coop_fits_image_cea)::imap, umap, qmap, hits, psmask
  type(coop_asy)::asy
  COOP_INT i, l
  type(coop_file) fp
  COOP_REAL, parameter::patchsize = 90.d0*coop_SI_arcmin
  COOP_UNKNOWN_STRING,parameter::output_dir = "ACTstacking/"
  COOP_REAL::mask_threshold
  COOP_REAL, parameter::smooth_scale = coop_SI_arcmin * fwhm_arcmin
  type(coop_healpix_maps)::hp, mask
  call coop_MPI_Init()
  call imap%open(Ifile)
  call qmap%open(Qfile)
  call umap%open(Ufile)
  if(coop_file_exists(hitsFile))then
     call hits%open(Hitsfile)
  else
     write(*,*) "Hits file "//trim(hitsfile)//" is not found; skipping..."     
     hits = imap
     hits%image = 1.
  endif
  if(coop_file_exists(PSfile))then
     call psmask%open(PSFile)
  else
     write(*,*) "Cluster mask file "//trim(psfile)//" is not found; skipping..."
     psmask = imap
     psmask%image = 1.
  endif
  print*,"# of pixels (I, Q, U, mask, ps):", imap%npix, qmap%npix, umap%npix, hits%npix, psmask%npix
  
  write(*,*) "max values:", maxval(abs(imap%image)), maxval(abs(qmap%image)), maxval(abs(umap%image))
  print*,"============comparing <T> in sources and global <T> ====="
  print*, "<T> in sources:", sum(imap%image*(1.-psmask%image)*hits%image)/sum((1.-psmask%image)*hits%image)
  print*, "global <T>:", sum(imap%image*hits%image)/sum(hits%image)  
  print*,"============masking sources===================="  
  hits%image = hits%image * psmask%image
  imap%image = imap%image*psmask%image
  qmap%image = qmap%image*psmask%image
  umap%image = umap%image*psmask%image
  write(*,*) "max values:", maxval(abs(imap%image)), maxval(abs(qmap%image)), maxval(abs(umap%image))  
  print*,"============regularizing===================="
  where(abs(imap%image) .gt. 300. .or. abs(qmap%image).gt.300 .or. abs(umap%image).gt. 300)
     hits%image = 0.
     imap%image = 0.
     qmap%image = 0.
     umap%image = 0.
  end where
  write(*,*) "max values:", maxval(abs(imap%image)), maxval(abs(qmap%image)), maxval(abs(umap%image))  
  print*,"=================================="  
  if(do_convert)then
     call hp%init(nside=2048, nmaps=3, genre="IQU", lmax=lmax)
     call mask%init(nside=2048, nmaps=1, genre="MASK", lmax=lmax)  
     call imap%convert2healpix(hp, 1, mask, hits=hits)
     call qmap%convert2healpix(hp, 2, mask, hits=hits)
     call umap%convert2healpix(hp, 3, mask, hits=hits)
     call hp%smooth(fwhm = smooth_scale, l_lower =lmin, l_upper = lmax)
     print*,"===== smoothed map min max ====="  
     print*, maxval(hp%map(:,1)), minval(hp%map(:,1))
     print*, maxval(hp%map(:,2)), minval(hp%map(:,2))
     print*, maxval(hp%map(:,3)), minval(hp%map(:,3))
     print*, "====  fsky = "//trim(coop_num2str(count(mask%map(:,1).gt.0.5)/dble(mask%npix)*100., "(F10.2)"))//"%======="
     
     print*,"=================================="  

     call hp%write(mapdir//"act_iqu_"//COOP_STR_OF(fwhm_arcmin)//"a_l"//COOP_STR_OF(lmin)//"-"//COOP_STR_OF(lmax)//".fits")
     call hp%write(mapdir//"act_qu_"//COOP_STR_OF(fwhm_arcmin)//"a_l"//COOP_STR_OF(lmin)//"-"//COOP_STR_OF(lmax)//".fits", index_list=(/ 2, 3/) )     
     call hp%get_QU()
     call hp%write(mapdir//"act_TQTUT_"//COOP_STR_OF(fwhm_arcmin)//"a_l"//COOP_STR_OF(lmin)//"-"//COOP_STR_OF(lmax)//".fits")
     call mask%write(mapdir//"act_mask.fits")
     stop
  endif
  mask_threshold = maxval(hits%image)*0.2
  where (hits%image .lt. mask_threshold)
     hits%image = 0.
  elsewhere
     hits%image = 1.
  end where
  call hits%get_flatmap(smooth_scale)
  imap%image = imap%image*hits%image
  qmap%image = qmap%image*hits%image
  umap%image = umap%image*hits%image    
  call imap%get_flatmap(smooth_scale)
  call qmap%get_flatmap(smooth_scale)
  call umap%get_flatmap(smooth_scale)
  call imap%get_QU(qmap, umap)
  call imap%smooth_flat(lmin = lmin, lmax = lmax)
  call qmap%smooth_flat(lmin = lmin, lmax = lmax)
  call umap%smooth_flat(lmin = lmin, lmax = lmax)
  call coop_random_init()
  call imap%find_extrema(hits, "spots/act_Tmax.txt", "Tmax", patchsize, irepeat)
  call imap%stack2fig("spots/act_Tmax.txt", "T", patchsize, output_dir//"act_T_onTmax.txt", caption="$T$ on $T_{\max}$", label = "$T (\mu K)$", color_table = "Rainbow")
!$  call system("../utils/fasy.sh "//output_dir//"act_T_onTmax.txt")  
  call imap%stack2fig("spots/act_Tmax.txt", "Qr", patchsize, output_dir//"act_Qr_onTmax.txt", caption="$Q_r$ on $T_{\max}$", label = "$Q_r (\mu K)$", color_table = "Rainbow")
!$  call system("../utils/fasy.sh "//output_dir//"act_Qr_onTmax.txt")

  call imap%stack2fig("spots/act_Tmax.txt", "Q", patchsize, output_dir//"act_Q_onTmax.txt", caption="$Q$ on $T_{\max}$", label = "$Q (\mu K)$", color_table = "Rainbow")
!$  call system("../utils/fasy.sh "//output_dir//"act_Q_onTmax.txt")
  
  call coop_MPI_FInalize()

end program test
