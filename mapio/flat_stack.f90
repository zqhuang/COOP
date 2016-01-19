program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  logical,parameter::do_convert = .true.
  COOP_INT,parameter::lmin = 350
  COOP_INT,parameter::lmax = 2500
  COOP_INT,parameter::irepeat = 1
  COOP_REAL, parameter::reg_limit = 0.01
  COOP_UNKNOWN_STRING,parameter::mapdir = "act16/"
  COOP_UNKNOWN_STRING, parameter::postfix="7ar2"
  COOP_UNKNOWN_STRING,parameter::Ifile = mapdir//"dataCoadd_I_"//postfix//".fits"
  COOP_UNKNOWN_STRING,parameter::Qfile = mapdir//"dataCoadd_Q_"//postfix//".fits"
  COOP_UNKNOWN_STRING,parameter::Ufile = mapdir//"dataCoadd_U_"//postfix//".fits"
  COOP_UNKNOWN_STRING,parameter::Imaskfile = mapdir//"mask_"//postfix//".fits"
  type(coop_fits_image_cea)::imap, umap, qmap, imask, cutmask
  type(coop_asy)::asy
  COOP_INT ix, iy, i, l
  type(coop_file) fp
  COOP_REAL, parameter::patchsize = 90.d0*coop_SI_arcmin
  COOP_UNKNOWN_STRING,parameter::output_dir = "ACTstacking/"
  COOP_REAL::mask_threshold
  COOP_REAL, parameter::smooth_scale = coop_SI_arcmin * 5.d0
  type(coop_healpix_maps)::hp, mask
  call coop_MPI_Init()
  call imap%open(Ifile)
  call qmap%open(Qfile)
  call umap%open(Ufile)
  if(coop_file_exists(imaskFile))then
     call imask%open(Imaskfile)
  else
     imask = imap
     imask%image = 1.
  endif
  print*,"=================================="    
  print*,"# of pixels (I, Q, U, mask):", imap%npix, qmap%npix, umap%npix, imask%npix
  print*,"======== regularized ============"  
  where(abs(imap%image) .gt. 500. .or. abs(qmap%image).gt.500 .or. abs(umap%image).gt. 500)
     imask%image = 0.
     imap%image = 0.
     qmap%image = 0.
     umap%image = 0.
  end where
  print*,"=================================="  
  if(do_convert)then
     call hp%init(nside=2048, nmaps=3, genre="IQU", lmax=lmax)
     call mask%init(nside=2048, nmaps=1, genre="MASK", lmax=lmax)  
     call imap%convert2healpix(hp, 1, mask, hits=imask)
     call qmap%convert2healpix(hp, 2, mask, hits=imask)
     call umap%convert2healpix(hp, 3, mask, hits=imask)
     call hp%smooth(fwhm = 5.*coop_SI_arcmin, l_lower =lmin, l_upper = lmax)
     print*,"===== smoothed map min max ====="  
     print*, maxval(hp%map(:,1)), minval(hp%map(:,1))
     print*, maxval(hp%map(:,2)), minval(hp%map(:,2))
     print*, maxval(hp%map(:,3)), minval(hp%map(:,3))
     print*, "====  fsky = "//trim(coop_num2str(count(mask%map(:,1).gt.0.5)/dble(mask%npix)*100., "(F10.2)"))//"%======="
     print*,"=================================="  

     call hp%write(mapdir//"act_iqu_5a_l"//COOP_STR_OF(lmin)//"-"//COOP_STR_OF(lmax)//".fits")
     call hp%write(mapdir//"act_qu_5a_l"//COOP_STR_OF(lmin)//"-"//COOP_STR_OF(lmax)//".fits", index_list=(/ 2, 3/) )     
     call hp%get_QU()
     call hp%write(mapdir//"act_TQTUT_5a_l"//COOP_STR_OF(lmin)//"-"//COOP_STR_OF(lmax)//".fits")
     call mask%write(mapdir//"act_mask.fits")
     stop
  endif
  mask_threshold = maxval(imask%image)*0.2
  where (imask%image .lt. mask_threshold)
     imask%image = 0.
  elsewhere
     imask%image = 1.
  end where
  call imask%get_flatmap(smooth_scale)
  imap%image = imap%image*imask%image
  qmap%image = qmap%image*imask%image
  umap%image = umap%image*imask%image    
  call imap%get_flatmap(smooth_scale)
  call qmap%get_flatmap(smooth_scale)
  call umap%get_flatmap(smooth_scale)
  call imap%get_QU(qmap, umap)
  call imap%smooth_flat(lmin = lmin, lmax = lmax)
  call qmap%smooth_flat(lmin = lmin, lmax = lmax)
  call umap%smooth_flat(lmin = lmin, lmax = lmax)
  call coop_random_init()
  call imap%find_extrema(imask, "spots/act_Tmax.txt", "Tmax", patchsize, irepeat)
  call imap%stack2fig("spots/act_Tmax.txt", "T", patchsize, output_dir//"act_T_onTmax.txt", caption="$T$ on $T_{\max}$", label = "$T (\mu K)$", color_table = "Rainbow")
  call system("../utils/fasy.sh "//output_dir//"act_T_onTmax.txt")  
  call imap%stack2fig("spots/act_Tmax.txt", "Qr", patchsize, output_dir//"act_Qr_onTmax.txt", caption="$Q_r$ on $T_{\max}$", label = "$Q_r (\mu K)$", color_table = "Rainbow")
  call system("../utils/fasy.sh "//output_dir//"act_Qr_onTmax.txt")

  call imap%stack2fig("spots/act_Tmax.txt", "Q", patchsize, output_dir//"act_Q_onTmax.txt", caption="$Q$ on $T_{\max}$", label = "$Q (\mu K)$", color_table = "Rainbow")
  call system("../utils/fasy.sh "//output_dir//"act_Q_onTmax.txt")
  

!!simulations

!!$  call imap%simulate_flat(lmin = lmin,lmax = lmax, cls_file = "cls.dat")
!!$  call imap%find_extrema(cutmask, "spots/simu_I"//mapid//"_Tmax.txt", "Tmax", patchsize, irepeat)
!!
!!$  call imap%stack2fig("spots/simu_I"//mapid//"_Tmax.txt", "Qr", patchsize, output_dir//"simu_Qr_onTmax.txt", caption="$Q_r$ on $T_{\max}$", label = "$Q_r (\mu K)$", color_table = "Rainbow")

  call coop_MPI_FInalize()

end program test
