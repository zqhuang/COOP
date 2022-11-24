program test
  use coop_healpix_mod
  use coop_wrapper_utils
  implicit none

#include "constants.h"

#ifdef HAS_HEALPIX
  COOP_STRING :: spot_type = "T"
  COOP_STRING :: map_file =  "planck14/dx11_v2_smica_int_cmb_010a_1024.fits"
  COOP_STRING :: imask_file ="planck14/dx11_v2_common_int_mask_010a_1024.fits"
  COOP_STRING:: polmask_file = "planck14/dx11_v2_common_pol_mask_010a_1024.fits"
  COOP_UNKNOWN_STRING,parameter:: color_table = "Rainbow"
  COOP_REAL,parameter::r_degree = 2.d0
  COOP_REAL,parameter::r=2.d0*sin(r_degree*coop_SI_degree/2.d0)
  COOP_INT, parameter::n = 36
  COOP_REAL,parameter::dr = r/n
  COOP_UNKNOWN_STRING, parameter :: prefix = "stacked/"
  COOP_STRING fout,fout2, caption, fname, inline
  COOP_INT,parameter::mmax = 4
  integer i, m
  type(coop_healpix_maps) map, mask
  type(coop_healpix_patch) patch
  logical::do_mask
  type(coop_asy)::fp
  COOP_REAL::zmin = 1.1e31
  COOP_REAL::zmax  = -1.1e31
  COOP_REAL::zmin2 = 1.1e31
  COOP_REAL::zmax2  = -1.1e31

  if(iargc() .ge. 5)then
     map_file = trim(adjustl(coop_InputArgs(1)))
     spots_file = trim(adjustl(coop_InputArgs(2)))
     spot_type = trim(adjustl(coop_InputArgs(3)))     
     imask_file = trim(adjustl(coop_InputArgs(4)))
     polmask_file = trim(adjustl(coop_InputArgs(5)))
     if(iargc().ge. 7)then
        inline = trim(adjustl(coop_InputArgs(6)))
        read(inline, *) zmin
        inline = trim(adjustl(coop_InputArgs(7)))
        read(inline, *) zmax        
     endif
     if(iargc().ge.9)then
        inline = trim(adjustl(coop_InputArgs(8)))
        read(inline, *) zmin2
        inline = trim(adjustl(coop_InputArgs(9)))
        read(inline, *) zmax2               
     endif
  endif
  coop_healpix_IAU_headless_vector = .true.       
  if(.not. coop_file_exists(spots_file))then
     write(*,*) "Cannot find the file "//trim(adjustl(spots_file))
     stop
  endif
  if(.not. coop_file_exists(map_file))then
     write(*,*) "Cannot find the file "//trim(adjustl(map_file))
     stop
  endif
  select case(trim(spot_type))
  case("T", "I", "zeta")
     if(trim(imask_file) .ne. "")then
        do_mask = .true.
        call mask%read(trim(imask_file), nmaps_wanted = 1 )
     else
        do_mask = .false.
     endif
  case("E","B")
     if(trim(polmask_file) .ne. "")then
        call mask%read(trim(polmask_file), nmaps_wanted = 1)
        do_mask = .true.
     else
        do_mask = .false.
     endif
  case( "QrUr", "QU")
     if(index(map_file, "QTUT").eq.0 .and. index(map_file, "qzuz").eq.0)then
        if(trim(polmask_file) .ne. "")then
           call mask%read(trim(polmask_file), nmaps_wanted = 1)
           do_mask = .true.
        else
           do_mask = .false.
        endif
     else
        if(trim(imask_file) .ne. "")then
           do_mask = .true.
           call mask%read(trim(imask_file), nmaps_wanted = 1 )
        else
           do_mask = .false.
        endif
     endif
  case default
     stop "Unknown spot_type"
  end select
  call map%read(trim(map_file))
  if(do_mask)then
     if(mask%nside .ne.map%nside) stop "mask must have the same nside"
     do i=1, min(map%nmaps, 3)
        map%map(:, i) = map%map(:, i)*mask%map(:, 1)
     enddo
  endif
  if(maxval(abs(map%map(0:100,1))).lt. 1.e-2)then
     map%map = map%map*1.e6
  endif

  if(index(spots_file, "_Tmax_QTUTOrient") .gt. 0)then
     caption = "$T$ maxima, oriented"
  elseif(index(spots_file, "_TQULb").gt.0)then
     caption = "$T$ maxima, $\nabla^2$ oriented"
     str_2nd = "\nu_e^{\rm upper}"
  elseif(index(spots_file, "_TQULt").gt.0)then
     caption = "$T$ maxima, $\nabla^2$ oriented"
     str_2nd = "\nu_e"     
  elseif(index(spots_file, "_Tmax").gt.0)then
     caption = "$T$ maxima, random orientation"
  elseif(index(spots_file, "_zetamax_qzuzOrient").gt.0)then
     caption = "$\zeta$ maxima, oriented"     
  elseif(index(spots_file, "_zetamax").gt.0)then
     caption = "$\zeta$ maxima, random orientation"
  elseif(index(spots_file, "_zetamax_qzuzorient").gt.0)then
     caption = "$\zeta$ maxima, random orientation"     
  elseif(index(spots_file, "_Emax").gt.0)then
     caption = "$E$ maxima, random orientation"
  elseif(index(spots_file, "_Bmax").gt.0)then
     caption = "$B$ maxima, random orientation"
  elseif(index(spots_file, "_PTmax").gt.0)then
     caption = "$P_T$ maxima, oriented"
  elseif(index(spots_file, "_PZmax").gt.0)then
     caption = "$P_\zeta$ maxima, oriented"     
  elseif(index(spots_file, "_Pmax_").gt.0)then
     caption = "$P$ maxima, oriented"
  elseif(index(spots_file, "_Tmin_QTUTOrient") .gt. 0)then
     caption = "$T$ minima, oriented"
  elseif(index(spots_file, "_Tmin").gt.0)then
     caption = "$T$ minima, random orientation"
  elseif(index(spots_file, "_zetamin_qzuzOrient").gt.0)then
     caption = "$\zeta$ minima, oriented"     
  elseif(index(spots_file, "_zetamin").gt.0)then
     caption = "$\zeta$ minima, random orientation"
  elseif(index(spots_file, "_zetamin_qzuzorient").gt.0)then
     caption = "$\zeta$ minima, random orientation"     
  elseif(index(spots_file, "_Emin").gt.0)then
     caption = "$E$ minima, random orientation"
  elseif(index(spots_file, "_Bmin").gt.0)then
     caption = "$B$ minima, random orientation"
  elseif(index(spots_file, "_PTmin").gt.0)then
     caption = "$P_T$ minima, oriented"
  elseif(index(spots_file, "_PZmin").gt.0)then
     caption = "$P_\zeta$ minima, oriented"     
  elseif(index(spots_file, "_Pmin").gt.0)then
     caption = "$P$ minima, oriented"
  else
     stop "Unknown spots_file class"
  endif

  if(index(spots_file, "_threshold0pt1").gt.0)then
     caption = trim(caption)//", $\nu=0.1$"
  elseif(index(spots_file, "_threshold0pt2").gt.0)then
     caption = trim(caption)//", $\nu=0.2$"
  elseif(index(spots_file, "_threshold0pt3").gt.0)then
     caption = trim(caption)//", $\nu=0.3$"
  elseif(index(spots_file, "_threshold0pt4").gt.0)then
     caption = trim(caption)//", $\nu=0.4$"
  elseif(index(spots_file, "_threshold0pt5").gt.0)then
     caption = trim(caption)//", $\nu=0.5$"
  elseif(index(spots_file, "_threshold0pt6").gt.0)then
     caption = trim(caption)//", $\nu=0.6$"
  elseif(index(spots_file, "_threshold0pt7").gt.0)then
     caption = trim(caption)//", $\nu=0.7$"
  elseif(index(spots_file, "_threshold0pt8").gt.0)then
     caption = trim(caption)//", $\nu=0.8$"
  elseif(index(spots_file, "_threshold0pt9").gt.0)then
     caption = trim(caption)//", $\nu=0.9$"
  elseif(index(spots_file, "_threshold1pt5").gt.0)then
     caption = trim(caption)//", $\nu=1.5$"
  elseif(index(spots_file, "_threshold2pt5").gt.0)then
     caption = trim(caption)//", $\nu=2.5$"     
  elseif(index(spots_file, "_threshold3pt5").gt.0)then
     caption = trim(caption)//", $\nu=3.5$"
  elseif(index(spots_file, "_threshold0").gt.0)then
     caption = trim(caption)//", $\nu=0$"     
  elseif(index(spots_file, "_threshold1").gt.0)then
     caption = trim(caption)//", $\nu=1$"
  elseif(index(spots_file, "_threshold2").gt.0)then
     caption = trim(caption)//", $\nu=2$"
  elseif(index(spots_file, "_threshold3").gt.0)then
     caption = trim(caption)//", $\nu=3$"
  endif
  
  
  if(index(spots_file, "_2ndthreshold0pt1").gt.0)then
     caption = trim(caption)//", $"//trim(str_2nd)//"=0.1$"
  elseif(index(spots_file, "_2ndthreshold0pt2").gt.0)then
     caption = trim(caption)//", $"//trim(str_2nd)//"=0.2$"
  elseif(index(spots_file, "_2ndthreshold0pt3").gt.0)then
     caption = trim(caption)//", $"//trim(str_2nd)//"=0.3$"
  elseif(index(spots_file, "_2ndthreshold0pt4").gt.0)then
     caption = trim(caption)//", $"//trim(str_2nd)//"=0.4$"
  elseif(index(spots_file, "_2ndthreshold0pt5").gt.0)then
     caption = trim(caption)//", $"//trim(str_2nd)//"=0.5$"
  elseif(index(spots_file, "_2ndthreshold0pt6").gt.0)then
     caption = trim(caption)//", $"//trim(str_2nd)//"=0.6$"
  elseif(index(spots_file, "_2ndthreshold0pt7").gt.0)then
     caption = trim(caption)//", $"//trim(str_2nd)//"=0.7$"
  elseif(index(spots_file, "_2ndthreshold0pt8").gt.0)then
     caption = trim(caption)//", $"//trim(str_2nd)//"=0.8$"
  elseif(index(spots_file, "_2ndthreshold0pt9").gt.0)then
     caption = trim(caption)//", $"//trim(str_2nd)//"=0.9$"
  elseif(index(spots_file, "_2ndthreshold1pt5").gt.0)then
     caption = trim(caption)//", $"//trim(str_2nd)//"=1.5$"
  elseif(index(spots_file, "_2ndthreshold2pt5").gt.0)then
     caption = trim(caption)//", $"//trim(str_2nd)//"=2.5$"     
  elseif(index(spots_file, "_2ndthreshold3pt5").gt.0)then
     caption = trim(caption)//", $"//trim(str_2nd)//"=3.5$"
  elseif(index(spots_file, "_2ndthreshold0").gt.0)then
     caption = trim(caption)//", $"//trim(str_2nd)//"=0$"
  elseif(index(spots_file, "_2ndthreshold1").gt.0)then
     caption = trim(caption)//", $"//trim(str_2nd)//"=1$"
  elseif(index(spots_file, "_2ndthreshold2").gt.0)then
     caption = trim(caption)//", $"//trim(str_2nd)//"=2$"
  elseif(index(spots_file, "_2ndthreshold3").gt.0)then
     caption = trim(caption)//", $"//trim(str_2nd)//"=3$"
  endif



  call patch%init(trim(spot_type), n, dr, mmax = mmax)
  if(do_mask)then
     call map%stack(patch, trim(spots_file), mask, do_weight = .true.)
  else
     call map%stack(patch, trim(spots_file), do_weight = .true.)
  endif
  patch%caption = trim(coop_num2str(patch%nstack_raw))//" patches on "//trim(caption)
  patch%color_table = color_table
  fname = coop_file_name_of(trim(spots_file))
  fout2 = ""
  select case(trim(spot_type))
  case("QrUr")
     if(index(map_file, "QTUT") .eq. 0)then
        fout = prefix//"Qr_on_"//trim(fname)
        fout2 = prefix//"Ur_on_"//trim(fname)
     else
        patch%label(1) = "$Q_{T,r}(\mu K)$"
        patch%label(2) = "$U_{T,r}(\mu K)$"
        fout = prefix//"QTr_on_"//trim(fname)
        fout2 = prefix//"UTr_on_"//trim(fname)
     endif
  case("QU")
     if(index(map_file, "QTUT") .eq. 0)then
        fout = prefix//"Q_on_"//trim(fname)
        fout2 = prefix//"U_on_"//trim(fname)
     else
        fout = prefix//"QT_on_"//trim(fname)
        fout2 = prefix//"UT_on_"//trim(fname)
        patch%label(1) = "$Q_T(\mu K)$"
        patch%label(2) = "$U_T(\mu K)$"
     endif
  case("T", "E", "B", "I", "zeta") 
     fout = prefix//trim(spot_type)//"_on_"//trim(fname)
  end select
  if(zmin .lt. zmax)then
     patch%zmin = zmin
     patch%zmax = zmax
  endif
  call patch%plot(imap = 1, output =trim(fout))
  if(trim(fout2).ne."")then
     if(zmin2 .lt. zmax2)then
        patch%zmin = zmin2
        patch%zmax = zmax2
     endif
     call patch%plot(imap = 2, output =trim(fout2))
  endif
  write(*,*) "the output file is: "//trim(fout)
  call system("../utils/fasy.sh "//trim(fout))
  if(trim(fout2).ne."")then
     call system("../utils/fasy.sh "//trim(fout2))
  endif
  call patch%get_all_radial_profiles()
  do m = 0, 4, 2
     call fp%open(trim(coop_file_add_postfix(fout, "_m"//COOP_STR_OF(m))))
     call fp%init(xlabel="$r$", ylabel="radial profile")
     call coop_asy_curve(fp, patch%r, patch%fr(:, m/2, 1))
     call fp%close()
  enddo
#else
  print*, "You need to install Healpix"
#endif  
end program test



