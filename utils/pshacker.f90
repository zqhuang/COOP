MODULE ps_utils
  use file_io_utils
  USE general_utils
  use sort_utils
  IMPLICIT NONE
#include "utils.h"

  Integer(IB)::PS_Resolution_level = 0  !! 0, 1, 2, 3;  2 is usually sufficiently high resolution; 1 is not bad in many cases; sometimes you may want to use 0 for arxiv submission  
  real(dl)::ps_default_linewidth = 2.

  integer,PARAMETER::PS_N_CoorLabels=2**8
  integer,PARAMETER::PS_PathDepth=2**14
  UNKNOWN_STRING ,parameter::Default_FONT="Times-Roman"
  Integer(IB),parameter::Hermite_Interpolation_Option = INTERP_HERMITE
  Integer(IB),parameter::CubicSpline_Interpolation_Option = INTERP_SPLINE
  Integer(IB),parameter::Linear_Interpolation_Option = INTERP_LINEAR
  Integer(IB),parameter::Chebyshev_Interpolation_Option = INTERP_CHEBYSHEV


  TYPE PSBox
     integer XL,XR,YL,YR
  END TYPE PSBox

  Type PSFont
     STRING T
     integer S
  END Type PSFont

  Type PSColor
     Real(dl) R
     Real(dl) G
     Real(dl) B
  End Type PSColor

  TYPE PSPath
     real(dl) X(PS_PathDepth),Y(PS_PathDepth)
     integer N
  END TYPE PSPath

  Type PSPoint
     real(dl) x,y
  End Type PSPoint

  Type PSLine
     real(dl) xmin,ymin
     real(dl) xmax,ymax
  End Type PSLine

  Type PSGridPoint
     Integer(IB) x,y
     Integer(IB) Direction !! 1 down 2 right 3 up 4 left
     Real(dl) lambda,a,b
     logical used
  End Type PSGridPoint

  Type PSMark
     STRING str, pos
     SHORT_STRING color
     real(dl) size
  End type PSMark


  TYPE CurveSettings
     Integer(IB) FittingType !!0 no fitting, 1 fit, 2 fit and keep the original as a gray curve
     Integer(IB) FittingMethod !! 0 cubic spline, 1 chebyshev, 2+ chebyshev fit gaussianlike functions
     Integer(IB) FittingOrder
     Integer(IB) ProbRenorm 
     !! 0: usual normalization \int P(x) dx =1
     !! 1: normalization \int_{x_min}^{x_max} P(x) dx =1
     !! 2: normalization max P(x) = 1
     real(dl)SigmaRange
     integer NumBins
     Integer(IB) twodsmooth
     real(dl) LineWidth
     real(dl) dashlen(2)
     real(dl) dotlen(2)
     real(dl) dotdashlen(4)
     real(dl) xamp,yamp
  END TYPE Curvesettings

  TYPE PSCoor
     TYPE(PSBOX) Bound
     real(dl) OX,OY,XUR,YUR
     real(dl) xscale, yscale
     real(dl) CX(PS_N_COORLABELS),CY(PS_N_COORLABELS)
     integer LongTick,ShortTick
     integer CXTYPE(PS_N_COORLABELS),CYTYPE(PS_N_COORLABELS)
     !! 1 long tick; 2 long tick with number; 3 long tick with 10^(number)
     !! 4 short tick; 5 short tick with number;6 short tick with 10^(number)
     !! 7 nothing 
     integer XPOW,YPOW
     SHORT_STRING XFORM,YFORM
     STRING xlabel,ylabel
     logical UR !!if ur=.true., plot ticks on the upper and right axises
     logical YLabelRot !!if Ylabelrot=.true., rotate 90 degrees before printing ylabel
     logical XLog,YLog
     integer XTicks,YTicks
     logical AutoAdjust
     Logical FracBoundX,FracBoundY
     logical Locked
     Logical PlotXticks
     Logical PlotYticks
     integer(IB) YlabelX
  END TYPE PSCoor


  TYPE PSFrame
     STRING FileName
     integer FileUnit
     TYPE(PSBox)Bound
     TYPE(PSFont)Font
     TYPE(PSCoor)Coor
     type(pscolor) color
     TYPE(CurveSettings)CS
     STRING Dict
  END TYPE PSFrame


  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  INTERFACE PSPlotCurves
     MODULE PROCEDURE PSPlotXY, PSPlotXY_2, PSPlotXY_3, PSPlotXY_4
  END INTERFACE


  INTERFACE PSDefaultCoor
     MODULE PROCEDURE PSDefaultCoor_D,PSDefaultCoor_DARR, PSDefaultCoor_Darr1, PSDefaultCoor_Darr2
  END INTERFACE

  INTERFACE PSScale
     MODULE PROCEDURE PsScale_I,PsScale_D
  END INTERFACE

  INTERFACE PSSetLineWidth
     MODULE PROCEDURE PsSetLineWidth_I,PsSetLineWidth_D
  END INTERFACE

  INTERFACE PSSetRGBColor
     MODULE PROCEDURE PsSetRGBColor_D,PsSetRGBColor_STR,PsSetRGBColor_Color
  END INTERFACE

  INTERFACE PSCoorTrans
     MODULE PROCEDURE PSCoorTrans_D,PSCoorTrans_DD,PsCoorTrans_pt
  END INTERFACE

  INTERFACE PSMoveTo
     MODULE PROCEDURE PsMoveTo_D,PsMoveTo_I
  END INTERFACE

  INTERFACE PSLineTo
     MODULE PROCEDURE PSLineTo_D,PSLineTo_I
  END INTERFACE

  INTERFACE PSRMoveTo
     MODULE PROCEDURE PSRMoveTo_D,PSRMoveTo_I
  END INTERFACE

  INTERFACE PSRLineTo
     MODULE PROCEDURE PSRLineTo_D,PSRLineTo_I
  END INTERFACE

  INTERFACE PSCMoveTo
     MODULE PROCEDURE PsCMoveTo_D
  END INTERFACE


  INTERFACE PSPlotPoint
     MODULE PROCEDURE PSPlotPoint_D,PSPlotPoint_DV
  END INTERFACE

  INTERFACE PSPlotCurve
     MODULE PROCEDURE PSPlotCurve_D,PSPlotDenseCurve_D, PsPlotCurve_fromfile
  END INTERFACE

  INTERFACE PSPlotLine
     MODULE PROCEDURE PSPlotLine_D,PSPlotLine_Line,PSPlotLine_I
  END INTERFACE


  INTERFACE PsPlotBox
     MODULE PROCEDURE PsPlotBox,PsPlotBox_Color
  END INTERFACE

  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CONTAINS

  !!%%%%%%%%%%%%%%%%%%%%% PsStart %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !!%%%%%%%%%%%%%%%%%%%% PSDefaultFrame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine PSDefaultFrame(Frame,Width,Height,NoCoor)
    TYPE(PSFrame)Frame
    integer,optional::Width,Height
    Logical,Optional::Nocoor
    Frame%FileName="mypsplot.eps"
    Frame%FileUnit=New_File_Unit()
    Frame%Bound%XL=0
    Frame%Bound%YL=0
    IF(PRESENT(Width))then
       Frame%Bound%XR=Width
    ELSE
       Frame%Bound%XR=800
    endif
    IF(PRESENT(HEIGHT))then
       Frame%Bound%YR=Height
    ELSE
       Frame%Bound%YR=NINT(Frame%Bound%XR*0.75)
    endif
    Frame%Font%T="Times-Roman"
    Frame%Font%S=MAX(16,NINT(Frame%Bound%XR*0.036))
    frame%color%r=0.
    frame%color%g=0.
    frame%color%b=0.

    Frame%Coor%LongTick=10
    Frame%Coor%ShortTick=5
    
    If(Present(NoCoor))then
       Frame%Coor%Bound=Frame%Bound
       Frame%Coor%Locked=.true. 
       Frame%Coor%PlotXTicks = .false.
       Frame%COOR%PlotYticks = .false.
       Frame%Coor%AUTOADJUST = .false.
    else
       Frame%Coor%PlotXTicks = .true.
       Frame%COOR%PlotYticks = .true.
       Frame%Coor%LOCKED = .FALSE.
       Frame%Coor%Bound%XL=MAX(40,Frame%Font%S*5+4*Frame%Coor%LongTick)  !! Note that Coor%Bound is relative to Bound
       Frame%Coor%Bound%YL=MAX(25,Frame%Font%S*2+4*Frame%Coor%LongTick)
       Frame%Coor%Bound%XR=Frame%Bound%XR-Frame%Bound%XL-Frame%Font%S*2-Frame%Coor%LongTick*2
       Frame%Coor%Bound%YR=Frame%Bound%YR-Frame%Bound%YL-Frame%Font%S
       Frame%Coor%AUTOADJUST = .TRUE.
    endif

    Frame%Coor%OX=0._dl
    Frame%Coor%OY=0._dl
    Frame%Coor%XSCALE=1._dl
    Frame%Coor%YSCALE=1._dl
    Frame%Coor%CX=0._dl
    Frame%Coor%CY=0._dl
    Frame%Coor%CXTYPE=0
    Frame%Coor%CYTYPE=0
    Frame%Coor%XFORM=""
    Frame%Coor%YFORM=""
    Frame%Coor%FracBoundX=.false. !! use interger boundaries e.g. x from 0 to 5, y from -2 to 10 
    Frame%Coor%FracBoundY=.false.
    Frame%Coor%UR=.true. !! plot ticks on upper and right axis
    Frame%Coor%XPOW=0
    Frame%Coor%YPOW=0
    Frame%Coor%xlabel = "x"
    Frame%Coor%ylabel = "y"
    Frame%Coor%XLOG = .FALSE.
    Frame%Coor%YLOG = .FALSE.
    Frame%Coor%XTicks = 41
    Frame%Coor%YTicks = 41
    frame%coor%ylabelx = 0
    Frame%Dict = ""
    Frame%CS%ProbRenorm = 0 
    Frame%CS%FittingType = 0
    Frame%CS%FittingMethod = 0
    Frame%CS%SigmaRange = 5.
    Frame%CS%Numbins = 20
    Frame%CS%TwoDSmooth = 0
    Frame%CS%FittingOrder = 1
    Frame%CS%dashlen(1) = 10
    Frame%CS%dashlen(2) = 6
    Frame%CS%dotlen(1) = 3
    Frame%CS%dotlen(2) = 5
    Frame%CS%dotdashlen(1) = 10
    Frame%CS%dotdashlen(2) = 4
    Frame%CS%dotdashlen(3) = 3
    Frame%CS%dotdashlen(4) = 4
    Frame%CS%xamp=1.
    Frame%CS%yamp=1.
    Frame%CS%LineWidth=ps_default_linewidth
  END subroutine PSDefaultFrame

  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  !!%%%%%%%%%%%%%%%%%%%%%%%%%PSStart %%%%%%%%%%%%%%%%%%%%%%%%
  subroutine PSStart(Frame,BD)
    TYPE(PSFrame)Frame
    Type(PSBox),optional::BD
    If(Frame%CS%Numbins.le.0) Frame%CS%Numbins=32
    IF(Frame%FileUnit.le.0 .or. Frame%FileUnit.ge.100.or.trim(Frame%FileName).eq."")CALL PSDefaultFrame(Frame)
    IF(Frame%Coor%YLabel.ne."" .AND. Frame%Coor%AUTOADJUST)then
       Frame%Coor%Bound%XL=PSSlTex(Frame,Frame%Coor%YLabel)+Frame%Font%S*4+Frame%Coor%LongTick*2
       IF(Frame%Coor%Bound%XL.GT.Frame%Coor%Bound%XR/4)then
          Frame%Coor%Bound%XL=Frame%Font%S*4+Frame%Coor%LongTick*4
          Frame%Coor%YLabelRot=.TRUE.
       endif
    endif
    Frame%Coor%XUR=Frame%Coor%OX+Frame%Coor%XSCALE*(Frame%Coor%BOUND%XR-Frame%Coor%BOUND%XL)
    Frame%Coor%YUR=Frame%Coor%OY+Frame%Coor%YSCALE*(Frame%Coor%BOUND%YR-Frame%Coor%BOUND%YL)

    OPEN(Frame%FILEUNIT,FILE=trim(Frame%Filename),STATUS="Unknown",Form="formatted")
    WRITE(Frame%FileUnit,"(A)") "%!PS-Adobe-3.0 EPSF-3.0"
    IF(PRESENT(BD))then
       WRITE(Frame%FILEUNIT,"(8a)") "%%BoundingBox: "//trim(INT2STR(BD%XL/2))//" "  &
            //trim(INT2STR(BD%YL/2))//" "//trim(INT2STR((BD%XR+1)/2))//" "//trim(INT2STR((BD%YR+1)/2))
    ELSE
       WRITE(Frame%FILEUNIT,"(8a)") "%%BoundingBox: "//trim(INT2STR(Frame%BOUND%XL/2))//" "  &
            //trim(INT2STR(Frame%BOUND%YL/2))//" "//trim(INT2STR((Frame%BOUND%XR+1)/2))//" "//  &
            trim(INT2STR((Frame%BOUND%YR+1)/2))
    endif
    WRITE(Frame%FILEUNIT,"(A)")"%%Pages: 1"
    WRITE(Frame%FILEUNIT,"(A)")"%%Page: 1 1"
    call PsLoadDict(Frame) !! Load additional postscript dictionary. Sepcify Frame%Dict to be the path/filename.
    CALL PSLoadMacro(Frame) !! Intrinsic postscript macros that will be used in this package
    CALL PsScale(Frame,0.5d0,0.5d0) !! get higher resolution.
    CALL PSsetfont(Frame)
    CALL PSGSAVE(Frame)
    CALL PSOUTLINEBOX(Frame,Frame%BOUND)
    CALL PsPush(Frame,"clip")
    CALL PSTRANSLATE(Frame,Frame%BOUND%XL,Frame%BOUND%YL)
    Call PSSetLineWidth(Frame,Frame%CS%Linewidth)
    CALL PsNewPath(Frame)
  END subroutine PSStart


  subroutine PsDoubleFrames(Frame1, Frame2, filename, width, height, position)
    UNKNOWN_STRING   position, filename
    type(Psframe) Frame1, Frame2
    integer width,height
    call psdefaultframe(frame1)
    select case(position)
    case("LR")
       frame1%filename = filename
       frame1%bound%XL = 0
       Frame1%bound%yl = 0
       frame1%bound%xr = width
       frame1%bound%yr = height
       frame1%coor%bound%xl = nint(width*0.1)
       frame1%coor%bound%yl = nint(height*0.2)
       frame1%coor%bound%xr = nint(width*0.52)
       frame1%coor%bound%yr = nint(height*0.96)
       frame1%coor%autoadjust = .false.
       frame1%font%s = max(min(nint(width*0.05),nint(height*(0.05))),6)
       frame1%dict = "cal.dict"
       frame1%coor%ylabelrot = .true.
       frame2 = frame1
       frame2%coor%bound%xl= nint(width*0.54)
       frame2%coor%bound%xr = nint(width*0.96)
       frame2%coor%ylabel = ""
       frame2%coor%yform = "none"
       frame2%dict=""
       call psstart(frame1)
    case("UD")
       stop "psdoubleframe: this code yet to be written"
       frame2 = frame1
    case("SUB_UPPER_RIGHT")
       call psdefaultframe(frame1, width, height)
       frame1%dict = "cal.dict"
       frame1%filename = filename
       frame2 = frame1
       frame2%coor%bound%xl= nint(width*0.55)
       frame2%coor%bound%xr = frame1%coor%bound%xr - width/20
       frame2%coor%bound%yl= nint(height*0.6)
       frame2%coor%bound%yr = nint(height*0.92)
       frame2%coor%ylabel = ""
       frame2%coor%xlabel = ""
       call psstart(frame1)
    case("SUB_LOWER_RIGHT")
       call psdefaultframe(frame1, width, height)
       frame1%dict = "cal.dict"
       frame1%filename = filename
       frame2 = frame1
       frame2%coor%bound%xl= nint(width*0.58)
       frame2%coor%bound%xr = frame1%coor%bound%xr - width/20
       frame2%coor%bound%yl=  frame1%coor%bound%yl + frame1%font%s * 2
       frame2%coor%bound%yr = nint(height*0.52)
       frame2%coor%ylabel = ""
       frame2%coor%xlabel = ""
       call psstart(frame1)
    case default
       stop "unknown option in psdoubleframe"
    end select
  end subroutine PsDoubleFrames
  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !!%%%%%%%%%%%%%%%%%%%%%% PSEnd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Subroutine PSEnd(Frame)
    TYPE(PSFrame)Frame
    CALL PSGRESTORE(Frame)
    call pspush(frame,"showpage")
    call close_file_unit(Frame%fileunit)
  End Subroutine PSEnd
  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  !!%%%%%%%%%%%%%%%%%%%%%%%%%% shift the frame%%%%%%%%%%%%%%%%%%%%%%%%%
  Subroutine PSShift(Frame,X,Y)
    Type(PSFrame)Frame
    integer(IB) x,y
    Frame%Bound%XL=Frame%Bound%XL+X
    Frame%Bound%YL=Frame%Bound%YL+Y
    Frame%Bound%XR=Frame%Bound%XR+X
    Frame%Bound%YR=Frame%Bound%YR+Y
  End subroutine PSShift

  Subroutine PSRescale(Frame,xs,ys)
    Type(PSFrame)Frame
    real(dl) xs,ys
    Frame%Bound%XR=nint(Frame%Bound%XL+(Frame%Bound%XR-Frame%Bound%XL)*xs)
    Frame%Bound%yR=nint(Frame%Bound%yL+(Frame%Bound%yR-Frame%Bound%yL)*ys)
  end  Subroutine PSRescale
  !!+===================================================

  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine PS1DDistr(Frame, x, w)
    TYPE(PSFrame)Frame
    real(dl),dimension(:),INTENT(IN)::X
    real(dl),dimension(:),intent(IN),optional::w
    real(dl) SIGMA,DX,XMIN,XMAX,XBAR,YYMAX
    real(dl) XX(Frame%CS%Numbins),YY(Frame%CS%Numbins)
    integer I,N,K
    real(dl) tot
    REAL(dl) FRACINDX
    if(frame%CS%Numbins.le.0) stop "PS1DDistr: Numbins<=0"
    IF(PRESENT(W))then
       N=GETDIM("PS1DDistr",SIZE(X),SIZE(W))
    else
       N=SIZE(X)
    endif
    if(frame%coor%locked)then
       xmin=max(frame%coor%ox,minval(x))
       xmax=min(frame%coor%xur,maxval(x))
    elseif(frame%CS%sigmaRange.le.0.)then
       xmin=minval(x)
       xmax=maxval(x)
    else
       IF(PRESENT(W))then
          XBAR=SUM(X*W)/SUM(W)
          SIGMA=SQRT(SUM((X-XBAR)**2*W)/SUM(W))
       ELSE
          XBAR=SUM(X)/N
          SIGMA=SQRT(SUM((X-XBAR)**2)/N)
       endif
       XMIN=XBAR-Frame%CS%SigmaRange*SIGMA  
       XMAX=XBAR+Frame%CS%SigmaRange*SIGMA
       XMIN=MAX(MINVAL(X),XMIN) 
       XMAX=MIN(MAXVAL(X),XMAX)
    endif
    dx=(XMAX-XMIN)/(Frame%CS%Numbins-1)
    CALL FINDGEN(XX(1:Frame%CS%Numbins),XMIN,XMAX)
    YY=0.
    tot=0.
    DO I=1,N
       IF(X(I).LT.XMIN .OR. X(I).GT.XMAX)CYCLE
       FRACINDX=(X(I)-XMIN)/DX+1.
       K=INT(FRACINDX) !!1<=K<=M
       IF(K.ge.Frame%CS%Numbins)K=Frame%CS%Numbins-1 !1<=K<=M-1
       FRACINDX=FRACINDX-K
       IF(PRESENT(w))then
          YY(K)=YY(K)+W(I)*(1.-FRACINDX)
          YY(K+1)=YY(K+1)+W(I)*FRACINDX
          tot=tot+w(i)
       ELSE
          YY(K)=YY(K)+1.-FRACINDX
          YY(K+1)=YY(K+1)+FRACINDX
          tot=tot+1._dl
       endif
    enddo
    tot=tot+YY(1)+YY(Frame%CS%Numbins)
    YY(1)=YY(1)*2.
    YY(Frame%CS%Numbins)=YY(Frame%CS%Numbins)*2.    
    Select case(Frame%CS%ProbRenorm)
    case(1)
       if(present(w))then
          YY=YY/sum(w(1:N))/dx
       else
          YY=YY/N/dx
       endif
    case(2)
       YY=YY/maxval(YY)
    case default
       YY=YY/tot/dx
    end Select
    i=1
    K=frame%CS%Numbins
    If(Frame%CS%FittingMethod.eq.3)then
       YYMAX=MAXVAL(YY)*0.01
       DO WHILE(YY(I).LT.YYMAX.AND.I.LT.Frame%CS%NumBins)
          I=I+1
       enddo
       if(I.gt.1)I=I-1
       DO WHILE(YY(K).LT.YYMAX.and.k.gt.1)
          K=K-1
       enddo
       if(K.LT.Frame%CS%Numbins)K=K+1
    endif
    if(.not.Frame%Coor%Locked)then
       If(Frame%CS%ProbRenorm.ne.2)then
          Frame%Coor%PlotYTicks=.false. 
          CALL PSDefaultCoor(Frame,XX(I),XX(K),0._dl,MAXVAL(YY(I:K))*1.08D0)
       else
          CALL PSDefaultCoor(Frame,XX(I),XX(K),0._dl,1.1D0)
       endif
    endif
    IF(Frame%CS%FittingType.ge.1)then
       CALL PsSetGray(FRAME,0._dl)
       CALL PSPlotSmoothCurve(Frame,XX(I:K),YY(I:K),Frame%CS%FittingOrder)
       if(Frame%CS%FittingType.eq.2)then
          CALL PsSetGray(FRAME,0.65d0)
          CALL PSPLOTCURVE(FRAME,XX(I:K),YY(I:K))
       endif
    ELSE
       CALL PSPLOTCURVE(FRAME,XX(I:K),YY(I:K))
    endif
  END subroutine PS1DDistr

  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Subroutine PsPlotxy(X,Y,filename,xlabel,ylabel, xlog, ylog, smooth)
    real(dl),dimension(:),INTENT(IN)::X,Y
    UNKNOWN_STRING  ,optional::filename,xlabel,ylabel
    logical,optional::xlog,ylog
    integer,optional::smooth
    TYPE(PSFrame)Frame
    CALL PsDefaultFrame(Frame)
    if(present(xlog)) Frame%coor%xlog = xlog
    if(present(ylog)) Frame%coor%ylog = ylog
    IF(Present(FileName))then
       Frame%FileName = trim(FileName)
       call add_file_ext(Frame%FileName, "eps")
     endif
    IF(PRESENT(xlabel))Frame%Coor%xlabel=xlabel
    IF(PRESENT(ylabel))Frame%Coor%ylabel=ylabel
    CALL PsStart(Frame)
    CALL PsDefaultCoor(Frame,X,Y)
    if(present(smooth))then
       CALL PsPlotCurve(Frame,X,Y, smooth)
    else
       CALL PsPlotCurve(Frame,X,Y)
    endif
    CALL PsEnd(Frame)
  END Subroutine PsPlotxy

  Subroutine PsPlotxy_2(x1, y1, x2, y2, form1, form2, FileName, xlabel, ylabel, xlog, ylog)
    real(dl),dimension(:),INTENT(IN)::x1, x2, y1, y2
    real(dl) xmin, xmax, ymin, ymax, dx, dy
    UNKNOWN_STRING  ,optional::form1, form2, FileName, xlabel, ylabel
    TYPE(PSFrame)Frame
    logical,optional::xlog,ylog
    Call PsDefaultFrame(Frame, 800, 600)
    IF(PRESENT(filename))then
       Frame%filename=trim(filename)
       call add_file_ext(Frame%filename, "eps")
    end IF
    IF(PRESENT(xlabel))Frame%Coor%xlabel=trim(xlabel)
    IF(PRESENT(ylabel))Frame%Coor%ylabel=trim(ylabel)
    if(present(xlog)) Frame%coor%xlog = xlog
    if(present(ylog)) Frame%coor%ylog = ylog

    CALL PsStart(Frame)
    xmin = min(minval(x1),minval(x2))
    xmax = max(maxval(x1),maxval(x2))
    dx = (xmax-xmin)/25.
    ymin = min(minval(y1),minval(y2))
    ymax = max(maxval(y1), maxval(y2))
    dy = (ymax - ymin)/25.
    CALL PsDefaultCoor(Frame, xmin-dx, xmax+dx, ymin-dy , ymax+dy )
    if(present(form1))then
       call pssetcolor(Frame,trim(form1))
    else
       call pssetcolor(Frame,"solid_red")
    endif
    CALL PsPlotCurve(Frame,X1,Y1)
    if(present(form2))then
       call pssetcolor(Frame,trim(form2))
    else
       call pssetcolor(Frame,"dash_blue")
    endif
    CALL PsPlotCurve(Frame,X2,Y2)
    CALL PsEnd(Frame)
  END Subroutine PsPlotxy_2


  Subroutine PsPlotxy_3(x1, y1, x2, y2, x3, y3, form1, form2, form3, FileName, xlabel, ylabel, xlog, ylog)
    real(dl),dimension(:),INTENT(IN)::x1, x2, y1, y2, x3, y3
    UNKNOWN_STRING  ,optional::form1, form2, form3, FileName, xlabel, ylabel
    TYPE(PSFrame)Frame
    real(dl) xmin,ymin,xmax,ymax,dx,dy
    logical,optional::xlog,ylog
    CALL PsDefaultFrame(Frame, 800, 600)
    IF(PRESENT(filename))then
       Frame%filename=trim(filename)
       call add_file_ext(Frame%filename, "eps")       
    endif
    IF(PRESENT(xlabel))Frame%Coor%xlabel=trim(xlabel)
    IF(PRESENT(ylabel))Frame%Coor%ylabel=trim(ylabel)
    if(present(xlog)) Frame%coor%xlog = xlog
    if(present(ylog)) Frame%coor%ylog = ylog
    CALL PsStart(Frame)
    xmin = min(minval(x1),minval(x2),minval(x3))
    xmax = max(maxval(x1),maxval(x2),maxval(x3))
    dx = (xmax-xmin)/25.
    ymin = min(minval(y1),minval(y2),minval(y3))
    ymax = max(maxval(y1), maxval(y2),maxval(y3))
    dy = (ymax - ymin)/25.
    CALL PsDefaultCoor(Frame, xmin-dx, xmax+dx, ymin-dy , ymax+dy )
    if(present(form1))then
       call pssetcolor(Frame,trim(form1))
    else
       call pssetcolor(Frame,"solid_red")
    endif
    CALL PsPlotCurve(Frame,X1,Y1)
    if(present(form2))then
       call pssetcolor(Frame,trim(form2))
    else
       call pssetcolor(Frame,"dash_blue")
    endif
    CALL PsPlotCurve(Frame,X2,Y2)
    if(present(form3))then
       call pssetcolor(Frame,trim(form3))
    else
       call pssetcolor(Frame,"dot_green")
    endif
    CALL PsPlotCurve(Frame,X3,Y3)
    CALL PsEnd(Frame)
  END Subroutine PsPlotxy_3


  Subroutine PsPlotxy_4(x1, y1, x2, y2, x3, y3, x4, y4, form1, form2, form3, form4, FileName, xlabel, ylabel, xlog, ylog)
    real(dl),dimension(:),INTENT(IN)::x1, x2, y1, y2, x3, y3, x4, y4
    UNKNOWN_STRING  ,optional::form1, form2, form3, form4, FileName, xlabel, ylabel
    TYPE(PSFrame)Frame
    real(dl) xmin,ymin,xmax,ymax,dx,dy
    logical,optional::xlog,ylog
    CALL PsDefaultFrame(Frame, 800, 600)
    IF(PRESENT(filename))then
       Frame%filename=trim(filename)
       call add_file_ext(Frame%filename, "eps")       
    endif
    IF(PRESENT(xlabel))Frame%Coor%xlabel=trim(xlabel)
    IF(PRESENT(ylabel))Frame%Coor%ylabel=trim(ylabel)
    if(present(xlog)) Frame%coor%xlog = xlog
    if(present(ylog)) Frame%coor%ylog = ylog
    CALL PsStart(Frame)
    xmin = min(minval(x1),minval(x2),minval(x3),minval(x4))
    xmax = max(maxval(x1),maxval(x2),maxval(x3),maxval(x4))
    dx = (xmax-xmin)/25.
    ymin = min(minval(y1),minval(y2),minval(y3),minval(y4))
    ymax = max(maxval(y1), maxval(y2),maxval(y3),maxval(y4))
    dy = (ymax - ymin)/25.
    CALL PsDefaultCoor(Frame, xmin-dx, xmax+dx, ymin-dy , ymax+dy )
    if(present(form1))then
       call pssetcolor(Frame,trim(form1))
    else
       call pssetcolor(Frame,"solid_red")
    endif
    CALL PsPlotCurve(Frame,X1,Y1)
    if(present(form2))then
       call pssetcolor(Frame,trim(form2))
    else
       call pssetcolor(Frame,"dash_blue")
    endif
    CALL PsPlotCurve(Frame,X2,Y2)
    if(present(form3))then
       call pssetcolor(Frame,trim(form3))
    else
       call pssetcolor(Frame,"dot_green")
    endif
    CALL PsPlotCurve(Frame,X3,Y3)
    if(present(form4))then
       call pssetcolor(Frame,trim(form4))
    else
       call pssetcolor(Frame,"dotdash_violet")
    endif
    CALL PsPlotCurve(Frame,X4,Y4)
    CALL PsEnd(Frame)
  END Subroutine PsPlotxy_4

  SubRoutine PSPlotLogXY(X,Y,filename,xlabel,ylabel)
    real(dl),dimension(:),INTENT(IN)::X,Y
    UNKNOWN_STRING  ,optional::filename,xlabel,ylabel
    TYPE(PSFrame)Frame
    CALL PsDefaultFrame(Frame)
    Frame%Coor%xlog = .true.
    IF(PRESENT(filename))then
       Frame%filename=trim(filename)
       call add_file_ext(frame%filename,"eps")
    endif
    IF(PRESENT(xlabel))Frame%Coor%xlabel=xlabel
    IF(PRESENT(ylabel))Frame%Coor%ylabel=ylabel
    CALL PsStart(Frame)
    CALL PsDefaultCoor(Frame,log10(X),y)
    CALL PsPlotCurve(Frame,log10(X),Y)
    CALL PsEnd(Frame)
  End SubRoutine PSPlotLogXY

  SubRoutine PSPlotXLogY(X,Y,filename,xlabel,ylabel)
    real(dl),dimension(:),INTENT(IN)::X,Y
    UNKNOWN_STRING  ,optional::filename,xlabel,ylabel
    TYPE(PSFrame)Frame
    CALL PsDefaultFrame(Frame)
    Frame%Coor%Ylog = .true.
    IF(PRESENT(filename))then
       Frame%filename=trim(filename)
       call add_file_ext(frame%filename,"eps")
    endif
    IF(PRESENT(xlabel))Frame%Coor%xlabel=xlabel
    IF(PRESENT(ylabel))Frame%Coor%ylabel=ylabel
    CALL PsStart(Frame)
    CALL PsDefaultCoor(Frame,X,log10(Y))
    CALL PsPlotCurve(Frame,X,log10(Y))
    CALL PsEnd(Frame)
  End SubRoutine PSPlotXLogY

  SubRoutine PSPlotLogXLogY(X,Y,filename,xlabel,ylabel)
    real(dl),dimension(:),INTENT(IN)::X,Y
    UNKNOWN_STRING  ,optional::filename,xlabel,ylabel
    TYPE(PSFrame)Frame
    CALL PsDefaultFrame(Frame)
    Frame%Coor%Xlog = .true.
    Frame%Coor%Ylog = .true.
    IF(PRESENT(filename))then
       Frame%filename=trim(filename)
       call add_file_ext(frame%filename,"eps")
    endif
    IF(PRESENT(xlabel))Frame%Coor%xlabel=xlabel
    IF(PRESENT(ylabel))Frame%Coor%ylabel=ylabel
    CALL PsStart(Frame)
    CALL PsDefaultCoor(Frame,log10(X),log10(Y))
    CALL PsPlotCurve(Frame,log10(X),log10(Y))
    CALL PsEnd(Frame)
  End SubRoutine PSPlotLogXLogY

  SubRoutine PsPlotFunction(func,xstart,xend,FileName,xlabel,ylabel,rescale, xlog, ylog)
    integer(IB),parameter::N=512
    real(dl),external::func
    real(dl) xstart,xend
    UNKNOWN_STRING  ,optional::filename,xlabel,ylabel
    real(dl) x(N),y(N)
    integer(IB) i
    real(dl),optional::rescale
    logical,optional::xlog, ylog
    TYPE(PSFrame)Frame
    CALL PsDefaultFrame(Frame)
    IF(PRESENT(filename))then
       Frame%filename=trim(filename)
       call add_file_ext(frame%filename,"eps")
    endif
    IF(PRESENT(xlabel))Frame%Coor%xlabel=xlabel
    IF(PRESENT(ylabel))Frame%Coor%ylabel=ylabel
    if(present(xlog)) Frame%coor%xlog = xlog
    if(present(ylog)) Frame%coor%ylog = ylog
    CALL PsStart(Frame)


    if(Frame%coor%xlog)then
       call findgen(x,log10(xstart),log10(xend))
       do i=1,N
          y(i)=func(10._dl**x(i))
       enddo
    else
       call findgen(x,xstart,xend)
       do i=1,N
          y(i)=func(x(i))
       enddo
    endif
    if(Frame%coor%ylog) y = log10(abs(y))
    if(present(rescale))then
       y=y*rescale
    endif
    CALL PSDefaultCoor(Frame,X,Y)
    CALL PSPlotCurve(Frame,X,Y)
    CALL PsEnd(Frame)
  end SubRoutine PsPlotFunction

  SubRoutine PsPlotFunction2(func, xstart,xend, func2, xstart2, xend2, form1, form2, FileName,xlabel,ylabel,rescale, xlog, ylog)
    integer(IB),parameter::N=512
    real(dl),external::func, func2
    real(dl) xstart,xend, xstart2, xend2
    UNKNOWN_STRING  ,optional::filename,xlabel,ylabel, form1, form2
    real(dl) x(N),y(N), x2(N), y2(N)
    integer(IB) i
    real(dl),optional::rescale
    logical,optional::xlog, ylog
    TYPE(PSFrame)Frame
    CALL PsDefaultFrame(Frame)
    IF(PRESENT(filename))then
       Frame%filename=trim(filename)
       call add_file_ext(frame%filename,"eps")
    endif
    IF(PRESENT(xlabel))Frame%Coor%xlabel=xlabel
    IF(PRESENT(ylabel))Frame%Coor%ylabel=ylabel
    if(present(xlog)) Frame%coor%xlog = xlog
    if(present(ylog)) Frame%coor%ylog = ylog
    CALL PsStart(Frame)
    if(Frame%coor%xlog)then
       call findgen(x,log10((xstart)),log10((xend)))
       do i=1,N
          y(i)=func(10._dl**x(i))
       enddo
       call findgen(x2,log10((xstart2)),log10((xend2)))
       do i=1,N
          y2(i)=func2(10._dl**x2(i))
       enddo
    else
       call findgen(x,xstart,xend)
       do i=1,N
          y(i)=func(x(i))
       enddo
       call findgen(x2,xstart2,xend2)
       do i=1,N
          y2(i)=func2(x2(i))
       enddo
    endif
    if(Frame%coor%ylog)then
       if(any(y.le.0.d0) .or. any(y2.le.0.d0)) print*,"warning: plotting negative function in log scale, taking absolute values"
       y = log10(abs(y)+1.e-99_dl)
       y2 = log10(abs(y2)+1.e-99_dl)
    endif
    if(present(rescale))then
       y=y*rescale
       y2 = y2*rescale
    endif
    CALL PSDefaultCoor(Frame,min(minval(x),minval(x2)),max(maxval(x),maxval(x2)),min(minval(y),minval(y2)), max(maxval(y), maxval(y2)))
    if(present(form1))call pssetcolor(Frame, form1)
    CALL PSPlotCurve(Frame,X,Y)
    if(present(form2))then
       call pssetcolor(Frame, form2)
    else
       call psdash(Frame)
    endif
    CALL PSPlotCurve(Frame,X2,Y2)
    CALL PsEnd(Frame)
  end SubRoutine PsPlotFunction2





  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine PsLabelX(Frame,label)
    UNKNOWN_STRING  ,optional::label
    TYPE(PSFrame)Frame
    STRING xlabel
    IF(PRESENT(LABEL))then
       xlabel=Trim(label)
    ELSE
       xlabel=Trim(Frame%Coor%xlabel)
    endif
    CALL PsMoveTo(FRAME,(Frame%Coor%BOUND%XL+Frame%Coor%BOUND%XR)/2, &
         MAX(Frame%Coor%Bound%YL-Frame%Font%S*2-Frame%Coor%LongTick*2,Frame%Coor%LongTick))
    CALL PSPrtTexAtCenter(Frame,xlabel) 
  END subroutine PsLabelX

  subroutine PsLabelY(Frame,LABEL)
    UNKNOWN_STRING  ,optional::LABEL
    TYPE(PSFrame)Frame
    STRING  ylabel
    IF(PRESENT(LABEL))then
       ylabel=trim(LABEL)
    ELSE
       ylabel=trim(Frame%Coor%ylabel)
    endif
    IF(Frame%Coor%YlabelRot)then
       IF(Frame%Coor%ylabelx.gt.0)then
          CALL PsMoveTo(FRAME,Frame%coor%ylabelx,(Frame%Coor%BOUND%YL+Frame%Coor%BOUND%YR)/2)          
       ELSEIF(FRAME%COOR%CYTYPE(1).GT.7)then
          CALL PsMoveTo(FRAME,FRAME%COOR%BOUND%XL-FRAME%COOR%LONGTICK,(Frame%Coor%BOUND%YL+Frame%Coor%BOUND%YR)/2)
       ELSE
          CALL PSMoveto(Frame,MAX(Frame%Font%S*4/3, &
               Frame%Coor%BOUND%XL-5*Frame%Coor%LongTick- &
               Frame%Font%S*(LEN_trim(NUM2STR(Frame%Coor%CY(1),"f")))), &
               (Frame%Coor%BOUND%YL+Frame%Coor%BOUND%YR)/2 ) 
       endif
       CALL PSGSAVE(FRAME)
       CALL PSROTATE(FRAME,90)
       CALL PSPrtTexAtCenter(Frame,Frame%Coor%ylabel) 
       CALL PSGRESTORE(FRAME)
    ELSE
       IF(Frame%Coor%ylabelx.gt.0)then
          CALL PsMoveTo(FRAME,Frame%coor%ylabelx,(Frame%Coor%BOUND%YL+Frame%Coor%BOUND%YR)/2)  
       elseif(FRAME%COOR%CYTYPE(1).GT.7)then
          CALL PsMoveTo(FRAME,FRAME%COOR%BOUND%XL-FRAME%COOR%LONGTICK,(Frame%Coor%BOUND%YL+Frame%Coor%BOUND%YR)/2+Frame%Font%S)
       ELSE
          CALL PSMoveto(Frame,MAX(Frame%Font%S, &
               Frame%Coor%BOUND%XL-2*Frame%Coor%LongTick- &
               Frame%Font%S*LEN_trim(NUM2STR(Frame%Coor%CY(1),Frame%Coor%YFORM))), &
               (Frame%Coor%BOUND%YL+Frame%Coor%BOUND%YR)/2+Frame%Font%S )
       endif
       CALL PSPrtTexAtLeft(Frame,ylabel) 
    endif
  END subroutine PsLabelY

  subroutine PSMoveToRatio(Frame,X,Y,LINES)
    TYPE(PSFrame)Frame
    real(dl)X,Y
    integer,optional::LINES
    IF(PRESENT(LINES))then
       CALL PsMoveTo(Frame,NINT(Frame%Coor%BOUND%XL*(1.-X)+Frame%Coor%BOUND%XR*X), &
            NINT(Frame%Coor%BOUND%YL*(1.-Y)+Frame%Coor%BOUND%YR*Y)-(Frame%Font%S+10)*Lines)
       LINES=LINES+1
    ELSE
       CALL PsMoveTo(Frame,NINT(Frame%Coor%BOUND%XL*(1.-X)+Frame%Coor%BOUND%XR*X), &
            NINT(Frame%Coor%BOUND%YL*(1.-Y)+Frame%Coor%BOUND%YR*Y))
    endif
  END subroutine PSMoveToRatio

  subroutine PSDefaultCoor_D(Frame,XL,XR,YL,YR,notdrawnow)
    TYPE(PSFrame)Frame
    real(dl)XL,XR,YL,YR
    integer I
    logical,optional::notdrawnow
    IF(Frame%FileUnit.le.0 .or. Frame%FileUnit.ge.100.or.trim(Frame%FileName).eq."")CALL PSStart(Frame)
    IF(Frame%Coor%Locked)RETURN
    if(xr.le.xl)then
       xl = xl -0.5
       xr = xl +1
    endif
    if(yr.le.yl) then
       yl = yl-0.5
       yr =yl+1.
    endif

    CALL PSGetCoor(Frame%Coor%CX,Frame%Coor%CXTYPE,Frame%Coor%XPOW,Frame%Coor%XFORM,XL,XR,Frame%Coor%XLOG,Frame%Coor%XTicks)
    Select case(FRAME%COOR%XPOW)
    Case(3)
       Frame%Coor%CX=Frame%Coor%CX*1000.
       Frame%Coor%XPOW=0
    Case(2)
       Frame%Coor%CX=Frame%Coor%CX*100.
       Frame%Coor%XPOW=0
    case(1)
       Frame%Coor%CX=Frame%Coor%CX*10.
       Frame%Coor%XPOW=0
    case(-1)
       Frame%Coor%CX=Frame%Coor%CX/10.
       Frame%Coor%XPOW=0
       if(trim(Frame%Coor%xform) .eq. "")Frame%Coor%XFORM="f"
     case(-2)
       Frame%Coor%CX=Frame%Coor%CX/100.
       Frame%Coor%XPOW=0
       if(trim(Frame%Coor%xform) .eq. "")Frame%Coor%XFORM="f"
    case(-3)
       Frame%Coor%CX=Frame%Coor%CX/1000.
       Frame%Coor%XPOW=0
       if(trim(Frame%Coor%xform) .eq. "")Frame%Coor%XFORM="f"
    case(-4)
       Frame%Coor%CX=Frame%Coor%CX/10000.
       Frame%Coor%XPOW=0
       if(trim(Frame%Coor%xform) .eq. "")Frame%Coor%XFORM="f"
    END Select

    CALL PSGetCoor(Frame%Coor%CY,Frame%Coor%CYTYPE,Frame%Coor%YPOW,Frame%Coor%YFORM,YL,YR,Frame%Coor%YLOG,Frame%Coor%YTicks)

    Select case (Frame%coor%YPOW)
    case(3)
       Frame%Coor%CY=Frame%Coor%CY*1000.
       Frame%Coor%YPOW=0
    case(2)
       Frame%Coor%CY=Frame%Coor%CY*100.
       Frame%Coor%YPOW=0
    case(1)
       Frame%Coor%CY=Frame%Coor%CY*10.
       Frame%Coor%YPOW=0
    case(-1)
       Frame%Coor%CY=Frame%Coor%CY/10.
       Frame%Coor%YPOW=0
       if(trim(Frame%Coor%yform) .eq. "")Frame%Coor%YFORM="f"
    case(-2)
       Frame%Coor%CY=Frame%Coor%CY/100.
       Frame%Coor%YPOW=0
       if(trim(Frame%Coor%yform) .eq. "")Frame%Coor%YFORM="f"
    case(-3)
       Frame%Coor%CY=Frame%Coor%CY/1000.
       Frame%Coor%YPOW=0
       if(trim(Frame%Coor%yform) .eq. "")Frame%Coor%YFORM="f"
    case(-4)
       Frame%Coor%CY=Frame%Coor%CY/10000.
       Frame%Coor%YPOW=0
       if(trim(Frame%Coor%yform) .eq. "")Frame%Coor%YFORM="f"
    ENd select
    IF(Frame%Coor%FracBoundX)then
       Frame%Coor%OX=XL
       Frame%Coor%XUR=xR
    else
       Frame%Coor%OX=MIN(XL,Frame%Coor%CX(1)*10.**Frame%Coor%XPOW)
       I=0
       DO WHILE(Frame%Coor%CXTYPE(I+1).GT.0 .AND. Frame%Coor%CXTYPE(I+1).LE.7 .and. I.lt.PS_N_CoorLabels-1)
          I=I+1
       enddo
       if(i.gt.0)then
          Frame%Coor%XUR=MAX(XR,Frame%Coor%CX(I)*10.**Frame%Coor%XPOW)
       else
          Frame%Coor%XUR=XR
       endif
    endif
    Frame%Coor%XSCALE=(Frame%Coor%XUR-Frame%Coor%OX)/(Frame%Coor%BOUND%XR-Frame%Coor%BOUND%XL)
    IF(Frame%Coor%FracBoundY)then
       Frame%Coor%OY=YL
       Frame%Coor%YUR=YR
    else
       Frame%Coor%OY=MIN(YL,Frame%Coor%CY(1)*10.**Frame%Coor%YPOW)
       I=0
       DO WHILE(Frame%Coor%CYTYPE(I+1).GT.0 .AND. Frame%Coor%CYTYPE(I+1).LE.7 .and. I.lt.PS_N_CoorLabels-1)
          I=I+1
       enddo
       IF(I.GT.0)then
          Frame%Coor%YUR=MAX(YR,Frame%Coor%CY(I)*10.**Frame%Coor%YPOW)
       ELSE
          fRAME%COOR%yur=yr
       endif
    endif
    Frame%Coor%YSCALE=(Frame%Coor%YUR-Frame%Coor%OY)/(Frame%Coor%BOUND%YR-Frame%Coor%BOUND%YL)

    if(present(notdrawnow))then
       if(.not.notdrawnow)then
          CALL PSDrawCoor(Frame)
       endif
    else
       CALL PSDrawCoor(Frame)
    endif
    CALL PsLabelX(FRAME)
    CALL PsLabelY(FRAME)
  END subroutine PSDefaultCoor_D



  subroutine PSDefaultCoor_DARR(Frame,X,Y)
    TYPE(PSFrame)Frame
    real(dl),dimension(:),INTENT(IN)::X,Y
    real(dl) XL,XR,YL,YR
    XL=MINVAL(X)
    XR=MAXVAL(X)
    YL=MINVAL(Y)
    YR=MAXVAL(Y)
    CALL PSDefaultCoor_D(Frame,xl-(xr-xl)*0.05 , xr + (xr-xl)*0.05 ,YL-(yr-yl)*0.05,YR+(yr-yl)*0.05)
  END subroutine PSDefaultCoor_DARR

  subroutine PSDefaultCoor_DARR1(Frame,xmin,xmax,Y)
    TYPE(PSFrame)Frame
    real(dl),dimension(:),INTENT(IN)::Y
    real(dl) Xmin,xmax,YL,YR
    YL=MINVAL(Y)
    YR=MAXVAL(Y)
    CALL PSDefaultCoor_D(Frame,xmin, xmax ,YL-(yr-yl)*0.05,YR+(yr-yl)*0.05)
  END subroutine PSDefaultCoor_DARR1

  subroutine PSDefaultCoor_DARR2(Frame,X,Ymin,ymax)
    TYPE(PSFrame)Frame
    real(dl),dimension(:),INTENT(IN)::x
    real(dl) XL,XR,ymin,ymax
    XL=MINVAL(x)
    XR=MAXVAL(x)
    CALL PSDefaultCoor_D(Frame,xl-(xr-xl)*0.05 , xr + (xr-xl)*0.05 , ymin, ymax)
  END subroutine PSDefaultCoor_DARR2


  subroutine PSDrawCoor(Frame)
    TYPE(PSFrame)Frame
    integer I,J
    STRING CurrentType
    CALL PSGSAVE(Frame)
    CurrentType=trim(Frame%Font%T)
    Frame%Font%T="Symbol"
    WRITE(Frame%FileUnit,'(a)')trim(INT2STR(Frame%Font%S))//" gft"
    CALL PSOUTLINEBOX(Frame,Frame%Coor%BOUND)


    if(frame%coor%plotxticks)then
       DO I=1,PS_N_COORLABELS
          IF(Frame%Coor%CXTYPE(I).LE.1.e-10 .OR. Frame%Coor%CXTYPE(I).GE.7)then
             IF(Frame%Coor%CXTYPE(I).EQ.7)then
                CYCLE
             ELSE
                EXIT
             endif
          ELSE
             IF(Frame%Coor%FracBoundX.and.(Frame%Coor%CX(I).lt.Frame%Coor%OX.or.Frame%Coor%CX(I).gt.Frame%Coor%XUR))CYCLE
             CALL PSCMoveTo(Frame,Frame%Coor%CX(I)*10.**Frame%Coor%XPOW,Frame%Coor%OY)
             IF(Frame%Coor%CXTYPE(I).LE.3)then
                CALL PSRLineTo(Frame,0,Frame%Coor%LongTick)
             ELSE
                CALL PSRLineTo(Frame,0,Frame%Coor%ShortTick)
             endif
             IF(Frame%Coor%UR)then
                CALL PSCMoveTo(Frame,Frame%Coor%CX(I)*10.**Frame%Coor%XPOW,Frame%Coor%YUR)
                IF(Frame%Coor%CXTYPE(I).LE.3)then
                   CALL PSRLineTo(Frame,0,-Frame%Coor%LongTick)
                ELSE
                   CALL PSRLineTo(Frame,0,-Frame%Coor%ShortTick)
                endif
             endif

             if(trim(frame%coor%xform) .eq. "none") cycle
             J=MOD(Frame%Coor%CXTYPE(I),3)
             IF(J.ne.1)then
                CALL PSCMoveTo(Frame,Frame%Coor%CX(I)*10.**Frame%Coor%XPOW,Frame%Coor%OY)
                CALL PSRMOVETO(Frame,0,-Frame%Font%S-Frame%Coor%LongTick)
                IF(J.EQ.2)then
                   IF(Frame%Coor%XPOW.EQ.0 .OR. ABS(Frame%Coor%CX(I)).LT.1.D-5)then
                      IF(Frame%Coor%XFORM.EQ."" .OR. ABS(Frame%Coor%CX(I)).LT.1.D-5)then
                         WRITE(Frame%FileUnit,'(a)') "("// &
                              trim(INT2STR(NINT(Frame%Coor%CX(I))))//") prtcc"
                      ELSE
                         WRITE(Frame%FileUnit,'(a)') "("//  &
                              trim(NUM2STR(Frame%Coor%CX(I),Frame%Coor%XForm))//") prtcc"
                      endif
                   ELSE
                      WRITE(Frame%FileUnit,'(a)') "("//  &
                           trim(NUM2STR(Frame%Coor%CX(I),Frame%Coor%XFORM))//"\26410) prtcc"
                      CALL PSSUPS(Frame,trim(INT2STR(Frame%Coor%XPOW)))
                   endif
                ELSE
                   IF(ABS(Frame%Coor%CX(I)).GT.3.01)then
                      write(frame%fileunit,'(a)')"(10) prtcc"
                      if(NINT(Frame%Coor%CX(I)).ne.1)then
                         CALL PSSUPS(Frame,trim(INT2STR(NINT(Frame%Coor%CX(I)))))
                      endif
                   ELSE
                      Select case(Nint(Frame%Coor%CX(i)))
                      case(0)
                         write(frame%fileunit,'(a)') "(1) prtcc"
                      case(1)
                         write(frame%fileunit,'(a)') "(10) prtcc"
                      case(2)
                         write(frame%fileunit,'(a)') "(100) prtcc"                         
                      case(3)
                         write(frame%fileunit,'(a)') "(1000) prtcc"                         
                      case(-1)
                         write(frame%fileunit,'(a)') "(0.1) prtcc"                        
                      case(-2)
                         write(frame%fileunit,'(a)') "(0.01) prtcc"                        
                      case(-3)
                         write(frame%fileunit,'(a)') "(0.001) prtcc"                        
                      end Select
                   endif
                endif
             endif
          endif
       enddo
    endif
    if(frame%coor%plotyticks)then
       DO I=1,PS_N_COORLABELS
          IF(Frame%Coor%CYTYPE(I).LE.0 .OR. Frame%Coor%CYTYPE(I).GE. 7)then
             IF(Frame%Coor%CYTYPE(I).EQ.7)then
                CYCLE
             ELSE
                EXIT
             endif
          ELSE
             IF(Frame%Coor%FracBoundY.and.(Frame%Coor%CY(I).lt.Frame%Coor%OY.or.Frame%Coor%CY(I).gt.Frame%Coor%YUR))CYCLE
             CALL PSCMoveTo(Frame,Frame%Coor%OX,Frame%Coor%CY(I)*10.**Frame%Coor%YPOW)
             IF(Frame%Coor%CYTYPE(I).LE.3)then
                CALL PSRLineTo(Frame,Frame%Coor%LongTick,0)
             ELSE
                CALL PSRLineTo(Frame,Frame%Coor%ShortTick,0)
             endif
             IF(Frame%Coor%UR)then
                CALL PSCMoveTo(Frame,Frame%Coor%XUR,Frame%Coor%CY(I)*10.**Frame%Coor%YPOW)
                IF(Frame%Coor%CYTYPE(I).LE.3)then
                   CALL PSRLineTo(Frame,-Frame%Coor%LongTick,0)
                ELSE
                   CALL PSRLineTo(Frame,-Frame%Coor%ShortTick,0)
                endif
             endif
             if(trim(frame%coor%yform) .eq. "none") cycle
             J=MOD(Frame%Coor%CYTYPE(I),3)
             IF(J.ne.1)then
                CALL PSCMoveTo(Frame,Frame%Coor%OX,Frame%Coor%CY(I)*10.**Frame%Coor%YPOW)
                CALL PSRMOVETO(Frame,-Frame%Coor%LongTick,-2)
                IF(J.EQ.2)then
                   IF(Frame%Coor%YPOW.EQ.0 .OR. ABS(Frame%Coor%CY(I)).LT.1.D-3)then
                      write(frame%fileunit,'(a)') "("//  &
                           trim(NUM2STR(Frame%Coor%CY(I),Frame%Coor%YFORM))//") prtcr"
                   ELSE
                      CALL PSRMOVETO(Frame,-MIN(Frame%Coor%BOUND%XL/3,Frame%Font%S* &
                           LEN_trim(INT2STR(Frame%Coor%YPOW))/3+Frame%Coor%LongTick),0)
                      IF(Frame%Coor%YFORM.EQ."")then
                         write(frame%fileunit,'(a)') "("// &
                              trim(INT2STR(NINT(Frame%Coor%CY(I))))//"\26410) prtcr"
                      ELSE
                         write(frame%fileunit,'(a)') "("//  &
                              trim(NUM2STR(Frame%Coor%CY(I),Frame%Coor%YFORM))//"\26410) prtcr"
                      endif
                      CALL PSSUPS(Frame,trim(INT2STR(Frame%Coor%YPOW)))
                   endif
                ELSE
                   IF(ABS(Frame%Coor%CY(I)).GT.3.01d0)then
                      CALL PSRMOVETO(Frame,-MIN(Frame%Coor%BOUND%XL/3,Frame%Font%S* &
                           LEN_trim(INT2STR(NINT(Frame%Coor%CY(I))))/3+Frame%Coor%LongTick),0)
                      write(frame%fileunit,'(a)') "(10) prtcr"
                      if(NINT(Frame%Coor%CY(I)).ne.1)then
                         CALL PSSUPS(Frame,trim(INT2STR(NINT(Frame%Coor%CY(I)))))
                      endif
                   ELSE
                      Select case(Nint(Frame%Coor%cy(i)))
                      case(0)
                         write(frame%fileunit,'(a)') "(1) prtcr"
                      case(1)
                         write(frame%fileunit,'(a)') "(10) prtcr"
                      case(2)
                         write(frame%fileunit,'(a)') "(100) prtcr"                         
                      case(3)
                         write(frame%fileunit,'(a)') "(1000) prtcr"                         
                      case(-1)
                         write(frame%fileunit,'(a)') "(0.1) prtcr"                        
                      case(-2)
                         write(frame%fileunit,'(a)') "(0.01) prtcr"                        
                      case(-3)
                         write(frame%fileunit,'(a)') "(0.001) prtcr"                        
                      end Select
                   endif
                endif
             endif
          endif
       enddo
    endif
    CALL PSSTROKE(Frame)
    CALL PSsetfont(Frame,CurrentType)
    CALL PSGRESTORE(Frame)
  END subroutine PSDrawCoor

  subroutine PSREDRAWCOOR(FRAME)
    TYPE(PSFrame)Frame
    integer I
    CALL PSGSAVE(Frame)
    CALL PSOUTLINEBOX(Frame,Frame%Coor%BOUND)
    DO I=1,PS_N_COORLABELS
       IF(Frame%Coor%CXTYPE(I).LE.0 .OR. Frame%Coor%CXTYPE(I).GE. 7)then
          IF(Frame%Coor%CXType(I).EQ.7)then
             CYCLE
          ELSE
             EXIT
          endif
       ELSE
          IF(Frame%Coor%FracBoundX.and.(Frame%Coor%CX(I).lt.Frame%Coor%OX.or.Frame%Coor%CX(I).gt.Frame%Coor%XUR))CYCLE
          CALL PSCMoveTo(Frame,Frame%Coor%CX(I)*10.**Frame%Coor%XPOW,Frame%Coor%OY)
          IF(Frame%Coor%CXTYPE(I).LE.3)then
             CALL PSRLineTo(Frame,0,Frame%Coor%LongTick)
          ELSE
             CALL PSRLineTo(Frame,0,Frame%Coor%ShortTick)
          endif
          IF(Frame%Coor%UR)then
             CALL PSCMoveTo(Frame,Frame%Coor%CX(I)*10.**Frame%Coor%XPOW,Frame%Coor%YUR)
             IF(Frame%Coor%CXTYPE(I).LE.3)then
                CALL PSRLineTo(Frame,0,-Frame%Coor%LongTick)
             ELSE
                CALL PSRLineTo(Frame,0,-Frame%Coor%ShortTick)
             endif
          endif
       endif
    enddo
    DO I=1,PS_N_COORLABELS
       IF(Frame%Coor%CYTYPE(I).LE.0 .OR. Frame%Coor%CYTYPE(I).GE. 7)then
          IF(FRAME%COOR%CYTYPE(I).EQ.7)then
             CYCLE
          ELSE
             EXIT
          endif
       ELSE
          IF(Frame%Coor%FracBoundY.and.(Frame%Coor%CY(I).lt.Frame%Coor%OY.or.Frame%Coor%CY(I).gt.Frame%Coor%YUR))CYCLE
          CALL PSCMoveTo(Frame,Frame%Coor%OX,Frame%Coor%CY(I)*10.**Frame%Coor%YPOW)
          IF(Frame%Coor%CYTYPE(I).LE.3)then
             CALL PSRLineTo(Frame,Frame%Coor%LongTick,0)
          ELSE
             CALL PSRLineTo(Frame,Frame%Coor%ShortTick,0)
          endif
          IF(Frame%Coor%UR)then
             CALL PSCMoveTo(Frame,Frame%Coor%XUR,Frame%Coor%CY(I)*10.**Frame%Coor%YPOW)
             IF(Frame%Coor%CYTYPE(I).LE.3)then
                CALL PSRLineTo(Frame,-Frame%Coor%LongTick,0)
             ELSE
                CALL PSRLineTo(Frame,-Frame%Coor%ShortTick,0)
             endif
          endif
       endif
    enddo
    CALL PSSTROKE(Frame)
    CALL PSGRESTORE(Frame)
  END subroutine PSREDRAWCOOR
  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine PSNoCoorCurve(Frame,X,Y)
    type(psframe)frame
    real(dl),dimension(:),intent(in)::X,Y
    integer(IB) i,n
    N=GetDim("PSNoCoorCurve",SIZE(X),Size(Y))
    call psnewpath(frame)
    call psmoveto(frame,x(1),y(1))
    do i=2,N
       call pslineto(frame,x(i),y(i))
    enddo
    call psstroke(frame)
  End subroutine PSNoCoorCurve

  subroutine PSPlotCurve_D(Frame,X,Y, stroke_color,fill_color)
    UNKNOWN_STRING  ,optional::stroke_color, fill_color
    STRING str,laststr
    TYPE(PSFrame)Frame
    real(dl),dimension(:),INTENT(IN)::X,Y
    real(dl) X1,Y1,accuracy
    integer(IB) ix,iy,ix2,iy2
    integer N,i,ns
    type(PSPath)path
    logical stroked, filled
    N=GetDim("PSplotCurve",SIZE(X),Size(Y))
    if(n.le.1) return
    CALL PSGSAVE(Frame)
    CALL PSCLIPCOOR(FRAME)
    filled =  .false.
    stroked = .false.
50  CALL PsNewPath(Frame)
    if(ps_resolution_level.gt.0)then 
       path%n=0
       accuracy = min( abs(Frame%Coor%Bound%XR-Frame%Coor%Bound%XL), abs(Frame%Coor%Bound%YR-Frame%Coor%Bound%YL) )/10.**(ps_resolution_level+3) !!less than 0.1% is almost invisible by eye
       laststr=''
       do i=1,N
          call pscoortrans(Frame,x(i),y(i),x1,y1)
          str=trim(PS_Num_String(x1))//' '//trim(PS_Num_String(y1))
          if(trim(str).ne.trim(laststr))then
             call add2path(path,x1,y1,accuracy)
             laststr=trim(str)
          endif
       enddo
       write(frame%fileunit,'(a)') trim(PS_Num_String(path%x(1)))//' '//trim(PS_Num_String(path%y(1)))//' m'
       do i=path%n,2,-1
          write(frame%fileunit,'(a)') trim(PS_Num_String(path%x(i)))//' '//trim(PS_Num_String(path%y(i)))
       enddo
       write(frame%fileunit,'(a)') trim(int2str(path%n-1))//' {l} repeat'
    else
       ns=0
       call pscoortrans(frame,x(1),y(1),ix,iy)
       write(frame%fileunit,'(a)') trim(Int2Str(ix))//' '//trim(Int2Str(iy))//' m'
       call pscoortrans(frame,x(N),y(N),ix,iy)
       do i=N-1,1,-1
          call pscoortrans(frame,x(i),y(i),ix2,iy2)
          if(ix2.ne.ix .or. iy2.ne.iy)then
             write(frame%fileunit,'(a)') Trim(PSAbbrev(trim(Int2Str(ix-ix2))//' '//trim(Int2Str(iy-iy2))))
             ix=ix2
             iy=iy2
             ns=ns+1
          endif
       enddo
       write(frame%fileunit,'(a)') trim(int2str(ns))//' {r} repeat'
    endif

    if(.not. filled)then
       if(present(fill_color))then
          if(trim(fill_color).ne."")then
             call pssetcolor(Frame,fill_color)
             call PSClosePath(frame)
             call Psfill(frame)
             filled = .true.
             if(.not. present(stroke_color))then
                goto 50
             else
                if(trim(stroke_color) .ne. "") goto 50
             endif
          endif
       endif
    endif
    if(.not. stroked)then
       if(present(stroke_color))then
          if(trim(stroke_color).ne."")then
             call pssetcolor(frame, stroke_color)
             call PsStroke(Frame)
             stroked = .true.
          endif
       else
          call PsStroke(Frame)
          stroked = .true.
       endif
    endif
    CALL PSGRESTORE(Frame)
  END subroutine PSPlotCurve_D

  subroutine PSSortPlotData(Frame, filename, colx, coly)
    type(PSFrame) Frame
    UNKNOWN_STRING   filename
    integer colx, coly
    real(dl),dimension(:),allocatable::x, y, tmp
    integer n, i
    type(file_pointer)fp
    if(colx.lt.1 .or. coly.lt.1) call ReturnError("PSSortPlotData","negative column number", "stop")
    n = file_NumLines(filename)
    if(n.le.0)then
       write(*,*) "warning: cannot find data in "//trim(filename)
       return
    endif
    allocate(x(n),y(n), tmp(max(colx, coly)))
    fp = open_file(filename, "r")
    do i=1, n
       if(.not. file_readOneLine(fp, tmp)) call ReturnError("PSSortPlotData","data format error or not enough columns in "//trim(filename), "stop")
       x(i) = tmp(colx)
       y(i) = tmp(coly)
    enddo
    call close_file(fp)
    call quicksortAcc(x,y)
    call psplotcurve_d(frame,x,y)
    deallocate(x,y,tmp)
  End subroutine PSSortPlotData

  Subroutine PSPlotDenseCurve_D(Frame,X,Y,N,Interp_method) !!0 linear; 1 natural cubic spline; 2 chebyshev; 3 monotone cubic hermite
    type(PSFrame)Frame
    real(dl),dimension(:),intent(IN)::x,y
    integer(IB) N
    integer(IB),optional::Interp_method
    integer(IB) intm
    real(dl),dimension(:),allocatable::xbin,ybin,y2
    if(size(x).ne.size(y))stop "PSPlotDenseCurve: size of arrays do not agree."
    if(present(interp_method))then
       intm = interp_method
    else
       intm = 3
    endif
    select case(intm)
    case(Cubicspline_Interpolation_Option)
       allocate(xbin(N),ybin(N), y2(size(y)))
       call findgen(xbin,minval(x),maxval(x))
       call splines(x,y,y2)
       call splints(x,y,y2,xbin,ybin)
       call psplotcurve_D(Frame,xbin,ybin)
       deallocate(xbin,ybin,y2)
    case(Chebyshev_Interpolation_Option) !!In this option actually only the boundaries from x are used. I always assume y(:) are the values on the nodal points.
       if(size(y).gt.15)stop "Chebyshev Interpolation beyond 15 nodal points is not recommended. You may want to use cubic spline or hermite interpolations."
       allocate(xbin(N),ybin(N), y2(size(y)))       
       call findgen(xbin,minval(x),maxval(x))
       y2=matmul(Mat_ChebTrans(size(y)),y)  !!get Chebyshev coefficients       
       ybin = Chebfit_value( xbin(1), xbin(N), y2(1:size(y)), xbin(1:N))
       call psplotcurve_D(Frame,xbin,ybin)
       deallocate(xbin,ybin,y2)       
    case(Hermite_Interpolation_Option)
       allocate(xbin(N),ybin(N))       
       call findgen(xbin,minval(x),maxval(x))
       call Hermite_Interp_Monotone(size(x),x,y,xbin,ybin)
       call psplotcurve_D(Frame,xbin,ybin)
       deallocate(xbin,ybin) 
    case default
       call psplotcurve_d(frame,x,y)
    end select
  end Subroutine PSPlotDenseCurve_D


  subroutine PsPlotCurve_FromFile(frame,filename)
    STRING str,laststr
    TYPE(PSFrame)Frame
    real(dl) x,y,x1,y1,accuracy
    integer i
    type(file_pointer)fp
    UNKNOWN_STRING   filename
    type(PSPath)path
    fp = open_file(filename,"r")
    CALL PSGSAVE(Frame)
    CALL PSCLIPCOOR(FRAME)
    CALL PsNewPath(Frame)
    path%n=0
    accuracy = min( abs(Frame%Coor%Bound%XR-Frame%Coor%Bound%XL), abs(Frame%Coor%Bound%YR-Frame%Coor%Bound%YL) )/10.**(ps_resolution_level+3) !!less than 0.1% is almost invisible by eye
    laststr=''
    do 
       read(fp%unit,*, END=50) x, y
       call pscoortrans(Frame,x,y,x1,y1)
       str=trim(PS_Num_String(x1))//' '//trim(PS_Num_String(y1))
       if(trim(str).ne.trim(laststr))then
          call add2path(path,x1,y1,accuracy)
          laststr=trim(str)
       endif
    enddo
50  write(frame%fileunit,'(a)') trim(PS_Num_String(path%x(1)))//' '//trim(PS_Num_String(path%y(1)))//' m'
    do i=path%n,2,-1
       write(frame%fileunit,'(a)') trim(PS_Num_String(path%x(i)))//' '//trim(PS_Num_String(path%y(i)))
    enddo
    write(frame%fileunit,'(a)') trim(int2str(path%n-1))//' {l} repeat'
    CALL PSSTROKE(Frame)
    CALL PSGRESTORE(Frame)
    call close_file(fp)
  end subroutine PsPlotCurve_FromFile


  subroutine PSPlotSmoothCurve(Frame,X,Y,Smoothorder) 
    integer,PARAMETER::M=12
    TYPE(PSFRAME)Frame
    Integer(IB),optional::Smoothorder
    !! 1<=n<=10, n-th order chebyshev polynomial fitting; 0 spline;
    real(dl),dimension(:),Intent(IN)::X,Y
    real(dl),dimension(:),allocatable::Y2
    real(dl) XS(M),YS(M),XX(M*30),YY(M*30),a,b,C(20),tmp
    integer n
    n=GetDim("plotsmoothCurve",SIZE(X),size(y))
    If(Present(Smoothorder).and.Smoothorder.gt.0)then
       If(Smoothorder.gt.12)then
          Stop "PlotSmoothCurve (using chebyshev): n>12 smooth order not supported."
       EndIf
       a=Minval(x)
       b=maxval(x)
       tmp=(b-a)/30.
       CALL FINDGEN(XX,a,b)
       a=a-tmp
       b=b+tmp
       Call ChebFit(X,Y,a,b,C(1:Smoothorder))
       YY=ChebFit_value(a,b,c(1:smoothorder),XX)
    else !! cubic spline fitting
       CALL FINDGEN(XS,X(1),X(N))
       CALL FINDGEN(XX,X(1),X(N))
       allocate(Y2(N))
       CALL SPLINES(X(1:N),Y(1:N),Y2(1:N))
       CALL SPLINTS(X(1:N),Y(1:N),Y2(1:N),XS,YS)
       CALL SPLINTS(X(1:N),Y(1:N),Y2(1:N),XX,YY)
       deallocate(Y2)
    endif
    CALL PSPLOTCURVE(FRAME,XX,YY)
  END subroutine PSPlotSmoothCurve


  !!=============PLOTYEB=====================!!
  Subroutine PSPlotYbar(Frame,X,Y,YEB)
    Type(PSFrame)Frame
    Real(Dl),INTENT(IN)::X,Y,YEB
    Integer(IB) ix,iy,iyeb
    Call PSGsave(Frame)
    CALL PSCOORTRANS(Frame,X,Y,IX,IY)
    IYEB=NINT(YEB/Frame%Coor%YSCALE/2.)
    write(frame%fileunit,'(a)') trim(INT2STR(IYEB))// &
         " "//trim(INT2STR(IX))//" "//trim(INT2STR(IY))//" psysbar"
    Call PSSTRoke(Frame)
    Call PSGRESTORE(Frame)
  End subroutine PSPlotYbar
    


  subroutine PSPLOTYEB(Frame,X,Y,YEB)
    TYPE(PSFrame)Frame
    real(dl),dimension(:),INTENT(IN)::X,Y,YEB
    integer I,N,IX,IY,IYEB
    N=SIZE(X)
    CALL PSGSAVE(Frame)
    CALL PsNewPath(Frame)
    CALL PSOUTLINEBOX(Frame,Frame%Coor%BOUND)
    CALL PsPush(Frame,"clip")
    CALL PsNewPath(Frame)
    DO I=1,N
       CALL PSCOORTRANS(Frame,X(I),Y(I),IX,IY)
       IYEB=NINT(YEB(I)/Frame%Coor%YSCALE)
       IF(IYEB.GE.4)then
          write(frame%fileunit,'(a)') trim(INT2STR(IYEB))//  &
               " "//trim(INT2STR(IX))//" "//trim(INT2STR(IY))//" psybar"
       ELSE
          write(frame%fileunit,'(a)') trim(INT2STR(IYEB))//  &
               " "//trim(INT2STR(IX))//" "//trim(INT2STR(IY))//" psysbar"
       endif
    enddo
    CALL PSSTROKE(Frame)
    CALL PSGRESTORE(Frame)
  END subroutine PSPLOTYEB


  subroutine PSPLOTXEB(Frame,X,Y,XEB)
    TYPE(PSFrame)Frame
    real(dl),dimension(:),INTENT(IN)::X,Y,XEB
    integer I,N,IX,IY,IXEB
    N=SIZE(X)
    CALL PSGSAVE(Frame)
    CALL PsNewPath(Frame)
    CALL PSOUTLINEBOX(Frame,Frame%Coor%BOUND)
    CALL PsPush(Frame,"clip")
    CALL PsNewPath(Frame)
    DO I=1,N
       CALL PSCOORTRANS(Frame,X(I),Y(I),IX,IY)
       IXEB=nint(XEB(I)/Frame%Coor%XSCALE)
       IF(IXEB.GE.4)then
          write(frame%fileunit,'(a)') trim(INT2STR(IXEB))// &
               " "//trim(INT2STR(IX))//" "//trim(INT2STR(IY))//" psxbar"
       ELSE
          write(frame%fileunit,'(a)') trim(INT2STR(IXEB)), &
               " "//trim(INT2STR(IX))//" "//trim(INT2STR(IY))//" psxsbar"
       endif
    enddo
    CALL PSSTROKE(Frame)
    CALL PSGRESTORE(Frame)
  END subroutine PSPLOTXEB


  subroutine PSPLOTCEB(Frame,X,Y,XEB,YEB)
    TYPE(PSFrame)Frame
    real(dl),dimension(:),INTENT(IN)::X,Y,XEB,YEB
    integer I,N,IX,IY,IXEB,IYEB
    N=SIZE(X)
    CALL PSGSAVE(Frame)
    CALL PsNewPath(Frame)
    CALL PSOUTLINEBOX(Frame,Frame%Coor%BOUND)
    CALL PsPush(Frame,"clip")
    CALL PsNewPath(Frame)
    DO I=1,N
       CALL PSCOORTRANS(Frame,X(I),Y(I),IX,IY)
       IXEB=nint(XEB(I)/Frame%Coor%XSCALE)
       IYEB=nint(YEB(I)/Frame%Coor%YSCALE)
       IF(IXEB.GE.4 .AND. IYEB.GE.4)then
          write(frame%fileunit,'(a)') trim(INT2STR(IXEB))//  &
               " "//trim(INT2STR(IYEB))//" "//trim(INT2STR(IX))//" "//trim(INT2STR(IY))//" pscbar"
       ELSE
          write(frame%fileunit,'(a)') trim(INT2STR(IXEB))//  &
               " "//trim(INT2STR(IYEB))//" "//trim(INT2STR(IX))//" "//trim(INT2STR(IY))//" pscsbar"
       endif
    enddo
    CALL PSSTROKE(Frame)
    CALL PSGRESTORE(Frame)
  END subroutine PSPLOTCEB

  subroutine PSDrawArrow(Frame,X,Y,LENGTH,A,HEADL,HEADA)
    TYPE(PSFrame)Frame
    integer X,Y,LENGTH,A,HEADL,HEADA
    CALL PSGSAVE(Frame)
    CALL PsNewPath(Frame)
    CALL PSTRANSLATE(Frame,X,Y)
    CALL PSRotate(Frame,A)
    CALL PsMoveTo(Frame,0,0)
    CALL PsLineTO(Frame,LENGTH,0)
    CALL PSTRANSLATE(Frame,LENGTH,0)
    CALL PSROTATE(Frame,180-HEADA)
    CALL PsLineTO(Frame,HEADL,0)
    CALL PsMoveTo(Frame,0,0)
    CALL PSROTATE(Frame,HEADA*2)
    CALL PsLineTO(Frame,HEADL,0)
    CALL PSSTROKE(Frame)
    CALL PSGRESTORE(Frame)
  END subroutine PSDrawArrow

  subroutine PSDrawSolidArrow(Frame,X,Y,LENGTH,A,HEADL,HEADA)
    TYPE(PSFrame)Frame
    integer X,Y,LENGTH,A,HEADL,HEADA
    CALL PSGSAVE(Frame)
    CALL PsNewPath(Frame)
    CALL PSTRANSLATE(Frame,X,Y)
    CALL PSROTATE(Frame,A)
    CALL PsMoveTo(Frame,0,0)
    CALL PsLineTO(Frame,LENGTH,0)
    CALL PSSTROKE(Frame)
    CALL PSTRANSLATE(Frame,LENGTH,0)
    CALL PSROTATE(Frame,180-HEADA)
    CALL PsNewPath(Frame)
    CALL PsMoveTo(Frame,HEADL,0)
    CALL PsLineTO(Frame,0,0)
    CALL PSROTATE(Frame,HEADA*2)
    CALL PsLineTO(Frame,HEADL,0)
    CALL PSCLOSEPATH(Frame)
    CALL PSFILL(Frame)
    CALL PSGRESTORE(Frame)
  END subroutine PSDrawSolidArrow


  subroutine PSPlotArrow(Frame,X,Y,LENGTH,DX,DY,HEADL,HEADA)
    TYPE(PSFrame)Frame
    integer LENGTH,HEADL,A,HEADA
    real(dl) DX,DY
    real(dl) X,Y
    integer IX,IY
    CALL PSCOORTRANS(Frame,X,Y,IX,IY)
    A=NINT(ATAN2(DY/Frame%Coor%YSCALE,DX/Frame%Coor%XSCALE)*180._dl/const_PI)
    CALL PSDrawArrow(Frame,IX,IY,LENGTH,A,HEADL,HEADA)
  END subroutine PSPlotArrow


  subroutine PSPlotSolidArrow(Frame,X,Y,LENGTH,DX,DY,HEADL,HEADA)
    TYPE(PSFrame)Frame
    integer LENGTH,HEADL,A,HEADA
    real(dl) DX,DY
    real(dl) X,Y
    integer IX,IY
    CALL PSCOORTRANS(Frame,X,Y,IX,IY)
    A=NINT(ATAN2(DY/Frame%Coor%YSCALE,DX/Frame%Coor%XSCALE)*180./const_PI)
    CALL PSDrawSolidArrow(Frame,IX,IY,LENGTH,A,HEADL,HEADA)
  END subroutine PSPlotSolidArrow

  Subroutine PSPlotLine_line(Frame,Line)
    Type(PSFrame)Frame
    Type(PSLine) Line
    call PSPlotLine_D(Frame,Line%xmin,Line%ymin,Line%xmax,Line%ymax)
  ENd Subroutine PSPlotLine_line

  subroutine PSPlotLine_D(Frame,X1,Y1,X2,Y2)
    TYPE(PSFrame)Frame
    real(dl) X1,Y1,X2,Y2
    CALL PSGSAVE(Frame)
    CALL PsNewPath(Frame)
    CALL PsCMoveTo(Frame,X1,Y1)
    CALL PSCLINETO(Frame,X2,Y2)
    CALL PSSTROKE(Frame)
    CALL PSGRESTORE(Frame)
  END subroutine PSPlotLine_D


  subroutine PSPlotLine_I(Frame,X,color)
    TYPE(PSFrame)Frame
    integer x
    UNKNOWN_STRING  ,optional::color
    CALL PSGSAVE(Frame)
    if(present(color)) call pssetrgbcolor(frame,trim(color))
    CALL PSCURRENTPOINT(FRAME)
    CALL PsNewPath(Frame)
    CALL PsPush(Frame,"m")
    CALL PSRMOVETO(FRAME,0,FRAME%FONT%S/2-2)
    CALL PsRLineTo(Frame,X,0)
    CALL PSSTROKE(Frame)
    CALL PSGRESTORE(Frame)
    CALL PSRMOVETO(FRAME,X+6,0)
  END subroutine PSPlotLine_I
  
  Subroutine PSPlotVerticalLine(Frame,x)
    type(PSframe)Frame
    real(dl) x
    call PSPlotLine_D(Frame,x,Frame%Coor%OY+Frame%Coor%Yscale,x,Frame%Coor%YUR-Frame%Coor%Yscale)
  end Subroutine PSPlotVerticalLine


  Subroutine PSPlotHorizontalLine(Frame,y)
    type(PSframe)Frame
    real(dl) y
    call PSPlotLine_D(Frame,Frame%Coor%OX+Frame%Coor%Xscale,y,Frame%Coor%XUR-Frame%Coor%Xscale,y)
  end Subroutine PSPlotHorizontalLine


  subroutine PSPlotBox(Frame,WIDTH,HEIGHT)
    TYPE(PSFRAME)FRAME
    integer WIDTH
    integer,optional::HEIGHT
    integer H
    CALL PSRMOVETO(FRAME,WIDTH,0)
    CALL PSGSAVE(FRAME)
    CALL PSCURRENTPOINT(FRAME)
    CALL PsNewPath(FRAME)
    CALL PsPush(FRAME,"m")
    IF(PRESENT(HEIGHT))then
       H=HEIGHT
    ELSE
       H=CEILING(FRAME%FONT%S*0.85)
    endif
    CALL PsRLineTo(FRAME,0,H)
    CALL PsRLineTo(FRAME,-WIDTH,0)
    CALL PsRLineTo(FRAME,0,-H)
    CALL PSCLOSEPATH(FRAME)
    CALL PSFILL(FRAME)
    CALL PSGRESTORE(FRAME)
    CALL PSRMOVETO(FRAME,4,0)
  END subroutine PSPlotBox

  subroutine PSPlotBox_Color(Frame,Color,WIDTH,HEIGHT)
    TYPE(PSFRAME)FRAME
    integer WIDTH
    integer,optional::HEIGHT
    integer H
    UNKNOWN_STRING  COLOR
    CALL PSRMOVETO(FRAME,WIDTH,0)
    CALL PSGSAVE(FRAME)
    CALL PsSetRGBColor(FRAME,COLOR)
    CALL PSCURRENTPOINT(FRAME)
    CALL PsNewPath(FRAME)
    CALL PsPush(FRAME,"m")
    IF(PRESENT(HEIGHT))then
       H=HEIGHT
    ELSE
       H=CEILING(FRAME%FONT%S*0.85)
    endif
    CALL PsRLineTo(FRAME,0,H)
    CALL PsRLineTo(FRAME,-WIDTH,0)
    CALL PsRLineTo(FRAME,0,-H)
    CALL PSCLOSEPATH(FRAME)
    CALL PSFILL(FRAME)
    CALL PSGRESTORE(FRAME)
    CALL PSRMOVETO(FRAME,4,0)
  END subroutine PSPlotBox_Color

  Subroutine PSPlotEcllipse(Frame,cx,cy,covmat, CL, stroke_color,fill_color) 
    !!plot equal likelihood for a given confidence level;
    !!the likelihood is Gaussian with covariance matrix = covmat
    !!(cx, cy) is the best-fit point
    UNKNOWN_STRING   fill_color, stroke_color
    integer(IB),parameter::N=3600
    type(PSFrame)Frame
    real(dl) covmat(2,2), invcov(2,2), vec(2)
    real(dl) r, theta, chisq, x(N), y(N), cx, cy, CL
    integer(IB) i    
    invcov =  covmat
    if(covmat(1,2)**2 - covmat(1,1)*covmat(2,2) .ge. 0.) then
       print*,"PSPlotEcllipse: covmat not positive definite"
       return
    endif
    call matsym_inverse(invcov)
    if(CL.ge. 1._dl .or. CL.le.0._dl) stop "PSPlotEcllipse: 0< CL < 1 not satisfied"
    chisq = (-2._dl * log(1._dl - CL))
    do i=1, N
       theta = i*const_2pi/N
       vec=(/ cos(theta), sin(theta) /)
       r = sqrt(chisq/(dot_product(vec,matmul(invcov, vec))))
       x(i) = r*vec(1) + cx
       y(i) = r*vec(2) + cy
    enddo
    call psplotcurve(frame, x, y, stroke_color, fill_color)
  End Subroutine PSPlotEcllipse


  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine PSPrt(Frame,CHR,FontSize)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  CHR
    integer,optional::FontSize
    integer CurrentSize
    CurrentSize=Frame%Font%S
    IF(PRESENT(FontSize))then
       CALL PSsetfont(Frame,Frame%Font%T,FontSize)
    endif
    IF(CHR.ne."")then
       write(frame%fileunit,'(a)') "("//CHR//") show"
    endif
    CALL PSsetfont(Frame,Frame%Font%T,CurrentSize)
  END subroutine PSPrt

  subroutine PSPrtInt(Frame,N)
    TYPE(PSFrame)Frame
    integer N
    write(frame%fileunit,'(a)')"("//trim(INT2STR(N))//") show"
  END subroutine PSPrtInt

  subroutine PSPrtS(Frame,CHR,FontSize)
    TYPE(PSFrame) Frame
    UNKNOWN_STRING  CHR
    integer,optional::FontSize
    IF(PRESENT(FontSize))then
       CALL PSSymbolSwitch(Frame,FontSize)
    ELSE
       CALL PSSymbolSwitch(Frame)
    endif
    IF(CHR.ne."")then
       write(frame%fileunit,'(a)') "("//CHR//") show"
    endif
    CALL PSSymbolSwitch(Frame)
  END subroutine PSPrtS

  subroutine PSPrtSI(Frame,CHR,FontSize)
    TYPE(PSFrame) Frame
    UNKNOWN_STRING  CHR
    integer,optional::FontSize
    IF(PRESENT(FontSize))then
       CALL PSSymbolItalicSwitch(Frame,FontSize)
    ELSE
       CALL PSSymbolItalicSwitch(Frame)
    endif
    IF(CHR.ne."")then
       write(frame%fileunit,'(a)') "("//CHR//") show"
    endif
    CALL PSSymbolItalicSwitch(Frame)
  END subroutine PSPrtSI

  subroutine PSPrtN(Frame,CHR,FontSize)
    TYPE(PSFrame) Frame
    UNKNOWN_STRING  CHR
    integer,optional::FontSize
    IF(PRESENT(FontSize))then
       CALL PSNormalSwitch(Frame,FontSize)
    ELSE
       CALL PSNormalSwitch(Frame)
    endif
    IF(CHR.ne."")then
       write(frame%fileunit,'(a)') "("//CHR//") show"
    endif
    CALL PSNormalSwitch(Frame)
  END subroutine PSPrtN

  subroutine PSPrtCal(Frame,Chr,FontSize)
    TYPE(PSFrame) Frame
    UNKNOWN_STRING  CHR
    integer,optional::FontSize
    If(trim(Frame%dict).eq."")then
       print*,"warning: the cal font is not loaded"
       return
    endif
    IF(PRESENT(FontSize))then
       CALL PSCalSwitch(Frame,FontSize)
    ELSE
       CALL PSCalSwitch(Frame)
    endif
    IF(CHR.ne."")then
       write(frame%fileunit,'(a)') "("//CHR//") show"
    endif
    CALL PSCalSwitch(Frame)
  end subroutine PSPrtCal

  subroutine PSPrtI(Frame,CHR,FontSize)
    TYPE(PSFrame) Frame
    UNKNOWN_STRING  CHR
    integer,optional::FontSize
    IF(PRESENT(FontSize))then
       CALL PSItalicSwitch(Frame,FontSize)
    ELSE
       CALL PSItalicSwitch(Frame)
    endif
    IF(CHR.ne."")then
       write(frame%fileunit,'(a)') "("//CHR//") show"
    endif
    CALL PSItalicSwitch(Frame)
  END subroutine PSPrtI

  subroutine PSPrtB(Frame,CHR,FontSize)
    TYPE(PSFrame) Frame
    UNKNOWN_STRING  CHR
    integer,optional::FontSize
    IF(PRESENT(FontSize))then
       CALL PSBoldSwitch(Frame,FontSize)
    ELSE
       CALL PSBoldSwitch(Frame)
    endif
    IF(CHR.ne."")then
       write(frame%fileunit,'(a)') "("//CHR//") show"
    endif
    CALL PSBoldSwitch(Frame)
  END subroutine PSPrtB

  subroutine PSPrtBI(Frame,CHR,FontSize)
    TYPE(PSFrame) Frame
    UNKNOWN_STRING  CHR
    integer,optional::FontSize
    IF(PRESENT(FontSize))then
       CALL PSBoldItalicSwitch(Frame,FontSize)
    ELSE
       CALL PSBoldItalicSwitch(Frame)
    endif
    IF(CHR.ne."")then
       write(frame%fileunit,'(a)') "("//CHR//") show"
    endif
    CALL PSBoldItalicSwitch(Frame)
  END subroutine PSPrtBI

  subroutine PSShowSub(Frame,CHR)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  CHR
    CALL PSRMOVETO(Frame,0,-4)
    write(frame%fileunit,'(a)') "("//CHR//") show"
    CALL PSRMOVETO(Frame,0,4)
  END subroutine PSShowSub

  subroutine PSSub(Frame,CHR)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  CHR
    integer CurrentSize
    CurrentSize=Frame%Font%S
    CALL PSsetfont(Frame,Frame%Font%T,MAX(Frame%Font%S/2+2,12))
    CALL PSSHOWSUB(Frame,CHR)
    CALL PSsetfont(Frame,Frame%Font%T,CurrentSize)
  END subroutine PSSub

  subroutine PSSubS(Frame,CHR)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  CHR
    CALL PSSymbolSwitch(Frame,MAX(Frame%Font%S/2+2,12))
    CALL PSSHOWSUB(Frame,CHR)
    CALL PSSymbolSwitch(Frame)
  END subroutine PSSubS

  subroutine PSSubI(Frame,CHR)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  CHR
    CALL PSItalicSwitch(Frame,MAX(Frame%Font%S/2+2,12))
    CALL PSSHOWSUB(Frame,CHR)
    CALL PSItalicSwitch(Frame)
  END subroutine PSSubI

  subroutine PSSubB(Frame,CHR)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  CHR
    CALL PSBoldSwitch(Frame,MAX(Frame%Font%S/2+2,12))
    CALL PSSHOWSUB(Frame,CHR)
    CALL PSBoldSwitch(Frame)
  END subroutine PSSubB

  subroutine PSSubBI(Frame,CHR)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  CHR
    CALL PSBoldItalicSwitch(Frame,MAX(Frame%Font%S/2+2,12))
    CALL PSSHOWSUB(Frame,CHR)
    CALL PSBoldItalicSwitch(Frame)
  END subroutine PSSubBI

  subroutine PSSup(Frame,CHR)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  CHR
    integer CurrentSize
    CurrentSize=Frame%Font%S
    CALL PSsetfont(Frame,Frame%Font%T,MAX(Frame%Font%S/2+2,12))
    CALL PSRMOVETO(Frame,0,CurrentSize-Frame%Font%S+2)
    write(frame%fileunit,'(a)') "("//CHR//") show"
    CALL PSRMOVETO(Frame,0,Frame%Font%S-CurrentSize-2)
    CALL PSsetfont(Frame,Frame%Font%T,CurrentSize)
  END subroutine PSSup

  subroutine PSSupS(Frame,CHR)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  CHR
    CALL PSSymbolSwitch(Frame)
    CALL PSSup(Frame,CHR)
    CALL PSSymbolSwitch(Frame)
  END subroutine PSSupS

  subroutine PSSupI(Frame,CHR)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  CHR
    CALL PSItalicSwitch(Frame)
    CALL PSSup(Frame,CHR)
    CALL PSItalicSwitch(Frame)
  END subroutine PSSupI

  subroutine PSSupB(Frame,CHR)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  CHR
    CALL PSBoldSwitch(Frame)
    CALL PSSup(Frame,CHR)
    CALL PSBoldSwitch(Frame)
  END subroutine PSSupB

  subroutine PSSupBI(Frame,CHR)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  CHR
    CALL PSBoldItalicSwitch(Frame)
    CALL PSSup(Frame,CHR)
    CALL PSBoldItalicSwitch(Frame)
  END subroutine PSSupBI

  subroutine PSSupSub(Frame,SUPCHR,SUBCHR)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING SUPCHR,SUBCHR
    integer CurrentSize
    CurrentSize=Frame%Font%S
    CALL PSsetfont(Frame,Frame%Font%T,MAX(Frame%Font%S/2+2,12))
    CALL PSRMOVETO(Frame,0,CurrentSize-Frame%Font%S+3)
    write(frame%fileunit,'(a)') "("//SUPCHR//") prtcl"
    CALL PSRMOVETO(Frame,0,Frame%Font%S-CurrentSize-6)
    write(frame%fileunit,'(a)') "("//SUBCHR//") prtcl"
    CALL PsPush(Frame,"larger 3 rmoveto")
    CALL PSsetfont(Frame,Frame%Font%T,CurrentSize)
  END subroutine PSSupSub

  subroutine PSSupSSubS(Frame,SUPCHR,SUBCHR)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  SUPCHR,SUBCHR
    CALL PSSymbolSwitch(Frame)
    CALL PSSupSub(Frame,SUPCHR,SUBCHR)
    CALL PSSymbolSwitch(Frame)
  END subroutine PSSupSSubS

  subroutine PSSupSubS(Frame,SUPCHR,SUBCHR)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  SUPCHR,SUBCHR
    integer CurrentSize
    CurrentSize=Frame%Font%S
    CALL PSsetfont(Frame,Frame%Font%T,MAX(Frame%Font%S/2+2,12))
    CALL PSRMOVETO(Frame,0,CurrentSize-Frame%Font%S+3)
    write(frame%fileunit,'(a)') "("//SUPCHR//") prtcl"
    CALL PSRMOVETO(Frame,0,Frame%Font%S-CurrentSize-6)
    CALL PSSymbolSwitch(Frame)
    write(frame%fileunit,'(a)') "("//SUBCHR//") prtcl"
    CALL PSSymbolSwitch(Frame)
    CALL PsPush(Frame,"larger 3 rmoveto")
    CALL PSsetfont(Frame,Frame%Font%T,CurrentSize)
  END subroutine PSSupSubS

  subroutine PSSupSSub(Frame,SUPCHR,SUBCHR)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  SUPCHR,SUBCHR
    integer CurrentSize
    CurrentSize=Frame%Font%S
    CALL PSsetfont(Frame,Frame%Font%T,MAX(Frame%Font%S/2+2,12))
    CALL PSSymbolSwitch(Frame)
    CALL PSRMOVETO(Frame,0,CurrentSize-Frame%Font%S+3)
    write(frame%fileunit,'(a)') "("//SUPCHR//") prtcl"
    CALL PSSymbolSwitch(Frame)
    CALL PSRMOVETO(Frame,0,Frame%Font%S-CurrentSize-6)
    write(frame%fileunit,'(a)') "("//SUBCHR//") prtcl"
    CALL PsPush(Frame,"larger 3 rmoveto")
    CALL PSsetfont(Frame,Frame%Font%T,CurrentSize)
  END subroutine PSSupSSub

  subroutine PSSupSubI(Frame,SUPCHR,SUBCHR)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  SUPCHR,SUBCHR
    integer CurrentSize
    CurrentSize=Frame%Font%S
    CALL PSsetfont(Frame,Frame%Font%T,MAX(Frame%Font%S/2+2,12))
    CALL PSRMOVETO(Frame,0,CurrentSize-Frame%Font%S+3)
    write(frame%fileunit,'(a)') "("//SUPCHR//") prtcl"
    CALL PSRMOVETO(Frame,0,Frame%Font%S-CurrentSize-6)
    CALL PSItalicSwitch(Frame)
    write(frame%fileunit,'(a)') "("//SUBCHR//") prtcl"
    CALL PSItalicSwitch(Frame)
    CALL PsPush(Frame,"larger 3 rmoveto")
    CALL PSsetfont(Frame,Frame%Font%T,CurrentSize)
  END subroutine PSSupSubI

  subroutine PSSupSubB(Frame,SUPCHR,SUBCHR)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  SUPCHR,SUBCHR
    integer CurrentSize
    CurrentSize=Frame%Font%S
    CALL PSsetfont(Frame,Frame%Font%T,MAX(Frame%Font%S/2+2,12))
    CALL PSRMOVETO(Frame,0,CurrentSize-Frame%Font%S+3)
    write(frame%fileunit,'(a)') "("//SUPCHR//") prtcl"
    CALL PSRMOVETO(Frame,0,Frame%Font%S-CurrentSize-6)
    CALL PSBoldSwitch(Frame)
    write(frame%fileunit,'(a)') "("//SUBCHR//") prtcl"
    CALL PSBoldSwitch(Frame)
    CALL PsPush(Frame,"larger 3 rmoveto")
    CALL PSsetfont(Frame,Frame%Font%T,CurrentSize)
  END subroutine PSSupSubB

  subroutine PSSupSubBI(Frame,SUPCHR,SUBCHR)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  SUPCHR,SUBCHR
    integer CurrentSize
    CurrentSize=Frame%Font%S
    CALL PSsetfont(Frame,Frame%Font%T,MAX(Frame%Font%S/2+2,12))
    CALL PSRMOVETO(Frame,0,CurrentSize-Frame%Font%S+3)
    write(frame%fileunit,'(a)') "("//SUPCHR//") prtcl"
    CALL PSRMOVETO(Frame,0,Frame%Font%S-CurrentSize-6)
    CALL PSBoldItalicSwitch(Frame)
    write(frame%fileunit,'(a)') "("//SUBCHR//") prtcl"
    CALL PSBoldItalicSwitch(Frame)
    CALL PsPush(Frame,"larger 3 rmoveto")
    CALL PSsetfont(Frame,Frame%Font%T,CurrentSize)
  END subroutine PSSupSubBI

  subroutine PSFrac(Frame,Denominator,Numerator)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  DENOMINATOR,NUMERATOR
    integer CurrentSize
    CurrentSize=Frame%Font%S
    CALL PSRMOVETO(Frame,1,CurrentSize/3)
  !  CALL PSsetfont(Frame,Frame%Font%T,MAX(CurrentSize/2+6,10))
    write(frame%fileunit,'(a)')"("//DENOMINATOR//") ("//NUMERATOR//") fracm"
    CALL PSRMOVETO(Frame,1,-CurrentSize/3)
    CALL PSsetfont(Frame,Frame%Font%T,CurrentSize)
  END subroutine PSFrac

  subroutine PSFracSOS(Frame,Denominator,Numerator)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  DENOMINATOR,NUMERATOR
    CALL PSSymbolSwitch(Frame)
    CALL PSFrac(Frame,Denominator,Numerator)
    CALL PSSymbolSwitch(Frame)
  END subroutine PSFracSOS

  subroutine PSFracOS(Frame,Denominator,Numerator)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  DENOMINATOR,NUMERATOR
    integer CurrentSize
    SHORT_STRING CurrentType
    CurrentType=trim(Frame%Font%T)
    CurrentSize=Frame%Font%S
    CALL PSRMOVETO(Frame,1,CurrentSize/3)
  !  CALL PSsetfont(Frame,Frame%Font%T,MAX(CurrentSize/2+6,10))
    write(frame%fileunit,'(a)')trim(INT2STR(Frame%Font%S))//  &
         "("//DENOMINATOR//") ("//NUMERATOR//") fracms"
    CALL PSRMOVETO(Frame,1,-CurrentSize/3)
    CALL PSsetfont(Frame,CurrentType,CurrentSize)
  END subroutine PSFracOS

  subroutine PSFracON(Frame,Denominator,Numerator)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  DENOMINATOR,NUMERATOR
    integer CurrentSize
    SHORT_STRING  CurrentType
    CurrentType=trim(Frame%Font%T)
    CurrentSize=Frame%Font%S
    CALL PSRMOVETO(Frame,1,CurrentSize/3)
  !  CALL PSsetfont(Frame,Frame%Font%T,MAX(CurrentSize/2+6,10))
    write(frame%fileunit,'(a)')trim(INT2STR(Frame%Font%S))//  &
         "("//DENOMINATOR//") ("//NUMERATOR//") fracmn"
    CALL PSRMOVETO(Frame,1,-CurrentSize/3)
    CALL PSsetfont(Frame,CurrentType,CurrentSize)
  END subroutine PSFracON

  subroutine PSFracOI(Frame,Denominator,Numerator)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  DENOMINATOR,NUMERATOR
    integer CurrentSize
    SHORT_STRING  CurrentType
    CurrentType=trim(Frame%Font%T)
    CurrentSize=Frame%Font%S
    CALL PSRMOVETO(Frame,1,CurrentSize/3)
  !  CALL PSsetfont(Frame,Frame%Font%T,MAX(CurrentSize/2+6,10))
    write(frame%fileunit,'(a)')trim(INT2STR(Frame%Font%S))//  &
         "("//DENOMINATOR//") ("//NUMERATOR//") fracmi"
    CALL PSRMOVETO(Frame,1,-CurrentSize/3)
    CALL PSsetfont(Frame,CurrentType,CurrentSize)
  END subroutine PSFracOI

  subroutine PSFracOB(Frame,Denominator,Numerator)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  DENOMINATOR,NUMERATOR
    integer CurrentSize
    SHORT_STRING  CurrentType
    CurrentType=trim(Frame%Font%T)
    CurrentSize=Frame%Font%S
    CALL PSRMOVETO(Frame,1,CurrentSize/3)
  !  CALL PSsetfont(Frame,Frame%Font%T,MAX(CurrentSize/2+6,10))
    write(frame%fileunit,'(a)')trim(INT2STR(Frame%Font%S))//  &
         "("//DENOMINATOR//") ("//NUMERATOR//") fracmb"
    CALL PSRMOVETO(Frame,1,-CurrentSize/3)
    CALL PSsetfont(Frame,CurrentType,CurrentSize)
  END subroutine PSFracOB

  subroutine PSFracOBI(Frame,Denominator,Numerator)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  DENOMINATOR,NUMERATOR
    integer CurrentSize
    SHORT_STRING  CurrentType
    CurrentType=trim(Frame%Font%T)
    CurrentSize=Frame%Font%S
    CALL PSRMOVETO(Frame,1,CurrentSize/3)
  !  CALL PSsetfont(Frame,Frame%Font%T,MAX(CurrentSize/2+6,10))
    write(frame%fileunit,'(a)')trim(INT2STR(Frame%Font%S))//  &
         "("//DENOMINATOR//") ("//NUMERATOR//") fracmbi"
    CALL PSRMOVETO(Frame,1,-CurrentSize/3)
    CALL PSsetfont(Frame,CurrentType,CurrentSize)
  END subroutine PSFracOBI

  subroutine PSBar(Frame,CHR)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  CHR
    CALL PSRMOVETO(Frame,2,0)
    write(frame%fileunit,'(a)')"("//CHR//") psbar"
    CALL PSRMOVETO(Frame,2,0)
  END subroutine PSBar

  subroutine PSSqrt(Frame,CHR,N)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  CHR
    integer,optional::N
    IF(PRESENT(N))then
       CALL PSRMOVETO(Frame,0,Frame%Font%S/2+2)
       CALL PSSYMBOLSWITCH(Frame,MAX(Frame%Font%S/2-2,10))
       write(frame%fileunit,'(a)')"("//trim(INT2STR(N))//") show"
       CALL PSSYMBOLSWITCH(Frame)
       CALL PSRMOVETO(Frame,-Frame%Font%S/4,-Frame%Font%S/2-2)
    endif
    CALL PSPRTS(Frame,"\326")
    write(frame%fileunit,'(a)')"("//CHR//") psbar"
  END subroutine PSSqrt

  RECURSIVE subroutine PSPrtTex(Frame,CHR,FontSize)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  CHR
    integer,optional::FontSize
    integer StrLen,CurrentFS,I,J,FS,RM,RestoreSize,FRU,FRD,lu, ld
    SHORT_STRING  CurrentFT
    logical SQRTRESTORE,UDRESTORE
    lu = 0
    ld = 0
    StrLen=Len_Trim(CHR)
    IF(StrLen.EQ.0)RETURN
    SQRTRESTORE=.FALSE.
    CurrentFS=Frame%Font%S
    CurrentFT=Frame%Font%T
    IF(PRESENT(FontSize))then
       CALL PSsetfont(Frame,trim(Default_Font),FontSize)
       FS=FontSize
    ELSE
       CALL PSsetfont(Frame,trim(Default_Font))
       FS=CURRENTFS
    endif
    I=SCAN(CHR(1:StrLen),"\^_{")
    IF(I.EQ.0)then
       CALL PSPrtSingleTex(Frame,CHR(1:StrLen))
       RETURN
    ELSEIF(I.GT.1)then
       CALL PSPrtSingleTex(Frame,CHR(1:I-1))
    endif
    RM=0
    UDRESTORE=.FALSE.
    do
       if(I.GT.STRLEN)exit
       select case(CHR(I:I))
       case("^")
          CALL PSRMOVETO(FRAME,0,FS/2+RM)
          CALL PSSETFONT(FRAME,trim(Default_Font),MAX(FS*2/3,12))
          IF(CHR(I+1:I+1).ne."{")then
             IF(CHR(I+1:I+1).ne.Const_backslash)then
                j = 1
             ELSE
                J=VERIFY(CHR(I+3:STRLEN),"abcdefghijklmnopqrstuvwxyz")
                IF(J.EQ.0)J=STRLEN-I-1
                J=J+1
             endif
             CALL PSPRTSingleTex(FRAME,CHR(I+1:I+J))
             lu=PSSlSingleTex(FRAME,CHR(I+1:I+J))
             I=I+J+1
          ELSE
             J=NEXTBRACKET(CHR(I+2:STRLEN))
             IF(J.EQ.0)STOP "PSPrtTex: Syntax error."
             CALL PSPRTTEX(FRAME,CHR(I+2:I+J))
             lu=PSSlTex(FRAME,CHR(I+2:I+J))
             I=I+J+2
          endif
          CALL PSRMOVETO(FRAME,0,-FS/2-RM)
          IF(UDRESTORE)then
             CALL PSRMOVETO(FRAME,MAX(lu,ld)-lu,0)
             UDRESTORE=.FALSE.
          ELSEIF(I.LE.STRLEN .AND. CHR(I:I).EQ."_")then
             CALL PSRMOVETO(FRAME,-LU,0)
             UDRESTORE=.TRUE.
          endif
          CALL PSSETFONT(FRAME,trim(Default_Font),FS)
          CALL PSRMOVETO(FRAME,2,0)
          RM=0
       case("_")
          call PSRMOVETO(FRAME,0,-6)
          call PSSETFONT(FRAME,trim(Default_Font),MAX(FS*2/3,12))
          IF(CHR(I+1:I+1).ne."{")then
             IF(chr(i+1:i+1) .ne. Const_backslash)then
                j = 1
             else
                J=verify(CHR(I+3:STRLEN),"abcdefghijklmnopqrstuvwxyz")
                if(J.eq.0)J=STRLEN-I-1
                J=J+1
             endif
             CALL PSPRTSingleTex(FRAME,CHR(I+1:I+J))
             LD=PSSlSingleTex(FRAME,CHR(I+1:I+J))
             I=I+J+1
          else
             J=NEXTBRACKET(CHR(I+2:STRLEN))
             IF(J.EQ.0)STOP "PSPrtTex: Syntax error."
             CALL PSPRTTEX(FRAME,CHR(I+2:I+J))
             LD=PSSlTEX(FRAME,CHR(I+2:I+J))
             I=I+J+2
          endif
          call PSRMOVETO(FRAME,0,6)
          IF(UDRESTORE)then
             CALL PSRMOVETO(FRAME,MAX(LU,LD)-LD,0)
             UDRESTORE=.FALSE.
          ELSEIF(I.LE.STRLEN .AND. CHR(I:I).EQ."^")then
             CALL PSRMOVETO(FRAME,-LD,0)
             UDRESTORE=.TRUE.
          endif
          CALL PSSETFONT(FRAME,trim(Default_Font),FS)
          CALL PSRMOVETO(FRAME,2,0)
       case("{")
          J=NEXTBRACKET(CHR(I+1:STRLEN))
          IF(J.EQ.0)STOP "PSPrtTex: Syntax error."
          CALL PSPRTTEX(FRAME,CHR(I+1:I+J-1))
          I=I+J+1
          IF(SQRTRESTORE)then
             CALL PSGSAVE(FRAME)
             CALL PSCURRENTPOINT(Frame)
             CALL PSNewPath(FRAME)
             CALL PsPush(FRAME,"m")
             CALL PSRMOVETO(FRAME,0,RestoreSize*11/12)
             CALL PsPush(FRAME,"l")
             CALL PSSTROKE(FRAME)
             CALL PSGRESTORE(FRAME)
             CALL PSSETFONT(FRAME,FRAME%FONT%T,RESTORESIZE)
             SQRTRESTORE=.FALSE.
          endif
          RM=0
       case(Const_backslash)
          RM=0
          IF(CHR(I+1:I+1).eq." ")then
             call PSprtN(Frame,"\040")
             I=I+2
             cycle
          endif
          J=VERIFY(CHR(I+2:StrLen),"abcdefghijklmnopqrstuvwxyz")
          if(J.EQ.0)then
             call PSPrtSingleTex(Frame,CHR(I:StrLen))
             return
          endif
          Select case(trim(chr(i:i+j)))
          case("\sqrt")
             IF(CHR(I+5:I+5).EQ."{")then
                CALL PSPrtS(Frame,"\326")
                I=I+5
             ELSEIF(CHR(I+5:I+5).EQ."[")then
                J=SCAN(CHR(I+6:STRLEN),"]")
                CALL PSRMOVETO(FRAME,FRAME%FONT%S/4,FRAME%FONT%S*2/3)
                RestoreSize=Frame%Font%S
                CALL PSSETFONT(FRAME,FRAME%Font%T,MAX(Frame%Font%S/3+2,10))
                CALL PSPRTTEX(FRAME,CHR(I+6:I+J+4))
                CALL PSSETFONT(FRAME,FRAME%FONT%T,RESTORESIZE)
                CALL PSRMOVETO(FRAME,-FRAME%FONT%S/4,-FRAME%FONT%S*2/3)
                CALL PSPrtS(Frame,"\326")
                I=I+J+6
             else
                write(*,*) "PSPrtTex: Syntax error."
                stop
             endif
             RestoreSize=Frame%Font%S
             CALL PSRMOVETO(FRAME,0,Frame%Font%S*11/12)
             CALL PSCURRENTPOINT(FRAME)
             CALL PSRMOVETO(FRAME,0,-Frame%Font%S*11/12)
             CALL PSSETFONT(FRAME,FRAME%Font%T,Frame%Font%S*11/12)
             SQRTRESTORE=.TRUE.
          case("\frac")
             RestoreSize=Frame%Font%S
             CALL PSSETFONT(FRAME,FRAME%FONT%T,MAX(FRAME%FONT%S*3/4,10))
             IF(CHR(I+5:I+5).ne."{")STOP "PSPrtTex: Syntax Error."
             J=NEXTBRACKET(CHR(I+6:STRLEN))
             IF(J.EQ.0)STOP "PSPrtTex: Syntax Error."
             FRU=I+J+5
             LU=PSSlTex(Frame,CHR(I+6:FRU-1))+3
             IF(CHR(FRU+1:FRU+1).ne."{")STOP "PSPrtTex: Syntax Error."
             J=NEXTBRACKET(CHR(FRU+2:STRLEN))
             IF(J.EQ.0)STOP "PSPrtTex: Syntax Error."
             FRD=FRU+J+1
             LD=PSSlTex(Frame,CHR(FRU+2:FRD-1))+3
             CALL PSRMOVETO(FRAME,MAX(LU,LD),0)
             CALL PSRMOVETO(FRAME,0,RestoreSize/3)
             CALL PSGSAVE(FRAME)
             CALL PSCANCELDASH(FRAME)
             CALL PSCURRENTPOINT(FRAME)
             CALL PsNewPath(FRAME)
             CALL PsPush(FRAME,"m")
             CALL PsRLineTo(FRAME,-MAX(LU,LD),0)
             CALL PSSTROKE(FRAME)
             CALL PSGRESTORE(FRAME)
             CALL PSGSAVE(FRAME)
             CALL PSRMOVETO(FRAME,-(MAX(LU,LD)+LU)/2,6)
             CALL PSPRTTEX(FRAME,CHR(I+6:FRU-1))
             CALL PSGRESTORE(FRAME)
             CALL PSGSAVE(FRAME)
             if(scan(CHR(FRU+2:FRD-1),"ABCDEFGHIJKLMNOPQRSTUVWXYZ"//Const_backslash).EQ.0)then
                call psrmoveto(FRAME,-(MAX(LU,LD)+LD)/2,-FRAME%FONT%S)
             else
                call psrmoveto(FRAME,-(MAX(LU,LD)+LD)/2,-FRAME%FONT%S-2)
             endif
             CALL PSPRTTEX(FRAME,CHR(FRU+2:FRD-1))
             CALL PSGRESTORE(FRAME)
             CALL PSRMOVETO(FRAME,0,-RestoreSize/3)
             CALL PSSETFONT(FRAME,Frame%Font%T,RestoreSize)
             I=FRD+1
          case("\cal","\Cal")
             IF(CHR(I+4:I+4).ne."{")STOP "PSPrtTex:Syntax Error: only support \cal{}."
             J=NEXTBRACKET(CHR(I+5:STRLEN))
             if(J.eq.0)stop "PSPrtTex Syntax error: only support \cal{}."
             CALL PSPRTCAL(Frame,CHR(I+5:I+J+3))
             I=I+J+5
          case("\tilde")
             if(Chr(I+6:I+6) .ne. "{") stop "psprttex: syntex error, only support \tildeP{}"
             J=NEXTBRACKET(CHR(I+7:STRLEN))
             if(J.eq.0)stop "PSPrtTex Syntax error: only support \tilde{}."
             call pscurrentpoint(frame)
             call psRmoveto(Frame, 0, Frame%Font%S*2/3)
             call psprtN(Frame, "~")
             call pspush(frame,"m")
             CALL PSPrtTex(Frame,CHR(I+7:I+J+5))
             I=I+J+7
          case("\textit")
             IF(CHR(I+7:I+7).ne."{")STOP "PSPrtTex: Syntax Error: only support \textit{}."
             J=NEXTBRACKET(CHR(I+8:STRLEN))
             IF(J.EQ.0)STOP "PSPrtTex: Syntax Error: only support \textit{}."
             CALL PSPrtI(FRAME,CHR(I+8:I+J+6))
             I=I+J+8
          case("\text")
             IF(CHR(I+5:I+5).ne."{")STOP "PSPrtTex: Syntax Error: only support \text{}."
             J=NEXTBRACKET(CHR(I+6:STRLEN))
             IF(J.EQ.0)STOP "PSPrtTex: Syntax Error: only support \text{}."
             CALL PSPRTN(FRAME,CHR(I+6:I+J+4))
             I=I+J+6
          case("\textbf")
             IF(CHR(I+7:I+7).ne."{")STOP "PSPrtTex: Syntax Error: only support \textbf{}."
             J=NEXTBRACKET(CHR(I+8:STRLEN))
             IF(J.EQ.0)STOP "PSPrtTex: Syntax Error: only support \textbf{}."
             CALL PSPRTB(FRAME,CHR(I+8:I+J+6))
             I=I+J+8
          case default
             CALL PSPrtSingleTex(Frame,CHR(I:I+J))
             IF(CHR(I:I+J).EQ."\int")RM=Frame%Font%S/6
             I=I+J+1
          end select
       case default
          J=scan(CHR(I+1:STRLEN),"^_{"//Const_backslash)
          IF(J.EQ.0)then
             call PSPrtSingleTex(Frame,CHR(I:StrLen))
             return
          endif
          call PSPRTSINGLETEX(FRAME,CHR(I:I+J-1))
          I=I+J
       end select
    enddo
    call PSsetfont(Frame,CurrentFT,CurrentFS)
  END subroutine PSPrtTex

  subroutine PSPrtTexAtCenter(Frame,CHR,FontSize)
    TYPE(PSFRAME)FRAME
    UNKNOWN_STRING  CHR
    integer,optional::FontSize
    integer CurrentFS
    SHORT_STRING  CurrentFT
    CurrentFS=Frame%Font%S
    CurrentFT=Frame%Font%T
    IF(PRESENT(FontSize))then
       CALL PSsetfont(Frame,trim(Default_Font),FontSize)
    ELSE
       CALL PSsetfont(Frame,trim(Default_Font))
    endif
    CALL PSRMOVETO(FRAME,-PSSlTex(Frame,CHR)/2,0)
    CALL PSPrtTex(Frame,CHR)
    CALL PSsetfont(Frame,CurrentFT,CurrentFS)
  END subroutine PSPrtTexAtCenter

  subroutine PSPrtTexAtLeft(Frame,CHR,FontSize)
    TYPE(PSFRAME)FRAME
    UNKNOWN_STRING  CHR
    integer,optional::FontSize
    integer CurrentFS
    SHORT_STRING  CurrentFT
    CurrentFS=Frame%Font%S
    CurrentFT=Frame%Font%T
    IF(PRESENT(FontSize))then
       CALL PSsetfont(Frame,trim(Default_Font),FontSize)
    ELSE
       CALL PSsetfont(Frame,trim(Default_Font))
    endif
    CALL PSRMOVETO(FRAME,-PSSlTex(Frame,CHR),0)
    CALL PSPrtTex(Frame,CHR)
    CALL PSsetfont(Frame,CurrentFT,CurrentFS)
  END subroutine PSPrtTexAtLeft

  subroutine PSPrtSingle(Frame,C)
    TYPE(PSFrame)Frame
    character C
    integer ASCII
    Select case(C)
    Case("(",")")
       CALL PSPRTN(Frame,Const_backslash//C)
    case("=","-","+","[","]")
       CALL PSPRTN(Frame," "//C//" ")
    case("a")
       call PSPrtI(Frame,C)
    case default
       ASCII=ICHAR(C)
       IF((ASCII.GE.42 .AND. ASCII.LE.63).OR.(ASCII.GE.123 .AND. ASCII.LE.126))then
          CALL PSPRTS(Frame,C)
       ELSEIF(ASCII.GE.33 .AND. ASCII.LE.122)then
          CALL PSPRT(Frame,C)
       endif
    END select
  END subroutine PSPrtSingle

  subroutine PSPrtSingleTex(Frame,CHR)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  CHR
    integer StrLen, I
    StrLen=Len_Trim(CHR)
    IF(StrLen.EQ.0)RETURN
    IF(CHR(1:1).ne."\")then
       DO I=1,StrLen
          CALL PSPrtSingle(Frame,CHR(I:I))
       enddo
       return
    ELSEIF(StrLen.EQ.1)then
       CALL PSPrt(Frame,"\040")
       return
    Endif
    Select case(CHR(2:Strlen))
    case("\")
       CALL PSPrtN(Frame,"\134")
    case("^","_")
       CALL PSPrtN(Frame,CHR(2:2))
    case("omega")
       CALL PSPrtSI(Frame,"w")
    case("Omega")
       CALL PSPrtSI(Frame,"W")
    case("psi")
       CALL PSPrtSI(Frame,"y")
    case("Psi")
       CALL PSPrtSI(Frame,"Y")
    case("phi")
       CALL PSPrtSI(Frame,"f")
    case("Phi")
       CALL PSPrtSI(Frame,"F")
    case("theta")
       CALL PSPrtSI(Frame,"q")
    case("Theta")
       CALL PSPrtSI(Frame,"Q")
    case("eta")
       CALL PSPrtSI(Frame,"h")
    case("Upsilon")
       CALL PSPrtSI(Frame,"\241")
    case("varpi")
       CALL PSPrtSI(Frame,"v")
    case("varsigma")
       CALL PSPrtSI(Frame,"V")
    case("vartheta")
       CALL PSPrtSI(Frame,"J")
    case("varphi")
       CALL PSPrtSI(Frame,"j")
    case("exists")
       CALL PSPrts(Frame,"\044")
    case("forall")
       CALL PSPrts(Frame,"\042")
    case("cong")
       CALL PSPrts(Frame,"\100")
    case("le")
       CALL PSPrts(Frame,"\243")
    case("infty")
       CALL PSPrts(Frame,"\245")
    case("clubsuit")
       CALL PSPrts(Frame,"\247")
    case("diamondsuit")
       CALL PSPrts(Frame,"\250")
    case("heartsuit")
       CALL PSPrts(Frame,"\251")
    case("spadesuit")
       CALL PSPrts(Frame,"\252")
    case("leftrightarrow")
       CALL PSPrts(Frame,"\253")
    case("leftarrow")
       CALL PSPrts(Frame,"\254")
    case("uparrow")
       CALL PSPrts(Frame,"\255")
    case("rightarrow")
       CALL PSPrts(Frame,"\256")
    case("downarrow")
       CALL PSPrts(Frame,"\257")
    case("pm")
       CALL PSPrts(Frame,"\261")
    case("ge")
       CALL PSPrts(Frame,"\263")
    case("times")
       CALL PSPrts(Frame,"\264")
    case("propto")
       CALL PSPrts(Frame,"\265")
    case("partial")
       CALL PSPrts(Frame,"\266")
    case("bullet")
       CALL PSPrts(Frame,"\267")
    case("div")
       CALL PSPrts(Frame,"\270")
    case("neq")
       CALL PSPrts(Frame,"\271")
    case("equiv")
       CALL PSPrts(Frame,"\272")
    case("approx")
       CALL PSPrts(Frame,"\273")
    case("aleph")
       CALL PSPrts(Frame,"\300")
    case("Re")
       CALL PSPrts(Frame,"\302")
    case("otimes")
       CALL PSPrts(Frame,"\304")
    case("oplus")
       CALL PSPrts(Frame,"\305")
    case("emptyset")
       CALL PSPrts(Frame,"\306")
    case("cap")
       CALL PSPrts(Frame,"\307")
    case("cup")
       CALL PSPrts(Frame,"\310")
    case("supset")
       CALL PSPrts(Frame,"\311")
    case("supseteq")
       CALL PSPrts(Frame,"\312")
    case("subset")
       CALL PSPrts(Frame,"\314")
    case("subseteq")
       CALL PSPrts(Frame,"\315")
    case("in")
       CALL PSPrts(Frame,"\316")
    case("ni")
       CALL PSPrts(Frame,"\317")
    case("nabla")
       CALL PSPrts(Frame,"\321")
    case("prod")
       CALL PSPrts(Frame,"\325")
    case("cdot")
       CALL PSPrts(Frame,"\327")
    case("Leftrightarrow")
       CALL PSPrts(Frame,"\333")
    case("Leftarrow")
       CALL PSPrts(Frame,"\334")
    case("Uparrow")
       CALL PSPrts(Frame,"\335")
    case("Rightarrow")
       CALL PSPrts(Frame,"\336")
    case("Downarrow")
       CALL PSPrts(Frame,"\337")
    case("langle")
       CALL PSPrts(Frame,"\341")
    case("Sum")
       CALL PSPrts(Frame,"\345")
    case("rangle")
       CALL PSPrts(Frame,"\361")
    case("int")
       CALL PSPrts(Frame,"\362")
    case("hbar")
       CALL PSPRT(FRAME,"h")
       CALL PSGSAVE(FRAME)
       CALL PSCURRENTPOINT(FRAME)
       CALL PsNewPath(FRAME)
       CALL PsPush(FRAME,"m")
       CALL PSRMOVETO(FRAME,0,FRAME%FONT%S*3/5)
       CALL PsRLineTo(FRAME,-FRAME%FONT%S/2,-FRAME%FONT%S/12)
       CALL PSSTROKE(FRAME)
       CALL PSGRESTORE(FRAME)
    case("sin","cos","tan","cot","sec","csc","exp","ln","log","lg","sinh","cosh","tanh","coth","arcsin","arccos","arctan","lim","max","min" &
         ,"ker","dim","det","gcd","limsup","liminf","sup","inf","arg","hom")
       CALL PSPrtN(Frame,CHR(2:StrLen))
    case default
       CALL PSPrts(Frame,CHR(2:2))
    end select
  END subroutine PSPrtSingleTex


  RECURSIVE function PSSlTex(Frame,CHR) RESULT(SL)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  CHR
    integer I,J,SL,StrLen,FRD,FRU,LD,LU
    StrLen=LEN_trim(CHR)
    I=SCAN(CHR(1:StrLen),"\^_{")
    SL=0
    IF(I.EQ.0)then
       SL=PSSlSingleTex(Frame,CHR(1:StrLen))
       RETURN
    ELSEIF(I.GT.1)then
       SL=PSSlSingleTex(Frame,CHR(1:I-1))
    endif
    DO
       IF(I.GT.STRLEN)EXIT
       IF(CHR(I:I).EQ."^")then
          IF(CHR(I+1:I+1).ne."{")then
             IF(CHR(I+1:I+1).ne.Const_backslash)then
                J=1
             ELSE
                J=VERIFY(CHR(I+3:STRLEN),"abcdefghijklmnopqrstuvwxyz")
                IF(J.EQ.0)J=STRLEN-I-1
                J=J+1
             endif
             SL=SL+PSSlSingleTex(FRAME,CHR(I+1:I+J))/2
             I=I+J+1
          ELSE
             J=NEXTBRACKET(CHR(I+2:STRLEN))
             IF(J.EQ.0)STOP "PSSlTex: Syntax error."
             SL=SL+PSSlTex(FRAME,CHR(I+2:I+J))/2
             I=I+J+2
          endif
          IF(CHR(I:I).EQ."_") SL=SL-FRAME%FONT%S/3
       ELSEIF(CHR(I:I).EQ."_")then
          IF(CHR(I+1:I+1).ne."{")then
             IF(CHR(I+1:I+1).ne.Const_backslash)then
                J=1
             ELSE
                J=VERIFY(CHR(I+3:STRLEN),"abcdefghijklmnopqrstuvwxyz")
                IF(J.EQ.0)J=STRLEN-I-1
                J=J+1
             endif
             SL=SL+PSSlSingleTex(FRAME,CHR(I+1:I+J))/2
             I=I+J+1
          ELSE
             J=NEXTBRACKET(CHR(I+2:STRLEN))
             IF(J.EQ.0)STOP "PSPrtTex: Syntax error."
             SL=SL+PSSlTex(FRAME,CHR(I+2:I+J))/2
             I=I+J+2
          endif
          IF(CHR(I:I).EQ."^") SL=SL-FRAME%FONT%S/3
       ELSEIF(CHR(I:I).EQ."{")then
          J=NEXTBRACKET(CHR(I+1:STRLEN))
          IF(J.EQ.0)STOP "PSSlTex: Syntax error."
          SL=SL+PSSlTex(FRAME,CHR(I+1:I+J-1))
          I=I+J+1
       ELSEIF(CHR(I:I).EQ.Const_backslash)then
          J=VERIFY(CHR(I+2:StrLen),"abcdefghijklmnopqrstuvwxyz")
          IF(J.EQ.0)then
             SL=SL+PSSlSingleTex(Frame,CHR(I:StrLen))
             RETURN
          endif
          select case(Chr(I:I+J))
          case("\sqrt")
             IF(CHR(I+5:I+5).EQ."{")then
                SL=SL+(Frame%Font%S+2)*2/3
                I=I+5
             ELSEIF(CHR(I+5:I+5).EQ."[")then
                J=SCAN(CHR(I+6:STRLEN),"]")
                SL=SL+PSSlTEX(FRAME,CHR(I+6:I+J+4))/3
                SL=SL+(Frame%Font%S+2)*2/3
                I=I+J+6
             ELSE
                STOP "PSSlTex: Syntax error."
             endif
          case("\frac")
             IF(CHR(I+5:I+5).ne."{")STOP "PSSlTex: Syntax Error."
             J=NEXTBRACKET(CHR(I+6:STRLEN))
             IF(J.EQ.0)STOP "PSSlTex: Syntax Error."
             FRU=I+J+5
             LU=PSSlTex(Frame,CHR(I+6:FRU-1))
             IF(CHR(FRU+1:FRU+1).ne."{")STOP "PSSlTex: Syntax Error."
             J=NEXTBRACKET(CHR(FRU+2:STRLEN))
             IF(J.EQ.0)STOP "PSSlTex: Syntax Error."
             FRD=FRU+J+1
             LD=PSSlTex(Frame,CHR(FRU+2:FRD-1))
             SL=SL+MAX(LU,LD)
             I=FRD+1
          case("\textit","\textbf")
             IF(CHR(I+7:I+7).ne."{")STOP "PSPrtTex: Syntax Error."
             J=NEXTBRACKET(CHR(I+8:STRLEN))
             IF(J.EQ.0)STOP "PSPrtTex: Syntax Error."
             SL=SL+PSSlSingleTex(FRAME,CHR(I+8:I+J+6))
             I=I+J+8
          case("\text")
             IF(CHR(I+5:I+5).ne."{")STOP "PSPrtTex: Syntax Error."
             J=NEXTBRACKET(CHR(I+6:STRLEN))
             IF(J.EQ.0)STOP "PSPrtTex: Syntax Error."
             SL=SL+PSSlSingleTex(FRAME,CHR(I+6:I+J+4))
             I=I+J+6
          case("\tilde")
             if(Chr(I+6:I+6).ne. "{") stop "PsprtTex: Syntax error."
             J=nextbracket(CHR(I+7:STRLEN))
             IF(J.EQ.0)STOP "PSPrtTex: Syntax Error."
             SL=SL+PSSlSingleTex(FRAME,CHR(I+7:I+J+5))
             I=I+J+7
          case default
             SL=SL+PSSlSingleTex(Frame,CHR(I:I+J))
             I=I+J+1
          end select
       ELSE
          J=SCAN(CHR(I+1:STRLEN),"^_{"//Const_backslash)
          IF(J.EQ.0)then
             SL=SL+PSSlSingleTex(Frame,CHR(I:StrLen))
             RETURN
          endif
          SL=SL+PSSlSingleTex(FRAME,CHR(I:I+J-1))
          I=I+J
       endif
    enddo
  END function PSSlTex

  function PSSlSingle(Frame,C)
    TYPE(PSFrame)Frame
    character C
    integer ASCII
    integer PSSlSingle
    IF(C.EQ."=")then
       PSSlSingle=(Frame%Font%S+1)*3/2
    ELSE
       ASCII=ICHAR(C)
       IF(ASCII.GE.65 .AND. ASCII.LE.90)then
          PSSlSingle=(Frame%Font%S+2)*2/3
       ELSEIF(ASCII.GE.33 .AND. ASCII.LE.122)then
          PSSlSingle=(Frame%Font%S+1)/2
       ELSE
          PSSlSingle=0
       endif
    endif
  END function PSSlSingle

  function PSSlSingleTex(Frame,CHR)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING  CHR
    integer PSSlSingleTex
    integer StrLen,I
    PSSlSingleTex=0
    StrLen=Len_Trim(CHR)
    IF(StrLen.EQ.0)RETURN
    IF(CHR(1:1).ne.Const_backslash)then
       DO I=1,StrLen
          PSSlSingleTex=PSSlSingleTex+PSSlSingle(Frame,CHR(I:I))
       enddo
       RETURN
    ELSEIF(StrLen.EQ.1)then
       PSSlSingleTex=CEILING(Frame%Font%S*0.55+2)
    ELSEIF(CHR(2:StrLen).EQ.Const_backslash)then
       PSSlSingleTex=Frame%Font%S/2
    ELSEIF(CHR(2:StrLen).EQ."exists" .OR. CHR(2:StrLen).EQ."forall" &
         .OR. CHR(2:StrLen).EQ."cong" .OR. CHR(2:StrLen).EQ."le" &
         .OR. CHR(2:StrLen).EQ."infty" .OR. CHR(2:StrLen).EQ."pm" &
         .OR. CHR(2:StrLen).EQ."ge" .OR. CHR(2:StrLen).EQ."times" &
         .OR. CHR(2:StrLen).EQ."propto" .OR. CHR(2:StrLen).EQ."div" &
         .OR. CHR(2:StrLen).EQ."neq" .OR. CHR(2:StrLen).EQ."equiv" &
         .OR. CHR(2:StrLen).EQ."approx" .OR. CHR(2:StrLen).EQ."aleph" &
         .OR. CHR(2:StrLen).EQ."otimes" .OR. CHR(2:StrLen).EQ."oplus" &
         .OR. CHR(2:StrLen).EQ."supset" .OR. CHR(2:StrLen).EQ."supseteq" &
         .OR. CHR(2:StrLen).EQ."subset" .OR. CHR(2:StrLen).EQ."subseteq" &
         .OR. CHR(2:StrLen).EQ."prod" .OR. CHR(2:StrLen).EQ."nabla")then
       PSSlSingleTex=(Frame%Font%S)*2/3+2
    ELSEIF(CHR(2:StrLen).EQ."leftrightarrow" .OR. CHR(2:StrLen).EQ."leftarrow" &
         .OR. CHR(2:StrLen).EQ."rightarrow" .OR. CHR(2:StrLen).EQ."Leftrightarrow" &
         .OR. CHR(2:StrLen).EQ."Leftarrow" .OR. CHR(2:StrLen).EQ."Rightarrow" &
         )then
       PSSlSingleTex=Frame%Font%S
    ELSEIF(CHR(2:StrLen).EQ."sin" .OR. CHR(2:StrLen).EQ."cos" &
         .OR. CHR(2:StrLen).EQ."tan" .OR. CHR(2:StrLen).EQ."cot" &
         .OR. CHR(2:StrLen).EQ."sec" .OR. CHR(2:StrLen).EQ."csc" &
         .OR. CHR(2:StrLen).EQ."exp" .OR. CHR(2:StrLen).EQ."ln" &
         .OR. CHR(2:StrLen).EQ."log" .OR. CHR(2:StrLen).EQ."lg" &
         .OR. CHR(2:StrLen).EQ."sinh" .OR. CHR(2:StrLen).EQ."cosh" &
         .OR. CHR(2:StrLen).EQ."tanh" .OR. CHR(2:StrLen).EQ."coth" &
         .OR. CHR(2:StrLen).EQ."arcsin" .OR. CHR(2:StrLen).EQ."arccos" &
         .OR. CHR(2:StrLen).EQ."arctan" .OR. CHR(2:StrLen).EQ."lim" &
         .OR. CHR(2:StrLen).EQ."max" .OR. CHR(2:StrLen).EQ."min" &
         .OR. CHR(2:StrLen).EQ."ker" .OR. CHR(2:StrLen).EQ."dim" &
         .OR. CHR(2:StrLen).EQ."det" .OR. CHR(2:StrLen).EQ."gcd" &
         .OR. CHR(2:StrLen).EQ."limsup" .OR. CHR(2:StrLen).EQ."liminf" &
         .OR. CHR(2:StrLen).EQ."sup" .OR. CHR(2:StrLen).EQ."inf" &
         .OR. CHR(2:StrLen).EQ."arg" .OR. CHR(2:StrLen).EQ."hom" )then
       PSSlSingleTex=(Frame%Font%S+1)/2*(StrLen-1)
    else
       PSSlSingleTex=PSSlSingle(FRAME,CHR(2:2))
    endif
  END function PSSlSingleTex

  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine PSSetGray(Frame,GRAY)  !!0 black, 1 white
    TYPE(PSFrame)Frame
    real(dl) GRAY
    write(frame%fileunit,'(a)')trim(RGBSTR(GRAY))//" setgray"
  END subroutine PSSetGray

  Subroutine PSDot(Frame)
    type(PSFrame)Frame
    write(frame%fileunit,'(a)')"["//trim(num2str(Frame%CS%dotlen(1)*sqrt(frame%CS%xamp**2+frame%CS%yamp**2),"f"))//" "//trim(num2str(frame%CS%dotlen(2)*sqrt(frame%CS%xamp**2+frame%CS%yamp**2),"f"))//"] 0 setdash"
  end subroutine PSDot

  subroutine PSDash(Frame,N,Ndot)
    TYPE(PSFrame)Frame
    integer,optional:: N
    integer(IB),optional::Ndot
    if(present(N))then
       if(present(Ndot))then
          write(frame%fileunit,'(a)')"["//trim(INT2STR(N))//" "//trim(INT2STR(Ndot))//"] 0 setdash"
       else
          write(frame%fileunit,'(a)')"["//trim(INT2STR(N))//" "//trim(INT2STR(N))//"] 0 setdash"
       endif
    else
       write(frame%fileunit,'(a)')"["//trim(num2str(Frame%CS%dashlen(1),"f"))//" "//trim(Num2str(frame%CS%dashlen(2),"f"))//"] 0 setdash"
    endif
  END subroutine PSDash

  subroutine PSDotDash(Frame,N,NDOT)
    TYPE(PSFrame)Frame
    integer,optional::N
    integer,optional::NDOT
    if(present(N))then
       IF(PRESENT(NDOT))then
          write(frame%fileunit,'(a)')"["//trim(INT2STR(N))//" "//trim(INT2STR(N/2))//" "//trim(INT2STR(NDOT))//" "//trim(INT2STR(N/2))//"] 0 setdash"
       ELSE
          write(frame%fileunit,'(a)')"["//trim(INT2STR(N))//" "//trim(INT2STR(N/2))//" 1 "//trim(INT2STR(N/2))//"] 0 setdash"
       endif
    else
       write(frame%fileunit,'(a)')"["//trim(num2str(Frame%CS%dotdashlen(1),"f"))//" "//trim(num2str(Frame%CS%dotdashlen(2),"f"))//" "// &
            trim(num2str(Frame%CS%dotdashlen(3),"f"))//" "//trim(num2str(Frame%CS%dotdashlen(4),"f"))//"] 0 setdash"
    end if
  END subroutine PSDotDash

  subroutine PSCancelDash(Frame)
    TYPE(PSFrame)Frame
    write(frame%fileunit,'(a)') "[1 0] 0 setdash"
  END subroutine PSCancelDash


  

  subroutine PSSetRGBColor_D(Frame,R,G,B)
    TYPE(PSFrame)Frame
    real(dl) RED,GREEN,BLUE,R,G,B
    Red=Min(Max(R,0._dl),1._dl)
    Green=Min(Max(G,0._dl),1._dl)
    Blue=Min(Max(B,0._dl),1._dl)
    WRITE(Frame%fileunit,"(a)")Trim(RGBStr(RED))//" "//Trim(RGBStr(GREEN))//" "//Trim(RGBStr(BLUE))//" setrgbcolor"
    frame%color%r=red
    frame%color%g=green
    frame%color%b=blue
  END subroutine PSSetRGBColor_D

  Subroutine PSSetRGBColor_color(Frame,Color,dimfac)
    TYPE(PSFRAME)FRAME
    Type(PSColor) Color
    real(dl),optional::dimfac
    real(dl) dimming
    if(present(dimfac))then
       dimming=dimfac
    else
       dimming=0._dl
    endif
    call PsSetRGBColor_D(Frame,Color%R+dimming,Color%G+dimming,Color%B+dimming)
  End Subroutine PSSetRGBColor_color

  subroutine PSSetRGBColor_Str(Frame,Color,dimfac)
    TYPE(PSFRAME)FRAME
    UNKNOWN_STRING  Color
    real(dl) R,G,B
    Real(dl),optional::dimfac
    CALL PSColor2RGB(Color,R,G,B)
    if(present(dimfac))then
       R=R+dimfac
       G=G+dimfac
       B=B+dimfac
    endif
    call PSSetRGBColor_D(Frame,R,G,B)
  END subroutine PSSetRGBColor_Str


  Subroutine PSSetColor(Frame,Color)
    Type(PSFrame) Frame
    UNKNOWN_STRING   Color
    Integer(IB) istart,iend
    real(dl) R,G,B
    if(trim(color).eq."") return
    iend=Len_Trim(Color)
    istart=1
    if(iend.ge.4)then
       if(Color(1:4) .eq. "dot_")then
          call PSDot(Frame)
          istart = 5
       elseif(iend.ge.5)then
          if(Color(1:5) .eq. "dash_")then
             call PsDash(Frame)
             istart=6
          elseif(Color(1:6) .eq. "solid_")then
             call PScanceldash(Frame)
             istart=7
          elseif(iend.ge.8)then
             if(Color(1:8) .eq. "dotdash_")then
                call PSDotDash(Frame)
                istart = 9
             endif
          endif
       endif
    endif
    Call PScolor2RGB(Color(istart:iend),R,G,B)
    call PSSetRGBColor_D(Frame,R,G,B)
  End Subroutine PSSetColor
  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  !!%%%%%%%%%%%%%%%%%%%%% MOVETO LINETO ETC. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine PSCMoveTo_D(Frame,X,Y)
    real(dl) X,Y
    TYPE(PSFrame)Frame
    real(dl) IX,IY
    IX=(X-Frame%Coor%OX)/Frame%Coor%XSCALE+Frame%Coor%Bound%XL
    IY=(Y-Frame%Coor%OY)/Frame%Coor%YSCALE+Frame%Coor%Bound%YL
    call psmoveto(frame,ix,iy)
  END subroutine PSCMoveTo_D

  subroutine PSCMoveTo_S(Frame,X,Y)
    real(dl) X,Y
    TYPE(PSFrame)Frame
    call pscmoveto_D(frame,dble(x),dble(y))
  END subroutine PSCMoveTo_S

  subroutine PSCLineTo(Frame,X,Y)
    real(dl) X,Y
    TYPE(PSFrame)Frame
    real(dl) IX,IY
    CALL PSCOORTRANS(Frame,X,Y,IX,IY)
    CALL PsLineTO(Frame,IX,IY)
  END subroutine PSCLineTo



  subroutine PSCRMoveTo(Frame,X,Y)
    real(dl) X,Y
    TYPE(PSFrame)Frame
    real(dl) IX,IY
    IX=X/Frame%Coor%XSCALE
    IY=Y/Frame%Coor%YSCALE
    CALL PSRMOVETO(Frame,IX,IY)
  END subroutine PSCRMoveTo


  subroutine PSCRLineTo(Frame,X,Y)
    real(dl) X,Y
    TYPE(PSFrame)Frame
    real(dl) IX,IY
    IX=X/Frame%Coor%XSCALE
    IY=Y/Frame%Coor%YSCALE
    CALL PsRLineTo(Frame,IX,IY)
  END subroutine PSCRLineTo

  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine PSIMAGE(Frame,X,Y,GRAY,SZX,SZY)
    TYPE(PSFrame)Frame
    real(dl),dimension(:),INTENT(IN)::X,Y,GRAY
    integer SZX,SZY  !! SIZE OF IMAGE
    integer WD,HT,N,IX,IY,I,J,IGRAY1,IGRAY2
    real(dl) XMIN,XMAX,YMIN,YMAX,XWD,YHT
    integer IGRAY(8*SZX,8*SZY)
    integer IZ1(8*SZX,8*SZY),IZ2(8*SZX,8*SZY)
    WD=SZX*8
    HT=SZY*8
    N=SIZE(X)
    XMIN=MINVAL(X)
    XMAX=MAXVAL(X)
    YMIN=MINVAL(Y)
    YMAX=MAXVAL(Y)
    XWD=XMAX-XMIN
    YHT=YMAX-YMIN
    CALL PSCOORTRANS(Frame,XMIN,YMIN,IX,IY)
    CALL PSGSAVE(Frame)
    CALL PSTRANSLATE(Frame,IX,IY)
    CALL PsScale(Frame,(XMAX-XMIN)/Frame%Coor%XSCALE,(YMAX-YMIN)/Frame%Coor%YSCALE)
    IGRAY=255
    DO I=1,N
       IX=CEILING((X(I)-XMIN)/XWD*WD)
       IF(IX.LE.0)IX=1
       IF(IX.GT.WD)IX=WD
       IY=CEILING((Y(I)-YMIN)/YHT*HT)
       IF(IY.EQ.0)IY=1
       IF(IY.GT.HT)IY=HT
       IGRAY(IX,IY)=MIN(NINT(GRAY(I)*255.),IGRAY(IX,IY))
    enddo
    DO I=1,WD
       DO J=1,HT
          IGRAY1=ISHFT(IGRAY(I,J),-4)
          IGRAY2=IGRAY(I,J)-ISHFT(IGRAY1,4)
          IZ1(I,J)=IGRAY1
          IZ2(I,J)=IGRAY2
       enddo
    enddo
    write(frame%fileunit,'(a)') trim(INT2STR(WD))//" "// &
         trim(INT2STR(HT))//" 8 ["//trim(INT2STR(WD))//" 0 0 "//trim(INT2STR(HT))//" 0 0]"
    CALL PsPush(Frame,"{<")
    DO J=1,HT
       DO I=1,WD
          WRITE(Frame%FILEUNIT,"(2Z1)")IZ1(I,J),IZ2(I,J)
       enddo
    enddo
    CALL PsPush(Frame,">} image")
    CALL PSGRESTORE(Frame)
  END subroutine PSIMAGE



  subroutine PSHEXIMAGE_D(Frame,X,Y,GRAY,SZX,SZY)
    TYPE(PSFrame)Frame
    real(dl),dimension(:),INTENT(IN)::X,Y
    integer SZX,SZY
    integer WD,HT,N,IX,IY,I,J
    real(dl) XMIN,XMAX,YMIN,YMAX,XWD,YHT,GRAY
    integer IGRAY(8*SZX,8*SZY)
    WD=SZX*8
    HT=SZY*8
    N=SIZE(X)
    XMIN=MINVAL(X)
    XMAX=MAXVAL(X)
    YMIN=MINVAL(Y)
    YMAX=MAXVAL(Y)
    XWD=XMAX-XMIN
    YHT=YMAX-YMIN
    CALL PSCOORTRANS(Frame,XMIN,YMIN,IX,IY)
    CALL PSGSAVE(Frame)
    CALL PSTRANSLATE(Frame,IX,IY)
    CALL PsScale(Frame,(XMAX-XMIN)/Frame%Coor%XSCALE,(YMAX-YMIN)/Frame%Coor%YSCALE)
    IGRAY=15
    DO I=1,N
       IX=CEILING((X(I)-XMIN)/XWD*WD)
       IF(IX.LE.0)IX=1
       IF(IX.GT.WD)IX=WD
       IY=CEILING((Y(I)-YMIN)/YHT*HT)
       IF(IY.EQ.0)IY=1
       IF(IY.GT.HT)IY=HT
       IGRAY(IX,IY)=MIN(IGRAY(IX,IY),NINT(GRAY(I)*15.))
    enddo
    write(frame%fileunit,'(a)') trim(INT2STR(WD))//" ", &
         trim(INT2STR(HT))," 4 ["//trim(INT2STR(WD))//" 0 0 "//trim(INT2STR(HT))//" 0 0]"
    CALL PsPush(Frame,"{<")
    DO J=1,HT
       I=1
       DO WHILE(I+15.LT.WD)
          WRITE(Frame%FILEUNIT,"(16Z1)")IGRAY(I:I+15,J)
          I=I+16
       enddo
       WRITE(Frame%FILEUNIT,"(16Z1)")IGRAY(I:WD,J)
    enddo
    CALL PsPush(Frame,">} image")
    CALL PSGRESTORE(Frame)
  END subroutine PSHEXIMAGE_D


  subroutine PSBITIMAGE_D(Frame,X,Y,SZX,SZY)
    TYPE(PSFrame)Frame
    real(dl),dimension(:),INTENT(IN)::X,Y
    integer SZX,SZY
    integer XN,YN   !! 4*XN x 4*YN GRID
    integer WD,HT,N,IX,IY,I,J,K
    real(dl) XMIN,XMAX,YMIN,YMAX,XWD,YHT
    logical BMAP(8*SZX,8*SZY)
    integer BITMAP(2*SZX,8*SZY)
    XN=SZX*2
    YN=SZY*2
    WD=SZX*8
    HT=SZY*8
    N=SIZE(X)
    XMIN=MINVAL(X)
    XMAX=MAXVAL(X)
    YMIN=MINVAL(Y)
    YMAX=MAXVAL(Y)
    XWD=XMAX-XMIN
    YHT=YMAX-YMIN
    CALL PSCOORTRANS(Frame,XMIN,YMIN,IX,IY)
    CALL PSGSAVE(Frame)
    CALL PSTRANSLATE(Frame,IX,IY)
    CALL PsScale(Frame,(XMAX-XMIN)/Frame%Coor%XSCALE,(YMAX-YMIN)/Frame%Coor%YSCALE)
    BMAP=.FALSE.
    DO I=1,N
       IX=CEILING((X(I)-XMIN)/XWD*WD)
       IF(IX.LE.0)IX=1
       IF(IX.GT.WD)IX=WD
       IY=CEILING((Y(I)-YMIN)/YHT*HT)
       IF(IY.EQ.0)IY=1
       IF(IY.GT.HT)IY=HT
       BMAP(IX,IY)=.TRUE.
    enddo
    BITMAP=15
    DO I=1,XN
       DO J=1,HT
          K=I*4
          IF(BMAP(K-3,J))BITMAP(I,J)=BITMAP(I,J)-8
          IF(BMAP(K-2,J))BITMAP(I,J)=BITMAP(I,J)-4
          IF(BMAP(K-1,J))BITMAP(I,J)=BITMAP(I,J)-2
          IF(BMAP(K,J))BITMAP(I,J)=BITMAP(I,J)-1
       enddo
    enddo
    write(frame%fileunit,'(a)')trim(INT2STR(WD))//" "// &
         trim(INT2STR(HT))," 1 ["//trim(INT2STR(WD))//" 0 0 "//trim(INT2STR(HT))//" 0 0]"
    CALL PsPush(Frame,"{<")
    DO J=1,HT
       I=1
       DO WHILE(I+15.LT.XN)
          WRITE(Frame%FILEUNIT,"(16Z1)")BITMAP(I:I+15,J)
          I=I+16
       enddo
       WRITE(Frame%FILEUNIT,"(16Z1)")BITMAP(I:XN,J)
    enddo
    CALL PsPush(Frame,">} image")
    CALL PSGRESTORE(Frame)
  END subroutine PSBITIMAGE_D

 

  subroutine ReadBmpHead(BMPFILE,WIDTH,HEIGHT,HEADLENGTH,BIT)
    UNKNOWN_STRING  BMPFILE
    integer HEADLENGTH,WIDTH,HEIGHT,BIT
    integer I,BN(54)
    character BMPDATA(54)
    OPEN(tmp_file_unit,FILE=BMPFILE,STATUS='OLD',FORM='unformatted',ACCESS='SEQUENTIAL',RECL=1)
    DO I=1,30
       READ(tmp_file_unit)BMPDATA(I)
    enddo
    BN=IACHAR(BMPDATA)
    HEADLENGTH=ISHFT((ISHFT((ISHFT(BN(14),8)+BN(13)),8)+BN(12)),8)+BN(11)
    WIDTH=ISHFT((ISHFT((ISHFT(BN(22),8)+BN(21)),8)+BN(20)),8)+BN(19)
    HEIGHT=ISHFT((ISHFT((ISHFT(BN(26),8)+BN(25)),8)+BN(24)),8)+BN(23)
    BIT=BN(30)*256+BN(29)
    CLOSE(tmp_file_unit)
  END subroutine ReadBmpHead

  subroutine ReadBmpData(BMPFILE,BMPDATA,WIDTH,HEIGHT,HEADLENGTH,BIT)
    UNKNOWN_STRING  BMPFILE
    integer HEADLENGTH,BIT
    integer WIDTH,HEIGHT,PERROW,I,J,PT,PTT
    integer BMPDATA(HEIGHT,WIDTH)
    character TEMP
    integer MOV,BMOV
    PERROW=((WIDTH*BIT-1)/32+1)*4-((WIDTH*BIT-1)/8+1)
    BMPDATA=0
    OPEN(tmp_file_unit,FILE=BMPFILE,STATUS='OLD',FORM='unformatted',ACCESS='SEQUENTIAL',RECL=1)
    DO I=1,HEADLENGTH
       READ(TMP_FILE_UNIT)TEMP
    enddo
    DO I=1,HEIGHT
       PT=1
       DO J=1,WIDTH
          PTT=0
          DO WHILE(PTT.LT.BIT)
             IF(PT.GE.1)then
                READ(TMP_FILE_UNIT)TEMP
                MOV=IACHAR(TEMP)
                PT=-7
             endif
             BMOV=ISHFT(MOV,PT)
             BMPDATA(I,J)=ISHFT(BMPDATA(I,J),1)+BMOV
             PTT=PTT+1
             MOV=ISHFT(ISHFTC(MOV,8+PT,8),-PT-8)
             PT=PT+1
          enddo
       enddo
       DO J=1,PERROW
          READ(TMP_FILE_UNIT)TEMP
       enddo
    enddo
    CLOSE(TMP_FILE_UNIT)
  END subroutine ReadBmpData

  subroutine BMP2PS(Frame,BMPFILE)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING BMPFILE
    integer,dimension(:,:),allocatable::BMAP
    integer WIDTH,HEIGHT,HEADLENGTH,BIT,LENGTH,TIMES
    integer I,J,IZ1,IZ2,IZ3,IZ4,IZ5,IZ6,X1,X2,X3
    CALL READBMPHEAD(BMPFILE,WIDTH,HEIGHT,HEADLENGTH,BIT)
    allocate(BMAP(HEIGHT,WIDTH))
    CALL READBMPDATA(BMPFILE,BMAP,WIDTH,HEIGHT,HEADLENGTH,BIT)
    CALL PSGSAVE(Frame)
    IF(BIT.EQ.24)then
       write(frame%fileunit,*) trim(INT2STR(WIDTH))," ",trim(INT2STR(HEIGHT)),  &
            " scale ",trim(INT2STR(WIDTH))," ",trim(INT2STR(HEIGHT))," 8 [",  &
            trim(INT2STR(WIDTH))," 0 0 ",trim(INT2STR(HEIGHT))," 0 0]"
       CALL PsPush(Frame,"{<")
       DO J=1,HEIGHT
          DO I=1,WIDTH
             X1=ISHFT(BMAP(J,I),-16)
             X2=ISHFT(ISHFTC(BMAP(J,I),8,24),-16)
             X3=ISHFT(ISHFTC(BMAP(J,I),16,24),-16)
             IZ1=ISHFT(X1,-4)
             IZ2=ISHFT(ISHFTC(X1,4,8),-4)
             IZ3=ISHFT(X2,-4)
             IZ4=ISHFT(ISHFTC(X2,4,8),-4)
             IZ5=ISHFT(X3,-4)
             IZ6=ISHFT(ISHFTC(X3,4,8),-4)
             WRITE(Frame%FILEUNIT,"(6Z1)") IZ5,IZ6,IZ3,IZ4,IZ1,IZ2
          enddo
       enddo
       CALL PsPush(Frame,">} false 3 colorimage")
    ELSEIF(BIT.EQ.1)then
       CALL PsSetLineWidth(Frame,1.05d0)
       TIMES=0
       LENGTH=0
       DO J=1,WIDTH
          DO I=1,HEIGHT-1
             IF(BMAP(I,J).eq.0)then
                LENGTH=LENGTH-1
             ELSEIF(LENGTH.LT.0)then
                write(frame%fileunit,*)trim(INT2STR(LENGTH)), &
                     trim(INT2STR(I-1))//trim(INT2STR(J))
                TIMES=TIMES+1
                LENGTH=0
             endif
          enddo
          IF(BMAP(HEIGHT,J).eq.0)then
             LENGTH=LENGTH-1
             write(frame%fileunit,*)trim(INT2STR(LENGTH))// &
                  trim(INT2STR(HEIGHT)),trim(INT2STR(J))
             TIMES=TIMES+1
             LENGTH=0
          ELSEIF(LENGTH.LT.0)then
             write(frame%fileunit,*)trim(INT2STR(LENGTH)), &
                  trim(INT2STR(HEIGHT-1)),trim(INT2STR(J))
             TIMES=TIMES+1
             LENGTH=0
          endif
       enddo
       write(frame%fileunit,'(a)')trim(INT2STR(TIMES)) 
       CALL PsPush(Frame,"{newpath 0.5 sub exch moveto 0 exch rlineto stroke} repeat")
       CALL PsSetGray(Frame,0._dl)
    ELSE
       WRITE(*,*)" Your bmp file is a ",trim(INT2STR(BIT)), &
            " bit file. BMP2PS can only convert 1bit or 24bit bmp files." 
    endif
    CALL PSGRESTORE(Frame)
    deallocate(BMAP)
  END subroutine BMP2PS

  subroutine BlurBMP2PS(Frame,BMPFILE,BLUR)
    TYPE(PSFrame)Frame
    UNKNOWN_STRING BMPFILE
    real(dl) BLUR !! 0<BLUR<1, RECOMMENDED BLUR=0.05~0.2
    integer,PARAMETER::NSIZE=6
    integer,PARAMETER::GSIZE=2**NSIZE
    integer,dimension(:,:),allocatable::BMPDATA
    real(dl),dimension(:,:),allocatable::R
    real(dl),dimension(:,:),allocatable::G
    real(dl),dimension(:,:),allocatable::B
    real(dl),dimension(:,:),allocatable::BGR
    real(dl),dimension(:,:),allocatable::BGG
    real(dl),dimension(:,:),allocatable::BGB
    real(dl),dimension(:,:),allocatable::BGRCOPY
    real(dl),dimension(:,:),allocatable::BGGCOPY
    real(dl),dimension(:,:),allocatable::BGBCOPY
    integer WIDTH,HEIGHT,HEADLENGTH,BIT
    integer WD,HT,RT,GT,BT
    integer I,J,COUNT,CI,CJ
    real(dl) RAVE,GAVE,BAVE,TEMP
    real(dl) SQAVE,NPXL,RBAR,GBAR,BBAR
    integer CGS
    real(dl) CR,CG,CB,TCR,TCG,TCB
    TYPE(PSBOX)BOX
    logical CHANGECOLOR
    CALL READBMPHEAD(BMPFILE,WIDTH,HEIGHT,HEADLENGTH,BIT)
    IF(BIT.ne.24)then
       WRITE(*,*)" Your bmp file is a ",trim(INT2STR(BIT)),  &
            "bit file. BLURBMP2PS can only convert 24bit bmp files." 
       stop
    endif
    CALL PSGSAVE(Frame)
    CALL PSTRANSLATE(Frame,Frame%BOUND%XL,Frame%BOUND%YL)
    TEMP=MIN(REAL(Frame%BOUND%XR-Frame%BOUND%XL)/WIDTH,REAL(Frame%BOUND%YR-Frame%BOUND%YL)/HEIGHT)
    CALL PsScale(Frame,TEMP,TEMP)
    BOX%XL=0
    BOX%YL=0
    BOX%XR=WIDTH
    BOX%YR=HEIGHT
    CALL PSOUTLINEBOX(Frame,BOX)
    CALL PsPush(Frame,"clip")
    HT=(HEIGHT+GSIZE-1)/GSIZE*GSIZE
    WD=(WIDTH+GSIZE-1)/GSIZE*GSIZE
    if(ps_resolution_level.eq.0) stop 'Error in BlurBMP2PS: ps_resolution_level too low'
    CALL PsPush(Frame,"/q {setrgbcolor newpath .5 add moveto 0 rlineto stroke} def")
    CALL PsPush(Frame,"/t {setrgbcolor newpath .5 add moveto 1 0 rlineto stroke} def")
    CALL PsPush(Frame,"/s {setrgbcolor newpath .5 add moveto 2 0 rlineto stroke} def")
    CALL PsPush(Frame,"/u {newpath .5 add moveto 0 rlineto stroke} def")
    CALL PsPush(Frame,"/v {newpath .5 add moveto 1 0 rlineto stroke} def")
    CALL PsPush(Frame,"/w {newpath .5 add moveto 2 0 rlineto stroke} def")
    allocate(BMPDATA(HEIGHT,WIDTH))
    allocate(R(HT,WD))
    allocate(G(HT,WD))
    allocate(B(HT,WD))
    allocate(BGR(HT,WD))
    allocate(BGG(HT,WD))
    allocate(BGB(HT,WD))
    allocate(BGRCOPY(HT,WD))
    allocate(BGGCOPY(HT,WD))
    allocate(BGBCOPY(HT,WD))
    CALL READBMPDATA(BMPFILE,BMPDATA,WIDTH,HEIGHT,HEADLENGTH,BIT)
    DO J=1,HEIGHT
       DO I=1,WIDTH
          BT=ISHFT(BMPDATA(J,I),-16)
          GT=ISHFT(ISHFTC(BMPDATA(J,I),8,24),-16)
          RT=ISHFT(ISHFTC(BMPDATA(J,I),16,24),-16)
          R(J,I)=RT/255.
          G(J,I)=GT/255.
          B(J,I)=BT/255.
       enddo
    enddo
    RAVE=SUM(R(1:HEIGHT,1:WIDTH))/WIDTH/HEIGHT
    GAVE=SUM(G(1:HEIGHT,1:WIDTH))/WIDTH/HEIGHT
    BAVE=SUM(B(1:HEIGHT,1:WIDTH))/WIDTH/HEIGHT
    TEMP=0.
    COUNT=0
    DO I=1,HEIGHT
       DO J=1,WIDTH
          IF(ABS(R(I,J)-RAVE).LT.BLUR)then
             TEMP=TEMP+R(I,J)
             COUNT=COUNT+1
          endif
       enddo
    enddo
    IF(COUNT.ne.0)RAVE=TEMP/COUNT
    TEMP=0.
    COUNT=0
    DO I=1,HEIGHT
       DO J=1,WIDTH
          IF(ABS(G(I,J)-GAVE).LT.BLUR)then
             TEMP=TEMP+G(I,J)
             COUNT=COUNT+1
          endif
       enddo
    enddo
    IF(COUNT.ne.0)GAVE=TEMP/COUNT
    TEMP=0.
    COUNT=0
    DO I=1,HEIGHT
       DO J=1,WIDTH
          IF(ABS(B(I,J)-BAVE).LT.BLUR)then
             TEMP=TEMP+B(I,J)
             COUNT=COUNT+1
          endif
       enddo
    enddo
    IF(COUNT.ne.0)BAVE=REAL(TEMP)/COUNT
    DO J=HEIGHT+1,HT
       R(J,:)=RAVE
       G(J,:)=GAVE
       B(J,:)=BAVE
    enddo
    DO I=WIDTH+1,WD
       R(:,I)=RAVE
       G(:,I)=GAVE
       B(:,I)=BAVE
    enddo
    RAVE=RAVE
    GAVE=GAVE
    BAVE=BAVE
    CALL PsPush(Frame,"0 setlinecap")
    WRITE(Frame%FILEUNIT,"(I5,A21)")HT," setlinewidth newpath"
    CALL PsMoveTo(Frame,0,HT/2)
    CALL PsRLineTo(Frame,WD,0)
    CALL PsSetRGBColor(Frame,RAVE,GAVE,BAVE)
    CALL PSSTROKE(Frame)
    CALL PsPush(Frame,"1 setlinewidth")
    CGS=GSIZE
    SQAVE=BLUR**2
    BGRCOPY=RAVE
    BGGCOPY=GAVE
    BGBCOPY=BAVE
    CR=RAVE
    CG=GAVE
    CB=BAVE
    CI=0
    CJ=0
    CHANGECOLOR=.TRUE.
    DO WHILE(CGS.GE.1)
       CALL PSGSAVE(Frame)
       WRITE(Frame%FILEUNIT,"(I3,A1,I3,A6)")CGS," ",CGS," scale"
       NPXL=REAL(CGS*CGS)
       COUNT=0
       DO I=1,HT/CGS
          CI=I
          DO J=1,WD/CGS
             BGR(I,J)=SUM(R((I-1)*CGS+1:I*CGS,(J-1)*CGS+1:J*CGS))/NPXL
             BGG(I,J)=SUM(G((I-1)*CGS+1:I*CGS,(J-1)*CGS+1:J*CGS))/NPXL
             BGB(I,J)=SUM(B((I-1)*CGS+1:I*CGS,(J-1)*CGS+1:J*CGS))/NPXL
             RBAR=SUM((BGR(I,J)-R((I-1)*CGS+1:I*CGS,(J-1)*CGS+1:J*CGS))**2)/NPXL
             GBAR=SUM((BGG(I,J)-G((I-1)*CGS+1:I*CGS,(J-1)*CGS+1:J*CGS))**2)/NPXL
             BBAR=SUM((BGB(I,J)-B((I-1)*CGS+1:I*CGS,(J-1)*CGS+1:J*CGS))**2)/NPXL
             IF((ABS(BGR(I,J)-BGRCOPY((I+1)/2,(J+1)/2)).LT.BLUR  &
                  .AND. ABS(BGG(I,J)-BGGCOPY((I+1)/2,(J+1)/2)).LT.BLUR  &
                  .AND. ABS(BGB(I,J)-BGBCOPY((I+1)/2,(J+1)/2)).LT.BLUR)  &
                  .OR.  RBAR.GT.SQAVE .OR. GBAR.GT.SQAVE .OR. BBAR.GT.SQAVE) then
                BGR(I,J)=BGRCOPY((I+1)/2,(J+1)/2)
                BGG(I,J)=BGGCOPY((I+1)/2,(J+1)/2)
                BGB(I,J)=BGBCOPY((I+1)/2,(J+1)/2)
                CALL DRAWBMPBLOCK(Frame,COUNT,CI,CJ,CR,CG,CB,CHANGECOLOR)
             ELSE
                TCR=NINT(BGR(I,J)*1000.)/1000.
                TCG=NINT(BGG(I,J)*1000.)/1000.
                TCB=NINT(BGB(I,J)*1000.)/1000.
                IF(COUNT.GT.0)then
                   IF(ABS(TCR-CR).GE.0.001 .OR. ABS(TCG-CG).GE.0.001 .OR. ABS(TCB-CB).GE.0.001)then 
                      CALL DRAWBMPBLOCK(Frame,COUNT,CI,CJ,CR,CG,CB,CHANGECOLOR)
                      CR=TCR
                      CG=TCG
                      CB=TCB
                      CJ=J
                      CHANGECOLOR=.TRUE.
                   endif
                   COUNT=COUNT+1
                ELSE
                   IF(ABS(TCR-CR).GE.0.001 .OR. ABS(TCG-CG).GE.0.001 .OR. ABS(TCB-CB).GE.0.001)then
                      CHANGECOLOR=.TRUE.
                   ELSE
                      CHANGECOLOR=.FALSE.
                   endif
                   CR=TCR
                   CG=TCG
                   CB=TCB
                   CJ=J
                   COUNT=1
                endif
             endif
          enddo
          CALL DRAWBMPBLOCK(Frame,COUNT,CI,CJ,CR,CG,CB,CHANGECOLOR)
       enddo
       BGRCOPY(1:HT/CGS,1:WD/CGS)=BGR(1:HT/CGS,1:WD/CGS)
       BGGCOPY(1:HT/CGS,1:WD/CGS)=BGG(1:HT/CGS,1:WD/CGS)
       BGBCOPY(1:HT/CGS,1:WD/CGS)=BGB(1:HT/CGS,1:WD/CGS)
       CGS=CGS/2
       CALL PSGRESTORE(Frame)
    enddo
    CALL PSGRESTORE(Frame)
    deallocate(BGR)
    deallocate(BGG)
    deallocate(BGB)
    deallocate(BGRCOPY)
    deallocate(BGGCOPY)
    deallocate(BGBCOPY)
    deallocate(BMPDATA)
    deallocate(R)
    deallocate(G)
    deallocate(B)
  END subroutine BlurBMP2PS

  subroutine drawbmpblock(Frame,COUNT,CI,CJ,CR,CG,CB,CHANGECOLOR)
    TYPE(PSFrame)Frame
    integer COUNT,CI,CJ
    logical CHANGECOLOR
    real(dl) CR,CG,CB
    IF(COUNT.GT.0)then
       IF(CHANGECOLOR)then
          IF(COUNT.EQ.1)then
             write(frame%fileunit,'(a)') trim(INT2STR(CJ-1))//" "//trim(INT2STR(CI-1))//trim(prtrgb(CR,CG,CB))//" t" 
          ELSEIF(COUNT.EQ.2)then
             write(frame%fileunit,'(a)') trim(INT2STR(CJ-1))//" "//trim(INT2STR(CI-1))//trim(prtrgb(CR,CG,CB))//" s" 
          ELSE
             write(frame%fileunit,'(a)') trim(INT2STR(COUNT))//" "//trim(INT2STR(CJ-1))//" "//trim(INT2STR(CI-1))//trim(prtrgb(CR,CG,CB))//" q" 
          endif
       ELSE
          IF(COUNT.EQ.1)then
             write(frame%fileunit,'(a)') trim(INT2STR(CJ-1))//" "//trim(INT2STR(CI-1))//" v" 
          ELSEIF(COUNT.EQ.2)then
             write(frame%fileunit,'(a)') trim(INT2STR(CJ-1))//" "//trim(INT2STR(CI-1))//" w" 
          ELSE
             write(frame%fileunit,'(a)') trim(INT2STR(COUNT))//" "//trim(INT2STR(CJ-1))//" "//trim(INT2STR(CI-1))//" u" 
          endif
       endif
       COUNT=0
    endif
  END subroutine drawbmpblock


  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine PSPush(Frame,CHR)  !!push an element into the stack
    TYPE(PSFrame)Frame
    UNKNOWN_STRING CHR
    WRITE(Frame%FileUnit,"(A)")Trim(CHR)
  END subroutine PSPush

  subroutine PSPop(Frame) !! discard top element
    TYPE(PSFrame)Frame
    CALL PsPush(FRAME,"pop")
  END subroutine PSPop

  subroutine PSClear(Frame) !!discard all element
    TYPE(PSFrame)Frame
    CALL PsPush(FRAME,"clear")
  END subroutine PSClear

  subroutine PSLoadDict(Frame)
    TYPE(PSFrame)Frame
    LONG_STRING InLine
    If(Trim(Frame%Dict).eq."")return
    open(tmp_file_unit,file=trim(Frame%Dict),status="old",form="formatted",err=200)
    print*,"loading dictionary file ",trim(frame%dict)
    do 
       read(tmp_file_unit,"(A)",end=100) InLine
       call pspush(Frame,InLine)
    enddo
100 close(tmp_file_unit)
    return
200 print*,"Dictonary file ",Trim(Frame%Dict)," does not exist. Importing failed."
  END subroutine PSLoadDict


  subroutine PSLoadMacro(Frame)
    TYPE(PSFrame)Frame
    CALL PsPush(FRAME,"/ft {findfont exch scalefont setfont} def")
    CALL PsPush(Frame,"/nft {/Times-Roman ft} def")
    CALL PsPush(Frame,"/gft {/Symbol ft} def")
    CALL PsPush(FRAME,"/calft {/CMSY10 ft} def")
    CALL PsPush(Frame,"/gift {/Symbol ft} def")
    !! some printer does not support symbol-italic, so we just use symbol here.
!!    CALL PsPush(Frame,"/ift {/Times-Roman ft} def")
    CALL PsPush(Frame,"/ift {/Times-Italic ft} def")
    CALL PsPush(Frame,"/bft {/Times-Bold ft} def")
    CALL PsPush(Frame,"/bift {/Times-BoldItalic ft} def")
    CALL PsPush(Frame,"/m {moveto} def")
    CALL PsPush(Frame,"/s {rmoveto} def")
    CALL PsPush(Frame,"/l {lineto} def")
    CALL PsPush(Frame,"/r {rlineto} def")
    CALL PsPush(Frame,"/sl {stringwidth pop} def")
    CALL PsPush(Frame,"/larger {2 copy gt {pop} {exch pop} ifelse} def")
    CALL PsPush(Frame,"/prtcl {dup sl exch show dup -1 mul 0 s} def")
    CALL PsPush(Frame,"/prtcc {dup sl -2 div 0 s show} def")
    CALL PsPush(Frame,"/prtcr {dup sl -1 mul 0 s show} def")
    CALL PsPush(Frame,"/prtnl {20 string cvs show} def")
    CALL PsPush(Frame,"/prtnc {20 string cvs prtcc} def")
    CALL PsPush(Frame,"/prtnr {20 string cvs prtcr} def")
    CALL PsPush(Frame,"/frac1st {2 copy sl exch sl exch 2 copy larger 2 add dup currentpoint")
    CALL PsPush(Frame,"gsave newpath m 0 r stroke grestore exch} def")
    CALL PsPush(Frame,"/frac2nd {copy sub 2 div 6 s pop show pop} def")
    CALL PsPush(Frame,"/frac3rd {2 copy add -2 div -10 s pop 4 copy sub -2 div (W) sl")
    CALL PsPush(Frame,"-1 mul s pop show sub -2 div (W) sl 1 add s pop pop} def")
    CALL PsPush(Frame,"/prtfrac {2 copy sl exch sl exch 2 copy larger 2 add dup 1 (o) sl 2 div")
    CALL PsPush(Frame,"s currentpoint gsave m newpath 0 r stroke grestore exch 5 frac2nd")
    CALL PsPush(Frame,"2 copy add -2 div -4 s pop 4 copy sub -2 div (W) sl -1 mul s pop")
    CALL PsPush(Frame,"show sub -2 div (W) sl 1 add (o) sl 2 div sub s pop pop} def")
    CALL PsPush(Frame,"/fracm {frac1st 5 frac2nd frac3rd} def")
    CALL PsPush(Frame,"/fracms {frac1st 6 frac2nd gft frac3rd pop} def")
    CALL PsPush(Frame,"/fracmn {frac1st 6 frac2nd nft frac3rd pop} def")
    CALL PsPush(Frame,"/fracmi {frac1st 6 frac2nd ift frac3rd pop} def")
    CALL PsPush(Frame,"/fracmb {frac1st 6 frac2nd bft frac3rd pop} def")
    CALL PsPush(Frame,"/fracmbi {frac1st 6 frac2nd bift frac3rd pop} def")
    CALL PsPush(Frame,"/psbar {dup currentpoint gsave newpath m -2 (W) sl s sl 4 add")
    CALL PsPush(Frame,"0 r stroke grestore show} def")
    CALL PsPush(Frame,"/inch {144 mul} def")  
    CALL PsPush(Frame,"/pscircle {3 copy pop m 0 360 arc} def") !!stack: x,y,r
    CALL PsPush(Frame,"/psytip {-3 0 s 6 0 r -3 0 s} def")
    CALL PsPush(Frame,"/psyeb {dup 0 exch s psytip dup -2. mul 0 exch r psytip 0 exch s} def")
    CALL PsPush(Frame,"/psysbar {m psyeb} def")
    CALL PsPush(Frame,"/psybar {2 copy 3 pscircle psysbar} def") !!stack: sigma x y
    CALL PsPush(Frame,"/psxtip {0 -3 s 0 6 r 0 -3 s} def")
    CALL PsPush(Frame,"/psxeb {dup 0 s psxtip dup -2. mul 0 r psxtip 0 s} def")
    CALL PsPush(Frame,"/psxsbar {m psxeb} def")
    CALL PsPush(Frame,"/psxbar {2 copy 3 pscircle psxsbar} def") !!stack: sigma,x,y
    CALL PsPush(Frame,"/pscsbar {m psyeb psxeb} def") 
    CALL PsPush(Frame,"/pscbar {2 copy 3 pscircle  pscsbar} def") !!stack: sigmax,sigmay,x,y
    CALL PsPush(FRAME,"/pscross {dup 0 exch s dup -2 mul 0 exch")
    CALL PsPush(FRAME,"r 0 exch s dup -1 mul 0 s 2 mul 0 r} def") 
    if(ps_resolution_level.eq.0)then
       call pspush(Frame,"/a {-1 0} def")
       call pspush(Frame,"/d {1 0} def")
       call pspush(Frame,"/x {0 -1} def")
       call pspush(Frame,"/w {0 1} def")
       call pspush(Frame,"/q {-1 1} def")
       call pspush(Frame,"/e {1 1} def")
       call pspush(Frame,"/z {-1 -1} def")
       call pspush(Frame,"/c {1 -1} def")
    endif
  END subroutine PSLoadMacro

  subroutine PSPlotPoint_D(FRAME,X,Y,PointSIZE,STYLE)
    TYPE(PSFRAME)FRAME
    real(dl)X,Y
    integer,optional::PointSIZE
    UNKNOWN_STRING  ,optional::STYLE
    integer IX,IY,S,i
    type(pscolor) color
    CALL PSGSAVE(FRAME)
    CALL PsNewPath(FRAME)
    IF(PRESENT(PointSIZE))then
       S=MAX(PointSIZE,1)
    ELSE
       S=3
    endif
    IF(PRESENT(STYLE))then
       Select case(trim(STYLE))
       case("circle","Circle","CIRCLE")
          CALL PSCOORTRANS(FRAME,X,Y,IX,IY)
          CALL PsPush(FRAME,trim(INT2STR(IX))//" "//trim(INT2STR(IY))//" "//trim(INT2STR(S))//" pscircle")
          CALL PSSTROKE(FRAME)
       case("solidcircle","SolidCircle","SOLIDCIRCLE","Solidcircle")
          CALL PSCOORTRANS(FRAME,X,Y,IX,IY)
          CALL PsPush(FRAME,trim(INT2STR(IX))//" "//trim(INT2STR(IY))//" "//trim(INT2STR(S))//" pscircle")
          CALL PSCLOSEPATH(FRAME)
          CALL PSFILL(FRAME)

       case("ball","Ball","BALL")
          CALL PSCOORTRANS(FRAME,X,Y,IX,IY)
          color=frame%color
          do i=1,5
             call pssetrgbcolor(frame,color,(i-1)*(1.-min(color%r,color%g,color%b))/5._dl)
             CALL PsPush(FRAME,trim(INT2STR(IX))//" "//trim(INT2STR(IY))//" "//trim(Num2str(s/5._dl*(6-i)))//" pscircle")
             CALL PSCLOSEPATH(FRAME)
             CALL PSFILL(FRAME)
             if(i.lt.5)call psnewpath(frame)
          enddo
          call pssetrgbcolor(frame,color)
       case("box","Box","BOX")
          CALL PsCMoveTo(FRAME,X,Y)
          CALL PSRMOVETO(FRAME,-S,-S)
          CALL PsRLineTo(FRAME,2*S,0)
          CALL PsRLineTo(FRAME,0,2*S)
          CALL PsRLineTo(FRAME,-2*S,0)
          CALL PSCLOSEPATH(FRAME)
          CALL PSSTROKE(FRAME)
       case("solidbox","SolidBox","SOLIDBOX","Solidbox")
          CALL PsCMoveTo(FRAME,X,Y)
          CALL PSRMOVETO(FRAME,-S,-S)
          CALL PsRLineTo(FRAME,2*S,0)
          CALL PsRLineTo(FRAME,0,2*S)
          CALL PsRLineTo(FRAME,-2*S,0)
          CALL PSCLOSEPATH(FRAME)
          CALL PSFILL(FRAME)
       case("triangle","Triangle","TRIANGLE")
          CALL PsCMoveTo(FRAME,X,Y)
          CALL PSRMOVETO(FRAME,-INT(S*1.732),-S)
          CALL PsRLineTo(FRAME,2*INT(S*1.732),0)
          CALL PsRLineTo(FRAME,-INT(S*1.732),3*S)
          CALL PSCLOSEPATH(FRAME)
          CALL PSSTROKE(FRAME)
       case("solidtriangle","SolidTriangle","Solidtriangle","SOLIDTRIANGLE")
          CALL PsCMoveTo(FRAME,X,Y)
          CALL PSRMOVETO(FRAME,-INT(S*1.732),-S)
          CALL PsRLineTo(FRAME,2*INT(S*1.732),0)
          CALL PsRLineTo(FRAME,-INT(S*1.732),3*S)
          CALL PSCLOSEPATH(FRAME)
          CALL PSFILL(FRAME)
       case("star","STAR","Star")
          CALL PsCMoveTo(FRAME,X,Y)
          CALL PSRMOVETO(FRAME,-Frame%Font%S/2, - Frame%Font%S/2)
          call psprtN(Frame,"*")
       case default
          CALL PsCMoveTo(FRAME,X,Y)
          CALL PsPush(FRAME,trim(INT2STR(S)))
          CALL PsPush(FRAME,"dup pscross")
          CALL PSSTROKE(FRAME)
       END Select
    ELSE
       CALL PsCMoveTo(FRAME,X,Y)
       CALL PsPush(FRAME,trim(INT2STR(S)))
       CALL PsPush(FRAME,"dup pscross")
       CALL PSSTROKE(FRAME)
    endif
    CALL PSGRESTORE(FRAME)
  END subroutine PSPlotPoint_D


  subroutine PSPLOTPOINT_DV(FRAME,X,Y,PointSIZE,STYLE)
    TYPE(PSFRAME)FRAME
    real(dl),dimension(:),intent(in)::X,Y
    integer,optional::PointSIZE
    UNKNOWN_STRING  ,optional::STYLE
    integer(IB) N,i
    N=GetDim("PSPlotpoint",size(x),size(y))
    IF(PRESENT(PointSIZE))then
       IF(PRESENT(STYLE))then
          do i=1,N
             CALL PSPLOTPOINT_D(FRAME,x(i),y(i),PointSIZE,STYLE)
          enddo
       ELSE
          do i=1,N
             CALL PSPLOTPOINT_D(FRAME,x(i),y(i),PointSIZE)
          enddo
       endif
    ELSE
       do i=1,N
          CALL PSPLOTPOINT_D(FRAME,x(i),y(i))
       enddo
    endif
  END subroutine PSPLOTPOINT_DV



  subroutine PSsetfont(Frame,Font,SCALE)
    !!Font="Times-Roman","Symbol","Times-Bold","Times-Italic","Times-BoldItalic","Helvetica",etc
    TYPE(PSFrame)Frame
    UNKNOWN_STRING ,optional::Font
    integer,optional::SCALE
    IF(PRESENT(Font).AND.PRESENT(SCALE))then
       IF(trim(Frame%Font%T).EQ.Font .AND. Frame%Font%S.EQ.SCALE)RETURN
    endif
    IF(PRESENT(Font))Frame%Font%T=Font
    IF(PRESENT(SCALE))Frame%Font%S=SCALE
    select case(Trim(Frame%Font%T))
    case("Times-Roman")
       write(frame%fileunit,'(a)')trim(INT2STR(Frame%Font%S))//" nft"
    case("Symbol")
       write(frame%fileunit,'(a)')trim(INT2STR(Frame%Font%S))//" gft"
    case("Symbol-Italic")
       write(frame%fileunit,'(a)')trim(INT2STR(Frame%Font%S))//" gift"
    case("Times-Italic")
       write(frame%fileunit,'(a)')trim(INT2STR(Frame%Font%S))//" ift"
    case("Times-Bold")
       write(frame%fileunit,'(a)')trim(INT2STR(Frame%Font%S))//" bft"
    case("Times-BoldItalic")
       write(frame%fileunit,'(a)')trim(INT2STR(Frame%Font%S))//" bift"
    case("Cal")
       write(frame%fileunit,'(a)')trim(Int2Str(Frame%Font%S))//" calft"
    case default
       write(frame%fileunit,'(a)')"/",trim(Frame%Font%T)//" findfont"
       write(frame%fileunit,'(a)')trim(INT2STR(Frame%Font%S))//" scalefont setfont"
    END select
  END subroutine PSsetfont

  subroutine PSSymbolSwitch(Frame,FSIZE)
    TYPE(PSFrame)Frame
    integer,optional::FSIZE
    logical,Save::RESTORE=.FALSE.
    SHORT_STRING,Save::FontType="Times-Roman"
    integer,Save::FontSize=24
    IF(.not.RESTORE)then
       IF(trim(Frame%Font%T).EQ."Symbol")then
          IF(PRESENT(FSize))then
             IF(FSIZE.EQ.Frame%Font%S)RETURN
          ELSE
             RETURN
          endif
       endif
       FontType=trim(Frame%Font%T)
       FontSize=Frame%Font%S
       IF(PRESENT(FSize))then
          write(frame%fileunit,'(a)')trim(INT2STR(FSIZE))//" gft"
          Frame%Font%S=FSIZE
       ELSE
          write(frame%fileunit,'(a)')trim(INT2STR(FontSize))//" gft"
       endif
       Frame%Font%T="Symbol"
       RESTORE=.TRUE.
    ELSE
       CALL PSsetfont(Frame,FontType,FontSize)
       RESTORE=.FALSE.
    endif
  END subroutine PSSymbolSwitch

  subroutine PSSymbolItalicSwitch(Frame,FSIZE)
    TYPE(PSFrame)Frame
    integer,optional::FSIZE
    logical,SAVE::RESTORE=.FALSE.
    SHORT_STRING ,SAVE::FontType="Times-Roman"
    integer,SAVE::FontSize=24
    IF(.not.RESTORE)then
       IF(trim(Frame%Font%T).EQ."Symbol-Italic")then
          IF(PRESENT(FSize))then
             IF(FSIZE.EQ.Frame%Font%S)RETURN
          ELSE
             RETURN
          endif
       endif
       FontType=trim(Frame%Font%T)
       FontSize=Frame%Font%S
       IF(PRESENT(FSize))then
          write(frame%fileunit,'(a)')trim(INT2STR(FSIZE))//" gift"
          Frame%Font%S=FSIZE
       ELSE
          write(frame%fileunit,'(a)')trim(INT2STR(FontSize))//" gift"
       endif
       Frame%Font%T="Symbol-Italic"
       RESTORE=.TRUE.
    ELSE
       CALL PSsetfont(Frame,FontType,FontSize)
       RESTORE=.FALSE.
    endif
  END subroutine PSSymbolItalicSwitch

  subroutine PSItalicSwitch(Frame,FSIZE)
    TYPE(PSFrame)Frame
    integer,optional::FSIZE
    logical,SAVE::RESTORE=.FALSE.
    SHORT_STRING ,SAVE::FontType="Times-Roman"
    integer,SAVE::FontSize=24
    IF(.not.RESTORE)then
       IF(trim(Frame%Font%T).EQ."Times-Italic")then
          IF(PRESENT(FSize))then
             IF(FSIZE.EQ.Frame%Font%S)RETURN
          ELSE
             RETURN
          endif
       endif
       FontType=trim(Frame%Font%T)
       FontSize=Frame%Font%S
       IF(PRESENT(FSize))then
          write(frame%fileunit,'(a)')trim(INT2STR(FSIZE))//" ift"
          Frame%Font%S=FSIZE
       ELSE
          write(frame%fileunit,'(a)')trim(INT2STR(FontSize))," ift"
       endif
       Frame%Font%T="Times-Italic"
       RESTORE=.TRUE.
    ELSE
       CALL PSsetfont(Frame,FontType,FontSize)
       RESTORE=.FALSE.
    endif
  END subroutine PSItalicSwitch

  subroutine PSBoldSwitch(Frame,FSize)
    TYPE(PSFrame)Frame
    integer,optional::FSIZE
    logical,SAVE::RESTORE=.FALSE.
    SHORT_STRING ,SAVE::FontType="Times-Roman"
    integer,SAVE::FontSize=24
    IF(.not.RESTORE)then
       IF(trim(Frame%Font%T).EQ."Times-Bold")then
          IF(PRESENT(FSize))then
             IF(FSIZE.EQ.Frame%Font%S)RETURN
          ELSE
             RETURN
          endif
       endif
       FontType=trim(Frame%Font%T)
       FontSize=Frame%Font%S
       IF(PRESENT(FSize))then
          write(frame%fileunit,'(a)')trim(INT2STR(FSIZE))," bft"
          Frame%Font%S=FSIZE
       ELSE
          write(frame%fileunit,'(a)')trim(INT2STR(FontSize))," bft"
       endif
       Frame%Font%T="Times-Bold"
       RESTORE=.TRUE.
    ELSE
       CALL PSsetfont(Frame,FontType,FontSize)
       RESTORE=.FALSE.
    endif
  END subroutine PSBoldSwitch

  subroutine PSBoldItalicSwitch(Frame,FSize)
    TYPE(PSFrame)Frame
    integer,optional::FSIZE
    logical,SAVE::RESTORE=.FALSE.
    SHORT_STRING ,SAVE::FontType="Times-Roman"
    integer,SAVE::FontSize=24
    IF(.not.RESTORE)then
       IF(trim(Frame%Font%T).EQ."Times-BoldItalic")then
          IF(PRESENT(FSize))then
             IF(FSIZE.EQ.Frame%Font%S)RETURN
          ELSE
             RETURN
          endif
       endif
       FontType=trim(Frame%Font%T)
       FontSize=Frame%Font%S
       IF(PRESENT(FSize))then
          write(frame%fileunit,'(a)')trim(INT2STR(FSIZE))," bift"
          Frame%Font%S=FSIZE
       ELSE
          write(frame%fileunit,'(a)')trim(INT2STR(FontSize))//" bift"
       endif
       Frame%Font%T="Times-BoldItalic"
       RESTORE=.TRUE.
    ELSE
       CALL PSsetfont(Frame,FontType,FontSize)
       RESTORE=.FALSE.
    endif
  END subroutine PSBoldItalicSwitch

  subroutine PSNormalSwitch(Frame,FSIZE)
    TYPE(PSFrame)Frame
    integer,optional::FSIZE
    logical,SAVE::RESTORE=.FALSE.
    SHORT_STRING ,SAVE::FontType="Times-Roman"
    integer,SAVE::FontSize=24
    IF(.not.RESTORE)then
       IF(trim(Frame%Font%T).EQ."Times-Roman")then
          IF(PRESENT(FSize))then
             IF(FSIZE.EQ.Frame%Font%S)RETURN
          ELSE
             RETURN
          endif
       endif
       FontType=trim(Frame%Font%T)
       FontSize=Frame%Font%S
       IF(PRESENT(FSize))then
          write(frame%fileunit,'(a)')trim(INT2STR(FSIZE))//" nft"
          Frame%Font%S=FSIZE
       ELSE
          write(frame%fileunit,'(a)')trim(INT2STR(FontSize))//" nft"
       endif
       Frame%Font%T="Times-Roman"
       RESTORE=.TRUE.
    ELSE
       CALL PSsetfont(Frame,FontType,FontSize)
       RESTORE=.FALSE.
    endif
  END subroutine PSNormalSwitch

  subroutine PSCalSwitch(Frame,FSIZE)
    TYPE(PSFrame)Frame
    integer,optional::FSIZE
    logical,SAVE::RESTORE=.FALSE.
    SHORT_STRING ,SAVE::FontType="Times-Roman"
    integer,SAVE::FontSize=24
    IF(.not.RESTORE)then
     !  IF(trim(Frame%Font%T).EQ."Times-Roman")then
     !     IF(PRESENT(FSize))then
     !        IF(FSIZE.EQ.Frame%Font%S)RETURN
     !     ELSE
     !        RETURN
     !     endif
     !  endif
       FontType=trim(Frame%Font%T)
       FontSize=Frame%Font%S
       IF(PRESENT(FSize))then
          write(frame%fileunit,'(a)')trim(INT2STR(FSIZE))//" calft"
          Frame%Font%S=FSIZE
       ELSE
          write(frame%fileunit,'(a)')trim(INT2STR(FontSize))//" calft"
       endif
       Frame%Font%T="Cal"
       RESTORE=.TRUE.
    ELSE
       CALL PSsetfont(Frame,FontType,FontSize)
       RESTORE=.FALSE.
    endif
  END subroutine PSCalSwitch


  subroutine PSTranslate(Frame,X,Y)
    TYPE(PSFrame)Frame
    integer X,Y
    write(frame%fileunit,'(a)') trim(INT2STR(X))//" "//trim(INT2STR(Y))//" translate"
  END subroutine PSTranslate

  subroutine PSMoveTo_I(Frame,X,Y)
    TYPE(PSFrame)Frame
    integer X,Y
    write(frame%fileunit,'(a)') trim(INT2STR(X))//" "//trim(INT2STR(Y))//" m"
  END subroutine PSMoveTo_I

  subroutine PSMoveTo_D(Frame,X,Y)
    TYPE(PSFrame)Frame
    real(dl) X,Y
    write(frame%fileunit,'(a)') trim(PS_Num_String(X))//" "//trim(PS_Num_String(Y))//" m"
  END subroutine PSMoveTo_D


  subroutine PSRMoveTo_I(Frame,X,Y)
    TYPE(PSFrame)Frame
    integer X,Y
    write(frame%fileunit,'(a)') trim(INT2STR(X))//" "//trim(INT2STR(Y))//" s"
  END subroutine PSRMoveTo_I

  subroutine PSRMoveTo_D(Frame,X,Y)
    TYPE(PSFrame)Frame
    Real(Dl) X,Y
    write(frame%fileunit,'(a)') trim(PS_Num_String(X))//" "//trim(PS_Num_String(Y))//" s"
  END subroutine PSRMoveTo_D


  subroutine PSLineTo_I(Frame,X,Y)
    TYPE(PSFrame)Frame
    integer X,Y
    write(frame%fileunit,'(a)') trim(INT2STR(X))//" "//trim(INT2STR(Y))//" l"
  END subroutine PSLineTo_I

  subroutine PSLineTo_D(Frame,X,Y)
    TYPE(PSFrame)Frame
    Real(Dl) X,Y
    write(frame%fileunit,'(a)') trim(PS_Num_string(x))//" "//trim(PS_Num_String(Y))//" l"
  END subroutine PSLineTo_D



  subroutine PSRLineTo_I(Frame,X,Y)
    TYPE(PSFrame)Frame
    integer X,Y
    write(frame%fileunit,'(a)') trim(INT2STR(X))//" "//trim(INT2STR(Y))//" r"
  END subroutine PSRLineTo_I

  subroutine PSRLineTo_D(Frame,X,Y)
    TYPE(PSFrame)Frame
    Real(Dl) X,Y
    write(frame%fileunit,'(a)') trim(ps_num_string(x))//" ",trim(ps_num_string(Y))//" r"
  END subroutine PSRLineTo_D


  subroutine PSRotate(Frame,A)
    TYPE(PSFrame)Frame
    integer A
    write(frame%fileunit,'(a)') trim(INT2STR(A))//" rotate"
  END subroutine PSRotate

  subroutine PSNewPath(Frame)
    TYPE(PSFrame)Frame
    CALL PsPush(Frame,"newpath")
  END subroutine PSNewPath

  subroutine PSClosePath(Frame)
    TYPE(PSFrame)Frame
    CALL PsPush(Frame,"closepath")
  END subroutine PSClosePath

  subroutine PSStroke(Frame)
    TYPE(PSFrame)Frame
    CALL PsPush(Frame,"stroke")
  END subroutine PSStroke

  subroutine PSFill(Frame)
    TYPE(PSFrame)Frame
    CALL PsPush(Frame,"fill")
  END subroutine PSFill

  subroutine PSSetLinewidth_I(Frame,WIDTH)
    TYPE(PSFrame)Frame
    integer,optional::WIDTH
    If(Present(Width))Then
       Frame%CS%LineWidth=Width
    endif
    write(frame%fileunit,'(a)')trim(Num2Str(Frame%CS%LineWidth))//" setlinewidth"

  END subroutine PSSetLinewidth_I

 

  subroutine PSSetLineWidth_D(Frame,WIDTH)
    TYPE(PSFrame)Frame
    real(dl) WIDTH
    WRITE(Frame%FILEUNIT,"(F9.4,A)")WIDTH," setlinewidth"
  END subroutine PSSetLineWidth_D

  subroutine PSSetLineJoin(Frame,N)
    !!N=0 Miter joins, 1 round joins, 2 Bevel joins
    TYPE(PSFrame)Frame
    integer N
    WRITE(Frame%FILEUNIT,"(I1,A)")N," setlinejoin"
  END subroutine PSSetLineJoin

  subroutine PSSetLineCap(Frame,N)
    !!N=0 butt caps, 1 round caps, 2 projecting caps
    TYPE(PSFrame)Frame
    integer N
    WRITE(Frame%FILEUNIT,"(I1,A)")N," setlinecap"
  END subroutine PSSetLineCap

  subroutine PSgSave(Frame)
    TYPE(PSFrame)Frame
    CALL PsPush(Frame,"gsave")
  END subroutine PSgSave

  subroutine PSgRestore(Frame)
    TYPE(PSFrame)Frame
    CALL PsPush(Frame,"grestore")
  END subroutine PSgRestore

  subroutine PSOutlineBox(Frame,Box)
    TYPE(PSFrame)Frame
    TYPE(PSBox)Box
    CALL PsNewPath(Frame)
    CALL PsMoveto(Frame,Box%XL,Box%YL)
    CALL PsLineto(Frame,Box%XR,Box%YL)
    CALL PsLineto(Frame,Box%XR,Box%YR)
    CALL PsLineto(Frame,Box%XL,Box%YR)
    CALL PSClosePath(Frame)
  END subroutine PSOutlineBox

  subroutine PSScale_I(Frame,XSCALE,YSCALE)
    TYPE(PSFrame)Frame
    integer XSCALE,YSCALE
    write(frame%fileunit,'(a)')trim(INT2STR(XSCALE))//" "//trim(INT2STR(YSCALE))//" scale"
    Frame%CS%xamp=xscale
    Frame%CS%yamp=yscale
  END subroutine PSScale_I



  subroutine PSScale_D(Frame,XSCALE,YSCALE)
    TYPE(PSFrame)Frame
    real(dl) XSCALE,YSCALE
    WRITE(Frame%FILEUNIT,"(a)")Trim(Num2Str(XSCALE,'(F9.3)'))//" "//trim(Num2Str(YSCALE,'(F9.3)'))//" scale"
    Frame%CS%xamp=xscale
    Frame%CS%yamp=yscale
  END subroutine PSScale_D

  subroutine PSDrawBox(Frame,BOX)
    TYPE(PSFrame)Frame
    TYPE(PSBOX)BOX
    CALL PSGSAVE(Frame)
    CALL PsNewPath(Frame)
    CALL PSOutlineBox(Frame,BOX)
    CALL PSSTROKE(Frame)
    CALL PSGRESTORE(Frame)
  END subroutine PSDrawBox

  subroutine PSDrawSolidBox(Frame,BOX)
    TYPE(PSFrame)Frame
    TYPE(PSBOX)BOX
    CALL PSGSAVE(Frame)
    CALL PsNewPath(Frame)
    CALL PSOUTLINEBOX(Frame,BOX)
    CALL PSFILL(Frame)
    CALL PSGRESTORE(Frame)
  END subroutine PSDrawSolidBox

  Subroutine PSColorTex(Frame,Box,text,offset,bgcolor)
    Type(PSframe)Frame
    Type(PSBox) Box
    UNKNOWN_STRING   Text,bgcolor
    Type(PSPoint) Offset
    Call PSGSave(Frame)
    call pssetrgbcolor(Frame,bgcolor)
    call PSDrawSolidBox(Frame,Box)
    call psgrestore(Frame)
    call psgsave(Frame)
    call psmoveto(Frame,Box%xl+offset%x,Box%Yl+offset%y)
    call psprtTex(Frame,Text)
    call psgrestore(frame)
  End Subroutine PSColorTex

  Subroutine PSShiftBox(Box,xshift,yshift)
    Type(PSBox) Box
    Integer(IB) xshift,yshift
    Box%xl=Box%xl+xshift
    Box%xr=Box%xr+xshift
    Box%yl=Box%yl+yshift
    Box%yr=Box%yr+yshift
  End Subroutine PSShiftBox

  Subroutine PSRescaleBox(Box1,Box2,xfac,yfac)
    Type(PSBox) Box1,Box2
    real(dl) xfac
    real(dl),optional::yfac
    Box2%xl=Box1%xl
    Box2%yl=Box1%yl
    Box2%xr=nint((Box1%xr-Box1%xl)*xfac)+Box2%xl
    if(present(yfac))then
       Box2%yr=nint((Box1%yr-Box1%yl)*yfac)+Box2%yl
    else
       Box2%yr=nint((Box1%yr-Box1%yl)*xfac)+Box2%yl
    endif
  End Subroutine PSRescaleBox


  subroutine PSArc(Frame,X,Y,RADIUS,ASTART,AEND)
    TYPE(PSFrame)Frame
    integer X,Y,RADIUS,ASTART,AEND
    write(frame%fileunit,'(a)')trim(INT2STR(X))//" "//trim(INT2STR(Y))//" "  &
         //trim(INT2STR(RADIUS))//" "//trim(INT2STR(ASTART))//" "//trim(INT2STR(AEND))//" arc"
  END subroutine PSArc

  subroutine PSDrawArc(Frame,X,Y,RADIUS,ASTART,AEND)
    TYPE(PSFrame)Frame
    integer X,Y,RADIUS,ASTART,AEND
    CALL PSGSAVE(Frame)
    CALL PsNewPath(Frame)
    CALL PSARC(Frame,X,Y,RADIUS,ASTART,AEND)
    CALL PSSTROKE(Frame)
    CALL PSGRESTORE(Frame)
  END subroutine PSDrawArc

  subroutine PSDrawSolidArc(Frame,X,Y,RADIUS,ASTART,AEND)
    TYPE(PSFrame)Frame
    integer X,Y,RADIUS,ASTART,AEND
    CALL PSGSAVE(Frame)
    CALL PsNewPath(Frame)
    CALL PsMoveTo(Frame,X,Y)
    CALL PSARC(Frame,X,Y,RADIUS,ASTART,AEND)
    CALL PSCLOSEPATH(Frame)
    CALL PSFILL(Frame)
    CALL PSGRESTORE(Frame)
  END subroutine PSDrawSolidArc

  subroutine PSStartProc(Frame)
    TYPE(PSFrame)Frame
    CALL PsPush(Frame,"{")
  END subroutine PSStartProc

  subroutine PSEndProc(Frame)
    TYPE(PSFrame)Frame
    CALL PsPush(Frame,"}")
  END subroutine PSEndProc

  subroutine PSRepeatProc(Frame,N)
    TYPE(PSFrame)Frame
    integer N
    write(frame%fileunit,'(a)')trim(INT2STR(N))//" exch repeat"
  END subroutine PSRepeatProc

  subroutine PSCurrentPoint(Frame)
    TYPE(PSFrame)Frame
    CALL PsPush(Frame,"currentpoint")
  END subroutine PSCurrentPoint

  subroutine PSCLIPCOOR(Frame)
    TYPE(PSFRAME)FRAME
    CALL PsNewPath(FRAME)
    CALL PSOutLineBox(FRAME,Frame%Coor%Bound)
    CALL PsPush(FRAME,"clip")
  END subroutine PSCLIPCOOR


  subroutine PsClipBox(Frame,xmin,ymin,xmax,ymax)
    TYPE(PSFRAME)FRAME
    real(dl) xmin,xmax,ymin,ymax
    CALL PsNewPath(FRAME)
    CALL PScmoveto(FRAME,xmin,ymin)
    call psclineto(frame,xmax,ymin)
    call psclineto(frame,xmax,ymax)
    call psclineto(frame,xmin,ymax)
    call psclosepath(frame)
    CALL PsPush(FRAME,"clip")
  END subroutine PSCLIPBOX
  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !! converting cooridinates
  Subroutine PSCoorTrans_DD(Frame,X,Y,XX,YY)
    Real(Dl) X,y
    Real(dl),optional::xx,yy
    type(PSFrame)Frame
    if(present(xx).and.present(yy))then
       xx=(X-Frame%Coor%OX)/Frame%Coor%XSCALE+Frame%Coor%Bound%XL
       yy=(Y-Frame%Coor%OY)/Frame%Coor%YSCALE+Frame%Coor%Bound%YL
    else
       x=(X-Frame%Coor%OX)/Frame%Coor%XSCALE+Frame%Coor%Bound%XL
       y=(Y-Frame%Coor%OY)/Frame%Coor%YSCALE+Frame%Coor%Bound%YL
    endif
    
  End subroutine PSCoorTrans_DD

  Subroutine PSCoorTrans_pt(Frame,pt)
    Type(Psframe)Frame
    Type(PSPoint)pt
    pt%x=(pt%X-Frame%Coor%OX)/Frame%Coor%XSCALE+Frame%Coor%Bound%XL
    pt%y=(pt%Y-Frame%Coor%OY)/Frame%Coor%YSCALE+Frame%Coor%Bound%YL
  End Subroutine PSCoorTrans_pt


  subroutine PSCoorTrans_D(Frame,X,Y,IX,IY)
    real(dl) X,Y
    integer IX,IY
    TYPE(PSFrame)Frame
    IX=NINT((X-Frame%Coor%OX)/Frame%Coor%XSCALE)+Frame%Coor%Bound%XL
    IY=NINT((Y-Frame%Coor%OY)/Frame%Coor%YSCALE)+Frame%Coor%Bound%YL
  END subroutine PSCoorTrans_D


  subroutine PSGetCoor(CX,CXTYPE,INDX,FORM,XMIN,XMAX,LLOG,Ticks)
    real(dl) CX(PS_N_COORLABELS)
    integer CXTYPE(PS_N_COORLABELS)
    real(dl) XSTART,XEND,TMP
    real(dl) XMIN,XMAX
    real(dl) LOGARR(8)
    integer INDX
    integer Ticks
    real(dl) MaxTicks,MinTicks
    logical LLOG,CHANGETYPE
    SHORT_STRING FORM
    integer N,IXMIN,IXMAX,I,J
    DATA LOGARR /0.30103,0.47712,0.60206,0.69897,0.77815,0.84510,0.90309,0.95424/
    IF(XMAX.LE.XMIN)then
       PRINT*,XMIN,XMAX
       STOP "Error in PSGetCoor: X_max must be grater than X_min."
    endif
    CHANGETYPE=ALL(CXTYPE.EQ.0)
    INDX=0
    MaxTicks=DBLE(Ticks)
    MinTicks=MaxTicks/10._dl
    IF(.not.LLOG)then
       TMP=XMAX-XMIN
       DO WHILE(TMP.Ge.MaxTicks)
          TMP=TMP/10._dl
          INDX=INDX+1
       enddo
       DO WHILE(TMP.Lt.MinTicks)
          TMP=TMP*10._dl
          INDX=INDX-1
       enddo
    endif
    XSTART=XMIN*10.**(-INDX)
    XEND=XMAX*10.**(-INDX)
    IF(XEND.LE.XSTART .OR. ABS(XEND).GT.1.D9)then
       PRINT*,XMIN,XMAX
       STOP "Error in PSGetCoor: X_max must be grater than X_min."
    endif
    IXMIN=GoodFloor(XSTART)
    IXMAX=GoodCeiling(XEND)
    IF(LLOG .OR. XSTART.LT.IXMIN+0.5)then
       XSTART=IXMIN
    ELSE
       IXMIN=IXMIN+1
    endif
    IF(LLOG .OR. XEND.GT.IXMAX-0.5) then
       XEND=IXMAX
    ELSE
       IXMAX=IXMAX-1
    endif
    N=IXMAX-IXMIN+1
    IF(LLOG)then
       DO I=0,N-2
          CX(9*I+1)=IXMIN+I
          IF(CHANGETYPE)CXTYPE(9*I+1)=3
          DO J=1,8
             CX(9*I+1+J)=CX(9*I+1)+LOGARR(J)
          enddo
          IF(CHANGETYPE)CXTYPE(9*I+2:9*I+9)=4
       enddo
       CX(9*N-8)=IXMAX
       IF(CHANGETYPE)CXTYPE(9*N-8)=3
       RETURN
    endif
    IF(N.LE.1)then
       CX(1)=IXMIN-0.5D0
       CX(2)=IXMIN
       CX(3)=IXMIN+0.5D0
       IF(CHANGETYPE)CXTYPE(1:3)=2
       if(trim(form).eq."")FORM="f"
       RETURN
    endif
    IF(N.EQ.4)then
       DO I=0,N-1
          CX(I+1)=IXMIN+I
       enddo
       IF(CHANGETYPE)CXTYPE(1:N)=2
       RETURN
    endif
    IF(N.GE.5 .AND. N.LT.11)then
       DO I=0,N-1
          CX(I+1)=IXMIN+I
          IF(CHANGETYPE)then
             IF(MOD(NINT(CX(I+1)),2).EQ.0)then
                CXTYPE(I+1)=2
             ELSE
                CXTYPE(I+1)=4
             endif
          endif
       enddo
       RETURN
    endif
    IF(N.GE.11.AND.N.LE.20)then
       DO I=0,N-1
          CX(I+1)=IXMIN+I
          IF(CHANGETYPE)then
             IF(MOD(NINT(CX(I+1)),5).EQ.0)then
                CXTYPE(I+1)=2
             ELSE
                CXTYPE(I+1)=4
             endif
          endif
       enddo
       RETURN
    endif
    IF(N.GE.21)then
       DO I=0,N-1
          CX(I+1)=IXMIN+I
          IF(CHANGETYPE)then
             IF(MOD(NINT(CX(I+1)),10).EQ.0)then
                CXTYPE(I+1)=2
             ELSEIF(MOD(NINT(CX(I+1)),2).EQ.0)then
                CXTYPE(I+1)=4
             ELSE
                CXTYPE(I+1)=7
             endif
             CX(I+1)=CX(I+1)/10.
          endif
       enddo
       IF(CHANGETYPE)INDX=INDX+1
       RETURN
    endif
    I=1
    IF(IXMIN-XSTART.GE.0.35D0)then
       CX(I)=IXMIN-0.5D0
       IF(CHANGETYPE)then
          IF(N.EQ.2)then
             CXTYPE(I)=2
          ELSE
             CXTYPE(I)=4
          endif
       endif
       I=I+1
    ELSE
       CX(I)=IXMIN
       IF(CHANGETYPE)CXTYPE(I)=2
       I=I+1
    endif
    DO WHILE(CX(I-1)+0.35 .LT. XEND)
       CX(I)=CX(I-1)+0.5D0
       IF(CHANGETYPE)then
          IF(N.EQ.2)then
             CXTYPE(I)=2
          ELSE
             IF(CXTYPE(I-1).EQ.4)then
                CXTYPE(I)=2
             ELSE
                CXTYPE(I)=4
             endif
          endif
       endif
       I=I+1
    enddo
    IF(N.EQ.2 .and. trim(form).eq."") FORM="f"
    IF(N.LE.0) STOP "Unknown Error in PSGetCoor."
  END subroutine PSGetCoor

  !!string library
  function rgbstr(R)
    real(dl),INTENT(IN)::R
    SHORT_STRING rgbstr
    integer K
    K=NINT(1000*R)
    IF(K.EQ.0)then
       rgbstr="0"
    ELSEIF(K.LT.10)then
       rgbstr=".00"//trim(INT2STR(K))
    ELSEIF(K.LT.100)then
       rgbstr=".0"//trim(INT2STR(K))
    ELSEIF(K.LT.1000)then
       rgbstr="."//trim(INT2STR(K))
    ELSE
       rgbstr="1"
    endif
  END function rgbstr



  function prtrgb(R,G,B)
    real(dl),INTENT(IN)::R,G,B
    SHORT_STRING prtrgb
    prtrgb = " " // trim(rgbstr(R)) // " " // trim(rgbstr(G)) // " "// trim(rgbstr(B))
  END function prtrgb

  function GoodFloor(X)
    real(dl) X
    integer GoodFloor
    GoodFloor=FLOOR(X)
    IF(X-GoodFloor-1.GT.-1.D-6)then
       GoodFloor=GoodFloor+1
    endif
  END function GoodFloor

  function GoodCeiling(X)
    real(dl) X
    integer GoodCeiling
    GoodCeiling=CEILING(X)
    IF(X-GoodCeiling+1.LT.1.D-6)then
       GoodCeiling=GoodCeiling-1
    endif
  END function GoodCeiling

  function AnyInsideTriangle(X1,Y1,X2,Y2,X3,Y3,X,Y,N)
    logical AnyInsideTriangle
    integer N,I
    integer X1,Y1,X2,Y2,X3,Y3,X(N),Y(N)
    integer C1(3),C2(3),C3(3)
    C1(1)=Y2-Y3
    C1(2)=X3-X2
    C1(3)=X2*Y3-Y2*X3
    IF((Y1-Y2)*(X3-X2)-(Y3-Y2)*(X1-X2).LT.0)C1=-C1
    C2(1)=Y3-Y1
    C2(2)=X1-X3
    C2(3)=X3*Y1-Y3*X1
    IF((Y2-Y3)*(X1-X3)-(Y1-Y3)*(X2-X3).LT.0)C2=-C2
    C3(1)=Y1-Y2
    C3(2)=X2-X1
    C3(3)=X1*Y2-Y1*X2
    IF((Y3-Y1)*(X2-X1)-(Y2-Y1)*(X3-X1).LT.0)C3=-C3
    DO I=1,N
       IF(C1(1)*X(I)+C1(2)*Y(I)+C1(3) .LT. 0)CYCLE
       IF(C2(1)*X(I)+C2(2)*Y(I)+C2(3) .LT. 0)CYCLE
       IF(C3(1)*X(I)+C3(2)*Y(I)+C3(3) .LT. 0)CYCLE
       AnyInsideTriangle=.TRUE.
       RETURN
    enddo
    AnyInsideTriangle=.FALSE.
  END function AnyInsideTriangle

  function InsideTriangle(X1,Y1,X2,Y2,X3,Y3,X,Y)
    logical InsideTriangle
    integer X1,Y1,X2,Y2,X3,Y3,X,Y
    InsideTriangle=SAME_SIGN((Y3-Y1)*(X2-X1)-(Y2-Y1)*(X3-X1),(Y-Y1)*(X2-X1)-(Y2-Y1)*(X-X1)) .AND. &
         SAME_SIGN((Y1-Y2)*(X3-X2)-(Y3-Y2)*(X1-X2),(Y-Y2)*(X3-X2)-(Y3-Y2)*(X-X2)) .AND. &
         SAME_SIGN((Y2-Y3)*(X1-X3)-(Y1-Y3)*(X2-X3),(Y-Y3)*(X1-X3)-(Y1-Y3)*(X-X3))
  END function InsideTriangle

  function same_sign(X,Y)
    integer X,Y
    logical same_sign
    same_sign=(X.GE.0 .AND. Y.GE.0) .OR. (X.LE.0 .AND. Y.LE.0)
  END function same_sign

  function TwiceArea(X1,Y1,X2,Y2,X3,Y3)
    !! positive iff 1,2,3 conterclockwise
    integer TwiceArea
    integer X1,Y1,X2,Y2,X3,Y3
    TwiceArea=X1*(Y2-Y3)+X2*(Y3-Y1)+X3*(Y1-Y2)
  END function TwiceArea

  function NextBracket(CHR)
    UNKNOWN_STRING  CHR
    integer I,STRLEN,COUNT
    integer NextBracket
    logical IGNORE
    STRLEN=len(CHR)
    COUNT=0
    I=1
    IGNORE=.FALSE.
    DO I=1,STRLEN
       IF(CHR(I:I) .eq. Const_backslash)then
          IGNORE=.True.
          CYCLE
       ELSE
          IGNORE=.False.
       endif
       IF(.not.IGNORE)then
          IF(CHR(I:I).EQ."{")then
             COUNT=COUNT+1
          ELSEIF(CHR(I:I).EQ."}")then
             COUNT=COUNT-1
             IF(COUNT.LT.0)then
                NextBracket=I
                RETURN
             endif
          endif
       endif
    enddo
    NextBracket=0
  END function NextBracket

  subroutine PSColor2RGB(ColorName,R,G,B)
    UNKNOWN_STRING  ColorName
    REAL(dl) R,G,B
    integer I
    SHORT_STRING Color
    if(Trim(ColorName).eq.'')then
       R=0.
       G=0.
       B=0.
       return
    endif
    Color=trim(COLORNAME)
    IF(color(1:4).EQ."rgb[" .OR. Color(1:4).EQ."RGB[")then
       I=SCAN(COLOR,"]")
       READ(COLOR(5:I-1),*) R,G,B
       RETURN
    endif
    COLOR=""
    DO I=1,MIN(LEN_trim(ColorName),20)
       IF(ICHAR(ColorName(I:I)).GE.65 .AND. ICHAR(ColorName(I:I)).LE.90 )then
          Color=Trim(Color)//CHAR(ICHAR(ColorName(I:I))+32)
       ELSEIF(ICHAR(ColorName(I:I)).GE.97 .AND. ICHAR(ColorName(I:I)).LE.122)then
          Color=Trim(Color)//ColorName(I:I)
       endif
    enddo
    Select case(trim(COLOR))
    case("")
       R=0.
       G=0.
       B=0.
       RETURN
    case("red")
       R=1.
       G=0.
       B=0.
       RETURN
    case("green")
       R=0.
       G=1.
       B=0.
       RETURN
    case("blue")
       R=0.
       G=0.
       B=1.
       RETURN
    case("yellow")
       R=1.
       G=1.
       B=0.
       RETURN
    case("magenta")
       R=1.
       G=0.
       B=1.
       RETURN
    case("cyan","turquoise")
       R=0.
       G=1.
       B=1.
       RETURN
    case("black")
       R=0.
       G=0.
       B=0.
       RETURN
    case("white")
       R=1.
       G=1.
       B=1.
       RETURN
    case("gray","grey")
       R=0.5
       G=0.5
       B=0.5
       RETURN
    case("orange")
       R=1.
       G=0.5
       B=0.2
       RETURN
    case("violet")
       R=0.55
       G=0.22
       B=0.79
       RETURN
    case("brown")
       R=0.5
       G=0.25
       B=0.
       RETURN
    case("pink")
       R=0.98
       G=0.69
       B=0.73
       RETURN
    case("gold")
       R=0.83
       G=0.63
       B=0.09
       RETURN
    case("purple")
       R=0.56
       G=0.21
       B=0.94
       RETURN
    case("maroon")
       R=0.51
       G=0.02
       B=0.25
       RETURN
    case("slateblue")
       R=0.21
       G=0.45
       B=0.78
       RETURN
    case("skyblue")
       R=0.24
       G=0.6
       B=1.
       RETURN
    case("lawngreen","grassgreen")
       R=0.53
       G=0.97
       B=0.10
       RETURN
    case("darkgray","darkgrey")
       R=0.2
       G=0.2
       B=0.2
       RETURN
    case("lightgray","lightgrey")
       R=0.8
       G=0.8
       B=0.8
       RETURN
    case("darkblue")
       R=0.
       G=0.
       B=0.5
       RETURN
    case("darkgreen")
       R=0.
       G=0.5
       B=0.
       RETURN
    case("darkred")
       R=0.5
       G=0.
       B=0.
       RETURN
    case DEFAULT
       Print*, "Unknown color:",trim(ColorName)
       R=0.
       G=0.
       B=0.
    END select
  END subroutine PSColor2RGB

  Subroutine PSDelay(Frame,n)
    Type(PSFrame) Frame
    integer(IB) n
    call pscurrentpoint(Frame)
    call pspush(Frame,Trim(Num2Str(n*1000))//" {0 0 m} repeat m")
  End Subroutine PSDelay

  Subroutine PSUpdateTex(Frame,rx,ry,width,height,Str)
    Type(PSFrame)Frame
    UNKNOWN_STRING   Str
    real(dl) rx,ry
    integer(IB) width,height
    call psgsave(Frame)
    call psnewpath(Frame)
    call PSMoveToRatio(Frame,rx,ry)
    call psrlineto(Frame,width,0)
    call psrlineto(Frame,0,height)
    call psrlineto(Frame,-width,0)
    call psclosePath(Frame)
    call Pssetgray(frame,1._dl)
    call psfill(Frame)
    call pssetgray(frame,0._dl)
    call psmovetoratio(Frame,rx,ry)
    call psrmoveto(Frame,2,2)
    call psprttex(Frame,str)
    call psgrestore(Frame)
  End Subroutine PSUpdateTex


  Subroutine PSPlotContour(Frame,h,nx,ny,xmin,xmax,ymin,ymax,hcut)
    real(dl),parameter::epsilon=1.d-20  !!make sure it is small enough
    Integer(IB),parameter::Nmid=5
    Type(Psframe) frame
    Type(PSGridPoint) gp(4096)
    Integer(IB),dimension(:),allocatable::istart
    integer(IB) nx,ny
    real(dl) h(nx,ny)
    real(dl) xmin,xmax,ymin,ymax,hcut,hmax,hmin,dh
    real(dl)dx,dy
    integer(IB) ns,nump
    real(dl) xmid(Nmid),ymid(Nmid),s(Nmid)
    integer(IB) i,j,k,ib,it,im
    Type(PSGridPoint) gNext
    type(PSPath) Path
    real(dl) accuracy

    accuracy = min( abs(Frame%Coor%XUR-Frame%Coor%OX), abs(Frame%Coor%YUR-Frame%Coor%OY) )/10.**(ps_resolution_level +3)
    
    hmax=maxval(h)
    hmin=minval(h)
    dh=(hmax-hmin)*epsilon
    if(hcut.ge.hmax.or. hcut.le.hmin)return  !!nothing to plot
    do i=1,nx
       do j=1,ny
          if(abs(h(i,j)-hcut).lt.dh)h(i,j)=hcut+sign(dh,h(i,j)-hcut)  !!make sure no grid precisely hit hcut
       enddo
    enddo
    allocate(istart(nx))
    dx=(xmax-xmin)/(nx-1)
    dy=(ymax-ymin)/(ny-1)
    istart(1)=1
    nump=0

#define HERE h(i,j)
#define RIGHT h(i+1,j)
#define UP h(i,j+1)
#define LEFT h(i-1,j)
#define UPRIGHT h(i+1,j+1)
#define UPLEFT h(i-1,j+1)
#define LOWERLEFT h(i-1,j-1)
#define LOWER h(i,j-1)
    do i=2,nx-1
       istart(i)=nump+1
       do j=1,ny-1
          if(HERE.gt.hcut)then
             if(hcut.gt.RIGHT)then
                nump=nump+1
                gp(nump)%x=i
                gp(nump)%y=j
                gp(nump)%direction=2
                gp(nump)%lambda=(HERE-hcut)/(HERE-RIGHT)
                gp(nump)%a=-(RIGHT+UP-HERE-UPRIGHT)/(HERE-RIGHT)
                gp(nump)%b=(UP-HERE)/(HERE-hcut)+gp(nump)%a
             endif
             if(hcut.gt. UP)then
                nump=nump+1
                gp(nump)%x=i
                gp(nump)%y=j
                gp(nump)%direction=3
                gp(nump)%lambda=(HERE-hcut)/(HERE-UP)
                gp(nump)%a=-(LEFT+UP-HERE-UPLEFT)/(HERE-UP)
                gp(nump)%b=(LEFT-HERE)/(HERE-hcut)+gp(nump)%a
             endif
             if(hcut.gt.LEFT)then
                nump=nump+1
                gp(nump)%x=i
                gp(nump)%y=j
                gp(nump)%direction=4
                gp(nump)%lambda=(HERE-hcut)/(HERE-LEFT)
                gp(nump)%a=-(LEFT+LOWER-LOWERLEFT-HERE)/(HERE-LEFT)
                gp(nump)%b=(LOWER-HERE)/(HERE-hcut)+gp(nump)%a
             endif
          elseif(hcut.lt.UP)then
             nump=nump+1
             gp(nump)%x=i
             gp(nump)%y=j+1
             gp(nump)%direction=1
             gp(nump)%lambda=(UP-hcut)/(UP-HERE)
             gp(nump)%a=-(UPRIGHT+HERE-UP-RIGHT)/(UP-HERE)
             gp(nump)%b=(UPRIGHT-UP)/(UP-hcut)+gp(nump)%a
          endif
       enddo
    enddo

#undef HERE
#undef RIGHT
#undef UP
#undef LEFT
#undef UPRIGHT
#undef UPLEFT
#undef LOWER
#undef LOWERLEFT
    istart(nx)=nump+1

    gp(1:nump)%used=.false.

    ns=0
!    do j=1,nump
!       write(*,"(4I5)") j,gp(j)%x,gp(j)%y,gp(j)%direction
!    enddo
    k=1
   

    do while(ns.lt.nump) !!not done ,need to find a new contour
       do while(gp(k)%used)
          k=k+1
       !   if(k.gt.nump) stop "k >nump"
       !   print*,k
       enddo
       !!now construct a contour from k
       call ClearPath(path)
       i=k
       xmid = 1.d0
       ymid = 1.d0
       do
          if(gp(i)%b.gt.0._dl)then
             if(gp(i)%a + gp(i)%b/(1._dl/gp(i)%lambda-1._dl).ge.1._dl )then
                call findgen(s,0._dl,1._dl/(gp(i)%a + gp(i)%b/(1._dl/gp(i)%lambda-1._dl)))
                select case(gp(i)%direction)
                case(1)
                   gnext%x=gp(i)%x+1
                   gnext%y=gp(i)%y-1
                   gnext%direction=4
                case(2)
                   gnext%x=gp(i)%x+1
                   gnext%y=gp(i)%y+1
                   gnext%direction=1
                case(3)
                   gnext%x=gp(i)%x-1
                   gnext%y=gp(i)%y+1
                   gnext%direction=2
                case(4)
                   gnext%x=gp(i)%x-1
                   gnext%y=gp(i)%y-1
                   gnext%direction=3
                end select
             else
                call findgen(s,0._dl,1._dl)
                select case(gp(i)%direction)
                case(1)
                   gnext%x=gp(i)%x+1
                   gnext%y=gp(i)%y
                case(2)
                   gnext%x=gp(i)%x
                   gnext%y=gp(i)%y+1
                case(3)
                   gnext%x=gp(i)%x-1
                   gnext%y=gp(i)%y
                case(4)
                   gnext%x=gp(i)%x
                   gnext%y=gp(i)%y-1
                end select    
                gnext%direction=gp(i)%direction            
             endif
          else
             if(gp(i)%a-gp(i)%b.ge.1._dl)then
                call findgen(s,0._dl,1._dl/(gp(i)%a-gp(i)%b))
                gnext%x=gp(i)%x
                gnext%y=gp(i)%y
                gnext%direction=mod(gp(i)%direction,4)+1
             else
                call findgen(s,0._dl,1._dl)
                select case(gp(i)%direction)
                case(1)
                   gnext%x=gp(i)%x+1
                   gnext%y=gp(i)%y
                case(2)
                   gnext%x=gp(i)%x
                   gnext%y=gp(i)%y+1
                case(3)
                   gnext%x=gp(i)%x-1
                   gnext%y=gp(i)%y
                case(4)
                   gnext%x=gp(i)%x
                   gnext%y=gp(i)%y-1
                end select
                gnext%direction=gp(i)%direction
             endif
          endif
          select case(gp(i)%direction)
          case(1)
             xmid(1)=gp(i)%x
             ymid(1)=gp(i)%y-gp(i)%lambda
             do j=2,nmid-1
                xmid(j)=gp(i)%x+s(j)
                ymid(j)=gp(i)%y-gp(i)%lambda*(1._dl+gp(i)%b/(1._dl/s(j) -gp(i)%a))
             enddo
          case(2)
             xmid(1)=gp(i)%x+gp(i)%lambda
             ymid(1)=gp(i)%y
             do j=2,nmid-1
                xmid(j)=gp(i)%x+gp(i)%lambda*(1._dl+gp(i)%b/(1._dl/s(j) -gp(i)%a))
                ymid(j)=gp(i)%y+s(j)
             enddo
          case(3)
             xmid(1)=gp(i)%x
             ymid(1)=gp(i)%y+gp(i)%lambda
             do j=2,nmid-1
                xmid(j)=gp(i)%x-s(j)
                ymid(j)=gp(i)%y+gp(i)%lambda*(1._dl+gp(i)%b/(1._dl/s(j) -gp(i)%a))
             enddo
          case(4)
             xmid(1)=gp(i)%x-gp(i)%lambda
             ymid(1)=gp(i)%y
             do j=2,nmid-1
                xmid(j)=gp(i)%x-gp(i)%lambda*(1._dl+gp(i)%b/(1._dl/s(j) -gp(i)%a))
                ymid(j)=gp(i)%y-s(j)
             enddo
          end select


          ib=istart(gnext%x)
          it=istart(gnext%x+1)

          do while(it.gt.ib+1)
             im=(it+ib)/2
             if(gp(im)%y.gt.gnext%y .or. (gp(im)%y.eq.gnext%y .and. gp(im)%direction.gt.gnext%direction))then

                it=im
             else
                ib=im
             endif
          enddo
         if(gp(ib)%x.ne.gnext%x .or. gp(ib)%y.ne.gnext%y .or. gp(ib)%direction .ne. gnext%direction) then
            print*,gnext%x,gnext%y,gnext%direction
            print*,gp(ib)%x,gp(ib)%y,gp(ib)%direction
            stop "Unknown error in Psplotcontour"
         endif
         do j=1,nmid-1
            call Add2path(path,xmin+dx*(xmid(j)-1._dl),ymin+dy*(ymid(j)-1._dl),accuracy)
         enddo
         gp(i)%used=.true. !!clear i
         ns=ns+1  !!count cleared points
         i=ib !!next point
         if(gp(ib)%used)exit
      enddo
       if(path%N.gt.0)then
          if(Frame%cs%twodsmooth.gt.1 .and. path%n.gt.Frame%cs%twodsmooth)then
             call PSRunBsplinePath(Frame,path,Frame%cs%twodsmooth)
          else
             call PsRunPath(Frame,path,.true.)
          endif
       Endif
    enddo
    deallocate(istart)
  
  End Subroutine PSPlotContour


  Subroutine ClearPath(path)
    Type(PSPath) Path
    Path%N=0
  End Subroutine ClearPath


  subroutine PsRunPath(Frame,PATH,closepath)
    TYPE(PSFrame)Frame
    TYPE(PSPATH)PATH
    real(dl) xx,yy
    integer I,ns,ix,iy,ix2,iy2
    logical,optional::closepath
    STRING str,laststr
    if(path%N.le.0) return
    CALL pscmoveto(Frame,PATH%X(1),PATH%Y(1))
    ns=0
    if(ps_resolution_level.gt.1)then
       laststr=''
       DO I=PATH%N,2,-1
          call psCoorTrans(Frame,path%x(i),path%y(i),xx,yy)
          str=trim(PS_Num_STRing(xx))//" "//trim(PS_Num_STRing(yy))
          if(trim(str).ne.trim(laststr))then
             write(frame%fileunit,'(a)') trim(str)
             laststr=trim(str)
             ns=ns+1
          endif
       enddo
       write(frame%fileunit,'(a)')trim(INT2STR(ns))//" {l} repeat"
    else
       call pscoortrans(frame,path%x(1),path%y(1),ix,iy)
       write(frame%fileunit,'(a)') trim(Int2Str(ix))//' '//trim(Int2Str(iy))//' m'
       call pscoortrans(frame,path%x(path%N),path%y(path%N),ix,iy)
       do i=path%N-1,1,-1
          call pscoortrans(frame,path%x(i),path%y(i),ix2,iy2)
          if(ix2.ne.ix .or. iy2.ne.iy)then
             write(frame%fileunit,'(a)') Trim(PSAbbrev(trim(Int2Str(ix-ix2))//' '//trim(Int2Str(iy-iy2))))
             ix=ix2
             iy=iy2
             ns=ns+1
          endif
       enddo
       write(frame%fileunit,'(a)') trim(int2str(ns))//' {r} repeat'
    endif
    if(present(closepath))then
       if(closepath)call psclineto(Frame,Path%x(1),path%y(1))
    endif
  END subroutine PsRunPath

  Subroutine PSRunBSplinePath(Frame,path,bspline_order)
    Integer(IB),parameter::Nb=11
    Type(PSFrame)Frame
    Type(PSPath)Path
    Integer(IB),intent(in)::bspline_order
    Integer(IB) i,j,ns,ix,iy,ix2,iy2
    real(dl) t(nb)
    Real(dl) xx,yy
    STRING str,laststr
    real(dl) accuracy
    accuracy = min( abs(Frame%Coor%XUR-Frame%Coor%OX), abs(Frame%Coor%YUR-Frame%Coor%OY))/10.**(ps_resolution_level +3)
    if(path%N.eq.0) return
    do i=1,bspline_order
       CALL add2path(path,path%x(i),path%y(i),0._dl)   !!wrap data
    enddo
    call findgen(t,0._dl,1._dl)
    i=bspline_order+1
    xx=UBSplineSeg(bspline_order,0._dl,path%x(i-bspline_order:i))
    yy=UBSplineSeg(bspline_order,0._dl,path%y(i-bspline_order:i))
    call pscmoveto(Frame,xx,yy)
    ns=0
    laststr=''
    if(ps_resolution_level.ge.1)then
       do i=path%n, bspline_order+1, -1
          do j=nb,2,-1
             xx=UBSplineSeg(bspline_order,t(j),path%x(i-bspline_order:i))
             yy=UBSplineSeg(bspline_order,t(j),path%y(i-bspline_order:i))
             
             call PSCoorTrans(Frame,xx,yy)
          
             str=trim(PS_Num_STRing(xx))//" "//trim(PS_Num_String(yy))
             if(trim(str).ne.trim(laststr))then
                laststr=trim(str)
                write(frame%fileunit,'(a)')trim(str)
                ns=ns+1
             endif
          enddo
       enddo
       WRITE(Frame%FileUnit,'(a)')trim(Int2Str(ns))//" {l} repeat"
    else
       i=path%n
       xx=UBSplineSeg(bspline_order,1._dl,path%x(i-bspline_order:i))
       yy=UBSplineSeg(bspline_order,1._dl,path%y(i-bspline_order:i))
       call PSCoorTrans(Frame,xx,yy,ix,iy)

       do i=path%n, bspline_order+1, -1
          do j=nb,2,-1
             xx=UBSplineSeg(bspline_order,t(j),path%x(i-bspline_order:i))
             yy=UBSplineSeg(bspline_order,t(j),path%y(i-bspline_order:i))
             call PSCoorTrans(Frame,xx,yy,ix2,iy2)
             if(ix2.ne.ix .or. iy2.ne.iy)then
                write(frame%fileunit,'(a)') Trim(PSAbbrev(trim(Int2Str(ix-ix2))//' '//trim(Int2Str(iy-iy2))))
                ix=ix2
                iy=iy2
                ns=ns+1
             endif
          enddo
       enddo
       WRITE(Frame%FileUnit,'(a)')trim(Int2Str(ns))//" {r} repeat"
    endif
    call psclosepath(frame)
    path%n=path%n-bspline_order
  End Subroutine PSRunBSplinePath


  Function PSInterpPoints(p1,p2,lambda) !!return p1(1-lambda)+p2*lambda
    type(PSPoint) p1,p2
    real(dl) lambda
    Type(PSPoint) PSInterpPoints
    PSInterpPoints%x=p1%x*(1._dl-lambda)+p2%x*lambda
    PSInterpPoints%y=p1%y*(1._dl-lambda)+p2%y*lambda
  End Function PSInterpPoints


  Subroutine Add2Path(path,x,y,accuracy)
    type(PsPath)path
    real(dl) x,y,delta
    real(dl) accuracy
    if(path%n.le.1)then
       path%n=path%n+1
       if(path%n .gt. ps_pathdepth) stop "path over flow"
       path%x(path%n)=x
       path%y(path%n)=y
    else
       if(abs(x-path%x(path%n))+abs(y-path%y(path%n)).lt.accuracy)return
       delta = path%x(path%n-2)*path%y(path%n-1)+path%x(path%n-1)*y+x*path%y(path%n-2) &
            -path%x(path%n-2)*y - path%x(path%n-1)*path%y(path%n-2) -x*path%y(path%n-1)
       if(abs(delta).lt.accuracy**2*1.d-6)then !!a straightline
          path%x(path%n)=x
          path%y(path%n)=y          
       else
          path%n=path%n+1
          if(path%n .gt. ps_pathdepth) stop "path over flow"
          path%x(path%n)=x
          path%y(path%n)=y
       endif
    endif
  end Subroutine Add2Path



  Function PS_Num_String(x)
    SHORT_STRING PS_Num_String
    character head
    Real(dl) x,ax
    Integer(IB) i
    ax=abs(x)
    if(ax.lt.0.)then
       head='-'
    else
       head=''
    endif
    select case(ps_resolution_level)
    case(3)
       PS_Num_String=Trim(Num2Str(x,'(G15.5)'))
       return
    case(2)
       if(ax.ge.1._dl)then
          if(abs(nint(x)-x).le.1.d-4*ax)then
             PS_Num_String=Trim(Int2Str(nint(x))) 
             return
          endif
          if(ax.lt.10.)then  
             PS_Num_String=Trim(Num2Str(x,'(F12.3)'))
             PS_Num_String=PS_Num_String( 1 : verify(PS_Num_String,'0',BACK=.true.))
             return 
          elseif(ax.lt.100.)then
             PS_Num_String=Trim(Num2Str(x,'(F12.2)'))
             PS_Num_String=PS_Num_String( 1 : verify(PS_Num_String,'0',BACK=.true.))             
             return    
          else
             PS_Num_String=Trim(Num2Str(x,'(G15.4)'))
             return
          endif
       endif


       !!less than 1 
       i=nint(ax*1.d5)
       select case(i)
       case(0)
          PS_Num_String='0'
          return
       case(1:9)
          PS_Num_String=head//'.0000'//trim(Int2Str(i))
       case(10:99)
          PS_Num_String=head//'.000'//trim(Int2Str(i))
       case(100:999)
          PS_Num_String=head//'.00'//trim(Int2Str(i))
       case(1000:9999)
          PS_Num_String=head//'.0'//trim(Int2Str(i))
       case(10000:99999)
          PS_Num_String=head//'.'//trim(Int2Str(i))
       case default
          PS_Num_String=head//'1'
          return
       end select
       PS_Num_String=PS_Num_String( 1 : verify(PS_Num_String,'0',BACK=.true.))
       return
    case(1)
       if(ax.ge.1._dl)then

          if(abs(nint(x)-x).le.0.005*ax)then
             PS_Num_String=Trim(Int2Str(nint(x))) !!well assuming for a plot the 0.5% accuracy is irrelevant
             return
          endif
          if(ax.lt.10.)then  
             PS_Num_String=Trim(Num2Str(x,'(F12.2)'))
             PS_Num_String=PS_Num_String( 1 : verify(PS_Num_String,'0',BACK=.true.))
             return 
          elseif(ax.lt.100.)then
             PS_Num_String=Trim(Num2Str(x,'(F12.1)'))
             PS_Num_String=PS_Num_String( 1 : verify(PS_Num_String,'0',BACK=.true.))
             return    
          else
             PS_Num_String=Trim(Num2Str(x,'(G15.4)'))
             return
          endif
       endif


       !!less than 1 
       i=nint(ax*1.d4)
       select case(i)
       case(0)
          PS_Num_String='0'
          return
       case(1:9)
          PS_Num_String=head//'.000'//trim(Int2Str(i))
       case(10:99)
          PS_Num_String=head//'.00'//trim(Int2Str(i))
       case(100:999)
          PS_Num_String=head//'.0'//trim(Int2Str(i))
       case(1000:9999)
          PS_Num_String=head//'.'//trim(Int2Str(i))
       case default
          PS_Num_String=head//'1'
          return
       end select
       PS_Num_String=PS_Num_String( 1 : verify(PS_Num_String,'0',BACK=.true.))
       return
    case(0)
       PS_Num_String=Trim(Int2Str(nint(x))) 
       return
    case default
       print*,"Unknown resolution level", PS_Resolution_level
       stop
    end select
    
  End Function PS_Num_String

  Function PSAbbrev(str)
    UNKNOWN_STRING  str
    SHORT_STRING PSAbbrev
    select case(trim(str))
    case('-1 -1')
        PSAbbrev='z'
     case('-1 0')
        PSAbbrev='a'
     case('-1 1')
        PSAbbrev='q'
     case('0 -1')
        PSAbbrev='x'
     case('0 1')
        PSAbbrev='w'
     case('1 -1')
        PSAbbrev='c'
     case('1 0')
        PSAbbrev='d'
     case('1 1')
        PSAbbrev='e'
     case default
        PSAbbrev=str
     end select
   end Function PSAbbrev


   Subroutine PSLabelMark(frame,mark)
     type(psframe)frame 
     type(PSmark) mark
     real(dl)::tmp(4)
     integer i
     call psgsave(frame)
     call pssetcolor(frame,mark%color)
     select case(trim(mark%str))
     case("\point","\cross","\triangle","\solidtriangle","\circle","\solidcircle","\box","\solidbox")
        read(mark%pos,*) tmp(1:2)
        if(mark%size.le.0)then
           i=3
        else
           i=nint(mark%size)
        endif
        call psplotpoint(frame,tmp(1),tmp(2),i,trim(mark%str(2:)))
     case("\line")
        read(mark%pos,*) tmp(1:4)
        if(mark%size.gt.0)then
           call pssetlinewidth(frame,mark%size)
        endif
        call psplotline(frame,tmp(1),tmp(2),tmp(3),tmp(4))
     case default
        read(mark%pos,*) tmp(1:2)
        call pscmoveto(frame,tmp(1),tmp(2))
        if(mark%size .gt. 0)then
           call pssetfont(frame,frame%font%t,nint(mark%size))
        endif
        call psprttex(frame,trim(mark%str))
     end select
     call psgrestore(frame)
   End Subroutine PSLabelMark

   subroutine psplotdatafile(frame,filename,color1,color2,color3,color4,color5,color6)
     type(psframe)frame
     UNKNOWN_STRING  filename
     UNKNOWN_STRING  ,optional::color1,color2,color3,color4,color5,color6
     type(file_pointer) fp
     real(dl),dimension(:),allocatable::x,y
     real(dl) lbd,ubd
     integer i,n
     fp = open_file(trim(filename),"r")
     read(fp%unit,*)n,lbd,ubd
     allocate(x(n),y(n))
     call findgen(x,lbd,ubd)
     i=1
     do 
        read(fp%unit,*,END=100) y
        select case(i)
        case(1)
           if(present(color1)) call pssetcolor(frame,color1)
        case(2)
           if(present(color2)) call pssetcolor(frame,color2)
        case(3)
           if(present(color3)) call pssetcolor(frame,color3)
        case(4)
           if(present(color4)) call pssetcolor(frame,color4)
        case(5)
           if(present(color5)) call pssetcolor(frame,color5)
        case(6)
           if(present(color6)) call pssetcolor(frame,color6)
        end select
        call psplotcurve(frame,x,y)
        i=i+1
     enddo
100  continue
     deallocate(x,y)
     call close_file(fp)
   end subroutine psplotdatafile


   subroutine psplottabfunc(id1, id2, id3, form1, form2, form3, xlog, ylog, filename, xlabel, ylabel, xmin, xmax)
     use tabfunc_utils
     integer id1
     integer, optional::id2, id3
     logical,optional::xlog, ylog
     UNKNOWN_STRING ,optional::filename, xlabel, ylabel, form1, form2, form3
     real(dl),optional::xmin, xmax
     real(dl) ymin, ymax, xl, xr
     type(PsFrame)frame
     integer,parameter::N=2048
     real(dl) x1(N), x2(N), x3(N),y1(N), y2(N), y3(N)
     if(present(xmin))then
        xl = xmin
     else
        xl = tabfunc_space(id1)%left
     endif
     if(present(xmax))then
        xr = xmax
     else
        xr = tabfunc_space(id1)%right
     endif
     call findgen(x1, max(xl, tabfunc_space(id1)%left), min(tabfunc_space(id1)%right, xr))
     y1 = tabfunc_eval(id1, x1)
     ymin = minval(y1)
     ymax = maxval(y1)
     if(present(id2))then
        if(.not. present(xmin))then
           xl = min(tabfunc_space(id2)%left, xl)
        endif
        if(.not. present(xmax))then
           xr = max(tabfunc_space(id2)%right, xr)
        endif
        call findgen(x2, max(tabfunc_space(id2)%left, xl), min(tabfunc_space(id2)%right, xr))
        y2 = tabfunc_eval(id2,x2)
        ymin = min(ymin, minval(y2))
        ymax = max(ymax, maxval(y2))
     endif
     if(present(id3))then
        if(.not. present(xmin))then
           xl = min(tabfunc_space(id3)%left, xl)
        endif
        if(.not. present(xmax))then
           xr = max(tabfunc_space(id3)%right, xr)
        endif
        call findgen(x3, max(tabfunc_space(id3)%left,xl), min(tabfunc_space(id3)%right,xr))
        y3 = tabfunc_eval(id3,x3)
        ymin = min(ymin, minval(y3))
        ymax = max(ymax, maxval(y3))
     endif
     call psdefaultframe(frame)
     if(present(xlog)) frame%coor%xlog = xlog
     if(present(ylog))frame%coor%ylog = ylog
     if(present(filename)) frame%filename = trim(filename)//".eps"
     if(present(xlabel)) frame%coor%xlabel = xlabel
     if(present(ylabel)) frame%coor%ylabel = ylabel
     call psstart(frame)
     call psdefaultcoor(frame, xl, xr, ymin, ymax)
     if(present(form1)) call pssetcolor(frame, form1)
     call psplotcurve(frame, x1, y1)
     if(present(id2))then
        if(present(form2))then
           call pssetcolor(frame, form2)
        else
           call psdash(frame)
        endif
        call psplotcurve(frame, x2, y2)
     endif
     if(present(id3))then
        if(present(form3))then
           call pssetcolor(frame, form2)
        else
           call psdot(frame)
        endif
        call psplotcurve(frame, x3, y3)
     endif
     call psend(frame)
   end subroutine psplottabfunc


 End module ps_utils

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!13 Standard PS Fonts
!!Times-Roman
!!Helvetica
!!Courier
!!Symbol
!!Times-Italic
!!Helvetica-Oblique
!!Courier-Oblique
!!Times-Bold
!!Helvetica-Bold
!!Courier-Bold
!!Times-BoldItalic
!!Helvetica-BoldOblique
!!Courier-BoldOblique
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

