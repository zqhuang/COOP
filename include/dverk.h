    dimension y(n),c(24), w(nw,9)
    external fcn
!c
!c     ******************************************************************
!c     * begin initialization, parameter checking, interrupt re-entries *
!c     ******************************************************************
!c
!c  ......abort if ind out of range 1 to 6
    if (ind.lt.1 .or. ind.gt.6) goto 500
!c
!c        cases - initial entry, normal re-entry, interrupt re-entries
    goto (5, 5, 45, 1111, 2222, 2222), ind
!c        case 1 - initial entry (ind .eq. 1 or 2)
!c  .........abort if n.gt.nw or tol.le.0
5   if (n.gt.nw .or. tol.le.0.d0) go to 500
    if (ind.eq. 2) go to 15
!c              initial entry without options (ind .eq. 1)
!c              set c(1) to c(9) equal to 0
    do k = 1, 9
       c(k) = 0.d0
    enddo
10 continue
    go to 35
15  continue
!c              initial entry with options (ind .eq. 2)
!c              make c(1) to c(9) non-negative
    do  k = 1, 9
       c(k) = dabs(c(k))
    enddo
20  continue
!c              make floor values non-negative if they are to be used
    if (c(1).ne.4.d0 .and. c(1).ne.5.d0) go to 30
    do  k = 1, n
       c(k+30) = dabs(c(k+30))
    enddo
25  continue
30  continue
35  continue
!c           initialize rreb, dwarf, prev xend, flag, counts
    c(10) = 2.d0**(-56)
    c(11) = 1.d-35
!c           set previous xend initially to initial value of x
    c(20) = x
    do  k = 21, 24
       c(k) = 0.d0
    enddo
40  continue
    go to 50
!c        case 2 - normal re-entry (ind .eq. 3)
!c  .........abort if xend reached, and either x changed or xend not
45  if (c(21).ne.0.d0 .and. & 
         (x.ne.c(20) .or. xend.eq.c(20))) go to 500
!c           re-initialize flag
    c(21) = 0.d0
    go to 50
!c        case 3 - re-entry following an interrupt (ind .eq. 4 to 6)
!c           transfer control to the appropriate re-entry point..........
!c           this has already been handled by the computed go to        .
!c        end cases                                                     v
50  continue
!c
!c     end initialization, etc.
!c
!c     ******************************************************************
!c     * loop through the following 4 stages, once for each trial  step *
!c     * until the occurrence of one of the following                   *
!c     *    (a) the normal return (with ind .eq. 3) on reaching xend in *
!c     *        stage 4                                                 *
!c     *    (b) an error return (with ind .lt. 0) in stage 1 or stage 4 *
!c     *    (c) an interrupt return (with ind  .eq.  4,  5  or  6),  if *
!c     *        requested, in stage 1 or stage 4                        *
!c     ******************************************************************
!c
99999 continue
!c
!c        ***************************************************************
!c        * stage 1 - prepare - do calculations of  hmin,  hmax,  etc., *
!c        * and some parameter  checking,  and  end  up  with  suitable *
!c        * values of hmag, xtrial and htrial in preparation for taking *
!c        * an integration step.                                        *
!c        ***************************************************************
!c
!c***********error return (with ind=-1) if no of fcn evals too great
    if (c(7).eq.0.d0 .or. c(24).lt.c(7)) go to 100
    ind = -1
    return
100 continue
!c
!c           calculate slope (adding 1 to no of fcn evals) if ind .ne. 6
    if (ind .eq. 6) go to 105
    call fcn(n, x, y, w(1,1) DVERK_ARGUMENTS)
    c(24) = c(24) + 1.d0
105 continue
!c
!c           calculate hmin - use default unless value prescribed
    c(13) = c(3)
    if (c(3) .ne. 0.d0) go to 165
!c              calculate default value of hmin
!c              first calculate weighted norm y - c(12) - as specified
!c              by the error control indicator c(1)
    temp = 0.d0
    if (c(1) .ne. 1.d0) go to 115
!c                 absolute error control - weights are 1
    do k = 1, n
       temp = dmax1(temp, dabs(y(k)))
    enddo
110 continue
    c(12) = temp
    go to 160
115 if (c(1) .ne. 2.d0) go to 120
!c                 relative error control - weights are 1/dabs(y(k)) so
!c                 weighted norm y is 1
    c(12) = 1.d0
    go to 160
120 if (c(1) .ne. 3.d0) go to 130
!c                 weights are 1/max(c(2),abs(y(k)))
    do  k = 1, n
       temp = dmax1(temp, dabs(y(k))/c(2))
    enddo
125 continue
    c(12) = dmin1(temp, 1.d0)
    go to 160
130 if (c(1) .ne. 4.d0) go to 140
!c                 weights are 1/max(c(k+30),abs(y(k)))
    do k = 1, n
       temp = dmax1(temp, dabs(y(k))/c(k+30))
    enddo
135 continue
    c(12) = dmin1(temp, 1.d0)
    go to 160
140 if (c(1) .ne. 5.d0) go to 150
!c                 weights are 1/c(k+30)
    do k = 1, n
       temp = dmax1(temp, dabs(y(k))/c(k+30))
    enddo
145 continue
    c(12) = temp
    go to 160
150 continue
!c                 default case - weights are 1/max(1,abs(y(k)))
    do  k = 1, n
       temp = dmax1(temp, dabs(y(k)))
    enddo
155 continue
    c(12) = dmin1(temp, 1.d0)
160 continue
    c(13) = 10.d0*dmax1(c(11),c(10)*dmax1(c(12)/tol,dabs(x)))
165 continue
 !c
!c           calculate scale - use default unless value prescribed
    c(15) = c(5)
    if (c(5) .eq. 0.d0) c(15) = 1.d0
!c
!c           calculate hmax - consider 4 cases
!c           case 1 both hmax and scale prescribed
    if (c(6).ne.0.d0 .and. c(5).ne.0.d0) & 
         c(16) = dmin1(c(6), 2.d0/c(5))
!c           case 2 - hmax prescribed, but scale not
    if (c(6).ne.0.d0 .and. c(5).eq.0.d0) c(16) = c(6)
!c           case 3 - hmax not prescribed, but scale is
    if (c(6).eq.0.d0 .and. c(5).ne.0.d0) c(16) = 2.d0/c(5)
!c           case 4 - neither hmax nor scale is provided
    if (c(6).eq.0.d0 .and. c(5).eq.0.d0) c(16) = 2.d0
!c
!c***********error return (with ind=-2) if hmin .gt. hmax
    if (c(13) .le. c(16)) go to 170
    ind = -2
    return
170 continue
!c
!c           calculate preliminary hmag - consider 3 cases
    if (ind .gt. 2) go to 175
!c           case 1 - initial entry - use prescribed value of hstart, if
!c              any, else default
    c(14) = c(4)
    if (c(4) .eq. 0.d0) c(14) = c(16)*tol**(1./6.)
    go to 185
175 if (c(23) .gt. 1.d0) go to 180
!c           case 2 - after a successful step, or at most  one  failure,
!c              use min(2, .9*(tol/est)**(1/6))*hmag, but avoid possible
!c              overflow. then avoid reduction by more than half.
    temp = 2.d0*c(14)
    if (tol .lt. (2.d0/.9d0)**6*c(19)) & 
         temp = .9d0*(tol/c(19))**(1./6.)*c(14)
    c(14) = dmax1(temp, .5d0*c(14))
    go to 185
180 continue
!c           case 3 - after two or more successive failures
    c(14) = .5d0*c(14)
185 continue
!c
!c           check against hmax
    c(14) = dmin1(c(14), c(16))
!c
!c           check against hmin
    c(14) = dmax1(c(14), c(13))
!c
!c***********interrupt no 1 (with ind=4) if requested
    if (c(8) .eq. 0.d0) go to 1111
    ind = 4
    return
!c           resume here on re-entry with ind .eq. 4   ........re-entry..
1111 continue
!c
!c           calculate hmag, xtrial - depending on preliminary hmag, xend
    if (c(14) .ge. dabs(xend - x)) go to 190
!c              do not step more than half way to xend
    c(14) = dmin1(c(14), .5d0*dabs(xend - x))
    c(17) = x + dsign(c(14), xend - x)
    go to 195
190 continue
!c              hit xend exactly
    c(14) = dabs(xend - x)
    c(17) = xend
195 continue
!c
!c           calculate htrial
    c(18) = c(17) - x
!c
!c        end stage 1
!c
!c        ***************************************************************
!c        * stage 2 - calculate ytrial (adding 7 to no of  fcn  evals). *
!c        * w(*,2), ... w(*,8)  hold  intermediate  results  needed  in *
!c        * stage 3. w(*,9) is temporary storage until finally it holds *
!c        * ytrial.                                                     *
!c        ***************************************************************
!c
    temp = c(18)/1398169080000.d0
!c
    do k = 1, n
       w(k,9) = y(k) + temp*w(k,1)*233028180000.d0
    enddo
200 continue
    call fcn(n, x + c(18)/6.d0, w(1,9), w(1,2) DVERK_ARGUMENTS)
!c
    do k = 1, n
       w(k,9) = y(k) + temp*(   w(k,1)*74569017600.d0 & 
            + w(k,2)*298276070400.d0  )
    enddo
205 continue
    call fcn(n, x + c(18)*(4.d0/15.d0), w(1,9), w(1,3) DVERK_ARGUMENTS)
!c
    do  k = 1, n
       w(k,9) = y(k) + temp*(   w(k,1)*1165140900000.d0 & 
            - w(k,2)*3728450880000.d0 & 
            + w(k,3)*3495422700000.d0 )
    enddo
210 continue
    call fcn(n, x + c(18)*(2.d0/3.d0), w(1,9), w(1,4) DVERK_ARGUMENTS)
!c
    do k = 1, n
       w(k,9) = y(k) + temp*( - w(k,1)*3604654659375.d0 & 
            + w(k,2)*12816549900000.d0 & 
            - w(k,3)*9284716546875.d0 & 
            + w(k,4)*1237962206250.d0 )
    enddo
    continue
    call fcn(n, x + c(18)*(5.d0/6.d0), w(1,9), w(1,5) DVERK_ARGUMENTS)
!c
    do k = 1, n
       w(k,9) = y(k) + temp*(   w(k,1)*3355605792000.d0 & 
            - w(k,2)*11185352640000.d0 & 
            + w(k,3)*9172628850000.d0 & 
            - w(k,4)*427218330000.d0 & 
            + w(k,5)*482505408000.d0  )
    enddo
220 continue
    call fcn(n, x + c(18), w(1,9), w(1,6) DVERK_ARGUMENTS)
!c
    do k = 1, n
       w(k,9) = y(k) + temp*( - w(k,1)*770204740536.d0 & 
            + w(k,2)*2311639545600.d0 & 
            - w(k,3)*1322092233000.d0 & 
            - w(k,4)*453006781920.d0 & 
            + w(k,5)*326875481856.d0  )
    enddo
225    continue
    call fcn(n, x + c(18)/15.d0, w(1,9), w(1,7) DVERK_ARGUMENTS)
!c
    do k = 1, n
       w(k,9) = y(k) + temp*(   w(k,1)*2845924389000.d0 & 
            - w(k,2)*9754668000000.d0 & 
            + w(k,3)*7897110375000.d0 & 
            - w(k,4)*192082660000.d0 & 
            + w(k,5)*400298976000.d0 & 
            + w(k,7)*201586000000.d0  )
    enddo
230 continue
    call fcn(n, x + c(18), w(1,9), w(1,8) DVERK_ARGUMENTS)
!c
!c           calculate ytrial, the extrapolated approximation and store
!c              in w(*,9)
    do k = 1, n
       w(k,9) = y(k) + temp*(   w(k,1)*104862681000.d0 & 
            + w(k,3)*545186250000.d0 & 
            + w(k,4)*446637345000.d0 & 
            + w(k,5)*188806464000.d0 & 
            + w(k,7)*15076875000.d0 & 
            + w(k,8)*97599465000.d0   )
    enddo
235    continue
!c
!c           add 7 to the no of fcn evals
    c(24) = c(24) + 7.d0
!c
!c        end stage 2
!c
!c        ***************************************************************
!c        * stage 3 - calculate the error estimate est. first calculate *
!c        * the  unweighted  absolute  error  estimate vector (per unit *
!c        * step) for the unextrapolated approximation and store it  in *
!c        * w(*,2).  then  calculate the weighted max norm of w(*,2) as *
!c        * specified by the error  control  indicator  c(1).  finally, *
!c        * modify  this result to produce est, the error estimate (per *
!c        * unit step) for the extrapolated approximation ytrial.       *
!c        ***************************************************************
!c
!c           calculate the unweighted absolute error estimate vector
    do k = 1, n
       w(k,2) = (   w(k,1)*8738556750.d0 & 
            + w(k,3)*9735468750.d0 & 
            - w(k,4)*9709507500.d0 & 
            + w(k,5)*8582112000.d0 & 
            + w(k,6)*95329710000.d0 & 
            - w(k,7)*15076875000.d0 & 
            - w(k,8)*97599465000.d0)/1398169080000.d0
    enddo
300 continue
!c
!c           calculate the weighted max norm of w(*,2) as specified by
!c           the error control indicator c(1)
    temp = 0.d0
    if (c(1) .ne. 1.d0) go to 310
!c              absolute error control
    do  k = 1, n
       temp = dmax1(temp,dabs(w(k,2)))
    enddo
305 continue
    go to 360
310 if (c(1) .ne. 2.d0) go to 320
!c              relative error control
    do k = 1, n
       temp = dmax1(temp, dabs(w(k,2)/y(k)))
    enddo
315 continue
    go to 360
320 if (c(1) .ne. 3.d0) go to 330
!c              weights are 1/max(c(2),abs(y(k)))
    do k = 1, n
       temp = dmax1(temp, dabs(w(k,2)) & 
            / dmax1(c(2), dabs(y(k))) )
    enddo
325 continue
    go to 360
330 if (c(1) .ne. 4.d0) go to 340
!c              weights are 1/max(c(k+30),abs(y(k)))
    do k = 1, n
       temp = dmax1(temp, dabs(w(k,2)) & 
            / dmax1(c(k+30), dabs(y(k))) )
    enddo
335 continue
    go to 360
340 if (c(1) .ne. 5.d0) go to 350
!c              weights are 1/c(k+30)
    do k = 1, n
       temp = dmax1(temp, dabs(w(k,2)/c(k+30)))
    enddo
345 continue
    go to 360
350 continue
!c              default case - weights are 1/max(1,abs(y(k)))
    do k = 1, n
       temp = dmax1(temp, dabs(w(k,2)) & 
            / dmax1(1.d0, dabs(y(k))) )
    enddo
355 continue
360 continue
!c
!c           calculate est - (the weighted max norm of w(*,2))*hmag*scale
!c              - est is intended to be a measure of the error  per  unit
!c              step in ytrial
    c(19) = temp*c(14)*c(15)
!c
!c        end stage 3
!c
!c        ***************************************************************
!c        * stage 4 - make decisions.                                   *
!c        ***************************************************************
!c
!c           set ind=5 if step acceptable, else set ind=6
    ind = 5
    if (c(19) .gt. tol) ind = 6
!c
!c***********interrupt no 2 if requested
    if (c(9) .eq. 0.d0) go to 2222
    return
!c           resume here on re-entry with ind .eq. 5 or 6   ...re-entry..
2222 continue
!c
    if (ind .eq. 6) go to 410
!c              step accepted (ind .eq. 5), so update x, y from xtrial,
!c                 ytrial, add 1 to the no of successful steps, and set
!c                 the no of successive failures to zero
    x = c(17)
    do k = 1, n
       y(k) = w(k,9)
    enddo
400 continue
    c(22) = c(22) + 1.d0
    c(23) = 0.d0
!c**************return(with ind=3, xend saved, flag set) if x .eq. xend
    if (x .ne. xend) go to 405
    ind = 3
    c(20) = xend
    c(21) = 1.d0
    return
405 continue
    go to 420
410 continue
!c              step not accepted (ind .eq. 6), so add 1 to the no of
!c                 successive failures
    c(23) = c(23) + 1.d0
!c**************error return (with ind=-3) if hmag .le. hmin
    if (c(14) .gt. c(13)) go to 415
    ind = -3
    return
415 continue
420 continue
!c
!c        end stage 4
!c
    go to 99999
!c     end loop
!c
!c  begin abort action
500 continue
!c
    write(*,*) 'dverk crashes with ind =',  ind
    write(*,*) 'x = ', x
    call fcn(n, x, y, w(1, 1) DVERK_ARGUMENTS)
    write(*, "(2A16)") "  y  ",  " dy/dx  "
    do i=1, min(n, 100)
       write(*, "(2E16.7)"), y(i), w(i, 1)
    enddo
    stop
!c
!c  end abort action
!c
