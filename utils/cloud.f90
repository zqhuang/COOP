module cloud_utils
  use basic_utils
  use file_io_utils
  use general_utils
  implicit none

  !!for 3d cloud interpolation 
#define CLOUD_DIM 3

#if CLOUD_DIM == 3
#define CLOUD_INDS i1, i2, i3
#define CLOUD_LOOP  do i3 = 1, cloud_np; do i2 = 1,cloud_np; do i1 = 1, cloud_np
#define CLOUD_ARR_INDS i(1), i(2), i(3)
#define CLOUD_REC (((i3-1)*cloud_np+i2-1)*cloud_np + i1)
#define CLOUD_ENDLOOP enddo; enddo; enddo
#elif CLOUD_DIM == 2
#define CLOUD_INDS i1, i2
#define CLOUD_LOOP  do i2=1,cloud_np; do i1 = 1, cloud_np
#define CLOUD_ARR_INDS i(1), i(2)
#define CLOUD_REC ((i2-1)*cloud_np + i1)
#define CLOUD_ENDLOOP enddo; enddo
#endif


#define CLOUD_PT c%g(CLOUD_ARR_INDS)
#define CLOUD_POINT c%g(CLOUD_INDS)
#define CLOUD_DSQ(x, j) sum((x(1:CLOUD_DIM) - CLOUD_PT%p(1:CLOUD_DIM,j))**2)
#define CLOUD_DISTSQ(x, j) sum((x(1:CLOUD_DIM) - CLOUD_POINT%p(1:CLOUD_DIM,j))**2)

  integer,parameter::cloud_np = 32
  integer,parameter::cloud_ncolor = 3
  integer,parameter::cloud_capacity = 32
  integer,parameter::cloud_recl = 2 * ib + (cloud_ncolor + CLOUD_DIM)*cloud_capacity* dl
  real(dl),parameter::cloud_min_d = 1.e-5_dl
  real(dl),parameter::cloud_min_dsq = cloud_min_d**2

  integer::cloud_nprocs = 1
  integer::cloud_id = 0

  type cloud_grid
     logical updated
     integer n
     integer ex
     real(dl),dimension(CLOUD_DIM, cloud_capacity):: p
     real(dl),dimension(cloud_ncolor, cloud_capacity):: c
  end type cloud_grid

  type cloud_body
     integer total_n, total_ex
#if CLOUD_DIM == 3
     type(cloud_grid),dimension(cloud_np, cloud_np, cloud_np):: g
#elif CLOUD_DIM == 2
     type(cloud_grid),dimension(cloud_np, cloud_np):: g
#endif

  end type cloud_body


contains

  subroutine cloud_initialize(c)
    type(cloud_body)c
    c%g%n = 0
    c%g%ex = cloud_capacity
    c%g%updated = .false.
    c%total_n = 0
    c%total_ex = cloud_capacity * cloud_np**3
  end subroutine cloud_initialize

  subroutine cloud_get_expected(c, xgen)
    type(cloud_body)c
    external xgen
    real(dl) x(CLOUD_DIM)
    integer i(CLOUD_DIM), CLOUD_INDS, ntry
    c%g%ex = 0
    CLOUD_LOOP
    i = (/ CLOUD_INDS /)
    ntry = cloud_capacity*8
    call xgen(x, i, ntry)
    if(ntry .eq. 0)then
       CLOUD_PT%ex  = 0
    else
       CLOUD_PT%ex = (cloud_capacity*9 - 1 - ntry)/(cloud_capacity*8-ntry)
    endif
    CLOUD_ENDLOOP
    c%total_ex = sum(c%g%ex)
  end subroutine cloud_get_expected

  subroutine cloud_locate(x, i)
    real(dl) x(CLOUD_DIM)
    integer i(CLOUD_DIM)
    i = min(max(ceiling(x*cloud_np), 1), cloud_np)
  end subroutine cloud_locate

  subroutine cloud_mindsq_p2g(x, i, dsq)
    integer i(CLOUD_DIM)
    real(dl) x(CLOUD_DIM)
    real(dl) dsq, y1(CLOUD_DIM), y2(CLOUD_DIM)
    if(any(i.le.0) .or. any(i.gt. cloud_np))then
       dsq = 100.d0
       return
    endif
    y1 = x - (i-1.d0)/cloud_np
    y2 = x - dble(i)/cloud_np
    dsq = sum(min(abs(y1), abs(y2))**2, mask = (y1.gt.0 .and. y2.gt.0) .or. (y1.lt.0 .and. y2.lt.0) )
  end subroutine cloud_mindsq_p2g

    
  subroutine cloud_maxdsq_p2g(x, i, dsq)
    integer i(CLOUD_DIM)
    real(dl) x(CLOUD_DIM)
    real(dl) dsq, y1(CLOUD_DIM), y2(CLOUD_DIM)
    y1 = abs(x - (i-1.d0)/cloud_np)
    y2 = abs(x - dble(i)/cloud_np)
    dsq = sum(max(y1, y2)**2)
  end subroutine cloud_maxdsq_p2g


  subroutine cloud_find_nearest(c, x, n, p, d, dsq, nf)
    type(cloud_body) c
    real(dl),intent(IN)::x(CLOUD_DIM)
    integer,intent(IN)::n
    integer,intent(OUT)::nf
    real(dl),intent(OUT):: dsq(n), p(CLOUD_DIM, n), d(cloud_ncolor, n)

    real(dl)  maxdsq, d2, mindsq(3**CLOUD_DIM)
    integer i(CLOUD_DIM), j, k, imax
    integer loc(CLOUD_DIM)
    integer, dimension(CLOUD_DIM, 3**CLOUD_DIM),parameter:: step = reshape( (/ &
#if CLOUD_DIM == 3
         0, 0, 0,   0, 0, 1,   0, 1, 0,   1, 0, 0, &
         0, 0,-1,   0,-1, 0,  -1, 0, 0,   1, 1, 0, &
         1, 0, 1,   0, 1, 1,  -1,-1, 0,  -1, 0,-1,  &
         0,-1,-1,   -1, 1, 0,  1, -1, 0, -1, 0, 1,  &
         1, 0, -1,   0, -1, 1,  0, 1, -1,  1, 1, 1,  &
         -1,-1,-1,   1, 1, -1,  1, -1, 1,  -1, 1, 1, &
         -1, -1, 1,  -1, 1, -1,  1, -1, -1 &
#else
         0, 0,    0, 1,  1, 0,  0,-1,  -1, 0, &
         1, 1,    -1,-1,  1, -1, -1, 1 &
#endif
         /), (/ CLOUD_DIM, 3**CLOUD_DIM /) )

    call cloud_locate(x, loc)
    do k=1, 3**CLOUD_DIM
       call cloud_mindsq_p2g(x, loc+step(:, k), mindsq(k))
    enddo
    nf = 0
    maxdsq = 0
    do k = 1, 3**CLOUD_DIM
       if(mindsq(k) .gt. 1.d0 .or. (nf.eq.n .and. mindsq(k) .le. maxdsq) ) cycle
       i = loc+step(:,k)
       do j = 1, CLOUD_PT%n
          d2 = CLOUD_DSQ(x, j)
          if(nf.lt.n)then
             nf = nf + 1
             p(:, nf) = CLOUD_PT%p(:,j)
             d(:, nf) = CLOUD_PT%c(:,j)
             dsq(nf) = d2
             if(d2 .gt. maxdsq)then
                maxdsq = d2
                imax = nf
             endif
          else
             if(d2 .lt. maxdsq)then
                p(:, imax) = CLOUD_PT%p(:,j)
                d(:, imax) = CLOUD_PT%c(:,j)
                dsq(imax) = d2
                imax = imaxloc(dsq)
                maxdsq = dsq(imax)
             endif
          endif
       enddo
    enddo
   end subroutine cloud_find_nearest


  function cloud_sparseness(c, x) result(s)
    type(cloud_body) c
    real(dl),intent(IN):: x(CLOUD_DIM)
    real(dl) p(CLOUD_DIM, 4), d(cloud_ncolor, 4), dsq(4), s
    integer nf
    call cloud_find_nearest(c, x, 4, p, d, dsq, nf)
    if(nf .eq. 0)then
       s = 1.d0
    else
       s = product(dsq(1:nf))**0.125d0
    endif
  end function cloud_sparseness



  function cloud_pick_sparse(c, x, xgen, filter) result(success)
    integer,parameter::max_nsteps = 256
    integer,parameter::max_nhit = 8
    type(cloud_body) c
    integer CLOUD_INDS, i(CLOUD_DIM)
    real(dl) x(CLOUD_DIM), xtest(CLOUD_DIM)
    external xgen
    external filter
    logical filter
    real(dl) s, maxs
    integer nsteps, nhit, ntry
    real(dl) rat, r
    logical success
    maxs = 0.d0
    nhit = 0
    nsteps = 0
    rat = cloud_capacity
    i = 0
    !!pick out the lowest ratio
    CLOUD_LOOP
    if(filter(CLOUD_INDS) .and. CLOUD_POINT%ex .gt. 0 .and. CLOUD_POINT%n .lt. cloud_capacity)then
       r = dble(CLOUD_POINT%n)/CLOUD_POINT%ex
       if(r .lt. rat)then
          rat = r
          i = (/ CLOUD_INDS /)
       endif
    endif
    CLOUD_ENDLOOP
    if(all(i.eq.0))then
       success = .false.
       return
    endif
    do 
       ntry = cloud_capacity * 64
       call xgen(xtest, i, ntry)
       if(ntry.eq.0)then
          success = (nhit .gt. 0)          
          return 
       endif
       s = cloud_sparseness(c, xtest)
       if(s .gt. maxs)then
          x = xtest
          maxs = s
          nhit = nhit + 1
       endif
       nsteps = nsteps + 1
       if(nsteps .ge. max_nsteps .or. nhit .ge. max_nhit)exit
    enddo
    success = (nhit .gt. 0)
  end function cloud_pick_sparse

  function cloud_insert_sym(c, x, d) result(success)
    logical success
    type(cloud_body) c
    real(dl) x(CLOUD_DIM), d(cloud_ncolor)
#if CLOUD_DIM == 3
    logical s(6)
    s(1) = cloud_insert(c, x, d)
    s(2) = cloud_insert(c, (/ x(2), x(3), x(1) /), d)
    s(3) = cloud_insert(c, (/ x(3), x(1), x(2) /), d)
    s(4) = cloud_insert(c, (/ x(3), x(2), x(1) /), d)
    s(5) = cloud_insert(c, (/ x(1), x(3), x(2) /), d)
    s(6) = cloud_insert(c, (/ x(2), x(1), x(3) /), d)
#elif CLOUD_DIM == 2
    logical s(2)
    s(1) = cloud_insert(c, x, d)
    s(2) = cloud_insert(c, (/ x(2), x(1) /), d)
#endif
    success  = any(s)
  end function cloud_insert_sym

  !!-1 < x1, x2, x3 < 1
  function cloud_insert(c, x, d) result(success)
    logical success
    type(cloud_body) c
    real(dl) x(CLOUD_DIM), d(cloud_ncolor)
    integer i(CLOUD_DIM), j
    call cloud_locate(x, i)
    if(CLOUD_PT%n .eq. cloud_capacity)then
       success = .false.
       return
    endif
    do j = 1, CLOUD_PT%n
       if( CLOUD_DSQ(x, j) .lt. cloud_min_dsq)then
          success = .false.
          return
       endif
    enddo
    CLOUD_PT%n = CLOUD_PT%n + 1
    CLOUD_PT%p(:, CLOUD_PT%n) = x
    CLOUD_PT%c(:, CLOUD_PT%n) = d
    c%total_n  = c%total_n + 1
    success = .true.
    CLOUD_PT%updated = .true.
    if(CLOUD_PT%ex .lt. 1) CLOUD_PT%ex = 1
  end function cloud_insert


  subroutine cloud_dump(c, fname)
    type(file_pointer) fp
    character(LEN=*) fname
    type(cloud_body) c
    integer CLOUD_INDS, j
    fp = open_file(trim(fname), 'b', recl = cloud_recl ) 
    j = 0
    CLOUD_LOOP
    j = j + 1
    write(fp%unit, rec = j) CLOUD_POINT%ex, CLOUD_POINT%n, CLOUD_POINT%p, CLOUD_POINT%c
    CLOUD_ENDLOOP
    call close_file(fp)
  end subroutine cloud_dump


  function cloud_load(c, fname) result(loaded)
    type(file_pointer) fp
    character(LEN=*) fname
    type(cloud_body) c
    integer CLOUD_INDS, j
    logical loaded
    if(file_exists(fname))then
       fp = open_file(trim(fname), 'br', recl = cloud_recl)
       j = 0
       CLOUD_LOOP
       j = j + 1
       read(fp%unit, rec = j) CLOUD_POINT%ex, CLOUD_POINT%n, CLOUD_POINT%p, CLOUD_POINT%c
       CLOUD_ENDLOOP
       call close_file(fp)
       c%total_n = sum(c%g%n)
       c%total_ex = sum(c%g%ex)
       c%g%updated = .false.
       loaded = .true.
    else
       write(*,*) "warning: "//trim(fname)//" does not exist, initializing a null file"
       call cloud_initialize(c)
       loaded = .false.
    endif
  end function cloud_load

  subroutine cloud_update(c, fname, filter)
    integer CLOUD_INDS
    type(cloud_body) c
    type(file_pointer) fp
    character(LEN=*) fname
    external filter
    logical filter
    fp = open_file(trim(fname), 'b', recl = cloud_recl)
    CLOUD_LOOP
    if(filter(CLOUD_INDS) .and. CLOUD_POINT%updated)then
       write(fp%unit, rec = CLOUD_REC) CLOUD_POINT%ex, CLOUD_POINT%n, CLOUD_POINT%p, CLOUD_POINT%c
       CLOUD_POINT%updated = .false.
    endif
    CLOUD_ENDLOOP
    call close_file(fp)
  end subroutine cloud_update

  function cloud_no_filter(CLOUD_INDS) result(f)
    integer CLOUD_INDS
    logical f
    f = .true.
    return
  end function cloud_no_filter

  function cloud_default_filter(CLOUD_INDS) result(f)
    integer CLOUD_INDS
    logical f
#if CLOUD_DIM == 3
    f = ( mod(i1+i2+i3, cloud_nprocs) .eq. cloud_id )
#elif CLOUD_DIM == 2
    f = ( mod(i1+i2, cloud_nprocs) .eq. cloud_id )
#endif
  end function cloud_default_filter


  subroutine cloud_triangle_xgen(x, i, ntry) 
    real(dl) x(CLOUD_DIM)
    integer i(CLOUD_DIM), ntry
100 call random_number(x)
    x = (i-1+x)/cloud_np
    ntry = ntry - 1
#if CLOUD_DIM == 3
    if (x(1)+x(2) .lt. x(3) .or. x(2)+x(3).lt. x(1) .or. x(1)+x(3) .lt. x(2))then
       if(ntry .eq. 0)then
          x = -1.d0
          return
       endif
       goto 100
    endif
#endif
  end subroutine cloud_triangle_xgen



  subroutine cloud_default_xgen(x, i, ntry) 
    real(dl) x(CLOUD_DIM)
    integer i(CLOUD_DIM)
    integer ntry
    call random_number(x)
    x = (i-1+x)/cloud_np
    ntry = ntry - 1
  end subroutine cloud_default_xgen

  subroutine cloud_print(c, fname)
    type(cloud_body) c
    integer CLOUD_INDS, j
    character(LEN=1024) fmt
    type(file_pointer) fp
    character(LEN=*) fname
    fmt = '('//trim(num2str(CLOUD_DIM))//'F9.5, '//trim(num2str(cloud_ncolor))//'E14.5)'
    fp = open_file(fname, "w")
    CLOUD_LOOP
    do j = 1, CLOUD_POINT%n
       write(fp%unit,trim(fmt))CLOUD_POINT%p(:, j), CLOUD_POINT%c(:,j)
    enddo
    CLOUD_ENDLOOP
    call close_file(fp)
  end subroutine cloud_print

  subroutine cloud_interpolation(c, x, y)
#if CLOUD_DIM == 3
    integer,parameter::nd = 10
#elif CLOUD_DIM == 2
    integer,parameter::nd = 6
#endif
    integer,parameter::n = nd * 2
    type(cloud_body) c
    real(dl) coef(nd, cloud_ncolor)
    real(dl) x(CLOUD_DIM), y(cloud_ncolor)
    real(dl) p(CLOUD_DIM, n), d(cloud_ncolor, n), dsq(n), w(n), wtot
    integer nf, i
    real(dl),dimension(:,:), allocatable:: fx, ys
    call cloud_find_nearest(c, x, n, p, d, dsq, nf)
    select case(nf)
    case(0)
       y = 0.d0  !!default 0
    case(1) 
       y = d(:, 1) 
    case(2:n-1)  
       w(1:nf) = 1.d0/(dsq(1:nf) + cloud_min_dsq)
       wtot = sum(w(1:nf))
       w(1:nf) = w(1:nf) / wtot
       y = d(:,1)*w(1)
       do i=2, nf
          y = y + d(:, i)*w(i)
       enddo
    case default
       allocate(fx(nf, nd), ys(nf, cloud_ncolor))
       fx(:, 1) = 1.d0
       fx(:, 2) = p(1, 1:nf)
       fx(:, 3) = p(2, 1:nf)
#if CLOUD_DIM == 3
       fx(:, 4) = p(3, 1:nf)
       fx(:, 5) = p(1, 1:nf) ** 2
       fx(:, 6) = p(2, 1:nf) ** 2
       fx(:, 7) = p(3, 1:nf) ** 2
       fx(:, 8) = p(2, 1:nf) * p(3, 1:nf)
       fx(:, 9) = p(1, 1:nf) * p(3, 1:nf)
       fx(:, 10) = p(1, 1:nf) * p(2, 1:nf)
#elif CLOUD_DIM == 2
       fx(:, 4) = p(1, 1:nf)**2
       fx(:, 5) = p(2, 1:nf)**2
       fx(:, 6) = p(1, 1:nf) * p(2, 1:nf)
#endif
       dsq(1:nf) = 1.d0/sqrt(dsq(1:nf)+cloud_min_dsq)
       ys =  transpose(d(1:cloud_ncolor, 1:nf))
       do i=1, nf
          fx(i, :) = fx(i, :)*dsq(i)
          ys(i, :) = ys(i, :)*dsq(i)
       enddo
       call svd_least_square_all(nf, nd, cloud_ncolor, fx, ys, coef)  !!this is small matrix operation, I don't call lapack
#if CLOUD_DIM == 3
       y = matmul( (/ 1.d0, x(1), x(2), x(3), x(1)**2, x(2)**2, x(3)**2, x(2)*x(3), x(1)*x(3), x(1)*x(2) /), coef)
#elif CLOUD_DIM == 2
       y = matmul( (/ 1.d0, x(1), x(2), x(1)**2, x(2)**2, x(1)*x(2) /), coef)
#endif
       deallocate(fx, ys)
    end select

  end subroutine cloud_interpolation



  subroutine cloud_interpolation_sym(c, x, y)
#if CLOUD_DIM == 3
    integer,parameter::nd = 7
#elif CLOUD_DIM == 2
    integer,parameter::nd = 6
#endif
    integer,parameter::n = nd * 2
    type(cloud_body) c
    real(dl) coef(nd, cloud_ncolor)
    real(dl) x(CLOUD_DIM), y(cloud_ncolor)
    real(dl) p(CLOUD_DIM, n), d(cloud_ncolor, n), dsq(n), w(n), wtot
    integer nf, i
    real(dl),dimension(:,:), allocatable:: fx, ys
    call cloud_find_nearest(c, x, n, p, d, dsq, nf)
    select case(nf)
    case(0)
       y = 0.d0  !!default 0
    case(1) 
       y = d(:, 1) 
    case(2:n-1)  
       w(1:nf) = 1.d0/(dsq(1:nf) + cloud_min_dsq)
       wtot = sum(w(1:nf))
       w(1:nf) = w(1:nf) / wtot
       y = d(:,1)*w(1)
       do i=2, nf
          y = y + d(:, i)*w(i)
       enddo
    case default
       allocate(fx(nf, nd), ys(nf, cloud_ncolor))
       fx(:, 1) = 1.d0
#if CLOUD_DIM == 3
       fx(:, 2) = p(1, 1:nf) +  p(2, 1:nf) +  p(3, 1:nf)
       fx(:, 3) = p(1, 1:nf)**2 +  p(2, 1:nf)**2 +  p(3, 1:nf)**2
       fx(:, 4) = ( p(1, 1:nf) + p(3,1:nf) ) *  p(2, 1:nf) +  p(3, 1:nf)*p(1, 1:nf)
       fx(:, 5) = p(1, 1:nf)**3 + p(2, 1:nf)**3 + p(3, 1:nf) **3
       fx(:, 6) = p(1, 1:nf) * p(2, 1:nf) * p(3, 1:nf)
       fx(:, 7) = p(1, 1:nf) * ( p(2, 1:nf) * (p(1, 1:nf) + p(2, 1:nf)) &
            + p(3, 1:nf) * (p(1, 1:nf) + p(3, 1:nf)) ) &
            + p(3, 1:nf) * p(2, 1:nf) * (p(3, 1:nf) + p(2, 1:nf)) 
#elif CLOUD_DIM == 2
       fx(:, 2) = p(1, 1:nf)+p(2, 1:nf)
       fx(:, 3) = p(1, 1:nf)**2 + p(2, 1:nf)**2
       fx(:, 4) = p(1, 1:nf)*p(2,1:nf)
       fx(:, 5) = p(1, 1:nf)**3 + p(2, 1:nf) **3
       fx(:, 6) = p(1, 1:nf) *  p(2, 1:nf) * (p(1, 1:nf) + p(2, 1:nf))
#endif
       dsq(1:nf) = 1.d0/sqrt(dsq(1:nf)+cloud_min_dsq)
       ys =  transpose(d(1:cloud_ncolor, 1:nf))
       do i=1, nf
          fx(i, :) = fx(i, :)*dsq(i)
          ys(i, :) = ys(i, :)*dsq(i)
       enddo
       call svd_least_square_all(nf, nd, cloud_ncolor, fx, ys, coef)  !!this is small matrix operation, I don't call lapack
#if CLOUD_DIM == 3
       y = matmul( (/ 1.d0, sum(x),  sum(x**2) , (x(2)+x(1))*x(3)+x(1)*x(2), sum(x**3), product(x), x(1)*(x(2)*(x(1)+x(2))+x(3)*((x(1)+x(3)))) + x(2)*x(3)*(x(2)+x(3)) /), coef)
#elif CLOUD_DIM == 2
       y = matmul( (/ 1.d0, sum(x), sum(x**2), product(x), sum(x**3), x(1)*x(2)*(x(1)+x(2)) /), coef)
#endif
       deallocate(fx, ys)
    end select

  end subroutine cloud_interpolation_sym



end module cloud_utils
