program crack2048
  use wrap_utils
  implicit none
#include "utils.h"
  integer,parameter:: risky_level = 4  !! 0 <= risky_level <= 10, recommended 3-6
  integer,parameter:: think_depth = 8
  integer, parameter::n=4
  integer, parameter::up = 1
  integer, parameter::down = 2
  integer, parameter::left = 3
  integer, parameter::right = 4
  character(LEN=5), dimension(4)::names = (/ "up   ", "down ", "left ", "right" /)
  logical::save_and_exit = .false.
  integer mygrids(n, n), previous(n,n), i, j, k
  STRING inline
  type(list_integer) il
  real gain(4)
  real,parameter:: prob4 = 0.21
  integer loc(1)
  type(file_pointer) fp
  call global_initialization(nobessel = .true.)
  if(risky_level .lt. 0 .or. risky_level .gt. 10) stop "risky_level must be set between 0 and 10"
50 fp = open_file("2048.ini","r")
  do i=1, n
     read(fp%unit, *) mygrids(i,:)
  enddo
  call close_file(fp)
100 continue
  print*,"+++++++++++++++++++++"
  do i =  1, 4
     gain(i)= mcmc_gain(mygrids, i, think_depth)
     write(*,"(A6, G12.4)") trim(names(i)), gain(i)
  enddo
  loc = maxloc(gain)
  write(*,*) "*** Move "//trim(names(loc(1)))//"! ****"
  previous = mygrids
  call move(mygrids, loc(1))
  call show(previous, mygrids)
  write(*,*) "entropy = ", entropy(mygrids)
  write(*,*) "Randomly fill:"
200 read(*,"(A)") inline
  inline = adjustl(inline)
  select case(trim(inline))
  case("1", "11")
     i = 1
     j = 1
  case("2","22")
     i = 1
     j = 2
  case("3","33")
     i = 1
     j = 3
  case("4", "44")
     i = 1
     j = 4
  case("q","qq")
     i = 2
     j = 1
  case("w","ww")
     i = 2
     j = 2
  case("e","ee")
     i = 2
     j = 3
  case("r","rr")
     i = 2
     j = 4
  case("a","aa")
     i = 3
     j = 1
  case("s","ss")
     i = 3
     j = 2
  case("d","dd")
     i = 3
     j = 3
  case("f","ff")
     i = 3
     j = 4
  case("z", "zz")
     i = 4
     j = 1
  case("x","xx")
     i = 4
     j = 2
  case("c", "cc")
     i = 4
     j = 3
  case("v", "vv")
     i = 4
     j = 4
  case("reload")
     goto 50
  case("modify")
     write(*,*) "Enter location and value (x, y, value) (enter value < -1 to finish):"
     read(*,*) i, j, k
     do while(k .ge. 0) 
        if(k.eq. 0 .and. i.ge. 1 .and. j.ge.1 .and. i.le.n .and. j.le.n)then
           mygrids(i, j) = k
           call show(mygrids)
        elseif(is_integer(log2(k*1.d0)).and. i.ge. 1 .and. j.ge.1 .and. i.le.n .and. j.le.n)then
           mygrids(i, j) = k
           call show(mygrids)
        else
           write(*,*) "Wrong input values"
        endif
        write(*,*) "Enter location and value (x, y, value) (enter value < -1 to finish):"
        read(*,*) i, j, k
     enddo
     goto 200
  case("save")
     save_and_exit = .true.
     goto 200
  case default
     write(*,*) "wrong input"
     goto 200
  end select
  if(mygrids(i,j).ne.0)then
     write(*,*) "wrong location"
     goto 200
  endif
  mygrids(i,j)= len_trim(inline)*2
  if(save_and_exit)then
     fp = open_file("2048.ini","w")
     do i=1, n
        write(fp%unit, *) mygrids(i,:)
     enddo
     call close_file(fp)     
     write(*,*) "data saved"
  else
     goto 100
  endif

contains


  function entropy(f)
    integer f(n,n)
    integer i, j
    real entropy
    entropy = (sum(abs(f(1:n-1, :)-f(2:n,:)))   + &
         sum(abs( f(:, 1:n-1) - f(:, 2:n) )))/(256.*(risky_level+1))
  end function entropy

  recursive function gain_in_direction(fin, direction, nsteps) result(gain)
    integer fin(n,n), fm(n,n), direction, nsteps, dir, gmax
    integer g(4), loc(2), prevloc(2)
    integer  gain
    gain =  merged(fin, fm, direction)
    if(gain .lt. 0)return
    if(nsteps.eq.1)then
       prevloc = maxloc(fin)
       loc = maxloc(fm)
       if(fin(prevloc(1) , prevloc(2)) .lt. fm(loc(1), loc(2)))then
          gain = gain + 0.2   !bonus for getting a larger number
       endif
       gain = gain + ((loc(1)-1)*(loc(1)-n)+(loc(2)-1)*(loc(2)-n)-(prevloc(1)-1)*(prevloc(1)-n)-(prevloc(2)-1)*(prevloc(2)-n))/(real(n)**2) !!punish for moving largest numbers to the center
       gain = gain + entropy(fin) - entropy(fm) !punish for chaos
       return
    endif
    if(nsteps .gt. 5)then  !!see if greedy algorithm works
       do dir = 1, 4
          g(dir) = gain_in_direction(fm, dir, 4)
       enddo
       if(all(g .lt. 0))then
          gain =  gain + (risky_level - 10)*5
          return
       else
          gmax = maxval(g) 
          if(gmax .gt. 5.6)then
             gain = gain + gmax/4.*(nsteps-1)
             return
          endif
       endif
    endif
    do dir = 1, 4
       g(dir) = gain_in_direction(fm, dir, nsteps-1)
    enddo
    if(all(g .lt. 0))then
       gain =  gain + (risky_level - 10)*5
    else
       gain = gain + maxval(g)
    endif
  end function gain_in_direction

  function mcmc_gain(fin, dir, nsteps) result(gain)
    integer,parameter::mcmc_steps = 32
    integer dir, nsteps, fin(n,n), g(0:mcmc_steps)
    integer i
    real gain
    g(0) = gain_in_direction(fin, dir, nsteps)
    if(g(0) .lt. -0.5)then
       gain = g(0)
       return
    endif
    !$omp parallel do
    do i=1, mcmc_steps
       g(i) = gain_in_direction(fin, dir, nsteps)
    enddo
    !$omp end parallel do
    gain = (sum(g)/(mcmc_steps+1.)*risky_level + minval(g)*(10.-risky_level))/10. !!add some weight to the worst case scenario 
  end function mcmc_gain


  function merged(fin, fm, direction)
    integer fin(n,n), fm(n,n), j(n),  direction
    logical success
    integer merged
    merged = 0
    select case(direction)
    case(up)
       call move_up(fin, fm, success, j, merged)
       if(success)call fill_up(fm, j)
    case(down)
       call move_down(fin, fm, success, j, merged)
       if(success)call fill_down(fm, j)
    case(left)
       call move_left(fin, fm, success, j, merged)
       if(success)call fill_left(fm, j)
    case(right)
       call move_right(fin, fm, success, j, merged)
       if(success)call fill_right(fm, j)
    case default
       merged = -9999
    end select
    if(.not. success)then
       merged = -9999
    endif
  end function merged





  subroutine fill_up(fin, j)
    integer fin(n,n), j(n)
    integer i, k, sumj
    real r
    k = random_index(sum(n-j))
    i = 1
    sumj = n - j(1)
    do while(sumj .lt. k)
       i = i + 1
       sumj = sumj + n - j(i)
    enddo
    sumj = n - (sumj - k)
    call random_number(r)
    if(r.lt. prob4)then
       fin(sumj, i) = 4
    else
       fin(sumj, i) = 2
    endif
  end subroutine fill_up

  subroutine fill_down(fin, j)
    integer fin(n,n), j(n)
    integer i, k, sumj
    real r
    k = random_index(sum(n-j))
    i = 1
    sumj = n - j(1)
    do while(sumj .lt. k)
       i = i + 1
       sumj = sumj + n - j(i)
    enddo
    sumj =  (sumj - k) + 1
    call random_number(r)
    if(r.lt. prob4)then
       fin(sumj, i) = 4
    else
       fin(sumj, i) = 2
    endif
  end subroutine fill_down


  subroutine fill_left(fin, j)
    integer fin(n,n), j(n)
    integer i, k, sumj
    real r
    k = random_index(sum(n-j))
    i = 1
    sumj = n - j(1)
    do while(sumj .lt. k)
       i = i + 1
       sumj = sumj + n - j(i)
    enddo
    sumj =  n - (sumj - k)
    call random_number(r)
    if(r.lt. prob4)then
       fin(i, sumj) = 4
    else
       fin(i, sumj) = 2
    endif
  end subroutine fill_left

  subroutine fill_right(fin, j)
    integer fin(n,n), j(n)
    integer i, k, sumj
    real r
    k = random_index(sum(n-j))
    i = 1
    sumj = n - j(1)
    do while(sumj .lt. k)
       i = i + 1
       sumj = sumj + n - j(i)
    enddo
    sumj =  n - (sumj - k)
    call random_number(r)
    if(r.lt. prob4)then
       fin(i, sumj) = 4
    else
       fin(i, sumj) = 2
    endif
  end subroutine fill_right


  subroutine move(fin, direction)
    integer fin(n,n), fmoved(n,n), direction, j(n), nzeros
    logical success
    select case(direction)
    case(up)
       call move_up(fin, fmoved, success, j, nzeros)
    case(down)
       call move_down(fin, fmoved, success, j, nzeros)
    case(left)
       call move_left(fin, fmoved, success, j, nzeros)
    case(right)
       call move_right(fin, fmoved, success, j, nzeros)
    end select
    if(success)then
       fin = fmoved
    else
       stop "Cannot move!"
    endif
  end subroutine move


  subroutine move_up(fin, fmoved, success, j, nzeros)
    logical success
    integer nzeros, i, j(n)
    integer fmoved(n, n), fin(n,n)
    success = .false.
    do i=1, n
       call collapse(fin(:, i), fmoved(:, i), j(i), success, nzeros)
    enddo
  end subroutine move_up


  subroutine move_down(fin, fmoved, success, j, nzeros)
    logical success
    integer nzeros, i, j(n)
    integer fmoved(n, n), fin(n,n)
    success = .false.
    do i=1, n
       call collapse(fin(n:1:-1, i), fmoved(n:1:-1, i), j(i), success, nzeros)
    enddo
  end subroutine move_down

  subroutine move_left(fin, fmoved,  success, j, nzeros)
    logical success
    integer nzeros, i, j(n)
    integer fin(n,n), fmoved(n, n)
    success = .false.
    do i=1, n
       call collapse(fin(i,:), fmoved(i,:), j(i), success, nzeros)
    enddo
  end subroutine move_left


  subroutine move_right(fin, fmoved, success, j,  nzeros)
    logical success
    integer nzeros, i, j(n)
    integer fin(n,n), fmoved(n, n)
    success = .false.
    do i=1, n
       call collapse(fin(i,n:1:-1), fmoved(i,n:1:-1), j(i), success, nzeros)
    enddo
  end subroutine move_right


  function strnum(i)
    integer i
    character(LEN=6) strnum
    if(i.eq.0)then
       strnum = "      "
    else
       strnum = trim(num2str(i))
       strnum = adjustr(strnum)
    endif
  end function strnum

  function strnumarr(i)
    integer i(n)
    character(LEN=6*n) strnumarr
    integer j
    do j=1, n
       strnumarr(j*6-5:j*6) = strnum(i(j))
    enddo
  end function strnumarr

  subroutine show(fin, f2)
    integer fin(n,n)
    integer,optional::f2(n,n)
    integer i, j
    print*
    if(present(f2))then
       do i=1,n
          write(*, "(A)") "|"//strnumarr(fin(i,:))//"|     --->    |"//strnumarr(f2(i,:))//"|"
       enddo
    else
       do i=1,n
          write(*, "(A)") "|"//strnumarr(fin(i,:))//"|"
       enddo
    endif
    print*
  end subroutine show

  subroutine collapse(a, c, j, changed, nzeros) 
    integer,intent(IN):: a(n)
    integer,intent(OUT):: c(n)
    integer  i, inext, j, nzeros
    logical changed
    i = 1
    j = 0
    c = 0
10  if(i .gt. n)then
       return
    elseif(a(i).eq.0)then
       i = i + 1
       goto 10
    endif
20  inext = i+1
30  if(inext .gt. n)then
       j = j + 1
       c(j) = a(i)
       changed = j.ne.i .or. changed
       return
    elseif(a(inext).eq.0)then
       inext = inext + 1
       goto 30
    endif

    if(a(inext).eq.a(i))then
       j = j + 1
       c(j) = a(inext)*2
       changed = j.ne.inext .or. changed
       nzeros = nzeros + 1
       i = inext + 1
       goto 10
    else
       j = j + 1
       c(j) = a(i) 
       changed = j.ne.i .or. changed
       i = inext
       goto 20
    endif
  end subroutine collapse


end program crack2048
