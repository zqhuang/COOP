module coop_weiqi_mod
  use coop_wrapper_typedef
#include "constants.h"
  COOP_INT, parameter::coop_weiqi_max_nmoves = 4096
  COOP_INT, parameter::coop_weiqi_board_size = 19
  COOP_INT, parameter::coop_weiqi_max_score = 999999
  COOP_INT, parameter::coop_weiqi_max_depth = 20
  COOP_INT, parameter::coop_weiqi_max_nstones = coop_weiqi_board_size**2
  COOP_INT, parameter::coop_weiqi_max_ngroups = coop_weiqi_max_nstones
  COOP_INT, parameter::coop_weiqi_black = 1
  COOP_INT, parameter::coop_weiqi_white = -1
  COOP_INT, parameter::coop_weiqi_pass = -1  
  COOP_INT, parameter::coop_weiqi_capacity_increment = 16
  character,parameter::coop_weiqi_symbol(-1:1) = (/ "w", "o", "b" /)

  type coop_weiqi_list
     COOP_INT::n = 0
     COOP_INT::capacity = 0
     COOP_INT,dimension(:), allocatable::l
   contains
     procedure::expand => coop_weiqi_list_expand
     procedure::insert => coop_weiqi_list_insert
     procedure::check_insert => coop_weiqi_list_check_insert     
     procedure::delete => coop_weiqi_list_delete
     procedure::is_element => coop_weiqi_list_is_element
     procedure::free => coop_weiqi_list_free
  end type coop_weiqi_list
  
  
  type coop_weiqi_board
     COOP_INT,dimension(0:coop_weiqi_max_nstones-1)::stones, prevstones, nn, gid
     COOP_REAL,dimension(0:coop_weiqi_max_nstones-1):: prior
     COOP_INT,dimension(4, 0:coop_weiqi_max_nstones-1)::nbs     
     COOP_INT::turn = coop_weiqi_black, prevmove = -1
     COOP_INT::maxid(-1:1) = 0     
     logical::gused(-coop_weiqi_max_ngroups:coop_weiqi_max_ngroups) = .false.
     type(coop_weiqi_list),dimension(-coop_weiqi_max_ngroups:coop_weiqi_max_ngroups)::groups, libs
     type(coop_weiqi_list)::init_stones(-1:1)
     COOP_INT::tolive, tokill
     logical::check_ko = .true.
   contains
     procedure::init => coop_weiqi_board_init
     procedure::spread => coop_weiqi_board_spread
     procedure::move => coop_weiqi_board_move
     procedure::quick_move => coop_weiqi_board_quick_move     
     procedure::is_allowed_move => coop_weiqi_board_is_allowed_move     
     procedure::merge_groups => coop_weiqi_board_merge_groups
     procedure::print => coop_weiqi_board_print
     procedure::group_add_stone => coop_weiqi_board_group_add_stone
     procedure::is_alive => coop_weiqi_board_is_alive
     procedure::score => coop_weiqi_board_score
     procedure::get_prior => coop_weiqi_board_get_prior
  end type coop_weiqi_board


  
  
contains

  subroutine coop_weiqi_point_to_coor(pt, x, y)
    COOP_INT::pt, x, y
    x = pt/coop_weiqi_board_size 
    y = pt - x*coop_weiqi_board_size
  end subroutine coop_weiqi_point_to_coor

  subroutine coop_weiqi_coor_to_point(x,y, pt)
    COOP_INT::pt, x, y
    pt = y + x*coop_weiqi_board_size
  end subroutine coop_weiqi_coor_to_point
  

  subroutine coop_weiqi_neighbors(pt, nbs, nn)
    COOP_INT::pt, nbs(:), nn, x, y
    call coop_weiqi_point_to_coor(pt, x, y)
    if(x .eq. 0)then
       if(y .eq. 0)then
          nn = 2
          call coop_weiqi_coor_to_point(0, 1, nbs(1))
          call coop_weiqi_coor_to_point(1, 0, nbs(2))          
       elseif(y .eq. coop_weiqi_board_size-1)then
          nn = 2
          call coop_weiqi_coor_to_point(0, coop_weiqi_board_size-2, nbs(1))
          call coop_weiqi_coor_to_point(1, coop_weiqi_board_size-1, nbs(2))          
       else
          nn = 3
          call coop_weiqi_coor_to_point(0, y-1, nbs(1))
          call coop_weiqi_coor_to_point(0, y+1, nbs(2))          
          call coop_weiqi_coor_to_point(1, y, nbs(3))
       endif
    elseif(x.eq.coop_weiqi_board_size-1)then
       if(y .eq. 0)then
          nn = 2
          call coop_weiqi_coor_to_point(coop_weiqi_board_size-2, 0, nbs(1))                    
          call coop_weiqi_coor_to_point(coop_weiqi_board_size-1, 1, nbs(2))
       elseif(y .eq. coop_weiqi_board_size-1)then
          nn = 2
          call coop_weiqi_coor_to_point(coop_weiqi_board_size-2, coop_weiqi_board_size-1, nbs(1))                    
          call coop_weiqi_coor_to_point(coop_weiqi_board_size-1, coop_weiqi_board_size-2, nbs(2))
       else
          nn = 3
          call coop_weiqi_coor_to_point(coop_weiqi_board_size-2, y, nbs(1))
          call coop_weiqi_coor_to_point(coop_weiqi_board_size-1, y-1, nbs(2))
          call coop_weiqi_coor_to_point(coop_weiqi_board_size-1, y+1, nbs(3))          
       endif       
    else
       if(y .eq. 0)then
          nn = 3
          call coop_weiqi_coor_to_point(x-1, 0, nbs(1))
          call coop_weiqi_coor_to_point(x, 1, nbs(2))
          call coop_weiqi_coor_to_point(x+1, 0, nbs(3))
       elseif(y .eq. coop_weiqi_board_size-1)then
          nn = 3
          call coop_weiqi_coor_to_point(x-1, coop_weiqi_board_size-1, nbs(1))
          call coop_weiqi_coor_to_point(x, coop_weiqi_board_size-2, nbs(2))
          call coop_weiqi_coor_to_point(x+1, coop_weiqi_board_size-1, nbs(3))
       else
          nn = 4
          call coop_weiqi_coor_to_point(x-1, y, nbs(1))
          call coop_weiqi_coor_to_point(x, y-1, nbs(2))
          call coop_weiqi_coor_to_point(x, y+1, nbs(3))                   
          call coop_weiqi_coor_to_point(x+1, y, nbs(4))                   
       endif       
    endif
  end subroutine coop_weiqi_neighbors

  subroutine coop_weiqi_list_expand(this)
    class(coop_weiqi_list)::this
    COOP_INT,dimension(:),allocatable::copy
    if(allocated(this%l))then
       allocate(copy(this%n))
       copy = this%l(1:this%n)
       deallocate(this%l)
       this%capacity = this%capacity + coop_weiqi_capacity_increment
       allocate(this%l(this%capacity))
       this%l(1:this%n) = copy
       deallocate(copy)
    else
       this%capacity = coop_weiqi_capacity_increment
       allocate(this%l(this%capacity))
    endif
  end subroutine coop_weiqi_list_expand

  subroutine coop_weiqi_list_insert(this, pt)
    class(coop_weiqi_list)::this
    COOP_INT::pt
    if(this%n .ge. this%capacity) call this%expand()
    this%n = this%n + 1
    this%l(this%n) = pt
  end subroutine coop_weiqi_list_insert

  subroutine coop_weiqi_list_check_insert(this, pt)
    class(coop_weiqi_list)::this
    COOP_INT::pt
    if(this%is_element(pt))return
    if(this%n .ge. this%capacity) call this%expand()
    this%n = this%n + 1
    this%l(this%n) = pt
  end subroutine coop_weiqi_list_check_insert
  

  !!I do not do sorted list for two reasons
  !!1. the groups are usually small, typically 2-10 stones
  !!2. less operation on memory copying  
  subroutine coop_weiqi_list_delete(this, pt)
    class(coop_weiqi_list)::this
    COOP_INT,intent(IN)::pt
    COOP_INT::i
    do i = 1, this%n
       if(this%l(i).eq.pt)then
          this%l(i) = this%l(this%n)
          this%n = this%n - 1 
          exit
       endif
    enddo
  end subroutine coop_weiqi_list_delete

  function coop_weiqi_list_is_element(this, pt) result(isele)
    class(coop_weiqi_list)::this
    COOP_INT,intent(IN)::pt
    logical:: isele
    isele  = any(this%l(1:this%n) .eq. pt)
  end function coop_weiqi_list_is_element
  
  subroutine coop_weiqi_list_free(this)
    class(coop_weiqi_list)::this
    this%n = 0
    this%capacity = 0
    if(allocated(this%l))deallocate(this%l)
  end subroutine coop_weiqi_list_free


  subroutine coop_weiqi_board_init(this, blist, wlist, turn, tolive)
    class(coop_weiqi_board)::this
    COOP_INT::i, j, k
    COOP_INT,dimension(:),optional::blist, wlist
    COOP_INT, optional::turn, tolive
    COOP_INT::id
    do i = -this%maxid(-1), this%maxid(1)
       call this%groups(i)%free()
       call this%libs(i)%free()
    enddo
    call this%init_stones(1)%free()
    call this%init_stones(-1)%free()
    this%gused = .false.    
    this%stones = 0
    this%gid = 0
    !$omp parallel do
    do i=0, coop_weiqi_max_nstones-1
       call coop_weiqi_neighbors(i, this%nbs(:,i), this%nn(i))
    enddo
    !$omp end parallel do
    if(present(turn))then
       this%turn = turn
    else
       this%turn = coop_weiqi_black
    endif
    if(present(tolive))then
       this%tolive = tolive
    else
       this%tolive = coop_weiqi_white
    endif
    this%tokill = - this%tolive
    if(present(blist))then
       this%stones(blist) = coop_weiqi_black
       id = 0
       do i = 1, size(blist)
          call this%init_stones(coop_weiqi_black)%insert(blist(i))
          if(this%gid(blist(i)) .eq. 0)then !!start a new id
             id = id + coop_weiqi_black
             this%gid(blist(i)) = id
             call this%spread(blist(i))
          endif
       enddo
       this%maxid(1) = id
    endif
    if(present(wlist))then       
       this%stones(wlist) = coop_weiqi_white
       id = 0
       do i = 1, size(wlist)
          call this%init_stones(coop_weiqi_white)%insert(wlist(i))
          if(this%gid(wlist(i)) .eq. 0)then !!start a new id
             id = id +  coop_weiqi_white
             this%gid(wlist(i)) = id
             call this%spread(wlist(i))
          endif
       enddo
       this%maxid(-1) = id
    endif
    do i = 0, coop_weiqi_max_nstones-1
       if(this%gid(i).ne.0)then
          call this%group_add_stone(this%gid(i), i)
       endif
    enddo
    this%prevstones = 0
    this%prevmove = -1
  end subroutine coop_weiqi_board_init

  recursive subroutine coop_weiqi_board_spread(this, pt)
    class(coop_weiqi_board)::this
    COOP_INT,intent(IN)::pt
    COOP_INT::i, color, id
    color = this%stones(pt)
    id = this%gid(pt)
    do i = 1, this%nn(pt)
       if( this%stones(this%nbs(i, pt)).eq. color .and. this%gid(this%nbs(i, pt)) .ne. id ) then
          this%gid(this%nbs(i, pt)) = id
          call coop_weiqi_board_spread(this, this%nbs(i, pt))
       endif
    enddo
  end subroutine coop_weiqi_board_spread


  function coop_weiqi_board_is_allowed_move(this, pt) result(valid)
    class(coop_weiqi_board)::this
    COOP_INT,dimension(0:coop_weiqi_max_nstones-1)::stones
    COOP_INT,intent(IN)::pt
    COOP_INT::i, id
    logical::invalid
    logical::valid
    !!put on the stone
    if(pt .eq. coop_weiqi_pass)then
       valid = .true.
       return
    endif
    if(this%stones(pt) .ne. 0)then
       valid = .false.
       return
    endif
    invalid = all(this%stones(this%nbs(1:this%nn(pt), pt)).ne. 0)
    if(invalid)then
       do i = 1, this%nn(pt)
          id = this%gid(this%nbs(i, pt))       
          if(this%stones(this%nbs(i, pt)).eq.this%turn .and. this%libs(id)%n.gt.1)then  
             invalid = .false.
             exit
          elseif( this%stones(this%nbs(i, pt)).eq.-this%turn .and. this%libs(id)%n .le. 1 )then 
             invalid = .false.
             exit
          endif
       enddo
    endif
    if(invalid)then
       valid = .false.
       return
    endif
    if(this%check_ko .and. this%prevmove .ge. 0)then
       call this%quick_move(pt, stones)
       valid = any(stones .ne. this%prevstones)
    else
       valid = .true.
       return
    endif
  end function coop_weiqi_board_is_allowed_move

  subroutine coop_weiqi_board_move(this, move, valid)
    class(coop_weiqi_board)::this
    COOP_INT::move
    COOP_INT::pt
    COOP_INT::i, id, j, k, ss
    COOP_INT::lib_added(4)
    logical, optional::valid
    pt = move
    if(pt .eq. coop_weiqi_pass)then
       this%turn = - this%turn
       if(present(valid))then
          valid = .true.
       endif
       return
    endif
    if(this%is_allowed_move(pt))then
       if(present(valid))then
          valid = .true.
       endif
    else
       if(present(valid))then
          valid = .false.
          return          
       else
          call this%print()
          write(*,"(A5, I5, A15)") "move:", pt, " is invalid"
          stop
       endif
    endif
    this%prevmove = pt
    this%prevstones = this%stones

    !!now change stuff
    this%stones(pt) = this%turn

    !!update group id     
    this%gid(pt) = 0 !!null list
    if(any(this%stones(this%nbs(1:this%nn(pt), pt)).eq. this%turn))then  !!combine groups
       do i = 1, this%nn(pt)
          if(this%stones(this%nbs(i, pt)).eq. this%turn)then
             if(this%groups(this%gid(this%nbs(i, pt)))%n .gt. this%groups(this%gid(pt))%n)then
                this%gid(pt) = this%gid(this%nbs(i, pt))
             endif
          endif
       enddo
       call this%group_add_stone(this%gid(pt), pt)
       do i = 1, this%nn(pt)
          if(this%stones(this%nbs(i, pt)).eq. this%turn)then  !!merge into this gourp
             id = this%gid(this%nbs(i, pt)) 
             call this%merge_groups(this%gid(pt), id, pt)
          endif
       enddo
    else  !!an independent group
       do i = this%turn, this%maxid(this%turn)+this%turn, this%turn

          if(.not.this%gused(i))then
             call this%group_add_stone(i, pt)
             if(abs(i) .gt. abs(this%maxid(this%turn))) &
                  this%maxid(this%turn) = i
             exit
          endif
       enddo
    endif
    !!remove dead stones
    do i = 1, this%nn(pt)
       if(this%stones(this%nbs(i, pt)).eq. -this%turn)then
          id = this%gid(this%nbs(i, pt))
          call this%libs(id)%delete(pt)
          if(this%libs(id)%n .eq. 0)then
             do j = 1, this%groups(id)%n
                ss = this%groups(id)%l(j)
                lib_added = 0
                do k = 1, this%nn(ss)
                   if(this%stones(this%nbs(k, ss)) .eq. this%turn)then  !!if there is a black group around                      
                      if(k.eq.1 .or. (.not. any(lib_added(1:k-1).eq. this%gid(this%nbs(k, ss)))))then
                         call this%libs( this%gid(this%nbs(k, ss) ) )%check_insert( ss )
                         lib_added(k) = this%gid(this%nbs(k, ss))                         
                      endif
                   endif

                   this%gid(ss) = 0
                   this%stones(ss) = 0 !!delte the dead stone
                   
                enddo
             enddo
             this%gused(id) = .false.
             this%groups(id)%n = 0
          endif
       endif
    enddo
    if(present(valid))then
       valid = (this%libs(this%gid(pt))%n .gt. 0)
       if(.not.valid)then
          call this%print(want_groups = .true.)
          write(*,"(A5, I5, A15)") "move:", pt, " is invalid"
          stop "This should never happen!"
       endif
    endif
    !!change turn
    this%turn = - this%turn
  end subroutine coop_weiqi_board_move


  subroutine coop_weiqi_board_quick_move(this, move, stones)
    class(coop_weiqi_board)::this
    COOP_INT, dimension(0:coop_weiqi_max_nstones-1)::stones
    COOP_INT::move
    COOP_INT::pt
    COOP_INT::i, id
    pt = move
    stones = this%stones
    stones(pt) = this%turn
    !!remove dead stones
    do i = 1, this%nn(pt)
       if(stones(this%nbs(i, pt)).eq. -this%turn)then
          id = this%gid(this%nbs(i, pt))
          if(this%libs(id)%n .le. 1) &
               stones(this%groups(id)%l(1:this%groups(id)%n)) = 0 
       endif
    enddo
  end subroutine coop_weiqi_board_quick_move
  


  subroutine coop_weiqi_board_print(this, want_groups, want_prior)
    class(coop_weiqi_board)::this
    logical,optional::want_groups, want_prior
    COOP_INT::x,y, pt(0:coop_weiqi_board_size - 1)
    write(*,"(A)") "-----------------------------------------"    
    write(*,"(A)") "stones:"    
    do y=0, coop_weiqi_board_size - 1  
       do x = 0, coop_weiqi_board_size - 1
          call coop_weiqi_coor_to_point(x, y, pt(x))
       enddo
       write(*, "(19A4)") coop_weiqi_symbol(this%stones(pt))
    enddo
    if(present(want_groups))then
       if(want_groups)then
          write(*,"(A)") "groups:"
          do y=0, coop_weiqi_board_size - 1  
             do x = 0, coop_weiqi_board_size - 1
                call coop_weiqi_coor_to_point(x, y, pt(x))
             enddo
             write(*, "(19I4)") this%gid(pt)
          enddo
       endif
    endif
    if(present(want_prior))then
       if(want_prior)then
          call this%get_prior()    
          write(*,"(A)") "prior:"
          do y=0, coop_weiqi_board_size - 1  
             do x = 0, coop_weiqi_board_size - 1
                call coop_weiqi_coor_to_point(x, y, pt(x))
             enddo
             write(*, "(19F6.2)") this%prior(pt)
          enddo
       endif
    endif
    write(*,"(A)") "-----------------------------------------"        
  end subroutine coop_weiqi_board_print


  subroutine coop_weiqi_board_merge_groups(this, id1, id2, connect_stone)
    class(coop_weiqi_board)::this
    COOP_INT::id1, id2, i, j, connect_stone
    if(id2 .eq. id1) return
    this%gid(this%groups(id2)%l(1:this%groups(id2)%n)) = id1
    do i = 1, this%groups(id2)%n
       call this%groups(id1)%insert(this%groups(id2)%l(i))
    enddo
    do i = 1, this%libs(id2)%n
       if(this%libs(id2)%l(i) .ne. connect_stone) call this%libs(id2)%check_insert(this%libs(id2)%l(i))
    enddo
    this%gused(id2) = .false.
    this%groups(id2)%n = 0
    this%libs(id2)%n = 0
    this%gused(id1) = this%groups(id1)%n .gt. 0
  end subroutine coop_weiqi_board_merge_groups

  subroutine coop_weiqi_board_group_add_stone(this, id, stone)
    class(coop_weiqi_board)::this
    COOP_INT::id, stone, i
    call this%groups(id)%insert(stone)
    call this%libs(id)%delete(stone)
    do i = 1, this%nn(stone)
       if(this%stones(this%nbs(i, stone)).eq.0) &
            call this%libs(id)%check_insert(this%nbs(i, stone))
    enddo
    this%gid(stone) = id
    this%gused(id) = this%groups(id)%n.gt.0
  end subroutine coop_weiqi_board_group_add_stone


  function coop_weiqi_board_score(this) result(score)
    class(coop_weiqi_board)::this
    COOP_REAL::score
    COOP_INT::id, j
    if( all(this%stones(this%init_stones(this%tolive)%l(1:this%init_stones(this%tolive)%n)) .ne. this%tolive)) then !!killed
       score = coop_weiqi_max_score * this%tokill
    elseif(this%is_alive(this%tolive))then  !!live
       score = coop_weiqi_max_score * this%tolive
    else  !!not sure, need to keep playing
       score = 0
       do id =  this%tolive, this%maxid(this%tolive), this%tolive
          if(this%gused(id))then
             score = score + this%libs(id)%n - 2
          endif
       enddo
       score = score * this%tolive
    endif
  end function coop_weiqi_board_score

  subroutine coop_weiqi_board_get_prior(this)
    class(coop_weiqi_board)::this
    type(coop_weiqi_board)::copy
    logical valid
    COOP_INT::id, i, j, j1, k1, j2, k2, j3, k3, j4, k4, x, y, xx, yy
    COOP_REAL::p(0: coop_weiqi_max_nstones-1, -1:1)
    p = 0.d0
    do id = this%maxid(-1), -1
       if(this%gused(id))then
          do j1 = 1, this%libs(id)%n
             k1 = this%libs(id)%l(j1)
             p(k1, -1) = p(k1, 1) + 0.6d0
             do j2 = 1, this%nn(k1)
                k2 = this%nbs(j2, k1)
                if(this%stones(k2).eq. 0)then
                   p(k2, -1) = p(k2, 1) + 0.4d0
                   do j3 = 1, this%nn(k2)
                      k3 = this%nbs(j3, k2)
                      if(this%stones(k3).eq.0)p(k3, -1) = p(k3, -1) + 0.3d0
                      do j4 = 1, this%nn(k3)
                         k4 = this%nbs(j4, k3)
                         if(this%stones(k4).eq.0)p(k4, -1) = p(k4, -1) + 0.2d0
                      enddo
                   enddo
                endif
             enddo
          enddo
       endif       
    enddo
    do id=1, this%maxid(1)
       if(this%gused(id))then
          do j1 = 1, this%libs(id)%n
             k1 = this%libs(id)%l(j1)
             p(k1, 1) = p(k1, 1) + 0.6d0
             do j2 = 1, this%nn(k1)
                k2 = this%nbs(j2, k1)
                if(this%stones(k2).eq.0)then
                   p(k2, 1) = p(k2, 1) + 0.4d0
                   do j3 = 1, this%nn(k2)
                      k3 = this%nbs(j3, k2)
                      if(this%stones(k3).eq.0) p(k3, 1) = p(k3, 1) + 0.3d0
                      do j4 = 1, this%nn(k3)
                         k4 = this%nbs(j4, k3)
                         if(this%stones(k4).eq.0)p(k4, 1) = p(k4, 1) + 0.2d0
                      enddo
                   enddo
                endif
             enddo
          enddo
       endif       
    enddo
    this%prior = p(:,1)*p(:,-1)
    do i = 0, coop_weiqi_max_nstones-1
       if(this%stones(i) .eq. 0 .and. this%prior(i).gt. 0.d0)then
          if(this%is_allowed_move(i))then
             call coop_weiqi_point_to_coor(i, x, y)       
             xx = coop_weiqi_board_size/2 - abs(x - coop_weiqi_board_size/2)
             yy = coop_weiqi_board_size/2 - abs(y - coop_weiqi_board_size/2)
             this%prior(i) = this%prior(i)*(1.+ exp(- (min(xx, yy) + 1.)**2-0.5d0) + exp(-(xx**2+yy**2)/32.d0-0.5d0))
          else
             this%prior(i) = 0.d0
          endif
       endif
    enddo

    this%prior = this%prior/max(maxval(this%prior), 1.d-99)
  end subroutine coop_weiqi_board_get_prior

  function coop_weiqi_board_is_alive(this, color) result(alive)
    class(coop_weiqi_board)::this
    COOP_INT::color, id, i, j
    logical::alive, valid, anyvalid
    type(coop_weiqi_board)::copy
    select type(this)
    type is(coop_weiqi_board)
       copy = this
       copy%check_ko = .false.       
    class default
       stop "coop_weiqi: class copy failed"       
    end select

100 anyvalid = .false.
    do id = color, copy%maxid(color), color
       if(copy%gused(id))then
          do i=1, copy%libs(id)%n
             copy%turn = - color
             call copy%move(copy%libs(id)%l(i), valid)
             anyvalid = valid .or. anyvalid
          enddo
       endif
    enddo
    if(anyvalid) goto 100
    alive = any(copy%gused( color:copy%maxid(color):color ))
  end function coop_weiqi_board_is_alive


end module coop_weiqi_mod
