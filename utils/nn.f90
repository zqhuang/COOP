Module coop_nn_mod
  use coop_random_mod
  use coop_list_mod
  use coop_matrix_mod
  use coop_file_mod
  use coop_wrapper_typedef
  use coop_MPI_mod
  use coop_cholesky_mod
  implicit none
#include "constants.h"  

#define ACTIVATE_FUNCTION 3
  !! 0 = sigmoid
  !! 1 = ReLU (Rectified Linear Unit)
  !! 2 = tanh
  !! 3 = modified ReLU
 

  type coop_nn_connection
     COOP_INT:: in = 0
     COOP_INT:: out = 0
     COOP_REAL :: weight = 0.d0
     COOP_REAL :: best_weight = 0.d0     
     COOP_REAL :: dw = 0.d0
  end type coop_nn_connection

  type coop_nn_layer
     COOP_INT:: nin =0
     COOP_INT:: nc=0
     COOP_INT:: nout = 0
     COOP_REAL,dimension(:),allocatable::v
     COOP_REAL,dimension(:),allocatable::dv
     type(coop_nn_connection),dimension(:),allocatable::c
     COOP_REAL,dimension(:),allocatable::bias
     COOP_REAL,dimension(:),allocatable::best_bias     
     COOP_REAL,dimension(:),allocatable::db
   contains
     procedure::init=>coop_nn_layer_init
     procedure::full_init=>coop_nn_layer_full_init     
     procedure::free=>coop_nn_layer_free     
     procedure::propogate => coop_nn_layer_propogate
     procedure::chain_derv => coop_nn_layer_chain_derv
     procedure::init_connections => coop_nn_layer_init_connections
  end type coop_nn_layer

  type coop_nn
     COOP_INT::nlayers, nin, nout
     COOP_REAL::best_Err = 1.d99     
     type(coop_nn_layer),dimension(:),allocatable::layers
     COOP_REAL,dimension(:),allocatable::true_out
   contains
     procedure::init=>coop_nn_init
     procedure::full_init=>coop_nn_full_init     
     procedure::free=>coop_nn_free
     procedure::fp => coop_nn_fp
     procedure::bp => coop_nn_bp
     procedure::rand => coop_nn_rand
     procedure::Err => coop_nn_Err
     procedure::walk => coop_nn_walk
  end type coop_nn

contains


  subroutine  coop_nn_layer_init_connections(this, in, out)
    class(coop_nn_layer)::this
    COOP_INT,dimension(:),optional::in, out
    COOP_INT::i, j
    if(present(in) .and. present(out))then
       if(size(in) .ne. this%nc .or. size(out) .ne. this%nc) call coop_return_error("coop_nn_layer_init_connections","# of connections do not match", "stop")

       do i = 1, this%nc
          this%c(i)%in = in(i)
          this%c(i)%out = out(i)
       enddo
    else
       if(this%nc .ne. this%nin*this%nout) call coop_return_error("coop_nn_layer_init_connections", "# of connections do not match", "stop")
       do i=0, this%nin-1
          do j=0, this%nout-1
             this%c(i*this%nout + j + 1)%in = i+1
             this%c(i*this%nout + j + 1)%out = j+1 
          enddo
       enddo
    endif
  end subroutine coop_nn_layer_init_connections

  subroutine coop_nn_layer_free(this)
    class(coop_nn_layer)::this
    COOP_DEALLOC(this%v)
    COOP_DEALLOC(this%c)
    COOP_DEALLOC(this%bias)
    COOP_DEALLOC(this%best_bias)    
    COOP_DEALLOC(this%db)
!!$  COOP_DEALLOC(this%last_db)    
    COOP_DEALLOC(this%dv)    
    this%nin = 0
    this%nc = 0
    this%nout = 0
  end subroutine coop_nn_layer_free



  subroutine coop_nn_layer_init(this, nin, nout, nc)
    class(coop_nn_layer)::this
    COOP_INT::nin, nout, nc
    if(this%nc .ne. nc)then
       COOP_DEALLOC(this%c)
       this%nc = nc
       if(nc.gt.0) allocate(this%c(nc))
    endif
    if(this%nin .ne. nin)then
       COOP_DEALLOC(this%v)
       COOP_DEALLOC(this%dv)
       this%nin = nin
       allocate(this%v(nin))
       allocate(this%dv(nin))
    endif
    if(this%nout .ne. nout)then
       COOP_DEALLOC(this%bias)
       COOP_DEALLOC(this%best_bias)       
       COOP_DEALLOC(this%db)
!!$       COOP_DEALLOC(this%last_db)       
       this%nout = nout
       if(nout .gt. 0)then
          allocate(this%bias(nout))
          allocate(this%best_bias(nout))          
          allocate(this%db(nout))
!!$          allocate(this%last_db(nout))          
       endif
    endif
  end subroutine coop_nn_layer_init



  subroutine coop_nn_layer_full_init(this, nin, nout)
    class(coop_nn_layer)::this
    COOP_INT::nin, nout
    if(this%nin .ne. nin)then
       COOP_DEALLOC(this%v)
       COOP_DEALLOC(this%dv)
       this%nin = nin
       allocate(this%v(nin))
       allocate(this%dv(nin))
    endif
    if(this%nout .ne. nout)then
       COOP_DEALLOC(this%bias)
       COOP_DEALLOC(this%best_bias)       
       COOP_DEALLOC(this%db)
!!$       COOP_DEALLOC(this%last_db)       
       this%nout = nout
       if(nout .gt. 0)then
          allocate(this%bias(nout))
          allocate(this%best_bias(nout))          
          allocate(this%db(nout))
!!$          allocate(this%last_db(nout))          
       endif
    endif
    this%nc = this%nin*this%nout
    if(this%nc .gt. 0)then
       allocate(this%c(this%nc))
       call this%init_connections()
    endif
  end subroutine coop_nn_layer_full_init


  subroutine coop_nn_layer_propogate(this, next)
    class(coop_nn_layer)::this
    class(coop_nn_layer)::next
    COOP_INT::i
    if(next%nin .ne. this%nout) call coop_return_error("coop_nn_layer_propogate","layers do not fit", "stop") !!quick check
    next%v(1:this%nout) = this%bias(1:this%nout)
    do i=1, this%nc
       next%v(this%c(i)%out) = next%v(this%c(i)%out) + this%v(this%c(i)%in)*this%c(i)%weight
    enddo
    do i=1, this%nout
       next%v(i) = coop_nn_activate(next%v(i))
    enddo
  end subroutine coop_nn_layer_propogate


  subroutine coop_nn_layer_chain_derv(this, next)
    class(coop_nn_layer)::this
    class(coop_nn_layer)::next
    COOP_INT::i, j
!!$    this%last_db = this%db
!!$    this%c%last_dw = this%c%dw    
    do i=1, this%nout
       this%db(i) = next%dv(i)*coop_nn_activate_derv(next%v(i))
    enddo
    this%dv = 0.d0
    do i=1, this%nc
       this%dv(this%c(i)%in) = this%dv(this%c(i)%in) + next%dv(this%c(i)%out)  *  coop_nn_activate_derv(next%v(this%c(i)%out)) * this%c(i)%weight
       this%c(i)%dw = next%dv(this%c(i)%out) * coop_nn_activate_derv(next%v(this%c(i)%out)) * this%v(this%c(i)%in)
    enddo
  end subroutine coop_nn_layer_chain_derv


  subroutine coop_nn_rand(this)
    class(coop_nn)::this
    COOP_INT::i,j
    do i=1, this%nlayers
       call random_number(this%layers(i)%bias)
       do j=1, this%layers(i)%nc
          call random_number(this%layers(i)%c(j)%weight)
       enddo
    enddo
  end subroutine coop_nn_rand

  
  subroutine coop_nn_fp(this)
    class(coop_nn)::this
    COOP_INT::i
    do i=2, this%nlayers
       call this%layers(i-1)%propogate(this%layers(i))
    enddo
  end subroutine coop_nn_fp

  subroutine coop_nn_bp(this)
    class(coop_nn)::this
    COOP_INT::i, j
    this%layers(this%nlayers)%dv = this%layers(this%nlayers)%v - this%true_out
    do i=this%nlayers-1, 1, -1
       call this%layers(i)%chain_derv(this%layers(i+1))
    enddo
  end subroutine coop_nn_bp
  
  subroutine coop_nn_free(this)
    class(coop_nn)::this
    COOP_INT::i
    do i=1, this%nlayers
       call this%layers(i)%free()
    enddo
    COOP_DEALLOC(this%layers)
    this%nlayers = 0
    this%nin = 0
    this%nout = 0
  end subroutine coop_nn_free

  function coop_nn_Err(this) result(Err)
    class(coop_nn)::this    
    COOP_REAL::Err
    Err = sum((this%true_out - this%layers(this%nlayers)%v)**2)/2.d0
  end function coop_nn_Err


  subroutine coop_nn_Walk(this, step)
    class(coop_nn)::this
    COOP_REAL::step, s, sumd2, err
    COOP_INT::i, j
    call this%fp()
    err = this%Err()
    if(err .lt. this%best_Err)then
       this%best_Err = err
       do i = 1, this%nlayers-1
          this%layers(i)%best_bias = this%layers(i)%bias
          this%layers(i)%c%best_weight = this%layers(i)%c%weight
       enddo
    elseif(exp(this%best_Err-err)  .lt. coop_random_unit() )then
       do i = 1, this%nlayers-1
          this%layers(i)%bias = this%layers(i)%best_bias
          this%layers(i)%c%weight = this%layers(i)%c%best_weight
       enddo
       call this%fp()
    endif
    call this%bp()

    sumd2 = 0.d0
    do i = 1, this%nlayers-1
       sumd2 = sumd2 + sum(this%layers(i)%db**2) + sum(this%layers(i)%c%dw**2)
    enddo
    
    s = step*err/max(sumd2, 1.d-99)
    
    do i = 1, this%nlayers-1
       this%layers(i)%bias = this%layers(i)%bias - s* coop_random_Gaussian() * this%layers(i)%db 
       this%layers(i)%c%weight = this%layers(i)%c%weight - s* coop_random_Gaussian() * this%layers(i)%c%dw 
    enddo
    
  end subroutine coop_nn_Walk
  
  
  subroutine coop_nn_init(this, nlayers, nin, nout)
    class(coop_nn)::this
    COOP_INT::nlayers, nin, nout
    call this%free()
    this%nlayers = nlayers
    allocate(this%layers(nlayers))
    this%nin = nin
    this%nout = nout
    allocate(this%true_out(nout))    
  end subroutine coop_nn_init


  subroutine coop_nn_full_init(this, n)
    class(coop_nn)::this
    COOP_INT,dimension(:)::n
    COOP_INT::i
    call this%free()
    this%nlayers = size(n)
    allocate(this%layers(this%nlayers))
    this%nin = n(1)
    this%nout = n(this%nlayers)
    allocate(this%true_out(this%nout))
    do i=1, this%nlayers-1
       call this%layers(i)%full_init(n(i), n(i+1))
    enddo
    call this%layers(this%nlayers)%full_init(n(this%nlayers), 0)
    call this%rand()
  end subroutine coop_nn_full_init
  
  

#if ACTIVATE_FUNCTION == 0 
  function coop_nn_activate(x) result(f)
    COOP_REAL::x, f
    f = 1.d0/(1.d0+exp(-x))
  end function coop_nn_activate

  function coop_nn_activate_derv(f) result(df)
    COOP_REAL::f,df
    df  = f*(1.d0-f)
  end function coop_nn_activate_derv

#elif ACTIVATE_FUNCTION == 1

  function coop_nn_activate(x) result(f)
    COOP_REAL::x, f
    f = max(0.d0, x)
  end function coop_nn_activate

  function coop_nn_activate_derv(f) result(df)
    COOP_REAL::f,df
    if(f.gt.0.d0)then
       df = 1.d0
    else
       df = 0.d0
    endif
  end function coop_nn_activate_derv

#elif ACTIVATE_FUNCTION == 2

  function coop_nn_activate(x) result(f)
    COOP_REAL::x, f
    f = tanh(x)
  end function coop_nn_activate

  function coop_nn_activate_derv(f) result(df)
    COOP_REAL::f,df
    df = 1.d0-f**2
  end function coop_nn_activate_derv

#elif ACTIVATE_FUNCTION == 3

  function coop_nn_activate(x) result(f)
    COOP_REAL::x, f
    if(x.gt.0.d0)then
       f = x
    else
       f = x/10.d0
    endif
  end function coop_nn_activate

  function coop_nn_activate_derv(f) result(df)
    COOP_REAL::f,df
    if(f.gt.0.d0)then
       df = 1.d0
    else
       df = 0.1d0
    endif
  end function coop_nn_activate_derv
  
#endif  

  !!==============================


  
end module coop_nn_mod
    

