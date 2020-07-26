module coop_hnn_mod
  use coop_healpix_mod
  use coop_wrapper_firstorder
  implicit none
  
#include "constants.h"  

#define ACTIVATE_FUNCTION 2
  !! 0 = sigmoid
  !! 1 = ReLU (Rectified Linear Unit)
  !! 2 = tanh
  !! 3 = modified ReLU
  !!all maps are converted to nested ordering

  
  interface coop_hnn_activate
     module procedure coop_hnn_activate_s, coop_hnn_activate_v, coop_hnn_activate_vv
  end interface coop_hnn_activate

  interface coop_hnn_activate_derv
     module procedure coop_hnn_activate_derv_s, coop_hnn_activate_derv_v, coop_hnn_activate_derv_vv
  end interface coop_hnn_activate_derv
  
  type coop_hnn_connection
     COOP_INT:: in = 0
     COOP_INT:: out = 0
     COOP_REAL :: weight = 0.d0
     COOP_REAL :: best_weight = 0.d0     
     COOP_REAL :: dw = 0.d0
  end type coop_hnn_connection

  type coop_hnn_layer
     COOP_INT:: nin =0
     COOP_INT:: nc=0
     COOP_INT:: nout = 0
     COOP_INT:: nside = 0
     type(coop_healpix_maps)::v, dv
     type(coop_hnn_connection),dimension(:),allocatable::c
     COOP_REAL,dimension(:),allocatable::bias
     COOP_REAL,dimension(:),allocatable::best_bias     
     COOP_REAL,dimension(:),allocatable::db
   contains
     procedure::init=>coop_hnn_layer_init
     procedure::free=>coop_hnn_layer_free     
     procedure::propogate => coop_hnn_layer_propogate
     procedure::chain_derv => coop_hnn_layer_chain_derv
     procedure::init_connections => coop_hnn_layer_init_connections
  end type coop_hnn_layer

  type coop_hnn
     COOP_INT::nlayers, nin, nout
     COOP_REAL::best_Err = 1.d99     
     type(coop_hnn_layer),dimension(:),allocatable::layers
     COOP_REAL,dimension(:),allocatable::true_out
   contains
     procedure::init=>coop_hnn_init
     procedure::free=>coop_hnn_free
     procedure::fp => coop_hnn_fp
     procedure::bp => coop_hnn_bp
     procedure::rand => coop_hnn_rand
     procedure::Err => coop_hnn_Err
     procedure::walk => coop_hnn_walk
  end type coop_hnn

contains


  subroutine  coop_hnn_layer_init_connections(this, in, out)
    class(coop_hnn_layer)::this
    COOP_INT,dimension(:),optional::in, out
    COOP_INT::i, j
    if(present(in) .and. present(out))then
       if(size(in) .ne. this%nc .or. size(out) .ne. this%nc) call coop_return_error("coop_hnn_layer_init_connections","# of connections do not match", "stop")

       do i = 1, this%nc
          this%c(i)%in = in(i)
          this%c(i)%out = out(i)
       enddo
    else
       if(this%nc .ne. this%nin*this%nout) call coop_return_error("coop_hnn_layer_init_connections", "# of connections do not match", "stop")
       do i=0, this%nin-1
          do j=0, this%nout-1
             this%c(i*this%nout + j + 1)%in = i+1
             this%c(i*this%nout + j + 1)%out = j+1 
          enddo
       enddo
    endif
  end subroutine coop_hnn_layer_init_connections

  subroutine coop_hnn_layer_free(this)
    class(coop_hnn_layer)::this
    call this%v%free()
    call this%dv%free()
    COOP_DEALLOC(this%c)
    COOP_DEALLOC(this%bias)
    COOP_DEALLOC(this%best_bias)    
    COOP_DEALLOC(this%db)
    this%nin = 0
    this%nc = 0
    this%nout = 0
  end subroutine coop_hnn_layer_free



  subroutine coop_hnn_layer_init(this, nin, nout, nside, nc)
    class(coop_hnn_layer)::this
    COOP_INT::nin, nout, nside
    COOP_INT,optional::nc
    if(this%nin .ne. nin .or. this%nside .ne. nside)then
       call this%v%free
       call this%dv%free
       this%nin = nin
       this%nside = nside
       call this%v%init(nside = nside, nmaps = nin, genre="UNKNOWN", nested = .true.)
       this%dv = this%v
    endif
    if(this%nout .ne. nout)then
       COOP_DEALLOC(this%bias)
       COOP_DEALLOC(this%best_bias)       
       COOP_DEALLOC(this%db)
       this%nout = nout
       if(nout .gt. 0)then
          allocate(this%bias(nout))
          allocate(this%best_bias(nout))          
          allocate(this%db(nout))
       endif
    endif
    if(present(nc))then
       if(this%nc .ne. nc)then
          COOP_DEALLOC(this%c)
          this%nc = nc
          if(nc.gt.0) allocate(this%c(nc))
       endif
    else
       this%nc = this%nin*this%nout
       if(this%nc .gt. 0)then
          allocate(this%c(this%nc))
          call this%init_connections()
       endif
    endif
  end subroutine coop_hnn_layer_init




  subroutine coop_hnn_layer_propogate(this, next)
    class(coop_hnn_layer)::this
    class(coop_hnn_layer)::next
    COOP_INT::i, j
    if(next%nin .ne. this%nout) call coop_return_error("coop_hnn_layer_propogate","layers cannot be connected", "stop") !!quick check
    if(next%nside .eq. this%nside .and. this%nc .gt. 0)then
       do i=1, this%nout
          next%v%map(:, i) = this%bias(i)
       enddo
       do i=1, this%nc
          next%v%map(:, this%c(i)%out) = next%v%map(:, this%c(i)%out) + this%v%map(:, this%c(i)%in)*this%c(i)%weight
       enddo
       next%v%map = coop_hnn_activate(next%v%map)
    elseif(next%nside .lt. this%nside .and. this%nin .eq. next%nin .and. this%nc .eq. 0)then !!pooling
       call coop_healpix_maps_max_udgrade(this%v, next%v)
    else
       call coop_return_error("coop_hnn_layer_propogate","layers cannot be connected", "stop") !!quick check
    endif
  end subroutine coop_hnn_layer_propogate


  subroutine coop_hnn_layer_chain_derv(this, next)
    class(coop_hnn_layer)::this
    class(coop_hnn_layer)::next
    COOP_INT::i, j, div
    if(next%nin .ne. this%nout) call coop_return_error("coop_hnn_layer_chain_derv","layers cannot be connected", "stop") !!quick check
    if(next%nside .eq. this%nside .and. this%nc .gt. 0)then
       do i=1, this%nout
          this%db(i) =  sum(next%dv%map(:,i)*coop_hnn_activate_derv(next%v%map(:,i)))
       enddo
       this%dv%map = 0.
       do i=1, this%nc
          this%dv%map(:, this%c(i)%in) = this%dv%map(:, this%c(i)%in) + next%dv%map(:, this%c(i)%out)  *  coop_hnn_activate_derv(next%v%map(:, this%c(i)%out)) * this%c(i)%weight
          this%c(i)%dw = sum(next%dv%map(:, this%c(i)%out) * coop_hnn_activate_derv(next%v%map(:, this%c(i)%out)) * this%v%map(:, this%c(i)%in))
       enddo
    elseif(next%nside .lt. this%nside .and. this%nin .eq. next%nin .and. this%nc .eq. 0)then !!pooling
       this%db = 0.
       this%dv%map = 0.
       div=(this%nside/next%nside)**2
       do j=1, this%dv%nmaps
          do i=0, next%dv%npix-1
             this%dv%map(i*div+coop_maxloc(this%dv%map(i*div:(i+1)*div-1, j))-1, j) = next%dv%map(i, j)
          enddo
       enddo
    else
       call coop_return_error("coop_hnn_layer_chain_derv","layers cannot be connected", "stop") !!quick check
    endif
       
  end subroutine coop_hnn_layer_chain_derv


  subroutine coop_hnn_rand(this)
    class(coop_hnn)::this
    COOP_INT::i,j
    do i=1, this%nlayers
       if(this%layers(i)%nout .gt. 0) call random_number(this%layers(i)%bias)
       do j=1, this%layers(i)%nc
          call random_number(this%layers(i)%c(j)%weight)
       enddo
    enddo
  end subroutine coop_hnn_rand

  
  subroutine coop_hnn_fp(this)
    class(coop_hnn)::this
    COOP_INT::i
    do i=2, this%nlayers
       call this%layers(i-1)%propogate(this%layers(i))
    enddo
  end subroutine coop_hnn_fp

  subroutine coop_hnn_bp(this)
    class(coop_hnn)::this
    COOP_INT::i, j
    do i=1, this%layers(this%nlayers)%dv%nmaps
       this%layers(this%nlayers)%dv%map(:, i) = sum(this%layers(this%nlayers)%v%map(:, i)) - this%true_out(i)
    enddo
    do i=this%nlayers-1, 1, -1
       call this%layers(i)%chain_derv(this%layers(i+1))
    enddo
  end subroutine coop_hnn_bp
  
  subroutine coop_hnn_free(this)
    class(coop_hnn)::this
    COOP_INT::i
    do i=1, this%nlayers
       call this%layers(i)%free()
    enddo
    COOP_DEALLOC(this%layers)
    this%nlayers = 0
    this%nin = 0
    this%nout = 0
  end subroutine coop_hnn_free

  function coop_hnn_Err(this) result(Err)
    class(coop_hnn)::this    
    COOP_REAL::Err
    COOP_INT::i
    Err = 0.
    do i=1, this%layers(this%nlayers)%v%nmaps
       Err = Err + (this%true_out(i) - sum(this%layers(this%nlayers)%v%map(:,i)))**2
    enddo
    Err  =  Err/2.d0
  end function coop_hnn_Err


  subroutine coop_hnn_Walk(this, step)
    class(coop_hnn)::this
    COOP_REAL::step, s, sumd2, err
    COOP_INT::i, j
    call this%fp()
    err = this%Err()
    if(err .lt. this%best_Err)then
       this%best_Err = err
       do i = 1, this%nlayers-1
          if(this%layers(i)%nc .gt. 0) then                    
             this%layers(i)%best_bias = this%layers(i)%bias
             if(this%layers(i)%nc .gt. 0)  this%layers(i)%c%best_weight = this%layers(i)%c%weight
          endif
       enddo
    elseif(exp(this%best_Err-err)  .lt. coop_random_unit() )then
       do i = 1, this%nlayers-1
          if(this%layers(i)%nc .gt. 0) then          
             this%layers(i)%bias = this%layers(i)%best_bias
             this%layers(i)%c%weight = this%layers(i)%c%best_weight
          endif
       enddo
       call this%fp()
    endif
    call this%bp()

    sumd2 = 0.d0
    do i = 1, this%nlayers-1
       if(this%layers(i)%nc .gt. 0)then
          sumd2 = sumd2 + sum(this%layers(i)%db**2) + sum(this%layers(i)%c%dw**2)
       endif
    enddo
    
    s = step*err/max(sumd2, 1.d-99)
    
    do i = 1, this%nlayers-1
       if(this%layers(i)%nc .gt. 0)then
          do j=1, this%layers(i)%nout
             this%layers(i)%bias(j) = this%layers(i)%bias(j) - s* coop_random_Gaussian() * this%layers(i)%db(j)
          enddo
          do j = 1, this%layers(i)%nc
             this%layers(i)%c(j)%weight = this%layers(i)%c(j)%weight - s* coop_random_Gaussian() * this%layers(i)%c(j)%dw
          enddo
       endif
    enddo
    
  end subroutine coop_hnn_Walk
  
  

  subroutine coop_hnn_init(this, nmaps, nside, input, delta_ell)
    class(coop_hnn)::this
    COOP_INT,dimension(:)::nmaps, nside
    COOP_INT::i
    COOP_INT::delta_ell, l
    type(coop_healpix_maps)::input
    call this%free()
    this%nlayers = coop_getdim("hnn_init", size(nmaps), size(nside))
    allocate(this%layers(this%nlayers))
    this%nin = nmaps(1)
    this%nout = nmaps(this%nlayers)
    allocate(this%true_out(this%nout))
    do i=1, this%nlayers-1
       if(nside(i) .eq. nside(i+1))then
          call this%layers(i)%init(nmaps(i), nmaps(i+1), nside(i))
       elseif(nmaps(i) .eq. nmaps(i+1) .and. nside(i) .gt. nside(i+1))then
          call this%layers(i)%init(nmaps(i), nmaps(i+1), nside(i), 0)
       else
          call coop_return_error("hnn_init: nside / nmaps do not match")
       endif
    enddo
    call this%layers(this%nlayers)%init(nmaps(this%nlayers), 0, nside(this%nlayers))
    call this%rand()
    call input%allocate_alms(lmax = delta_ell * nmaps(1)-1)
    call input%map2alm(lmax = delta_ell * nmaps(1)-1)
    call this%layers(1)%v%allocate_alms(lmax = delta_ell * nmaps(1)-1)    
    do i=1, nmaps(1)
       do l= (i-1)*delta_ell, i*delta_ell-1
          this%layers(1)%v%alm(l,0:l, i) = input%alm(l, 0:l, 1)
       enddo
       call this%layers(1)%v%alm2map()
    enddo
  end subroutine coop_hnn_init
  
  

#if ACTIVATE_FUNCTION == 0 
  function coop_hnn_activate_s(x) result(f)
    COOP_SINGLE::x, f
    f = 1.d0/(1.d0+exp(-x))
  end function coop_hnn_activate_s

  function coop_hnn_activate_derv_s(f) result(df)
    COOP_SINGLE::f,df
    df  = f*(1.d0-f)
  end function coop_hnn_activate_derv_s

#elif ACTIVATE_FUNCTION == 1

  function coop_hnn_activate_s(x) result(f)
    COOP_SINGLE::x, f
    f = max(0.d0, x)
  end function coop_hnn_activate_s

  function coop_hnn_activate_derv_s(f) result(df)
    COOP_SINGLE::f,df
    if(f.gt.0.d0)then
       df = 1.d0
    else
       df = 0.d0
    endif
  end function coop_hnn_activate_derv_s

#elif ACTIVATE_FUNCTION == 2

  function coop_hnn_activate_s(x) result(f)
    COOP_SINGLE::x, f
    f = tanh(x)
  end function coop_hnn_activate_s

  function coop_hnn_activate_derv_s(f) result(df)
    COOP_SINGLE::f,df
    df = 1.d0-f**2
  end function coop_hnn_activate_derv_s

#elif ACTIVATE_FUNCTION == 3

  function coop_hnn_activate_s(x) result(f)
    COOP_SINGLE::x, f
    if(x.gt.0.d0)then
       f = x
    else
       f = x/10.d0
    endif
  end function coop_hnn_activate_s

  function coop_hnn_activate_derv_s(f) result(df)
    COOP_SINGLE::f,df
    if(f.gt.0.d0)then
       df = 1.d0
    else
       df = 0.1d0
    endif
  end function coop_hnn_activate_derv_s
  
#endif

  function coop_hnn_activate_v(x) result(f)
    COOP_SINGLE,dimension(:)::x
    COOP_SINGLE:: f(size(x))
    COOP_INT::i, n
    n = size(x)
    !$omp parallel do
    do i=1, n
       f(i) = coop_hnn_activate_s(x(i))
    enddo
    !$omp end parallel do
  end function coop_hnn_activate_v

  function coop_hnn_activate_derv_v(f) result(df)
    COOP_SINGLE,dimension(:)::f
    COOP_SINGLE:: df(size(f))
    COOP_INT::i, n
    n = size(f)
    !$omp parallel do
    do i=1, n
       df(i) = coop_hnn_activate_derv_s(f(i))
    enddo
    !$omp end parallel do
  end function coop_hnn_activate_derv_v

  function coop_hnn_activate_vv(x) result(f)
    COOP_SINGLE,dimension(:,:)::x
    COOP_SINGLE:: f(size(x,1), size(x, 2))
    COOP_INT::i, j,  n1, n2
    n1 = size(x, 1)
    n2 = size(x, 2)
    !$omp parallel do private(i, j)
    do j=1, n2
       do i=1, n1
          f(i, j) = coop_hnn_activate_s(x(i, j))
       enddo
    enddo
    !$omp end parallel do
  end function coop_hnn_activate_vv

  function coop_hnn_activate_derv_vv(f) result(df)
    COOP_SINGLE,dimension(:,:)::f
    COOP_SINGLE:: df(size(f,1), size(f,2))
    COOP_INT::i,j, n1, n2
    n1 = size(f,1)
    n2 = size(f,2)
    !$omp parallel do private(i, j)
    do j=1, n2
       do i=1, n1
          df(i, j) = coop_hnn_activate_derv_s(f(i, j))
       enddo
    enddo
    !$omp end parallel do
  end function coop_hnn_activate_derv_vv
  
  !!==============================


  
end module coop_hnn_mod
    

