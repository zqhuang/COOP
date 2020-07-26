module coop_stoinf_mod
  use coop_wrapper_utils  
  use coop_lattice_fields_mod
  implicit none

#include "constants.h"

  COOP_INT,parameter::stoinf_max_num_chidren = 64

  type stoinf_patch
     COOP_INT::coor(3) =  (/ 0, 0, 0 /) !!the relative coordinate in the parent node, if exists
     COOP_INT::nside=0
     COOP_INT::nflds=0
     COOP_INT::num_children=0
     COOP_REAL::lna = 0.d0
     COOP_REAL, dimension(:,:,:,:), allocatable :: f !!field values
     COOP_REAL, dimension(:,:,:,:), allocatable :: pi !!field momenta (\dot\phi a^3)
     type(stoinf_patch),pointer:: children(:) => null()
   contains
     procedure::init => stoinf_patch_init
     procedure::free => stoinf_patch_free
     procedure::evolve => stoinf_patch_evolve
     procedure::add_children=> stoinf_patch_add_children     
  end type stoinf_patch

contains


  subroutine stoinf_patch_free(this)
    class(stoinf_patch)::this
    COOP_INT::i
    if(associated(this%children))then
       do i=1, this%num_children
          call this%children(i)%free()
       enddo
       nullify(this%children)
    endif
    if(allocated(this%f))deallocate(this%f)
    if(allocated(this%pi))deallocate(this%pi)    
    this%nflds = 0
    this%nside = 0
    this%num_children = 0
  end subroutine stoinf_patch_free
  
  
  subroutine stoinf_patch_init(this, nside, nflds)
    class(stoinf_patch)::this
    COOP_INT::nside, nflds
    call this%free()
    this%nside = nside
    this%nflds = nflds
    allocate(this%f(this%nside, this%nside, this%nside, this%nflds))
    allocate(this%pi(this%nside, this%nside, this%nside, this%nflds))    
  end subroutine stoinf_patch_init

  subroutine stoinf_patch_add_children(this, num_children, nside, coors)
    class(stoinf_patch)::this
    COOP_INT::nside, num_children
    COOP_INT::coors(3,num_children)
    COOP_INT::i, j
    if(any(coors .lt. 1) .or. any(coors .gt. this%nside))then
       write(*,*) minval(coors), maxval(coors)
       write(*,*) "nside = ", this%nside
       write(*,*) "fatal error  in add_children: coordinates exceed the boundary."
       stop
    endif
    if(associated(this%children))then
       do i=1, this%num_children
          call this%children(i)%free()
       enddo
       nullify(this%children)
    endif
    allocate(this%children(num_children))
    this%num_children = num_children
    do i=1, this%num_children
       this%children(i)%coor = coors(1:3, i)
       call this%children(i)%init(nside = nside, nflds = this%nflds)
       do j=1, this%nflds
          this%children(i)%f(:,:,:, j) = this%f(coors(1,i),coors(2,i),coors(3,i), j)
          this%children(i)%pi(:,:,:,j) = this%pi(coors(1,i),coors(2,i),coors(3,i), j)          
       enddo
    enddo
  end subroutine stoinf_patch_add_children


  subroutine stoinf_patch_evolve(this, dN, noise_norm, noise_mean)
    class(stoinf_patch)::this
    COOP_REAL::dN, noise_mean(this%nflds), noise_norm(this%nflds)
    COOP_REAL::noise(this%nside, this%nside, this%nside, this%nflds )
    integer::fld, i, j, k
    COOP_REAL::k1(2*this%nflds), k2(2*this%nflds), k3(2*this%nflds), k4(2*this%nflds), sqrtdN, local_eps, local_H
    !!-------evaluate golbal Hubble ------------------
!!$    COOP_REAL::global_H, global_eps, Etot, local_H, n3
!!$    n3 = dble(this%nside)**3
!!$    global_eps = sum(this%pi**2)/2.d0/n3/coop_lattice_Mpsq
!!$    Etot = 0.d0
!!$    !$omp parallel do private(i,j,k) reduction(+:Etot)
!!$    do k=1, this%nside
!!$       do j=1, this%nside
!!$          do i=1, this%nside
!!$             Etot = Etot + coop_lattice_fields_V(this%f(:, i, j, k)) 
!!$          enddo
!!$       enddo
!!$    enddo
!!$    !$omp end parallel do
!!$    global_H = sqrt(Etot/n3/(3.d0-global_eps)/coop_lattice_Mpsq)  !!global H is used to generate the white noise
    !!===================================================================================
    do fld = 1, this%nflds
       call coop_generate_white_noise(n=this%nside, output= noise(:, :, :, fld), norm = noise_norm(fld) , mean = noise_mean(fld))
    enddo
    sqrtdN = sqrt(dN)
    do k=1, this%nside
       do j=1, this%nside
          do i=1, this%nside
             !!Runge Kutta 4th
             k1(1:this%nflds) = this%pi(i, j, k, :) + noise(i, j, k, :)/sqrtdN
             k1(this%nflds+1, 2*this%nflds) = -(3.d0-loal_eps)* this%pi(i, j, k, :)- coop_lattice_fields_dVdphi(this%f(i,j,k,:))/local_H**2
             
             this%pi(fld, i,j,k) = this%pi(fld, i, j, k) -  coop_lattice_fields_dVdphi(this%f(:,i,j,k))*exp(3.d0*this%lna)*dN/local_H
          enddo
       enddo
    enddo
  end subroutine stoinf_patch_evolve

  


end module coop_stoinf_mod
