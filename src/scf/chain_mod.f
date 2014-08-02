!-----------------------------------------------------------------------
! PSCF - Polymer Self-Consistent Field Theory
! Copyright (2007) David C. Morse
! email: morse@cems.umn.edu
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation. A copy of this license is included in
! the LICENSE file in the top-level PSCF directory. 
!****m scf/chain_mod ---------------------------------------------------
! MODULE
!   chain_mod
! PURPOSE
!   Define derived type chain_grid_type, and allocate arrays of this type.
!   Subroutines to make and destroy objects of chain_type.
! SOURCE
!-----------------------------------------------------------------------
module chain_mod
   use const_mod
   use fft_mod, only : fft_plan
   implicit none

   private
   public :: chain_grid_type
   public :: null_chain_grid
   public :: make_chain_grid
   public :: destroy_chain_grid
   !***

   !****t chain_mod/chain_grid_type  -----------------------------------
   !  TYPE 
   !     chain_grid_type
   !  PURPOSE
   !     Data structures defining discretization of s for a chain.
   !     Pointers to qf(r,s) and qr(r,s) functions for a chain
   !     Pointers to del_qf and del_qr functions for perturbation theory
   !  SOURCE
   !-------------------------------------------------------------------
   type chain_grid_type
      integer,    pointer     :: block_bgn(:)     ! 1st element of block
      real(long), pointer     :: block_ds(:)      ! step size for block
      real(long), pointer     :: qf(:,:,:,:)      ! function qf(x,y,x,s) 
      real(long), pointer     :: qr(:,:,:,:)      ! function qr(x,y,z,s) 
      real(long), pointer     :: rho(:,:,:,:)     ! rho(x,y,z,s,block)
      type(fft_plan)          :: plan             ! fft plan, see fft_mod
      real(long)              :: bigQ             ! chain partition func.
      complex(long), pointer  :: del_qf(:,:,:,:)  ! perturbation in qf
      complex(long), pointer  :: del_qr(:,:,:,:)  ! perturbation in qr
      complex(long)           :: delQ             ! perturbation in bigQ
   end type
   !***

contains

   !------------------------------------------------------------------
   !****p chain_mod/make_chain_grid
   ! PURPOSE
   !    Nullify the pointers to ensure they have the "dissociated"
   !    status.
   ! SOURCE
   !------------------------------------------------------------------
   subroutine null_chain_grid(chain)
   implicit none

   type(chain_grid_type), intent(INOUT)  :: chain
   !***

   nullify( chain%block_bgn, &
            chain%block_ds,  &
            chain%qf,        &
            chain%qr,        &
            chain%rho,       &
            chain%del_qf,    &
            chain%del_qr)

   end subroutine null_chain_grid
   !------------------------------------------------------------------

   !------------------------------------------------------------------
   !****p chain_mod/make_chain_grid
   ! SUBROUTINE
   !    make_chains(chain,plan,nblk,blk_length,ds,allocate_q)
   ! PURPOSE
   !    allocate memory for a single chain_grid_type variable
   !    initiate chain(:)%block_bgn, chain(:)%block_ds
   ! ARGUMENTS
   !    plan       = grid dimensions and FFT plan, see fft_mod
   !    nblk       = # of blocks of a single chain
   !    blk_length = block lengths
   !    ds         = segment length used to discretize the block
   !    allocate_q = true if qf and qr need to be allocated
   ! COMMENTS
   !    The # of segments for each block need to be even, because
   !    of the use of Simpson's rule in density/stress calculation
   ! SOURCE
   !------------------------------------------------------------------
   subroutine &
       make_chain_grid(chain,plan,nblk,blk_length,ds0,allocate_q,perturb,order)
   implicit none
   type(chain_grid_type), intent(INOUT)  :: chain
   type(fft_plan),        intent(IN)     :: plan
   integer,               intent(IN)     :: nblk
   real(long),            intent(IN)     :: blk_length(:)
   real(long),            intent(IN)     :: ds0
   logical,               intent(IN)     :: allocate_q
   logical,optional,      intent(IN)     :: perturb
   integer,optional,      intent(IN)     :: order
   !***
   real(long)            :: ds
   integer               :: iblk         ! index to block
   integer               :: bgns         ! segments counter
   integer               :: nx,ny,nz,i   ! loop indices
   integer               :: error        ! allocation error-message

   nx=plan%n(1)-1
   ny=plan%n(2)-1
   nz=plan%n(3)-1
   chain%plan=plan

   if ( .not. associated(chain%block_bgn) ) then
      allocate(chain%block_bgn(nblk+1),STAT=error)
      if (error /= 0) stop 'chain%block_bgn allocation error'
   end if

   if ( .not. (associated(chain%block_ds)) ) then
      allocate(chain%block_ds(nblk),STAT=error)
      if (error /= 0) stop 'chain%block_ds allocation error'
   end if

   if ( .not. associated(chain%rho) ) then
      allocate(chain%rho(0:nx,0:ny,0:nz,nblk),STAT=error)
      if (error /= 0) stop 'chain%rho allocation error!'
   end if

   bgns=1
   chain%block_bgn(1)=bgns

   ds = ds0

   do i=1, nblk    ! loop over blocks
      iblk=int(blk_length(i)/ds/2.0_long+0.5_long)

      if (iblk == 0) then
         iblk = 1
         chain%block_ds(i)=blk_length(i)/2.0_long
      else
         chain%block_ds(i)=blk_length(i)/dble(iblk)/2.0_long
      endif

      if ( present(order) ) then
        chain%block_ds(i) = chain%block_ds(i) / 2.0_long**order
        iblk = iblk * 2**order
      end if

      bgns = bgns + iblk * 2
      chain%block_bgn(i+1) = bgns
   end do

   if ( allocate_q ) then
      allocate(chain%qf(0:nx,0:ny,0:nz,bgns),STAT=error)
      if (error /= 0) stop "chain%qf allocation error!"

      allocate(chain%qr(0:nx,0:ny,0:nz,bgns),STAT=error)
      if (error /= 0) stop "chain%qr allocation error!"
   end if

   if ( (present(perturb)) .and. perturb ) then
      allocate(chain%del_qf(0:nx,0:ny,0:nz,bgns),STAT=error)
      if (error /= 0) stop "chain%qf allocation error!"
      
      allocate(chain%del_qr(0:nx,0:ny,0:nz,bgns),STAT=error)
      if (error /= 0) stop "chain%qr allocation error!"
    end if

   end subroutine make_chain_grid
   !==================================================================


   !-------------------------------------------------------------------
   !****p chain_mod/destroy_chain_grid
   ! SUBROUTINE
   !    destroy_chain_grid(chain)
   ! PURPOSE
   !    Deallocate memories use by chain%...
   ! ARGUMENTS
   !    chain = the chain_grid_type to be deallocated
   ! SOURCE
   !-------------------------------------------------------------------
   subroutine destroy_chain_grid(chain)
   implicit none
   type(chain_grid_type) :: chain
   !***
   integer  :: error     ! deallocation error msg

   if ( associated(chain%block_bgn) ) then
      deallocate(chain%block_bgn,STAT=error)
      if (error /= 0) stop "chain%block_bgn deallocation error!"
   endif

   if ( associated(chain%block_ds) ) then
      deallocate(chain%block_ds,STAT=error)
      if (error /= 0) stop "chain%block_ds deallocation error!"
   endif

   if ( associated(chain%rho) ) then
      deallocate(chain%rho,STAT=error)
      if (error /= 0) stop "chain%rho deallocation error!"
   endif

   if ( associated(chain%qf) ) then
      deallocate(chain%qf,STAT=error)
      if (error /= 0) stop "chain%qf deallocation error!"
   endif

   if ( associated(chain%qr) ) then
      deallocate(chain%qr,STAT=error)
      if (error /= 0) stop "chain%qr deallocation error!"
   endif

   if ( associated(chain%del_qf) ) then
      deallocate(chain%del_qf,STAT=error)
      if (error /= 0) stop "chain%qf deallocation error!"
   endif

   if ( associated(chain%del_qr) ) then
      deallocate(chain%del_qr,STAT=error)
      if (error /= 0) stop "chain%qr deallocation error!"
   endif

   end subroutine destroy_chain_grid
   !====================================================================

end module chain_mod
