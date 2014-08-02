!-----------------------------------------------------------------------
! PSCF - Polymer Self-Consistent Field Theory
! Copyright (2007) David C. Morse
! email: morse@cems.umn.edu
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation. A copy of this license is included in
! the LICENSE file in the top-level PSCF directory. 
!-----------------------------------------------------------------------
!****m  scf/grid_mod
! PURPOSE
!    Declare dimensions ngrid(3) of fft grid
!    Declare global data structures that are defined on an fft grid
!    Routines to allocate, deallocate and manipulate grid data
!    Functions G_to_fft, G_to_bz, norm to manipulate wavevectors
! SOURCE
!-----------------------------------------------------------------------
module grid_mod
   use const_mod
   !use string_mod
   implicit none

   private

   ! Public data structures
   public :: ngrid          ! dimensions ngrid(:)=(N1,N2,N3) of grid
   public :: rho_grid       ! rho on r-grid   (r,N_monomer)
   public :: omega_grid     ! omega on r-grid (r,N_monomer)
   public :: ksq_grid       ! k**2 on k-grid 

   ! Public procedures
   public :: input_grid
   public :: allocate_grid
   public :: deallocate_grid
   public :: make_ksq
   public :: max_Gabs
   public :: G_to_fft
   public :: G_to_bz
   public :: Greal_to_bz
   public :: norm

   integer                  :: ngrid(3)
   real(long), ALLOCATABLE  :: rho_grid(:,:,:,:)
   real(long), ALLOCATABLE  :: omega_grid(:,:,:,:)
   real(long), ALLOCATABLE  :: ksq_grid(:,:,:)
   !***

   !---------------------------------------------------------------
   !****v grid_mod/ngrid
   ! VARIABLE
   !   integer ngrid(3) :: Number of grid points in each direction
   !***
   !---------------------------------------------------------------
   !****v grid_mod/rho_grid
   ! VARIABLE
   !   real(long) rho_grid(0:,0:,0:,N_monomer)
   !              density of each monomer type
   !***
   !---------------------------------------------------------------
   !****v grid_mod/omega_grid
   ! VARIABLE
   !   real(long) omega_grid(0:,0:,0:,N_monomer)
   !              potential for each monomer type
   !***
   !---------------------------------------------------------------
   !****v grid_mod/ksq_grid
   ! VARIABLE
   !   real(long) ksq_grid(0:,0:,0:)
   !              wave vector magnitude squared, for FFT use
   !***
   !---------------------------------------------------------------

contains

   !--------------------------------------------------------------   
   !****p grid_mod/input_grid
   ! SUBROUTINE
   !   input_grid
   ! PURPOSE
   !   Input integer array ngrid(1:dim) from input script
   ! SOURCE
   !---------------------------------------------------------------
   subroutine input_grid
   !***
   use io_mod, only : input

   call input(ngrid(:),dim,f='A')
   if (dim < 3) then
      ngrid(dim+1:3) = 1
   end if

   end subroutine input_grid
   !==============================================================


   !--------------------------------------------------------------   
   !****p grid_mod/allocate_grid
   ! SUBROUTINE
   !   allocate_grid
   ! PURPOSE
   !   Allocate memory needed by the grid module data structures
   !   rho_grid, omega_grid, ksq_grid
   ! SOURCE
   !---------------------------------------------------------------
   subroutine allocate_grid(N_monomer)
   implicit none
   integer, intent(IN) :: N_monomer
   !***

   integer  :: nx, ny, nz
   integer  :: error

   nx=ngrid(1)-1
   ny=ngrid(2)-1
   nz=ngrid(3)-1

   ALLOCATE(rho_grid(0:nx,0:ny,0:nz,N_monomer),STAT=error)
   if (error /= 0) stop "rho_grid allocation error!"

   ALLOCATE(omega_grid(0:nx,0:ny,0:nz,N_monomer),STAT=error)
   if (error /= 0) stop "rho_grid allocation error!"

   ALLOCATE(ksq_grid(0:(nx+1)/2,0:ny,0:nz),STAT=error)
   if (error /= 0) stop "rho_grid allocation error!"

   end subroutine allocate_grid
   !==============================================================


   !--------------------------------------------------------------   
   !****p grid_mod/make_ksq
   ! SUBROUTINE
   !   make_ksq
   ! PURPOSE
   !   Calculate values of k**2 and store in ksq_grid(:,:,:)
   ! COMMENT
   !   calculates norm(G_to_bz(G)) for each wavevector G
   ! SOURCE
   !--------------------------------------------------------------   
   subroutine make_ksq(G_basis)
   implicit none
   real(long)  :: G_basis(:,:)
   !***
   integer :: G(3), Gbz(3)
   integer :: i1, i2, i3

   do i1=0, ngrid(1)/2
      G(1)=i1
      do i2=0, ngrid(2)-1
         G(2)=i2
         do i3=0, ngrid(3)-1
            G(3)=i3
            Gbz=G_to_bz(G)
            ksq_grid(i1,i2,i3)=norm(Gbz,G_basis)
         enddo
      enddo
   enddo

   end subroutine make_ksq
   !==============================================================


   !--------------------------------------------------------------   
   !****p  grid_mod/deallocate_grid
   ! SUBROUTINE
   !    destroy_grid
   ! PURPOSE
   !    deallocate(return) the memory used by grid_data
   ! SOURCE
   !--------------------------------------------------------------   
   subroutine deallocate_grid
   implicit none
   !***
   integer     :: error      ! deallocation error index

   DEALLOCATE(ksq_grid,STAT=error)
   if (error /= 0) stop "rho_grid deallocation error!"

   DEALLOCATE(rho_grid,STAT=error)
   if (error /= 0) stop "rho_grid deallocation error!"

   DEALLOCATE(omega_grid,STAT=error)
   if (error /= 0) stop "rho_grid deallocation error!"

   end subroutine deallocate_grid
   !==============================================================


   !--------------------------------------------------------------   
   !****p  grid_mod/G_to_fft
   ! SUBROUTINE
   !   G_to_fft(G)
   ! PURPOSE
   !   Shift any reciprocal wave vector G = [G(1),..,G(dim)] to 
   !   an equivalent wavevector with indices 0 <= G(i) < ngrid(i) 
   ! SOURCE
   !--------------------------------------------------------------   
   function G_to_fft(G)
   implicit none
   integer             :: G_to_fft(3)
   integer, intent(IN) :: G(:)
   !***
   integer             :: i
   G_to_fft=0
   do i=1,dim
      G_to_fft(i)=modulo( G(i), ngrid(i) )
   end do
   end function G_to_fft
   !==============================================================


   !------------------------------------------------------------------
   !****p  grid_mod/G_to_bz
   ! SUBROUTINE
   !    G_to_bz
   ! PURPOSE
   !   Shift any reciprocal wave vector of the crystal to the first 
   !   Brillouin zone of the grid unit cell, i.e., to the equivalent
   !   wavevector with the smallest absolute magnitude. 
   ! ARGUMENTS
   !   integer G(3) - array containing integer indices of wavevector
   ! SOURCE
   !------------------------------------------------------------------
   function G_to_bz(G)
   use unit_cell_mod, only : G_basis
   implicit none
   integer             :: G_to_bz(3)
   integer, intent(IN) :: G(:)         ! wave vector index
   !***

   ! local variables
   integer, parameter :: mup=1, mdn=-1
   integer            :: G_min(3),G_try(3)
   integer            :: i1,i2,i3
   real(long)         :: Gsq_min

   G_try   = 0
   G_min   = 0
   Gsq_min = 1.0E+10
   select case(dim)
   case(1)
      do i1=mup,mdn,-1
         G_try(1)=G(1)+i1*ngrid(1)
         call choose
      enddo
   case(2)
      do i1=mup,mdn,-1
         G_try(1)=G(1)+i1*ngrid(1)
         do i2=mup,mdn,-1
            G_try(2)=G(2)+i2*ngrid(2) 
            call choose
         enddo
      enddo
   case(3)
      do i1=mup,mdn,-1
         G_try(1)=G(1)+i1*ngrid(1)
         do i2=mup,mdn,-1
            G_try(2)=G(2)+i2*ngrid(2)
            do i3=mup,mdn,-1
               G_try(3)=G(3)+i3*ngrid(3)
               call choose
            enddo
         enddo
      enddo
   end select
   G_to_bz=G_min

   contains ! internal subroutine choose

     !-----------------------------------------
      subroutine choose
      real(long), parameter :: delta=1.0E-8
      real(long)            :: Gsq
      Gsq=norm(G_try,G_basis)
      if (Gsq < Gsq_min - delta) then
         G_min(1:dim) = G_try(1:dim)
         Gsq_min      = Gsq 
      endif
      end subroutine choose
     !-----------------------------------------

   end function G_to_bz
   !==============================================================

   function Greal_to_bz(Greal)
   use unit_cell_mod, only : G_basis
   implicit none
   real(long)                :: Greal_to_bz(3)
   real(long),intent(IN)     :: Greal(:)         ! wave vector index
!   real(long),intent(IN)  :: G_basis(:,:)

!  local variables
   integer, parameter :: mup=1, mdn=-1
   real(long)            :: G_min(3),G_try(3)
   integer            :: i1,i2,i3
   real(long)         :: Gsq_min

   G_try   = 0
   G_min   = 0
   Gsq_min = 1.0E+10

   select case(dim)
   case(1)
      do i1=mup,mdn,-1
         G_try(1)=Greal(1) + i1 * ngrid(1)
         call choose
      enddo
   case(2)
      do i1=mup,mdn,-1
         G_try(1)=Greal(1)+ i1*ngrid(1)
         do i2=mup,mdn,-1
            G_try(2)=Greal(2)+i2*ngrid(2) 
            call choose
         enddo
      enddo
   case(3)
      do i1=mup,mdn,-1
         G_try(1)=Greal(1)+i1*ngrid(1)
         do i2=mup,mdn,-1
            G_try(2)=Greal(2)+i2*ngrid(2)
            do i3=mup,mdn,-1
               G_try(3)=Greal(3)+i3*ngrid(3)
               call choose
            enddo
         enddo
      enddo
   end select
   Greal_to_bz=G_min

   contains      ! internal subroutine
     !-----------------------------------------
      subroutine choose
        use group_mod,only: operator(.dot.)
      real(long), parameter :: delta=1.0E-8
      real(long)            :: Gsq
      real(long)            :: v(3)
      integer               ::i
      
      v = 0.0_long
      do i = 1,dim
         v(i) = G_try(:) .dot. G_basis(:,i)
      end do
      Gsq = (v .dot. v)
!      Gsq=norm(G_try,G_basis)
      if (Gsq < Gsq_min - delta) then
         G_min(1:dim) = G_try(1:dim)
         Gsq_min      = Gsq 
      endif
      end subroutine choose
     !-----------------------------------------

    end function Greal_to_bz
! ==============================================================




   !--------------------------------------------------------------   
   !****p  grid_mod/max_Gabs
   ! SUBROUTINE
   !   max_Gabs
   ! PURPOSE
   !   calculate the max magniture of vectors in FFT grid
   ! SOURCE
   !--------------------------------------------------------------   
   function max_Gabs(G_basis)
   implicit none
   real(long)            :: max_Gabs
   real(long),intent(IN) :: G_basis(:,:)
   !***
   real(long)            :: length,tmp
   real(long),parameter  :: delta=1.0E-10
   integer               :: G(3),Gbz(3)
   integer               :: i1,i2,i3
   length=1.0E-15
   do i3=0,ngrid(3)-1
      G(3)=i3
      do i2=0,ngrid(2)-1
         G(2)=i2
         do i1=0,ngrid(1)-1
            G(1)=i1
            Gbz=G_to_bz(G)
            tmp=norm(Gbz,G_basis)
            if(tmp > length+delta) length=tmp
         end do
      end do
   end do
   max_Gabs=sqrt(length)
   end function max_Gabs
   !==============================================================


   !--------------------------------------------------------------
   !****p  grid_mod/norm
   ! SUBROUTINE
   !   norm(G,G_basis)
   ! PURPOSE
   !   Calculate the squared magnitude of vector G
   ! SOURCE
   !--------------------------------------------------------------
   function norm(G,G_basis)
   use group_mod,     only : operator(.dot.)
   implicit none
   real(long)             :: norm
   integer,intent(IN)     :: G(:)
   real(long),intent(IN)  :: G_basis(:,:)
   !***
   real(long)             :: Gvec(3)

   Gvec=G.dot.G_basis
   norm=Gvec.dot.Gvec
   end function norm
   !=================================================================


end module grid_mod
