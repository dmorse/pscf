!-----------------------------------------------------------------------
! PSCF - Polymer Self-Consistent Field Theory
! Copyright (2007) David C. Morse
! email: morse@cems.umn.edu
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation. A copy of this license is included in
! the LICENSE file in the top-level PSCF directory. 
!----------------------------------------------------------------------
!****m scf/response_step_mod
! PURPOSE
!   Integrator for the first-order perturbation theory for q(r,s).
!   Propagates both q(r,s) and delq(r,s)
! SOURCE
!----------------------------------------------------------------------
module response_step_mod
   use const_mod
   use chemistry_mod
   use fft_mod
 
   private
   public :: response_step_startup
   public :: fft_step
   public :: pspropagate
   public :: alloc_exparrays
   public :: calc_exparrays
   public :: calc_exp_grid
   public :: ps_propagate
   !***

   ! Private module arrays
   real(long), allocatable     :: exp_ksq_new(:,:,:,:)
   real(long), allocatable     :: exp_ksq(:,:,:,:)
   real(long), allocatable     :: exp_omega(:,:,:,:)
   real(long), allocatable     :: q_step(:,:,:,:)
   complex(long), allocatable  :: delq_step(:,:,:,:)
   complex(long), allocatable  :: fin_step(:,:,:)
   complex(long), allocatable  :: fout_step(:,:,:)
 

contains

   !-----------------------------------------------------------------
   !****p response_step_mod/response_step_startup
   ! SUBROUTINE
   !   response_step_startup(ngrid,order)
   ! PURPOSE
   !   Allocate the module variable exponential arrays
   ! ARGUMENTS
   !   integer ngrid  - number of grid points
   !   integer order  - order of Richardson extrapolation
   ! SOURCE
   !-----------------------------------------------------------------
   subroutine response_step_startup(Ngrid, order)
   integer,intent(IN)  :: ngrid(3)
   integer,intent(IN)  :: order
   !***
   integer             :: err
   integer             :: nx,ny,nz,nstep
     
   nx = ngrid(1)-1
   ny = ngrid(2)-1
   nz = ngrid(3)-1
   nstep = 2 ** order
 
   if (.not.(allocated(q_step))) allocate (q_step(0:nx,0:ny,0:nz,nstep+1),stat=err)
   if (err /= 0) stop "Error allocating q_step in propagate"
 
   if (.not.(allocated(delq_step))) allocate (delq_step(0:nx,0:ny,0:nz,nstep+1),stat=err)
   if (err /= 0) stop "Error allocating delq_step in propagate"
 
   if (.not.(allocated(fin_step))) allocate (fin_step(0:nx,0:ny,0:nz),stat=err)
   if (err /= 0) stop "Error allocating fin_step in propagate"
 
   if (.not.(allocated(fout_step))) allocate (fout_step(0:nx,0:ny,0:nz),stat=err)
   if (err /= 0) stop "Error allocating fout_step in propagate"
 
   if (.not.(allocated(exp_ksq_new))) allocate(exp_ksq_new(0:nx,0:ny,0:nz,0:order+1),stat=err)
   if (err /= 0) stop "Error allocating exp_ksq_new in alloc_exparray"
 
   if (.not.(allocated(exp_omega))) allocate(exp_omega(0:nx,0:ny,0:nz,0:order+1),stat=err)
   if (err /= 0) stop "Error allocating exp_omega in alloc_exparray"
   
   if (.not.(allocated(exp_ksq))) allocate(exp_ksq(0:(nx+1)/2,0:ny,0:nz,0:order+1),stat=err)
   if (err /= 0) stop "Error allocating exp_ksq in alloc_exparray"
 
   end subroutine response_step_startup
   !==================================================================
 
 
   !-----------------------------------------------------------------
   !****p response_step_mod/ps_propagate
   ! SUBROUTINE
   !   ps_propagate(plan,planc,ds,delqin,q0,exp_gpdotr,order,delta_omega, &
   !       monomer,p_mon,delqout,dsloc,qout)
   ! PURPOSE
   !   Propagate q and delq by solving the inhomog pde. 
   !   call  fft_step and fftc_step for unperturbed and perturbation
   !   partition functions.
   ! ARGUMENTS
   !   fft_plan planc      -  see fft_mod for details
   !   complex  delqin     -  delta q at beginning of the block
   !   complex  delqout    -  Output delta q array
   !   complex  exp_gpdotr -  exponential array of G.r for perturbation at G+k
   !    complex  f_inhomo  -  = d_w * q_0(1:2**order)
   !    integer  p_mon     -  index of the perturbed monomer
   !    integer  monomer   -  index of the monomer of the concerned chain site
   !                          over which the propagation is being implemented
   !    integer  order     -  order of extrapolation
   !    real     dsloc     -  contour step size for this block
   ! SOURCE
   !-----------------------------------------------------------------
   subroutine ps_propagate(delqin, delqout, p_mon, exp_gpdotr, &
                          f_inhomo, monomer, planc, order, dsloc) 
     
   complex(long), intent(IN)  :: delqin(:,:,:)
   complex(long), intent(OUT) :: delqout(:,:,:)
   complex(long), intent(IN)  :: exp_gpdotr(:,:,:)
   real(long),    intent(IN)  :: f_inhomo(:,:,:,:)
   integer,       intent(IN)  :: p_mon
   integer,       intent(IN)  :: monomer
   type(fft_plan)             :: planc
   real(long),    intent(IN)  :: dsloc
   integer,       intent(IN)  :: order
   !***
 
   integer  :: i
 
   !------------------------------------------------
   ! Keep using Amit's definition of arrays for now
   ! change the dimension later after the code has
   ! been successfully modified.
   !------------------------------------------------
   delq_step(:,:,:,1) = delqin(:,:,:)
 
   do i = 1, 2**order
      call fftc_step(planc,order,delq_step(:,:,:,i),delq_step(:,:,:,i+1))
 
      if (monomer .eq. p_mon) then
         fin_step(:,:,:) = f_inhomo(:,:,:,i) * exp_gpdotr(:,:,:)
 
         call fftc_step(planc,order+1,fin_step,fout_step)
 
         delq_step(:,:,:,i+1) = delq_step(:,:,:,i+1)         &
                   - fout_step(:,:,:) * dsloc
      end if
   end do
  
   delqout(:,:,:) = delq_step(:,:,:,2**order+1)
 
   end subroutine ps_propagate
   !==================================================================
 
 
   !-----------------------------------------------------------------
   !****P response_step_mod/fft_step
   ! SUBROUTINE
   !   fft_step(plan,order,q_in,q_out)
   ! PURPOSE
   !   solve for unperturbed partition function of one step
   !   call FFT and inverse FFT
   ! ARGUMENTS
   !   fft_plan plan  -  see fft_mod for details
   !   integer order  -  order of extrapolation
   !   real q_in      -  input unperturbed partition function
   !   real q_out     -  output unperturbed partition function
   ! SOURCE
   !-----------------------------------------------------------------
   subroutine fft_step(plan,order,q_in,q_out)
   implicit none
   type(fft_plan), intent(IN)   :: plan
   integer,        intent(IN)   :: order
   real(long),     intent(IN)   :: q_in(0:,0:,0:)
   real(long),    intent(OUT)   :: q_out(0:,0:,0:)
   !***    
 
   ! local variables
   complex(long)  :: qk(0:plan%n(1)/2,0:plan%n(2)-1,0:plan%n(3)-1) 
   real(long)     :: qr(0:plan%n(1)-1,0:plan%n(2)-1,0:plan%n(3)-1)
   real(long)     :: r_npoints
   
   r_npoints=dble(plan%n(1)*plan%n(2)*plan%n(3))
   
   qr(:,:,:)=exp_omega(:,:,:,order)*q_in(:,:,:)   ! 1/3 update
   
   call fft(plan,qr,qk)
   qk(:,:,:)=exp_ksq(:,:,:,order)*qk(:,:,:)       ! 2/3 update
   
   call ifft(plan,qk,qr)
   q_out(:,:,:)=exp_omega(:,:,:,order)*qr(:,:,:)/r_npoints ! 3/3 update
 
   end subroutine fft_step
   !====================================================================
 
 
   !--------------------------------------------------------------------
   !****p response_step_mod/fftc_step
   ! SUBROUTINE
   !   fftc_step(plan,order,q_in,q_out)
   ! PURPOSE
   !   solve for perturbed partition function of one step
   !   call FFT and inverse FFT
   ! ARGUMENTS
   !   fft_plan plan  -  see fft_mod for details
   !   integer order  -  extrapolation order
   !   real q_in      -  input partition
   !   real q_out     -  output partition
   ! SOURCE
   !-----------------------------------------------------------------
   subroutine fftc_step(plan,order,q_in,q_out)
   implicit none
   type(fft_plan), intent(IN)   :: plan
   integer,        intent(IN)   :: order
   complex(long),  intent(IN)   :: q_in(0:,0:,0:)
   complex(long), intent(OUT)   :: q_out(0:,0:,0:)
   !***
   !local variables
   complex(long)  :: qk(0:plan%n(1)-1,0:plan%n(2)-1,0:plan%n(3)-1) 
   complex(long)     :: qr(0:plan%n(1)-1,0:plan%n(2)-1,0:plan%n(3)-1)

   qr(:,:,:)=exp_omega(:,:,:,order)*q_in(:,:,:)   ! 1/3 update

   call fftc(1,plan,qr,qk)
   qk(:,:,:)=exp_ksq_new(:,:,:,order)*qk(:,:,:)   ! 2/3 update

   call fftc(-1,plan,qk,qr)
   q_out(:,:,:) = qr(:,:,:)*exp_omega(:,:,:,order)
   
   end subroutine fftc_step
   !====================================================================


   !--------------------------------------------------------------------
   !****p response_step_mod/calc_exp_grid
   ! SUBROUTINE
   !   calc_exp_grid(ksq_grid, ksq_grid_pert, omega_grid, 
   !                 monomer, ds, order, dsarray, pertb)
   ! PURPOSE
   !   Calculate the module variable exponential arrays
   ! ARGUMENTS
   !   real      ksq_grid  -  |G|^2 for reciprocal lattice vectors
   !   real    omega_grid  -  omega fields at grid points
   !   real ksq_grid_pert  -  |k+G|^2 on grids 
   !   integer    monomer  -  monomer type for block
   !   real            ds  -  chain step
   !   integer      order  -  order of extrapolation
   !   real       dsarray  -  output array that has stepsizes as elements
   !   logical      pertb  -  calculate the perturbed exponential if true
   ! SOURCE
   !----------------------------------------------------------------------
   subroutine calc_exp_grid(ksq_grid, ksq_grid_pert, omega_grid, &
                            monomer, ds, order, dsarray, pertb)
   implicit none
 
   real(long),  intent(IN) :: ksq_grid(:,:,:)
   real(long),  intent(IN) :: ksq_grid_pert(:,:,:)
   real(long),  intent(IN) :: omega_grid(:,:,:,:)
   integer,     intent(IN) :: monomer
   real(long),  intent(IN) :: ds
   integer,     intent(IN) :: order
   real(long), intent(out) :: dsarray(:)
   logical,     intent(IN) :: pertb
   !***
 
   real(long)              :: dsloc, bloc
   integer                 :: i
 
   if (order > 2) stop "Error: perturbation calculation not coded for order > 2"
 
   dsloc = ds
   do i = 1, order+2
      dsarray(i) = dsloc
      bloc = kuhn(monomer)
      exp_ksq(:,:,:,i-1)     = exp(-ksq_grid(:,:,:)*dsloc*bloc**2/6.0_long)
      exp_omega(:,:,:,i-1)   = exp(-omega_grid(:,:,:,monomer)*dsloc/2.0_long)
      if ( pertb ) then
         exp_ksq_new(:,:,:,i-1) = exp(-ksq_grid_pert(:,:,:)*dsloc*bloc**2/6.0_long)
      end if
      dsloc = dsloc/2.0_long
   end do
 
   end subroutine calc_exp_grid
   !==========================================================================
  
end module response_step_mod
