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
!****m  scf/step_mod
! MODULE
!   step_mod
! PURPOSE
!   Implements pseudo-spectral algorithm for integration of the modified
!   diffusion equation. The algorithm combines the operator-splitting 
!   method of Rasmussen and Kaloskas with Richardson extrapolation to 
!   obtain an algorithm with errors of O(ds**4). 
!
!   Subroutine init_step must be called once to allocate the FFT arrays
!   used by the module. Subroutine make_propg must be called once at the 
!   beginning of each block of a block copolymer, to set variables that
!   are used throughout that block. Subroutine step_exp implements one
!   'time' step of the integration algorithm. 
! SOURCE
!-----------------------------------------------------------------------
module step_mod
   use const_mod
   use fft_mod
   implicit none

   private

   public :: init_step      ! allocate array needed by module
   public :: make_propg     ! evaluate exp(-ds*b^2*nabla/6), exp(-ds*omega/2)
   public :: step_exp       ! evaluate on integration step 
   !***

   real(long), allocatable    :: exp_omega1(:,:,:)
   real(long), allocatable    :: exp_omega2(:,:,:)
   real(long), allocatable    :: exp_ksq1(:,:,:)
   real(long), allocatable    :: exp_ksq2(:,:,:)
   real(long), allocatable    :: q1(:,:,:)
   real(long), allocatable    :: q2(:,:,:)
   real(long), allocatable    :: qr(:,:,:)
   complex(long), allocatable :: qk(:,:,:)

contains

   !----------------------------------------------------------------
   !****p step_mod/init_step
   ! SUBROUTINE
   !   init_step(n)
   ! PURPOSE
   !   allocate memory to module arrays
   ! ARGUMENTS
   !   integer n(3)  - grid dimensions
   ! SOURCE
   !----------------------------------------------------------------
   subroutine init_step(n)
   implicit none

   integer, intent(IN) :: n(3)
   !***

   integer :: error

   ALLOCATE(exp_omega1(0:n(1)-1,0:n(2)-1,0:n(3)-1), STAT=error)
   if (error /= 0) stop "exp_omega1 allocation error"

   ALLOCATE(exp_omega2(0:n(1)-1,0:n(2)-1,0:n(3)-1), STAT=error)
   if (error /= 0) stop "exp_omega2 allocation error"

   ALLOCATE(exp_ksq1(0:n(1)/2,0:n(2)-1,0:n(3)-1), STAT=error)
   if (error /= 0) stop "exp_ksq1 allocation error"

   ALLOCATE(exp_ksq2(0:n(1)/2,0:n(2)-1,0:n(3)-1), STAT=error)
   if (error /= 0) stop "exp_ksq2 allocation error"

   ALLOCATE(q1(0:n(1)-1,0:n(2)-1,0:n(3)-1), STAT=error)
   if (error /= 0) stop "q1 allocation error"

   ALLOCATE(q2(0:n(1)-1,0:n(2)-1,0:n(3)-1), STAT=error)
   if (error /= 0) stop "q2 allocation error"

   ALLOCATE(qr(0:n(1)-1,0:n(2)-1,0:n(3)-1), STAT=error)
   if (error /= 0) stop "qr allocation error"

   ALLOCATE(qk(0:n(1)/2,0:n(2)-1,0:n(3)-1), STAT=error)
   if (error /= 0) stop "qk allocation error"

   end subroutine init_step
   !================================================================


   !----------------------------------------------------------------
   !****p step_mod/make_propg
   ! SUBROUTINE
   !   make_propg(ds,ksq,omega)
   ! PURPOSE
   !   evaluate exp_ksq, exp_omega's
   ! ARGUMENTS
   !   real(long) ds     - step size
   !   real(long) b      - statistical segment length
   !   real(long) ksq    - k^2 on grid
   !   real(long) omega  - omega field on grid
   ! SOURCE
   !----------------------------------------------------------------
   subroutine make_propg(ds,b,ksq,omega)
   implicit none

   real(long), intent(IN) :: ds
   real(long), intent(IN) :: b
   real(long), intent(IN) :: ksq(0:,0:,0:)
   real(long), intent(IN) :: omega(0:,0:,0:)
   !***

   real(long)  :: lap_coeff, pot_coeff

   lap_coeff = ds * b**2 /  6.0_long
   pot_coeff = ds / 2.0_long

   exp_ksq1 = exp( - lap_coeff * ksq )
   exp_ksq2 = exp( - lap_coeff / 2.0_long * ksq )

   exp_omega1 = exp( - pot_coeff * omega )
   exp_omega2 = exp( - pot_coeff / 2.0_long * omega )

   end subroutine make_propg
   !================================================================


   !----------------------------------------------------------------
   !****p step_mod/step_exp
   ! SUBROUTINE
   !   step_exp
   ! PURPOSE
   !   Calculate one step in pseudo-spectral algorithm for 
   !   integrating the modified diffusion equation.
   ! ARGUMENTS
   !   real q_in      -  input  q(r,s)
   !   real q_out     -  output q(r,s+-ds)
   !   fft_plan plan  -  see fft_mod for details
   ! SOURCE
   !----------------------------------------------------------------
   subroutine step_exp(q_in, q_out, plan)
   implicit none

   real(long), intent(IN)     :: q_in(0:,0:,0:)
   real(long), intent(OUT)    :: q_out(0:,0:,0:)
   type(fft_plan), intent(IN) :: plan
   !***

   ! local variables
   real(long)     :: r_npoints
   r_npoints = dble(plan%n(1)*plan%n(2)*plan%n(3))

   qr=exp_omega1*q_in            ! qr(r) = exp{-omega(r)*ds/2)*q_in(r)
   call fft(plan,qr,qk)          ! qk    = F[qr]
   qk=exp_ksq1*qk                ! qk(k) = exp(-ds*(k*b)**2/6)*qk(k)
   call ifft(plan,qk,qr)         ! qr    = F^{-1}[qk]
   q1=exp_omega1*qr/r_npoints    ! q1    = exp^{-omega(r)*ds/2)*qr(r)

   qr=exp_omega2*q_in            ! qr(r) = exp^{-omega(r)*ds/4)*q_in(r)
   call fft(plan,qr,qk)          ! qk    = F[qr]
   qk=exp_ksq2*qk                ! qk(k) = exp(-ds*(k*b)**2/12)*qk(k)
   call ifft(plan,qk,qr)         ! qr    = F^{-1}[qk]
   qr=exp_omega1*qr/r_npoints    ! q2    = exp^{-omega(r)*ds/2)*qr(r)
   call fft(plan,qr,qk)          ! qk    = F[qr]
   qk=exp_ksq2*qk                ! qk(k) = exp(-ds*(k*b)**2/12)*qk(k)
   call ifft(plan,qk,qr)         ! qr    = F^{-1}[qk]
   q2=exp_omega2*qr/r_npoints    ! q2    = exp^{-omega(r)*ds/4)*qr(r)

   q_out = ( 4.0_long * q2 - q1 ) / 3.0_long

   end subroutine step_exp
   !================================================================

end module step_mod
