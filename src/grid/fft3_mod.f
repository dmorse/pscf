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
!****m scf/fft_mod
! PURPOSE
!    Derived type fft_plan
!    Wrapper subroutines for fftw-3
! COMMENTS
!    Consider use even/odd fftw, and r-to-r fftw in the future
! SOURCE
!-----------------------------------------------------------------------
module fft_mod
   use const_mod
   implicit none

   PRIVATE 
   PUBLIC :: fft_plan
   PUBLIC :: create_fft_plan    ! initialize an fft_plan
   PUBLIC :: fftc               ! complex FFT for 1, 2, or 3D
   PUBLIC :: fft                ! Forward FFT for 1, 2, or 3D
   PUBLIC :: ifft               ! Inverse FFT for 1, 2, or 3D
   !***

   ! Parameters required by fftw3 supplied in fftw3.f
   integer, parameter :: FFTW_ESTIMATE=64
   integer, parameter :: FFTW_FORWARD=-1 ,FFTW_BACKWARD=1

   !-------------------------------------------------------------------
   !****t fft_mod/fft_plan
   ! TYPE
   !    fft_plan 
   ! PURPOSE
   !    Contains grid dimensions for FFT grid and integer pointers to
   !    the "plan" structures used by the FFTW package
   ! SOURCE
   !-------------------------------------------------------------------
   type fft_plan
      integer    ::  n(3)   ! grid dimensions, 0 for unused dimensions
      integer*8  ::  f      ! fftw plan object for forward transform
      integer*8  ::  r      ! fftw plan object for inverse transform
   end type fft_plan
   !***

contains

   !-------------------------------------------------------------------
   !****p fft_mod/create_fft_plan
   ! SUBROUTINE
   !    create_fft_plan
   ! PURPOSE
   !    Creates an fft_plan object for grids with dimensions 
   !    ngrid(1),..,ngrid(dim)
   !-------------------------------------------------------------------
   subroutine create_fft_plan(ngrid,plan,fft_c2c)
   integer,intent(IN)             :: ngrid(3) ! dimensions of grid
   type(fft_plan),intent(OUT)     :: plan
   logical, optional, intent(IN)  :: fft_c2c
   !***
   plan%n=ngrid
   end subroutine create_fft_plan
   !===================================================================


   !-------------------------------------------------------------------
   !****p fft_mod/fft
   ! SUBROUTINE
   !     fft(plan,in,out)
   ! PURPOSE
   !    Calculates forward fft of in, returns result in out.
   !    Wrapper for 1, 2, & 3 dimensional real -> complex transforms
   ! ARGUMENTS
   !    plan    - fft plan object
   !    in, out - real(long) 3D arrays 
   ! COMMENT
   !    in and out are dimensioned 0:ngrid(i)-1 for all i <= dim, 
   !    and 0:0 for any unused dimensions with dim < i <= 3
   !-------------------------------------------------------------------
   subroutine fft(plan,in,out)
   type(fft_plan),intent(IN)   :: plan
   real(long), intent(IN)      :: in(0:,0:,0:)
   complex(long), intent(OUT)  :: out(0:,0:,0:)
   !***
   call dfftw_plan_dft_r2c_3d(plan%f,plan%n(1),plan%n(2),plan%n(3),&!
                              in,out,FFTW_ESTIMATE)
   call dfftw_execute(plan%f)
   call dfftw_destroy_plan(plan%f)
   end subroutine fft
   !===================================================================


   !-------------------------------------------------------------------
   !****p fft_mod/ifft
   ! SUBROUTINE
   !     ifft(plan,in,out)
   ! PURPOSE
   !    Calculates inverse fft of real array in, returns in out.
   !    Wrapper for 1, 2, & 3 dimensional complex -> real transforms
   ! ARGUMENTS
   !    plan - fft plan object
   !    in   - complex(long) 3D input array
   !    out  - real(long) 3D input array
   ! COMMENT
   !    in and out are dimensioned 0:ngrid(i)-1 for all i <= dim, 
   !    and 0:0 for any unused dimensions with dim < i <= 3
   !-------------------------------------------------------------------
   subroutine ifft(plan,in,out)
   type(fft_plan),intent(IN)   :: plan
   complex(long), intent(IN)   :: in(0:,0:,0:)
   real(long), intent(OUT)     :: out(0:,0:,0:)
   !***
   call dfftw_plan_dft_c2r_3d(plan%r,plan%n(1),plan%n(2),plan%n(3),&!
                              in,out,FFTW_ESTIMATE)
   call dfftw_execute(plan%r)
   call dfftw_destroy_plan(plan%r)
   end subroutine ifft
   !===================================================================


   !-------------------------------------------------------------------
   !****p fft_mod/fftc
   ! SUBROUTINE
   !     fftc(direction,plan,in,out)
   ! PURPOSE
   !    Calculates forward fft of in, returns result in out.
   !    Wrapper for 1, 2, & 3 dimensional real -> complex transforms
   ! ARGUMENTS
   !    plan    - fft plan object
   !    in, out - real(long) 3D arrays 
   ! COMMENT
   !    in and out are dimensioned 0:ngrid(i)-1 for all i <= dim, 
   !    and 0:0 for any unused dimensions with dim < i <= 3
   !-------------------------------------------------------------------
   subroutine fftc(direction,plan,in,out)
   integer,intent(IN)          :: direction
   type(fft_plan),intent(IN)   :: plan
   complex(long), intent(IN)   :: in(0:,0:,0:)
   complex(long), intent(OUT)  :: out(0:,0:,0:)
   !***
   if (direction == 1) then
      call dfftw_plan_dft_3d(plan%f,plan%n(1),plan%n(2),plan%n(3),&!
             in,out,FFTW_FORWARD,FFTW_ESTIMATE)
   else
      call dfftw_plan_dft_3d(plan%f,plan%n(1),plan%n(2),plan%n(3),&!
             in,out,FFTW_BACKWARD,FFTW_ESTIMATE)
   end if
   call dfftw_execute(plan%f)
   call dfftw_destroy_plan(plan%f)

   if (direction == +1) out = out/dcmplx( dble(plan%n(1)*plan%n(2)*plan%n(3)) , 0.0_long)
   end subroutine fftc
   !===================================================================


end module fft_mod
