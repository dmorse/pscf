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
!    Fortran 90 wrapper subroutines for fftw2 package
! SOURCE
!-----------------------------------------------------------------------
module fft_mod
   use const_mod
   implicit none

   PRIVATE 
   PUBLIC :: fft_plan
   PUBLIC :: create_fft_plan    ! Routine to create an fft_plan
   PUBLIC :: destroy_fft_plan   ! Routine to destroy an fft_plan
   PUBLIC :: fft                ! Forward FFT for 1, 2, or 3D
   PUBLIC :: ifft               ! Inverse FFT for 1, 2, or 3D
   PUBLIC :: fftc               ! Complex FFT for 1, 2, or 3D
   !***

   ! Parameters required by fftw
   integer, parameter :: FFTW_FORWARD=-1 ,FFTW_BACKWARD=1
   integer, parameter :: FFTW_REAL_TO_COMPLEX=-1, FFTW_COMPLEX_TO_REAL=1
   integer, parameter :: FFTW_ESTIMATE=0, FFTW_MEASURE=1
   integer, parameter :: FFTW_OUT_OF_PLACE=0, FFTW_IN_PLACE=8
   integer, parameter :: FFTW_USE_WISDOM=16
   integer, parameter :: FFTW_THREADSAFE=128

   !-------------------------------------------------------------------
   !****t fft_mod/fft_plan
   ! TYPE
   !    fft_plan 
   ! PURPOSE
   !    Contains dimensions for FFT grid, and integer pointers to the 
   !    "plan" structures used by the FFTW package
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
   if ( present(fft_c2c) .and. (fft_c2c .eq. .true.) ) then

      call fftwnd_f77_create_plan_(plan%f,dim,plan%n(1:dim),&
         FFTW_FORWARD,FFTW_ESTIMATE)
      call fftwnd_f77_create_plan_(plan%r,dim,plan%n(1:dim),&
         FFTW_BACKWARD,FFTW_ESTIMATE)

   else

      call rfftwnd_f77_create_plan_(plan%f,dim,plan%n(1:dim),&
                        FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE)
      call rfftwnd_f77_create_plan_(plan%r,dim,plan%n(1:dim),&
                        FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE)
   end if
   end subroutine create_fft_plan
   !===================================================================


   !-------------------------------------------------------------------
   !****p fft_mod/destroy_fft_plan
   ! SUBROUTINE
   !    destroy_fft_plan
   ! PURPOSE
   !    destroys argument plan of type fft_plan
   !-------------------------------------------------------------------
   subroutine destroy_fft_plan(plan)
   type(fft_plan),intent(IN) :: plan
   !***
   call rfftwnd_f77_destroy_plan_(plan%f)
   call rfftwnd_f77_destroy_plan_(plan%r)
   end subroutine destroy_fft_plan
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
   select case(dim)
   case (1)
      call rfftwnd_f77_real_to_complex_(plan%f,1,&!
                    in(:,0,0),1,0,out(:,0,0),1,0)
   case (2)
      call rfftwnd_f77_real_to_complex_(plan%f,1,&!
                    in(:,:,0),1,0,out(:,:,0),1,0)
   case (3)
      call rfftwnd_f77_real_to_complex_(plan%f,1,&!
                    in,1,0,out,1,0)
   case default
      STOP "dim < 1 or > 3 when performing FFT"
   end select 
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
   select case(dim)
   case (1)
      call rfftwnd_f77_complex_to_real_(plan%r,1,&!
                    in(:,0,0),1,0,out(:,0,0),1,0)
   case (2)
      call rfftwnd_f77_complex_to_real_(plan%r,1,&!
                    in(:,:,0),1,0,out(:,:,0),1,0)
   case (3)
      call rfftwnd_f77_complex_to_real_(plan%r,1,&!
                    in,1,0,out,1,0)
   case default
      STOP "dim < 1 or > 3 when IFFT"
   end select 
   end subroutine ifft
   !===================================================================

   !-------------------------------------------------------------------
   !****p fft_mod/fftc
   ! SUBROUTINE
   !     fftc(plan,in,out)
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
   ! SOURCE
   !-------------------------------------------------------------------
   subroutine fftc(direction,plan,in,out)
   integer,intent(IN)          :: direction
   type(fft_plan),intent(IN)   :: plan
   complex(long), intent(IN)   :: in(0:,0:,0:)
   complex(long), intent(OUT)  :: out(0:,0:,0:)
   !***
   integer*8                   :: lplan
   if (direction == +1) lplan = plan%f
   if (direction == -1) lplan = plan%r
   select case(dim)
   case (1)
      call fftwnd_f77_one_(lplan, in(:,0,0), out(:,0,0) )
   case (2)
      call fftwnd_f77_one_(lplan,in(:,:,0), out(:,:,0) )
   case (3)
      call fftwnd_f77_one_(lplan, in, out)
   case default
      write(6,*) "dimension = ",dim," incorrect!"
      STOP
   end select
 
   if (direction == +1) out = out / dble( plan%n(1) * plan%n(2) * plan%n(3) )
 
   end subroutine fftc
  !===================================================================


end module fft_mod
