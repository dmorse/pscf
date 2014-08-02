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
!****m scf/extrapolate_mod
! PURPOSE
!    Richardson extrapolation - extrapolate to ds = 0
! SOURCE
!----------------------------------------------------------------------
module extrapolate_mod
   use const_mod
   implicit none
 
   private
   public :: extrapolate_real
   public :: extrapolate_complex
   !***  

contains

   !-----------------------------------------------------------------
   !****p extrapolate_mod/extrapolate_real
   ! SUBROUTINE
   !    extrapolate_real(order,xarray,yarray,extrap_yvalue)
   ! PURPOSE
   !    Calculate the extrapolated value of a real array
   ! ARGUMENTS
   !    integer order         -  order of extrapolation
   !    real    xarray        -  the array of x 
   !    real    yarray        -  the array of y that depends on x
   !    real    extrap_value  -  the y value extrapolated to x = 0 case
   ! SOURCE
   !-----------------------------------------------------------------
   subroutine extrapolate_real(order,xarray,yarray,extrap_yvalue)
   integer,intent(IN)     :: order
   real(long),intent(IN)  :: xarray(order+1)
   real(long),intent(IN)  :: yarray(order+1)
   real(long),intent(out) :: extrap_yvalue
   !***  
   
   integer                :: i
   integer                :: j
   real(long)             :: num
   real(long)             :: den
   real(long)             :: a

   !----------------------------------------------------------    
   ! Using Lagrange's formula for an nth order polynomial that 
   ! passes through n points, for extrapolation
   !----------------------------------------------------------    
   !    extrap_yvalue = 0.0_long        
   !    do i = 1,order+1
   !       num = yarray(i)
   !       a   = xarray(i)
   !       den = 1.0_long
   !       do j = 1,order+1
   !          if (j /= i) then
   !             num = num*xarray(j) * (-1.0_long)
   !             den = den *(a-xarray(j))
   !          end if
   !       end do
   !       extrap_yvalue = extrap_yvalue + (num/den)
   !    end do
   !----------------------------------------------------------    
    
   if (order == 1) then
      extrap_yvalue = 4.0_long*yarray(2)-yarray(1)
      extrap_yvalue = extrap_yvalue/3.0_long
   end if
   if (order == 2) then
      ! extrap_yvalue = -6*yarray(1)
      ! extrap_yvalue = extrap_yvalue + 80.0_long*yarray(2)
      ! extrap_yvalue = extrap_yvalue + 256.0_long*yarray(3)
      ! extrap_yvalue = extrap_yvalue/330.0_long
      extrap_yvalue = 64.0_long*yarray(3)
      extrap_yvalue = extrap_yvalue - 20.0_long*yarray(2)
      extrap_yvalue = extrap_yvalue + yarray(1)
      extrap_yvalue = extrap_yvalue/45.0_long
   end if
    
   end subroutine extrapolate_real
  
   !-----------------------------------------------------------------
   !****p extrapolate_mod/extrapolate_complex
   ! SUBROUTINE
   !    extrapolate_complex(order,xarray,yarray,extrap_yvalue)
   ! PURPOSE
   !    Calculate the extrapolated value of a complex array
   ! ARGUMENTS
   !    integer    order         -  order of extrapolation
   !    complex    xarray        -  the array of x 
   !    complex    yarray        -  the array of y that depends on x
   !    complex    extrap_value  -  the y value extrapolated to x = 0 case
   ! SOURCE
   !-----------------------------------------------------------------
   subroutine extrapolate_complex(order,xarray,yarray,extrap_yvalue)
   integer,intent(IN)     :: order
   real(long),intent(IN)  :: xarray(order+1)
   complex(long),intent(IN)  :: yarray(order+1)
   complex(long),intent(out) :: extrap_yvalue
   !***
   integer                :: i
   integer                :: j
   complex(long)          :: num
   complex(long)          :: den
   real(long)             :: a

   !----------------------------------------------------------    
   ! Using Lagrange's formula for an nth order polynomial that 
   ! passes through n points, for extrapolation
   !----------------------------------------------------------        
   !    extrap_yvalue = (0.0_long,0.0_long)    
   !    do i = 1,order+1
   !       num = yarray(i)
   !       a   = xarray(i)
   !       den = 1.0_long
   !       do j = 1,order+1
   !          if (j /= i) then
   !             num = num*xarray(j) * (-1.0_long)
   !             den = den *(a-xarray(j))
   !          end if
   !       end do
   !       extrap_yvalue = extrap_yvalue + (num/den)
   !    end do
   !----------------------------------------------------------      

   if (order == 1) then
      extrap_yvalue = 4.0_long*yarray(2)-yarray(1)
      extrap_yvalue = extrap_yvalue/3.0_long
   end if
   if (order == 2) then

      ! extrap_yvalue = -6*yarray(1)
      ! extrap_yvalue = extrap_yvalue + 80.0_long*yarray(2)
      ! extrap_yvalue = extrap_yvalue + 256.0_long*yarray(3)
      ! extrap_yvalue = extrap_yvalue/330.0_long

      extrap_yvalue = 64.0_long*yarray(3)
      extrap_yvalue = extrap_yvalue - 20.0_long*yarray(2)
      extrap_yvalue = extrap_yvalue + yarray(1)
      extrap_yvalue = extrap_yvalue/45.0_long

   end if
    
   end subroutine extrapolate_complex
   !===================================================================
  
end module extrapolate_mod
