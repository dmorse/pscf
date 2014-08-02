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
!****m* scf/string_mod
! MODULE
!   string_mod - string utilities
! PURPOSE
!   Provides character(int_len) function int_string
! SOURCE
!----------------------------------------------------------------------
module string_mod
   implicit none

   private
   public :: int_string  !  returns character representation of integer
   public :: int_len     !  parameter, # characters returned by int_string

   ! Parameter definition
   integer, parameter :: int_len = 10  
   !***

contains

   !-------------------------------------------------------------------
   !****p* string_mod/int_string
   ! FUNCTION
   !   int_string(n,[l]) - returns character representation of integer
   ! RETURN
   !   character(int_len) representation of the integer n,
   !   with no leading white space, padding by trailing spaces
   ! ARGUMENTS
   !   integer n - integer to be converted to a string
   !   integer l - (optional) trimmed length of resulting string
   ! SOURCE
   !-------------------------------------------------------------------
   character(int_len) function int_string(n,l)
   integer, intent(IN)            :: n ! integer input value
   integer, intent(OUT), optional :: l ! trimmed length string
   !***

   integer :: i, j, k, m 
   logical :: negative

   if (n < 0) then
     i        = -n
     negative = .true.
   else
     i        =  n
     negative = .false.
   endif
   int_string = ' '
   m = 0
   j = 1
10 k = mod(i,j*10)
   i = i - k 
   k = k/j
   if ( (m + 1) <= int_len ) then
      int_string = adjustl(int_char(k)//trim(adjustl(int_string)))
      m = m + 1
   else
      write(6,*) 'Error: Argument n=',n,' too large in int_string'
      stop
   endif
   if (i > 0) then
      j = j*10
      go to 10
   endif
   if (negative) then
      if ( (m + 1) <= int_len) then
         int_string = adjustl('-'//trim(adjustl(int_string)))
         m = m + 1
      else
         write(6,*) 'Error: Argument n=',n,' too large in int_string'
         stop
      endif
   endif

   ! Pad with trailing spaces
   if (m < int_len) then
      do k=m+1, int_len
         int_string(k:k) = ' '
      end do
   endif

   if (present(l)) then
      l = m 
   endif
   end function int_string
   !==============================================================


   !-------------------------------------------------------------------
   !****ip string_mod/int_char
   ! FUNCTION
   !   int_char(n)  
   ! RETURN
   !   char*1 representation of integer n=0-9
   ! SOURCE
   !--------------------------------------------------------------
   character(1) function int_char(n)
   integer :: n
   !***
   select case(n)
   case(0)
     int_char = '0'
   case(1)
     int_char = '1'
   case(2)
     int_char = '2'
   case(3)
     int_char = '3'
   case(4)
     int_char = '4'
   case(5)
     int_char = '5'
   case(6)
     int_char = '6'
   case(7)
     int_char = '7'
   case(8)
     int_char = '8'
   case(9)
     int_char = '9'
   case default
     write(6,*) 'Error: Illegal argument n=',n,' in int_char'
   end select 
   end function int_char
   !==============================================================

end module string_mod
