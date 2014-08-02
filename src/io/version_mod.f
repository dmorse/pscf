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
!****m* scf/version_mod
! MODULE
!   version_mod 
! PURPOSE
!   Defines derived type version_type, which represents a file format
!   or program version number, with a major version and minor version.
!   Defines several associated subroutines and functions
! SOURCE
!-----------------------------------------------------------------------
module version_mod
implicit none

   private
   public :: version_type
   public :: input_version, output_version
   public :: version_is, version_ge
   !***

   !-------------------------------------------------------------------
   !****t* version_mod/version_type
   ! TYPE
   !    version_type
   ! PURPOSE
   !    Represents a file format or program version number
   ! VARIABLE
   !    major - major version number 
   !    minor - minor version number 
   ! SOURCE
   !-------------------------------------------------------------------
   type version_type
      integer :: major
      integer :: minor
   end type version_type
   !***

contains

   !-------------------------------------------------------------------
   !****p* version_mod/input_version
   ! SUBROUTINE
   !    input version
   ! PURPOSE
   !    Read major and minor numbers of format from 1st line of file
   ! ARGUMENTS
   !    version - version_type version object
   !    unit    - file unit number for reading
   ! SOURCE
   !-------------------------------------------------------------------
   subroutine input_version(version, unit)
   type(version_type), intent(OUT)  :: version
   integer, intent(IN)              :: unit
   character(80) comment
   read(unit,*) comment, version%major, version%minor
   !***
   end subroutine input_version
   !===================================================================


   !-------------------------------------------------------------------
   !****p* version_mod/output_version
   ! SUBROUTINE
   !    output_version
   ! PURPOSE
   !    Read major and minor numbers of format to 1st line of file
   ! ARGUMENTS
   !    version - version_type version object
   !    unit    - file unit number for reading
   ! SOURCE
   !-------------------------------------------------------------------
   subroutine output_version(version, unit)
   type(version_type), intent(IN)  :: version
   integer, intent(IN)             :: unit
   write(unit,FMT='(A7,I3,I3)') 'format', version%major, version%minor
   !***
   end subroutine output_version
   !===================================================================


   !-------------------------------------------------------------------
   !****p* version_mod/version_is
   ! FUNCTION
   !    version_is(version, major, minor)
   ! PURPOSE
   !    Returns true if version has specified value, false otherwise
   ! ARGUMENTS
   !    version - version_type version object
   !    major   - major version number for comparison
   !    minor   - minor version number for comparison
   ! SOURCE
   !-------------------------------------------------------------------
   logical function version_is(version, major, minor)
   type(version_type), intent(IN)  :: version
   integer, intent(IN)             :: major, minor

   if ( version%major == major .and. version%minor == minor ) then
      version_is  = .TRUE.
   else
      version_is  = .FALSE.
   endif
   !***
   end function version_is
   !===================================================================


   !-------------------------------------------------------------------
   !****p* version_mod/version_ge
   ! FUNCTION
   !    version_ge(version, major, minor)
   ! PURPOSE
   !    Returns true if version is greater than or equal to a 
   !    specified value, false otherwise.
   ! ARGUMENTS
   !    version - version_type version object
   !    major   - major version number for comparison
   !    minor   - minor version number for comparison
   ! SOURCE
   !-------------------------------------------------------------------
   logical function version_ge(version, major, minor)
   type(version_type), intent(IN)  :: version
   integer, intent(IN)             :: major, minor
   if ((version%major >= major).and.(version%minor >= minor)) then
      version_ge  = .TRUE.
   else
      version_ge  = .FALSE.
   endif
   !***
   end function version_ge
   !===================================================================

end module version_mod
