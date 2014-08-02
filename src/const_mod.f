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
!****m* scf/const_mod
! MODULE
!   const_mod
! PURPOSE
!   Define integers variables used by most other modules
! SOURCE
!-----------------------------------------------------------------------
module const_mod 

   public
   integer            :: dim  ! = dimensionality of space
   integer, parameter :: long = selected_real_kind(13)

end module const_mod
!***

!****v const_mod/dim ------------------------------------------
! VARIABLE
!   integer    dim    = dimensionality of space
!*** ----------------------------------------------------------
!****v const_mod/long 
! VARIABLE
!   integer    long   = selected_real_kind(13)
!                     = double precision type kind 
! COMMENT
!   All real variables in scf are declared real(long), which
!   is a portable double precision type.
!*** ----------------------------------------------------------
