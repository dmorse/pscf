! fortran_dialect=elf
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
!****m* scf/io_mod
! MODULE
!   io_mod - generic subroutine interfaces for file io
!
! PURPOSE
!
!   The generic subroutine interfaces 'input' and 'output' provide 
!   a standard interface for reading parameters from and writing 
!   to file. Similar interfaces are used for integer, real, logical,
!   and character(*) data, for 1D arrays of integers, real, or
!   character(*) data, and 2D arrays of integers or real data.
!
!   Input and output styles can be chosen from several that allow
!   each data item to be accompanied by a comment string, namely:
!   comment on line above data, comment before or after data on the
!   same line, or no comment. Subroutines with names beginning with
!   set_... allow the user to set values of private module variables
!   that determine default input and output formats. These default 
!   values may also be overridden by the use of optional arguments 
!   to input and output routines
!
! PUBLIC PROCEDURES
!
!   input  - generic interface for input subroutines 
!   output - generic interface for output subroutine  
!
!   set_io_units   - set input and output file unit numbers
!   set_com_style  - set comment style (same for input and output)
!   set_output_fmt - set output format (field widths, float format)
!   set_echo       - choose whether to echo input
!   set_com_use    - choose whether to echo input comments
!
! AUTHOR
!   David Morse (12/2003)
!
!*** -------------------------------------------------------------------
!****p* io_mod/input
! SUBROUTINE
!
!   input  - generic interface for reading data from file
!   output - generic interface for writing data to file
!
! SYNOPSIS
!
!   Scalar integer, real(long), character(*) data
!
!       call input(data,[c,i,o,e,f,u])
!       call output(data,[c,o,e,f])     ! integer, real, and logical
!       call output(data,[c,o,e,f,l])   ! character(*)
!
!       data = variable to be read and/or written to file
!              integer, real(long) 
!
!   Vectors (1D arrays) integer, real, or character data:
!
!       input(data,n,[c,s,i,o,e,f,u])
!       output(data,n,[c,s,o,e,f])      ! integer or real
!       output(data,n,[c,o,e,f,l])      ! character(*)
!
!       data(:) = 1D array to be to be read from file
!                 integer or real(long) 
!
!       n       = logical dimension of data(1:n)
!                 integer
!
!   Matrices (2D arrays) of integer or real data (input & output)
!
!       input(data,m,n,[c,s,i,o,e,f,u])
!       output(data,m,n,[c,s,o,e,f])
!
!       data(:,:) = 2D array to be to be read from file
!                   integer or real(long) 
!
!       m, n      = logical dimensions of data(1:m,1:n)
!                   integers
!
! ARGUMENTS
!       
!       data = scalar, vector, or matrix data, as discussed above
!
!       c    = comment string for output (or echoed input)
!              character(*) (optional)
! 
!       i, o = input and output file unit numbers
!              integer (optional)
! 
!       e    = echo flag: e=1 echo input, e=0 no echo
!              integer (optional)
!
!       f    = comment style flag (see discussion)
!              character(1) (optional)
!
!       u    = comment usage flag (see discussion)
!              character(1) (optional) - only for input
!
!       l    = field width for character string output
!              integer(optional) - only for character output
!
!       s    = vector or matrix io format flag - see discussion below
!              character(1) (optional) - only for vector or matrix data
!
!***
!---------------------------------------------------------------------
!****p* io_mod/output
! SUBROUTINE
!
!   output - generic interface for writing data to file
!
! SYNOPSIS
!
!   Scalar integer, real(long), character(*) data
!
!       call output(data,[c,o,e,f])     ! integer, real, and logical
!       call output(data,[c,o,e,f,l])   ! character(*)
!
!   Vectors (1D arrays) integer, real, or character data:
!
!       output(data,n,[c,s,o,e,f])      ! integer or real
!       output(data,n,[c,o,e,f,l])      ! character(*)
!
!   Matrices (2D arrays) of integer or real data (input & output)
!
!       output(data,m,n,[c,s,o,e,f])
!
! ARGUMENTS
!  
!   See arguments of subroutine input
!
!***
!-------------------------------------------------------------------
!****c io_mod/io_format
! COMMENT
!
!   Comment styles:
!
!   The module defines four styles in which comments strings
!   can be associated with inputs and output data. These 
!   styles are associated with four possible values of a 
!   character(1) variable:
!
!   'N' -> None  - no comment 
!   'A' -> Above - comment on separate line above data
!   'L' -> Left  - comment to the left of data on the same line
!   'R' -> Right - comment to the right of data on the same line
!
!   The 'N' and 'A' style are defined for any type of data,
!   while only 'N' and 'A' are defined for matrix data.
!
!   A default comment style is given by the value of one the 
!   three global character(1) variables, which must have one of 
!   above values: 
!
!      scalar_com_style for scalar   (int, real, or char) 
!      vector_com_style for vectors  (int or real)
!      matrix_com_style for matrices (int or real)
!
!   These default comment style variables may modified by 
!   a subroutine call:
!   
!      call set_com_style([s],[v],[m])
!
!   in which optional character(1) variables s, v, and m hold
!   the desired values of the scalar, vector, and matrix 
!   comment styles, respectively. 
!
!   The default comment style may be overridden by passing 
!   input or output the optional character(1) argument f with 
!   one of the above values.
!
!-----------------------------------------------------------------
! Input format:
!   
!   All input routines use the default format read(iunit,*) for
!   both comments and data. As a result of this:
!
!   1) The spacing of data and comments within an input record 
!      is irrelevant.
!
!   2) Input comment strings must consist of single character(*)
!      tokens, i.e., they must either be strings with no spaces or 
!      other delimiters, or surrounded by quotation marks. 
! 
!   Rule (2) must be obeyed in the 'L' comment style in order to
!   to allow the comment to be distinguished from the input data.
!   In the 'A' and 'R' styles, if a comment string contains spaces, 
!   the data is read corectly, but only the first word of comment 
!   is actually read. 
!
!-----------------------------------------------------------------
! Output formats:
!
!   The output format is controlled by the following global
!   variables, whose values may be modified by calling 
!   set_output_format:
!
!   com_width  = width of comment output field        (integer)
!   data_width = width of data output field           (integer)
!   frac_width = # digits after decimal for reals #s  (integer)
!   fmt_ef     = format for reals = 'E','F', or 'ES'  (character*2)
!
!   The field width data_width is used for scalar integer, 
!   real, and character data.  Initial values are given in 
!   the declarations of these variables.
!
!   Row vectors are output on a single record by repeating the 
!   format string for the corresponding scalar, preceded by a
!   comment line in the 'A' comment style.
!  
!   Matrices are output as a sequence of row vectors, 
!   preceded by a comment line in the 'A' comment style.
!   If the symmetry flag s='A' or s='L', then only the
!   lower diagonal or below diagonal sector that is read
!   on input is output (see below)
!   
!-----------------------------------------------------------------
! Vector io format flag:
!
!   s = 'R'    Row vector - all data in one record
!              (default if argument s is absent)
!
!              Format (for comment style 'N' or 'A')
!
!                   data(1)  data(2) data(3) ... data(N)
!
!   s = 'C'    Column vector -each element on a separate line
!
!              Format (for comment style 'N' or 'A')
!
!                   data(1)  
!                   data(2) 
!                    ...
!                   data(N)
!
!-----------------------------------------------------------------
! Matrix io format flag:
!
!   s = 'N'    Normal or No symmetry
!              (default if argument s is absent)
!
!              Format for 3 x 3 matrix
!
!                data(1,1) data(1,2) data(1,3) 
!                data(2,1) data(2,2) data(2,3)
!                data(3,1) data(3,2) data(3,3)
!
!   s = 'S'    Symmetric matrix
!              data(i,j) = data(j,i)
!              read only data(i,j) for j <= i, i=1,..,m
!
!              Format for 3 x 3 matrix:
!
!                data(1,1) 
!                data(2,1) data(2,2) 
!                data(3,1) data(3,2) data(3,3)
!
!
!   s = 'L'    Symmetric matrix with zero diagonal elements
!              data(i,j) = data(j,i)
!              data(i,i) = 0
!              read only data(i,j) for j < i, i=2,..,m
!
!              Format for 3 x 3 matrix:
!
!                data(2,1) 
!                data(3,1) data(3,2) 
!
!-----------------------------------------------------------------
! Echoing: 
!
!   Input data and a comment are printed to an output file if:
!
!   1) The global variable default_echo has the value 
!      default_echo = 1, and the optional argument e is absent
!
!   2) the optional argument e is present, and e=1.
!
!   Echoed output is output using the same comment style as that
!   used for input. 
!
!   The default_echo variable may be reset by the subroutine call
!
!      call set_echo(e)
!
!   where e is the desired integer default value (e=1 for echoing, 
!   e=0 for no echoing)
!
!-----------------------------------------------------------------
! Comment Usage:
!
!   The module defines three possible treatments of the
!   comments that are read from the input file with data,
!   which are associated with four possible values of a
!   character(1) com_use variable:
!
!   'K' -> Keep    - Return the input comment as argument c 
!                    on output, if c is present, and use it 
!                    in any echoed output.
!
!   'R' -> Replace - Replace the input comment by the input 
!                    value of argument c in any echoed output.
!
!   'C' -> Check   - Check that the input comment matches the
!                    input value of argument c, write error
!                    message if they do not match (not yet
!                    implemented)
!
!   The choice of one of these actions is determined either by
!   by the value of a global argument*1 variable default_com_use,
!   or the value of the character(1) argument 'u', if present,
!   both of which must take on one of the above values.
!
!***
!-------------------------------------------------------------------
module io_mod 
   use const_mod
   use string_mod
   implicit none

   PRIVATE
   PUBLIC :: input, output 
   PUBLIC :: set_io_units, set_com_style, set_output_fmt
   PUBLIC :: set_echo, set_com_use


   !-------------------------------------------------------------
   ! Global variables - declare and set initial values 

   ! Default input and output file units :
   integer     :: default_iunit = 5 ! input file unit
   integer     :: default_ounit = 6 ! output file unit

   ! Default echo flag (echo 1, don't if 0)
   integer     :: default_echo  = 1

   ! Default comment styles for scalar, vector, and matrix data:
   character(1) :: scalar_com_style = 'A' ! data comment
   character(1) :: vector_com_style = 'A' ! comment / data
   character(1) :: matrix_com_style = 'A' ! comment / data

   ! Default com_use flag
   character(1) :: default_com_use = 'R'  ! replace

   ! Default vector and matrix symmetry flags:
   character(1) :: default_vec_sym = 'R'  ! row vector
   character(1) :: default_mat_sym = 'N'  ! normal (full matrix)

   ! Integer field widths for output of comments and data:
   integer :: com_width   = 20  ! comment field width
   integer :: data_width  = 20  ! data field width (total)
   integer :: frac_width  = 12  ! # digits after decimal

   ! Choice of format descriptor (E or F) for real numbers
   character(2) :: fmt_ef = 'ES'     ! must = F, E, or ES

   ! Output format strings for scalar variables:
   character(3) :: fmt_c  = 'A20'     ! comment format
   character(3) :: fmt_i  = 'I20'     ! integer format
   character(7) :: fmt_r  = 'ES20.10' ! real format
   character(3) :: fmt_l  = 'L20'     ! logical format

   ! The initial values of the format strings 
   ! fmt_c, fmt_i, fmt_r fmt_l, may be modified, but should
   ! be consistent with initial values of the integer field 
   ! widths com_width , data_width, frac_width. Also, the
   ! floating point io descriptor E, F, or ES in fmt_r should
   ! be that given in the variable fmt_ef. Calls to the
   ! subroutine set_output_format allow the integer field 
   ! widths and character format strings to be reset in
   ! a consistent manner.

   ! Global character length parameters (used in declarations)
   integer, parameter :: com_len  = 50  ! comment variables
   integer, parameter :: fmt_len  = 25  ! io format strings

   ! End global variables and parameters
   !-------------------------------------------------------------

   interface input
      module procedure input_int
      module procedure input_real
      module procedure input_char
      module procedure input_logic
      module procedure input_int_vec
      module procedure input_real_vec
      module procedure input_char_vec
      module procedure input_int_mat
      module procedure input_real_mat
   end interface

   interface output
      module procedure output_int
      module procedure output_real
      module procedure output_char
      module procedure output_logic
      module procedure output_int_vec
      module procedure output_real_vec
      module procedure output_char_vec
      module procedure output_int_mat
      module procedure output_real_mat
   end interface

contains


   !------------------------------------------------------------
   !****p* io_mod/set_io_units
   ! SUBROUTINE
   !   set_io_units   - set input and output file unit numbers
   ! SOURCE
   !-------------------------------------------------------------------
   subroutine set_io_units(i,o)
   integer, intent(IN), optional :: i  ! default input unit #
   integer, intent(IN), optional :: o  ! default output unit #
   !***

   if (present(i)) then
      default_iunit = i
   endif
   if (present(o)) then
      default_ounit = o
   endif
   end subroutine set_io_units
   !-----------------------------------------------------------


   !-------------------------------------------------------------------
   !****p* io_mod/set_echo
   ! SUBROUTINE
   !   set_echo(e)    - choose whether to echo input
   ! SOURCE
   !-------------------------------------------------------------------
   subroutine set_echo(e)
   integer, intent(IN) :: e
   !***
   select case(e)
   case(0,1)
      default_echo = e
   case default
      write(6,*) 'Error: Incorrect argument',e,'to set_echo'
   end select
   end subroutine set_echo
   !------------------------------------------------------------------


   !------------------------------------------------------------------
   !****p* io_mod/set_com_style
   ! SUBROUTINE
   !   set_com_style(s,v,m)  - set comment style for input and output
   ! PURPOSE
   !    Allows new values to be set for the default comment
   !    styles for scalars, vectors, and matrices. All 
   !    arguments are optional character(1) variables which
   !    must (if present) must have a valid value 'N', 'A',
   !    'L' or 'R for a comment 
   ! SOURCE
   !------------------------------------------------------------------
   subroutine set_com_style(s,v,m)
   character(1), intent(IN), optional :: s  ! scalar_com_style
   character(1), intent(IN), optional :: v  ! vector_com_style
   character(1), intent(IN), optional :: m  ! matrix_com_style
   !***

   if (present(s)) then
      if ((s=='N').or.(s=='A').or.(s=='L').or.(s=='R')) then
         scalar_com_style = s
      else
         write(6,*)'Error: Illegal scalar_com_style =', s
      endif 
   endif
   if (present(v)) then
      if ((v=='N').or.(v=='A').or.(v=='L').or.(v=='R')) then
         vector_com_style = v
      else
         write(6,*)'Error: Illegal vector_com_style =', v
      endif 
   endif
   if (present(m)) then
      if ((m=='N').or.(m=='A')) then
         matrix_com_style = m
      else
         write(6,*)'Error: Illegal matrix_com_style =', m
      endif 
   endif
   end subroutine set_com_style
   !-------------------------------------------------------


   !-------------------------------------------------------------------
   !****p* io_mod/set_com_use
   ! SUBROUTINE
   !   set_com_use(u) - choose whether to echo input comments
   ! PURPOSE
   !    Reset default comment usage variable default_com_use 
   !    to input value of argument u
   !-------------------------------------------------------------------
   subroutine set_com_use(u)
   character(1), intent(IN) :: u  ! new value of default_com_use
   !***
   if ((u=='K').or.(u=='R')) then
      default_com_use = u
   else if (u=='C') then
      write(6,*)'Error checking not yet implemented, ignored' 
   else
      write(6,*)'Error: Illegal scalar_com_use =', u
   endif 
   end subroutine set_com_use
   !-------------------------------------------------------



   !-------------------------------------------------------------------
   !****p* io_mod/set_output_fmt
   ! SUBROUTINE
   !   set_output_fmt - set output format (field widths, float format)
   !
   ! PURPOSE
   ! Routine allows new values to be set for the integer
   ! variables com_width , data_width, and/or frac_width,
   ! which determine field widths for comment and data 
   ! fields, and for the character(2) fmt_ef, which is the
   ! format specifier 'E' or 'F' used to output real numbers.
   !
   ! Also resets the format strings fmt_c, fmt_i, and/or 
   ! fmt_r, as needed, so as to agree with the new values
   ! of the field widths and fmt_ef.
   !
   ! Subroutine arguments (all are optional arguments):
   !
   !  c = com_width  = width of comment field
   !  d = data_width = width of scalar data field
   !  f = frac_width = # of digits after decimal in floats
   !  e = fmt_ef     = 'E', 'F', or 'ES' format style for floats
   !
   ! SOURCE
   !--------------------------------------------------------
   subroutine set_output_fmt(c,d,f,e)
   integer, intent(IN), optional      :: c  ! com_width 
   integer, intent(IN), optional      :: d  ! data_width
   integer, intent(IN), optional      :: f  ! frac_width
   character(2), intent(IN), optional :: e  ! fmt_ef
   !***

   if (present(c)) then
      if (c <= com_len) then
         com_width  = c
         fmt_c = 'A' // trim(int_string(com_width ))
      else
         write(6,*) 'Error: Input com_width > com_len'
      endif
   endif
   if (present(d)) then
      data_width = d
      fmt_i = 'I' // trim(int_string(data_width))
      fmt_l = 'L' // trim(int_string(data_width))
   endif
   if (present(f)) then
      frac_width = f
   endif
   if (present(e)) then
      if ((trim(e) =='E').or.(trim(e) == 'F').or.e =='SE' ) then
         fmt_ef = e
      else
         write(6,*) 'Error: Illegal value of fmt_ef=',e
      endif
   endif
   if (present(d).or.present(f).or.present(e)) then
      fmt_r = trim(fmt_ef) // trim(int_string(data_width))
      fmt_r = trim(fmt_r) // '.' // trim(int_string(frac_width))
   endif
   end subroutine set_output_fmt
   !-------------------------------------------------------


   ! Input subroutines (generic interface input)
  

   subroutine input_int(data,c,i,o,e,f,u)
   !--------------------------------------------------------------
   integer, intent(OUT)                   :: data
   character(*), optional                 :: c     ! comment
   integer, intent(IN), optional          :: i     ! input unit #
   integer, intent(IN), optional          :: o     ! output unit #
   integer, intent(IN), optional          :: e     ! echo flag
   character(1), intent(IN), optional     :: f     ! comment style
   character(1), intent(IN), optional     :: u     ! comment usage

   ! Local variables
   integer            :: iunit, ounit, echo
   character(1)       :: com_style, com_use
   character(com_len) :: comment

   ! Set defaults 
   if (present(i)) then
      iunit = i
   else
      iunit = default_iunit
   endif
   if (present(o)) then
      ounit = o
   else
      ounit = default_ounit   
   endif
   if (present(e)) then
      if ((e==0).or.(e==1)) then
         echo = e
      else
         write(6,*) 'Error: Illegal value of e=',e
      endif
   else
      echo = default_echo
   endif
   if (present(f)) then
      if ((f=='N').or.(f=='A').or.(f=='L').or.(f=='R')) then
         com_style = f
      else
         write(ounit,*) 'Error: incorrect f = ' , f , ' in in_int'
         com_style = scalar_com_style
      endif
   else
      com_style = scalar_com_style
   endif
   if (present(u)) then
      com_use = u
   else
      com_use = default_com_use
   endif
   if (.not.present(c)) com_use = 'K'

   ! Read data
   select case(com_style)
   case('N') ! No comment
      read(iunit,*) data
   case('A') ! comment above
      read(iunit,*) comment
      read(iunit,*) data
   case('L') ! comment to left of data
      read(iunit,*) comment, data
   case('R') ! comment to right of data
      read(iunit,*) data, comment
   end select

   if (com_style /= 'N') then
      comment = adjustl(comment)
      if ((com_use == 'K').and.present(c)) then
         c = comment
      endif
   endif

   if (echo == 1) then ! echo in to ounit
      if ((com_use == 'R').and.present(c)) then
         comment = adjustl(c)
      endif
      call output(data,comment,ounit,com_style)
   endif

   end subroutine input_int
   !--------------------------------------------------------------


   subroutine input_real(data,c,i,o,e,f,u)
   !--------------------------------------------------------------
   real(long), intent(OUT)               :: data  ! 
   character(*), optional                :: c     ! comment
   integer, intent(IN), optional         :: i     ! input unit #
   integer, intent(IN), optional         :: o     ! output unit #
   integer, intent(IN), optional         :: e     ! echo flag
   character(1), intent(IN), optional    :: f     ! comment style
   character(1), intent(IN), optional    :: u     ! comment usage

   ! Local variables
   integer            :: iunit, ounit, echo
   character(1)       :: com_style, com_use
   character(com_len) :: comment

   ! Set defaults 
   if (present(i)) then
      iunit = i
   else
      iunit = default_iunit
   endif
   if (present(o)) then
      ounit = o    
   else
      ounit = default_ounit    
   endif
   if (present(e)) then
      if ((e==0).or.(e==1)) then
         echo = e
      else
         write(6,*) 'Error: Illegal value of e=',e
      endif
   else
      echo = default_echo
   endif
   if (present(f)) then
      if ((f=='N').or.(f=='A').or.(f=='L').or.(f=='R')) then
         com_style = f
      else
         write(ounit,*) 'Error: incorrect f = ' , f , ' in in_real'
         com_style = scalar_com_style
      endif
   else
      com_style = scalar_com_style
   endif
   if (present(u)) then
      com_use = u
   else
      com_use = default_com_use
   endif

   ! Read data
   select case(com_style)
   case('N') ! No comment
      read(iunit,*) data
   case('A') ! comment above
      read(iunit,*) comment
      read(iunit,*) data
   case('L') ! comment to left of data
      read(iunit,*) comment, data
   case('R') ! comment to right of data
      read(iunit,*) data, comment
   end select

   if (com_style /= 'N') then
      comment = adjustl(comment)
      if ((com_use == 'K').and.present(c)) then
         c = comment
      endif
   endif

   if (echo == 1) then 
      if ((com_use == 'R').and.present(c)) then
         comment = adjustl(c)
      endif
      call output(data,comment,ounit,com_style)
   endif

   end subroutine input_real
   !--------------------------------------------------------------


   subroutine input_char(data,c,i,o,e,f,u)
   !--------------------------------------------------------------
   character(*), intent(OUT)             :: data  ! input string 
   character(*), optional                :: c     ! comment
   integer, intent(IN), optional         :: i     ! input unit #
   integer, intent(IN), optional         :: o     ! output unit #
   integer, intent(IN), optional         :: e     ! echo flag
   character(1), intent(IN), optional    :: f     ! comment style
   character(1), intent(IN), optional    :: u     ! comment usage

   ! Local variables
   integer            :: iunit, ounit, echo
   character(1)       :: com_style, com_use
   character(com_len) :: comment

   ! Set defaults 
   if (present(i)) then
      iunit = i
   else
      iunit = default_iunit
   endif
   if (present(o)) then
      ounit = o    
   else
      ounit = default_ounit    
   endif
   if (present(e)) then
      if ((e==0).or.(e==1)) then
         echo = e
      else
         write(6,*) 'Error: Illegal value of e=',e
      endif
   else
      echo = default_echo
   endif
   if (present(f)) then
      if ((f=='N').or.(f=='A').or.(f=='L').or.(f=='R')) then
         com_style = f
      else
         write(ounit,*) 'Error: incorrect f = ' , f , ' in in_char'
         com_style = scalar_com_style
      endif
   else
      com_style = scalar_com_style
   endif
   if (present(u)) then
      com_use = u
   else
      com_use = default_com_use
   endif
   if (.not.present(c)) com_use = 'K'

   ! Read data
   select case(com_style)
   case('N') ! No comment
      read(iunit,*) data
   case('A') ! comment above
      read(iunit,*) comment
      read(iunit,*) data
   case('L') ! comment to left of data
      read(iunit,*) comment, data
   case('R') ! comment to right of data
      read(iunit,*) data, comment
   end select

   if (com_style /= 'N') then
      comment = adjustl(comment)
      if ((com_use == 'K').and.present(c)) then
         c = comment
      endif
   endif

   if (echo == 1) then 
      if ((com_use == 'R').and.present(c)) then
         comment = adjustl(c)
      endif
      call output(trim(data),comment,ounit,com_style)
   endif

   end subroutine input_char
   !--------------------------------------------------------------


   subroutine input_logic(data,c,i,o,e,f,u)
   !--------------------------------------------------------------
   logical, intent(OUT)                   :: data
   character(*), optional                 :: c     ! comment
   integer, intent(IN), optional          :: i     ! input unit #
   integer, intent(IN), optional          :: o     ! output unit #
   integer, intent(IN), optional          :: e     ! echo flag
   character(1), intent(IN), optional     :: f     ! comment style
   character(1), intent(IN), optional     :: u     ! comment usage

   ! Local variables
   integer            :: iunit, ounit, echo
   character(1)       :: com_style, com_use
   character(com_len) :: comment

   ! Set defaults 
   if (present(i)) then
      iunit = i
   else
      iunit = default_iunit
   endif
   if (present(o)) then
      ounit = o
   else
      ounit = default_ounit   
   endif
   if (present(e)) then
      if ((e==0).or.(e==1)) then
         echo = e
      else
         write(6,*) 'Error: Illegal value of e=',e
      endif
   else
      echo = default_echo
   endif
   if (present(f)) then
      if ((f=='N').or.(f=='A').or.(f=='L').or.(f=='R')) then
         com_style = f
      else
         write(ounit,*) 'Error: incorrect f = ' , f , ' in in_int'
         com_style = scalar_com_style
      endif
   else
      com_style = scalar_com_style
   endif
   if (present(u)) then
      com_use = u
   else
      com_use = default_com_use
   endif
   if (.not.present(c)) com_use = 'K'

   ! Read data
   select case(com_style)
   case('N') ! No comment
      read(iunit,*) data
   case('A') ! comment above
      read(iunit,*) comment
      read(iunit,*) data
   case('L') ! comment to left of data
      read(iunit,*) comment, data
   case('R') ! comment to right of data
      read(iunit,*) data, comment
   end select

   if (com_style /= 'N') then
      comment = adjustl(comment)
      if ((com_use == 'K').and.present(c)) then
         c = comment
      endif
   endif

   if (echo == 1) then ! echo in to ounit
      if ((com_use == 'R').and.present(c)) then
         comment = adjustl(c)
      endif
      call output(data,comment,ounit,com_style)
   endif

   end subroutine input_logic
   !--------------------------------------------------------------


   subroutine input_int_vec(data,n,c,s,i,o,e,f,u)
   !--------------------------------------------------------------
   integer, intent(OUT)                  :: data(:) ! in array 
   integer, intent(IN)                   :: n       ! dimension
   character(*), optional                :: c       ! comment
   character(1), intent(IN), optional    :: s       ! column or row
   integer, intent(IN), optional         :: i       ! input unit #
   integer, intent(IN), optional         :: o       ! output unit #
   integer, intent(IN), optional         :: e       ! echo flag
   character(1), intent(IN), optional    :: f       ! comment style
   character(1), intent(IN), optional    :: u       ! comment usage

   ! Local variables
   integer            :: iunit, ounit, echo, j
   character(1)       :: com_style, com_use, symmetry
   character(com_len) :: comment

   ! Set defaults 
   if (present(s)) then
      if ((s=='R').or.(s=='C')) then
         symmetry = s
      else
         write(6,*) 'Error: Illegal s=',s,' in input_int_vec'
         symmetry = default_vec_sym
      endif
   else
      symmetry = default_vec_sym
   endif
   if (present(i)) then
      iunit = i
   else
      iunit = default_iunit
   endif
   if (present(o)) then
      ounit = o    
   else
      ounit = default_ounit    
   endif
   if (present(e)) then
      if ((e==0).or.(e==1)) then
         echo = e
      else
         write(6,*) 'Error: Illegal value of e=',e
      endif
   else
      echo = default_echo
   endif
   if (present(f)) then
      if ((f=='N').or.(f=='A').or.(f=='L').or.(f=='R')) then
         com_style = f
      else
         write(ounit,*) 'Error: incorrect f = ' , f , ' in in_int_vec'
         com_style = vector_com_style
      endif
   else
      com_style = vector_com_style
   endif
   if (present(u)) then
      com_use = u
   else
      com_use = default_com_use
   endif
   if (.not.present(c)) com_use = 'K'

   ! Comment line
   if (com_style=='A') then
      read(iunit,*) comment
      comment = adjustl(comment)
      if ((com_use == 'K').and.present(c)) then
         c = comment
      endif
   endif

   ! Read data
   select case(symmetry)
   case('R')
      read(iunit,*) data(1:n)
   case('C')
      do j=1, n
         read(iunit,*) data(j)
      enddo
   end select

   ! Echo
   if (echo == 1) then 
      if ((com_use == 'R').and.present(c)) then
         comment = adjustl(c)
      endif
      call output(data,n,comment,symmetry,ounit,com_style)
   endif

   end subroutine input_int_vec
   !--------------------------------------------------------------


   subroutine input_real_vec(data,n,c,s,i,o,e,f,u)
   !--------------------------------------------------------------
   real(long), intent(OUT)               :: data(:) ! input array 
   integer, intent(IN)                   :: n       ! dimension
   character(*), optional                :: c       ! output comment
   character(1), intent(IN), optional    :: s       ! column or row
   integer, intent(IN), optional         :: i       ! input unit #
   integer, intent(IN), optional         :: o       ! output unit #
   integer, intent(IN), optional         :: e       ! echo flag
   character(1), intent(IN), optional    :: f       ! comment style
   character(1), intent(IN), optional    :: u       ! comment usage

   ! Local variables
   integer            :: iunit, ounit, echo, j
   character(1)       :: com_style, com_use, symmetry
   character(com_len) :: comment

   ! Set defaults 
   if (present(s)) then
      if ((s=='R').or.(s=='C')) then
         symmetry = s
      else
         write(6,*) 'Error: Illegal s=',s,' in input_real_vec'
         symmetry = default_vec_sym
      endif
   else
      symmetry = default_vec_sym
   endif
   if (present(i)) then
      iunit = i
   else
      iunit = default_iunit
   endif
   if (present(o)) then
      ounit = o    
   else
      ounit = default_ounit    
   endif
   if (present(e)) then
      if ((e==0).or.(e==1)) then
         echo = e
      else
         write(6,*) 'Error: Illegal value of e=',e
      endif
   else
      echo = default_echo
   endif
   if (present(f)) then
      if ((f=='N').or.(f=='A').or.(f=='L').or.(f=='R')) then
         com_style = f
      else
         write(ounit,*) 'Error: incorrect f = ' , f , ' in in_real_vec'
         com_style = vector_com_style
      endif
   else
      com_style = vector_com_style
   endif
   if (present(u)) then
      com_use = u
   else
      com_use = default_com_use
   endif
   if (.not.present(c)) com_use = 'K'

   ! Comment line
   if (com_style=='A') then
      read(iunit,*) comment
      comment = adjustl(comment)
      if ((com_use == 'K').and.present(c)) then
         c = comment
      endif
   endif

   ! Read data
   select case(symmetry)
   case('R')
      read(iunit,*) data(1:n)
   case('C')
      do j=1, n
         read(iunit,*) data(j)
      enddo
   end select

   ! Echo
   if (echo == 1) then ! echo in to ounit
      if ((com_use == 'R').and.present(c)) then
         comment = adjustl(c)
      endif
      call output(data,n,comment,symmetry,ounit,com_style)
   endif

   end subroutine input_real_vec
   !--------------------------------------------------------------


   subroutine input_char_vec(data,n,c,s,i,o,e,f,u)
   !--------------------------------------------------------------
   character(*), intent(OUT)             :: data(:) ! in array 
   integer, intent(IN)                   :: n       ! dimension
   character(*), optional                :: c       ! comment
   character(1), intent(IN), optional    :: s       ! column or row
   integer, intent(IN), optional         :: i       ! input unit #
   integer, intent(IN), optional         :: o       ! output unit #
   integer, intent(IN), optional         :: e       ! echo flag
   character(1), intent(IN), optional    :: f       ! comment style
   character(1), intent(IN), optional    :: u       ! comment usage

   ! Local variables
   integer            :: iunit, ounit, echo, j
   character(1)       :: com_style, com_use, symmetry
   character(com_len) :: comment

   ! Set defaults 
   if (present(s)) then
      if ((s=='R').or.(s=='C')) then
         symmetry = s
      else
         write(6,*) 'Error: Illegal s=',s,' in input_int_vec'
         symmetry = default_vec_sym
      endif
   else
      symmetry = default_vec_sym
   endif
   if (present(i)) then
      iunit = i
   else
      iunit = default_iunit
   endif
   if (present(o)) then
      ounit = o    
   else
      ounit = default_ounit    
   endif
   if (present(e)) then
      if ((e==0).or.(e==1)) then
         echo = e
      else
         write(6,*) 'Error: Illegal value of e=',e
      endif
   else
      echo = default_echo
   endif
   if (present(f)) then
      if ((f=='N').or.(f=='A').or.(f=='L').or.(f=='R')) then
         com_style = f
      else
         write(ounit,*) 'Error: incorrect f = ' , f , ' in in_int_vec'
         com_style = vector_com_style
      endif
   else
      com_style = vector_com_style
   endif
   if (present(u)) then
      com_use = u
   else
      com_use = default_com_use
   endif
   if (.not.present(c)) com_use = 'K'

   ! Comment line
   if (com_style=='A') then
      read(iunit,*) comment
      comment = adjustl(comment)
      if ((com_use == 'K').and.present(c)) then
         c = comment
      endif
   endif

   ! Read data
   select case(symmetry)
   case('R')
      read(iunit,*) data(1:n)
   case('C')
      do j=1, n
         read(iunit,*) data(j)
      enddo
   end select

   ! Echo
   if (echo == 1) then 
      if ((com_use == 'R').and.present(c)) then
         comment = adjustl(c)
      endif
      call output(data,n,comment,symmetry,ounit,com_style)
   endif

   end subroutine input_char_vec


   subroutine input_int_mat(data,m,n,c,s,i,o,e,f,u)
   !--------------------------------------------------------------
   integer, intent(OUT)                  :: data(:,:) ! in array 
   integer, intent(IN)                   :: m, n      ! dimensions
   character(*), optional                :: c         ! comment
   character(1), intent(IN), optional    :: s         ! symmetry flag
   integer, intent(IN), optional         :: i         ! input unit #
   integer, intent(IN), optional         :: o         ! output unit #
   integer, intent(IN), optional         :: e         ! echo flag
   character(1), intent(IN), optional    :: f         ! comment style
   character(1), intent(IN), optional    :: u         ! comment usage

   ! Local variables
   integer            :: iunit, ounit, echo, k, l
   character(1)       :: com_style, com_use, symmetry
   character(com_len) :: comment

   ! Set defaults 
   if (present(s)) then
      if ((s=='N').or.(s=='S').or.(s=='L')) then
         symmetry = s
      else
         write(6,*) 'Error: Illegal s=',s,' in input_int_mat'
         symmetry = default_mat_sym
      endif
   else
      symmetry = default_mat_sym
   endif
   if (present(i)) then
      iunit = i
   else
      iunit = default_iunit
   endif
   if (present(o)) then
      ounit = o   
   else
      ounit = default_ounit    
   endif
   if (present(e)) then
      if ((e==0).or.(e==1)) then
         echo = e
      else
         write(6,*) 'Error: Illegal value of e=',e
      endif
   else
      echo = default_echo
   endif
   if (present(f)) then
      if ((f=='N').or.(f=='A')) then
         com_style = f
      else
         write(ounit,*) 'Error: incorrect f = ' , f , ' in in_int_mat'
         com_style = matrix_com_style
      endif
   else
      com_style = matrix_com_style
   endif
   if (present(u)) then
      com_use = u
   else
      com_use = default_com_use
   endif
   if (.not.present(c)) com_use = 'K'

   ! Comment line
   if (com_style=='A') then
      read(iunit,*) comment
      comment = adjustl(comment)
      if ((com_use == 'K').and.present(c)) then
         c = comment
      endif
   endif

   ! Read data
   select case(symmetry)
   case ('N')
      do k=1, m
         read(iunit,*) data(k,1:n)
      enddo
   case ('S')
      if (m /= n) then
         write(6,*) &
         'Error: Symmetry flag S illegal for non-square matrix'
      endif
      do k=1, m
         read(iunit,*) data(k,1:k)
         do l=1,k
            data(l,k) = data(k,l)
         enddo
      enddo
   case ('L')
      if (m /= n) then
         write(6,*) &
         'Error: Symmetry flag L illegal for non-square matrix'
      endif
      if (m > 1) then
         do k=2, m
            read(iunit,*) data(k,1:k-1)
            do l=1, k-1
               data(l,k) = data(k,l)
            enddo
         enddo
         do k=1, m
            data(k,k) = 0
         enddo
      endif
   case default
      write(ounit,*) 'Error: Illegal symmetry=',S,' in in_int_mat'
      stop
   end select

   if (echo == 1) then 
      if ((com_use == 'R').and.present(c)) then
         comment = adjustl(c)
      endif
      if (present(s)) then
         call output(data,m,n,comment,s,ounit,com_style)
      else
         call output(data,m,n,comment,o=ounit,f=com_style)
      endif
   endif

   end subroutine input_int_mat
   !--------------------------------------------------------------


   subroutine input_real_mat(data,m,n,c,s,i,o,e,f,u)
   !--------------------------------------------------------------
   real(long), intent(OUT)               :: data(:,:) ! input array 
   integer, intent(IN)                   :: m, n      ! dimensions
   character(*), optional                :: c         ! comment
   character(1), intent(IN), optional    :: s         ! symmetry flag
   integer, intent(IN), optional         :: i         ! input unit #
   integer, intent(IN), optional         :: o         ! output unit #
   integer, intent(IN), optional         :: e         ! echo flag
   character(1), intent(IN), optional    :: f         ! comment style
   character(1), intent(IN), optional    :: u         ! comment usage

   ! Local variables
   integer            :: iunit, ounit, echo, k, l
   character(1)       :: com_style, com_use, symmetry
   character(com_len) :: comment

   ! Set defaults 
   if (present(s)) then
      if ((s=='N').or.(s=='S').or.(s=='L')) then
         symmetry = s
      else
         write(6,*) 'Error: Illegal s=',s,' in input_real_mat'
         symmetry = default_mat_sym
      endif
   else
      symmetry = default_mat_sym
   endif
   if (present(i)) then
      iunit = i
   else
      iunit = default_iunit
   endif
   if (present(o)) then
      ounit = o    
   else
      ounit = default_ounit    
   endif
   if (present(e)) then
      if ((e==0).or.(e==1)) then
         echo = e
      else
         write(6,*) 'Error: Illegal value of e=',e
      endif
   else
      echo = default_echo
   endif
   if (present(f)) then
      if ((f=='N').or.(f=='A')) then
         com_style = f
      else
         write(ounit,*) 'Error: incorrect f = ' , f , ' in in_real_mat'
         com_style = matrix_com_style
      endif
   else
      com_style = matrix_com_style
   endif
   if (present(u)) then
      com_use = u
   else
      com_use = default_com_use
   endif
   if (.not.present(c)) com_use = 'K'

   ! Comment line
   if (com_style=='A') then
      read(iunit,*) comment
      comment = adjustl(comment)
      if ((com_use == 'K').and.present(c)) then
         c = comment
      endif
   endif

   ! Read data
   select case(symmetry)
   case ('N')
      do k=1, m
         read(iunit,*) data(k,1:n)
      enddo
   case ('S')
      if (m /= n) then
         write(6,*) &
         'Error: Symmetry flag S illegal for non-square matrix'
      endif
      do k=1, m
         read(iunit,*) data(k,1:k)
         do l=1,k
            data(l,k) = data(k,l)
         enddo
      enddo
   case ('L')
      if (m /= n) then
         write(6,*) &
         'Error: Symmetry flag L illegal for non-square matrix'
      endif
      if (m > 1) then
         do k=2, m
            read(iunit,*) data(k,1:k-1)
            do l=1, k-1
               data(l,k) = data(k,l)
            enddo
         enddo
         do k=1, m
            data(k,k) = 0
         enddo
      endif
   case default
      write(ounit,*) 'Error: Illegal symmetry=',S,' in in_int_mat'
      stop
   end select

   ! Echo
   if (echo == 1) then 
      if ((com_use == 'R').and.present(c)) then
         comment = adjustl(c)
      endif
      if (present(s)) then
         call output(data,m,n,comment,s,ounit,com_style)
      else
         call output(data,m,n,comment,o=ounit,f=com_style)
      endif
   endif

   end subroutine input_real_mat
   !--------------------------------------------------------------


   ! Output subroutines 


   subroutine output_int(data,c,o,f)
   !--------------------------------------------------------------
   ! Output a single integer
   !--------------------------------------------------------------
   integer, intent(IN)                   :: data
   character(*) , intent(IN), optional   :: c
   integer, intent(IN), optional         :: o
   character(1), intent(IN), optional    :: f

   integer            :: ounit
   character(1)       :: com_style
   character(com_len) :: comment

   ! Set defaults
   if (present(o)) then
      ounit = o
   else
      ounit = default_ounit
   endif
   if (present(f)) then
      com_style = f      
   else
      com_style = scalar_com_style
   endif
   if (present(c)) then
      comment = adjustl(c)
   else
      comment = ' '
   endif

   select case(com_style)
   case('N')
      write(ounit,FMT=fmt_int(com_style)) data
   case('A','L')
      write(ounit,FMT=fmt_int(com_style)) comment, data
   case('R')
      write(ounit,FMT=fmt_int(com_style)) data, comment
   end select

   end subroutine output_int
   !--------------------------------------------------------------


   subroutine output_real(data,c,o,f)
   !--------------------------------------------------------------
   ! output a single real
   !--------------------------------------------------------------
   real(long), intent(IN)               :: data
   character(*), intent(IN), optional   :: c
   integer, intent(IN), optional        :: o
   character(1), intent(IN), optional   :: f

   integer            :: ounit
   character(1)       :: com_style
   character(com_len) :: comment

   ! Set defaults
   if (present(o)) then
      ounit = o
   else
      ounit = default_ounit
   endif
   if (present(f)) then
      com_style = f      
   else
      com_style = scalar_com_style
   endif
   if (present(c)) then
      comment = adjustl(c)
   else
      comment = ' '
   endif

   ! Output
   select case(com_style)
   case('N')
      write(ounit,FMT=fmt_real(com_style)) data
   case('A','L')
      write(ounit,FMT=fmt_real(com_style)) comment, data
   case('R')
      write(ounit,FMT=fmt_real(com_style)) data, comment
   end select
   end subroutine output_real
   !--------------------------------------------------------------


   subroutine output_char(data,c,o,f,l,j,q)
   !--------------------------------------------------------------
   character(*), intent(IN)             :: data
   character(*), intent(IN), optional   :: c  ! comment
   integer, intent(IN), optional        :: o  ! output file
   character(1), intent(IN), optional   :: f  ! format for comment
   integer, intent(IN), optional        :: l  ! length of string
   character(1), intent(IN), optional   :: j  ! justify = 'L' or 'R'
   logical, intent(IN), optional        :: q  ! quote = T or F

   integer            :: ounit, n_char, length
   character(1)       :: com_style, justify
   character(com_len) :: comment
   character(80)      :: buffer
   logical            :: quote

   ! Set defaults
   if (present(o)) then
      ounit = o
   else
      ounit = default_ounit
   endif
   if (present(f)) then
      com_style = f 
   else
      com_style = scalar_com_style
   endif
   if (present(c)) then
      comment = adjustl(c)
   else
      comment = ' '
   endif
   if (present(l)) then
      n_char = l
   else
      n_char = data_width
   endif
   if (present(j)) then
      justify = j
   else
      justify = 'R'
   endif
   if (present(q)) then
      quote = q
   else
      if (justify == 'R' .or. justify == 'r') then
          quote = .TRUE.
      else
          quote = .FALSE.
      endif
   endif

   ! Add quotes (or not)
   if (quote) then 
       buffer = adjustl("'" // trim(data)// "'")
   else
       buffer = adjustl(trim(data))
   endif

   ! Check length
   length = len(trim(buffer))
   if (length > n_char) then
      n_char = length
   endif

   ! Right justify data if required
   if (justify=='R'.or.justify=='r') then
      buffer(1:n_char) = adjustr(buffer(1:n_char))
   endif

   ! Output
   select case(com_style)
   case('N')
      write(ounit,FMT=fmt_char(com_style,l=n_char)) buffer(1:n_char)
   case('A','L')
      write(ounit,FMT=fmt_char(com_style,l=n_char)) comment, buffer(1:n_char)
   case('R')
      write(ounit,FMT=fmt_char(com_style,l=n_char)) buffer(1:n_char), comment
   end select

   end subroutine output_char
   !--------------------------------------------------------------


   subroutine output_logic(data,c,o,f)
   !--------------------------------------------------------------
   ! Output a single logical variable
   !--------------------------------------------------------------
   logical, intent(IN)                   :: data
   character(*) , intent(IN), optional   :: c
   integer, intent(IN), optional         :: o
   character(1), intent(IN), optional    :: f

   integer            :: ounit
   character(1)       :: com_style
   character(com_len) :: comment

   ! Set defaults
   if (present(o)) then
      ounit = o
   else
      ounit = default_ounit
   endif
   if (present(f)) then
      com_style = f      
   else
      com_style = scalar_com_style
   endif
   if (present(c)) then
      comment = adjustl(c)
   else
      comment = ' '
   endif

   select case(com_style)
   case('N')
      write(ounit,FMT=fmt_logic(com_style)) data
   case('A','L')
      write(ounit,FMT=fmt_logic(com_style)) comment, data
   case('R')
      write(ounit,FMT=fmt_logic(com_style)) data, comment
   end select

   end subroutine output_logic
   !--------------------------------------------------------------


   subroutine output_int_vec(data,n,c,s,o,f)
   !--------------------------------------------------------------
   ! output a 1D array (vector) of integers
   !--------------------------------------------------------------
   integer, intent(IN)                  :: data(:)
   integer, intent(IN)                  :: n 
   character(*), intent(IN), optional   :: c
   character(1), intent(IN), optional   :: s
   integer, intent(IN), optional        :: o
   character(1), intent(IN), optional   :: f

   integer            :: ounit, j
   character(1)       :: com_style, symmetry
   character(com_len) :: comment

   ! Set defaults
   if (present(s)) then
      if ((s=='R').or.(s=='C')) then
         symmetry = s
      else
         write(6,*) 'Error: Illegal s=',s,' in output_int_vec'
         symmetry = 'R'
      endif
   else
      symmetry = 'R'
   endif
   if (present(o)) then
      ounit = o
   else
      ounit = default_ounit
   endif
   if (present(f)) then
      com_style = f      
   else
      com_style = vector_com_style
   endif
   if (present(c)) then
      comment = adjustl(c)
   else
      comment = ' '
   endif

   ! Output
   select case(symmetry)
   case('R')
      select case(com_style)
      case('N')
         write(ounit,FMT=fmt_int('N',n)) data(1:n)
      case('A','L')
         write(ounit,FMT=fmt_int(com_style,n)) comment, data(1:n)
      case('R')
         write(ounit,FMT=fmt_int('R',n)) data(1:n), comment
      end select
   case('C')
      select case(com_style)
      case('N')
         do j=1, n
            write(ounit,FMT=fmt_int('N',1)) data(j)
         enddo
      case('A')
         write(ounit,FMT=fmt_int('A',1)) comment, data(1)
         if (n > 1) then
            do j=2, n
               write(ounit,FMT=fmt_int('N',1)) data(j)
            enddo
         endif
      case('L')
         write(ounit,FMT=fmt_int(com_style,1)) comment, data(1)
         if (n > 1) then
            do j=2, n
               write(ounit,FMT=fmt_int('L',1)) ' ', data(j)
            enddo
         endif
      case('R')
         write(ounit,FMT=fmt_int(com_style,1)) data(1), comment
         if (n > 1) then
            do j=2, n
               write(ounit,FMT=fmt_int('N',1)) data(j)
            enddo
         endif
      end select
   end select

   end subroutine output_int_vec
   !--------------------------------------------------------------


   subroutine output_real_vec(data,n,c,s,o,f)
   !--------------------------------------------------------------
   ! output a 1D array (vector) of reals
   !--------------------------------------------------------------
   real(long), intent(IN)               :: data(:)
   integer, intent(IN)                  :: n 
   character(*), intent(IN), optional   :: c
   character(1), intent(IN), optional   :: s
   integer, intent(IN), optional        :: o
   character(1), intent(IN), optional   :: f

   integer            :: ounit, j
   character(1)       :: com_style, symmetry
   character(com_len) :: comment

   ! Set defaults
   if (present(s)) then
      if ((s=='R').or.(s=='C')) then
         symmetry = s
      else
         write(6,*) 'Error: Illegal s=',s,' in output_int_vec'
         symmetry = 'R'
      endif
   else
      symmetry = 'R'
   endif
   if (present(o)) then
      ounit = o
   else
      ounit = default_ounit
   endif
   if (present(f)) then
      com_style = f      
   else
      com_style = vector_com_style
   endif
   if (present(c)) then
      comment = adjustl(c)
   else
      comment = ' '
   endif

   ! Output
   select case(symmetry)
   case('R')
      select case(com_style)
      case('N')
         write(ounit,FMT=fmt_real(com_style,n)) data(1:n)
      case('A','L')
         write(ounit,FMT=fmt_real(com_style,n)) comment, data(1:n)
      case('R')
         write(ounit,FMT=fmt_real(com_style,n)) data(1:n), comment
      end select
   case('C')
      select case(com_style)
      case('N')
         do j=1, n
            write(ounit,FMT=fmt_real('N',1)) data(j)
         enddo
      case('A')
         write(ounit,FMT=fmt_real(com_style,1)) comment, data(1)
         if (n > 1) then
            do j=2, n
               write(ounit,FMT=fmt_real('N',1)) data(j)
            enddo
         endif
      case('L')
         write(ounit,FMT=fmt_real(com_style,1)) comment, data(1)
         if (n > 1) then
            do j=2, n
               write(ounit,FMT=fmt_real('L',1)) ' ', data(j)
            enddo
         endif
      case('R')
         write(ounit,FMT=fmt_real(com_style,1)) data(1), comment
         if (n > 1) then
            do j=2, n
               write(ounit,FMT=fmt_real('N',1)) data(j)
            enddo
         endif
      end select
   end select

   end subroutine output_real_vec
   !--------------------------------------------------------------


   subroutine output_char_vec(data,n,c,s,o,f)
   !--------------------------------------------------------------
   ! output a 1D array (vector) of character(*) data
   !--------------------------------------------------------------
   character(*), intent(IN)             :: data(:)
   integer, intent(IN)                  :: n 
   character(*), intent(IN), optional   :: c
   character(1), intent(IN), optional   :: s
   integer, intent(IN), optional        :: o
   character(1), intent(IN), optional   :: f

   integer            :: ounit, j
   character(1)       :: com_style, symmetry
   character(com_len) :: comment
   character(80)      :: buffer(n), shift

   ! Set defaults
   if (present(s)) then
      if ((s=='R').or.(s=='C')) then
         symmetry = s
      else
         write(6,*) 'Error: Illegal s=',s,' in output_char_vec'
         symmetry = 'R'
      endif
   else
      symmetry = 'R'
   endif
   if (present(o)) then
      ounit = o
   else
      ounit = default_ounit
   endif
   if (present(f)) then
      com_style = f      
   else
      com_style = vector_com_style
   endif
   if (present(c)) then
      comment = adjustl(c)
   else
      comment = ' '
   endif

   ! Right justify data
   do j=1, n
      shift = adjustl( "'" // trim(adjustl(data(j))) // "'" )
      shift(1:data_width) = adjustr(shift(1:data_width))
      buffer(j) = trim(shift)
   enddo

   ! Output
   select case(symmetry)
   case('R')
      select case(com_style)
      case('N')
         write(ounit,FMT=fmt_char('N',n=n)) buffer(1:n)
      case('A','L')
         write(ounit,FMT=fmt_char(com_style,n=n)) comment, buffer(1:n)
      case('R')
         write(ounit,FMT=fmt_char('R',n=n)) buffer(1:n), comment
      end select
   case('C')
      select case(com_style)
      case('N')
         do j=1, n
            write(ounit,FMT=fmt_char('N',n=1)) buffer(j)
         enddo
      case('A')
         write(ounit,FMT=fmt_char('A',n=1)) comment, buffer(1)
         if (n > 1) then
            do j=2, n
               write(ounit,FMT=fmt_char('N',n=11)) buffer(j)
            enddo
         endif
      case('L')
         write(ounit,FMT=fmt_char(com_style,n=1)) comment, buffer(1)
         if (n > 1) then
            do j=2, n
               write(ounit,FMT=fmt_char('L',n=1)) ' ', buffer(j)
            enddo
         endif
      case('R')
         write(ounit,FMT=fmt_char(com_style,n=1)) buffer(1), comment
         if (n > 1) then
            do j=2, n
               write(ounit,FMT=fmt_char('N',n=1)) buffer(j)
            enddo
         endif
      end select
   end select

   end subroutine output_char_vec
   !--------------------------------------------------------------


   subroutine output_int_mat(data,m,n,c,s,o,f)
   !--------------------------------------------------------------
   ! output a matrix of integers
   !--------------------------------------------------------------
   integer, intent(IN)                  :: data(:,:)
   integer, intent(IN)                  :: m, n
   character(*), intent(IN), optional   :: c
   character(1), intent(IN), optional   :: s
   integer, intent(IN), optional        :: o
   character(1), intent(IN), optional   :: f

   integer            :: ounit, k
   character(1)       :: com_style
   character(com_len) :: comment

   ! Set defaults
   if (present(o)) then
      ounit = o
   else
      ounit = default_ounit
   endif
   if (present(f)) then
      com_style = f      
   else
      com_style = matrix_com_style
   endif
   if (present(c)) then
      comment = adjustl(c)
   else
      comment = ' '
   endif

   ! Output
   if (com_style=='A') write(ounit, FMT = '('//fmt_c//')' ) comment
   if ((s=='N').or.(.not.present(s))) then
      do k=1, m
         write(ounit,FMT=fmt_int('N',n)) data(k,1:n)
      enddo
   else if (s == 'S') then
      if (m /= n) then
         write(6,*) &
         'Error: Symmetry flag S illegal for non-square matrix'
      endif
      do k=1, m
         write(ounit,FMT=fmt_int('N',n)) data(k,1:k)
      enddo
   else if (s == 'L') then
      if (m /= n) then
         write(6,*) &
         'Error: Symmetry flag L illegal for non-square matrix'
      endif
      if (m > 1) then
         do k=2, m
            write(ounit,FMT=fmt_int('N',n)) data(k,1:k-1)
         enddo
      endif
   else
      write(ounit,*) 'Error: Illegal value s=',s,' in in_int_mat'
   endif

   end subroutine output_int_mat
   !--------------------------------------------------------------


   subroutine output_real_mat(data,m,n,c,s,o,f)
   !--------------------------------------------------------------
   ! output a matrix of reals
   !--------------------------------------------------------------
   real(long), intent(IN)               :: data(:,:)
   integer, intent(IN)                  :: m, n  
   character(*), intent(IN), optional   :: c
   character(1), intent(IN), optional   :: s
   integer, intent(IN), optional        :: o
   character(1), intent(IN), optional   :: f

   integer            :: ounit, k
   character(1)       :: com_style
   character(com_len) :: comment

   ! Set defaults
   if (present(o)) then
      ounit = o
   else
      ounit = default_ounit
   endif
   if (present(f)) then
      com_style = f      
   else
      com_style = matrix_com_style
   endif
   if (present(c)) then
      comment = adjustl(c)
   else
      comment = ' '
   endif

   ! Output
   if (com_style=='A') write(ounit, FMT = '('//fmt_c//')' ) comment
   if ((s=='N').or.(.not.present(s))) then
      do k=1, m
         write(ounit,FMT=fmt_real('N',n)) data(k,1:N)
      enddo
   else if (s == 'S') then
      if (m /= n) then
         write(6,*) &
         'Error: Symmetry flag S illegal for non-square matrix'
      endif
      do k=1, m
         write(ounit,FMT=fmt_real('N',n)) data(k,1:k)
      enddo
   else if (s == 'L') then
      if (m /= n) then
         write(6,*) &
         'Error: Symmetry flag L illegal for non-square matrix'
      endif
      if (m > 1) then
         do k=2, m
            write(ounit,FMT=fmt_real('N',n)) data(k,1:k-1)
         enddo
      endif
   else
      write(ounit,*) 'Error: Illegal value s=',s,' in in_real_mat'
   endif

   end subroutine output_real_mat
   !--------------------------------------------------------------


   !--------------------------------------------------------------
   ! Functions that return output format strings for 1D arrays of
   ! integers and reals, and for character strings (all private)
   !--------------------------------------------------------------


   character(fmt_len) function fmt_int(com_style,n)
   !--------------------------------------------------------------
   ! fmt_int = format string for array of n integers (default n=1)
   !--------------------------------------------------------------
   character(1), intent(IN)       :: com_style
   integer, intent(IN), optional  :: n
   character(int_len)             :: n_c

   if (present(n).and.(n>1)) then
      n_c = int_string(n)
   else
      n_c = ' '
   endif

   select case(com_style)
   case('N')
      fmt_int  = '(' // trim(n_c) // fmt_i // ')'
   case('A')
      fmt_int  = '(' // fmt_c // '/' // trim(n_c) // fmt_i // ')'
   case('L')
      fmt_int  = '(' // fmt_c // ',' // trim(n_c) // fmt_i // ')'
   case('R')
      fmt_int  = '(' // trim(n_c) // fmt_i // ',5X,' // fmt_c // ')'
   end select

   end function fmt_int
   !--------------------------------------------------------------
 

   character(fmt_len) function fmt_real(com_style,n)
   !--------------------------------------------------------------
   ! fmt_real = format string for array of n reals (default n=1)
   !--------------------------------------------------------------
   character(1), intent(IN)        :: com_style
   integer, intent(IN), optional  :: n
   character(int_len)             :: n_c

   if (present(n).and.(n>1)) then
      n_c = int_string(n)
   else
      n_c = ' '
   endif

   select case(com_style)
   case('N')
      fmt_real  = '(' // trim(n_c) // fmt_r // ')'
   case('A')
      fmt_real  = '(' // fmt_c // '/' // trim(n_c) // fmt_r // ')'
   case('L')
      fmt_real  = '(' // fmt_c // ',' // trim(n_c) // fmt_r // ')'
   case('R')
      fmt_real  = '(' // trim(n_c) // fmt_r // ',5X,' // fmt_c // ')'
   end select

   end function fmt_real
   !--------------------------------------------------------------


   character(fmt_len) function fmt_char(com_style,n,l)
   !--------------------------------------------------------------
   ! fmt_char = format string for character(n) string
   !-------------------------------------------------------------
   character(1), intent(IN)         :: com_style   
   integer,  intent(IN), optional  :: n    ! # strings/array
   integer,  intent(IN), optional  :: l    ! # characters/string
   character(int_len)              :: l_c, n_c

   ! l = # of characters in character variable
   if (present(l)) then
      l_c = int_string(l)
   else
      l_c = int_string(data_width)
   endif
   ! default strings of width data_width

   ! n = # of character(l) variables in array
   if (present(n).and.(n>1)) then
      n_c = int_string(n)
   else
      n_c = ' '
   endif
   ! default: scalar string variable, rather than array

   select case(com_style)
   case('N')
      fmt_char = '(' // trim(n_c) // 'A' // trim(l_c) // ')'
   case('A')
      fmt_char = '(' // fmt_c // '/' // trim(n_c) // 'A' // trim(l_c) // ')'
   case('L')
      fmt_char = '(' // fmt_c // ',' // trim(n_c) // 'A' // trim(l_c) // ')'
   case('R')
      fmt_char = '(' // trim(n_c) // 'A' // trim(l_c) // ',5X,' // fmt_c // ')'
   end select

   end function fmt_char
   !--------------------------------------------------------------


   character(fmt_len) function fmt_logic(com_style,n)
   !--------------------------------------------------------------
   ! fmt_logic = format string for array of n logical (default n=1)
   !--------------------------------------------------------------
   character(1), intent(IN)       :: com_style
   integer, intent(IN), optional  :: n
   character(int_len)             :: n_c

   if (present(n).and.(n>1)) then
      n_c = int_string(n)
   else
      n_c = ' '
   endif

   select case(com_style)
   case('N')
      fmt_logic = '(' // trim(n_c) // fmt_l // ')'
   case('A')
      fmt_logic = '(' // fmt_c // '/'// trim(n_c) // fmt_l // ')'
   case('L')
      fmt_logic = '(' // fmt_c // ',' // trim(n_c) // fmt_l // ')'
   case('R')
      fmt_logic = '(' // trim(n_c) // fmt_l // ',5X,' // fmt_c // ')'
   end select

   end function fmt_logic
   !--------------------------------------------------------------



end module io_mod
