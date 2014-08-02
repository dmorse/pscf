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
!****m scf/group_mod
! MODULE
!   group_mod
! PURPOSE
!   Define derived data types and basic operations 
!   for space group symmetry operations and groups
! AUTHOR
!    David Morse (2002)
! SOURCE
!-----------------------------------------------------------------------
module group_mod
   use const_mod,   only : dim, long
   use version_mod, only : version_type, input_version, output_version
   implicit none
   
   private

   ! Derived types
   public :: symmetry_type              ! space group symmetry
   public :: group_type                 ! space group

   ! Generic interfaces
   public :: operator(.dot.)            ! .dot. products for vectors,
                                        !    matrices and symmetry_type
   public :: inverse                    ! inversion of 2D and 3D real 
                                        !    matrices and symmetry_type
   public :: equal                      ! equality with tolerance for variety
                                        !    of data types, including symmetries

   ! Public procedures
   public :: read_group, output_group   ! io for groups
   public :: make_group                 ! complete and check space group
   !***

   !Private (make public for testing)
   !public :: max_order
   !public :: output_symmetry
   !public :: table_type, make_table
   !public :: valid_basis
   !public :: index_of_inverse
   !public :: make_identity, shift_translation
   !public :: is_sub_group

   ! Parameters (Private)
   integer, parameter      :: max_order = 192
   real(long), parameter   :: epsilon   = 1.0E-6_long
   character(9), parameter :: Cartesian = 'Cartesian'
   character(9), parameter :: Bravais   = 'Bravais  '
  
   !------------------------------------------------------------------
   !****t group_mod/symmetry_type
   ! TYPE
   !    symmetry_type
   ! VARIABLE
   !   character(9) basis   = 'Cartesian' or 'Bravais '
   !   real(long)   m(3,3)  = point group matrix
   !   real(long)   v(3)    = translation vector
   ! COMMENT
   ! Conventions for symmetry_type and related derived types:
   !
   !  a) The effect of a symmetry on a position vector is to take
   !    
   !     R -> m .dot. R + v
   !
   !     where m .dot. R represents contraction with first index of m
   !
   !  b) Point group matrix m operates on reciprocal G vectors by
   !     contraction with first index of m:  G -> G .dot. m
   ! 
   !  c) Symmetries can be expressed in either Cartesian or
   !     Bravais basis, as indicated by value of the character
   !     variable basis, which can have values equal to the string
   !     constants Cartesian = 'Cartesian' or Bravais='Bravais  '.
   !
   !  d) In the Bravais basis, position vectors are represented
   !     as coefficients in expansion in Bravais basis vectors,
   !     R_basis(:,i), i=1,..dim, and G vectors are represented 
   !     as (integer) coefficients in expansion in reciprocal 
   !     basis  vectors, G_basis(:,j), j=1,..,dim .
   !
   !  e) In the Bravais representation, elements of a point 
   !     group matrix m should be integers (though they are 
   !     stored as reals), and elements of the translation
   !     vector v should be low order fractions
   !       
   ! SOURCE
   !------------------------------------------------------------------
   type symmetry_type
      character(9) :: basis   ! must equal 'Cartesian' or 'Bravais '
      real(long)   :: m(3,3)  ! point group matrix
      real(long)   :: v(3)    ! translation vector
   end type symmetry_type
   !***
  

   !------------------------------------------------------------------
   !****t group_mod/group_type
   ! TYPE
   !    group_type
   ! VARIABLE
   !    order        = # of symmetry elements in group
   !    s(max_order) = array of symmetries
   ! COMMENT
   !    All symmetries in a group must have the same basis, i.e.,
   !    they must all be either in Bravais or Cartesian basis
   ! SOURCE
   !------------------------------------------------------------------
   type group_type
      integer             :: order          
      type(symmetry_type) :: s(max_order)  
   end type group_type
   !***
  

   !------------------------------------------------------------------
   !****it group_mod/table_type
   ! TYPE
   !    table_type - multipication table for symmetries in a group
   !
   ! VARIABLE
   !
   !    order      = order of symmetries
   !    s(i,j)     = product symmetries i .dot. symmetry j
   !    index(i,j) = index of s(i,j) in group
   !                 so that s(i,j) = group%s( index(i,j) )
   !         
   ! SOURCE
   !------------------------------------------------------------------
   type table_type
      integer             :: order
      integer             :: index(max_order,max_order) 
      type(symmetry_type) :: s(max_order,max_order)     
   end type table_type
   !***


   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! Generic Interfaces
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   !------------------------------------------------------------------
   !****p group_mod/dot
   ! OPERATOR
   !    .dot. - dot product generic interface
   !
   ! COMMENT
   !
   !  a) In names of the specific realizations the types of
   !     arguments are indicated with the shorthand:
   !
   !         integer          ivec(:)
   !         real(long)       vec(:)
   !         real(long)       mat(:,:)
   !         symmetry_type    sym
   !
   ! 
   !  b) When evaluating dot products, elements of the input
   !     arguments with indices > dim are ignored, and elements
   !     of any returned vector or matrix with indices > dim
   !     are padded with zeros. Because the return value of 
   !     the operator cannot be adjustable, the operator returns
   !     vectors and matrices (when appropriate) with dimension 3
   !
   ! SOURCE
   !------------------------------------------------------------------
   interface operator(.dot.)
      module procedure ivec_dot_ivec  ! integer 
      module procedure vec_dot_vec    ! real    
      module procedure ivec_dot_vec   ! real   
      module procedure vec_dot_ivec   ! real    
      module procedure mat_dot_vec    ! real(3) 
      module procedure ivec_dot_mat   ! real(3) 
      module procedure vec_dot_mat    ! real(3) 
      module procedure mat_dot_mat    ! real(3,3) 
      module procedure sym_dot_vec    ! real(3) = sym%m.dot.vec + sym%v
      module procedure vec_dot_sym    ! real(3) = vec.dot.sym%m
      module procedure ivec_dot_sym   ! integer(3) = ivec.dot.sym%m
      module procedure sym_dot_sym    ! symmetry_type
   end interface
   !***
   
   
   !------------------------------------------------------------------
   !****p group_mod/equal
   ! FUNCTION
   !    equal(a,b)
   ! RETURN
   !
   !    true if a nd b are equal to within a tolerance epsilon
   !
   !    The arguments a and b may be objects of type:
   !
   !        real(long)
   !        integer       dimension(3)
   !        real(long)    dimension(3)
   !        real(long)    dimension(3,3)
   !        symmetry_type
   !
   ! SOURCE
   !------------------------------------------------------------------
   interface equal
      module procedure real_equal
      module procedure ivector_equal
      module procedure r_vector_equal
      module procedure matrix_equal
      module procedure symmetry_equal
   end interface
   !***
  
   
   !------------------------------------------------------------------
   !****p group_mod/inverse
   ! FUNCTION
   !
   !    inverse(a) - generic interface for inversion
   !
   ! RETURN
   !
   !    inverse of argument a
   ! 
   !    The argument a may be:
   !
   !       real(long)    matrix_inverse(3,3) 
   !       symmetry_type symmetry_inverse    
   !
   ! SOURCE
   !------------------------------------------------------------------
   interface inverse
      module procedure matrix_inverse   ! (3,3) padded with zeros
      module procedure symmetry_inverse ! symmetry_type
   end interface
   !***


contains
   
   
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! Output Subroutines for Data Types & Input for Groups
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


   !------------------------------------------------------------------
   !****ip group_mod/output_scalar
   !  SUBROUTINE
   !      output_scalar(r,iunit)
   !  PURPOSE
   !      Write real scalar to iunit 
   !  SOURCE
   !------------------------------------------------------------------
   subroutine output_scalar(r,iunit)
   real(long), intent(IN)  :: r
   integer                 :: iunit
   !***
   write(iunit,FMT='(F12.3)') r
   write(iunit,*)
   end subroutine output_scalar
   !===================================================================


   !------------------------------------------------------------------
   !****ip group_mod/output_vector
   !  SUBROUTINE
   !      output_vector(v,iunit)
   !  PURPOSE
   !      Write real vector to iunit 
   !  SOURCE
   !------------------------------------------------------------------
   subroutine output_vector(v,iunit)
   real(long), intent(in)  :: v(3)
   integer, intent(in)     :: iunit
   !***
   integer                 :: i
   ! Check vector size
   if (size(v) < dim) then
      write(6,*) 'Error: Too small a vector in output_matrix'
   endif
   write(iunit,FMT='(3F12.3)') (v(i), i=1, dim)
   write(iunit,*)
   end subroutine output_vector
   !===================================================================


   !-------------------------------------------------------------------
   !****ip group_mod/output_matrix
   ! SUBROUTINE
   !    output_matrix(m,iunit)
   ! PURPOSE
   !    Write matrix m to file iunit 
   ! SOURCE
   !-------------------------------------------------------------------
   subroutine output_matrix(m,iunit)
   real(long), intent(in)  :: m(:,:)
   integer           :: iunit
   !***
   integer                   :: i, j 

   ! Check matrix size
   if ( (size(m,1) < dim).or.(size(m,2)< dim)) then
      write(6,*) 'Error: Too small a matrix in output_matrix'
   endif
   do i=1, dim
      write(iunit,FMT='(3F12.3)') (m(i,j), j=1, dim)
   enddo
   write(iunit,*)

   end subroutine output_matrix
   !===================================================================


   !-------------------------------------------------------------------
   !****p group_mod/output_symmetry
   ! SUBROUTINE
   !    output_symmetry(s,iunit)
   ! PURPOSE
   !    Write symmetry s to file iunit 
   ! SOURCE
   !-------------------------------------------------------------------
   subroutine output_symmetry(s,iunit)
   type(symmetry_type), intent(in):: s
   integer           :: iunit
   !***
   integer                   :: i, j

   do i=1, dim
     if (dim == 3) then
        write(iunit,FMT='(3F8.3,5X,F8.3)') &
             (s%m(i,j),j=1, dim),s%v(i)
     else if (dim == 2) then
        write(iunit,FMT='(2F8.3,5X,F8.3)') &
             (s%m(i,j),j=1, dim),s%v(i)
     else if (dim == 1) then
        write(iunit,FMT='(F8.3,5X,F8.3)') &
             (s%m(i,j),j=1, dim),s%v(i)
     else
        write(iunit,*) 'invalid dim in output_symmetry'
     endif
   enddo
   write(iunit,*)
   end subroutine output_symmetry
   !===================================================================


   !-------------------------------------------------------------------
   !****p group_mod/output_group
   ! SUBROUTINE
   !    output_group(g,iunit)
   ! PURPOSE
   !    Write group g to file iunit 
   ! SOURCE
   !-------------------------------------------------------------------
   subroutine output_group(g,iunit)
   type(group_type), intent(IN) :: g
   integer, intent(IN)          :: iunit
   !***

   type(version_type) :: version
   integer :: i

   version%major = 1
   version%minor = 0
   call output_version(version, iunit)
   write(iunit,FMT="(A9,' = basis')") g%s(1)%basis
   write(iunit,FMT="(i9,' = dimension')") dim
   write(iunit,FMT="(i9,' = # of symmetries')") g%order
   write(iunit,*)
   do i=1, g%order
      write(iunit,FMT='(i3)') i
      call output_symmetry(g%s(i),iunit)
   enddo

   end subroutine output_group
   !===================================================================


   !-------------------------------------------------------------------
   !****ip group_mod/output_table
   ! SUBROUTINE
   !    output_table(g,iunit)
   ! PURPOSE
   !    Write table t to file iunit 
   ! SOURCE
   !-------------------------------------------------------------------
   subroutine output_table(t,iunit)
   type(table_type), intent(IN) :: t
   integer, intent(IN)          :: iunit
   !***
   integer :: i,j

   do i=1, t%order
      write(iunit,FMT='(48I3)') (t%index(i,j),j=1,t%order)
   enddo

   end subroutine output_table
   !===================================================================

 
   !-------------------------------------------------------------------
   !****p group_mod/read_group
   ! SUBROUTINE
   !    read_group(g,iunit)
   ! PURPOSE
   !    Read group g from file unit iunit 
   ! SOURCE
   !-------------------------------------------------------------------
   subroutine read_group(g,iunit)

   type(group_type), intent(OUT) :: g
   integer, intent(IN)           :: iunit
   !***

   integer            :: i, j, k, i_input, dim_input
   character(9)       :: basis
   type(version_type) :: version

   ! Read file format
   call input_version(version,iunit)

   ! Read character string basis and check validity
   read(iunit,*) basis
   if ( ( basis /= Cartesian ).and.( basis /= Bravais  ) ) then
      write(6,*) 'Invalid value of s%basis in read_group'
      stop
   endif

   ! Read dimension, and check consistency
   read(iunit,*) dim_input
   if ( dim_input /= dim ) then
      write(6,*) 'Incompatible values of dim in read_group'
      write(6,*) 'dim_input = ', dim_input
      write(6,*) 'dim       = ', dim
      stop
   endif

   do i=1, max_order
      g%s(i)%m = 0.0_long
      g%s(i)%v = 0.0_long
   enddo
   read(iunit,*) g%order
   do i=1, g%order
      read(iunit,*)
      read(iunit,*) i_input
      do j=1, dim
         read(iunit,*) (g%s(i)%m(j,k),k=1, dim),g%s(i)%v(j)
      enddo
      g%s(i)%basis = basis
   enddo

   end subroutine read_group
   !===================================================================


   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! Definitions of Dot Product
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


   !------------------------------------------------------------------
   !****ip group_mod/ivec_dot_ivec
   ! FUNCTION
   !   integer ivec_dot_ivec(v1,v2)
   ! RETURN
   !   Dot product of 2 integer vectors
   ! SOURCE
   !------------------------------------------------------------------
   integer function ivec_dot_ivec(v1,v2)
   integer, intent(IN) :: v1(:), v2(:)
   !***
   integer     :: i
   ivec_dot_ivec = 0
   do i=1, dim
      ivec_dot_ivec = ivec_dot_ivec + v1(i)*v2(i)
   enddo
   end function ivec_dot_ivec
   !===================================================================


   !------------------------------------------------------------------
   !****ip group_mod/vec_dot_vec
   !  FUNCTION
   !      real(long) vec_dot_vec(v1,v2)
   !  RETURN
   !      Dot product of 2 real vectors
   ! SOURCE
   !------------------------------------------------------------------
   real(long) function vec_dot_vec(v1,v2)
   real(long),intent(IN), dimension(:) :: v1, v2
   !***
   integer     ::i
   vec_dot_vec = 0.0_long
   do i=1, dim
      vec_dot_vec = vec_dot_vec + v1(i)*v2(i)
   enddo
   end function vec_dot_vec
   !===================================================================


   !------------------------------------------------------------------
   !****ip group_mod/ivec_dot_vec
   ! FUNCTION
   !   real(long) ivec_dot_vec(v1,v2)
   ! RETURN
   !   Dot product of integer vector .dot. real vector
   ! SOURCE
   !------------------------------------------------------------------
   real(long) function ivec_dot_vec(v1,v2)
   integer, intent(IN), dimension(:)    :: v1
   real(long), intent(IN), dimension(:) :: v2
   !***
   integer     ::i
   ivec_dot_vec = 0.0_long
   do i=1, dim
      ivec_dot_vec = ivec_dot_vec + v1(i)*v2(i)
   enddo
   end function ivec_dot_vec
   !===================================================================


   !------------------------------------------------------------------
   !****ip group_mod/vec_dot_ivec
   !  FUNCTION
   !      real(long) vec_dot_ivec(v1,v2)
   !  RETURN
   !      Dot product of real vector .dot. integer vector
   ! SOURCE
   !------------------------------------------------------------------
   real(long) function vec_dot_ivec(v1,v2)
   real(long), intent(IN), dimension(:) :: v1
   integer, intent(IN), dimension(:)    :: v2
   !***
   vec_dot_ivec = ivec_dot_vec(v2,v1)
   end function vec_dot_ivec
   !===================================================================


   !------------------------------------------------------------------
   !****ip group_mod/mat_dot_vec
   ! FUNCTION
   !   real(long) mat_dot_vec(m,v)
   ! RETURN
   !   Dot product of real matrix m .dot. real column vector v
   !   Returns real array mat_dot_vec(3)
   ! SOURCE
   !------------------------------------------------------------------
   function mat_dot_vec(m,v)
   real(long)             :: mat_dot_vec(3)
   real(long),intent(IN)  :: v(:)
   real(long),intent(IN)  :: m(:,:)
   !***
   integer     ::i,j
   mat_dot_vec = 0.0_long
   do i=1, dim
      do j=1, dim
         mat_dot_vec(i) = mat_dot_vec(i) + m(i,j)*v(j)
      enddo
   enddo
   end function mat_dot_vec
   !===================================================================


   !------------------------------------------------------------------
   !****ip group_mod/vec_dot_mat
   ! FUNCTION
   !   real(long) vec_dot_mat(v,m)
   ! RETURN
   !   Dot product of real row vector v .dot. real matrix m
   !   Returns real array vec_dot_mat(3)
   ! SOURCE
   !------------------------------------------------------------------
   function vec_dot_mat(v,m)
   real(long)             :: vec_dot_mat(3)
   real(long),intent(IN)  :: v(:)
   real(long),intent(IN)  :: m(:,:)
   !***
   integer     ::i,j
   vec_dot_mat = 0.0_long
   do j=1, dim
      do i=1, dim
         vec_dot_mat(j) = vec_dot_mat(j) + v(i)*m(i,j)
      enddo
   enddo
   end function vec_dot_mat
   !===================================================================


   !------------------------------------------------------------------
   !****ip group_mod/ivec_dot_mat
   ! FUNCTION
   !   real(long) ivec_dot_mat(v,m)
   ! RETURN
   !   Dot product of integer row vector v .dot. real matrix m
   !   Returns real array ivec_dot_mat(3)
   ! SOURCE
   !------------------------------------------------------------------
   function ivec_dot_mat(v,m)
   real(long)             :: ivec_dot_mat(3)
   integer, intent(IN)    :: v(:)
   real(long), intent(IN) :: m(:,:)
   !***
   integer     ::i,j
   ivec_dot_mat = 0.0_long
   do j=1, dim
      do i=1, dim
         ivec_dot_mat(j) = ivec_dot_mat(j) + v(i)*m(i,j)
      enddo
   enddo
   end function ivec_dot_mat
   !===================================================================


   !------------------------------------------------------------------
   !****ip group_mod/mat_dot_mat
   ! FUNCTION
   !   real(long) mat_dot_mat(m1,m2)
   ! RETURN
   !   Dot product of real matrix m1 .dot. real matrix m2
   !   Returns real array mat_dot_mat(3,3)
   ! SOURCE
   !------------------------------------------------------------------
   function mat_dot_mat(m1,m2)
   real(long)             :: mat_dot_mat(3,3)
   real(long),intent(IN)  :: m1(:,:), m2(:,:)
   !***
   integer     ::i,j,k
   mat_dot_mat = 0.0_long
   do i=1, dim
     do j=1, dim
       do k=1, dim
         mat_dot_mat(i,j) = mat_dot_mat(i,j) + m1(i,k)*m2(k,j)
       enddo
     enddo
   enddo
   end function mat_dot_mat
   !===================================================================


   !------------------------------------------------------------------
   !****ip group_mod/sym_dot_vec
   ! FUNCTION
   !   real(long) sym_dot_vec(s,v)
   ! RETURN
   !   Result of symmetry element s acting on position vector v
   !   Returns real vector sym_dot_vec(3) = s%m .dot. v + s%v
   ! COMMENT
   !   Result is valid iff basis of symmetry and vector are same
   !   (i.e., both Cartesian or both Bravais bases)
   ! SOURCE
   !-------------------------------------------------------------------
   function sym_dot_vec(s,v)
   real(long), dimension(3)     :: sym_dot_vec
   type(symmetry_type), intent(IN):: s
   real(long),intent(IN)          :: v(:)
   !***
   integer:: i
   do i=1, dim
      sym_dot_vec = mat_dot_vec(s%m,v) + s%v
   enddo
   end function sym_dot_vec
   !===================================================================


   !------------------------------------------------------------------
   !****ip group_mod/vec_dot_sym
   ! FUNCTION
   !   real(long) vec_dot_sym(v,s)
   ! RETURN
   !   Result of symmetry s upon real Cartesian wavevector v
   !   Returns real array vec_dot_sym(3) = v .dot. s%m
   ! COMMENT
   !   Valid only if symmetry and wavector are in Cartesian basis
   ! SOURCE
   !-------------------------------------------------------------------
   function vec_dot_sym(v,s)
   real(long), dimension(3)             :: vec_dot_sym
   real(long), intent(IN), dimension(:) :: v
   type(symmetry_type), intent(IN)      :: s
   !***
   if ( s%basis == Cartesian  ) then
      vec_dot_sym = vec_dot_mat(v,s%m)
   else
      write(6,*) 'Improper use of vec_dot_sym - not a Cartesian symmetry'
      stop
   endif
   end function vec_dot_sym
   !===================================================================


   !------------------------------------------------------------------
   !****ip group_mod/ivec_dot_sym
   ! FUNCTION
   !   real(long) ivec_dot_sym(v,s)
   ! RETURN
   !   Result of symmetry s upon integer Cartesian wavevector v
   !   Returns integer array vec_dot_sym(3) = v .dot. s%m
   ! COMMENT
   !   Valid only if symmetry and wavector are in Bravais basis
   ! SOURCE
   !-------------------------------------------------------------------
   function ivec_dot_sym(ivec,s)
   integer, dimension(3)             :: ivec_dot_sym
   integer, intent(IN), dimension(:) :: ivec
   type(symmetry_type), intent(IN)   :: s
   !***
   if (s%basis == 'Bravais  ') then
      ivec_dot_sym = nint( ivec_dot_mat(ivec,s%m) )
   else
      write(6,*) 'Improper use of ivec_dot_sym - not a Bravais symmetry'
      stop
   endif
   end function ivec_dot_sym
   !===================================================================


   !------------------------------------------------------------------
   !****ip group_mod/sym_dot_sym
   ! FUNCTION
   !   type(symmetry_type) sym_dot_sym(v,s)
   ! RETURN
   !   Dot product of two symmetries
   ! SOURCE
   !-------------------------------------------------------------------
   type(symmetry_type) function sym_dot_sym(s1,s2)
   type(symmetry_type), intent(IN) :: s1,s2
   !***
   type(symmetry_type)             :: s3 
   ! Variable s3 added to compile on Regatta (6/2004)
   if (s1%basis == s2%basis) then
      s3%m     = mat_dot_mat(s1%m,s2%m)
      s3%v     = mat_dot_vec(s1%m,s2%v) + s1%v
      s3%basis = s1%basis
   else
      write(6,*) 'Incompatible bases in sym_dot_sym'
      stop
   endif
   sym_dot_sym = s3
   end function sym_dot_sym
   !===================================================================


   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! Definitions of inverse
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


   !-------------------------------------------------------------------
   !****ip group_mod/matrix_inverse
   ! FUNCTION
   !   matrix_inverse(m)
   !   real(long) dimension(3,3) :: matrix_inverse(3,3)
   ! RETURN
   !   Matrix inverse for dim=1, dim=2 or dim=3
   ! SOURCE
   !-------------------------------------------------------------------
   function matrix_inverse(m)
   real(long) :: matrix_inverse(3,3)
   real(long), intent(IN):: m(3,3)
   !***
   real(long) :: cofac(3,3)
   real(long):: det
   integer :: i, j
   matrix_inverse = 0.0_long
   if (dim == 3) then
      cofac(1,1) = m(2,2)*m(3,3) - m(3,2)*m(2,3)
      cofac(1,2) = m(2,1)*m(3,3) - m(3,1)*m(2,3)
      cofac(1,3) = m(2,1)*m(3,2) - m(3,1)*m(2,2)
      cofac(2,1) = m(1,2)*m(3,3) - m(3,2)*m(1,3)
      cofac(2,2) = m(1,1)*m(3,3) - m(3,1)*m(1,3)
      cofac(2,3) = m(1,1)*m(3,2) - m(3,1)*m(1,2)
      cofac(3,1) = m(1,2)*m(2,3) - m(2,2)*m(1,3)
      cofac(3,2) = m(1,1)*m(2,3) - m(2,1)*m(1,3)
      cofac(3,3) = m(1,1)*m(2,2) - m(2,1)*m(1,2)
      cofac(1,2) = -cofac(1,2)
      cofac(2,3) = -cofac(2,3)
      cofac(2,1) = -cofac(2,1)
      cofac(3,2) = -cofac(3,2)
      det = 0.0_long
      do i=1, 3
         det = det + m(1,i)*cofac(1,i) 
      enddo
      do i=1, 3
         do j=1, 3
            matrix_inverse(i,j) = cofac(j,i)/det
         enddo
      enddo
   else if (dim == 2) then
      cofac(1,1) =  m(2,2)
      cofac(2,2) =  m(1,1)
      cofac(1,2) = -m(2,1)
      cofac(2,1) = -m(1,2)
      det = m(1,1)*m(2,2) - m(1,2)*m(2,1)
      do i=1, 2
         do j=1, 2
            matrix_inverse(i,j) = cofac(j,i)/det
         enddo
      enddo
   else if (dim == 1) then
      matrix_inverse(1,1) = 1.0_long/m(1,1)
   else
      write(6,*) 'Invalid dim in matrix_inverse'
      stop
   endif
   end function matrix_inverse
   !===================================================================


   !-------------------------------------------------------------------
   !****ip group_mod/symmetry_inverse
   ! FUNCTION
   !   type(symmetry_type) symmetry_inverse(m)
   ! RETURN
   !   Inverse of symmetry operation s
   ! SOURCE
   !-------------------------------------------------------------------
   type(symmetry_type) function symmetry_inverse(s)
   type(symmetry_type), intent(IN) :: s
   !***
   real(long)                      :: inv_m(3,3)
   ! Inverse stored in variable inv_m to compile on Regatta (6/2004)

   symmetry_inverse%basis = s%basis
   inv_m = matrix_inverse(s%m) 
   symmetry_inverse%m = inv_m
   symmetry_inverse%v = (-1.0d0)*(symmetry_inverse%m .dot. s%v)
   end function symmetry_inverse
   !===================================================================


   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! Definitions of logical function equal(A,B)
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


   !-------------------------------------------------------------------
   !****ip group_mod/real_equal
   ! FUNCTION
   !   logical real_equal(r1,r2)
   ! RETURN
   !   Compares two real numbers r1 and r2
   !   True if |r1 - r2| <= epsilon
   ! SOURCE
   !-------------------------------------------------------------------
   logical function real_equal(r1,r2)
   real(long), intent(IN) :: r1, r2
   !***
   real_equal = .TRUE.
   if (abs(r1-r2) > epsilon ) then
      real_equal = .FALSE.
   endif
   end function real_equal
   !===================================================================


   !-------------------------------------------------------------------
   !****ip group_mod/ivector_equal
   ! FUNCTION
   !   logical ivector_equal(v1,v2)
   ! RETURN
   !   Compares two integer vectors of logical dimension dim 
   !   True if v1 = v2
   ! SOURCE
   !-------------------------------------------------------------------
   logical function ivector_equal(v1,v2)
   integer, intent(IN) :: v1(:), v2(:)
   !***
   integer             :: i
   ivector_equal = .TRUE.
   do i=1, dim 
      if ( v1(i) /= v2(i) ) then
         ivector_equal = .FALSE.
         return
      endif
   enddo
   end function ivector_equal
   !===================================================================


   !-------------------------------------------------------------------
   !****ip group_mod/r_vector_equal
   ! FUNCTION
   !   logical r_vector_equal(v1,v2)
   ! RETURN
   !   Compares two real(long) vectors of logical dimension dim 
   !   True if | v1(i) - v2(i) | < epsilon for i=1,...,dim  
   ! SOURCE
   !-------------------------------------------------------------------
   logical function r_vector_equal(v1,v2)
   real(long), intent(IN) :: v1(:), v2(:)
   !***
   integer :: i
   r_vector_equal = .TRUE.
   do i=1, dim 
      if (abs(v1(i)-v2(i)) > epsilon ) then
         r_vector_equal = .FALSE.
         return
      endif
   enddo
   end function r_vector_equal
   !===================================================================


   !-------------------------------------------------------------------
   !****ip group_mod/matrix_equal
   ! FUNCTION
   !   logical matrix_equal(m1,m2)
   ! RETURN
   !   Compares two real matrices of logical dimension dim x dim
   !   True if | m1(i,j) - m2(i,j) | < epsilon for i,j = 1,...,dim  
   ! SOURCE
   !-------------------------------------------------------------------
   logical function matrix_equal(m1,m2)
   real(long), intent(IN) :: m1(:,:), m2(:,:)
   !***
   integer :: i, j
   matrix_equal = .TRUE.
   do i=1, dim 
      do j=1, dim
         if ( abs(m1(i,j)-m2(i,j)) > epsilon ) then
            matrix_equal = .FALSE.
            return
         endif
      enddo
   enddo
   end function matrix_equal
   !===================================================================


   !-------------------------------------------------------------------
   !****ip group_mod/symmetry_equal
   ! FUNCTION
   !   logical symmetry_equal(s1,s2,R_basis,G_basis)
   ! RETURN
   !   Compares two symmetry elements. True if matrices s1%m and s2%m 
   !   are equal and translations s1%v and s2%v are equal modulo a 
   !   lattice translation
   ! SOURCE
   !-------------------------------------------------------------------
   logical function symmetry_equal(s1,s2,R_basis,G_basis)
   type(symmetry_type), intent(IN) :: s1, s2
   real(long), intent(IN) :: R_basis(:,:), G_basis(:,:)
   !***

   !  Local variables
   real(long) :: v1(3), v2(3)

   ! Check that both symmetries are defined in same basis
   if (s1%basis /= s2%basis) then
      write(6,*) 'Incompatible bases in symmetry_equal'
      stop
   endif

   symmetry_equal = .TRUE.

   ! Check equality of point group operations
   if (.not.matrix_equal(s1%m,s2%m)) then
       symmetry_equal = .FALSE.
       return
   endif

   ! Check equality of translations
   call shift_translation(s1,v1,R_basis,G_basis)
   call shift_translation(s2,v2,R_basis,G_basis)
   if (.not.equal(v1,v2)) then
       symmetry_equal = .FALSE.
       return
   endif

   end function symmetry_equal
   !===================================================================


   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! Operations on symmetries
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


   type(symmetry_type) function make_identity(basis)
   !--------------------------------------------------------------
   !  Constructor for identity symmetry element
   !--------------------------------------------------------------
   character(9) :: basis
   integer:: i
   make_identity%basis = basis
   make_identity%m     = 0.0_long
   make_identity%v     = 0.0_long
   do i=1, dim
      make_identity%m(i,i) = 1.0_long
   enddo
   end function make_identity
   !===================================================================


   type(symmetry_type) function make_inversion(basis)
   !-------------------------------------------------------------------
   !  Constructor for inversion symmetry element
   !-------------------------------------------------------------------
   character(9) :: basis
   integer:: i
   make_inversion%basis = basis
   make_inversion%m     = 0.0_long
   make_inversion%v     = 0.0_long
   do i=1, dim
      make_inversion%m(i,i) = -1.0_long
   enddo
   end function make_inversion
   !===================================================================


   subroutine shift_translation(s,v,R_basis,G_basis)
   !-------------------------------------------------------------------
   !  Shifts position or translation vector so as to lie in the
   !  first unit cell of a periodic structure 
   !-------------------------------------------------------------------
   type(symmetry_type), intent(IN) :: s
   real(long), intent(IN)    :: R_basis(:,:), G_basis(:,:)
   real(long), intent(OUT)   :: v(3)

   !  Local variables
   real(long)                :: coeff(3), twopi
   integer                   :: i, j

   if ( s%basis == Cartesian ) then 
   
      twopi = 4.0*acos(0.0)
      coeff = 0.0_long
      do i=1, dim
         ! coeff(i) = (s%v .dot. G_basis(i,:)) /twopi
         do j=1, dim
            coeff(i) = coeff(i) + G_basis(i,j)*s%v(j)
         enddo
         coeff(i) = coeff(i)/twopi
         coeff(i) = modulo(coeff(i),1.0_long)
      enddo
      v = 0.0_long
      do i=1, dim
         do j=1, dim
            v(i) = v(i) + coeff(j)*R_basis(j,i)
         enddo
      enddo
   
   else if (s%basis == Bravais  ) then 
   
      do i=1, dim
         v(i) = modulo(s%v(i),1.0_long)
      enddo
   
   else
   
      write(6,*) 'Invalid value of s%basis in shift_translation'
      stop
   
   endif
   
   end subroutine shift_translation
   !===================================================================


   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! Operations on groups
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


   type(group_type) function copy_of_group(in_group)
   !-------------------------------------------------------------------
   !  Returns a copy of a group, checking for consistency of bases
   !-------------------------------------------------------------------
   type(group_type), intent(IN) :: in_group
   integer      :: i
   character(9) :: basis
   basis = in_group%s(1)%basis
   do i=1, in_group%order
      copy_of_group%s(i)%m = in_group%s(i)%m
      copy_of_group%s(i)%v = in_group%s(i)%v
      if ( basis == in_group%s(i)%basis ) then
         copy_of_group%s(i)%basis = basis
      else
         write(6,*) 'Error: Inconsistent bases on input to copy_of_group'
         stop
      endif
   enddo
   end function copy_of_group
   !===================================================================


   type(group_type) function Brav_group(Cart_group,R_basis,G_basis)
   !-------------------------------------------------------------------
   !  Takes input Cartesian group, returns corresponding Bravais group
   !-------------------------------------------------------------------
   type(group_type), intent(IN) :: Cart_group
   real(long),  intent(IN)      :: R_basis(:,:), G_basis(:,:)
   integer     :: i
   real(long)  :: m(3,3), twopi

   twopi = 4.0_long*acos(0.0_long) 
   if ( Cart_group%s(1)%basis /= Cartesian) then
      write(6,*) 'Error: Input group not Cartesian in Brav_group'
      stop
   endif
   Brav_group%order = Cart_group%order
   do i = 1, Cart_group%order
      m = G_basis .dot. Cart_group%s(i)%m
      Brav_group%s(i)%m = ( m .dot. transpose(R_basis) )/twopi
      Brav_group%s(i)%v = ( G_basis .dot. Cart_group%s(i)%v )/twopi
      if ( Cart_group%s(i)%basis == Cartesian ) then
         Brav_group%s(i)%basis = Bravais
      else
         write(6,*) 'Error: Inconsistent bases on input to Brav_group'
         stop
      endif
   enddo
   end function Brav_group
   !===================================================================


   type(group_type) function Cart_group(Brav_group,R_basis,G_basis)
   !-------------------------------------------------------------------
   !  Takes input Bravais group, returns corresponding Cartesian group
   !-------------------------------------------------------------------
   type(group_type), intent(IN) :: Brav_group
   real(long),  intent(IN)      :: R_basis(:,:), G_basis(:,:)
   integer     :: i
   real(long)  :: m(3,3), twopi

   twopi = 4.0_long*acos(0.0_long)
   if ( Brav_group%s(1)%basis /= Bravais) then
      write(6,*) 'Error: Input group not Bravais in Cart_group'
      stop
   endif
   Cart_group%order = Brav_group%order
   do i = 1, Brav_group%order
      m = Brav_group%s(i)%m .dot. G_basis
      Cart_group%s(i)%m = ( transpose(R_basis) .dot. m )/twopi
      Cart_group%s(i)%v = Brav_group%s(i)%v .dot. R_basis
      if ( Brav_group%s(i)%basis == Bravais ) then
         Cart_group%s(i)%basis = Cartesian
      else
         write(6,*) 'Error: Inconsistent bases on input to Cart_group'
         stop
      endif
   enddo
   end function Cart_group
   !===================================================================


   subroutine make_table(g,table,R_basis,G_basis)
   !-------------------------------------------------------------------
   !  Routine to construct Cayley table for group
   !  Note: table%index(i,j) = 0 if element not in group
   !-------------------------------------------------------------------
   type(group_type), intent(IN) :: g
   real(long), intent(IN) :: R_basis(:,:), G_basis(:,:)
   type(table_type), intent(OUT):: table

   real(long) :: v(3)
   integer    :: i, j, k

   table%order = g%order
   do i=1, g%order
      do j=1, g%order

   !        Make symmetry element (i,j) of table
         table%s(i,j) = g%s(i) .dot. g%s(j) 

   !        Shift translation of element (i j)
         call shift_translation(table%s(i,j),v,R_basis,G_basis)
         table%s(i,j)%v = v

   !        Find index of symmetry element (i,j) in group
         table%index(i,j)  = 0
         do k=1, g%order
           if (equal(table%s(i,j),g%s(k),R_basis,G_basis)) then
              table%index(i,j) = k
           endif
         enddo

      enddo
   enddo
   end subroutine make_table
   !===================================================================

 
   !-------------------------------------------------------------------
   !****p* group_mod/make_group
   ! SUBROUTINE
   !    make_group(g,R_basis,G_basis)
   ! PURPOSE
   !    Construct complete group from incomplete proto-group, after 
   !    checking that all elements of proto-group have same basis
   ! ARGUMENTS
   !    group   - group, incomplete on input, complete on output
   !    R_basis - Bravais lattice basis vectors
   !    G_basis - reciprocal lattice basis vectors
   ! SOURCE
   !-------------------------------------------------------------------
   subroutine make_group(g,R_basis,G_basis)
   type(group_type), intent(INOUT) :: g
   real(long), intent(IN) :: R_basis(:,:), G_basis(:,:)
   !***
 
   type(symmetry_type) :: identity
   type(table_type)    :: table
   real(long)          :: v(3)
   integer             :: i,j,k,l, iterate
   logical             :: is_a_group, redundant
   character(9)        :: basis
 
   ! Check that all symmetry elements are in same, valid basis
   basis = g%s(1)%basis
   if (( basis /= Cartesian ).and.( basis /= Bravais )) then
       write(6,*) 'Invalid basis in make_group'
       stop
   endif
   do i=1, g%order
      if (g%s(i)%basis /= basis) then
         write(6,*) 'Incompatible bases for symmetries in make_group'
         stop
      endif
   enddo

   ! Check for validity of basis, consistency with group
   if (.not.valid_basis(g,R_basis,G_basis) ) then
      write(6,*) 'Invalid lattice basis in make_group'
      stop
   endif
 
   ! Check for presence of identity, add if not present
   identity = make_identity(basis)
   if ( in_group(identity,g,R_basis,G_basis) == 0) then
      g%order = g%order + 1
      g%s(g%order) = identity
      g%s(g%order)%basis = basis
   endif
 
   Do iterate=1, 20
 
     !write(6,*) "iteration #", iterate
     is_a_group = .true.
 
     ! Add missing inverses of existing elements to group
     k = g%order
     do i=1,k
        if ( index_of_inverse(g,i,R_basis,G_basis) == 0) then
           g%order = g%order + 1
           g%s(g%order) = inverse(g%s(i))
           call shift_translation(g%s(g%order),v,R_basis,G_basis)
           g%s(g%order)%v = v
           is_a_group = .false.
        endif
     enddo
 
     ! Add products of existing elements to group
     call make_table(g,table,R_basis,G_basis)
     k = g%order
     do i=1, k
        do j=1, k
          if ( table%index(i,j) == 0) then
              if ( in_group(table%s(i,j),g,R_basis,G_basis) == 0) then
                 g%order = g%order + 1
                 g%s(g%order) = table%s(i,j)
                 g%s(g%order)%basis = basis
                 is_a_group = .false.
              endif
           endif
        enddo
     enddo
  
     if (is_a_group) return
 
   enddo
   write(6,*) 'Iteration failed in make_group'
   stop
   
   end subroutine make_group
   !===================================================================


   logical function commute_group( A , g )
   !-----------------------------------------------------------------
   ! Returns true if matrix A and the group g commute, otherwise false
   !-----------------------------------------------------------------
   real(long), intent(IN)       :: A(:,:)
   type(group_type), intent(IN) :: g

   logical, dimension(g%order) :: test
   integer :: i

   do i = 1,g%order
      test(i) = equal( A .dot. g%s(i)%m, g%s(i)%m .dot. A )
   enddo
   commute_group = all(test)

   end function commute_group
   !===================================================================


   subroutine deformation_subgroup( def, g, R_basis, G_basis )
   !-------------------------------------------------------------------
   ! Determine subset of a group that commutes with deformation matrix
   ! On input,  g is the original group
   ! On output, g is equal to this subgroup
   !-------------------------------------------------------------------
   !Dummy Variables
   real(long), intent(IN), dimension(3,3) :: def ! deformation matrix
   real(long), intent(IN), dimension(3,3) :: R_basis, G_basis
   type(group_type), intent(INOUT)        :: g   ! group
   
   integer          :: i
   type(group_type) :: subg
   
   ! Check which elements of group commute with deformation 
   subg%order = 0
   do i = 1,g%order
      if ( equal(def .dot. g%s(i)%m ,g%s(i)%m .dot. def ) ) then
         subg%order         = subg%order + 1
         subg%s(subg%order) = g%s(i)
      endif
   enddo
   g = subg
   call make_group( g, R_basis, G_basis )
 
   end subroutine deformation_subgroup
   !===================================================================


   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! Operations that search for an element in a group
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


   integer function in_group(s,g,R_basis,G_basis)
   !-------------------------------------------------------------------
   ! Returns i if symmetry s is element i of group g
   ! Returns 0 if symmetry s is not in group g
   !-------------------------------------------------------------------
   type(symmetry_type), intent(IN) :: s
   type(group_type), intent(IN)    :: g
   real(long), intent(IN) :: R_basis(:,:), G_basis(:,:)
   integer :: i
   in_group = 0
   do i=1, g%order
      if ( equal(s,g%s(i),R_basis,G_basis) ) then
         in_group = i
         return
      endif
   enddo
   end function in_group
   !===================================================================
   
 
   integer function index_of_inverse(g,i,R_basis,G_basis)
   !-------------------------------------------------------------------
   ! Function returns index of first appearance of inverse of 
   ! symmetry s(i) of group g, if inverse(g%s(i)) is in g, and 
   ! returns 0 if inverse(g%s(i)) is not in g
   !-------------------------------------------------------------------
   type(group_type), intent(IN) :: g
   real(long), intent(IN) :: R_basis(:,:), G_basis(:,:)
   integer,     intent(IN) :: i
   integer :: j
   index_of_inverse = 0
   do j=1, g%order
      if ( equal(inverse(g%s(i)),g%s(j),R_basis,G_basis) ) then
          index_of_inverse = j
      endif
   enddo
   end function index_of_inverse
   !===================================================================
 
 
   logical function includes_inversion(g,R_basis,G_basis)
   !-------------------------------------------------------------------
   ! True if the group includes inversion, i.e., is centrosymmetric
   !-------------------------------------------------------------------
   type(group_type), intent(IN) :: g
   real(long),  intent(IN) :: R_basis(:,:), G_basis(:,:)
   type(symmetry_type) :: inversion
   inversion = make_inversion(g%s(1)%basis)
   includes_inversion = .true.
   if ( in_group(inversion,g,R_basis,G_basis) == 0 ) then
      includes_inversion = .false.
   endif
   end function includes_inversion
   !===================================================================


   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   !  Procedures that act on R_basis and G_basis
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



   logical function valid_basis(g,R_basis,G_basis)
   !-------------------------------------------------------------------
   !  Check if R_basis and G_basis form a valid set of basis vectors
   !-------------------------------------------------------------------
   type(group_type), intent(IN) :: g
   real(long),  intent(IN)      :: R_basis(:,:), G_basis(:,:)

   ! Local variables
   integer      :: i, j
   character(9) :: basis
   real(long)   :: m1(3,3), m2(3,3), twopi
   character(3) :: fmt_string

   type(group_type) :: g_Brav, g_Cart
   twopi = 4.0_long*acos(0.0_long) 


   valid_basis = .true.

   ! Check if R_basis and G_basis satisfy bi-orthogonality condition
   m1 = R_basis .dot. Transpose(G_basis) / twopi
   m2 = 0.0_long
   do i=1, dim
      m2(i,i) = 1.0_long
   enddo
   if (.not.equal(m1,m2)) then 
      valid_basis = .false.
      write(6,*) 'R_basis and G_basis are not bi-orthogonal'
      write(6,*)
      call output_matrix(m1,6) 
      write(6,*)
      call output_matrix(m2,6) 
      return
   endif

   ! Check rotation matrices contain integers in Bravais representation
   basis = g%s(1)%basis
   if ( basis == Bravais) then
      g_Brav = copy_of_group(g)
   else if ( basis == Cartesian ) then
      g_Brav = Brav_group(g,R_basis,G_basis)
   else
      write(6,*) 'Invalid group basis in valid_basis'
      stop
   endif
   m2 = 0.0_long
   do i=1, g%order
      m1 = g%s(i)%m - 1.0_long * nint( g_Brav%s(i)%m )
      if (.not.equal(m1,m2)) then
         valid_basis = .false.
         write(6,*) 'Elements of rotation', i, ' are not integers'
         return
      endif
   enddo

   ! Check rotation matrices are orthogonal in Cartesian representation
   basis = g%s(1)%basis
   if ( basis == Cartesian) then
      g_Cart = copy_of_group(g)
   else if ( basis == Bravais ) then
      g_Cart = Cart_group(g,R_basis,G_basis)
   else
      write(6,*) 'Invalid group basis in valid_basis'
      stop
   endif

   m2 = 0.0_long
   do i=1, dim
      m2(i,i) = 1.0_long
   enddo
   do i=1, g%order
      m1 = transpose(g_Cart%s(i)%m).dot.g_Cart%s(i)%m
      if (.not.equal(m1,m2)) then
         valid_basis = .false.
         write(6,*) 'Rotation ', i, ' is not an orthogonal matrix'
         write(6,'(3f12.4)') g%s(i)%m
         write(6,*)
         write(6,'(3f12.4)') g_Cart%s(i)%m
         write(6,*)
         write(unit=fmt_string,fmt='(i3)') dim
         write(6,'('//fmt_string//'f12.4)') R_basis(:dim,:dim)/R_Basis(1,1)
         write(6,'('//fmt_string//'f12.4)') R_basis
         return
      endif
   enddo

   end function valid_basis
   !===================================================================
 

 
   function lattice_type(g)
   !-------------------------------------------------------------------
   ! Determine the lattice type (cubic, trigonal, etc) given the group
   !-------------------------------------------------------------------
   character(len = 30) :: lattice_type
   type(group_type), intent(IN) :: g ! group

   ! Local Variables
   integer :: i

   real(long), dimension(dim,dim) :: def  ! deformation
   real(long), parameter:: epsilon = 0.05 ! small number for deformation
   real(long) :: theta ! angle
   real(long) :: pi

   pi = 4* atan(1.0_long)

   ! If group order is >= 48, cubic
   if (g%order >= 48 ) then
      lattice_type = 'cubic'
      return
   endif
    
   ! Check if triclinic
   ! If we can deform gamma, then it's triclinic
     def = 0.0_long
   do i = 1,dim
      def(i,i) = 1.0_long
   enddo
   def(1,2) = epsilon
   def(2,1) = epsilon
   if ( commute_group(def,g) ) then
      lattice_type = 'triclinic'
      return
   endif

   ! Check if Monoclinic or orthorhombic ( a and b independent )
   def =  0.0_long
   do i = 1,dim
      def(i,i) = 1.0_long
   enddo
   def(1,1) = def(1,1) + epsilon
    
   if ( commute_group(def,g) ) then
      ! Lattice is orthorhombic or monoclinic as a,b are independent
      ! Check if monoclinic by looking at deformation of angle beta
      def = 0.0_long
      do i = 1,dim
         def(i,i) = 1.0_long
      enddo
      def(1,3) = epsilon
      def(3,1) = epsilon

      if ( commute_group(def,g) ) then
         lattice_type = 'monoclinic'
         return
      else
         lattice_type = 'orthogonal'
         return
      endif
   endif

   ! Check if Tetragonal or Hexagonal ( change a and b together)
   def=0.0_long
   do i = 1,dim-1
      def(i,i) = 1.0_long + epsilon
   enddo
   def(dim,dim) = 1.0_long
    
   if ( commute_group(def,g) ) then
      ! Lattice is tetragonal or hexagonal
      ! Hexagonal point group has three-fold symmetry around the c-axis.  
      ! Check for it
      def = 0.0_long
      theta = 2.0_long * pi/3.0_long

      def(3,3) = 1.0_long
      def(:,1) = (/  cos(theta), sin(theta), 0.0_long /)
      def(:,2) = (/ -sin(theta), cos(theta), 0.0_long /)
       
      if (commute_group(def,g)) then
         lattice_type = 'hexagonal'
         return
      else
         lattice_type = 'tetragonal'
         return
      endif

   endif

   !If it isn't any of the above, it must be trigonal
   lattice_type = 'trigonal'
   return

   end function lattice_type
   !===================================================================


   logical function is_sub_group(sub_group,super_group,R_basis,G_basis)
   !------------------------------------------------------------------
   ! Returns true if all of sub_groups elements are in super_group
   !------------------------------------------------------------------
   type(group_type) :: sub_group, super_group
   real(long) :: R_basis(dim,dim), G_basis(dim,dim)
   
   integer :: i
   
   is_sub_group = .true.
   do i = 1, sub_group%order
      if( in_group(sub_group%s(i),super_group,R_basis,G_basis).ne.0) then
         cycle
      else
         is_sub_group = .false.
         exit
      endif
   enddo
   end function is_sub_group
   !===================================================================


end module group_mod
