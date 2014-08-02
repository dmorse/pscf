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
!****m scf/unit_cell_mod
! PURPOSE
!    Define crystal unit cell and lattice basis vectors
! AUTHOR
!    David Morse (2002)
!    Chris Tyler (2002-2004)
!  SOURCE
!----------------------------------------------------------------------
module unit_cell_mod
   use const_mod     
   use io_mod
   implicit none

   private

   ! Public procedures
   public :: input_unit_cell     ! read  dim, crystal_system, cell_param
   public :: output_unit_cell     ! write dim, crystal_system, cell_param
   public :: standard_cell_param  ! return parameters a,b,c, alpha,beta,gamma
   public :: make_unit_cell       ! create lattice basis vectors, etc.
   public :: define_unit_cell     ! reset cell parameters 
   public :: make_G_basis         ! make G_basis from R_basis

   ! Public variables
   public :: crystal_system, N_cell_param, cell_param
   public :: R_basis, G_basis, dGG_basis

   character(30) :: crystal_system   ! type of crystal cell (cubic, etc.)
   integer       :: N_cell_param     ! # of unit cell parameters
   real(long)    :: cell_param(6)    ! unit cell parameters
   real(long)    :: R_basis(:,:)     ! (dim,dim) lattice bases a_i
   real(long)    :: G_basis(:,:)     ! (dim,dim) reciprocal bases b_i
   real(long)    :: dGG_basis(:,:,:) ! (dim,dim,6) derivatives of b_i.b_j
   !***

   allocatable   :: R_basis, G_basis, dGG_basis

   ! Private variables
   real(long), allocatable :: dR_basis(:,:,:)  ! derivatives of a_i
   real(long), allocatable :: dG_basis(:,:,:)  ! derivatives of b_i
   real(long), allocatable :: dRR_basis(:,:,:) ! derivatives of a_i.a_j


   !--------------------------------------------------------------------
   !****v* unit_cell_mod/crystal_system
   ! VARIABLE
   ! character(30) crystal_system = string identifying crystal system
   !
   ! Allowed values:      
   !    3D crystal systems (dim = 3):
   !        'cubic', 'tetragonal', 'orthorhombic', 'monoclinic',  
   !        'hexagonal', 'triclinic'
   !    2D crystal systems (dim = 2)
   !        'square', 'rectangular', 'rhombus', 'hexagonal', 'oblique'      
   !    1D crystal system (dim = 1):
   !           'lamellar'      
   !*** ----------------------------------------------------------------
   !****v* unit_cell_mod/N_cell_param
   ! VARIABLE
   ! integer N_cell_param  = # of unit cell parameters 
   !                         Different values needed for different crystal 
   !                         systems (e.g., 1 for cubic, 6 for triclinic)
   !*** ----------------------------------------------------------------
   !****v* unit_cell_mod/cell_param
   ! VARIABLE
   ! real(long) cell_param(6) - array of cell parameters
   !                            Only elements 1:N_cell_param are used
   !*** ----------------------------------------------------------------
   !****v* unit_cell_mod/R_basis
   ! VARIABLE
   ! real(long) R_basis(:,:) - dimension(dim,dim)
   !            R_basis(i,:) = a_i = Bravais lattice basis vector i
   !*** ----------------------------------------------------------------
   !****v* unit_cell_mod/G_basis
   ! VARIABLE
   ! real(long) G_basis(:,:) - dimension(dim,dim)
   !            G_basis(i,:) = b_i = reciprocal lattice basis vector i
   !*** ----------------------------------------------------------------
   !****v* unit_cell_mod/dGG_basis
   ! VARIABLE
   ! real(long)   dGG_basis(:,:,:) - dimension(dim,dim,6)
   !              dGG_basis(i,j,k) = d(b_i.dot.b_j)/d(cell_param(k))
   !              Needed in calculation of stress by perturbation theory
   !*** ----------------------------------------------------------------

contains

   !---------------------------------------------------------------
   !****p* unit_cell_mod/input_unit_cell
   ! SUBROUTINE
   !    input_unit_cell(i_unit, fmt_flag)
   ! PURPOSE
   !    Read data needed to construct unit cell from file i_unit.
   !    Inputs dim, crystal_system, N_cell_param, and cell_param 
   !    If necessary, allocates R_basis, G_basis, related arrays
   ! ARGUMENTS
   !    integer      i_unit   - unit # of input file
   !    character(1) fmt_flag - flag specifying input format 
   ! COMMENT
   !    Allowed values of fmt_flag:
   !       F -> formatted ascii input
   !       U -> unformatted input
   ! SOURCE
   !---------------------------------------------------------------
   subroutine input_unit_cell(i_unit,fmt_flag)
   integer, intent(IN)            :: i_unit
   character(len = 1), intent(IN) :: fmt_flag
   !***

   call set_io_units(i=i_unit,o=6)
   select case(fmt_flag)
   case('F') ! Input formatted
      call input(dim,'dim')
      call input(crystal_system,'crystal_system')
      call input(N_cell_param,'N_cell_param')
      call input(cell_param,N_cell_param,'cell_param')
   case('U')
      read(i_unit) dim
      read(i_unit) crystal_system
      read(i_unit) N_cell_param
      read(i_unit) cell_param
   case default
      print *, 'Illegal format specified in input_unit_cell'
      print *, 'fmt_flag = ', fmt_flag
      stop
   end select
   if (.not.allocated(R_basis))   allocate(R_basis(dim,dim))
   if (.not.allocated(G_basis))   allocate(G_basis(dim,dim))
   if (.not.allocated(dR_basis))  allocate(dR_basis(dim,dim,6))
   if (.not.allocated(dG_basis))  allocate(dG_basis(dim,dim,6))
   if (.not.allocated(dRR_basis)) allocate(dRR_basis(dim,dim,6))
   if (.not.allocated(dGG_basis)) allocate(dGG_basis(dim,dim,6))

   end subroutine input_unit_cell
   !---------------------------------------------------------------


   !---------------------------------------------------------------
   !****p* unit_cell_mod/output_unit_cell
   ! SUBROUTINE
   !    output_unit_cell(o_unit,fmt_flag)
   ! PURPOSE
   !    Write crystal_system, N_cell_param, and cell_param to file 
   ! ARGUMENTS
   !    integer      o_unit   - unit # of output file
   !    character(1) fmt_flag - flag specifying output format 
   ! COMMENT
   !    Allowed values of fmt_flag:
   !        F  -> formatted output
   !        U  -> unformatted output
   ! SOURCE
   !---------------------------------------------------------------
   subroutine output_unit_cell(o_unit,fmt_flag)
   integer, intent(IN)            :: o_unit
   character(len = 1), intent(IN) :: fmt_flag
   !***

   integer :: k
  
   !call set_io_units(o=o_unit)
   select case(fmt_flag)
   case('F') ! Formatted for input
      call output(dim,'dim',o=o_unit)
      call output(trim(crystal_system),'crystal_system',o=o_unit)
      call output(N_cell_param,'N_cell_param',o=o_unit)
      call output(cell_param,N_cell_param,'cell_param',o=o_unit)
   case('U')
      write(o_unit) dim
      write(o_unit) crystal_system
      write(o_unit) N_cell_param
      write(o_unit) cell_param
      write(o_unit) R_basis
      write(o_unit) G_basis
   case default
      print *, 'Invalid fmt_flag in output_unit_cell'
      print *, 'fmt_flag = ', fmt_flag
      stop
   end select

   end subroutine output_unit_cell
   !------------------------------------------------------------------


   !---------------------------------------------------------------
   !****p* unit_cell_mod/standard_cell_param
   ! FUNCTION
   !    standard_cell_param(cell_param)
   ! PURPOSE
   !    Returns array (a, b, c, alpha, beta, gamma) for 3-d systems
   !      a, b, c are lengths of the three Bravais basis vectors
   !      alpha is the angle beween b, c
   !      beta is the angle between a, c
   !      gamma is the angle between a, b
   ! RETURN
   !    standard_cell_param(1:3) = (a,b,c)
   !    standard_cell_param(4:6) = (alpha,beta,gamma)
   ! AUTHOR
   !    Chris Tyler
   ! SOURCE
   !---------------------------------------------------------------
   function standard_cell_param(cell_param)
   real(long), dimension(6), intent(IN) :: cell_param
   real(long), dimension(6) :: standard_cell_param
   !***

   real(long) :: a,b,c,alpha,beta,gamma
   
   if ( dim .ne. 3 ) then
      standard_cell_param = cell_param
      return
   endif

   select case(crystal_system)
   case('cubic')
      a = cell_param(1)
      b = cell_param(1)
      c = cell_param(1)
      alpha = 90.0
      beta  = 90.0
      gamma = 90.0
   case('tetragonal')
      a = cell_param(1)
      b = cell_param(1)
      c = cell_param(2)
      alpha = 90.0
      beta = 90.0
      gamma = 90.0
   case('orthorhombic')
      alpha = 90.0
      beta = 90.0 
      gamma = 90.0
      a = cell_param(1)
      b = cell_param(2)
      c = cell_param(3)
   case('hexagonal')
      a = cell_param(1)
      b = cell_param(1)
      c = cell_param(2)
      gamma = 120.0
      beta = 90.0
      alpha = 90.0
   case('trigonal')
      a = cell_param(1)
      b = cell_param(1)
      c = cell_param(1)
      alpha = cell_param(2) * 90/asin(1.0)
      beta = alpha
      gamma = alpha
   case('monoclinic')
      a = cell_param(1)
      b = cell_param(2)
      c = cell_param(3)
      alpha = 90.0
      beta = cell_param(4)
      gamma = 90.0
   case('triclinic')
      a = cell_param(1)
      b = cell_param(2)
      c = cell_param(3)
      alpha = cell_param(4)
      beta = cell_param(5)
      gamma = cell_param(6)
   case default
      a = 1
      b = 1
      c = 1
      alpha = 90
      beta = 90
      gamma = 90
   end select

   standard_cell_param(1) = a
   standard_cell_param(2) = b
   standard_cell_param(3) = c
   standard_cell_param(4) = alpha
   standard_cell_param(5) = beta
   standard_cell_param(6) = gamma

   end function standard_cell_param
   !---------------------------------------------------------------


   !---------------------------------------------------------------
   !****p unit_cell_mod/make_unit_cell
   ! SUBROUTINE
   !    make_unit_cell 
   ! PURPOSE
   !    Constructs Bravais and reciprocal lattice vectors, and
   !    related arrays, from knowledge of module input variables.
   ! COMMENT
   !    All inputs and outputs are module variables, rather than
   !    explicit parameters. 
   !    Inputs:  crystal_system, N_cell_param, and cell_param
   !    Outputs: R_basis, G_basis, dRR_basis, dGG_basis
   ! SOURCE
   !---------------------------------------------------------------
   subroutine make_unit_cell 
   !***

   integer     :: i,j,k,l,m
   real(long)  :: a, b, c, alpha, beta, gamma, twopi

   !C if ( size(cell_param) < N_cell_param ) then
   !C  print *,'Error: size(cell_param)<N_cell_param in make_unit_cell'
   !C  stop
   !C endif

   twopi = 4.0_long*acos(0.0_long)

   R_basis   = 0.0_long
   G_basis   = 0.0_long
   dR_basis  = 0.0_long
   dG_basis  = 0.0_long
   dRR_basis = 0.0_long
   dGG_basis = 0.0_long

   select case(dim)
   case(3)

      select case(trim(crystal_system))
      case('cubic')
         If (N_cell_param /= 1) then
            write(6,*) 'Incorrect N_cell_param'
            stop
         endif
         a = cell_param(1)
         R_basis(1,1) = a
         R_basis(2,2) = a
         R_basis(3,3) = a
         dR_basis(1,1,1) = 1.0_long
         dR_basis(2,2,1) = 1.0_long
         dR_basis(3,3,1) = 1.0_long
      case('tetragonal')
         If (N_cell_param /= 2) then
            write(6,*) 'Incorrect N_cell_param'
            stop
         endif
         a = cell_param(1)
         c = cell_param(2)
         R_basis(1,1) = a
         R_basis(2,2) = a
         R_basis(3,3) = c
         dR_basis(1,1,1) = 1.0_long
         dR_basis(2,2,1) = 1.0_long
         dR_basis(3,3,2) = 1.0_long
      case('orthorhombic')
         If (N_cell_param /= 3) then
            write(6,*) 'Incorrect N_cell_param'
            stop
         endif
         a = cell_param(1)
         b = cell_param(2)
         c = cell_param(3)
         R_basis(1,1) = a
         R_basis(2,2) = b
         R_basis(3,3) = c
         dR_basis(1,1,1) = 1.0_long
         dR_basis(2,2,2) = 1.0_long
         dR_basis(3,3,3) = 1.0_long
      case('hexagonal')
         ! Note: Unique axis is c axis
         If (N_cell_param /= 2) then
            write(6,*) 'Incorrect N_cell_param'
            stop
         endif
         a = cell_param(1)
         c = cell_param(2)
         R_basis(1,1) = a
         R_basis(2,1) = -0.5_long*a
         R_basis(2,2) = a * sqrt(0.75_long)
         R_basis(3,3) = c
         dR_basis(1,1,1) = 1.0_long
         dR_basis(2,1,1) = -0.5_long
         dR_basis(2,2,1) = sqrt(0.75_long)
         dR_basis(3,3,2) = 1.0_long
      case('trigonal')
         !For Rhombohedral axes, otherwise use Hexagonal
         If (N_cell_param /= 2) then
            write(6,*) 'Incorrect N_cell_param'
            stop
         endif
         a = cell_param(1)    ! length of one edge
         beta = cell_param(2) ! angle between edges
         ! gamma is angle of rotation from a-b plane up to c-axis
         gamma = cos(beta)/cos(beta* 0.5_long)
         gamma = acos(gamma)
         R_basis(1,1) = a
         R_basis(2,1) = a * cos(beta)
         R_basis(2,2) = a * sin(beta)
         R_basis(3,1) = a * cos(gamma) * cos(beta*0.5_long)
         R_basis(3,2) = a * cos(gamma) * sin(beta*0.5_long)
         R_basis(3,3) = a * sin(gamma)
         dR_basis(1,1,1) = 1.0_long
         dR_basis(2,1,1) = cos(beta)
         dR_basis(2,2,1) = sin(beta)
         dR_basis(3,1,1) = cos(gamma) * cos(beta*0.5_long)
         dR_basis(3,2,1) = cos(gamma) * sin(beta*0.5_long)
         dR_basis(3,3,1) = sin(gamma)
         dR_basis(2,1,2) = -a*sin(beta)
         dR_basis(2,2,2) = a*cos(beta)
         ! alpha =d gamma/ d beta
         alpha = 2._long* sin(beta) - cos(beta)* tan(beta*0.5_long) *0.5_long &
              / ( cos(beta*0.5) &
                 * sqrt( (1 + 2._long* cos(beta))* (tan(beta*0.5)**2))  )
         dR_basis(3,1,2) = a * (-0.5_long * cos(gamma) * sin(beta*0.5_long) - &
                                 sin(gamma) * alpha * cos(beta*0.5_long))
         dR_basis(3,2,2) = a * ( 0.5_long * cos(beta*0.5_long) * cos(gamma) - &
                                 sin(gamma) * sin(beta*0.5_long) * alpha )
         dR_basis(3,3,2) = 1 * cos(gamma) * alpha
      case('monoclinic')
         ! Note: Unique axis is b axis
         If (N_cell_param /= 4) then
            write(6,*) 'Incorrect N_cell_param'
            stop
         endif
         a    = cell_param(1)
         b    = cell_param(2)
         c    = cell_param(3)
         beta = cell_param(4)
         R_basis(1,1) = a
         R_basis(2,2) = b
         R_basis(3,1) = c*cos(beta)
         R_basis(3,3) = c*sin(beta)
         dR_basis(1,1,1) = 1.0_long
         dR_basis(2,2,2) = 1.0_long
         dR_basis(3,1,3) = cos(beta)
         dR_basis(3,3,3) = sin(beta)
         dR_basis(3,1,4) = -c * sin(beta)
         dR_basis(3,3,4) =  c * cos(beta)
      case('triclinic')
         If (N_cell_param /= 6) then
            write(6,*) 'Incorrect N_cell_param'
            stop
         endif
         a = cell_param(1)
         b = cell_param(2)
         c = cell_param(3)
         alpha = cell_param(4) ! angle between c and x-y-plane
         beta = cell_param(5)  ! angle between c and z-axis
         gamma = cell_param(6) ! angle between a and b
         R_basis(1,1) = a
         R_basis(2,1) = b * cos(gamma)
         R_basis(2,2) = b * sin(gamma)
         R_basis(3,1) = c * cos(alpha)*sin(beta)
         R_basis(3,2) = c * sin(alpha)*sin(beta)
         R_basis(3,3) = c * cos(beta)
         dR_basis(1,1,1) = 1.0_long
         dR_basis(2,1,2) = cos(gamma)
         dR_basis(2,1,2) = sin(gamma)
         dR_basis(3,1,3) = cos(alpha) * sin(beta)
         dR_basis(3,2,3) = sin(alpha) * sin(beta)
         dR_basis(3,3,3) = cos(beta)
         dR_basis(3,1,4) = - c * sin(alpha) * sin(beta)
         dR_basis(3,2,4) = c * cos(alpha) * sin(beta)
         dR_basis(3,3,5) = - c * sin(beta)
         dR_basis(2,1,6) = - b * sin(gamma)
         dR_basis(2,2,6) = b * cos(gamma)
      case('R-3m')
         !For Rhombohedral axes, otherwise use Hexagonal
         If (N_cell_param /= 2) then
            write(6,*) 'Incorrect N_cell_param'
            stop
         endif
         a = cell_param(1)    ! length of one edge
         beta = cell_param(2) ! angle between edges
         ! gamma is angle of rotation from a-b plane up to c-axis
         gamma = cos(beta)/cos(beta* 0.5_long)
         gamma = acos(gamma)
         R_basis(1,1) = a
         R_basis(2,1) = a * cos(beta)
         R_basis(2,2) = a * sin(beta)
         R_basis(3,1) = a * cos(gamma) * cos(beta*0.5_long)
         R_basis(3,2) = a * cos(gamma) * sin(beta*0.5_long)
         R_basis(3,3) = a * sin(gamma)
         dR_basis(1,1,1) = 1.0_long
         dR_basis(2,1,1) = cos(beta)
         dR_basis(2,2,1) = sin(beta)
         dR_basis(3,1,1) = cos(gamma) * cos(beta*0.5_long)
         dR_basis(3,2,1) = cos(gamma) * sin(beta*0.5_long)
         dR_basis(3,3,1) = sin(gamma)
         N_cell_param=1
      case('pnna')
         if (N_cell_param /= 1 ) then
            write(6,*) 'Incorrect N_cell_param'
            stop
         endif
         a = cell_param(1)
         R_basis(1,1) = 2.0_long * a
         R_basis(2,2) = sqrt(3.0_long) *a
         R_basis(3,3) = 1.0_long * a
         dR_basis(1,1,1) = 2.0_long
         dR_basis(2,2,1) = sqrt(3.0_long)
         dR_basis(3,3,1) = 1.0_long
      case('fddd1')
         if (N_cell_param /= 1 ) then
            write(6,*) 'Incorrect N_cell_param'
            stop
         endif
         a = cell_param(1)
         R_basis(1,1) = a 
         R_basis(2,2) = sqrt(3.0_long) * a
         R_basis(3,3) = sqrt(3.0_long) * 2.0_long * a
         dR_basis(1,1,1) = 1
         dR_basis(2,2,1) = sqrt(3.0_long)
         dR_basis(3,3,1) = 2 * sqrt(3.0_long)
      case('fddd2')
         if (N_cell_param /= 2 ) then
            write(6,*) 'Incorrect N_cell_param'
            stop
         endif
         a = cell_param(1)
         alpha = cell_param(2) * atan(1.0)/45.0_long
         R_basis(1,1) = 2.0_long * sin(alpha/2) * a
         R_Basis(2,2) = 2.0_long * cos(alpha/2) * a
         R_basis(3,3) = sqrt(3.0_long) * 2.0_long * a
         dR_basis(1,1,1) = 2.0_long * sin(alpha/2)
         dR_basis(2,2,1) = 2.0_long * cos(alpha/2)
         dR_basis(3,3,1) = 2.0_long * sqrt(3.0_long)
         dR_basis(1,1,2) = cos(alpha/2)
         dR_basis(2,2,2) = -sin(alpha/2)
      case default
         write(6,*) 'Unknown crystal system, dim=3'
         stop
      end select

   case(2)

      select case(trim(crystal_system))
      case('square')
         if (N_cell_param /= 1 ) then
            write(6,*) 'Incorrect N_cell_param'
            stop
         endif
         a = cell_param(1)
         R_basis(1,1) = a
         R_basis(2,2) = a
         dR_basis(1,1,1) = 1.0_long
         dR_basis(2,2,1) = 1.0_long
      case('rectangular')
         if (N_cell_param /=2) then
            write(6,*) 'Incorrect N_cell_param'
            stop
         endif
         a = cell_param(1)
         b = cell_param(2)
         R_basis(1,1) = a
         R_basis(2,2) = b
         dR_basis(1,1,1) = 1.0_long
         dR_basis(2,2,2) = 1.0_long
      case('hexagonal')
         if (N_cell_param /= 1) then
            write(6,*) 'Incorrect N_cell_param'
            stop
         endif
         a = cell_param(1)
         R_basis(1,1) = a
         R_basis(2,1) = -0.5_long*a
         R_basis(2,2) = a * sqrt(0.75_long)
         dR_basis(1,1,1) = 1.0_long
         dR_basis(2,1,1) = -0.5_long
         dR_basis(2,2,1) = sqrt(0.75_long)
      case('oblique')
         if (N_cell_param /=3 ) then
            write(6,*) 'Incorrect N_cell_param'
            stop
         endif
         a = cell_param(1)
         b = cell_param(2)
         gamma = cell_param(3)
         R_basis(1,1) = a
         R_basis(2,1) = b * cos(gamma)
         R_basis(2,2) = b * sin(gamma)
         dR_basis(1,1,1) = a
         dR_basis(2,1,2) = cos(gamma)
         dR_basis(2,2,2) = sin(gamma)
         dR_basis(2,1,3) = - b * sin(gamma)
         dR_basis(2,2,3) = b * cos(gamma)
      case default
         write(6,*) 'Unknown crystal system, dim=2'
         stop
      end select

   case(1) ! 1D crystal system - lamellar

      if ( trim(crystal_system) == 'lamellar') then
         a = cell_param(1)
         R_basis(1,1) = a
         dR_basis(1,1,1) = 1.0_long
      else
         write(6,*) 'Unknown crystal system, dim=2'
         write(6,*) 'Only 1D system is "lamellar"'
         stop
      endif

   case default
         write(6,*) 'Invalid dimension, dim=', dim
         stop
   end select

   ! Invert R_basis to make G_basis
   call make_G_basis(R_basis,G_basis)

   ! Calculate dG_basis
   do k=1, N_cell_param
      do i=1, dim
         do j=1, dim
            do l=1, dim
               do m=1, dim
                  dG_basis(i,j,k) = dG_basis(i,j,k) &
                      - G_basis(i,l)*dR_basis(m,l,k)*G_basis(m,j)
               enddo
            enddo
         enddo
      enddo
   enddo
   dG_basis = dG_basis/twopi

   ! Calculate dRR_basis, dGG_basis
   !  do k=1, N_cell_param
   !     do i=1, dim
   !        do j=1, dim
   !           do l=1, dim
   !              dRR_basis(i,j,k) = dRR_basis(i,j,k)  &
   !                               + R_basis(i,l)*dR_basis(j,l,k)
   !              dGG_basis(i,j,k) = dGG_basis(i,j,k)  &
   !                               + G_basis(i,l)*dG_basis(j,l,k)
   !           enddo
   !           dRR_basis(i,j,k) = dRR_basis(i,j,k) + dRR_basis(j,i,k)
   !           dGG_basis(i,j,k) = dGG_basis(i,j,k) + dGG_basis(j,i,k)
   !        enddo
   !     enddo
   !  enddo


   do k = 1,N_cell_param
      do i = 1,dim
         do j = 1,dim
            do l = 1,dim
               dRR_basis(i,j,k) = dRR_basis(i,j,k) &
                                + R_basis(i,l) * dR_basis(l,j,k) &
                                + R_basis(j,l) * dR_basis(l,i,k)
               dGG_basis(i,j,k) = dGG_basis(i,j,k) &
                                + G_basis(i,l) * dG_basis(j,l,k) &
                                + G_basis(j,l) * dG_basis(i,l,k)
            enddo
         enddo
      enddo
   enddo

   !dGG_basis = -dGG_basis/twopi

   end subroutine make_unit_cell
   !===================================================================
   

   !---------------------------------------------------------------
   !****p* unit_cell_mod/define_unit_cell
   ! SUBROUTINE
   !    define_unit_cell( mylattice, my_N, my_param )
   ! PURPOSE
   !    Modify crystal system and/or unit cell parameters
   ! ARGUMENTS
   !    mylattice - crystal system
   !    my_N      - number of cell parameters 
   !    my_param  - array of cell parameters
   ! AUTHOR
   !    Chris Tyler
   ! SOURCE
   !---------------------------------------------------------------
   subroutine define_unit_cell( mylattice, my_N, my_param )
   character(*), intent(IN) :: mylattice
   integer, intent(IN)      :: my_N
   real(long), intent(IN)   :: my_param(:)
   !***

   integer :: i 

   if ( size(my_param) < my_N ) then
      print *, 'Error: size(my_param) < my_N in define_unit_cell'
      stop
   endif

   crystal_system = mylattice
   N_cell_param = my_N

   cell_param = 0.0_long
   do i = 1,N_cell_param
      cell_param(i) = my_param(i)
   enddo

   allocate(R_basis(dim,dim))
   allocate(G_basis(dim,dim))
   allocate(dR_basis(dim,dim,6))
   allocate(dG_basis(dim,dim,6))
   allocate(dRR_basis(dim,dim,6))
   allocate(dGG_basis(dim,dim,6))

   end subroutine define_unit_cell
   !==================================================================

   
   !-------------------------------------------------------------------
   !****p* unit_cell_mod/make_G_basis
   ! SUBROUTINE 
   !    make_G_basis(R_basis,G_basis)
   ! PURPOSE
   !    Construct array G_basis of reciprocal lattice basis vectors 
   !    from input array R_basis of Bravais lattice basis vectors 
   ! SOURCE
   !-------------------------------------------------------------------
   subroutine make_G_basis(R_basis,G_basis)
   use group_mod, only : Inverse
   real(long),  intent(IN)  :: R_basis(:,:)   ! (dim,dim) Bravais
   real(long),  intent(OUT) :: G_basis(:,:)   ! (dim,dim) reciprocal
   !***

   real(long)               :: R_local(3,3), G_local(3,3), twopi
   integer                  :: i, j
   twopi = 4.0_long*acos(0.0_long)

   ! Check dimensions for R_basis and G_basis
   if ( ( size(R_basis,1) /= dim).or.( size(R_basis,2) /= dim ) ) then
      write(6,*) 'Error in make_G_basis: Incorrect dimensions for R_basis'
   endif
   if ((size(G_basis,1)/=dim).or.(size(G_basis,2)/=dim)) then
      write(6,*) 'Error in make_G_basis: Incorrect dimensions for G_basis'
   endif

   R_local = 0.0_long
   G_local = 0.0_long
   do i=1, dim
      do j=1, dim
         R_local(i,j) = R_basis(i,j)
      enddo
   enddo
   R_local = inverse(R_local)         
   G_local = twopi*Transpose(R_local) ! Line split to compile on Regatta
   do i=1, dim
      do j=1, dim
         G_basis(i,j) = G_local(i,j)
      enddo
   enddo

   end subroutine make_G_basis
   !===================================================================

end module unit_cell_mod
