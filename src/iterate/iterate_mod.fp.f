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
!****m scf/iterate_mod
! MODULE
!   iterate_mod
! PURPOSE
!   Iteration of SCF equations
! ALGORITHM
!   Newton/Broyden algorithm that uses an approximate calculation 
!   for the initial Jacobian and Broyden updates to the inverse 
!   Jacobian. The approximate initial Jacobian is calculatd from 
!   a corresponding approximate linear response matrix, which is 
!   calculated in response_pd_mod. See the docblock for N_cut for 
!   a discussion of approximate response function. 
!
!   Subroutine iterate_NR can either:
!     i)  Iterate omega in a fixed unit cell     (domain = F)
!     ii) Iterate both omega and unit_cell_param (domain = T)
!   depending upon value of logical module variable domain
! AUTHOR
!   Chris Tyler - multiblock melt version
!   Amit Ranjan - multi-component blend version (2003)
!   David Morse - modifications (1/2004)
!   Jian Qin    - Broyden and backtracking (2006-2007)
! SOURCE
!-----------------------------------------------------------------------
module iterate_mod
   use const_mod
   use io_mod
   use response_pd_mod
   implicit none

   private

   ! public procedures
   public :: input_iterate_param   ! input iteration parameters 
   public :: output_iterate_param  ! output iteration variables
   public :: iterate_NR_startup    ! allocates arrays needed by iterate_NR
   public :: iterate_NR            ! main Newton-Raphson iteration routine

   ! public variables
   public :: max_itr    ! maximum allowable number of iterations
   public :: error_max  ! error tolerance: stop when error < error_max
   public :: domain     ! If true, adjust unit_cell_dimensions
   public :: itr_algo   ! Iteration algorithm 
   public :: N_cut      ! dimension of cutoff response and Jacobian
   !***

   !--------------------------------------------------------------------
   !****v iterate_mod/max_itr 
   ! VARIABLE
   !   integer max_itr - maximum allowed number of iterations
   !*** ----------------------------------------------------------------
   !****v iterate_mod/error_max 
   ! VARIABLE
   !   real(long) error_max - maximum error in converged solution
   !                        - Stop when error < error_max
   !*** ----------------------------------------------------------------
   !****v iterate_mod/domain 
   ! VARIABLE
   !   logical domain - If .true.,  iterate variable unit_cell shape
   !                    If .false., iterate omega in fixed unit cell
   !*** ----------------------------------------------------------------
   !****v iterate_mod/iter_algo
   ! VARIABLE
   !   character(10) iter_algo - Iteration algorithm to be used
   !                             For now, must equal 'NR' (Newton-Raphson)
   ! PURPOSE
   !   This variable will eventually allow the user to select from more
   !   than one iteration algorithm. Thus far, however, the only algorithm
   !   that has been implemented is a quasi-Newton / Broyden algorithm, 
   !   for which the value of iter_algo should be 'NR'. Implementation of
   !   a more robust relaxation algorithm is planned. 
   !*** ----------------------------------------------------------------
   !****v iterate_mod/N_cut 
   ! VARIABLE
   !   integer N_cut  - dimension of cutoff response and Jacobian
   ! PURPOSE
   !   The approximate linear response matrix used to construct the
   !   initial Jacobian is obtained by numerically calculating the
   !   response to the first N_cut basis functions to obtain an
   !   N_cut x N_cut submatrix within the subspace spanned by these
   !   functions, and using the analytically calculated response of
   !   a homogenous gas to approximate the response to the remaining
   !   basis functions.
   !*** ----------------------------------------------------------------

   ! public variable declarations
   integer       ::  max_itr       ! maximum # of iterations
   real(long)    ::  error_max     ! error tolerance 
   logical       ::  domain        ! true if domain iteration
   character(10) ::  itr_algo      ! Iteraction algorithm string. For
                                   ! now, must be 'NR' = Newton-Raphson
   integer       ::  N_cut         ! dimension of cutoff response matrix

   ! Private variable declarations
   integer                 :: N_residual      ! number of residuals
   integer                 :: N_corner        ! cut-off dimension of Jacobian
   integer, allocatable    :: ipvt(:)         ! Pivots for LU decomp 
   integer                 :: lwork           ! workspace dimension of inverse
   real(long), allocatable :: work(:)         ! workspace for matrix inverse
   real(long), allocatable :: residual(:)     ! residual array (N_residual)
   real(long), allocatable :: delta(:)        ! change in variables (N_residual)
   real(long), allocatable :: Jacobian(:,:)   ! Jacobian matrix
   real(long), allocatable :: J_corner(:,:)   ! Jacobian matrix
   logical                 :: Jacobian_empty  ! if true, Jacobian is empty
   logical                 :: full_inversion  ! if true, full Jacobian is inverted
 
contains

   !--------------------------------------------------------------------
   !****p iterate_mod/input_iterate_param
   ! SUBROUTINE
   !   input_iterate_param
   ! PURPOSE
   !   Read parameters max_itr, error_max, domain, itr_algo, N_cut
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine input_iterate_param
   !***
   call input(max_itr,'max_itr')      ! max # of iterations
   call input(error_max,'error_max')  ! error tolerance
   call input(domain,'domain')        ! T -> domain iteration
   call input(itr_algo,'itr_algo')    ! Iteration algorithm
   call input(N_cut,'N_cut')
   end subroutine input_iterate_param
   !====================================================================


   !--------------------------------------------------------------------
   !****p iterate_mod/output_iterate_param
   ! SUBROUTINE
   !   output_iterate_param
   ! PURPOSE
   !   Read parameters max_itr, error_max, domain, N_cut
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine output_iterate_param
   !***
   call output(max_itr,'max_itr')      ! max # of iterations
   call output(error_max,'error_max')  ! error tolerance
   call output(domain,'domain')        ! T -> domain iteration
   call output(itr_algo,'itr_algo')    ! Iteration algorithm
   call output(N_cut,'N_cut')
   end subroutine output_iterate_param
   !====================================================================


   !--------------------------------------------------------------------
   !****p iterate_mod/iterate_NR_startup
   ! SUBROUTINE
   !    iterate_NR_startup(N)
   ! PURPOSE
   !   1) Set integer module variable N_residual = # of residuals
   !   2) Allocate arrays needed by iterate_NR
   ! ARGUMENTS
   !   integer N      = # of basis functions
   ! SOURCE
   !---------------------------------------------------------------------
   subroutine iterate_NR_startup(N)
   use unit_cell_mod, only : N_cell_param
   use chemistry_mod, only : N_monomer, ensemble
   integer, intent(IN)  :: N       ! # of basis functions
   !***

   integer              :: info    ! message variable for dgetri

   if (N_cut > N) then
       N_cut = N
   endif
   if (domain) then
      N_residual  = N_monomer*( N - 1 + ensemble ) + N_cell_param
      if( N_cut > 0 ) then
         N_corner = N_monomer*(N_cut- 1 + ensemble ) + N_cell_param
      else
         N_corner = N_cell_param
      endif
   else
      N_residual  = N_monomer*( N - 1 + ensemble )
      if( N_cut > 0 ) then
         N_corner = N_monomer*(N_cut- 1 + ensemble )
      else
         N_corner = 0
      endif
   endif

   if (allocated(ipvt)) deallocate(ipvt)
   if (allocated(work)) deallocate(work)
   if (allocated(residual)) deallocate(residual)
   if (allocated(delta)) deallocate(delta)
   if (allocated(Jacobian)) deallocate(Jacobian)
   if (allocated(J_corner)) deallocate(J_corner)

   allocate(ipvt(N_residual))
   allocate(work(N_residual))
   allocate(residual(N_residual))
   allocate(delta(N_residual))
   allocate(Jacobian(N_residual,N_residual))
   allocate(J_corner(N_corner,N_corner))

   call init_response_pd(N_monomer,N)

   ! estimate the opitimal workspace dimension
   lwork = -1
   call dgetri(N_residual,Jacobian,N_residual,ipvt,work,lwork,info)
   lwork = work(1)
   if (allocated(work)) deallocate(work)
   allocate(work(lwork))

   Jacobian_empty = .true.   ! reset to true after 1st calculation
   full_inversion = .false.  ! not invert full Jacobian, see invert_Jacobian

   end subroutine iterate_NR_startup
   !======================================================================


   !---------------------------------------------------------------------
   !****p iterate_mod/iterate_NR
   ! SUBROUTINE
   !    subroutine iterate_NR
   ! PURPOSE
   !    Newton-Raphson iteration of the SCFT, using numerically
   !    calculated Jacobian
   ! COMMENT
   !   1) Algorithm assumes that density_startup, iterate_NR_startup,
   !      and make_unit_cell have been called prior to iterate_NR.
   !
   !   2) Routine will simultaneously iterates omega field and unit 
   !      cell parameters if domain = .true., or iterate the omega
   !      field for a fixed unit cell if domain = .false.
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine iterate_NR(  &
              N,           &!  # of basis functions
              omega,       &!  chemical potential field (IN/OUT)
              itr,         &!  actual number of interations
              converge,    &!  = .true. if converged
              error,       &!  final error = max(residuals)
              rho,         &!  monomer density field
              f_Helmholtz, &!  Helmholtz free energy per monomer / kT
              !# ifdef DEVEL
              f_component, &!  different contribution to f
              overlap,     &!  overlap integrals
              !# endif
              pressure,    &!  pressure * monomer volume / kT
              stress       &!  d(free energy)/d(cell parameters)
              )
   !--------------------------------------------------------------------

   use io_mod, only        : output
   use unit_cell_mod, only : N_cell_param, cell_param, &
                             make_unit_cell, G_basis, dGG_basis
   use basis_mod, only     : make_dGsq
   use chemistry_mod, only : N_monomer, ensemble, mu_chain, mu_solvent, &
                             phi_solvent, phi_chain, N_chain, N_solvent 
   !# ifdef DEVEL
   use scf_mod, only       : density, scf_stress, set_omega_uniform, &
                             mu_phi, free_energy, divide_energy
   !# else
   use scf_mod, only       : density, scf_stress, set_omega_uniform, &
                             mu_phi, free_energy
   !# endif
   use grid_mod, only      : make_ksq

   ! Arguments
   integer, intent(IN)         :: N        
   real(long), intent(INOUT)   :: omega(:,:)  ! omega(N_monomer,N)
   integer, intent(OUT)        :: itr         
   logical, intent(OUT)        :: converge
   real(long), intent(OUT)     :: error
   real(long), intent(OUT)     :: rho(:,:)    ! rho(N_monomer,N)
   real(long), intent(OUT)     :: f_Helmholtz
   real(long), intent(OUT)     :: pressure
   real(long), intent(OUT)     :: stress(:)   ! stress(N_cell_param)
   !# ifdef DEVEL
   real(long), intent(OUT)     :: f_component(:)
   real(long), intent(OUT)     :: overlap(:,:)
   !# endif
   !***

   ! parameter
   real(long), parameter :: large = 100.0

   ! local Variables
   integer    :: i, k, alpha, M, info
   integer    :: error_index
   real(long) :: dGsq(N,N_cell_param), q(N_chain), q_solvent(N_solvent)
   real(long) :: old_error
   real(long) :: delta_norm, u(N_residual), v(N_residual), uv

   ! parameter and variables for backtracking line search
   ! backtracking algorithm is as described in Numerical Recipes
   real(long), parameter :: STPMX=100.0_long
   real(long)  :: ow(N_monomer,N)    ! old omega fields
   real(long)  :: ocell(N_cell_param)! old cell param
   real(long)  :: fold, fnew         ! f=x^2/2, to be optimized
   real(long)  :: p(N_residual)      ! update step direction
   real(long)  :: stpmax             ! maximum update step
   real(long)  :: slope              ! line search direction
   real(long)  :: r_old(N_residual)  ! old residual
   real(long)  :: lam, vsum

   real(long)  :: error1, error2
   real(long),parameter :: stress_rescale=1.0d+2

   error_index = N_residual
   if(domain) error_index = error_index - N_cell_param

   M = N - 1 + ensemble 

   !  Initiallization
   call density(N, omega, rho, q, q_solvent)
   call mu_phi(mu_chain,phi_chain,q,mu_solvent,phi_solvent,q_solvent)
   call make_correlation(N)
   do i = 1, N_cell_param
      call make_dGsq(dGsq(:,i), dGG_basis(:,:,i))
   enddo
  
   ! set the homogeneous coefficents
   if (ensemble == 0) call set_omega_uniform(omega)

   ! calculate the maximum step size
   vsum = 0.0_long
   do alpha = 1, N_monomer
      vsum = vsum + dot_product(omega(alpha,2-ensemble:), &!
                                omega(alpha,2-ensemble:))
   end do
   if(domain) vsum = vsum + dot_product(cell_param(1:N_cell_param),&!
                                        cell_param(1:N_cell_param))
   stpmax = STPMX * max( sqrt(vsum), dble(N_residual) )

   itr = 0

   ! Calculate density/stress and residual/error
   if (domain) stress = scf_stress(N, N_cell_param, dGsq)
   call make_residual(N, domain, rho, stress, omega, residual)
   !error    = maxval(abs(residual(1:error_index)))
   !error     = maxval(abs(residual))
   error1 = maxval(abs(residual(1:error_index)))
   error2 = maxval(abs(residual(error_index:N_residual)))
   error  = max(error1, error2*stress_rescale)

   ! Calculate thermodynamic properties
   call free_energy(N,rho,omega,phi_chain,mu_chain,phi_solvent,mu_solvent,f_Helmholtz,pressure)
   !# ifdef DEVEL
   call divide_energy(rho,omega,phi_chain,phi_solvent,Q,f_Helmholtz,f_component,overlap)
   !# endif

   iterative_loop : do
    
      write(6,*)
      write(6,"('Iteration ',i3)") itr

      ! Output to log file
      call output(error,                     'error       =',o=6,f='L')
      if (domain) then
         call output(stress,N_cell_param,    'stress      =',o=6,f='L')
         call output(cell_param,N_cell_param,'cell_param  =',o=6,f='L')
      endif
      call output(f_Helmholtz,               'f_Helmholtz =',o=6,f='L')
      if (ensemble == 1) then
         call output(pressure,               'pressure    =',o=6,f='L')
      endif 

      if ( error < error_max ) then
         converge = .true.
         exit iterative_loop
      else if ( itr > max_itr ) then
         converge = .false.
         exit iterative_loop
      else if ( error > large ) then
         converge = .false.
         exit iterative_loop
      endif

      ! Calculate the inverse of Initial Jacobian/Broyden
      if ( Jacobian_empty ) then

         write(6,"('Initializing Jacobian ...')")
         call Jacobian_response(N, N_cut, domain, omega, Jacobian)
         Jacobian_empty = .false.      
        
         write(6,"('Jacobian ...')")
         write(6,"('Inverting Jacobian ...')")
         if( .not. full_inversion ) then
            call invert_Jacobian(N, N_cut,domain,Jacobian)
         else
            ! inverting the full Jacobian
            call dgetrf(N_residual,N_residual,Jacobian,N_residual,ipvt,info)
            if(info/=0) stop "Full Jacobian LU factorization failed."
            call dgetri(N_residual,Jacobian,N_residual,ipvt,work,lwork,info)
            if(info/=0) stop "Full Jacobian inversion failed."
         end if

      end if

      ! Update Jacobian^(-1) with Broyden's method
      ! To save memory, from now on, Jacobian = Jacobian^(-1)
      if ( itr > 0 ) then
         delta_norm = dot_product(delta, delta)
         delta = delta / delta_norm
         r_old = residual - (1.0_long - lam) * r_old
         u = matmul(Jacobian, r_old)
         do i = 1, N_residual
            v(i) = dot_product(delta, Jacobian(:,i))
         end do
         uv = 1.0_long + dot_product(v, r_old)
         u = u / uv
         do i = 1, N_residual
            Jacobian(:,i) = Jacobian(:,i) - v(i) * u
         end do
      endif
      r_old = residual
  
      ! Solve Newtonian equation: Jacobian * delta = - residual
      p     = - matmul(Jacobian,residual)
      fold  = dot_product(residual, residual) * 0.5_long
      slope = - fold * 2.0_long

      ! Update potential and density fields with/without backtracking
        call update_with_linesearch
      ! call update_without_linesearch

      ! Calculate thermodynamic properties
      call mu_phi(mu_chain,phi_chain,q,mu_solvent,phi_solvent,q_solvent)
      call free_energy(N,rho,omega,phi_chain,mu_chain,phi_solvent,mu_solvent,f_Helmholtz,pressure)
      !# ifdef DEVEL
      call divide_energy(rho,omega,phi_chain,phi_solvent,Q,f_Helmholtz,f_component,overlap)
      !# endif

      itr = itr + 1

   end do iterative_loop

   ! If fixed unit cell, calculate residual stress before output
   if (.not.domain) stress = scf_stress(N, N_cell_param, dGsq)

   contains
      !----------------------------------------------------------
      ! Accepting Full Newton/Broyden step unconditionally
      !----------------------------------------------------------
      subroutine update_without_linesearch()
      implicit none

      vsum = sqrt(dot_product(p,p))
      if (vsum > stpmax) then
         write(6,*) "residual rescaled"
         p = p * stpmax / vsum
      end if

      lam   = 1.0_long
      delta = lam * p

      ! Update omega
      k = 0
      do i = 1, M
         do alpha = 1, N_monomer
            k = k + 1
            omega(alpha,i+1-ensemble) = omega(alpha,i+1-ensemble) +  delta(k)
         enddo
      enddo
   
      ! Update unit cell 
      if (domain) then
         do i = 1, N_cell_param
            cell_param(i) = cell_param(i) +  delta(M*N_monomer + i)
         enddo
         call make_unit_cell
         call make_correlation(N)
         call make_ksq(G_basis)
         do i = 1, N_cell_param
            call make_dGsq(dGsq(:,i), dGG_basis(:,:,i))
         enddo
      endif

      ! Update density/stress and residuals/error
      call density(N, omega, rho, q, q_solvent)
      if (domain) stress = scf_stress(N, N_cell_param, dGsq)

      call make_residual(N, domain, rho, stress, omega, residual)
!     error = maxval(abs(residual))
      error1 = maxval(abs(residual(1:error_index)))
      error2 = maxval(abs(residual(error_index:N_residual)))
      error  = max(error1, error2*stress_rescale)

      end subroutine update_without_linesearch
      !----------------------------------------------------------


      !----------------------------------------------------------
      ! back tracking line search subroutine
      ! instead of accepting Newton/Broyden step straightforwardly
      ! we look for the optimized step along N/B direction
      !----------------------------------------------------------
      subroutine update_with_linesearch()
      implicit none
      real(long), parameter :: ALF=1.0D-4, TOLX=1.0D-7
      real(long)  :: lam2               ! (0,1], controlling step size
      real(long)  :: lamin, lamax       ! Min and Max step size
      real(long)  :: tmplam             !
      real(long)  :: test, temp         ! auxillary variables
      real(long)  :: rhs1, rhs2, f2     !
      real(long)  :: a, b, disc         !

      vsum = sqrt(dot_product(p,p))
      if (vsum > stpmax) then
         write(6,*) "residual rescaled"
         p = p * stpmax / vsum
         slope = slope * stpmax / vsum
      end if

      ! compute minimum lam
      test = 0.0_long
      k = 0
      do i = 1, M
         do alpha = 1, N_monomer
            k = k + 1
            temp = abs(p(k))/max(abs(omega(alpha,i+1-ensemble)), 1.0_long)
            if( temp > test ) test = temp
         enddo
      enddo
      if (domain) then
         do i = 1, N_cell_param
            k = k + 1
            temp = abs(p(k))/max(abs(cell_param(i)), 1.0_long)
            if( temp > test ) test = temp
         enddo
      endif
      lamin = TOLX / test

      lam   = 1.0_long
      ow    = omega
      ocell = cell_param(1:N_cell_param)
      backtrack_loop : do
         delta = lam * p

         ! Update omega
         k = 0
         do i = 1, M
            do alpha = 1, N_monomer
               k = k + 1
               omega(alpha,i+1-ensemble) = ow(alpha,i+1-ensemble) + delta(k)
            enddo
         enddo
   
         ! Update unit cell 
         if (domain) then
            do i = 1, N_cell_param
               cell_param(i) = ocell(i) + delta(M*N_monomer + i)
            enddo
            call make_unit_cell
            call make_correlation(N)
            call make_ksq(G_basis)
            do i = 1, N_cell_param
               call make_dGsq(dGsq(:,i), dGG_basis(:,:,i))
            enddo
         endif

         ! Calculate trial density/stress and residuals/error
         call density(N, omega, rho, q, q_solvent)
         if (domain) stress = scf_stress(N, N_cell_param, dGsq)
         call make_residual(N, domain, rho, stress, omega, residual)
         ! error = maxval(abs(residual(1:error_index)))
         error = maxval(abs(residual))
         error1 = maxval(abs(residual(1:error_index)))
         error2 = maxval(abs(residual(error_index:N_residual)))
         error  = max(error1, error2*stress_rescale)

         ! backtracking tests
         fnew  = dot_product(residual, residual) * 0.5_long
         if ( lam < lamin ) then
            write(6,"('Minimum lambda in line search reached.')")
            exit backtrack_loop
         else if( fnew < fold + ALF*lam*slope ) then 
            exit backtrack_loop
         else
            if ( lam == 1.0_long ) then
               tmplam = -slope/(2.0_long*(fnew-fold-slope))
            else
               rhs1 = fnew-fold-lam *slope
               rhs2 = f2  -fold-lam2*slope
               a = (rhs1/lam**2-rhs2/lam2**2)/(lam-lam2)
               b = (-lam2*rhs1/lam**2+lam*rhs2/lam2**2)/(lam-lam2)
               if (a==0.0_long) then
                  tmplam = -slope/(2.0_long*b)
               else
                  disc = b*b-3.0_long*a*slope
                  if (disc < 0.0_long) then
                     tmplam = 0.5_long*lam
                  else if (b < 0.0_long) then
                     tmplam = (-b+sqrt(disc))/(3.0_long*a)
                  else
                     tmplam = -slope/(b+sqrt(disc))
                  endif
               endif
               if (tmplam > 0.5_long*lam) tmplam=0.5_long*lam
            endif
         endif

         lam2 = lam
         f2   = fnew
         lam  = max(tmplam, 0.1_long * lam)
      end do backtrack_loop

      end subroutine update_with_linesearch
      !----------------------------------------------------------
   end subroutine iterate_NR
   !==============================================================

   !--------------------------------------------------------------
   !****ip iterate_mod/invert_Jacobian
   ! SUBROUTINE
   !    subroutine invert_Jacobian(N,N_cut,domain,Jacobian)
   ! PURPOSE
   !    Invert the long wave length block (J_corner) of Jacobian
   ! ARGUMENTS
   !    N             = # of basis function
   !    N_cut         = # of long wave length (cut-off dimension)
   !    domain        = false, no stress related entries are considered
   !                    true, stress related entries are combined with
   !                    long wave length entries of Jacobian
   !    Jacobian(:,:) = self explanary
   ! COMMENT
   !
   !            a |  | b              a | b |
   !           ---\  |---            ---+---|
   !               \ |        ===>    c | d |
   !           _____\|___            -------\
   !            c |  | d                     \
   !
   !           (Jacobian)         (block diagonal)
   !
   !    Use the sparseness of Jacobian to construct inversion.
   !      1) invert J_corner=(a,b;c,d)
   !      2) shift entries of inversion of J_corner back to Jacobian
   !      3) the diagonal blocks between a and d are then inverted revursively
   !    
   !    Note:
   !      The block diagonal Jacobian is not explicitly needed. The internal
   !      subroutine "mapping" is designed to shift a,b,c,d between Jacobian
   !      and J_corner.
   ! SOURCE
   !--------------------------------------------------------------
   subroutine invert_Jacobian(N,N_cut,domain,Jacobian)
   use unit_cell_mod, only     :  N_cell_param
   use chemistry_mod, only     :  ensemble, N_monomer
   implicit none

   integer,    intent(IN)     ::  N     ! # basis function
   integer,    intent(IN)     ::  N_cut ! basis coeff cutoff
   logical,    intent(IN)     ::  domain
   real(long), intent(INOUT)  ::  Jacobian(:,:)
   !***

   real(long)  ::  diagonal(N_monomer,N_monomer)
   integer     ::  info, i

   if( domain ) then
      call mapping(Jacobian,N_residual,J_corner,N_corner,N_cell_param)
   else
      call mapping(Jacobian,N_residual,J_corner,N_corner,0)
   end if

   ! LU factorization and inversion of J_corner
   call dgetrf(N_corner,N_corner,J_corner,N_corner,ipvt(1:N_corner),info)
   if ( info /=0 ) stop "Jacobian corner LU-factor failed."
   call dgetri(N_corner,J_corner,N_corner,ipvt(1:N_corner),work,lwork,info)
   if ( info /= 0) stop "Jacobian corner inversion failed."

   if( domain ) then
      call mapping(J_corner,N_corner,Jacobian,N_residual,N_cell_param)
   else
      call mapping(J_corner,N_corner,Jacobian,N_residual,0)
   end if
 
   ! recursively invert the other diagonal blocks of Jacobian
   do i = max(1,N_cut) + ensemble, N + ensemble -1
!  do i = N_cut + ensemble, N + ensemble -1
      diagonal = Jacobian( (i-1)*N_monomer+1: i*N_monomer,    &
                           (i-1)*N_monomer+1: i*N_monomer)
      ! LU factorization and inversion of diagonal blocks
      call dgetrf(N_monomer,N_monomer,diagonal,N_monomer,ipvt(1:N_monomer),info)
      if ( info /=0 ) stop "Jacobian diagonal LU-factor failed."
      call dgetri(N_monomer,diagonal,N_monomer,ipvt(1:N_monomer),work,lwork,info)
      if ( info /= 0) stop "Jacobian diagonal inversion failed."

      Jacobian( (i-1)*N_monomer+1: i*N_monomer,             &
                (i-1)*N_monomer+1: i*N_monomer) = diagonal

   end do

   contains
      !---------------------------------------------------------
      subroutine mapping(A,na,B,nb,nc)
      implicit none
      real(long)  ::  A(:,:), B(:,:)
      integer     ::  na, nb, nc, ntmp

      ntmp = min(na,nb) - nc
      if( ntmp < 0 ) stop "Error while mapping Jacobian."
      B(1:ntmp,1:ntmp) = A(1:ntmp,1:ntmp)

      if ( nc > 0 ) then
         B(nb-nc+1:nb,1:ntmp)     = A(na-nc+1:na,1:ntmp)
         B(1:ntmp,nb-nc+1:nb)     = A(1:ntmp,na-nc+1:na) 
         B(nb-nc+1:nb,nb-nc+1:nb) = A(na-nc+1:na,na-nc+1:na) 
      end if

      if ( nb > na ) then
         B(ntmp+1:nb-nc,1:ntmp)     = 0.0_long
         B(1:ntmp,ntmp+1:nb-nc)     = 0.0_long
         B(nb-nc+1:nb,ntmp+1:nb-nc) = 0.0_long
         B(ntmp+1:nb-nc,nb-nc+1:nb) = 0.0_long
      end if

      end subroutine mapping
      !---------------------------------------------------------

   end subroutine invert_Jacobian
   !==============================================================

   !---------------------------------------------------------------------
   !****p iterate_mod/make_residual
   ! SUBROUTINE
   !    make_residual( N, domain, rho, stress, omega, residual )
   ! PURPOSE
   !    Calculate array residual(1:N_residual) 
   !    Uses input values of rho, omega, and stress arrays
   ! ARGUMENTS
   !   integer N              = # of basis functions
   !   logical domain         = T -> domain iteration, F -> fixed cell
   !   real(long) rho(:,:)    = monomer density field
   !   real(long) stress(:)   = stress components
   !   real(long) omega(:,:)  = chemical potential field
   !   real(long) residual(:) = residual array, dimension(1:N_residual)
   ! SOURCE
   !---------------------------------------------------------------------
   subroutine make_residual( N, domain, rho, stress, omega, residual )
   use unit_cell_mod, only: N_cell_param
   use chemistry_mod, only: ensemble, N_monomer, chi
   implicit none

   integer,    intent(IN) :: N
   logical,    intent(IN) :: domain
   real(long), intent(IN) :: rho(:,:)
   real(long), intent(IN) :: stress(:)
   real(long), intent(IN) :: omega(:,:)
   real(long), intent(OUT):: residual(:)
   !***

   ! Local variables
   integer    :: i, alpha, l
   real(long) :: xi(N)    ! Lagrange pressure field

   ! Lagrange field 
   xi = omega(1,:)
   do i = 1, N
      xi(i) = xi(i) - sum( chi(:,1)*rho(:,i) )
   enddo
 
   residual = 0.0_long
   ! Incompressiblity and omega residuals
   l = 0
   do i = 2-ensemble, N
      l = l + 1
      residual(l) = residual(l) + sum( rho(:,i) )

      do alpha = 2, N_monomer
         l = l + 1
         residual(l) = omega(alpha,i) - xi(i)
         residual(l) = residual(l) - sum( chi(:,alpha)*rho(:,i) )
      enddo
   enddo
   if (ensemble == 1) residual(1) = residual(1) - 1.0_long
  
   ! Stress residuals
   if (domain) residual(l+1:l+N_cell_param) = stress(:)

   end subroutine make_residual
   !====================================================================


   !--------------------------------------------------------------------
   !****p iterate_mod/Jacobian_response
   ! SUBROUTINE
   !    Jacobian_response(N, N_cut, domain, omega, Jacobian)
   ! PURPOSE
   !    Construct Jacobian from response functions
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine Jacobian_response(N, N_cut, domain, omega, Jacobian)
   use chemistry_mod, only : ensemble, N_monomer, chi
   use unit_cell_mod, only : N_cell_param
   implicit none

   
   integer, intent(IN)       :: N
   integer, intent(IN)       :: N_cut
   logical, intent(IN)       :: domain
   real(long), intent(IN)    :: omega(:,:)
   real(long), intent(OUT)   :: Jacobian(:,:)
   !***

   
   ! Local Variables
   real(long)  :: drho_domega(:,:,:,:) 
   real(long)  :: dxi_domega(:,:,:) 
   allocatable :: drho_domega, dxi_domega
   real(long)  :: drho_dcell(N_monomer, N, N_cell_param)
   real(long)  :: dstress_dcell(N_cell_param,N_cell_param)
   real(long)  :: dxi_dcell(N,N_cell_param) 

   
   ! indices and temporary variables needed ...
   integer    :: j, beta, l   ! omega field indices
   integer    :: i, alpha, k  ! residual indices
   !integer   :: n_cut        ! cut + 1, index
   integer    :: imin, imax   ! output indices
   integer    :: info   ! output indices


   allocate(drho_domega(N_monomer,N,N_monomer,N),stat=info)
   if( info /= 0) stop "drho_domega memory allocation failed in Jacobian_Response"

   allocate(dxi_domega(N,N_monomer,N),stat=info)
   if( info /= 0) stop "dxi_domega memory allocation failed in Jacobian_Response"

   ! drho_domega (density response)
   call response_pd_omega(N,N_cut,omega,drho_domega)
   
   ! dxi_domega (Lagrange response)
   do i = 1, N
      do alpha = 1, N_monomer
         do j = 1, N
            dxi_domega(j,alpha,i) = - sum( chi(:,1)*drho_domega(:,j,alpha,i) )
         enddo
      enddo
      dxi_domega(i,1,i) = dxi_domega(i,1,i) + 1.0_long
   enddo 

   ! Updating Jacobian, outmost loop is omega
   Jacobian = 0.0_long
   do j = 2 - ensemble, N
      do beta = 1, N_monomer
         l = beta + N_monomer*(j-2+ensemble)
         ! incompressibility and scf omega response
         do i = 2 - ensemble, N
            k = 1 + N_monomer*(i-2+ensemble)
            Jacobian(k,l) = sum( drho_domega(:,i,beta,j) )
            !if ( k == l) then
            !   print *, k
            !   print *, Jacobian(k,k)
            !end if  


            do alpha = 2, N_monomer
               k = alpha + N_monomer*(i-2+ensemble)
               Jacobian(k,l)= - sum( chi(:,alpha)*drho_domega(:,i,beta,j) )&!
                              - dxi_domega(i,beta,j)
            enddo
         enddo
         if (beta > 1) Jacobian(l,l) = Jacobian(l,l) + 1.0_long
         !Jacobian(l,l) = Jacobian(l,l) + 1.0_long
         !print *, beta,l
         !print *, Jacobian(l,l)
         
      enddo
   enddo

   if ( domain ) then

      call response_pd_cell(N,omega,drho_dcell,dstress_dcell )

      ! dxi_dcell (Lagrangian response)
      do j = 1, N_cell_param
         do i = 1, N
            dxi_dcell(i,j) = - sum ( chi(:,1)*drho_dcell(:,i,j) )
         enddo
      enddo

      ! Stress Response to omega, by transpostion of drho_dcell
      do i = 1, N_cell_param
         k = i + N_monomer*(N-1+ensemble)
         do j = 2 - ensemble, N
            do alpha = 1, N_monomer
               l = alpha + N_monomer*(j-2+ensemble) 
               !Jacobian(k,l) = chain_length * drho_dcell(alpha,j,i)
               Jacobian(k,l) = drho_dcell(alpha,j,i)
            enddo
         enddo
      enddo

      do j = 1, N_cell_param
         l = j + N_monomer * (N-1+ensemble)

         do i = 2 - ensemble, N
            ! Incompressibilty Response
            k = 1 + N_monomer * (i-2+ensemble)
            Jacobian(k,l) = sum( drho_dcell(:,i,j) )

            ! Omega Response
            do alpha = 2, N_monomer
               k = alpha + N_monomer*(i-2+ensemble)
               Jacobian(k,l) = - sum( chi(:,alpha)*drho_dcell(:,i,j) )&!
                               - dxi_dcell(i,j)
            enddo
         enddo

         ! Stress Response
         Jacobian(k+1:k+N_cell_param,l) = dstress_dcell(:,j)
      enddo

   end if

   if( allocated(drho_domega) ) deallocate(drho_domega)
   if( allocated(dxi_domega) ) deallocate(dxi_domega)

   end subroutine Jacobian_response
   !===================================================================

end module iterate_mod
