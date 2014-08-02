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
!****m scf/sweep_mod
! MODULE
!    sweep_mod
! PURPOSE
!    Solve SCF repeatedly for a range of parameters, using
!    zeroth, first, or second order continuation of solutions
! AUTHOR
!    David Morse (1/2004)
!    based on earlier version for multiblock melts by Chris Tyler
! COMMENT
!    Module uses 1st or 2nd order contination of scf equations
!    at equally spaced values of a continuation parameter s along 
!    a line in parameter space along which each of the chemistry
!    input parameters kuhn, block_length, chi or Temperature 
!    (depending upon the value of chi_fmt), and phi or mu 
!    (depending upon the ensemble) has a constant increment per 
!    unit change in s. For solutions of a specified sequence of 
!    unit cell parameters, with domain = F, (rather than iteration 
!    of each set of chemical parameters to stress free solution, 
!    with domain = T), the elements of cell_param also have a 
!    constant increment per unit change in s.  A single parameter 
!    may be varied by setting all but one of the increments to zero. 
!     
!    The required increments are read from file by subroutine
!    input_increments and are stored in private module variables
!    d_kuhn, d_block_length, d_chi or d_Temperature, etc.
!   
!    The order of continuation is set at compile time by the value
!    of the private module parameter N_hist, which gives the maximum
!    number of solutions stored in a history stack. Compile with
!    N_hist = 3 (current solution and two previous solutions) for
!    second order continuation, and N_hist = 2 (current and one
!    previous solution) for first order continuation.
!
! SOURCE
!-------------------------------------------------------------------
module sweep_mod
   use const_mod
   use unit_cell_mod, only : N_cell_param
   use chemistry_mod
   implicit none

   private

   ! Public procedures
   public :: input_increments     ! Read nonzero increments from file
   public :: output_increments    ! Write increments to file
   public :: increment_parameters ! Increment parameters
   public :: history_setup        ! Allocate history stack arrays
   public :: update_history       ! Push new solution onto history stack
   public :: continuation         ! Predict next solution
   !***

   ! Private module variables

   ! Increment variables 
   real(long)  :: d_kuhn(:)           ! d(Kuhn length) of (monomer)
   real(long)  :: d_chi(:,:)          ! d(Flory chi) of (monomer,monomer)
   real(long)  :: d_temperature       ! d(absolute temperature)
   real(long)  :: d_block_length(:,:) ! d(# monomers) in (block,chain)
   real(long)  :: d_phi_chain(:)      ! d(volume fraction) (chain)
   real(long)  :: d_mu_chain(:)       ! d(chemical potential) (chain)
   real(long)  :: d_phi_solvent(:)    ! d(volume fraction) (solvent)
   real(long)  :: d_mu_solvent(:)     ! d(chemical potential) (solvent)
   real(long)  :: d_cell_param(6)     ! d(unit cell parameters)
   real(long)  :: d_solvent_size(:)   ! d(solvent size) (solvent)   

   ALLOCATABLE :: d_kuhn, d_chi, d_block_length, d_phi_chain,d_phi_solvent, d_mu_chain, d_mu_solvent, d_solvent_size

   ! Logical flags for nonzero increments
   logical :: sweep_kuhn          ! true if d_kuhn /= 0
   logical :: sweep_chi           ! true if d_chi or d_temperature /= 0
   logical :: sweep_block_length  ! true if d_block_length /=0
   logical :: sweep_composition   ! true if d_phi or d_mu /= 0
   logical :: sweep_cell_param    ! true if d_cell_param /=0
   logical :: sweep_solvent_size  ! true if d_solvent_size /=0

   ! History stack variables 
   integer     :: depth_hist
   real(long)  :: s_hist(:)
   real(long)  :: omega_hist(:,:,:)
   real(long)  :: cell_param_hist(:,:)

   ALLOCATABLE :: s_hist, omega_hist, cell_param_hist

   ! History stack depth parameter
   integer, parameter :: N_hist = 3

contains


   !------------------------------------------------------------------
   !****p sweep_mod/input_increments
   ! SUBROUTINE
   !    input_increments(i_unit,fmt_flag,domain)
   ! PURPOSE
   !    Allocate arrays for chemistry increment variables
   !    Read values of increment module variables from unit i_unit
   ! NOTE
   !    Default imput unit in module io_mod is set to i_unit
   ! SOURCE
   !------------------------------------------------------------------
   subroutine input_increments(i_unit,fmt_flag,domain)
   use io_mod
   integer,      intent(IN) :: i_unit        ! input unit #
   character(1), intent(IN) :: fmt_flag      ! F = formatted, U = unformatted
   logical,      intent(IN) :: domain        ! true if domain iteration
   !***

   integer                 :: i, j          ! loop indices
   character(len = 100)    :: comment_line  

   ! Allocate arrays and zero all increments
   allocate(d_kuhn(N_monomer))
   d_kuhn = 0.0_long
   select case(interaction_type)
   case('chi')
      allocate(d_chi(N_monomer,N_monomer))
      d_chi = 0.0_long
   case('chi_T')
      d_temperature = 0.0_long
   end select
   allocate(d_block_length(N_blk_max,N_chain))
   d_block_length = 0.0_long
   if (ensemble == 0) then
      allocate(d_phi_chain(N_chain))
      allocate(d_phi_solvent(N_solvent))
      d_phi_chain  = 0.0_long
      d_phi_solvent= 0.0_long
   else 
      allocate(d_mu_chain(N_chain))
      d_mu_chain = 0.0_long
      allocate(d_mu_solvent(N_solvent))
      d_phi_solvent=0.0_long 
   endif

   allocate(d_solvent_size(N_solvent))
   d_solvent_size = 0.0_long  

   select case(fmt_flag)
   case('F','f')

      call set_io_units(i=i_unit)
      call input(d_kuhn,N_monomer,'d_kuhn',s='R')
      select case(interaction_type)
      case('chi')
         call input(d_chi,N_monomer,N_monomer,'d_chi',s='L')
      case('chi_T')
         call input(d_temperature,'d_temperature')
      end select
      read(i_unit,*)comment_line
      call output(trim(comment_line),f='N',j='L')
      do j = 1, N_chain
         call input(d_block_length(:,j),N_block(j),f='N')
      enddo
      select case (ensemble)
      case (0)
         call input(d_phi_chain,N_chain,'d_phi_chain',s='C')
         call input(d_phi_solvent,N_solvent,'d_phi_solvent',s='C')
      case (1)  ! Grand canonical ensemble 
         call input(d_mu_chain, N_chain,'d_mu_chain', s='C')
         call input(d_mu_solvent, N_solvent,'d_mu_solvent',s='C')
      end select
      call input(d_solvent_size,N_solvent,'d_solvent_size',s='C')
      if (.not.domain) then
         call input(d_cell_param, N_cell_param, 'd_cell_param', s='R')
         sweep_cell_param   = .true.
      else
         sweep_cell_param   = .false.
      endif

      sweep_kuhn         = .true.
      sweep_chi          = .true.
      sweep_block_length = .true.
      sweep_composition  = .true.
      sweep_solvent_size = .true.

   case('N','n')

      sweep_kuhn         = .false.
      sweep_chi          = .false.
      sweep_block_length = .false.
      sweep_composition  = .false.
      sweep_cell_param   = .false.
      sweep_solvent_size = .false.

      call set_io_units(i=i_unit)
      do
         
         ! Read comment line
         read(i_unit,*) comment_line
         call output(trim(comment_line),f='N',j='L')

         select case(trim(comment_line))
         case('d_kuhn')
            call input(d_kuhn,N_monomer,'d_kuhn',s='R',f='N')
            sweep_kuhn = .true.
         case('d_chi')
            if (interaction_type=='chi') then
               call input(d_chi,N_monomer,N_monomer,'d_chi',s='L',f='N')
               sweep_chi = .true.
            else
               write(6,*) &
                     'Error: Increment chi when interaction_type /= chi'
               stop 
            endif 
         case('d_temperature')
            if (interaction_type=='chi_T') then
               call input(d_temperature,'d_temperature',f='N')
               sweep_chi = .true.
            else
               write(6,*) &
               'Error: Increment T when interaction_type /= chi_T'
               stop ! Error handling
            endif 
         case('d_block_length')
            do j = 1, N_chain
               call input(d_block_length(:,j),N_block(j),f='N')
            enddo
            sweep_block_length = .true.
         case('d_phi_chain')
            if (ensemble == 0) then 
               call input(d_phi_chain,N_chain,'d_phi_chain',s='C',f='N')
               sweep_composition = .true.
            else
               write(6,*) 'Error: Increment phi_chain when ensemble /= 0'
               stop
            endif 
         case('d_phi_solvent')
            if (ensemble == 0) then
               call input(d_phi_solvent,N_solvent,'d_phi_solvent',s='C',f='N')
               sweep_composition = .true.
            else
               write(6,*) 'Error: Increment phi_solvent when ensemble /=0'
               stop
            endif
         case('d_mu_chain')
            if (ensemble == 1) then 
               call input(d_mu_chain,N_chain,'d_mu_chain',s='C',f='N')
               sweep_composition = .true.
            else
               write(6,*) 'Error: Increment mu_chain when ensemble neq 1'
               stop 
            endif 
         case('d_mu_solvent')
            if (ensemble == 1) then
               call input(d_mu_solvent,N_solvent,'d_mu_solvent',s='C',f='N')
               sweep_composition = .true.
            else
               write(6,*) 'Error: Increment mu_solvent when ensemble neq 1'
               stop
            endif
         case('d_solvent_size')
            call input(d_solvent_size,N_solvent,'d_solvent_size',s='C',f='N')
            sweep_solvent_size = .true.
         case('d_cell_param')
            if (.not.domain) then
               call input &
                  (d_cell_param,N_cell_param,'d_cell_param',s='R',f='N')
               sweep_cell_param = .true.
            else
               write(6,*) 'Error: Increment cell_param when domain == T'
               stop 
            endif
         case('end_increments')
            exit
         case default
            write(6,*) 'Error: Invalid comment line in input_increments'
            stop
         end select

      enddo

   case default

      write(6,*) 'Error: Illegal fmt_flag in input_increments'
      stop

   end select

   end subroutine input_increments
   !=================================================================


   !------------------------------------------------------------------
   !****p sweep_mod/output_increments
   ! SUBROUTINE
   !    output_increments(o_unit,fmt_flag,domain)
   ! PURPOSE
   !    Write values of increments to unit o_unit 
   !    in format required by input_increments
   ! NOTE
   !    Default output unit in module io_mod set to o_unit
   ! SOURCE
   !------------------------------------------------------------------
   subroutine output_increments(o_unit,fmt_flag,domain)
   use io_mod
   integer,      intent(IN) :: o_unit        ! input unit #
   character(1), intent(IN) :: fmt_flag      ! F = formatted, U = unformatted
   logical,      intent(IN) :: domain        ! true if domain iteration
   !***

   integer                 :: i, j          ! loop indices


   ! Set default output unit to o_unit
   call set_io_units(o=o_unit)

   select case(fmt_flag)
   case('F','f') ! Old format - all increments

      call output(d_kuhn,N_monomer,'d_kuhn',s='R')
      select case(interaction_type)
      case('chi')
         call output(d_chi,N_monomer,N_monomer,'d_chi',s='L')
      case('chi_T')
         call output(d_temperature,'d_temperature')
      end select
      call output('d_block_length',f='N',j='L')
      do j = 1, N_chain
         call output(d_block_length(:,j),N_block(j),f='N')
      enddo
      select case (ensemble)
      case (0)
         call output(d_phi_chain,N_chain,'d_phi_chain',s='C')
         call output(d_phi_solvent,N_solvent,'d_phi_solvent',s='C') 
      case (1)  ! Grand canonical ensemble 
         call output(d_mu_chain,N_chain,'d_mu_chain', s='C')
         call output(d_mu_solvent,N_solvent,'d_mu_solvent',s='C')
      end select
      call output(d_solvent_size,N_solvent,'d_solvent_size',s='C')
      if (sweep_cell_param) then
         call output(d_cell_param, N_cell_param, 'd_cell_param', s='R')
      endif

   case('N','n') ! New format - selected increments

      if (sweep_kuhn) then
         call output(d_kuhn,N_monomer,'d_kuhn',s='R')
      endif
      if (sweep_chi) then
         select case(interaction_type)
         case('chi')
            call output(d_chi,N_monomer,N_monomer,'d_chi',s='L')
         case('chi_T')
            call output(d_temperature,'d_temperature')
         end select
      endif
      if (sweep_block_length) then
         call output('d_block_length',f='N',j='L')
         do j = 1, N_chain
            call output(d_block_length(:,j),N_block(j),f='N')
         enddo
      endif
      if (sweep_composition) then
         select case (ensemble)
         case (0)
            call output(d_phi_chain,N_chain,'d_phi_chain',s='C')
            call output(d_phi_solvent,N_solvent,'d_phi_solvent',s='C')
         case (1)  ! Grand canonical ensemble 
            call output(d_mu_chain, N_chain,'d_mu_chain', s='C')
            call output(d_mu_solvent, N_solvent,'d_mu_solvent', s='C')
         end select
      endif
      if (sweep_solvent_size) then
         call output(d_solvent_size,N_solvent,'d_solvent_size',s='C')
      endif
      if (sweep_cell_param) then
         call output(d_cell_param, N_cell_param, 'd_cell_param', s='R')
      endif
      call output('end_increments',f='N',j='L')

   case default

      write(6,*) 'Error: Illegal fmt_flag in output_increments'
      stop

   end select

   end subroutine output_increments
   !=================================================================


   !-----------------------------------------------------------------
   !****p sweep_mod/increment_parameters
   !SUBROUTINE
   !   increment_parameters(step,domain,cell_param)
   ! PURPOSE
   !   Increment chemisty parameters 
   !   if (.not.domain), also increment cell_param
   ! ARGUMENTS
   !   step       - change in continuation parameter s
   !   domain     - F -> specified unit cell, T -> variable unit cell
   !   cell_param - cell parameters
   ! SOURCE
   !-----------------------------------------------------------------
   subroutine increment_parameters(step,domain,cell_param)
   real(long) :: step
   logical    :: domain
   real(long) :: cell_param(:)
   !***

   integer    :: i, j
   kuhn = kuhn + step*d_kuhn
   select case(interaction_type)
   case('chi')
      chi = chi + step*d_chi
   case('chi_T')
      temperature = temperature + step*d_temperature
      chi = chi_A/temperature + chi_B
   end select
   do j=1, N_chain
      block_length(1:N_block(j),j) = block_length(1:N_block(j),j) &
                       + step*d_block_length(1:N_block(j),j)
   enddo
   do j = 1, N_chain
      chain_length(j) = 0.0_long
      do i=1, N_block(j)
         chain_length(j)  = chain_length(j) + block_length(i,j)
      enddo
   enddo
   select case (ensemble)
   case (0) 
      if( N_chain > 0 ) phi_chain   = phi_chain   + step*d_phi_chain
      if( N_solvent > 0 ) phi_solvent = phi_solvent + step*d_phi_solvent
   case (1)  
      if( N_chain > 0 ) mu_chain   = mu_chain   + step*d_mu_chain
      if( N_solvent > 0 ) mu_solvent = mu_solvent + step*d_mu_solvent
   end select
   if( N_solvent > 0 ) solvent_size = solvent_size + step*d_solvent_size
   if (.not.domain) then
      cell_param = cell_param + step*d_cell_param
   endif

   end subroutine increment_parameters
   !=================================================================


   !-----------------------------------------------------------------
   !****p sweep_mod/history_setup
   ! SUBROUTINE
   !    history_setup
   ! PURPOSE
   !    Allocates private history stack module variables:
   !       s_hist(N_hist)                      = prior values of s
   !       omega_hist(N_monomer,N_star,N_hist) = prior omega fields
   !       cell_param_hist(6,N_hist)           = prior cell parameters
   ! SOURCE
   !-----------------------------------------------------------------
   subroutine history_setup
   use basis_mod, only : N_star
   !***
   allocate(s_hist(N_hist))
   allocate(omega_hist(N_monomer,N_star,N_hist))
   allocate(cell_param_hist(6,N_hist))
   depth_hist      = 0
   omega_hist      = 0.0_long
   cell_param_hist = 0.0_long
   end subroutine history_setup
   !=================================================================



   !-----------------------------------------------------------------
   !****p sweep_mod/update_history
   ! SUBROUTINE
   !    update_history(s,omega,cell_param,domain)
   ! PURPOSE
   !    Push new values of s, omega, & cell_param onto history stacks.
   !    On output:
   !       s_hist(1)            = s
   !       omega_hist(:,:,1)    = omega
   !       cell_param_hist(:,1) = cell_param
   !    The cell_param stack is updated only if domain = T
   ! SOURCE
   !-----------------------------------------------------------------
   subroutine update_history(s,omega,cell_param,domain)
   use basis_mod, only : N_star
   real(long), intent(IN) :: s
   real(long), intent(IN) :: omega(:,:)
   real(long), intent(IN) :: cell_param(:)
   logical,    intent(IN) :: domain
   !***

   ! Local 
   integer :: i

   ! Calculate new depth for stack
   depth_hist = min(depth_hist+1,N_hist) 

   ! Shift current stack 
   if (depth_hist > 1) then
      do i = depth_hist - 1, 1, -1
         s_hist(i+1)            = s_hist(i)
         omega_hist(:,:,i+1)    = omega_hist(:,:,i)
         if (domain) then
            cell_param_hist(1:N_cell_param,i+1) &
               = cell_param_hist(1:N_cell_param,i)
         endif
      enddo
   endif

   ! Add input values to top of stack
   s_hist(1) = s 
   omega_hist(:,:,1) = omega
   if (domain) then
      cell_param_hist(1:N_cell_param,1) = cell_param(1:N_cell_param)
   endif

   end subroutine update_history
   !=================================================================


   !------------------------------------------------------------------
   !****p sweep_mod/continuation
   ! SUBROUTINE
   !    continuation(step,domain,omega,cell_param)
   ! PURPOSE
   !    Predict new values of omega and (if domain = T) cell_param
   !    for a change in s of size step. Uses the history stack of
   !    previous values of omega, cell parameters, and s.
   !
   !    Order of continuation depends on the value of the private
   !    variable depth_hist = current depth of stack:
   !       If depth_hist = 1 -> zeroth order continuation
   !       If depth_hist = 2 -> first order continuation
   !       If depth_hist > 2 -> second order continuation
   ! SOURCE
   !------------------------------------------------------------------
   subroutine continuation(step,domain,omega,cell_param)
   real(long), intent(IN)    :: step         ! change in s
   logical,    intent(IN)    :: domain       ! adjust cell_param
   real(long), intent(INOUT) :: omega(:,:)   ! chemical potential field
   real(long), intent(INOUT) :: cell_param(:)! unit cell parameters
   !***

   ! Local variables
   real(long) :: ds2, ds3, A2, A3, B2, B3, C1, C2, C3

   if (depth_hist > 2) then ! 2nd order continuation
      
      ds2 = s_hist(2) - s_hist(1)
      ds3 = s_hist(3) - s_hist(1)
      A2  = step*(ds3/ds2)/(ds3-ds2)
      A3  = step*(ds2/ds3)/(ds2-ds3)
      B2  = (step**2)/(ds2*(ds2-ds3))
      B3  = (step**2)/(ds3*(ds3-ds2))
      C2  = A2 + B2
      C3  = A3 + B3
      C1  = 1.0_long - C2 - C3
      omega = C1*omega_hist(:,:,1) &
            + C2*omega_hist(:,:,2) + C3*omega_hist(:,:,3)
      if (domain) then
         cell_param(1:N_cell_param) = &
              C1*cell_param_hist(1:N_cell_param,1) &
            + C2*cell_param_hist(1:N_cell_param,2) &
            + C3*cell_param_hist(1:N_cell_param,3) 
      endif

   else if (depth_hist == 2) then ! 1st order continuation

      ds2 = s_hist(2) - s_hist(1)
      C2  = step/ds2
      C1  = 1.0_long - C2
      omega = C1*omega_hist(:,:,1) + C2*omega_hist(:,:,2)
      if (domain) then
         cell_param(1:N_cell_param) = &
               C1*cell_param_hist(1:N_cell_param,1) &
             + C2*cell_param_hist(1:N_cell_param,2)
      endif

   else if (depth_hist == 1) then ! Zeroth order continuation

      omega = omega_hist(:,:,1)
      if (domain) then
         cell_param(1:N_cell_param) = &
               cell_param_hist(1:N_cell_param,1) 
      endif

   else

      write(6,*) 'Error in continuation'
      stop

   endif
   end subroutine continuation      
   !==================================================================

   end module sweep_mod
