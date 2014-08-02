!-----------------------------------------------------------------------
! PSCF - Polymer Self-Consistent Field Theory
! Copyright (2007) David C. Morse
! email: morse@cems.umn.edu
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation. A copy of this license is included in
! the LICENSE file in the top-level PSCF directory. 
!--------------------------------------------------------------------
!****m* scf/chemistry_mod
! PURPOSE
!    Defines module variables that describe the chemistry of a polymer 
!    mixture, including information about monomers, interactions, chain 
!    structure, solvent molecules, and composition of the mixture.
!
!    Subroutines to input and output chemistry data
! AUTHOR
!    Chris Tyler           (wrote multiblock melt scf module)
!    Amit Ranjan           (generalized to polymer mixtures)
!    David Morse           (separated chemistry_mod from scf_mod)
!    Raghuram Thiagarajan  (added small molecule solvent)
! SOURCE
!--------------------------------------------------------------------
   module chemistry_mod
   use const_mod
   implicit none

   private

   ! Public procedures
   public :: input_chemistry   ! read chemistry data from file
   public :: input_monomers, input_chains, input_solvents
   public :: input_composition, input_interaction
   public :: output_chemistry  ! write chemistry data to file
   public :: output_monomers, output_chains, output_solvents
   public :: output_composition, output_interaction
   public :: rescale_vref      ! rescale monomer reference volume

   ! Public variables
   public :: N_monomer, kuhn
   public :: N_chain, N_block, N_blk_max 
   public :: block_length, block_monomer, chain_length
   public :: N_solvent, solvent_monomer, solvent_size
   public :: ensemble, phi_chain, phi_solvent, mu_chain, mu_solvent
   public :: interaction_type, chi, chi_A, chi_B, temperature

   ! Monomer properties
   integer    :: N_monomer = 0     ! # of monomer types
   real(long) :: kuhn(:)           ! (N_monomer) Kuhn lengths

   ! Polymer properties
   integer    :: N_chain = 0        ! # of species types of polymer
   integer    :: N_blk_max          ! maximum # of blocks in any species
   integer    :: N_block(:)         ! (N_chain) # of blocks in species
   integer    :: block_monomer(:,:) ! (N_blk_max,N_chain) monomer type 
   real(long) :: block_length(:,:)  ! (N_blk_max,N_chain) # monomers 
   real(long) :: chain_length(:)    ! (N_chain) # monomers in polymer
   
   ! Solvent properties
   integer    :: N_solvent = 0      ! # of solvent types
   integer    :: solvent_monomer(:) ! (N_solvent) indices of solvents
   real(long) :: solvent_size(:)    ! (N_solvent) solvent volume / vref

   ! Statistical ensemble
   integer    :: ensemble           ! 0 = canonical, 1 = grand canonical

   ! Mixture composition
   real(long) :: phi_chain(:)   ! (N_chain) volume fractions of chains
   real(long) :: mu_chain(:)    ! (N_chain) chemical potentials of chains
   real(long) :: phi_solvent(:) ! (N_solvent) volume fractions of solvents
   real(long) :: mu_solvent(:)  ! (N_solvent) chemical potentials of solvents
  
   ! Interaction parameters
   real(long)   :: chi(:,:)         ! (N_monomer,N_monomer) Flory chi 
   character(20):: interaction_type ! 'chi' -> bare chi, 'chi_T' -> chi(T)
   real(long)   :: chi_A(:,:)       ! chi(T) = chi_A/T + chi_B
   real(long)   :: chi_B(:,:)
   real(long)   :: temperature      ! absolute temperature
   !***

   allocatable :: kuhn, chi, chi_A, chi_B, solvent_monomer, solvent_size
   allocatable :: N_block, block_monomer, block_length, chain_length
   allocatable :: phi_chain, mu_chain, phi_solvent, mu_solvent

!----------------------------------------------------------------------
!****v* chemistry_mod/N_monomer
! VARIABLE
!    integer      N_monomer     = # of monomer types
!*** ------------------------------------------------------------------
!****v* chemistry_mod/kuhn
! VARIABLE
!    real(long)   kuhn(:)       -  dimension(N_monomer) 
!                 kuhn(i)       =  Kuhn length for monomer i
!*** ------------------------------------------------------------------
!****v* chemistry_mod/chi
! VARIABLE
!    real(long)   chi(:,:)      - dimension(N_monomer,N_monomer) 
!                 chi(i,j)      =  Flory chi for monomers i and j
!*** ------------------------------------------------------------------
!****v* chemistry_mod/N_chain
! VARIABLE
!    integer      N_chain        = # of species types of polymers
!*** ------------------------------------------------------------------
!****v* chemistry_mod/N_blk_max
! VARIABLE
!    integer      N_blk_max     = maximum # of blocks in any species
!                               = max(N_block)
!*** ------------------------------------------------------------------
!****v* chemistry_mod/N_block
! VARIABLE
!    integer      N_block(:)    = (N_chain) # of blocks in species
!                 N_block(i)    = # of blocks in species i
!*** ------------------------------------------------------------------
!****v* chemistry_mod/block_monomer
! VARIABLE
!    integer block_monomer(:,:) - dimension(N_blk_max,N_chain) 
!            block_monomer(i,j) - monomer type of block i, species j
!*** ------------------------------------------------------------------
!****v* chemistry_mod/block_length
! VARIABLE
!   real(long) block_length(:,:) - dimension(N_blk_max,N_chain) 
!              block_length(i,j) = # monomers in block i, species j
!*** ------------------------------------------------------------------
!****v* chemistry_mod/chain_length
! VARIABLE
!    real(long) chain_length(:) -  dimension(N_chain)
!               chain_length(i) =  total # monomers in species i
!                               =  sum(block_length(1:N_block(i))
!*** ------------------------------------------------------------------
!****v* chemistry_mod/N_solvent
! VARIABLE
!    integer      N_solvent     = # of monomer types
!*** ------------------------------------------------------------------
!****v* chemistry_mod/solvent_monomer
! VARIABLE
!    integer      solvent_monomer    - dimension(N_solvent)
!                 solvent_monomer(i) = type of solvent monomer i (string)
!*** ------------------------------------------------------------------
!****v* chemistry_mod/solvent_size
! VARIABLE
!    real(long)   solvent_size    - dimension(N_solvent)
!                 solvent_size(i) = size of the solvent type i
!*** ------------------------------------------------------------------
!****v* chemistry_mod/ensemble
! VARIABLE
!    integer      ensemble   = index for statistical ensemble
!                              0 = canonical, 1 = grand canonical
!*** ------------------------------------------------------------------
!****v* chemistry_mod/phi_chain
! VARIABLE
!    real(long)   phi_chain(:)     - dimension(N_chain) 
!                 phi_chain(i)     = volume fraction of chain molecular species i
!*** ------------------------------------------------------------------
!****v* chemistry_mod/phi_solvent
! VARIABLE
!    real(long)   phi_solvent(:)    - dimension(N_solvent)
!                 phi_solvent(i)    = volume fraction of solvent molecular species i      
!***-------------------------------------------------------------------
!****v* chemistry_mod/mu_chain
! VARIABLE
!    real(long)   mu_chain(:)      - dimension(N_chain)
!                 mu_chain(i)      = chemical potential of chain molecular species i
!                              units kT = 1
!***-------------------------------------------------------------------
!****v* chemistry_mod/mu_solvent
! VARIABLE
!    real(long)   mu_solvent(:)     - dimension(N_solvent) 
!                 mu_solvent(i)     = chemical potential of solvent molecular species i
!                              units kT = 1
!*** ------------------------------------------------------------------
!****v* chemistry_mod/interaction_type
! VARIABLE
!    character(20) interaction_type -  method of specifying interaction
!                                   chi   -> specify bare chi, 
!                                   chi_T -> chi(T) = chi_A/T + chi_B
!*** ------------------------------------------------------------------
!****v* chemistry_mod/chi_A
! VARIABLE
!    real(long)   chi_A(:,:)    - dimension(N_monomer,N_monomer)
!                 chi(i,j;T)    = chi_A(i,j)/T + chi_B(i,j)
!*** ------------------------------------------------------------------
!****v* chemistry_mod/chi_B
! VARIABLE
!    real(long)   chi_B(:,:)    - dimension(N_monomer,N_monomer)
!                 chi(i,j;T)    = chi_A(i,j)/T + chi_B(i,j)
!*** ------------------------------------------------------------------
!****v* chemistry_mod/temperature
! VARIABLE
!    real(long)   temperature   = absolute temperature
!*** ------------------------------------------------------------------
!

contains

   !------------------------------------------------------------------
   !****p* chemistry_mod/input_chemistry
   ! SUBROUTINE
   !    input_chemistry(i_unit,fmt_flag)
   ! PURPOSE
   !    Read information about system chemistry and composition.
   !    Allocate public allocatable arrays if necessary
   ! SOURCE
   !------------------------------------------------------------------
   subroutine input_chemistry(i_unit,fmt_flag)
   use io_mod
   integer,      intent(IN) :: i_unit        ! input file unit #
   character(1), intent(IN) :: fmt_flag      ! F = formatted, U = un...
   !***

   integer                 :: i, j          ! loop indices
   character(len = 100)    :: comment_line  
   real(long)              :: tot_vol_frac

   select case(fmt_flag)
   case('F')

      ! Set io defaults
      call set_io_units(i=i_unit,o=6)
      call set_echo(1)                  ! echo inputs
      call set_com_style('A','A','A')   ! comment above data
      call set_com_use('R')             ! replace comment in echo

      ! Monomer Properties
      call input(N_monomer,'N_monomer')
      allocate(kuhn(N_monomer))
      allocate(chi(N_monomer,N_monomer))
      allocate(chi_A(N_monomer,N_monomer))
      allocate(chi_B(N_monomer,N_monomer))
      call input(kuhn,N_monomer,'kuhn',s='R')

      ! Flory-Huggins chi 
      call input(interaction_type,'interaction_type')
      select case(interaction_type)
      case('chi')
         call input(chi,N_monomer,N_monomer,'chi',s='L')
      case('chi_T')
         call input(chi_A,N_monomer,N_monomer,'chi_A',s='L')
         call input(chi_B,N_monomer,N_monomer,'chi_B',s='L')
         call input(temperature,'temperature')
         chi = chi_A/temperature + chi_B
      case default
         write(6,*) 'Error: Illegal interaction_type in input_chemistry'
         stop
      end select
      
      ! Polymer Properties
      call input(N_chain, 'N_chain')
      if (N_chain > 0) then
         allocate(phi_chain(N_chain))
         allocate(mu_chain(N_chain)) 
         allocate(N_block(N_chain))
         call input(N_block,N_chain,'N_block',s='C')
         N_blk_max = 0
         do i = 1, N_chain
            if (N_block(i) > N_blk_max) then
               N_blk_max = N_block(i)
            endif
         enddo
         allocate(block_monomer(N_blk_max,N_chain))
         allocate(block_length(N_blk_max,N_chain))
         read(i_unit,*) comment_line
         call output(trim(comment_line),f='N',j='L')
         do j = 1, N_chain
            call input(block_monomer(:,j),N_block(j),f='N')
         enddo
         read(i_unit,*) comment_line
         call output(trim(comment_line),f='N',j='L')
         do j = 1, N_chain
            call input(block_length(:,j),N_block(j),f='N')
         enddo
         allocate(chain_length(N_chain))
         do j = 1, N_chain
            chain_length(j) = 0.0_long
            do i=1, N_block(j)
               chain_length(j)  = chain_length(j) + block_length(i,j)
            enddo
         enddo
      endif
      
      ! Solvent properties
      call input(N_solvent,'N_solvent')
      if (N_solvent > 0) then
         allocate(phi_solvent(N_solvent))
         allocate(mu_solvent(N_solvent))    
         allocate(solvent_monomer(N_solvent))
         allocate(solvent_size(N_solvent))
         call input(solvent_monomer,N_solvent,'solvent_monomer',s='R')
         call input(solvent_size,N_solvent,'solvent_size',s='R')
      endif
       
      ! Select ensemble
      call input(ensemble,'ensemble')
      select case (ensemble)
      case (0) ! Canonical ensemble 
       
         tot_vol_frac = 0.0_long
         if (N_chain > 0) then
            call input(phi_chain,N_chain,'phi_chain',s='C')
            tot_vol_frac = tot_vol_frac + sum(phi_chain)
         endif    
         if (N_solvent > 0) then
            call input(phi_solvent,N_solvent,'phi_solvent',s='C')
            tot_vol_frac = tot_vol_frac + sum(phi_solvent)
         endif    
         if ( abs(tot_vol_frac - 1.0_long) < 0.02) then
            if ( N_chain > 0 ) phi_chain = phi_chain/tot_vol_frac
            if ( N_solvent > 0 ) phi_solvent = phi_solvent/tot_vol_frac 
         else
            write(6,*) 'Error: sum(phi) /= 1 in read chemistry'
            stop
         endif

      case (1)  ! Grand canonical ensemble 

         if (N_chain > 0) then
            call input(mu_chain,N_chain,'mu_chain',s='C')
         endif    
         if (N_solvent > 0) then
            call input(mu_solvent,N_solvent,'mu_solvent',s='C')
         endif    

      case default

         write(6,*) 'Invalid ensemble =',ensemble,' in input_chemistry'
         stop

      end select

   case('U')

      ! Monomer properties
      read(i_unit) N_monomer
      allocate(kuhn(N_monomer))
      allocate(chi(N_monomer,N_monomer))
      allocate(chi_A(N_monomer,N_monomer))
      allocate(chi_B(N_monomer,N_monomer))
      read(i_unit) kuhn
      ! Chi parameters
      read(i_unit) interaction_type
      select case(trim(interaction_type))
      case('chi')
         read(i_unit) chi
      case('chi_T')
         read(i_unit) chi_A
         read(i_unit) chi_B
         read(i_unit) temperature
         chi = chi_A/temperature + chi_B
      case default
         write(6,*) 'Error: Illegal interaction_type in input_chemistry'
         stop
      end select

      ! Polymer properties
      read(i_unit) N_chain
      if (N_chain > 0) then
         allocate(phi_chain(N_chain))
         allocate(mu_chain(N_chain))
         allocate(N_block(N_chain))
         read(i_unit) N_block
         N_blk_max = 0
         do i = 1, N_chain
            if (N_block(i) .gt. N_blk_max) then
               N_blk_max = N_block(i)
            endif
         enddo
         allocate(block_monomer(N_blk_max,N_chain))
         allocate(block_length(N_blk_max,N_chain))
         do j=1, N_chain
            read(i_unit) block_monomer(1:N_block(j),j)
         enddo
         do j=1, N_chain
            read(i_unit) block_length(1:N_block(j),j)
         enddo
      endif
      ! Solvent properties
      read(i_unit) N_solvent 
      if (N_solvent > 0) then
         allocate(phi_solvent(N_solvent))
         allocate(mu_solvent(N_solvent))
         allocate(solvent_monomer(N_solvent))
         allocate(solvent_size(N_solvent))
         do j=1, N_solvent
            read(i_unit) solvent_monomer(j)
         enddo
         do j=1, N_solvent
            read(i_unit) solvent_size(j)
         enddo
      endif
 
      ! Select ensemble
      read(i_unit) ensemble
      select case (ensemble)
      case (0) ! Canonical 

         tot_vol_frac = 0.0_long
         if (N_chain > 0) then
            read(i_unit) phi_chain
            tot_vol_frac = tot_vol_frac + sum(phi_chain)
         endif
         if (N_solvent > 0) then
            read(i_unit) phi_solvent
            tot_vol_frac = tot_vol_frac + sum(phi_solvent)
         endif
         if ( abs(tot_vol_frac - 1.0_long) < 0.02) then
            if( N_chain > 0 ) phi_chain   = phi_chain/tot_vol_frac
            if( N_solvent > 0) phi_solvent = phi_solvent/tot_vol_frac
         else
            write(6,*) 'Error: sum(phi) /= 1 in read chemistry'
            stop
         endif

      case (1)  ! Grand canonical 

         if (N_chain > 0) then
            read(i_unit) mu_chain
         endif
         if (N_solvent > 0) then
            read(i_unit) mu_solvent 
         endif

      case default

         write(6,*) 'Error: Invalid ensemble =',ensemble,' in input_chemistry'
         stop

      end select

   case default

      write(6,*) 'Error: Illegal fmt_flag in input_chemistry'
      stop

   end select

   ! Check block_monomer indices
   do j=1, N_chain
      do i=1, N_block(j)
         if ( block_monomer(i,j) <= 0) then
            write(6,&
            "('Error: block_monomer(',i2,',',i2,') <=0 ')") i,j
            stop
         endif
         if ( block_monomer(i,j) > N_monomer ) then
            write(6,&
            "('Error: block_monomer(',i2,',',i2,') > N_monomer(',i2,')')") &
            i,j,N_monomer
            stop
         endif
      enddo
   enddo
   
   end subroutine input_chemistry
   !===================================================================


   !------------------------------------------------------------------
   !****p* chemistry_mod/input_monomers
   ! SUBROUTINE
   !    input_monomers(i_unit,fmt_flag)
   ! PURPOSE
   !    Read N_monomer and kuhn
   !    Allocate chi, chi_A and chi_B arrays
   ! SOURCE
   !------------------------------------------------------------------
   subroutine input_monomers(i_unit,fmt_flag)
   use io_mod
   integer,      intent(IN) :: i_unit   ! input file unit #
   character(1), intent(IN) :: fmt_flag ! F = formatted, U = unformatted
   !***

   select case(fmt_flag)
   case('F')

      ! Set io defaults
      call set_io_units(i=i_unit,o=6)
      call set_echo(1)                  ! echo inputs
      call set_com_style('A','A','A')   ! comment above data
      call set_com_use('R')             ! replace comment in echo

      ! Monomer Properties
      call input(N_monomer,'N_monomer')
      allocate(kuhn(N_monomer))
      call input(kuhn,N_monomer,'kuhn',s='R')

   case('U')

      ! Monomer properties
      read(i_unit) N_monomer
      allocate(kuhn(N_monomer))
      read(i_unit) kuhn

   case default

      write(6,*) 'Error: Illegal fmt_flag in input_monomers'
      stop

   end select

   allocate(chi(N_monomer,N_monomer))
   allocate(chi_A(N_monomer,N_monomer))
   allocate(chi_B(N_monomer,N_monomer))

   end subroutine input_monomers
   !===================================================================


   !--------------------------------------------------------------------
   !****p* chemistry_mod/input_chains
   ! SUBROUTINE
   !    input_chains(i_unit,fmt_flag)
   ! PURPOSE
   !    Input information about linear block copolymers and homopolymers
   !    Read N_chain, N_block, block_monomer, block_length
   !    Allocate arrays with dimensions N_chain and N_block
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine input_chains(i_unit,fmt_flag)
   use io_mod
   integer,      intent(IN) :: i_unit   ! input file unit #
   character(1), intent(IN) :: fmt_flag ! F = formatted, U = unformatted
   !***

   integer              :: i, j         ! loop indices
   character(len = 100) :: comment_line  

   select case(fmt_flag)
   case('F')

      ! Set io defaults
      call set_io_units(i=i_unit,o=6)
      call set_echo(1)                  ! echo inputs
      call set_com_style('A','A','A')   ! comment above data
      call set_com_use('R')             ! replace comment in echo

      call input(N_chain, 'N_chain')
      if (N_chain > 0) then
         allocate(phi_chain(N_chain))
         allocate(mu_chain(N_chain)) 
         allocate(N_block(N_chain))
         call input(N_block,N_chain,'N_block',s='C')
         N_blk_max = 0
         do i = 1, N_chain
            if (N_block(i) > N_blk_max) then
               N_blk_max = N_block(i)
            endif
         enddo
         allocate(block_monomer(N_blk_max,N_chain))
         allocate(block_length(N_blk_max,N_chain))
         read(i_unit,*) comment_line
         call output(trim(comment_line),f='N',j='L')
         do j = 1, N_chain
            call input(block_monomer(:,j),N_block(j),f='N')
         enddo
         read(i_unit,*) comment_line
         call output(trim(comment_line),f='N',j='L')
         do j = 1, N_chain
            call input(block_length(:,j),N_block(j),f='N')
         enddo
         allocate(chain_length(N_chain))
         do j = 1, N_chain
            chain_length(j) = 0.0_long
            do i=1, N_block(j)
               chain_length(j)  = chain_length(j) + block_length(i,j)
            enddo
         enddo
      endif
      
   case('U')

      ! Polymer properties
      read(i_unit) N_chain
      if (N_chain > 0) then
         allocate(phi_chain(N_chain))
         allocate(mu_chain(N_chain))
         allocate(N_block(N_chain))
         read(i_unit) N_block
         N_blk_max = 0
         do i = 1, N_chain
            if (N_block(i) .gt. N_blk_max) then
               N_blk_max = N_block(i)
            endif
         enddo
         allocate(block_monomer(N_blk_max,N_chain))
         allocate(block_length(N_blk_max,N_chain))
         do j=1, N_chain
            read(i_unit) block_monomer(1:N_block(j),j)
         enddo
         do j=1, N_chain
            read(i_unit) block_length(1:N_block(j),j)
         enddo
      endif

   case default

      write(6,*) 'Error: Illegal fmt_flag in input_chains'
      stop

   end select

   ! Check block_monomer indices
   do j=1, N_chain
      do i=1, N_block(j)
         if ( block_monomer(i,j) <= 0) then
            write(6,&
            "('Error: block_monomer(',i2,',',i2,') <=0 ')") i,j
            stop
         endif
         if ( block_monomer(i,j) > N_monomer ) then
            write(6,&
            "('Error: block_monomer(',i2,',',i2,') > N_monomer(',i2,')')") &
            i,j,N_monomer
            stop
         endif
      enddo
   enddo
   
   end subroutine input_chains
   !===================================================================


   !------------------------------------------------------------------
   !****p* chemistry_mod/input_solvents
   ! SUBROUTINE
   !    input_solvents(i_unit,fmt_flag)
   ! PURPOSE
   !    Read information about solvent molecules.
   !    Read N_solvents, solvent_monomer, solvent_size
   ! SOURCE
   !------------------------------------------------------------------
   subroutine input_solvents(i_unit,fmt_flag)
   use io_mod
   integer,      intent(IN) :: i_unit   ! input file unit #
   character(1), intent(IN) :: fmt_flag ! F = formatted, U = unformatted
   !***

   integer :: j

   select case(fmt_flag)
   case('F')

      ! Set io defaults
      call set_io_units(i=i_unit,o=6)
      call set_echo(1)                  ! echo inputs
      call set_com_style('A','A','A')   ! comment above data
      call set_com_use('R')             ! replace comment in echo

      ! Solvent properties
      call input(N_solvent,'N_solvent')
      if (N_solvent > 0) then
         allocate(phi_solvent(N_solvent))
         allocate(mu_solvent(N_solvent))    
         allocate(solvent_monomer(N_solvent))
         allocate(solvent_size(N_solvent))
         call input(solvent_monomer,N_solvent,'solvent_monomer',s='R')
         call input(solvent_size,N_solvent,'solvent_size',s='R')
      endif
       
   case('U')

      ! Solvent properties
      read(i_unit) N_solvent 
      if (N_solvent > 0) then
         allocate(phi_solvent(N_solvent))
         allocate(mu_solvent(N_solvent))
         allocate(solvent_monomer(N_solvent))
         allocate(solvent_size(N_solvent))
         do j=1, N_solvent
            read(i_unit) solvent_monomer(j)
         enddo
         do j=1, N_solvent
            read(i_unit) solvent_size(j)
         enddo
      endif
 
   case default

      write(6,*) 'Error: Illegal fmt_flag in input_solvents'
      stop

   end select
   
   end subroutine input_solvents
   !===================================================================


   !------------------------------------------------------------------
   !****p* chemistry_mod/input_composition
   ! SUBROUTINE
   !    input_composition(i_unit,fmt_flag)
   ! PURPOSE
   !    Read information about system composition.
   !    Allocate public allocatable arrays if necessary
   ! SOURCE
   !------------------------------------------------------------------
   subroutine input_composition(i_unit,fmt_flag)
   use io_mod
   integer,      intent(IN) :: i_unit   ! input file unit #
   character(1), intent(IN) :: fmt_flag ! F = formatted, U = unformatted
   !***

   integer                 :: i, j      ! loop indices
   real(long)              :: tot_vol_frac

   select case(fmt_flag)
   case('F')

      ! Set io defaults
      call set_io_units(i=i_unit,o=6)
      call set_echo(1)                  ! echo inputs
      call set_com_style('A','A','A')   ! comment above data
      call set_com_use('R')             ! replace comment in echo

      ! Select ensemble
      call input(ensemble,'ensemble')

      select case (ensemble)
      case (0) ! Canonical ensemble 
       
         tot_vol_frac = 0.0_long
         if (N_chain > 0) then
            call input(phi_chain,N_chain,'phi_chain',s='C')
            tot_vol_frac = tot_vol_frac + sum(phi_chain)
         endif    
         if (N_solvent > 0) then
            call input(phi_solvent,N_solvent,'phi_solvent',s='C')
            tot_vol_frac = tot_vol_frac + sum(phi_solvent)
         endif    
         if ( abs(tot_vol_frac - 1.0_long) < 0.02) then
            if ( N_chain > 0 ) phi_chain = phi_chain/tot_vol_frac
            if ( N_solvent > 0 ) phi_solvent = phi_solvent/tot_vol_frac 
         else
            write(6,*) 'Error: sum(phi) /= 1 in read chemistry'
            stop
         endif

      case (1)  ! Grand canonical ensemble 

         if (N_chain > 0) then
            call input(mu_chain,N_chain,'mu_chain',s='C')
         endif    
         if (N_solvent > 0) then
            call input(mu_solvent,N_solvent,'mu_solvent',s='C')
         endif    

      case default

         write(6,*) 'Invalid ensemble =',ensemble,' in input_composition'
         stop

      end select

   case('U')

      ! Select ensemble
      read(i_unit) ensemble

      select case (ensemble)
      case (0) ! Canonical 

         tot_vol_frac = 0.0_long
         if (N_chain > 0) then
            read(i_unit) phi_chain
            tot_vol_frac = tot_vol_frac + sum(phi_chain)
         endif
         if (N_solvent > 0) then
            read(i_unit) phi_solvent
            tot_vol_frac = tot_vol_frac + sum(phi_solvent)
         endif
         if ( abs(tot_vol_frac - 1.0_long) < 0.02) then
            if( N_chain > 0 ) phi_chain   = phi_chain/tot_vol_frac
            if( N_solvent > 0) phi_solvent = phi_solvent/tot_vol_frac
         else
            write(6,*) 'Error: sum(phi) /= 1 in read chemistry'
            stop
         endif

      case (1)  ! Grand canonical 

         if (N_chain > 0) then
            read(i_unit) mu_chain
         endif
         if (N_solvent > 0) then
            read(i_unit) mu_solvent 
         endif

      case default

         write(6,*) 'Error: Invalid ensemble =',ensemble,' in input_composition'
         stop

      end select

   case default

      write(6,*) 'Error: Illegal fmt_flag in input_composition'
      stop

   end select

   end subroutine input_composition
   !===================================================================


   !------------------------------------------------------------------
   !****p* chemistry_mod/input_interaction
   ! SUBROUTINE
   !    input_interaction(i_unit,fmt_flag)
   ! PURPOSE
   !    Read information about interaction free energy
   !    Allocate public allocatable arrays if necessary
   ! SOURCE
   !------------------------------------------------------------------
   subroutine input_interaction(i_unit,fmt_flag)
   use io_mod
   integer,      intent(IN) :: i_unit        ! input file unit #
   character(1), intent(IN) :: fmt_flag      ! F = formatted, U = un...
   !***

   select case(fmt_flag)
   case('F')

      ! Set io defaults
      call set_io_units(i=i_unit,o=6)
      call set_echo(1)                  ! echo inputs
      call set_com_style('A','A','A')   ! comment above data
      call set_com_use('R')             ! replace comment in echo

      ! Flory-Huggins chi 
      call input(interaction_type,'interaction_type')
      select case(trim(interaction_type))
      case('chi')
         call input(chi,N_monomer,N_monomer,'chi',s='L')
      case('chi_T')
         call input(chi_A,N_monomer,N_monomer,'chi_A',s='L')
         call input(chi_B,N_monomer,N_monomer,'chi_B',s='L')
         call input(temperature,'temperature')
         chi = chi_A/temperature + chi_B
      case default
         write(6,*) 'Error: Illegal interaction_type in input_interaction'
         stop
      end select
      
   case('U')

      ! Chi parameters
      read(i_unit) interaction_type
      select case(trim(interaction_type))
      case('chi')
         read(i_unit) chi
      case('chi_T')
         read(i_unit) chi_A
         read(i_unit) chi_B
         read(i_unit) temperature
         chi = chi_A/temperature + chi_B
      case default
         write(6,*) 'Error: Illegal interaction_type in input_interaction'
         stop
      end select

   case default

      write(6,*) 'Error: Illegal fmt_flag in input_interaction'
      stop

   end select

   end subroutine input_interaction
   !===================================================================


   !------------------------------------------------------------------
   !****p* chemistry_mod/output_chemistry
   ! SUBROUTINE
   !    output_chemistry(o_unit,fmt_flag)
   ! PURPOSE
   !    Write information about system chemistry and composition.
   ! SOURCE
   !------------------------------------------------------------------
   subroutine output_chemistry(o_unit,fmt_flag)
   use io_mod
   integer, intent(IN)          :: o_unit    ! input unit #
   character(len=1), intent(IN) :: fmt_flag  ! F = formatted, U = un...
   !***

   integer                 :: i, j    ! loop indices

   select case(fmt_flag)
   case('F') ! formatted for input

      call set_io_units(o=o_unit)
      call output(N_monomer,'N_monomer')
      call output(kuhn,N_monomer,'kuhn')
      call output(interaction_type,'interaction_type')
      select case (trim(interaction_type))
      case('chi')
         call output(chi,N_monomer,N_monomer,'chi','L')
      case('chi_T')
         call output(chi_A,N_monomer,N_monomer,'chi_A','L')
         call output(chi_B,N_monomer,N_monomer,'chi_B','L')
         call output(temperature,'temperature')
      end select

      ! Polymer properties
      call output(N_chain,'N_chain')
      if (N_chain > 0) then
         call output(N_block,N_chain,'N_block','C')
         call output('block_monomer',f='N',j='L')
         do j = 1,N_chain
            call output(block_monomer(:,j),N_block(j),f='N')
         enddo
         call output('block_length',f='N',j='L')
         do j = 1,N_chain
            call output(block_length(:,j),N_block(j),f='N')
         enddo
      endif
      
      ! Solvent properties
      call output(N_solvent,'N_solvent')
      if (N_solvent > 0) then
         call output(solvent_monomer,N_solvent,'solvent_monomer')
         call output(solvent_size,N_solvent,'solvent_size')
      endif

      ! Select ensemble
      call output(ensemble,'ensemble')
      select case (ensemble)
      case (0) ! Canonical ensemble 
         if (N_chain > 0) then
            call output(phi_chain,N_chain,'phi_chain','C')
         endif
         if (N_solvent > 0) then
            call output(phi_solvent,N_solvent,'phi_solvent','C')
         endif
      case (1)  ! Grand canonical ensemble 
         if (N_chain > 0) then
            call output(mu_chain,N_chain,'mu_chain','C')
         endif
         if (N_solvent > 0) then
            call output(mu_solvent,N_solvent,'mu_solvent','C')
         endif
      end select

   case('U')

      write(o_unit) N_monomer
      write(o_unit) kuhn
      write(o_unit) interaction_type
      select case (trim(interaction_type))
      case('chi')
         write(o_unit) chi
      case('chi_T')
         write(o_unit) chi_A
         write(o_unit) chi_B
         write(o_unit) temperature
      end select
      write(o_unit) N_chain
      if (N_chain > 0) then
         write(o_unit) N_block
         do j=1, N_chain
            write(o_unit) block_monomer(1:N_block(j),j)
         enddo
         do j=1, N_chain
            write(o_unit) block_length(1:N_block(j),j)
         enddo
      endif
      write(o_unit) N_solvent
      if (N_solvent > 0) then
         do j=1, N_solvent
            write(o_unit) solvent_monomer(j)
         enddo
         do j=1, N_solvent
            write(o_unit) solvent_size(j)
         enddo
      endif
      write(o_unit) ensemble
      select case (ensemble)
      case (0) ! Canonical ensemble 
         if (N_chain > 0) then
            write(o_unit) phi_chain
         endif
         if (N_solvent > 0) then
            write(o_unit) phi_solvent
         endif
      case (1)  ! Grand canonical ensemble 
         if (N_chain > 0) then
            write(o_unit) mu_chain
         endif
         if (N_solvent > 0) then
            write(o_unit) mu_solvent
         endif
      end select

   case default

      write(6,*) 'Invalid format flag in output_chemistry'

   end select

   end subroutine output_chemistry
   !==================================================================
 

   !------------------------------------------------------------------
   !****p* chemistry_mod/output_monomers
   ! SUBROUTINE
   !    output_monomers(o_unit,fmt_flag)
   ! PURPOSE
   !    Output N_monomer and kuhn
   !    Allocate chi, chi_A and chi_B arrays
   ! SOURCE
   !------------------------------------------------------------------
   subroutine output_monomers(o_unit,fmt_flag)
   use io_mod
   integer,      intent(IN) :: o_unit   ! output file unit #
   character(1), intent(IN) :: fmt_flag ! F = formatted, U = unformatted
   !***

   select case(fmt_flag)
   case('F')

      ! Set io defaults
      call set_io_units(o=o_unit)
      call set_com_style('A','A','A')   ! comment above data

      ! Monomer Properties
      call output(N_monomer,'N_monomer')
      call output(kuhn,N_monomer,'kuhn',s='R')

   case('U')

      ! Monomer properties
      write(o_unit) N_monomer
      write(o_unit) kuhn

   case default

      write(6,*) 'Error: Illegal fmt_flag in output_monomers'
      stop

   end select

   end subroutine output_monomers
   !===================================================================


   !--------------------------------------------------------------------
   !****p* chemistry_mod/output_chains
   ! SUBROUTINE
   !    output_chains(o_unit,fmt_flag)
   ! PURPOSE
   !    Input information about linear block copolymers and homopolymers
   !    Allocate related arrays
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine output_chains(o_unit,fmt_flag)
   use io_mod
   integer,      intent(IN) :: o_unit   ! output file unit #
   character(1), intent(IN) :: fmt_flag ! F = formatted, U = unformatted
   !***

   integer              :: i, j         ! loop indices

   select case(fmt_flag)
   case('F')

      ! Set io defaults
      call set_io_units(o=o_unit)
      call set_com_style('A','A','A')   ! comment above data

      call output(N_chain, 'N_chain')
      if (N_chain > 0) then
         call output(N_block,N_chain,'N_block',s='C')
         call output('block_monomer',f='N',j='L')
         do j = 1, N_chain
            call output(block_monomer(:,j),N_block(j),f='N')
         enddo
         call output('block_length',f='N',j='L')
         do j = 1, N_chain
            call output(block_length(:,j),N_block(j),f='N')
         enddo
      endif
      
   case('U')

      ! Polymer properties
      write(o_unit) N_chain
      if (N_chain > 0) then
         write(o_unit) N_block
         N_blk_max = 0
         do i = 1, N_chain
            if (N_block(i) .gt. N_blk_max) then
               N_blk_max = N_block(i)
            endif
         enddo
         do j=1, N_chain
            write(o_unit) block_monomer(1:N_block(j),j)
         enddo
         do j=1, N_chain
            write(o_unit) block_length(1:N_block(j),j)
         enddo
      endif

   case default

      write(6,*) 'Error: Illegal fmt_flag in output_chains'
      stop

   end select
   
   end subroutine output_chains
   !===================================================================


   !------------------------------------------------------------------
   !****p* chemistry_mod/output_solvents
   ! SUBROUTINE
   !    output_solvents(o_unit,fmt_flag)
   ! PURPOSE
   !    Output information about solvent identities and volumes
   !    Allocate public allocatable arrays if necessary
   ! SOURCE
   !------------------------------------------------------------------
   subroutine output_solvents(o_unit,fmt_flag)
   use io_mod
   integer,      intent(IN) :: o_unit   ! output file unit #
   character(1), intent(IN) :: fmt_flag ! F = formatted, U = unformatted
   !***

   integer :: j

   select case(fmt_flag)
   case('F')

      ! Set io defaults
      call set_io_units(o=o_unit)
      call set_com_style('A','A','A')   ! comment above data

      ! Solvent properties
      call output(N_solvent,'N_solvent')
      if (N_solvent > 0) then
         call output(solvent_monomer,N_solvent,'solvent_monomer',s='R')
         call output(solvent_size,N_solvent,'solvent_size',s='R')
      endif
       
   case('U')

      ! Solvent properties
      write(o_unit) N_solvent 
      if (N_solvent > 0) then
         do j=1, N_solvent
            write(o_unit) solvent_monomer(j)
         enddo
         do j=1, N_solvent
            write(o_unit) solvent_size(j)
         enddo
      endif
 
   case default

      write(6,*) 'Error: Illegal fmt_flag in output_solvents'
      stop

   end select
   
   end subroutine output_solvents
   !===================================================================


   !------------------------------------------------------------------
   !****p* chemistry_mod/output_composition
   ! SUBROUTINE
   !    output_composition(o_unit,fmt_flag)
   ! PURPOSE
   !    Output information about system composition.
   !    Allocate public allocatable arrays if necessary
   ! SOURCE
   !------------------------------------------------------------------
   subroutine output_composition(o_unit,fmt_flag)
   use io_mod
   integer,      intent(IN) :: o_unit   ! output file unit #
   character(1), intent(IN) :: fmt_flag ! F = formatted, U = unformatted
   !***

   integer                 :: i, j      ! loop indices

   select case(fmt_flag)
   case('F')

      ! Set io defaults
      call set_io_units(o=o_unit)
      call set_com_style('A','A','A')   ! comment above data

      ! Select ensemble
      call output(ensemble,'ensemble')

      select case (ensemble)
      case (0) ! Canonical ensemble 
       
         if (N_chain > 0) then
            call output(phi_chain,N_chain,'phi_chain',s='C')
         endif    
         if (N_solvent > 0) then
            call output(phi_solvent,N_solvent,'phi_solvent',s='C')
         endif    

      case (1)  ! Grand canonical ensemble 

         if (N_chain > 0) then
            call output(mu_chain,N_chain,'mu_chain',s='C')
         endif    
         if (N_solvent > 0) then
            call output(mu_solvent,N_solvent,'mu_solvent',s='C')
         endif    

      case default

         write(6,*) 'Invalid ensemble =',ensemble,' in output_composition'
         stop

      end select

   case('U')

      ! Select ensemble
      write(o_unit) ensemble

      select case (ensemble)
      case (0) ! Canonical 

         if (N_chain > 0) then
            write(o_unit) phi_chain
         endif
         if (N_solvent > 0) then
            write(o_unit) phi_solvent
         endif

      case (1)  ! Grand canonical 

         if (N_chain > 0) then
            write(o_unit) mu_chain
         endif
         if (N_solvent > 0) then
            write(o_unit) mu_solvent 
         endif

      case default

         write(6,*) 'Error: Invalid ensemble =',ensemble,' in output_composition'
         stop

      end select

   case default

      write(6,*) 'Error: Illegal fmt_flag in output_composition'
      stop

   end select

   end subroutine output_composition
   !===================================================================


   !------------------------------------------------------------------
   !****p* chemistry_mod/output_interaction
   ! SUBROUTINE
   !    output_interaction(o_unit,fmt_flag)
   ! PURPOSE
   !    Output information about interaction free energy
   !    Allocate public allocatable arrays if necessary
   ! SOURCE
   !------------------------------------------------------------------
   subroutine output_interaction(o_unit,fmt_flag)
   use io_mod
   integer,      intent(IN) :: o_unit        ! output file unit #
   character(1), intent(IN) :: fmt_flag      ! F = formatted, U = un...
   !***

   select case(fmt_flag)
   case('F')

      ! Set io defaults
      call set_io_units(o=o_unit)
      call set_com_style('A','A','A')   ! comment above data

      ! Flory-Huggins chi 
      call output(interaction_type,'interaction_type')
      select case(trim(interaction_type))
      case('chi')
         call output(chi,N_monomer,N_monomer,'chi',s='L')
      case('chi_T')
         call output(chi_A,N_monomer,N_monomer,'chi_A',s='L')
         call output(chi_B,N_monomer,N_monomer,'chi_B',s='L')
         call output(temperature,'temperature')
         chi = chi_A/temperature + chi_B
      case default
         write(6,*) 'Error: Illegal interaction_type in output_interaction'
         stop
      end select
      
   case('U')

      ! Chi parameters
      write(o_unit) interaction_type
      select case(trim(interaction_type))
      case('chi')
         write(o_unit) chi
      case('chi_T')
         write(o_unit) chi_A
         write(o_unit) chi_B
         write(o_unit) temperature
         chi = chi_A/temperature + chi_B
      case default
         write(6,*) 'Error: Illegal interaction_type in output_interaction'
         stop
      end select

   case default

      write(6,*) 'Error: Illegal fmt_flag in output_interaction'
      stop

   end select

   end subroutine output_interaction
   !===================================================================


   !------------------------------------------------------------------
   !****p* chemistry_mod/rescale_vref
   ! SUBROUTINE
   !    rescale_vref(scale)
   ! PURPOSE
   !    Rescale reference volume, such that:
   !       vref         -> vref/scale 
   !       chi          -> chi/scale
   !       kuhn         -> kuhn/sqrt(scale)
   !       block_length -> block_length*scale
   !    Omega is not rescaled by routine, but must be scaled as:
   !       omega        -> omega/scale
   ! SOURCE
   !------------------------------------------------------------------
   subroutine rescale_vref(scale)
   real(long) :: scale
   !***
   kuhn = kuhn/sqrt(scale)
   select case(trim(interaction_type))
   case('chi')
      chi = chi/scale
   case('chi_T')
      chi_A = chi_A/scale
      chi_B = chi_B/scale
   end select
   if (N_chain > 0) then
      block_length = block_length*scale
   endif
   if (N_solvent > 0) then
      solvent_size = solvent_size*scale
   endif
   end subroutine rescale_vref
   !===================================================================

 end module chemistry_mod
