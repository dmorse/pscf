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
!----------------------------------------------------------------------
!****e scf/pscf_pd
! PROGRAM
!   pscf_pd
! PURPOSE
!   Main program for polymer self-consistent field theory (PSCF) for 
!   spatially periodic (PD) microstructures. The program can treat 
!   incompressible mixtures of any number of linear multiblock 
!   copolymers, linear homopolymers, and small molecule solvents.
!
!   This program contains a main loop that reads an input script from
!   standard input. The loop, and the operation of the program, ends
!   when a line containing the word FINISH is encountered. The manual 
!   contains detailed description of the format of this script, but a
!   brief description is also given here:
!   
!   Input Script Format:
!   The first line of the input script must be of the "format i j", in 
!   which i and j are major and minor version numbers for the input 
!   script file format (e.g., "format 1 0" for v1.0). The rest of the 
!   input script consists of a series of blocks of input data. Each
!   block begins with a blank line followed by a line containing a 
!   capitalized operation flag (op_flag) string, such as 'CHEMISTRY', 
!   'UNIT_CELL', 'DISCRETIZATION', 'ITERATION', etc.  The remainder 
!   of each block (if any) is a series of variable values. Each 
!   variable is preceded by a comment line containing the name of 
!   the variable, as used in the program. One and two-dimensional
!   array variables may be input in any of several formats (e.g., all
!   values on a single line, or on multiple lines). See the manual
!   for the format of such variables.  Some operation flags, such 
!   as 'ITERATE' and 'SWEEP', trigger an action, beyond reading in 
!   variables. Reading of the script, and the program, ends when the 
!   'FINISH' flag is encountered.
!
!   Within each block, input variables are read using the input()
!   interface of the io_mod module. This is an overloaded interface
!   that can be used to read in a scalar or array of integer, real,
!   character, or boolean data. The main program uses the comment
!   style 'A' (for 'Above') of io_mod, in which the value of each
!   input variable is preceded by a line containing the name of
!   the variable. Within each block, variables must appear in a 
!   predetetermined format.  The required format is documented in
!   the manual, and in example scripts provided with the program.
!   All data is read using the fortran read(5,*) format, so that
!   the data format and use of white space is flexible, as long 
!   as each data value is of the expected data type. 
!     
! SOURCE
!----------------------------------------------------------------------
program pscf_pd
   use const_mod
   use string_mod
   use io_mod
   use version_mod
   use field_io_mod
   use chemistry_mod, only : input_chemistry, rescale_vref, &
                             input_monomers, input_chains, input_solvents,&
                             input_composition, input_interaction, &
                             output_monomers, output_chains, output_solvents,&
                             output_composition, output_interaction, &
                             N_monomer, N_chain, N_block, N_solvent
   use unit_cell_mod, only : input_unit_cell, output_unit_cell, &
                             N_cell_param, cell_param, &
                             make_unit_cell, R_basis, G_basis
   use group_mod,     only : output_group
   use grid_mod,      only : ngrid, input_grid, allocate_grid, make_ksq
   use basis_mod,     only : N_wave, N_star, group, &
                             make_basis, output_waves, release_basis
   use grid_basis_mod
   use chain_mod
   use scf_mod,       only : density_startup, density
   use iterate_mod,   only : input_iterate_param, output_iterate_param, &
                             iterate_NR_startup, iterate_NR, domain
   use sweep_mod
   use response_mod,  only : response_startup, response_sweep
   !# ifdef COMMENTS
   use spinodal_mod
   !use group_rep_mod
   !# endif
   implicit none

   ! SCFT variables
   real(long)      :: omega(:,:)      ! chemical potential field 
                                      ! omega(monomer,basis function)
   real(long)      :: rho(:,:)        ! monomer density field
                                      ! rho(monomer,basis function)
   real(long)      :: f_Helmholtz     ! free energy, units kT / monomer 
   real(long)      :: pressure        ! pressure * V_monomer / kT 
   real(long)      :: stress(:)       ! d(f_Helmholtz)/d(cell_param)
   !# ifdef DEVEL
   real(long)      :: f_component(4)  !
   real(long)      :: overlap(:,:)    ! overlap integrals
   allocatable     :: omega, rho, stress, overlap
   !# else
   allocatable     :: omega, rho, stress
   !# endif

   ! Input parameters (all others declared in modules)
   character(60)   :: group_name      ! name of crystal space group
   real(long)      :: chain_step      ! contour length step size
   real(long)      :: vref_scale      ! Factor to rescale reference volume

   ! Input and output file names
   character(60)   :: input_prefix    ! prefix for input omega file:
                                      !   input_prefix//omega
   character(60)   :: output_prefix   ! prefix for output files:
                                      !   output_prefix//out
                                      !   output_prefix//omega
                                      !   output_prefix//rho
                                      !   output_prefix//group
                                      !   output_prefix//waves
   character(60)   :: input_file      ! input file for FIELD_TO_GRID
   character(60)   :: output_file     ! output file for FIELD_TO_GRID

   ! Variables for iteration (fixed chemistry)
   integer         :: extr_order      ! extrapolation order
   integer         :: itr             ! iteration index
   real(long)      :: error           ! error = max(residual)
   logical         :: converge        ! true if converged

   ! Variables for sweep (sequence of parameters)
   integer         :: i, j            ! step indices
   real(long)      :: s_max           ! maximum value of variable s
   real(long)      :: s               ! continuation variable
   real(long)      :: step            ! actual step size
   real(long)      :: step_unit       ! unit step size 

   ! Operation selection string from input script
   character(60)   :: op_string      ! Operation selection string

   ! Logical operation flags - set true as each operation requested
   ! Note: These are listed in normal sequence within input file
   logical :: monomer_flag     = .FALSE. ! monomer data read
   logical :: chain_flag       = .FALSE. ! chain data read
   logical :: solvent_flag     = .FALSE. ! solvent data read
   logical :: composition_flag = .FALSE. ! composition data read
   logical :: interaction_flag = .FALSE. ! interaction data read
   logical :: unit_cell_flag   = .FALSE. ! unit_cell made
   logical :: discretize_flag  = .FALSE. ! grid and ds made
   logical :: prefix_flag      = .FALSE. ! io file prefixes read
   logical :: basis_flag       = .FALSE. ! symmetrized basis made
   logical :: omega_flag       = .FALSE. ! initial omega exists
   logical :: iterate_flag     = .FALSE. ! 1st iteration requested
   logical :: output_flag      = .FALSE. ! deferred iterate output
   logical :: sweep_flag       = .FALSE. ! sweep requested

   ! Timing variables
   real(long) :: start_time, basis_time, scf_time
   real(long) :: rpa_time

   ! File Unit numbers (parameters)
   integer, parameter :: out_unit   = 21 ! output summary
   integer, parameter :: field_unit = 22 ! omega and rho fields
   integer            :: ierr            ! error msg for file io

   ! File format version numbers
   type(version_type) :: version      ! input script format
   !------------------------------------------------------------------

   call cpu_time(start_time)

   ! Set defaults for parameter I/O - see io_mod
   call set_echo(1)                ! echo inputs to standard out
   call set_com_style('A','A','A') ! comments on line above data
   call set_com_use('R')           ! replace comment in echoed output
   call set_io_units(i=5,o=6)      ! set standard in and out units

   ! Read file format version from input script, echo to stdout
   call input_version(version, 5)
   call output_version(version, 6)

   ! Main operation loop - infinite loop
   op_loop : do 
 
      ! Read operation string from stdin 
      read(5,*)
      read(5,*) op_string
 
      ! Echo to stdout
      write(6,*)
      call output(trim(op_string),f="N",j="L",o=6)
 
      select case(trim(op_string))
 
      case ("CHEMISTRY")

         ! Input old format for chemistry, in which information about
         ! monomers, chains, solvents, composition, and interaction 
         ! are in one large block. Included temporarily for backwards 
         ! compatability. 
         monomer_flag     = .TRUE.
         chain_flag       = .TRUE.
         solvent_flag     = .TRUE.
         composition_flag = .TRUE.
         interaction_flag = .TRUE.
 
         ! Read chemistry data from stdin ( see chemistry_mod )
         call input_chemistry(5,'F')  

      case ("MONOMERS")

        ! Input N_monomer and kuhn array
        ! See chemistry_mod and users manual
 
         monomer_flag = .TRUE.
         call input_monomers(5,'F')

      case ("CHAINS")
 
        ! Input N_chains, block_monomer, block_length
        ! See chemistry_mod and users manual
 
         ! Check preconditions (must know N_monomer first)
         if (.not.monomer_flag) then
            write(6,*) "Error: Must read MONOMERS before CHAINS" 
            exit op_loop
         end if

         chain_flag  = .TRUE.
         call input_chains(5,'F')

      case ("SOLVENTS")
 
        ! Input N_solvents, solvent_monomer, solvent_size
        ! See chemistry_mod and users manual
 
         ! Check preconditions (must know N_monomer first)
         if (.not.monomer_flag) then
            write(6,*) "Error: Must read MONOMERS before SOLVENTS"
            exit op_loop
         end if
         solvent_flag  = .TRUE.
         call input_solvents(5,'F')

      case ("COMPOSITION")
 
        ! Input ensemble, phi_chain and phi_solvent, or mu_chain and mu_solvent
        ! See chemistry_mod and users manual

         ! Check preconditions (must know N_chain and N_solvent first)
         if ( .not. (chain_flag .OR. solvent_flag) ) then
            write(6,*) &
            "Error: Must read CHAINS and/or SOLVENTS before COMPOSITION"
            exit op_loop
         end if
         composition_flag = .TRUE.

         call input_composition(5,'F')

      case ("INTERACTION")

        ! Input Flory-Huggins chi interaction parameters.
        ! Input interaction_type and:
        !    chi (if interaction_type = 'chi'), or
        !    chi_a, chi_b, & temperature (if interaction_type = 'chi_T')
        ! See chemistry_mod and users manual

         ! Check preconditions (must know N_monomer first)
         if (.not.monomer_flag) then
            write(6,*) "Error: Must read MONOMERS before INTERACTION"
            exit op_loop
         end if
         interaction_flag = .TRUE.

         call input_interaction(5,'F')

      !# ifdef COMMENTS
      !case ('RPA') (Under reconstruction)
      !
      !   if (.not.composition_flag) then
      !      write(6,*) "Error: Must read COMPOSITION before RPA"
      !      exit op_loop
      !   end if
      !   if (.not.interaction_flag) then
      !      write(6,*) "Error: Must read INTERACTION before RPA"
      !      exit op_loop
      !   end if
      !
      !   call rpa_homo_startup
      !   call rpa_homo(5)
      !
      !   call cpu_time(rpa_time)
      !   rpa_time = rpa_time - start_time
      !   call output(rpa_time,'rpa_time',o=6)
      !# endif
 
      case ('UNIT_CELL')
 
         unit_cell_flag = .TRUE.
 
         ! Read unit cell parameters (see unit_cell_mod)
         call input_unit_cell(5,'F')       
 
         ! Construct initial unit cell (see unit_cell_mod)
         call make_unit_cell 
 
      case ('DISCRETIZATION')

         ! Read spatial and contour length discretization of PDE

         ! Check preconditions (Needs N_monomer and unit_cell parameters)
         if (.not.monomer_flag) then
            write(6,*) "Error: Must read MONOMERS before DISCRETIZATION"
            exit op_loop
         else if ( .not. unit_cell_flag ) then
            write(6,*) "Error: Must make UNIT_CELL before DISCRETIZATION"
            exit op_loop
         end if
         discretize_flag = .TRUE.

         ! Input ngrid = number of points in FFT grid in each direction
         call input_grid
         call allocate_grid(N_monomer)
 
         !call input(extr_order,'extr_order')
         extr_order = 1
         call input(chain_step,'chain_step')
 
      case ('FILE_PREFIXES')

         ! Input prefixes used to contruct input and output file names.
         ! Prefixes can include directories paths, if they end with a
         ! trailing directory separator '/'. The input omega file name 
         ! is constructed by appending 'omega' to input_prefix. Output
         ! file names are constructed by appending 'out', 'rho' and
         ! 'omega' to output_prefix. 

         prefix_flag = .TRUE.
         call input(input_prefix, 'input_prefix')  ! input  file prefix
         call input(output_prefix,'output_prefix') ! output file prefix
 
      case ('BASIS')

         ! Construct symmetry-adapated basis functions. See basis_mod

         ! Check preconditions (needs unit cell and grid)
         if ( .not. unit_cell_flag ) then
            write(6,*) "Error: Must make UNIT_CELL before BASIS"
            exit op_loop
         else if (.not.discretize_flag) then
            write(6,*) "Error: Must make DISCRETIZATION before BASIS"
            exit op_loop
         else if ( .not.prefix_flag ) then
            write(6,*) "Error: Must provide FILE_PREFIXES before BASIS"
            exit op_loop
         end if
         basis_flag = .TRUE.
 
         ! Read name of space group used to construct basis functions.
         ! The string group_name can be a space group symbol, a space 
         ! group number, or the name of a file containing the elements 
         ! of the group. See space_groups in module space_group_mod.
         call input(group_name,'group_name')
 
         ! Construct basis functions (see basis_mod)
         call make_basis(R_basis,G_basis,group_name,ngrid,grid_flag=.TRUE.)
 
         ! Output N_star (# of symmetrized basis functions) to stdout
         call output(N_star,'N_star',o=6)
      
         ! Allocate omega, rho, stress (internal routine)
         call allocate_scf_arrays             
      
      case ('RESCALE')
 
         ! Rescale monomer reference volume. Change values of kuhn, 
         ! chi, block_length, and solvent_size parameter arrays, and
         ! the omega field, to obtain an equivalent set of parameters 
         ! and fields.
 
         ! If omega field has not been read previously, read it now
         if (.not.omega_flag) then

            ! Check preconditions for reading omega
            if (.not.basis_flag) then
               write(6,*) "Error: Must make BASIS before RESCALE"
               exit op_loop
            else if ( .not.prefix_flag ) then
               write(6,*) &
                  "Error: Must provide FILE_PREFIXES before RESCALE"
               exit op_loop
            end if

            ! Read omega
            open(unit=field_unit,file=trim(input_prefix)//'omega',&
                              status='old',iostat=ierr)
            if (ierr/=0) stop "Error while opening omega file"
            call input_field(omega,field_unit)
            close(field_unit)
            omega_flag = .TRUE.

         end if
 
         ! Read in scale factor: vref -> vref/vref_scale
         call input(vref_scale,'vref_scale')

         ! Rescale kuhn, chi, block_length, solvent_size 
         ! See chemistry_mod
         call rescale_vref(vref_scale)

         ! Rescale omega field
         omega = omega/vref_scale
 
      case ('ITERATE')

         ! Iterate to convergence for one set of parameters
         ! See iterate_mod
         ! Check preconditions
         if (.not.composition_flag) then
            write(6,*) "Error: Must read COMPOSITION before ITERATION"
            exit op_loop
         else if (.not.interaction_flag) then
            write(6,*) "Error: Must read INTERACTION before ITERATION"
            exit op_loop
         else if (.not.unit_cell_flag) then
            write(6,*) "Error: Must make UNIT_CELL before ITERATION"
            exit op_loop
         else if (.not.discretize_flag) then
            write(6,*) "Error: Must make DISCRETIZATION before ITERATION"
            exit op_loop
         else if (.not.basis_flag) then
            write(6,*) "Error: Must make BASIS before ITERATION"
            exit op_loop
         else if ( .not.prefix_flag ) then
            write(6,*) "Error: Must read FILE_PREFIXES before ITERATION"
            exit op_loop
         end if

         ! Read omega file, if not read previously
         if (.not.omega_flag) then
            open(unit=field_unit,file=trim(input_prefix)//'omega',&
                              status='old',iostat=ierr)
            if (ierr/=0) stop "Error while opening omega source file."
            call input_field(omega,field_unit)
            close(field_unit)
            omega_flag = .TRUE.
         end if
 
         iterate_flag = .TRUE.
 
         call cpu_time(basis_time)
         basis_time = basis_time - start_time
         call cpu_time(start_time)
 
         ! Read parameters for iteration
         call input_iterate_param
 
         ! Allocate and initialize chain objects used in scf_mod
         ! Create fft_plan, which is saved in scf_mod as public variable
         call density_startup(ngrid,extr_order,chain_step,&
                                  update_chain=.false.)
 
         ! Allocate private arrays for Newton-Raphson iteration
         call iterate_NR_startup(N_star)
      
         write(6,FMT = "( / '************************************' / )" )
         ! Main Newton-Raphson iteration loop
         call iterate_NR(      &
                  N_star,      &! # of basis functions
                  omega,       &! chemical potential field (IN/OUT)
                  itr,         &! actual number of interations
                  converge,    &! = .TRUE. if converged
                  error,       &! final error = max(residuals)
                  rho,         &! monomer density field
                  f_Helmholtz, &! Helmholtz free energy per monomer/kT
                  !# ifdef DEVEL
                  f_component, &! free energy components
                  overlap,     &! overlap integrals
                  !# endif
                  pressure,    &! pressure * monomer volume / kT
                  stress       &! d(free energy)/d(cell parameters)
                       )

         ! Defer output to SWEEP, RESPONSE, or FINISH operation
         ! If output by SWEEP, '0.' will be added to output_prefix
         output_flag = .TRUE.
      
         call cpu_time(scf_time)
         scf_time = scf_time - start_time
 
      case ('SWEEP')

         ! Iterate to convergences for a sequence of sets of parameters.
         ! See sweep_mod.
 
         ! Must be preceded by successful ITERATE for first solution.
         if (.not.iterate_flag) then
                 write(6,*) &
                  "Error: Must call ITERATE (1st iteration) before SWEEP"
            exit op_loop
         else if (.not.converge) then
            write(6,*) "Error: 1st iteration failed to converge"
            exit op_loop
         end if
         sweep_flag = .TRUE.
 
         ! Read parameters needed by sweep
         call input(s_max,'s_max')           ! max(contour variable s)  
         call input_increments(5,'N',domain) ! see sweep_mod
        
         ! Output from first iteration, add '0.' to output_prefix
         if (output_flag) then
            call output_summary(trim(output_prefix)//'0.')
            call output_fields(trim(output_prefix)//'0.')
            output_flag = .FALSE.
         endif
 
         ! Initialize contour variable s = 0.0 -> s_max
         s = 0.0_long
   
         ! Initialize history arrays
         call history_setup
         call update_history(s,omega,cell_param,domain)
   
         ! Loop over sweep through parameters
         i         = 0
         step_unit = 1.0_long
         sweep_loop : do 
    
            if (i == 0) then
               step = 0.1*step_unit
            else if (i == 1) then
               if (converge) then
                  step = 0.9*step_unit
               else 
                  step = step_unit - s 
               end if
            else if (i > 1) then
               step = step_unit
            end if
   
            s = s + step
            call increment_parameters(step,domain,cell_param)
   
            write(6, FMT = "('************************************'/ )" )
            write(6, FMT = "('s =',f10.4)" ) s
            call cpu_time(start_time)
         
            ! 1st order continuation of omega and cell_param
            call continuation(step,domain,omega,cell_param)
   
            ! Reconstruct unit cell
            ! if (domain) then
               call make_unit_cell 
               call make_ksq(G_basis)
            ! end if
 
            ! Rebuild chains
            call density_startup(ngrid,extr_order,chain_step,&
                                 update_chain=.TRUE.)
 
            ! Main iteration routine
            call iterate_NR( &
                N_star,      &! # of basis functions
                omega,       &! chemical potential field (IN/OUT)
                itr,         &! actual number of interations
                converge,    &! = .TRUE. if converged
                error,       &! final error = max(residuals)
                rho,         &! monomer density field
                f_Helmholtz, &! Helmholtz free energy per monomer/kT
                !# ifdef DEVEL
                f_component, &! free energy components
                overlap,     &! overlap integrals
                !# endif
                pressure,    &! pressure * monomer volume/kT
                stress       &! d(free energy)/d(cell parameters)
                )
         
            if (converge) then
   
               i = i + 1
               call update_history(s,omega,cell_param,domain)
   
               call cpu_time(scf_time)
               scf_time = scf_time - start_time
  
               ! Output out, rho, and omega files if s is an integer 
               if ( abs(s-float(nint(s))) < 0.001_long ) then
                  j = nint(s)
                  ! Appending 'j.' to output_prefix in file names
                  call output_summary( &
                      trim(output_prefix)//trim(int_string(j))//'.' )
                  call output_fields( &
                      trim(output_prefix)//trim(int_string(j))//'.' )
               end if
               write(6,*)
   
            else if (step_unit > ( (1.0/16.0) + 0.001 ) ) then
   
               ! Backtrack to previous chemistry
               s = s - step
               call increment_parameters(-step,domain,cell_param)
    
               ! Halve step size
               step_unit = 0.5*step_unit
    
               write(6, FMT="( / 'Backtrack and halve step_unit' / )" )
               cycle
     
            else ! If not.converge and step_unit <= 1/16, stop.
     
               write(6,*)
               write(6,*) 'Failed to converge - stop program'
               stop
   
            end if
   
            if ( (s + 0.001) >= s_max ) exit sweep_loop
   
         end do sweep_loop 
         ! end loop over sweep through parameters
 
      !# ifdef COMMENTS
      !case ('REPRESENTATION')
      !   call group_rep('rep.input')
      !# endif
 
      case ('RESPONSE')

         ! Calculate and diagonalize linear SCF response functions for
         ! a periodic microstructure. See response_mod.

         ! Check preconditions (Must iterate to convergence first)
         if (.not. iterate_flag) then
            write(6,*) "Error: Must ITERATE before calculating RESPONSE"
            exit op_loop
         end if

         ! Deferred output from ITERATE, if necessary
         if (output_flag) then
            call output_summary(output_prefix)
            call output_fields(output_prefix)
            output_flag = .FALSE.
         endif

         call response_startup(Ngrid, chain_step, order=1)
         call response_sweep(Ngrid, output_prefix)
 
         call release_basis()
         call make_basis(R_basis,G_basis,group_name,Ngrid,grid_flag=.TRUE.)

      case ('OUTPUT_GROUP')

         ! Check preconditions (Needs group created in BASIS block)
         if (.not. basis_flag) then
            write(6,*) "Error: Must make BASIS before OUTPUT GROUP"
            exit op_loop
         end if
         open(unit=field_unit, &
              file=trim(output_prefix)//'group',status='replace')
         call output_group(group,field_unit)
         close(field_unit)

      case ('OUTPUT_WAVES')

         ! Check preconditions (Needs BASIS)
         if (.not. basis_flag) then
            write(6,*) "Error: Must make BASIS before OUTPUT WAVES"
            exit op_loop
         end if

         open(unit=field_unit,file=trim(output_prefix)//'waves',&
              status='replace')
         call output_waves(field_unit, group_name)
         close(field_unit)

      case ('FIELD_TO_GRID')

         ! Read representation of field as list of coefficients of
         ! symmetrized basis functions, output file containing field
         ! values at grid points.

         ! Check preconditions for FIELD_TO_GRID
         if ( .not. unit_cell_flag ) then
            write(6,*) "Error: Must make UNIT_CELL before FIELD_TO_GRID"
            exit op_loop
         else if (.not.discretize_flag) then
            write(6,*) &
                  "Error: Must make DISCRETIZATION before FIELD_TO_GRID"
            exit op_loop
         else if (.not.basis_flag) then
            write(6,*) "Error: Must make BASIS before FIELD_TO_GRID"
            exit op_loop
         end if
 
         ! Deferred output from ITERATE, if necessary
         if (output_flag) then
            call output_summary(output_prefix)
            call output_fields(output_prefix)
            output_flag = .FALSE.
         endif

         ! Read input and output file names from input script
         call input(input_file,'input_filename')
         call input(output_file,'output_filename')
        
         ! Read field (coefficients of basis functions) from input_file
         open(unit=field_unit,file=trim(input_prefix)//trim(input_file),status='old')
         call input_field(rho,field_unit)
         close(field_unit)
 
         ! Write values of field on a grid to output_file 
         open(unit=field_unit,file=trim(output_prefix)//trim(output_file),status='replace')
         call output_field_grid(rho,field_unit,ngrid)
         close(field_unit)
 
      case ('FINISH')

         ! Deferred output of ITERATE, if necessary
         if (output_flag) then
            call output_summary(output_prefix)
            call output_fields(output_prefix)
            output_flag = .FALSE.
         endif

         ! Stop execution
         exit op_loop
 
      case default
    
         write(6,*) 'Error: Invalid op_string'
         exit op_loop
         
      end select

   end do op_loop


contains ! internal subroutines of program pscf_pd


   subroutine allocate_scf_arrays
   !-----------------------------------------------------------------
   ! Allocate arrays omega, rho, stress
   !-----------------------------------------------------------------
   integer                 :: i        ! error for file opening
   allocate(omega(N_monomer,N_star), stat = i)
   if (i.ne.0 ) stop 'Error allocating omega'
   allocate(rho(N_monomer,N_star), stat = i)
   if (i.ne.0) stop 'Error allocating rho'
   allocate(stress(N_cell_param), stat = i )
   if (i.ne.0) stop 'Error allocating stress'
   !# ifdef DEVEL
   allocate(overlap(N_monomer,N_monomer),stat = i )
   if (i.ne.0) stop 'Error allocating overlap'
   !# endif
   end subroutine allocate_scf_arrays
   !=================================================================


   subroutine output_summary(prefix)
   !-----------------------------------------------------------------
   ! Writes the output summary to the file output_prefix//suffix
   !-----------------------------------------------------------------
   use chemistry_mod, only : output_chemistry, N_chain, N_solvent, &
                             ensemble,phi_chain,phi_solvent,mu_chain, &
                             mu_solvent, interaction_type, chi
   use scf_mod, only       : free_energy_FH
   character(*) :: prefix

   real(long)   :: f_homo  ! FH free energy, kT / monomer

   open(file=trim(prefix)//'out',unit=out_unit,status='replace')
   call set_io_units(o=out_unit)

   ! Set version number for *.out file
   version%major = 1
   version%minor = 0
   call output_version(version, out_unit)

   ! Output in format of input driver file 
   if (monomer_flag) then
      write(out_unit,*)
      call output('MONOMERS',f='N',j='L')
      call output_monomers(out_unit,'F')
   end if
   if (chain_flag) then
      write(out_unit,*)
      call output('CHAINS',f='N',j='L')
      call output_chains(out_unit,'F')
   end if
   if (solvent_flag) then
      write(out_unit,*)
      call output('SOLVENTS',f='N',j='L')
      call output_solvents(out_unit,'F')
   end if
   if (composition_flag) then
      write(out_unit,*)
      call output('COMPOSITION',f='N',j='L')
      call output_composition(out_unit,'F')
   end if
   if (interaction_flag) then
      write(out_unit,*)
      call output('INTERACTION',f='N',j='L')
      call output_interaction(out_unit,'F')
   end if
   if (unit_cell_flag) then
      write(out_unit,*)
      call output('UNIT_CELL',f='N',j='L')
      call output_unit_cell(out_unit,'F')
   end if
   if (discretize_flag) then
      write(out_unit,*)
      call output('DISCRETIZATION',f='N',j='L')
      call output(ngrid,dim,'ngrid')
      call output(chain_step,'chain_step')
      !call output(extr_order,'extr_order')
   end if
   if (prefix_flag) then
      write(out_unit,*)
      call output('FILE_PREFIXES',f='N',j='L')
      call output(trim(input_prefix), 'input_prefix')
      call output(trim(output_prefix),'output_prefix')
   end if
   if (basis_flag) then
      write(out_unit,*)
      call output('BASIS',f='N',j='L')
      call output(trim(group_name),'group_name')
   end if
   if (iterate_flag) then
      write(out_unit,*)
      call output('ITERATE',f='N',j='L')
      call output_iterate_param
   end if
   if (sweep_flag) then
      write(out_unit,*)
      call output('SWEEP',f='N',j='L')
      call output(s_max,'s_max')
      call output_increments(out_unit,'N',domain)
   end if
   write(out_unit,*)
   call output('FINISH',f='N',j='L')

   ! End output of input script, begin additional information

   ! Thermodynamics 
   write(out_unit,*)
   call output('THERMO',f='N',j='L')
   call output(f_Helmholtz,'f_Helmholtz')
   f_homo = free_energy_FH(phi_chain,phi_solvent)
   call output(f_homo,     'f_homo')
   call output(pressure,   'pressure')
   select case(ensemble)
   case(0) ! phi was output by output chemistry, so
      if( N_chain > 0 ) then
         call output(mu_chain,N_chain,'mu_chain',s='C')
      end if
      if( N_solvent > 0 ) then
         call output(mu_solvent,N_solvent,'mu_solvent',s='C')
      end if
   case(1) ! mu was output chemistry, so
      if( N_chain > 0) then
         call output(phi_chain,N_chain,'phi_chain',s='C')
      end if
      if ( N_solvent > 0) then
         call output(phi_solvent,N_solvent,'phi_solvent',s='C')
      end if 
   end select
   call output(stress,N_cell_param,'stress')

   ! Output chi if chi is input as chi = chi_A/T + B
   if (interaction_type =='chi_T') then
      call output(chi,N_monomer,N_monomer,'chi',s='L')
   end if
 
   !# ifdef DEVEL
   ! Decomposition of free energy
   write(out_unit,*)
   call output('DECOMPOSE',f='N',j='L')
   call output(overlap(1,2),'overlap_AB')
   if ( N_monomer == 3 ) then
      call output(overlap(2,3),'overlap_BC')
      call output(overlap(3,1),'overlap_CA')
   end if
   call output(f_component(1),'f_enthalpy')
   call output(f_component(2),'f_head')
   call output(f_component(3),'f_tail')
   call output(f_component(4),'f_excess')
   !# endif

   ! Timing and resource statistics
   write(out_unit,*)
   call output('STATISTICS',f='N',j='L')
   call output(N_star,'N_star')
   call output(error,'Final Error') 
   call output(itr,'Iterations') 
   call output(basis_time,'Basis Time') 
   call output(scf_time,'SCF Time') 

   close(out_unit)             ! close output_prefix//'out' file
   call set_io_units(o=6)      ! reset default echo unit to stdout

   end subroutine output_summary
   !============================================================


   subroutine output_fields(prefix)
   !-----------------------------------------------------------------
   ! Writes the output summary to the file output_prefix//suffix
   !-----------------------------------------------------------------
   character(*) :: prefix

   open(file=trim(prefix)//'omega', unit=field_unit,status='replace')
   call output_field(omega,field_unit,group_name)
   close(field_unit)
 
   open(file=trim(prefix)//'rho',unit=field_unit,status='replace')
   call output_field(rho,field_unit,group_name)
   close(field_unit)

   end subroutine output_fields
   !============================================================

end program pscf_pd
!***
