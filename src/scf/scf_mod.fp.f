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
!****m scf/scf_mod
! PURPOSE
!   Calculate monomer concentrations, free energies, and stresses.
!   Solve the modified diffusion equation for chains by a pseudo-spectral 
!   algorithm. Use a simple Boltzmann weight on a grid for solvent.
! AUTHOR 
!   Jian Qin - Implemented pseudo-spectral algorithm (2005-2006)
!   Raghuram Thiagarajan - Added small molecule solvent (2007)
! SOURCE
!----------------------------------------------------------------------
module scf_mod 
   use const_mod
   use chemistry_mod
   use fft_mod
   use grid_mod
   use grid_basis_mod 
   use chain_mod
   use step_mod
   implicit none

   private

   ! public procedures
   public:: density_startup   ! allocates arrays needed by density
   public:: density           ! scf calculation of rho & q
   public:: scf_stress        ! calculates d(free energy)/d(cell_param)
   public:: mu_phi            ! calculates mu from phi (canonical)
                              !  or phi from mu (grand canonical)
   public:: free_energy       ! calculates helmholtz free energy 
                              ! (optionally calculates pressure)
   public:: free_energy_FH    ! Flory-Huggins helmholtz free energy
   public:: set_omega_uniform ! sets k=0 component of omega (canonical)
   !# ifdef DEVEL
   public:: divide_energy     ! calculates different components of free energy
   !# endif

   ! public module variable 
   public:: plan              ! module variable, used in iterate_mod
   !***

   type(fft_plan)                     :: plan
   type(chain_grid_type),allocatable  :: chains(:)
   integer                            :: extrap_order

   !****v scf_mod/plan -------------------------------------------------
   ! VARIABLE
   !     type(fft_plan) plan - Plan of grid sizes etc. used for FFTs
   !                           (Public because its used in iterate_mod)
   !*** ----------------------------------------------------------------

contains

   !--------------------------------------------------------------------
   !****p scf_mod/density_startup
   ! SUBROUTINE
   !    subroutine density_startup(N_grids,extr_order,chain_step,update_chain)
   !
   ! PURPOSE
   !    Initialize FFT_plan, grid_data.
   !    Allocate or update/re-allocate memory for chains
   !
   ! ARGUMENTS
   !    N_grids      = grid dimensions
   !    extr_order   = Richardson extrapolation order
   !    chain_step   = the discretized chain segment length
   !    update_chain = true if simply the chain memory need to be re-allocated
   !
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine density_startup(N_grids, extr_order, chain_step, update_chain)
   implicit none

   integer, intent(IN)    :: N_grids(3) ! # of grid points in each direction
   integer, intent(IN)    :: extr_order
   real(long), intent(IN) :: chain_step
   logical, intent(IN)    :: update_chain
   !***

   integer :: i, nblk, error
   integer :: nx,ny,nz    
   if ( .NOT. update_chain ) then
      call create_fft_plan(N_grids,plan)

      if (N_chain > 0) then
         allocate(chains(N_chain),STAT=error)
         if(error /= 0) STOP "chains allocation error in scf_mod/density_startup"
         do i=1, N_chain
           nblk      = N_block(i)
           call null_chain_grid(chains(i))
           call make_chain_grid(chains(i),plan,nblk,&
                       block_length(1:nblk,i),chain_step,allocate_q=.TRUE.)
         end do
      end if
      call init_step(N_grids)
      extrap_order = extr_order   ! set up global variable
   else
      do i=1, N_chain
        nblk = N_block(i)
        call destroy_chain_grid(chains(i))
        call make_chain_grid(chains(i),plan,nblk,&
                    block_length(1:nblk,i),chain_step,allocate_q=.TRUE.)
      end do
   end if
   
   end subroutine density_startup
   !====================================================================


   !-----------------------------------------------------------------
   !****p scf_mod/density
   ! SUBROUTINE
   !    density(N,omega,rho,qout,q_solvent)
   !
   ! PURPOSE
   !    Main SCFT calculation. Solve the modified diffusion equation
   !    for all polymer species, and calculate monomer density field
   !    for all monomer types. 
   !
   ! ARGUMENTS
   !    N                    = # of basis functions
   !    omega(N_monomer,N)   = chemical potential
   !    rho(N_monomer,N)     = monomer density fields
   !    qout(N_chain)        = partition functions
   !    q_solvent(N_solvent) = partition functions of solvent molecules
   !
   ! COMMENT   
   !      density_startup should be called prior to density to
   !      allocate arrays used by density and scf_stress.
   !
   ! SOURCE
   !------------------------------------------------------------------
   subroutine density( &
       N,              & ! # of basis functions
       omega,          & ! (N_monomer,N)chemical potential field
       rho,            & ! (N_monomer,N) monomer density field
       qout,           & ! (N_chain) 1-chain partition functions
       q_solvent       & ! (N_solvent) solvent partition functions
        )
   implicit none

   integer,    intent(IN)            :: N
   real(long), intent(IN)            :: omega(:,:)
   real(long), intent(OUT)           :: rho(:,:)
   real(long), intent(OUT), optional :: qout(N_chain)
   real(long), intent(OUT), optional :: q_solvent(N_solvent)
   !***
   
   ! local variables
!  complex(long)  :: kgrid(0:plan%n(1)/2, &
!                          0:plan%n(2)-1, &
!                          0:plan%n(3)-1)
   complex(long),allocatable  :: kgrid(:,:,:)

   real(long)     :: rnodes       ! number of grid points
   real(long)     :: partion      ! partion of single chain
   real(long)     :: bigQ_solvent ! partition of solvent
   integer        :: i_chain      ! index to chain
   integer        :: i_blk        ! index to block
   integer        :: alpha        ! index to monomer
   integer        :: i            ! dummy variable
   real(long)     :: Ns           ! number of solvent molecules in a reference volume  
   integer        :: info
   
   allocate( kgrid(0:plan%n(1)/2, 0:plan%n(2)-1, 0:plan%n(3)-1), stat=info)
   if( info /= 0 ) stop "density/kgrid(:,:,:) allocation error"
  
   rnodes=dble( plan%n(1) * plan%n(2) * plan%n(3) )

   ! Transform omega fields onto a grid
   do alpha = 1, N_monomer
     call basis_to_kgrid(omega(alpha,:),kgrid)
     call ifft(plan,kgrid,omega_grid(:,:,:,alpha))
   end do
   
   ! loop over chains
   do i_chain = 1, N_chain
     call chain_density(i_chain,chains(i_chain),ksq_grid,omega_grid)
     if(present(qout)) qout(i_chain) = chains(i_chain)%bigQ
   end do
  
   ! takes into account solvent monomer densities
   rho_grid = 0.0_long
   do i=1, N_solvent
      alpha = solvent_monomer(i)
      Ns  = solvent_size(i)        ! No. of reference volumes in a solvent molecular volume
      CALL solvent_density(alpha,Ns,omega_grid,rho_grid,bigQ_solvent)
      if(present(q_solvent)) q_solvent(i) = bigQ_solvent
   end do
   
   if ( present(qout) .AND. present(q_solvent) .AND. ensemble == 1 ) then
      call mu_phi(mu_chain,phi_chain,qout,mu_solvent,phi_solvent,q_solvent)
   end if
 
   ! calculate monomer densities
   do i_chain = 1, N_chain
     do i_blk = 1, N_block(i_chain)
         alpha = block_monomer(i_blk,i_chain)
         rho_grid(:,:,:,alpha) = rho_grid(:,:,:,alpha) &
             + phi_chain(i_chain) * chains(i_chain)%rho(:,:,:,i_blk)
     end do
   end do
  
   ! project monomer densities onto basis functions
   do alpha=1, N_monomer
     call fft(plan,rho_grid(:,:,:,alpha),kgrid)
     call kgrid_to_basis(kgrid,rho(alpha,:))
     rho(alpha,:)=rho(alpha,:)/rnodes
   end do

   if( allocated(kgrid) ) deallocate(kgrid)

   end subroutine density
   !=======================================================
  

   !--------------------------------------------------------------------------
   !****p scf_mod/solvent_density
   ! SUBROUTINE 
   !    solvent_density(monomer,s_size,omega,rho_grid,bigQ_solvent)
   !
   ! PURPOSE
   !    to calculate the density profile of a  solvent specie
   !
   ! ARGUMENTS
   !    monomer      - monomer type of the solvent species
   !    s_size       - volume occupied by solvent molecule / reference volume
   !                   (volume in units where reference volume = 1)
   !    omega        - omega fields on grid, per reference volume
   !    rho_grid     - density fields on grid    
   !    bigQ_solvent - spatial average of Boltzmann factor exp(-s_size*omega)
   !
   ! SOURCE
   !--------------------------------------------------------------------------
   subroutine solvent_density(monomer,s_size,omega,rho_grid,bigQ_solvent)
   implicit none
   
   real(long),intent(IN)              :: s_size
   real(long),intent(IN)              :: omega(0:,0:,0:,:)
   integer,intent(IN)                 :: monomer
   real(long),intent(INOUT)           :: rho_grid(0:,0:,0:,:)
   real(long),intent(OUT)             :: bigQ_solvent          
   !***
   
   real(long):: rnodes

   integer   :: ix,iy,iz,i    ! loop indices
   integer   :: solvent       ! solvent species index in phi array
   integer   :: error
 
   rnodes = dble(ngrid(1) * ngrid(2) * ngrid(3))
  
   ! calculating bigQ_solvent
   bigQ_solvent = 0.0_long  
   do iz=0, ngrid(3)-1 
      do iy=0, ngrid(2)-1
         do ix=0, ngrid(1)-1
            
            bigQ_solvent = bigQ_solvent + EXP((-s_size)&
                                              * omega(ix,iy,iz,monomer))
          
         end do
      end do
   end do     

   bigQ_solvent = bigQ_solvent/dble(rnodes)
      
   ! calculating the index of the solvent in the phi array
   do i=1, N_solvent
      if (solvent_monomer(i)==monomer) then
         solvent = i
      end if
      if ( ensemble == 1 )   phi_solvent(solvent) = bigQ_solvent*exp(mu_solvent(solvent))
   end do
 
   rho_grid(:,:,:,monomer) = rho_grid(:,:,:,monomer) + phi_solvent(solvent) * &
                             EXP((-s_size) * omega(:,:,:,monomer))/bigQ_solvent
 
   end subroutine solvent_density
   !====================================================================
  

   !--------------------------------------------------------------------
   !****p scf_mod/chain_density
   ! SUBROUTINE
   !    chain_density(i_chain, chain, ksq, omega)
   !
   ! PURPOSE
   !    solve the PDE for a single chain
   !    evaluate the density for each block
   !
   ! ARGUMENTS
   !    i_chain - index to the chain 
   !    chain   - chain_grid_type, see chain_mod
   !    ksq     - k^2 on grid, initialized in grid_mod
   !    omega   - omega fields on grid
   ! SOURCE   
   !--------------------------------------------------------------------
   subroutine chain_density(i_chain, chain, ksq, omega)
   implicit none

   integer,intent(IN)                   :: i_chain
   type(chain_grid_type),intent(INOUT)  :: chain
   real(long),intent(IN)                :: ksq(0:,0:,0:)
   real(long),intent(IN)                :: omega(0:,0:,0:,:)
   !***

   integer   :: chain_end, i_blk
   integer   :: istep, ibgn, iend
   real(long):: ds, b
   integer   :: i, j, monomer
   integer   :: ix, iy, iz

   ! Calculate qf, by integratin forward from s=0
   chain%qf(:,:,:,1) = 1.0_long
   do i_blk = 1, N_block(i_chain)
      monomer = block_monomer(i_blk, i_chain)
      ds = chain%block_ds(i_blk)
      b  = kuhn( monomer )
      call make_propg(ds, b, ksq, omega(:,:,:,monomer) )

      do istep = chain%block_bgn(i_blk), chain%block_bgn(i_blk+1)-1
         call step_exp(chain%qf(:,:,:,istep), &
                       chain%qf(:,:,:,istep+1), plan)
      end do
   end do

   ! Calculate qr, by integrating backward from s = chain_end
   chain_end = chain%block_bgn(N_block(i_chain)+1)
   chain%qr(:,:,:,chain_end) = 1.0_long
   do i_blk = N_block(i_chain), 1, -1
      monomer = block_monomer(i_blk,i_chain)
      ds = chain%block_ds(i_blk)
      b  = kuhn( monomer )
      call make_propg(ds, b, ksq, omega(:,:,:,monomer) )
      do istep = chain%block_bgn(i_blk+1), chain%block_bgn(i_blk)+1, -1
         call step_exp(chain%qr(:,:,:,istep), &
                       chain%qr(:,:,:,istep-1), plan)
      end do
   end do

   ! Calculate single chain partition function chain%bigQ
   chain%bigQ = sum(chain%qf(:,:,:,chain_end)) &
          / dble(size(chain%qf(:,:,:,chain_end)))

   ! Calculate monomer concentration fields, using Simpson's rule
   ! to evaluate the integral \int ds qr(r,s)*qf(r,s)
   chain%rho = 0.0_long
   do i = 1, N_block(i_chain)
      ! Chain ends: Add qf(r,ibgn)*qr(r,ibgn) & qf(r,iend)*qr(r,iend)
      ibgn=chain%block_bgn(i)
      iend=chain%block_bgn(i+1)

!     chain%rho(:,:,:,i)=chain%qf(:,:,:,ibgn)*chain%qr(:,:,:,ibgn)
      do iz=0,ngrid(3)-1
      do iy=0,ngrid(2)-1
      do ix=0,ngrid(1)-1
        chain%rho(ix,iy,iz,i)=chain%qf(ix,iy,iz,ibgn)*chain%qr(ix,iy,iz,ibgn)
      end do
      end do
      end do

!     chain%rho(:,:,:,i)=chain%rho(:,:,:,i)+chain%qf(:,:,:,iend)*  &
!                               chain%qr(:,:,:,iend)
      do iz=0,ngrid(3)-1
      do iy=0,ngrid(2)-1
      do ix=0,ngrid(1)-1
        chain%rho(ix,iy,iz,i)=chain%rho(ix,iy,iz,i) + &
            chain%qf(ix,iy,iz,iend)*chain%qr(ix,iy,iz,iend)
      end do
      end do
      end do

      ! Odd indices: Sum values of qf(i)*qr(i)*4.0 with i odd
      do j=ibgn+1,iend-1,2
!        chain%rho(:,:,:,i)=chain%rho(:,:,:,i)+chain%qf(:,:,:,j)*  &
!                          chain%qr(:,:,:,j)*4.0_long
         do iz=0,ngrid(3)-1
         do iy=0,ngrid(2)-1
         do ix=0,ngrid(1)-1
           chain%rho(ix,iy,iz,i)=chain%rho(ix,iy,iz,i) + &
               chain%qf(ix,iy,iz,j)*chain%qr(ix,iy,iz,j)*4.0_long
         end do
         end do
         end do

      end do

      ! Even indices: Sum values of qf(i)*qr(i)*2.0 with i even
      do j=ibgn+2,iend-2,2
!        chain%rho(:,:,:,i)=chain%rho(:,:,:,i)+chain%qf(:,:,:,j)*  &
!                            chain%qr(:,:,:,j)*2.0_long
         do iz=0,ngrid(3)-1
         do iy=0,ngrid(2)-1
         do ix=0,ngrid(1)-1
           chain%rho(ix,iy,iz,i)=chain%rho(ix,iy,iz,i) + &
               chain%qf(ix,iy,iz,j)*chain%qr(ix,iy,iz,j)*2.0_long
         end do
         end do
         end do

      end do

      ! Multiply sum by ds/3
      chain%rho(:,:,:,i)=chain%rho(:,:,:,i)*chain%block_ds(i)/3.0_long  
   end do

   chain%rho=chain%rho/chain_length(i_chain)/chain%bigQ

   end subroutine chain_density
   !====================================================================


   !--------------------------------------------------------------------
   !****p scf_mod/scf_stress
   ! FUNCTION
   !    scf_stress(N, size_dGsq, dGsq )
   !
   ! RETURN
   !    real(long) array of dimension(size_dGsq) containing
   !    derivatives of free energy with respect to size_dGsq 
   !    cell parameters or deformations
   !
   ! ARGUMENTS
   !    N         = number of basis functions
   !    size_dGsq = number of cell parameters or deformations
   !    dGsq      = derivatives of |G|^2 w.r.t. cell parameters
   !                dGsq(i,j) = d |G(i)|**2 / d cell_param(j)
   ! COMMENT
   !    Requires previous call to density, because scf_stress
   !    uses module variables computed in density.
   !
   ! SOURCE
   !--------------------------------------------------------------------
   function scf_stress(N, size_dGsq, dGsq )
   implicit none

   integer,    intent(IN) :: N
   integer,    intent(IN) :: size_dGsq
   real(long), intent(IN) :: dGsq(:,:)
   !***

   real(long)  :: scf_stress(size_dGsq)

   ! ngrid(3) was obtained by association
   ! Local Variables

   real(long)      :: dQ(size_dGsq)    ! change in q
   real(long)      :: qf_basis(N),qr_basis(N),q_swp(N)
   !complex(long)   :: kgrid(0:ngrid(1)/2,0:ngrid(2)-1,0:ngrid(3)-1)
   complex(long),allocatable   :: kgrid(:,:,:)

   real(long)      :: rnodes, normal
   real(long)      :: ds0, ds, b
   real(long)      :: increment
   integer         :: i, alpha, beta   ! summation indices
   integer         :: monomer             ! monomer index
   integer         :: sp_index            ! species index
   integer         :: ibgn,iend
   integer         :: info

   allocate( kgrid(0:ngrid(1)/2, 0:ngrid(2)-1, 0:ngrid(3)-1), stat=info )
   if ( info /= 0 ) stop "scf_mod/scf_stress/kgrid(:,:,:) allocation error"

   ! number of grid points
   rnodes = dble( ngrid(1) * ngrid(2) * ngrid(3) )

   ! normal = rnodes  * &! normalization of bigQ, divided by volume
   normal = rnodes   *  &! fft normal of forward partition
            rnodes   *  &! fft normal of backward partition
            3.0_long *  &! normal simpson's rule
            6.0_long     ! b**2/6

   scf_stress = 0.0_long

   ! Loop over chain species
   do sp_index = 1, N_chain
      dQ = 0.0_long

      ! Loop over blocks
      do alpha = 1,  N_block(sp_index) 
         monomer = block_monomer(alpha,sp_index)
               b = kuhn(monomer)
             ds0 = chains(sp_index)%block_ds(alpha)

            ibgn = chains(sp_index)%block_bgn(alpha)
            iend = chains(sp_index)%block_bgn(alpha+1)

         do i = ibgn, iend
            ! rgrid=dcmplx( chains(sp_index)%qf(:,:,:,i), 0.0_long)
            call fft(plan, chains(sp_index)%qf(:,:,:,i), kgrid )
            call kgrid_to_basis( kgrid, qf_basis )

            ! rgrid=dcmplx( chains(sp_index)%qr(:,:,:,i), 0.0_long)
            call fft(plan, chains(sp_index)%qr(:,:,:,i), kgrid )
            call kgrid_to_basis( kgrid, qr_basis )

            ds = ds0
            if ( i/= ibgn .and. i/= iend) then
               if (modulo(i,2) == 0) then
                  ds = 4.0_long * ds
               else
                  ds = 2.0_long * ds
               end if
            end if   ! Simpson's rule quadrature

            do beta = 1, size_dGsq
               q_swp     = qr_basis * dGsq(:,beta)
               increment = dot_product(q_swp, qf_basis)
               increment = increment * b**2 * ds / normal
               dQ(beta)  = dQ(beta) - increment
            end do

         end do      ! loop over nodes of single block
      end do         ! loop over blocks


      ! Note the mixing rule
      ! stress(total) = \sum_alpha \phi_alpha \cdot~stress(\alpha)
      select case(ensemble)
      case (0)
         scf_stress = scf_stress - (dQ / chains(sp_index)%bigQ)*  &
                      phi_chain(sp_index)/chain_length(sp_index)
      case (1)
         scf_stress = scf_stress - (dQ / chains(sp_index)%bigQ)*  &
                      exp(mu_chain(sp_index))*chains(sp_index)%bigQ  / &
                      chain_length(sp_index)
      end select

   end do

   if ( allocated(kgrid) ) deallocate( kgrid )

   end function scf_stress
   !===================================================================


   !# ifdef DEVEL
   !-------------------------------------------------------------------
   !****p scf_mod/divide_energy
   ! SUBROUTINE
   !    divide_energy(rho, omega, phi_chain, phi_solvent, Q, f_comp, ovlap)
   ! PURPOSE   
   !    Divide free energy into components arising from binary
   !    interaction free energy and from chain entropy
   ! ARGUMENTS
   !    rho         = density fields
   !  omega         = potential fields
   !    phi_chain  = volume fraction of species (chain)
   !    phi_solvent = volume fraction of species (solvent)
   !      Q         = partion function of species (chain)
   !  f_tot         = total free energy
   ! f_comp         = components of free energy (see below)
   !  ovlap         = overlap integrals
   ! COMMENT
   !
   !    a) Components of f_comp array:
   !       f_comp(1) = overall interaction energy
   !       f_comp(2) = conformational energy of first block
   !       f_comp(3) = conformational energy of last  block
   !       f_comp(4) = junction translational energy (diblock)
   !
   !    b) Calculation of junction translational entropy is correct
   !       only for diblocks, for which there is only one junction
   !
   !    c) Components of overlap integral array ovlap can be used
   !       to divide interaction energy into components arising from
   !       interactions between specific pairs of monomer types.
   !
   ! SOURCE
   !----------------------------------------------------------------
   subroutine divide_energy(rho, omega, phi_chain, phi_solvent, Q, f_tot, f_comp, ovlap)
   implicit none
   real(long), intent(IN)  :: rho(:,:)        ! monomer vol. frac fields
   real(long), intent(IN)  :: omega(:,:)      ! chemical potential field
   real(long), intent(IN)  :: phi_chain(:)    ! molecule vol. frac of chain mol
   real(long), intent(IN)  :: phi_solvent(:)  ! molecule vol. frac of solvent mol
   real(long), intent(IN)  :: Q(:)            ! chain partition functions
   real(long), intent(IN)  :: f_tot           ! components of free energy
   real(long), intent(OUT) :: f_comp(:)       ! components of free energy
   real(long), intent(OUT) :: ovlap(:,:)      ! overlap integrals,N_monomer**2
   !***

   real(long) :: rnodes
   real(long) :: enthalpy    ! interaction energy
   real(long) :: fhead       ! head block energy
   real(long) :: ftail       ! tail block energy
   real(long) :: fjct        ! junction translational entropy
   real(long) :: ftmp        ! junction translational entropy (temporary)
   integer    :: alpha, beta ! monomer indices
   integer    :: i, nh, nt   ! loop indices
   integer    :: ix,iy,iz    ! loop indices

   rnodes = dble( ngrid(1) * ngrid(2) * ngrid(3) )

   enthalpy = 0.0_long
   ovlap    = 0.0_long
   do alpha = 1, N_monomer
      do beta = alpha+1, N_monomer
         ovlap(alpha,beta) = dot_product(rho(alpha,:),rho(beta,:))
         ovlap(beta,alpha) = ovlap(alpha,beta)
         enthalpy = enthalpy + chi(alpha,beta) * ovlap(alpha,beta)
      end do
   end do

   fhead = 0.0_long
   ftail = 0.0_long
   do i=1, N_chain
      nh = chains(i)%block_bgn(2)
      do iz = 0, ngrid(3)-1
      do iy = 0, ngrid(2)-1
      do ix = 0, ngrid(1)-1
         if ( chains(i)%qf(ix,iy,iz,nh) > 0.0_long .AND. &
              chains(i)%qr(ix,iy,iz,nh) > 0.0_long ) then
            fhead = fhead - chains(i)%qf(ix,iy,iz,nh)   &
                          * chains(i)%qr(ix,iy,iz,nh)   &
                          / Q(i)                        &
                          * log( chains(i)%qf(ix,iy,iz,nh) ) & 
                          * phi_chain(i) / chain_length(i)
         end if
      end do
      end do
      end do

      nt = chains(i)%block_bgn(N_block(i))
      do iz = 0, ngrid(3)-1
      do iy = 0, ngrid(2)-1
      do ix = 0, ngrid(1)-1
         if ( chains(i)%qf(ix,iy,iz,nt) > 0.0_long .AND. &
              chains(i)%qr(ix,iy,iz,nt) > 0.0_long ) then
            ftail = ftail - chains(i)%qf(ix,iy,iz,nt)   &
                          * chains(i)%qr(ix,iy,iz,nt)   &
                          / Q(i)                        &
                          * log( chains(i)%qr(ix,iy,iz,nt) ) &
                          * phi_chain(i) / chain_length(i)
         end if
      end do
      end do
      end do
   end do
   fhead = fhead / rnodes
   ftail = ftail / rnodes

   ! When monomer types in the middle blocks are different from either
   ! head or tail block, the subtraction below is correct.
   do i=1, N_chain
      beta=block_monomer(1,i)
      fhead = fhead - dot_product(omega(beta,:),rho(beta,:)) * phi_chain(i)

      beta=block_monomer(N_block(i),i)
      ftail = ftail - dot_product(omega(beta,:),rho(beta,:)) * phi_chain(i)
   end do

   ! --------------------------------------------
   ! The following block was used to calculate
   ! junction entropy contribution to free
   ! energy of diblocks. 
   ! Since it is not universal, we now instead
   ! calculate by subtraction, which can be
   ! interpreted by excess entropies for arbitrary
   ! molecular (linear) architecture.
   ! --------------------------------------------
   !## fjct = 0.0_long
   !## do i=1, N_chain
   !##    nh = chains(i)%block_bgn(2)
   !##    do iz = 0, ngrid(3)-1
   !##    do iy = 0, ngrid(2)-1
   !##    do ix = 0, ngrid(1)-1
   !##       if ( chains(i)%qf(ix,iy,iz,nh) > 0.0_long .AND. &
   !##            chains(i)%qr(ix,iy,iz,nh) > 0.0_long ) then
   !##          ftmp  = chains(i)%qf(ix,iy,iz,nh) &
   !##                * chains(i)%qr(ix,iy,iz,nh) &
   !##                / Q(i) 
   !##          fjct = fjct + ftmp * log( ftmp ) * phi_chain(i) / chain_length(i)
   !##       end if
   !##    end do
   !##    end do
   !##    end do
   !## end do
   !## fjct  = fjct  / rnodes
   ! --------------------------------------------
   fjct = 0.0_long
   do i=1, N_chain
      fjct = fjct + phi_chain(i) / chain_length(i)
   end do
   fjct = f_tot + fjct - enthalpy - fhead - ftail

   f_comp(1) = enthalpy   ! overall interaction energy
   f_comp(2) = fhead      ! conformational energy of first block
   f_comp(3) = ftail      ! conformational energy of last  block
   f_comp(4) = fjct       ! junction translational energy (diblock)

   end subroutine divide_energy
   !=============================================================
   !# endif


   !------------------------------------------------------------
   !****p scf_mod/set_omega_uniform
   ! SUBROUTINE
   !    set_omega_uniform(omega)
   ! PURPOSE
   !   Sets uniform (k=0) component of field omega to convention
   !      omega(:,1) = chi(:,:) .dot. phi_mon(:)
   !   corresponding to vanishing Lagrange multiplier field
   ! SOURCE
   !------------------------------------------------------------
   subroutine set_omega_uniform(omega)
   real(long), intent(INOUT) :: omega(:,:)
   !***

   integer    :: i, j, alpha, beta  ! loop indices
   real(long) :: phi_mon(N_monomer) ! average monomer vol. frac.

   phi_mon = 0.0_long
   do i = 1, N_chain
      do j = 1, N_block(i)
         alpha = block_monomer(j,i)
         phi_mon(alpha) = phi_mon(alpha) &
                        + phi_chain(i)*block_length(j,i)/chain_length(i)
      end do
   end do
   do i = 1, N_solvent
      alpha = solvent_monomer(i)
      phi_mon(alpha) = phi_mon(alpha) + phi_solvent(i)
   end do
   do alpha = 1, N_monomer
      omega(alpha,1) = 0.0_long
      do beta = 1, N_monomer
         omega(alpha,1) = omega(alpha,1) &
              + chi(alpha,beta) * phi_mon(beta)
      end do
   end do
   end subroutine set_omega_uniform
   !================================================================


   !-------------------------------------------------------------
   !****p scf_mod/mu_phi
   ! SUBROUTINE
   !    mu_phi(mu_chain,phi_chain,q,mu_solvent,phi_solvent,q_solvent)
   ! PURPOSE
   !    If ensemble = 0 (canonical), calculate mu from phi
   !    If ensemble = 1 (grand can), calculate phi from mu
   ! ARGUMENTS
   !    mu_chain(N_chain)     = chemical potentials of chain(units kT=1)
   !    phi_chain(N_chain)    = molecular volume fractions of chain
   !    q(N_chain)             = single chain partition functions
   !    mu_solvent(N_solvent)  = chemical potentials of solvent species
   !    phi_solvent(N_solvent) = molecular volume fractions of solvent species
   !    q_solvent(N_solvent)   = partition functions of solvent species
   !
   ! SOURCE
   !-------------------------------------------------------------
   subroutine mu_phi(mu_chain,phi_chain,q,mu_solvent,phi_solvent,q_solvent)
   real(long), intent(INOUT) :: mu_chain(N_chain)
   real(long), intent(INOUT) :: phi_chain(N_chain) 
   real(long), intent(IN)    :: q(N_chain)
   real(long), intent(INOUT) :: mu_solvent(N_solvent)
   real(long), intent(INOUT) :: phi_solvent(N_solvent)
   real(long), intent(IN)    :: q_solvent(N_solvent) 
   !***

   integer :: i
   select case(ensemble)
   case (0)
      do i = 1, N_chain
         mu_chain(i) = log( phi_chain(i) / q(i) )
      end do
      do i = 1, N_solvent
         mu_solvent(i) = log( phi_solvent(i) / q_solvent(i) )  
      end do
   case (1)
      do i = 1, N_chain
         phi_chain(i) = q(i)*exp(mu_chain(i))
      end do
      do i = 1, N_solvent
         phi_solvent(i) = q_solvent(i)*exp(mu_solvent(i))
      end do
   end select
   end subroutine mu_phi
   !================================================================


   !--------------------------------------------------------------------
   !****p scf_mod/free_energy
   ! SUBROUTINE
   !    free_energy( N, rho, omega, phi_chain, mu_chain, phi_solvent,
   !                 mu_solvent, f_Helmholtz, [pressure] )
   ! PURPOSE   
   !    Calculates Helmholtz free energy / monomer and (optionally)
   !    the pressure, given phi, mu, and omega and rho fields
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine free_energy(N, rho, omega, phi_chain, mu_chain, &
                          phi_solvent, mu_solvent, f_Helmholtz, pressure )
   integer, intent(IN)    :: N              ! # of basis functions
   real(long), intent(IN) :: rho(:,:)       ! monomer vol. frac fields
   real(long), intent(IN) :: omega(:,:)     ! chemical potential field
   real(long), intent(IN) :: phi_chain(:)   ! molecule vol. frac of chain species
   real(long), intent(IN) :: mu_chain(:)    ! chemical potential of chain species
   real(long), intent(IN) :: phi_solvent(:) ! molecule vol. fraction of solvent species 
   real(long), intent(IN) :: mu_solvent(:)  ! chemical potential of solvent species
   real(long), intent(OUT):: f_Helmholtz    ! free energy/monomer
   real(long), intent(OUT), optional :: pressure 
   !***
 
   integer :: i, alpha, beta ! loop indices

   f_Helmholtz = 0.0_long
   do i = 1, N_chain
      if ( phi_chain(i) > 1.0E-8 ) then
         f_Helmholtz = f_Helmholtz &
                     + phi_chain(i)*( mu_chain(i) - 1.0_long )/chain_length(i)
      end if
   end do
   do i=1, N_solvent
      if ( phi_solvent(i) > 1.0E-8) then
         f_Helmholtz = f_Helmholtz &
                     + phi_solvent(i)*( mu_solvent(i) - 1.0_long)/solvent_size(i)
      end if
   end do
   do i = 1, N
      do alpha = 1, N_monomer
         do beta = alpha+1, N_monomer
            f_Helmholtz = f_Helmholtz &
                        + rho(alpha,i)*chi(alpha,beta)*rho(beta,i)
         end do
         f_Helmholtz = f_Helmholtz - omega(alpha,i) * rho(alpha,i)
      end do
   end do
   
   if (present(pressure)) then
      pressure = -f_Helmholtz
      do i = 1, N_chain
         pressure = pressure + mu_chain(i)*phi_chain(i)/chain_length(i)
      end do
      do i = 1, N_solvent
         pressure = pressure + mu_solvent(i)*phi_solvent(i)/solvent_size(i)
      end do
   end if
 
   end subroutine free_energy
   !====================================================================


   !--------------------------------------------------------------------
   !****p scf_mod/free_energy_FH
   ! FUNCTION
   !    real(long) function free_energy_FH(phi_chain,phi_solvent)
   ! RETURN
   !    Flory-Huggins Helmholtz free energy per monomer, in units
   !    such that kT =1, for a homogeneous mixture of the specified 
   !    composition.
   ! ARGUMENTS
   !    phi_chain(N_chain)     = molecular volume fractions of chains
   !    phi_solvent(N_solvent) = molecular volume fractions of solvents
   ! SOURCE
   !--------------------------------------------------------------------
   real(long) function free_energy_FH(phi_chain,phi_solvent)
   real(long), intent(IN)           :: phi_chain(N_chain)
   real(long), intent(IN), optional :: phi_solvent(N_solvent)

   real(long)             :: rho(N_monomer)
   !***
   integer :: i, j, i_block, i_mon
   free_energy_FH = 0.0_long
   rho = 0.0_long

   do i = 1, N_chain
      if ( phi_chain(i) > 1.0E-8 ) then
           free_energy_FH = free_energy_FH + & 
                            (phi_chain(i)/chain_length(i))*(log(phi_chain(i))-1)
      end if
      do i_block = 1, N_block(i)
         i_mon = block_monomer(i_block,i)
         rho(i_mon) = rho(i_mon) & 
                    + phi_chain(i)*block_length(i_block,i)/chain_length(i)
      end do
   end do

   if (present(phi_solvent)) then
      do i=1, N_solvent
         if ( phi_solvent(i) > 1.0E-8 ) then
              free_energy_FH = free_energy_FH + &
                         (phi_solvent(i)/solvent_size(i))*(log(phi_solvent(i))-1)
         end if
         i_mon = solvent_monomer(i)
         rho(i_mon) = rho(i_mon) + phi_solvent(i)
      end do
   end if

   do i = 1, N_monomer - 1
      do j = i+1, N_monomer
         free_energy_FH = free_energy_FH + chi(i,j)*rho(i)*rho(j)
      end do
   end do
   end function free_energy_FH
   !=============================================================

end module scf_mod
