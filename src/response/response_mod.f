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
!****m scf/response_mod
! MODULE
!   response_mod
! PURPOSE
!   Calculate the linear self-consistent field or random phase 
!   approximation (RPA) for the linear susceptibility of a periodic 
!   phase. As in the RPA for a homogeneous liquid, the calculation
!   requires a calculation of the linear response function for a 
!   gas of non-interacting chains, from which we then calculate the 
!   corresponding response function for an incompressible liquid. 
!
!   Each calculation is carried out for a wavevector k in the first 
!   Brillouin zone (FBZ). For each k, we calculate a matrix response 
!   in which each row describes a response of the monomer concentrations 
!   of the Bloch form exp{ik.r}d_rho(r) to a specific perturbation of
!   the form exp{ik.r} d_omega(r), where d_rho(r) and d_omega(r) have 
!   the periodicity of the unperturbed lattice. The perturbations and 
!   response matrices are represented in a basis of symmetry adapted
!   basis functions, which is constructed from the space group labeled
!   by varial k_group.
!
!   Symmetry adapted basis functions, if used, are superpositions of 
!   reciprocal lattice vectors that are related by elements of a user
!   specified group. This group (which is specified by module variable 
!   k_group) must be a subgroup of the "little group" associated with
!   wavevector k.  The little group of k is the subgroup of elements 
!   of the space group of the unperturbed crystal that leaves a 
!   function e{ik.r} invariant. In the current implementation, each 
!   basis function must be invariant under the action of all of the 
!   elements of group k_group (i.e., must correspond to an irreducible 
!   representation of k_group in which each element is represented 
!   by the identity). 
!  
! COMMENT
!   The current implementation is limited to systems containing only 
!   two types of monomer (N_monomer=2), and unperturbed structures 
!   with inversion symmetry. This allows investigation of the stability
!   of all of the known equilibrium phases of diblock copolymer melts.
!
!   For unperturbed crystals with inversion symmetry, it may be shown
!   by considering the perturbation theory used to calculate the ideal
!   gas response that the elements of the response matrix S(G,G';k) for
!   an ideal gas are real in a basis of reciprocal vectors. The same
!   statement then also applies to the corresponding response matrix 
!   for an incompressible liquid. The response matrix remains real for 
!   any representation of the response in terms of basis functions that
!   are superpositions of plane waves in which the coefficient of each
!   plane wave is real. The basis functions used in this calculation 
!   have this property. The current implementation of the code is
!   restricted to unperturbed crystals with inversion symmetry by the
!   fact that, to save memory, the response functions Smm, Spm, Spp 
!   have been declared to be real, rather than complex arrays. The
!   generalization would be straightforward, and would approximately
!   double the amount of memory required by the code. 
!
! SOURCE
!-----------------------------------------------------------------------
module response_mod
   use const_mod
   use io_mod
   use field_io_mod
   use string_mod
   use chain_mod
   use chemistry_mod
   use fft_mod
   use extrapolate_mod
   use response_step_mod
   use group_mod,     only  : operator(.dot.)
   use grid_mod,      only  : ngrid, ksq_grid, omega_grid, rho_grid
   use grid_basis_mod 
   use unit_cell_mod, only  : G_basis,R_basis
   use basis_mod,     only  : which_wave, sort, N_wave, N_star, wave, &
                              make_basis, release_basis
   !use group_rep_mod
   implicit none

   private
   public  :: k_group
   public  :: response_startup
   public  :: response_sweep
   !***

   ! public variables
   character(60)     :: k_group

   !****v response_mod/k_group 
   ! VARIABLE
   !   character(60) k_group - subgroup used to generate basis functions.
   !                           Group name or name of file containing group
   ! PURPOSE
   !   Subgroup used to construct symmetrized basis functions. Must be
   !   a subgroup of the little group of k. All basis functions are
   !   invariant under the chosen subgroup. (We have not implemented
   !   construction of basis functions for other types of irreducible 
   !   representation).
   !*** ---------------------------------------------------------------

   ! ------------------------------------------------------------------
   ! private module parameters initialized in response_startup
   ! ------------------------------------------------------------------
   type(fft_plan)           :: plan
   type(fft_plan)           :: planc
   integer                  :: nx, ny, nz
   integer                  :: n_points
   integer                  :: n_basis
   integer                  :: n_eigval_out
   integer                  :: n_eigvec_out
   integer                  :: n_basis_out
   integer                  :: extrap_order

   ! -----------------------------------------------------------------
   ! private module arrays allocated in response_startup
   ! -----------------------------------------------------------------
   type(chain_grid_type),allocatable   :: pchains(:)
   type(chain_grid_type),allocatable   :: chains0(:,:)
 
   real(long),           allocatable   :: rho_mon(:,:,:,:)
   complex(long),        allocatable   :: delrho(:,:,:,:)
   complex(long),        allocatable   :: field_rgrid(:,:,:)
 
   real(long),           allocatable   :: q_interm   (:,:,:,:)
   complex(long),        allocatable   :: delq_interm(:,:,:,:)
 
   complex(long),        allocatable   :: delrho_chain(:,:,:,:,:)
   real(long),           allocatable   :: rho_chain(:,:,:,:,:)
 
   integer,              allocatable   :: label_of_gridvec(:,:,:)
   integer,              allocatable   :: gof(:,:)           
   integer,              allocatable   :: kof(:,:)           
   real(long),           allocatable   :: ksq_grid_new(:,:,:)
 
   real(long),           allocatable   :: smm(:)
   real(long),           allocatable   :: spm(:,:)
   real(long),           allocatable   :: spp(:,:)
 

contains

   !-------------------------------------------------------------------
   !****p response_mod/response_startup
   ! SUBROUTINE
   !   response_startup(N_grids, ds, order)
   ! PURPOSE
   !   (1) allocate memory needed by module variables
   !   (2) release the basis of the full space group
   !   (3) generate symmetrized basis function from k_group
   ! ARGUMENTS
   !   integer N_grids(3)    -  grid dimensions
   !   real(long)      ds    -  chain step
   !   integer order         -  extrapolation order w.r.t. chain steps
   ! COMMENT
   !   character* k_group    -  subgoup used to generate basis function
   ! SOURCE
   !-------------------------------------------------------------------
   subroutine response_startup(N_grids, ds, order)
   integer,           intent(IN)      :: N_grids(:)
   real(long),        intent(IN)      :: ds
   integer,           intent(IN)      :: order
   !***
   integer                            :: error
   integer                            :: i_order, i_chain, nblk
 
   call create_fft_plan(n_grids,planc,fft_c2c=.true.)  ! For complex FFT data
   call create_fft_plan(n_grids,plan)                  ! For real FFT data
 
   nx=planc%n(1)-1
   ny=planc%n(2)-1
   nz=planc%n(3)-1
   n_points = (nx+1)*(ny+1)*(nz+1)  

   call input(k_group,'k_group')
   if ((dim == 1) .and. (trim(k_group) .eq. 'P 1')) k_group = '1'

   call release_basis()
   call make_basis(R_basis,G_basis,k_group,N_grids,grid_flag=.true.)
   n_basis = n_star
       
   extrap_order = order
 
   !  chain and monomer information
   if (.not.(allocated(chains0))) then
      allocate(chains0(order+1,n_chain),stat = error)
      if (error /= 0) stop "Error in allocating chains0"
   endif
 
   if (.not.(allocated(pchains))) then
      allocate(pchains(n_chain),stat = error)
      if (error /= 0) stop "Error in allocating pchains"
   endif
 
   if (.not.(allocated(field_rgrid))) then
      allocate(field_rgrid(0:nx,0:ny,0:nz),stat=error)
      if (error /= 0) stop "Error allocating field_rgrid"
   endif
 
   if (.not.(allocated(rho_mon))) then
      allocate(rho_mon(0:nx,0:ny,0:nz,n_monomer),stat=error)
      if (error /= 0) stop "Error allocating rho_mon"
   endif
 
   if (.not.(allocated(delrho))) then
      allocate(delrho(0:nx,0:ny,0:nz,n_monomer),stat=error)
      if (error /= 0) stop "Error allocating delrho"
   endif
 
   if (.not.(allocated(q_interm))) then
      allocate(q_interm(0:nx,0:ny,0:nz,2**extrap_order),stat=error)
      if (error /= 0) stop "Error allocating q_interm"
   endif
 
   if (.not.(allocated(delq_interm))) then
      allocate(delq_interm(0:nx,0:ny,0:nz,extrap_order+1),stat=error)
      if (error /= 0) stop "Error allocating delq_interm"
   endif
 
   if (.not.(allocated(delrho_chain))) then
      allocate(delrho_chain(0:nx,0:ny,0:nz,n_monomer,n_chain),stat=error)
      if (error /= 0) stop "Error allocating delrho_chain"
   endif
 
   if (.not.(allocated(rho_chain))) then
      allocate(rho_chain(0:nx,0:ny,0:nz,n_monomer,n_chain),stat=error)
      if (error /= 0) stop "Error allocating rho_chain"
   endif
 
 
   ! grid and wave information
   if (.not.(allocated(ksq_grid_new))) then
      allocate(ksq_grid_new(0:n_grids(1)-1,0:n_grids(2)-1,0:n_grids(3)-1),&
               stat=error) 
      if (error /= 0) stop "Error allocating ksq_grid_new"
   endif
 
   if (.not.(allocated(label_of_gridvec))) then
      allocate(label_of_gridvec(0:nx,0:ny,0:nz),stat = error)
      if (error /= 0) stop "Error in allocating label_of_gridvec"
   endif
 
   if (.not.(allocated(gof))) then
      allocate(gof(n_points,3),stat = error)
      if (error /= 0) stop "Error in allocating gof"
   endif
 
   if (.not.(allocated(kof))) allocate(kof(n_points,3),stat = error)
   if (error /= 0) stop "Error in allocating kof"    
 
   !  reduced response function
   if (.not.(allocated(smm))) then
      allocate(smm((n_basis+1)*n_basis/2),stat=error)
      if (error /= 0) stop "Error allocating smm in main"
   endif
 
   if (.not.(allocated(spm))) then
      allocate(spm(n_basis,n_basis),stat=error)
      if (error /= 0) stop "Error allocating spm in main"
   endif
 
   if (.not.(allocated(spp))) then
      allocate(spp(n_basis,n_basis),stat=error)
      if (error /= 0) stop "Error allocating spp in main"
   endif

   do i_chain = 1, n_chain
      nblk = n_block(i_chain)
 
      do i_order = 0, order
         call null_chain_grid(chains0(i_order+1,i_chain))
         call make_chain_grid(chains0(i_order+1,i_chain),plan,nblk,&
              block_length(1:nblk,i_chain),ds,allocate_q=.true.,   &
              perturb =.false.,order=i_order)
      end do

      call null_chain_grid(pchains(i_chain))
      call make_chain_grid(pchains(i_chain),plan,nblk,          &
         block_length(1:nblk,i_chain),ds,allocate_q=.true.,perturb=.true.)

   end do

   call response_step_startup(N_grids, extrap_order)
   call gridvecs_and_pwlabels()
 
   end subroutine response_startup
   !=========================================================================


   !---------------------------------------------------------------------------
   !****p response_mod/response_sweep
   ! SUBROUTINE
   !    response_sweep(output_prefix)
   ! PURPOSE
   !    sweep over a range of k values and calculate the response
   !    functions.
   ! ARGUMENTS
   !    integer N_grids(3)      - grid dimensions
   !    character* ouput_prefix - output_directory
   ! SOURCE
   !---------------------------------------------------------------------------
   subroutine response_sweep(Ngrid,output_prefix)
   implicit none
   integer,           intent(IN)  :: Ngrid(3)
   character(len=60), intent(IN)  :: output_prefix
   !***

   real(long)  :: kvec(3)   ! Perturbation wave-vector in crystallog units
   real(long)  :: kvec0(3)  ! Initial perturbation wave-vector, crys units
   real(long)  :: dkvec(3)  ! Perturbation wave-vector step in crys units
   integer     :: nkstep    ! # wavevectors in sweep
   integer     :: kdim      ! dimension of kvec. Must be dim, or dim+1

   real(long)  :: kvec_sq, ktrans_sq
   real(long)  :: pi, time1, time2
   integer     :: i, j

   character(len=60) :: file_prefix

   pi = 4.0_long * atan(1.0_long) 
  
   call input(kdim,'kdim')
   if ( (kdim /= dim) .AND. (kdim /= dim+1) ) then
       stop "Error: In response_sweep, kdim must equal dim or dim+1"
   endif
   kvec0 = 0.0_long
   dkvec = 0.0_long
   call input(kvec0(1:kdim),kdim,f='A')
   call input(nkstep,'nkstep')
   if (nkstep > 0) then
      call input(dkvec(1:kdim),kdim,f='A')
   endif
   call input(n_eigval_out,'n_eigval_out')
   call input(n_eigvec_out,'n_eigvec_out')
   call input(n_basis_out,'n_basis_out')
      
   kvec   = kvec0
      
   write(6,FMT = "( / '************************************'  )" )
   do j = 1, nkstep+1

      write(6,*)
      call output(j, 'k step      = ', f='L')
      call output(kvec, kdim, 'kvec        = ', f='L')
      ktrans_sq = 0.0

      ! Treat transverse component, if any
      if (kdim == dim + 1) then 

         do i=1, dim
            ktrans_sq = g_basis(1,i)*g_basis(1,i) 
         enddo
         ktrans_sq = ktrans_sq*kvec(dim+1)*kvec(dim+1)

      end if

      ! Ideal gas response
      call cpu_time(time1)
      call drho_domega_id(Ngrid,kvec,ktrans_sq)
      call cpu_time(time2)
      call output(time2-time1, 'IGR time    = ', f='L')
      
      smm = -smm
      spm = -spm
      spp = -spp

      ! non-ideal gas response
      call cpu_time(time1)
      file_prefix = trim(output_prefix)//trim(int_string(j-1))
      call calc_drho_domega_full(chi, kvec, Ngrid, &
           eval_file=trim(file_prefix)//"."//trim("eval"),&
           evec_file=trim(file_prefix)//"."//trim("evec."))
      call cpu_time(time2)
      call output(time2-time1, 'NGR time    = ', f='L')

      do i=1, kdim
         kvec(i) = kvec(i) + dkvec(i)
      enddo
      
   end do

   end subroutine response_sweep
   !=========================================================================


   !-------------------------------------------------------------------
   !****ip response_mod/drho_domega_id
   ! SUBROUTINE
   !   drho_domega_id(n_grids,kvec,ktrans_sq)
   ! PURPOSE
   !   Calculation of ideal gas response function.
   !   Calls create_fft_plan, calc_ksqgrid_new, alloc_exparrays, calc_delq,
   !         calc_delrho, fftc,destroy_fftplan, and destroy_fftplanc.
   ! ARGUMENTS
   !   integer    N_grids     -  number of grid points
   !   real(long) kvec        -  perturbation wave-vector lying in the FBZ, 
   !                             in units of RL basis vectors
   !   real(long) ktrans_sq   -  square of the transverse component of the
   !                             perturbation wave-vector (for dim < 3)
   ! SOURCE
   !-----------------------------------------------------------------
   subroutine drho_domega_id(N_grids,kvec,ktrans_sq)
   use grid_mod, only: G_to_bz, Greal_to_bz
   integer,   intent(IN)      :: N_grids(3)
   real(long),intent(IN)      :: kvec(:)
   real(long),intent(IN)      :: ktrans_sq
   !***
   complex(long)  :: delrho_k(0:nx,0:ny,0:nz) 
   real(long)     :: delrho_basis(N_basis)
   real(long)     :: G(3)
   integer        :: sp_index,i,j,k,bgn,end,error
   integer        :: i1,j1,k1,i2,i3,j2,j3
   integer        :: grid_vector(3),G_grid(3)
   integer        :: Gp(3), smm_index
   integer        :: i_gp,p_mon
   integer        :: order
   integer        :: s_end

   smm = 0.0_long
   spm = 0.0_long
   spp = 0.0_long

   order = extrap_order

   !---------------------------------------------------
   ! calculate unperturbed partition functions
   !---------------------------------------------------
   call calc_q0(order)

   !---------------------------------------------------
   ! calculate perturbed ksq_grid
   !---------------------------------------------------
   call calc_ksqgrid_new(n_grids,kvec,ktrans_sq)

   !------------------------------------------------------
   ! start perturbing potential fields; loop over monomer
   !------------------------------------------------------
   do p_mon = 1, n_monomer

      !------------------------------------------------------
      ! loop over components of the potential fields
      !------------------------------------------------------
      do i_gp = 1, n_basis

            call calc_perturb_field(i_gp)

            call calc_delq(p_mon, order)
            call calc_delrho()

            do i = 1,n_monomer

               call fftc(1,planc,delrho(:,:,:,i),delrho_k)
               call kgrid_to_basis(delrho_k, delrho_basis, no_parity=.true.)

               do k1 = 1, N_basis
                  if ( k1 >= i_gp ) then
                     smm_index      = (n_basis+1)*(i_gp-1) - i_gp*(i_gp-1)/2 + (k1-i_gp+1)
                     smm(smm_index) = smm(smm_index) + (-1)**(i+p_mon)*delrho_basis(k1)
                  end if
                  spm(k1,i_gp) = spm(k1,i_gp) + (-1)**(1+p_mon)*delrho_basis(k1)
                  spp(k1,i_gp) = spp(k1,i_gp) + delrho_basis(k1)
               end do

            end do

      end do  ! potential components loop
   end do     ! monomer loop

   smm = smm/2.0_long
   spm = spm/2.0_long
   spp = spp/2.0_long    

   end subroutine drho_domega_id
   !====================================================================


   !--------------------------------------------------------------------
   !****ip response_mod/calc_drho_domega_full
   ! SUBROUTINE
   !   calc_drho_domega_full(chi, kvec, N_grids, evec_file, eval_file)
   ! PURPOSE
   !   Calculation of interacting incompressible gas response function.
   ! ARGUMENTS
   !   real chi         -  chi parameters
   !   real kvec        -  the perturbation wave-vector lying in the FBZ, 
   !                       in units of RL basis vectors
   !   integer N_grids  -  number of grid points
   ! SOURCE
   !------------------------------------------------------------------------
   subroutine calc_drho_domega_full(chi,kvec,N_grids,evec_file,eval_file)

   use grid_mod, only : G_to_fft

   real(long),intent(IN)                :: chi(:,:)
   real(long),intent(IN)                :: kvec(:)
   integer,intent(IN)                   :: N_grids(:) 
   character(len=*),intent(IN),optional :: evec_file, eval_file
   !***

   real(long),allocatable               :: resp_mod(:,:)
   real(long),allocatable               :: work(:)
   real(long)                           :: lwork(1)
   real(long)                           :: eval(n_basis)
   integer                              :: eindex(n_basis) 
   integer                              :: ipiv(n_basis), info
   integer                              :: i, j, k, ig, jg, smm_index
   integer                              :: trgr  ! truncated wave vector

   !---------------------------------------------------
   ! invert spp
   !---------------------------------------------------
   if( allocated(work) ) deallocate(work)
   call dsytrf('L',n_basis,spp,n_basis,ipiv,lwork,-1,info)
   allocate( work( max(N_basis, int(lwork(1)) ) ) )

   call dsytrf('L',n_basis,spp,n_basis,ipiv,work,int(lwork(1)),info)
   call dsytri('L',n_basis,spp,n_basis,ipiv,work(1:N_basis),info)
   if (info /= 0) then
      stop "failed to invert spp"
   end if

   do j=2,N_basis
      do i=1,j-1 
         spp(i,j) = spp(j,i)
      end do
   end do

   !---------------------------------------------------------------
   ! calculate: S(N_basis,N_basis) = smm - smp * spp_inv * spm
   !            S is stored in spp to save memory.
   !---------------------------------------------------------------
   do i = 1, N_basis
      work(1:N_basis) = spp(i,:)
      do j = 1,N_basis
         spp(i,j) = dot_product(work(1:N_basis),spm(:,j))
      end do
   end do

   do j = 1, N_basis
      work(1:N_basis) = spp(:,j)
      do i = 1,N_basis
         spp(i,j) = dot_product(work(1:N_basis),spm(:,i))
      end do
   end do

   ! spp = smm - spp
   do j = 1, N_basis
      do i = 1, N_basis
         if ( i >= j ) then
            smm_index = (N_basis+1)*(j-1) - j*(j-1)/2 + (i-j+1)
            spp(i,j)  = smm(smm_index) - spp(i,j)
         else
            spp(i,j)  = spp(j, i)
         end if
      end do
   end do

   !--------------------------------------------------------------------
   ! calculate the eigenvectors and eigenvalues of S
   !--------------------------------------------------------------------
   if( allocated(work) ) deallocate(work)
   call dsyev('V','L',N_basis,spp,N_basis,eval,lwork,-1,info )
   allocate(work( int(lwork(1)) ))
   call dsyev('V','L',N_basis,spp,N_basis,eval,work,int(lwork(1)),info )

   if ( info < 0 ) then
      write(6,*) -info,"-th argument of dsyev has an illegal value."
      stop
   else if ( info > 0 ) then
      write(6,*) "dsyev failed to converge:"
      write(6,*) info, "-th off-diagonal element of an intermediate tridiagonal"
      write(6,*) "form did not converge to zero."
      stop
   end if

   do i = 1,n_basis
      eindex(i) = i
   end do
   eval = 1.0_long / eval - chi(1,2)

   call sort(n_basis,eval,eindex)

   !--------------------------------------------------------------------
   ! divergent eigenvalue at k=0 caused by incompressibility condition
   !--------------------------------------------------------------------
   if (abs(eval(1)) > 1.0E05) then 
      n_eigval_out = min(n_eigval_out,N_basis-1)

      call output(eval(1), 'Large eval  = ', f='L')

      do i = 1, n_eigval_out
         eval(i) = eval(i+1)
         spp(:,i)  = spp(:,i+1)
         eindex(i) = eindex(i+1)
      end do
   else
      n_eigval_out = min(n_eigval_out, N_basis)
   end if

   if (present(eval_file)) then
      open(unit=35,file=eval_file,status='unknown')
      call output(kvec,3, 'q_vec   = ',o=35,f='L')
      call output(N_basis,'N_basis = ',o=35,f='L')
      write(35,*)
      do k = 1, n_eigval_out
         write(35,'(E20.9"  "I5)') eval(k), k
      end do
      close(35)
   end if

   !--------------------------------------------------------------------
   ! output the first few eigen pair; 6 is an arbitrary number
   !--------------------------------------------------------------------
   if (present(evec_file)) then

      if (n_eigvec_out > n_eigval_out) n_eigvec_out = n_eigval_out

      n_basis_out = min(n_basis_out, N_basis)
      allocate(resp_mod(n_monomer,n_basis_out), stat=info)
      if( info /= 0 ) stop "error while allocating resp_mod"

      do k = 1, n_eigvec_out

         ! Index of the corresponding eigenvlaue
         i = eindex(k) 

         open(unit=34, file=evec_file//trim(int_string(k)), status='unknown')

         if ( N_monomer == 2 ) then
            resp_mod(1,:) =  spp(1:n_basis_out,i)
            resp_mod(2,:) = -resp_mod(1,:)
         else
            stop "Response code works only in binary system." 
         end if

         call output_field(resp_mod,34,k_group,n_basis_out)

         close(34)

      end do

      if( allocated(resp_mod) ) deallocate( resp_mod )

   end if

   end subroutine calc_drho_domega_full
   !====================================================================


   !--------------------------------------------------------------------
   !****ip response_mod/calc_q0
   !
   ! SUBROUTINE
   !    calc_q0 (order)
   !
   ! PURPOSE
   !    Calculate the unperturbed forward and backward partition functions of
   !    each chain with on-the-fly strategy using Richardson extrapolation to
   !    the required order.
   !    
   ! ARGUMENTS
   !    integer order  -  order of extrapolation
   !
   ! COMMENTS
   !    The potential fields are synchronized with public module variables:
   !       grid_mod/ksq_grid   -- initiated  in "basis_mod/make_basis"
   !       grid_mod/omega_grid -- calculated in "iterate_mod/iterate_NR"
   !
   !    The value of the partition functions at the intermediate nodes for
   !    higher order chains are saved in order to calculate del_qf & del_qr
   !    of high order. On exit, the values defined on each chain is stored
   !    in the following format:
   !
   !         segment     :      --a-----------b-- 
   !         
   !         extrapolated:      --*-----------*--     pchains(:)
   !
   !         i_order  =  0      --+-----------^--     chains0(1,:)
   !                     1      --+-----^-----^--     chains0(2,:)
   !                     2      --+--^--^--^--^--     chains0(3,:)
   !                    ...
   !
   !    where ``*'' refers to the extrapolated value for segment a and b
   !    respectively; ``^'' refers to the value calculated with the pseudo
   !    -spectral propagator starting from a-th ``*''; ``+'' refers to the
   !    value propagated from the segment right before a.
   !
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine calc_q0(order)
   implicit none
  
   integer,intent(IN)          :: order
   !***
  
   real(long)                  :: ds, dsarray(order+1), qarray(order+1), extrap_swp
   integer                     :: monomer, blk_bgn, blk_end, n_seg, n_step, pace
   integer                     :: i_chain, i_blk, i_order
   integer                     :: i1, i2, i3
  
   !---------------------------------------------
   !  initialize unperturbed partion functions:
   !    chains0 stores the as-calculated values
   !    pchains stores the extrapolated values
   !---------------------------------------------
   do i_chain = 1,N_chain
      do i_order = 0, order
         blk_end = chains0(i_order+1,i_chain)%block_bgn(N_block(i_chain)+1)
         chains0(i_order+1,i_chain)%qf(:,:,:,:)       = 0.0_long
         chains0(i_order+1,i_chain)%qf(:,:,:,1)       = 1.0_long
         chains0(i_order+1,i_chain)%qr(:,:,:,:)       = 0.0_long
         chains0(i_order+1,i_chain)%qr(:,:,:,blk_end) = 1.0_long
      end do
  
      blk_end = pchains(i_chain)%block_bgn(N_block(i_chain)+1)
      pchains(i_chain)%qf(:,:,:,:)       = 0.0_long
      pchains(i_chain)%qf(:,:,:,1)       = 1.0_long
      pchains(i_chain)%qr(:,:,:,:)       = 0.0_long
      pchains(i_chain)%qr(:,:,:,blk_end) = 1.0_long
   end do
  
   !---------------------------------------------
   ! loop over chains
   !---------------------------------------------
   do i_chain = 1, N_chain
  
      !---------------------------------------------
      ! calculate qf; loop over blocks 
      !---------------------------------------------
      do i_blk = 1, N_block(i_chain)
  
         monomer = block_monomer(i_blk, i_chain)
         ds      = pchains(i_chain)%block_ds(i_blk)
  
         call calc_exp_grid(ksq_grid, ksq_grid_new, omega_grid, &
                    monomer, ds, order, dsarray, pertb=.false.)
  
         blk_bgn = pchains(i_chain)%block_bgn(i_blk)
         blk_end = pchains(i_chain)%block_bgn(i_blk+1)
  
        !---------------------------------------------
        ! loop over segments within a block
        !---------------------------------------------
         do n_seg = blk_bgn+1, blk_end
  
           !---------------------------------------------
           ! calculate qf up to the given order
           !---------------------------------------------
            do i_order = 0, order
               pace   = 2**i_order
               n_step = (n_seg-2)*pace+2
  
               call fft_step(plan, i_order,                          &
                       pchains(i_chain)%qf(:,:,:,n_seg-1),           &
                       chains0(i_order+1,i_chain)%qf(:,:,:,n_step) )
  
               do n_step = (n_seg-2)*pace + 3, (n_seg-1)*pace + 1
                  call fft_step(plan, i_order,                            &
                          chains0(i_order+1,i_chain)%qf(:,:,:,n_step-1),  &
                          chains0(i_order+1,i_chain)%qf(:,:,:,n_step) )
               end do
            end do
  
           !---------------------------------------------
           ! extrapolate qf
           !---------------------------------------------
            if ( order == 0 ) then
               pchains(i_chain)%qf(:,:,:,n_seg) = chains0(1,i_chain)%qf(:,:,:,n_seg)
            else
               do i1 = 0,nx
               do i2 = 0,ny
               do i3 = 0,nz
  
                  do i_order = 0, order
                    pace = 2**i_order 
                    qarray(i_order+1) = &
                      chains0(i_order+1, i_chain)%qf(i1,i2,i3,(n_seg-1)*pace+1)
                  end do
  
                  call extrapolate_real(order,dsarray,qarray,extrap_swp)
  
                  pchains(i_chain)%qf(i1,i2,i3,n_seg) = extrap_swp
  
               end do
               end do
               end do
            end if
  
         end do ! n_seg
  
      end do    ! i_blk
  
     !---------------------------------------------
     ! calculate qr; loop over blocks 
     !---------------------------------------------
      do i_blk = N_block(i_chain), 1, -1
  
         monomer = block_monomer(i_blk, i_chain)
         ds      = pchains(i_chain)%block_ds(i_blk)
  
         call calc_exp_grid(ksq_grid, ksq_grid_new, omega_grid, &
                    monomer, ds, order, dsarray, pertb=.false.)
  
         blk_bgn = pchains(i_chain)%block_bgn(i_blk+1)
         blk_end = pchains(i_chain)%block_bgn(i_blk)
  
        !---------------------------------------------
        ! loop over segments within a block
        !---------------------------------------------
         do n_seg = blk_bgn-1, blk_end, -1
  
           !---------------------------------------------
           ! calculate qr up to the given order
           !---------------------------------------------
            do i_order = 0, order
               pace   = 2**i_order
               n_step = n_seg*pace
  
               call fft_step(plan, i_order,                          &
                       pchains(i_chain)%qr(:,:,:,n_seg+1),           &
                       chains0(i_order+1,i_chain)%qr(:,:,:,n_step) )
  
               do n_step = n_seg*pace - 1, (n_seg-1)*pace + 1, -1
                  call fft_step(plan, i_order,                           &
                          chains0(i_order+1,i_chain)%qr(:,:,:,n_step+1), &
                          chains0(i_order+1,i_chain)%qr(:,:,:,n_step) )
               end do
            end do
  
           !---------------------------------------------
           ! extrapolate qr
           !---------------------------------------------
            if ( order == 0 ) then
               pchains(i_chain)%qr(:,:,:,n_seg) = chains0(1,i_chain)%qr(:,:,:,n_seg)
            else
               do i1 = 0,nx
               do i2 = 0,ny
               do i3 = 0,nz
  
                  do i_order = 0, order
                    pace = 2**i_order
                    qarray(i_order+1) = &
                      chains0(i_order+1, i_chain)%qr(i1,i2,i3,(n_seg-1)*pace+1)
                  end do
  
                  call extrapolate_real(order,dsarray,qarray,extrap_swp)
  
                  pchains(i_chain)%qr(i1,i2,i3,n_seg) = extrap_swp
  
               end do
               end do
               end do
            end if
  
         end do ! n_seg
  
      end do    ! i_block
  
      pchains(i_chain)%bigQ = sum( pchains(i_chain)%qr(:,:,:,1) ) / dble( n_points )

   end do       ! i_chain
  
   end subroutine calc_q0
   !=========================================================================
  
  
   !-------------------------------------------------------------------------
   !****ip response_mod/calc_delq
   ! SUBROUTINE
   !   calc_delq(p_mon, order)
   ! PURPOSE
   !   Calculate periodic part of delta_q(r,s) Bloch function for all chains.
   ! ARGUMENTS
   !   integer p_mon  -  monomer index of the perturbation
   !   integer order  -  order of extrapolation
   ! COMMENTS
   !   The subroutine contains 4 nest loops:
   !      chains: i_chain ==> 1 : N_chain
   !      blocks: i_block ==> 1 : N_block(i_chain)
   !                          use reverse order for qr
   !      steps:  n_seg   ==> blk_bgn : blk_end
   !      order:  i_order ==> 0 : order
   ! 
   !   within the inner-most loop
   !      q_interm(1:2^order) defines the inhomogeneous term
   !      containing unperturbed values and del_omega
   !
   !   Calls subroutines calc_exparrays, pspropagate, extrapolate_real,
   !   and extrapolate_complex.
   ! SOURCE
   !-----------------------------------------------------------------
   subroutine calc_delq(p_mon, order)
   implicit none
 
   integer,intent(IN)          :: p_mon
   integer,intent(IN)          :: order
   !***
 
   real(long)                  :: ds, dsarray(order+1)
   complex(long)               :: extrap_swp, qarray(order+1)
   integer                     :: monomer, blk_bgn, blk_end, n_seg, n_step, pace
   integer                     :: i_chain, i_blk, i_order
   integer                     :: i1, i2, i3, j
 
   !---------------------------------------------
   !  initialize partition function
   !---------------------------------------------
   do i_chain = 1, n_chain
      pchains(i_chain)%del_qf(:,:,:,:) = 0.0_long
      pchains(i_chain)%del_qr(:,:,:,:) = 0.0_long
   end do
 
   !---------------------------------------------
   !  loop over chains
   !---------------------------------------------
   do i_chain = 1, N_chain
 
      !---------------------------------------------
      !  calculate del_qf; loop over blocks 
      !---------------------------------------------
      do i_blk = 1, N_block(i_chain)
 
         monomer = block_monomer(i_blk, i_chain)
         ds      = pchains(i_chain)%block_ds(i_blk)
         blk_bgn = pchains(i_chain)%block_bgn(i_blk)
         blk_end = pchains(i_chain)%block_bgn(i_blk+1)
 
         !--------------------------------------------------
         ! if p_mon == monomer(head block), del_qf = 0
         !--------------------------------------------------
         if ( i_blk == 1 .and. monomer /= p_mon ) then
            pchains(i_chain)%del_qf(:,:,:,blk_bgn:blk_end) = dcmplx(0.0,0.0)
            goto 100
         end if
 
         !--------------------------------------------------
         !  evaluate the propagators; prepare for fft_step 
         !--------------------------------------------------
         call calc_exp_grid(ksq_grid, ksq_grid_new, omega_grid, &
                   monomer, ds, order, dsarray, pertb=.true.)
 
         !---------------------------------------------
         !  loop over segments within a block
         !---------------------------------------------
         do n_seg = blk_bgn+1, blk_end
 
            !---------------------------------------------
            !  calculate del_qf up to given order
            !---------------------------------------------
            do i_order = 0, order
              pace   = 2**i_order
              n_step = (n_seg-2)*pace + 1
 
              q_interm(:,:,:,1) =                                   &
                  chains0(i_order+1,i_chain)%qf(:,:,:,n_step+1)     &
                + pchains(i_chain)%qf(:,:,:,n_seg-1)               
               
              do j = 2, pace
                 q_interm(:,:,:,j) =                                   &
                     chains0(i_order+1,i_chain)%qf(:,:,:,n_step+j)     &
                   + chains0(i_order+1,i_chain)%qf(:,:,:,n_step+j-1)    
              end do
 
              q_interm(:,:,:,1:pace) = q_interm(:,:,:,1:pace) / 2.0_long
 
              call ps_propagate(pchains(i_chain)%del_qf(:,:,:,n_seg-1),  &
                                delq_interm(:,:,:,i_order+1) ,           &
                                p_mon, field_rgrid,                      &
                                q_interm(:,:,:,1:pace),                  &
                                monomer, planc,                          &
                                i_order, dsarray(i_order+1))              
 
            end do
 
            !---------------------------------------------
            ! extrapolate del_qf
            !---------------------------------------------
            if ( order == 0 ) then
              pchains(i_chain)%del_qf(:,:,:,n_seg) = delq_interm(:,:,:,1)
            else
              do i1 = 0,nx
              do i2 = 0,ny
              do i3 = 0,nz
 
                 do i_order = 0, order
                    qarray(i_order+1) = delq_interm(i1,i2,i3,i_order+1)
                 end do
 
                 call extrapolate_complex(order,dsarray,qarray,extrap_swp)
 
                 pchains(i_chain)%del_qf(i1,i2,i3,n_seg) = extrap_swp
 
              end do
              end do
              end do
            end if
 
         end do ! n_seg
 
 100     continue
 
      end do    ! i_blk
 
      !---------------------------------------------
      !  calculate del_qr; loop over blocks
      !---------------------------------------------
      do i_blk = N_block(i_chain), 1, -1
 
         monomer = block_monomer(i_blk, i_chain)
         ds      = pchains(i_chain)%block_ds(i_blk)
         blk_bgn = pchains(i_chain)%block_bgn(i_blk+1)
         blk_end = pchains(i_chain)%block_bgn(i_blk)
 
         !--------------------------------------------------
         ! if p_mon == monomer(end block), del_qr = 0
         !--------------------------------------------------
         if ( i_blk == N_block(i_chain) .and. monomer /= p_mon ) then
            pchains(i_chain)%del_qr(:,:,:,blk_end:blk_bgn) = dcmplx(0.0,0.0)
            goto 200
         end if
 
         !--------------------------------------------------
         !  evaluate the propagators; prepare for fft_step 
         !--------------------------------------------------
         call calc_exp_grid(ksq_grid, ksq_grid_new, omega_grid, &
                   monomer, ds, order, dsarray, pertb=.true.)
 
         !---------------------------------------------
         !  loop over segments within a block
         !---------------------------------------------
         do n_seg = blk_bgn-1, blk_end, -1
 
            !---------------------------------------------
            !  calculate del_qr up to given order
            !---------------------------------------------
            do i_order = 0, order
                pace   = 2**i_order
                n_step = n_seg*pace+1
   
                q_interm(:,:,:,1) =                               &
                    chains0(i_order+1,i_chain)%qr(:,:,:,n_step-1) &
                  + pchains(i_chain)%qr(:,:,:,n_seg+1)               
                 
                do j = 2, pace
                   q_interm(:,:,:,j) =                                 &
                       chains0(i_order+1,i_chain)%qr(:,:,:,n_step-j)   &
                     + chains0(i_order+1,i_chain)%qr(:,:,:,n_step-j+1)    
                end do
   
                q_interm(:,:,:,1:pace) = q_interm(:,:,:,1:pace) / 2.0_long
   
                call ps_propagate(pchains(i_chain)%del_qr(:,:,:,n_seg+1), &
                                  delq_interm(:,:,:,i_order+1),           &
                                  p_mon, field_rgrid,                     &
                                  q_interm(:,:,:,1:pace),                 &
                                  monomer, planc,                         &
                                  i_order, dsarray(i_order+1))              
            end do
 
            !---------------------------------------------
            !  extrapolate del_qr
            !---------------------------------------------
            if ( order > 0 ) then
               do i1 = 0,nx
               do i2 = 0,ny
               do i3 = 0,nz
  
                  do i_order = 0, order
                     qarray(i_order+1) = delq_interm(i1,i2,i3,i_order+1)
                  end do
  
                  call extrapolate_complex(order,dsarray,qarray,extrap_swp)
  
                  pchains(i_chain)%del_qr(i1,i2,i3,n_seg) = extrap_swp
  
               end do
               end do
               end do
            end if
 
         end do ! n_seg
 
 200     continue
 
      end do ! i_blk
 
      pchains(i_chain)%delQ = &
              sum( pchains(i_chain)%del_qr(:,:,:,1) ) / dble( n_points )
 
   end do ! i_chain
 
   end subroutine calc_delq
   !=========================================================================


   !--------------------------------------------------------------------
   !****ip response_mod/calc_delrho
   ! SUBROUTINE
   !   calc_delrho
   ! PURPOSE
   !   Calculate periodic part of delta_rho(r) Bloch function for each chain.
   ! ARGUMENTS
   !   No arguments. Uses the module variable pchains for calculation.
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine calc_delrho()
   !***

   integer  :: sp_index,iblk,mon,i,i_bgn,i_end
   integer  :: bgn,end,j,k
   integer  :: error

   delrho_chain = (0.0_long,0.0_long)
   rho_chain    = 0.0_long

   do sp_index = 1,n_chain
      do iblk = 1,n_block(sp_index)
         mon = block_monomer(iblk,sp_index)
         i_bgn = pchains(sp_index)%block_bgn(iblk)
         i_end = pchains(sp_index)%block_bgn(iblk+1)
         i = i_bgn

         delrho_chain(:,:,:,mon,sp_index) = delrho_chain(:,:,:,mon,sp_index)&
          + (pchains(sp_index)%qf(:,:,:,i)* pchains(sp_index)%del_qr(:,:,:,i)) &
          + (pchains(sp_index)%qr(:,:,:,i)*pchains(sp_index)%del_qf(:,:,:,i))

         rho_chain(:,:,:,mon,sp_index) = rho_chain(:,:,:,mon,sp_index) &
          + (pchains(sp_index)%qf(:,:,:,i)* pchains(sp_index)%qr(:,:,:,i))

         i = i_end

         delrho_chain(:,:,:,mon,sp_index) = delrho_chain(:,:,:,mon,sp_index) &
          + (pchains(sp_index)%qf(:,:,:,i)*pchains(sp_index)%del_qr(:,:,:,i)) &
          + (pchains(sp_index)%qr(:,:,:,i)*pchains(sp_index)%del_qf(:,:,:,i)) 

        rho_chain(:,:,:,mon,sp_index) = rho_chain(:,:,:,mon,sp_index) &
          + (pchains(sp_index)%qf(:,:,:,i)*pchains(sp_index)%qr(:,:,:,i))
         
         do i = i_bgn+1,i_end-1,2
           delrho_chain(:,:,:,mon,sp_index) = delrho_chain(:,:,:,mon,sp_index) &
            +((pchains(sp_index)%qf(:,:,:,i)*pchains(sp_index)%del_qr(:,:,:,i))&
            +(pchains(sp_index)%qr(:,:,:,i)*pchains(sp_index)%del_qf(:,:,:,i)))&
               * 4.0_long

            rho_chain(:,:,:,mon,sp_index) = rho_chain(:,:,:,mon,sp_index)  &
             + pchains(sp_index)%qf(:,:,:,i)*pchains(sp_index)%qr(:,:,:,i) &
               * 4.0_long
         end do
         do i = i_bgn+2,i_end-2,2
           delrho_chain(:,:,:,mon,sp_index) = delrho_chain(:,:,:,mon,sp_index)&
            +((pchains(sp_index)%qf(:,:,:,i)*pchains(sp_index)%del_qr(:,:,:,i))&
            +(pchains(sp_index)%qr(:,:,:,i)*pchains(sp_index)%del_qf(:,:,:,i)))&
            * 2.0_long

            rho_chain(:,:,:,mon,sp_index) = rho_chain(:,:,:,mon,sp_index) &
             + pchains(sp_index)%qf(:,:,:,i)*pchains(sp_index)%qr(:,:,:,i) &
             * 2.0_long
         end do
         delrho_chain(:,:,:,mon,sp_index) = delrho_chain(:,:,:,mon,sp_index) &
             * (pchains(sp_index)%block_ds(iblk))/3.0_long
         rho_chain(:,:,:,mon,sp_index) = rho_chain(:,:,:,mon,sp_index) &
             * (pchains(sp_index)%block_ds(iblk))/3.0_long
      end do
   end do
   
   delrho(:,:,:,:)  = 0.0_long
   rho_mon(:,:,:,:) = 0.0_long
   do i = 1,N_monomer
      do j = 1, n_chain
         delrho(:,:,:,i) = delrho(:,:,:,i) + delrho_chain(:,:,:,i,j)*&
              phi_chain(j)/(chain_length(j)*(pchains(j)%bigQ))

         rho_mon(:,:,:,i) = rho_mon(:,:,:,i) + rho_chain(:,:,:,i,j)*&
         phi_chain(j)/(chain_length(j)*(pchains(j)%bigQ))
      end do
   end do
   
   end subroutine calc_delrho
   !====================================================================

 
   !--------------------------------------------------------------------
   !****ip response_mod/calc_ksqgrid_new
   ! SUBROUTINE
   !     calc_ksqgrid_new(n_grids,kvec,ktrans_sq)
   ! PURPOSE
   !   Calculation of the k-square grid after perturbation
   ! ARGUMENTS
   !   integer n_grids  -  number of grid points
   !   real  kvec       -  the perturbation wave-vector lying in the FBZ, 
   !                       in units of RL basis vectors
   !   real ktrans_sq   -  square of tranverse component of perturbation
   !                       wave-vector 
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine calc_ksqgrid_new(n_grids,kvec,ktrans_sq)
   use grid_mod,only:Greal_to_Bz

   integer,intent(IN)    :: n_grids(3)
   real(long),intent(IN) :: kvec(:)
   real(long),intent(IN) :: ktrans_sq
   !***

   real(long)            :: temp_vec(3)
   real(long)            :: shifted_vec(3)
   real(long)            :: Gbz(3),Gbzl(3)
   integer               :: i1,i2,i3,k,i
   integer               :: error
   
   ksq_grid_new = 0.0_long
   shifted_vec = 0.0_long
   
   i = 0
   do i1=0,n_grids(1)-1
      temp_vec(1)=dble(i1)
      do i2=0,n_grids(2)-1
         temp_vec(2)=dble(i2)
         do i3=0,n_grids(3)-1
            temp_vec(3)=dble(i3)
            shifted_vec = temp_vec + kvec ! kvec is in units of RL basis.
            Gbz = Greal_to_bz(shifted_vec)
            do k = 1,dim
               Gbzl(k) = Gbz(:) .dot. G_basis(:,k)
            end do
            ksq_grid_new(i1,i2,i3) = (Gbzl .dot. Gbzl)
         enddo
         i = i+1
      enddo
   enddo
   
   ksq_grid_new = ksq_grid_new + ktrans_sq
   
   end subroutine calc_ksqgrid_new
   !================================================================


   !-----------------------------------------------------------------
   !****ip response_mod/gridvecs_and_pwlabels
   ! SUBROUTINE
   !    gridvecs_and_pwlabels()
   ! PURPOSE
   !    label the waves and shifted (to FBZ) ones
   ! SOURCE
   !-----------------------------------------------------------------
   subroutine gridvecs_and_pwlabels()
   use grid_mod,only:G_to_bz
   !***
   integer :: i,i1,i2,i3,ERROR
   integer :: Gp(3)
 
   i = 1
   do i1 = 0,nx
      Gp(1) = i1
      do i2 = 0,ny
         Gp(2) = i2
         do i3 = 0,nz
            Gp(3) = i3
            gof(i,:) = G_to_bz(Gp)
            kof(i,:) = Gp(:)
            label_of_gridvec(Gp(1),Gp(2),Gp(3)) = i
            i = i+1
         end do
      end do
   end do
   end subroutine gridvecs_and_pwlabels
   ! ================================================================


   !-----------------------------------------------------------------
   !****ip response_mod/cal_perturb_field
   ! SUBROUTINE
   !    calc_perturb_field(i_basis)
   ! PURPOSE
   !    map the perturbation field onto r-grid
   ! ARGUMENTS
   !    integer i_basis  - index to the perturbing component of w
   ! SOURCE
   !-----------------------------------------------------------------
   subroutine calc_perturb_field(i_basis)
   implicit none
   integer,intent(IN)               :: i_basis
   !***
 
   real(long)                       :: delta_w(N_basis)
   complex(long)                    :: kgrid(0:nx,0:ny,0:nz)
    
   delta_w          = 0.0_long
   delta_w(i_basis) = 1.0_long
 
   call basis_to_kgrid(delta_w, kgrid, karray_full=.true., no_parity=.true.)
   call fftc(-1, planc, kgrid, field_rgrid)
 
   end subroutine calc_perturb_field
   !====================================================================
 
 
   !--------------------------------------------------------------------
   !****ip response_mod/calc_pwfield_rgrid
   ! SUBROUTINE
   !   calc_pwfield_rgrid(Gp)
   ! PURPOSE
   !   map the perturbation on to grid
   ! ARGUMENTS
   !   integer  Gp(3)  -  wave index of omega field perturbation
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine calc_pwfield_rgrid(Gp)
   integer,intent(IN)               :: Gp(3)
   !***

   integer                          :: i1,i,j,k
   real(long)                       :: G(3)
   real(long)                       :: temp_vec(3)
   real(long)                       :: coordinate(3)
   real(long)                       :: foo
   
   do i1 = 1,dim
      G(i1)       = Gp(:) .dot. g_basis(:,i1)
   end do
   
   coordinate = 0.0_long
   do i = 0,nx
      temp_vec(1) = dble(i)/dble(nx+1)
      do j = 0,ny
         temp_vec(2) = dble(j)/dble(ny+1)
         do k = 0,nz
            temp_vec(3) = dble(k)/dble(nz+1)
            do i1 = 1,dim
               coordinate(i1) = temp_vec(:) .dot. r_basis(:,i1)
            end do
            foo = G .dot. coordinate
            field_rgrid(i,j,k) = exp(dcmplx(0,foo))
         end do
      end do
   end do
    
   end subroutine calc_pwfield_rgrid
   !================================================================

end module response_mod
