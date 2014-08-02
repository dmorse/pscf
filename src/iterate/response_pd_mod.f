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
!****m scf/response_pd_mod
! MODULE
!   response_pd_mod
! PURPOSE
!   Calculate response of periodic structure to periodic perturbation,
!   for use in construction of approximate Jacobian in iterate_mod.
!   Calculate approximate functional derivatives of monomer concentration 
!   field rho with respect to changes in the periodic potential field 
!   omega, i.e.,  derivatives of the coefficients in an expansion of 
!   rho with respect to coefficients in the corresponding expansion of 
!   omega. Also calculate derivatives of rho with respect to changes 
!   in unit cell dimensions. These derivatives can be to construct an 
!   approximate Jacobian for use in quasi-Newton iteration schemes. 
! SOURCE
!-----------------------------------------------------------------------
module response_pd_mod
   use const_mod
   use io_mod
   implicit none

   PRIVATE
   PUBLIC :: init_response_pd  ! allocate memory to array corrlt
   PUBLIC :: response_pd_omega ! calculate d_rho/d_omega
   PUBLIC :: response_pd_cell  ! calculate d_rho/d_cell_param
   PUBLIC :: make_correlation  ! block-block correlation of ideal chains
   PUBLIC :: corrlt            ! used by spinodal_mod

   real(long), allocatable :: corrlt(:,:,:)   ! correlation function
   !***

contains

   !--------------------------------------------------------------------
   !****p response_pd_mod/init_response_pd
   ! SUBROUTINE
   !    subroutine init_response_pd(N_monomer,N)
   ! PURPOSE
   !    Allocate memory for module array corrlt. This routine must be 
   !    called before make_correlation is called
   ! ARGUMENTS
   !    N_monomer = number monomers
   !    N         = number of basis functions
   ! SOURCE
   !-------------------------------------------------------------------
   subroutine init_response_pd(N_monomer,N)
   integer :: N_monomer, N
   !***   
   integer :: error

   if(.not.allocated(corrlt)) then
      allocate(corrlt(N_monomer,N_monomer,N),STAT=error)
   endif
   if(error /= 0) stop "Error allocating corrlt"

   end subroutine init_response_pd
   !=================================================================


   !-------------------------------------------------------------------
   !****p response_pd_mod/response_pd_omega
   ! SUBROUTINE
   !    response_pd_omega( N, cut, omega, drho_domega )
   ! PURPOSE
   !    Calculates matrix representation of response function d rho/d omega 
   ! ARGUMENTS
   !    N           - total number of symmetry-adapted basis functions
   !    cut         - number of basis functions treated numerically
   !    omega       - cofficients of monomer chemical potential field
   !    drho_domega - response matrix
   ! 
   !    drho_domega(alpha,i,beta,j) = d rho(alpha,i) / d omega(beta,j)
   ! COMMENT
   !    Long wavelength entries, with i,j in [1, cut], are calculated 
   !    numerically, by numerical differentiation. 
   !
   !    Short wavelength entries, i,j in (cut, N], are approximated by the 
   !    response of a homogeneous ideal gas, which is diagonal in k-space.
   ! SOURCE
   !-------------------------------------------------------------------
   subroutine response_pd_omega( N, cut, omega, drho_domega )
   use chemistry_mod, only : ensemble, N_monomer
   use unit_cell_mod, only : N_cell_param
   use scf_mod,       only : density

   integer, intent(IN)       :: N
   integer, intent(IN)       :: cut
   real(long), intent(IN)    :: omega(:,:)           !(N_monomer, N)
   real(long), intent(OUT)   :: drho_domega(:,:,:,:) !(N_monomer,N,N_monomer,N)
   !***

   ! Local Variables
   real(Long), dimension(N_monomer,N) :: w, p, rho
   real(long), parameter              :: increment = 1.0d-8
   real(long), parameter              :: x = 1.0_long / increment
   integer        :: ncut
   integer        :: i, alpha, j, beta
   
   ! Unperturbed density
   call density(N, omega, rho)
   drho_domega = 0.0_long
   

   if ( cut > 1 ) then
      do j = 2 - ensemble, cut
         do beta = 1, N_monomer   ! omega field loop
            w = omega
            w(beta,j) = w(beta,j) + increment
            call density(N, w, p)
 
            i = 2 - ensemble
            drho_domega(:,i:,beta,j) = x * ( p(:,i:) - rho(:,i:) )
         end do       
      end do         

      if ( cut < N ) then   ! fill the symmetric part
      do j = cut + 1, N
         do beta = 1, N_monomer
            do i = 2 - ensemble, cut
               do alpha = 1, N_monomer
                  drho_domega(alpha,i,beta,j) = drho_domega(beta,j,alpha,i)
               end do
            end do
         end do
      end do
      end if
      ncut = cut + 1
   else
     ncut = 2 - ensemble
   end if
   
   if ( ncut <= N ) then
     ! call modified_debye(N,ncut,rho,drho_domega)
       call debye_rsp(N,ncut,drho_domega)
   end if
   
   end subroutine response_pd_omega
   !=================================================================


   !-----------------------------------------------------------
   !****ip response_pd_mod/debye_rsp
   ! SUBROUTINE
   !    debye_rsp(N,ncut,drho_domega)
   ! PURPOSE
   !    Use debye function to approximate high frequency response
   !    The response matrix is block-diagonalized
   ! SOURCE
   !-----------------------------------------------------------
   subroutine debye_rsp(N,ncut,drho_domega)
   use chemistry_mod, only : N_monomer
   integer, intent(IN)       :: N
   integer, intent(IN)       :: ncut
   real(long), intent(INOUT) :: drho_domega(:,:,:,:)
   !***

   integer    :: beta, j

   do j = ncut, N
      do beta = 1, N_monomer
         drho_domega(:,j,beta,j) = - corrlt(:,beta,j)
      end do
   end do

   !print *, corrlt(:,1,:)
   !print *, corrlt(:,2,:)
   !print *, corrlt(:,3,:)
   !stop

   end subroutine debye_rsp
   !==============================================================


   !-----------------------------------------------------------
   !****ip response_pd_mod/modified_debye
   ! SUBROUTINE
   !    modified_debye(N,ncut,rho,drho_domega)
   ! PURPOSE
   !    Use density modified debye response function
   !    to approximate hight frequency response
   ! SOURCE
   !-----------------------------------------------------------
   subroutine modified_debye(N,ncut,rho,drho_domega)
   use chemistry_mod, only : N_monomer
   use unit_cell_mod, only : N_cell_param
   use SCF_mod,       only : plan
   use fft_mod
   use grid_basis_mod
   integer, intent(IN)       :: N
   integer, intent(IN)       :: ncut
   real(long), intent(IN)    :: rho(:,:)
   real(long), intent(INOUT) :: drho_domega(:,:,:,:)
   !***

   real(long)    :: rgrid(:,:,:), r_rho(:,:,:,:)
   complex(long) :: kgrid(:,:,:)
   allocatable   :: rgrid, r_rho, kgrid

   real(long)    :: kvec(N),drho(N)
   real(long)    :: rnodes
   integer       :: alpha, beta
   integer       :: i,j,info

   allocate( rgrid(0:plan%n(1)-1, 0:plan%n(2)-1, 0:plan%n(3)-1), stat=info )
   if( info /= 0 ) stop "modified_debye/rgrid(:,:,:) allocation error"

   allocate( kgrid(0:plan%n(1)/2, 0:plan%n(2)-1, 0:plan%n(3)-1), stat=info )
   if( info /= 0 ) stop "modified_debye/kgrid(:,:,:) allocation error"

   allocate( r_rho(0:plan%n(1)-1, 0:plan%n(2)-1, 0:plan%n(3)-1, N_monomer), stat=info )
   if( info /= 0 ) stop "modified_debye/r_rho(:,:,:,:) allocation error"

   rnodes = dble( plan%n(1) * plan%n(2) * plan%n(3) )

   do alpha = 1, N_monomer
      call basis_to_kgrid(rho(alpha,:), kgrid)
      call ifft(plan, kgrid, r_rho(:,:,:,alpha))
   end do

   do j = ncut, N
   do beta  = 1, N_monomer
      kvec = 0.0_long
      kvec(j) = 1.0_long
      call basis_to_kgrid(kvec, kgrid)
      call ifft(plan, kgrid, rgrid) 

      rgrid = rgrid * sqrt( r_rho(:,:,:,beta) / rho(beta,1) )

      call fft(plan, rgrid, kgrid)
      kgrid = kgrid / rnodes
      call kgrid_to_basis(kgrid, kvec)

      do alpha = 1, N_monomer
         drho = - corrlt(alpha,beta,:) * kvec
         call basis_to_kgrid(drho, kgrid)
         call ifft(plan, kgrid, rgrid)

         rgrid = rgrid * sqrt( r_rho(:,:,:,alpha) / rho(alpha,1) )

         call fft(plan, rgrid, kgrid)
         kgrid = kgrid / rnodes
         call kgrid_to_basis(kgrid, drho)

         drho_domega(alpha,ncut:N,beta,j) = drho(ncut:N)
      end do
   end do
   end do

   if( allocated(rgrid) ) deallocate( rgrid )
   if( allocated(kgrid) ) deallocate( kgrid )
   if( allocated(r_rho) ) deallocate( r_rho )

   end subroutine modified_debye
   !==============================================================


   !--------------------------------------------------------------
   !****p response_pd_mod/response_pd_cell
   ! SUBROUTINE
   !    response_pd_cell(N,omega,drho_dcell,dstress_dcell)
   ! PURPOSE
   !    drho_dcell(alpha,i,j) = d rho(alpha,i) / d cell_param(j)
   !    dstress_dcell(i,j)    = d stress(i) / d cell_param(j)
   ! SOURCE
   !--------------------------------------------------------------
   subroutine response_pd_cell(N,omega,drho_dcell,dstress_dcell)
   use chemistry_mod, only : ensemble, N_monomer
   use scf_mod,       only : density, scf_stress
   use grid_mod,      only : make_ksq
   use basis_mod,      only : make_dGsq
   use unit_cell_mod, only : N_cell_param, cell_param, &!
                             G_basis, dGG_basis, make_unit_cell

   integer, intent(IN)       :: N
   real(long), intent(IN)    :: omega(:,:)
   real(long), intent(OUT)   :: drho_dcell(:,:,:)
   real(long), intent(OUT)   :: dstress_dcell(:,:)
   !***

   ! Local variables
   real(long), parameter :: increment = 1.0d-8
   real(long), parameter :: x = 1.0_long / increment
   real(long), dimension(N_monomer,N)     :: p, rho
   real(long), dimension(N_cell_param)    :: old_param
   real(long), dimension(N_cell_param)    :: stress,new_stress
   real(long), dimension(N, N_cell_param) :: dGsq
   integer :: i, j, alpha, beta

   call make_unit_cell
   call make_ksq(G_basis)

   call density(N, omega, rho)
   do i = 1, N_cell_param
      call make_dGsq(dGsq(:,i), dGG_basis(:,:,i))
   end do
   stress = scf_stress(N, N_cell_param, dGsq )

   old_param = cell_param
   do i = 1, N_cell_param
      cell_param    = old_param
      cell_param(i) = old_param(i) + increment

      call make_unit_cell
      call make_ksq(G_basis)

      ! density response
      call density(N, omega, p)
      j = 2 - ensemble
      drho_dcell(:,j:,i) = x * ( p(:,j:) - rho(:,j:) )

      ! stress response
      do j = 1, N_cell_param
         call make_dGsq( dGsq(:,j), dGG_basis(:,:,j) )
      end do
      new_stress = scf_stress(N, N_cell_param, dGsq )
      dstress_dcell(:,i) = x * ( new_stress - stress )
   end do

   cell_param = old_param
   call make_unit_cell
   call make_ksq(G_basis)

   end subroutine response_pd_cell
   !===================================================================


   !-----------------------------------------------------------------
   !****p response_pd_mod/make_correlation
   ! SUBROUTINE
   !    make_correlation(N,[lambda])
   ! PURPOSE
   !    calculate correlation function of ideal chains melt
   ! SOURCE
   !-----------------------------------------------------------------
   subroutine make_correlation(N, lambda)
   use chemistry_mod
   use group_mod
   use basis_mod,     only : wave_of_star
   use scf_mod,     only : mu_phi, density
   use unit_cell_mod, only : G_basis
   integer        :: N           ! # of stars
   real(long), optional  :: lambda(:)   ! $k^2$
   !***

   ! local variables and indices
   integer      ::  i,  j,  k  ! block index
   real(Long), dimension(N_monomer,N) :: w, p, rho
   real(Long), dimension(N_chain) :: qout
   real(Long), dimension(N_solvent) :: q_solvent
   real(long)   :: xi, xj, xk  ! f * q * R_g**2
   real(long)   :: propg_ij    ! propogator between blk i & j
   real(long)   :: fi,fj       ! length fraction of block i & j
   integer      :: i_chain     ! chain index
   integer      :: alpha,beta  ! monomer index
   integer      :: gama        ! monomer index
   integer      :: i_star      ! basis index
   real(long)   :: qstar       ! wave vector ** 2
   real(long), parameter :: eps=1D-14
   real(long)   :: G_vec(3)

   corrlt = 0.0_long

   do i_star = 1, N
      if (present(lambda)) then
          qstar = lambda(i_star) 
      else
         G_vec = wave_of_star(:,i_star) .dot. G_basis
         qstar = G_vec .dot. G_vec
      endif
      qstar = qstar / 6.0_long

      do i_chain = 1, N_chain
         
         phi_chain(i_chain) = phi_chain(i_chain) * chain_length(i_chain)

         do i = 1, N_block(i_chain)
            alpha = block_monomer(i,i_chain)
            xi=block_length(i,i_chain) * Kuhn(alpha)**2 * qstar
            fi=block_length(i,i_chain) / chain_length(i_chain)
            corrlt(alpha,alpha,i_star)=corrlt(alpha,alpha,i_star) &!
                   + g1(xi) * fi**2 * phi_chain(i_chain)
         end do

         do i = 1, N_block(i_chain)
         do j = i+1, N_block(i_chain)
            alpha=block_monomer(i,i_chain)
            beta=block_monomer(j,i_chain)
            xi=block_length(i,i_chain) * Kuhn(alpha)**2 * qstar
            xj=block_length(j,i_chain) * Kuhn(beta)**2 * qstar
            fi=block_length(i,i_chain) / chain_length(i_chain)
            fj=block_length(j,i_chain) / chain_length(i_chain)

            propg_ij = 1.0_long
            do k = i+1, j-1
               gama=block_monomer(k,i_chain)
               xk=block_length(k,i_chain) * Kuhn(gama)**2 * qstar
               propg_ij = propg_ij * exp( -xk )
            end do

            if (alpha == beta) then
               corrlt(alpha,alpha,i_star)=corrlt(alpha,alpha,i_star)&!
                 + 2.0_long * g2(xi,xj) * propg_ij * fi * fj * phi_chain(i_chain)
            else
               corrlt(alpha,beta,i_star)=corrlt(alpha,beta,i_star)&!
                 + g2(xi,xj) * propg_ij * fi * fj * phi_chain(i_chain)
               corrlt(beta,alpha,i_star)=corrlt(beta,alpha,i_star)&!
                 + g2(xi,xj) * propg_ij * fi * fj * phi_chain(i_chain)
            end if
         end do
         end do

         if(chain_length(i_chain)>1.0D-12) then
            phi_chain(i_chain) = phi_chain(i_chain)/chain_length(i_chain)
         end if

      end do  ! loop over species
   end do     ! loop over stars

   contains

      !---------------------------------------------------------------
      function g1(x)       ! debye function, self correlation
      real(long)           :: g1,x
      g1 = 1.0_long
      if( x > eps ) g1 = g1 * 2.0_long * ( exp(-x) + x - 1.0_long ) /x/x
      return
      end function g1
      !---------------------------------------------------------------

      !---------------------------------------------------------------
      function g2(x1,x2)   ! block-pair correlation function
      real(long)           :: g2,x1,x2
      g2 = 1.0_long
      if ( x1 > eps ) g2 = g2 * ( 1.0_long - exp(-x1) ) / x1
      if ( x2 > eps ) g2 = g2 * ( 1.0_long - exp(-x2) ) / x2
      return
      end function g2
      !---------------------------------------------------------------

   end subroutine make_correlation
   !===================================================================

end module response_pd_mod
