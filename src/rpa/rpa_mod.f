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
!****m* scf/rpa_mod 
! MODULE 
!    rpa_mod
! PURPOSE
!    Provides procedures to calculate response/correlation functions 
!    for a homogenous polymer mixture at a specified wavenumber
! SOURCE
!-----------------------------------------------------------------------
module rpa_mod
use const_mod
use chemistry_mod
implicit none

private
public :: S_block   ! correlation function matrix for one chain
public :: S_chain   ! correlation function matrix for one chain
public :: S_ideal   ! ideal gas correlation function matrix
public :: S_inc     ! Incompressible correlation function matrix, chi=0
public :: S_rpa     ! RPA correlation function matrix, arbitrary chi
!***

contains

   !****p rpa_mod/S_block -----------------------------------
   ! FUNCTION
   !    S_block(i_chain, i_block, j_block, q)
   ! PURPOSE
   !    Return dimensionless cross correlation function for blocks
   !    i_block and j_block on a chain of species i_chain at a single
   !    wavenumber q
   !
   !    Convention: S_{ij}(k=0) = N_i N_j, where N_i is the number of i
   !    monomers in the chain of interest.
   !
   ! ARGUMENTS
   !    i_chain  - chain species index
   !    i_block  - block index
   !    j_block  - block index
   !    q        - wavenumber
   ! SOURCE
   !--------------------------------------------------------------------
   real(long) function S_block(i_chain, i_block, j_block, q)
   integer, intent(IN)     :: i_chain             ! chain species index
   integer, intent(IN)     :: i_block, j_block    ! block indices
   real(long), intent(IN)  :: q                   ! wavenumber
   !***
   real(long) :: li, lj, lk                ! block lengths
   real(long) :: xi, xj, xk                ! l*(qb)^2/6
   integer    :: i_mon, j_mon, k_mon       ! monomer indices
   integer    :: k_block, min, max         ! intermediate block indices
   real(long) :: eps = 1.0E-8

   if (i_block == j_block) then ! Self-correlation

      i_mon = block_monomer(i_block,i_chain)
      li    = block_length(i_block,i_chain)
      xi    = li * (kuhn(i_mon)*q)**2 /6.0_long
      if( xi > eps ) then
         S_block =  2.0_long * ( exp(-xi) + xi - 1.0_long )/(xi*xi)
      else
         S_block = 1.0_long
      endif
      S_block = S_block*li*li

   else ! Cross correlation

      S_block = 1.0_long

      ! Calculate factor from intermediate blocks (if any)
      if (abs(j_block-i_block) > 1) then
         if (j_block > i_block+1) then
            min = i_block + 1
            max = j_block - 1
         else if (i_block > j_block+1) then
            min = j_block + 1
            max = i_block - 1
         end if
         do k_block = min, max
            k_mon = block_monomer(k_block, i_chain)
            lk = block_length(k_block, i_chain)
            xk = lk*(kuhn(k_mon)*q)**2/6.0_long
            S_block = S_block * exp( -xk )
         end do
      endif

      ! Calculate factors from blocks i_block and j_block
      i_mon  = block_monomer(i_block, i_chain)
      j_mon  = block_monomer(j_block, i_chain)
      li = block_length(i_block, i_chain)
      lj = block_length(j_block, i_chain)
      xi = li * (kuhn(i_mon)*q)**2 / 6.0_long
      xj = lj * (kuhn(j_mon)*q)**2 / 6.0_long
      if ( xi > eps ) then 
         S_block = S_block * ( 1.0_long - exp(-xi) ) / xi
      endif
      if ( xj > eps ) then 
         S_block = S_block * ( 1.0_long - exp(-xj) ) / xj
      endif
      S_block = S_block*li*lj

   endif
   end function S_block
   !====================================================================


   !****p rpa_mod/S_chain  -----------------------------------
   ! SUBROUTINE
   !    S_chain(i_chain, q, S)
   ! PURPOSE
   !    Calculate dimensionless single-chain correlation function.
   !    Use convention such that S_{a,ij}(k=0) = N_{ai}N_{aj} 
   ! ARGUMENTS
   !    i_chain  - chain species index
   !    q        - wavenumber
   !    S(:,:)   - Correlation matrix (N_monomer x N_monomer)
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine S_chain(i_chain, q, S)
   integer, intent(IN)     :: i_chain ! chain species index
   real(long), intent(IN)  :: q       ! wavenumber
   real(long), intent(OUT) :: S(:,:)  ! correlation function
   !***
   integer :: i_block, j_block ! block indices
   integer :: i_mon,   j_mon   ! monomer type indices
   S = 0.0_long
   do i_block = 1, N_block(i_chain)
      i_mon = block_monomer(i_block, i_chain) 
      do j_block = 1, N_block(i_chain)
         j_mon = block_monomer(j_block, i_chain)
         S(i_mon, j_mon) = S(i_mon, j_mon) &
                         + S_block(i_chain, i_block, j_block, q)
      enddo
   enddo
   end subroutine S_chain
   !====================================================================


   !****p rpa_mod/S_ideal  -----------------------------------
   ! SUBROUTINE
   !    S_ideal
   ! PURPOSE
   !    Calculate correlation function for an ideal gas with compositon
   !    specified in chemistry. Includes contributions of all chain
   !    and solvent species.
   !     
   !    Convention: Fuction S_{ij}(k) is the dimensionless function:
   !
   !    S_{ij}(k) = v \int dr e^{ik.r} <\delta c_{i}(r)\delta c_{j}(0)>
   !
   !    where v is monomer reference volume.
   !
   ! ARGUMENTS
   !    q      - wavenumber
   !    S(:,:) - Correlation matrix (N_monomer x N_monomer)
   ! SOURCE
   !------------------------------------------------------------
   subroutine S_ideal(q, S, N)
   real(long), intent(IN)  :: q       ! wavenumber
   real(long), intent(OUT) :: S(:,:)  ! correlation function 
                                      ! (N_monomer,N_monomer)
   integer, intent(IN)     :: N       ! physical dimension of S
   !***
   real(long) :: Sc(N, N)
   integer    :: i_spec

   ! Chains
   S = 0.0_long
   do i_spec =1, N_chain
      call S_chain(i_spec, q, Sc)
      Sc = phi_chain(i_spec)*Sc
      Sc = Sc/chain_length(i_spec)
      S = S + Sc
   enddo

   ! Solvents
   do i_spec =1, N_solvent
      S = S + phi_solvent(i_spec)*solvent_size(i_spec)
   enddo

   end subroutine S_ideal
   !====================================================================
 
 
   !****p rpa_mod/S_inc  -------------------------------------
   ! SUBROUTINE
   !    S_inc
   ! PURPOSE
   !    Calculates RPA correlation function for an incompressible liquid
   !    with the composition specified in chemistry, but chi_{ij}=0.
   ! COMMENT
   !    Uses the convention:
   !    S_{ij}(k) = v \int dr e^{ik.r} <\delta c_{i}(r)\delta c_{j}(0)>
   !    The resulting matrix is singular, due to incompressibility.
   ! ARGUMENTS
   !    q      wavenumber
   !    S(:,:) correlation matrix (N_monomer x N_monomer)
   !    N      physical dimensions of matrix S
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine S_inc(q, S, N)
   real(long), intent(IN)  :: q       ! wavenumber
   real(long), intent(OUT) :: S(N,N)  ! correlation function matrix
                                      ! (N_monomer,N_monomer)
   integer, intent(IN)     :: N       ! physical dimension of S
   !***
   real(long) :: Sid(N,N), Sid_p(N), Sid_pp
   integer    :: i, j

   ! Calculate ideal gas correlation matrix
   call S_ideal(q,Sid,N)

   ! Calculate partial sums 
   Sid_pp = 0.0_long
   do i=1, N_monomer
      Sid_p(i) = 0.0_long
      do j=1, N_monomer
         Sid_p(i) = Sid_p(i) + Sid(i,j)
      enddo
      Sid_pp = Sid_pp + Sid_p(i)
   enddo

   ! Construct response of incompressible liquid
   do i=1, N_monomer
      do j=1, N_monomer
         S(i,j) = Sid(i,j) - Sid_p(i)*Sid_p(j)/Sid_pp
      enddo
   enddo

   end subroutine S_inc
   !====================================================================


   !****p rpa_mod/S_rpa  -------------------------------------
   ! SUBROUTINE
   !    S_rpa
   ! PURPOSE
   !    Return RPA correlation function for incompressible mixture
   !    specified in chemistry module. 
   ! COMMENT
   !    Uses the convention:
   !    S_{ij}(k) = v \int dr e^{ik.r} <\delta c_{i}(r)\delta c_{j}(0)>
   !    The resulting matrix is singular, due to incompressibility
   ! ARGUMENTS
   !    q      wavenumber
   !    S(:,:) correlation matrix (N_monomer x N_monomer)
   !    N      physical dimensions of matrix S
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine S_rpa(q, S, N)
   real(long), intent(IN)  :: q       ! wavenumber
   real(long), intent(OUT) :: S(:,:)  ! correlation function 
                                      ! (N_monomer,N_monomer)
   integer, intent(IN)     :: N       ! physical dimension of S
   !***
   real(long) :: P(N,N), Sinc(N,N)
   integer    :: i, j, k

   ! Construct response of incompressible liquid
   call S_inc(q, Sinc, N)

   ! Construct matrix P = I - Sinc*chi
   do i=1, N_monomer
      do j=1, N_monomer
         if (i == j) then
            P(i,j) = 1.0_long
         else
            P(i,j) = 0.0_long
         endif
         do k=1, N_monomer
            P(i,j) = P(i,j) + Sinc(i,k)*chi(k,j)
         enddo
      enddo
   enddo

   ! Invert P (P is replaced by its inverse upon return)
   call invert(P,N) 

   ! Evaluate matrix product S = P^{-1}*Sinc
   do i=1, N_monomer
      do j=1, N_monomer
         S(i,j) = 0.0_long 
         do k=1, N_monomer
            S(i,j) = S(i,j) + P(i,k)*Sinc(k,j) 
         enddo
      enddo
   enddo

   end subroutine S_rpa
   !====================================================================


   !--------------------------------------------------------------------
   ! Private routine for N x N matrix inversion. 
   ! Matrix A is overwritten by its inverse on output
   ! Analytic expression is used for N=2.
   !--------------------------------------------------------------------
   subroutine invert(A, N)
   real(long), intent(INOUT) :: A(:,:)  ! Matrix (N_monomer,N_monomer)
   integer, intent(IN)       :: N       ! physical dimension of S

   real(long)                :: B(N,N)  ! Work space
   real(long)                :: det
   integer                   :: i, j

   ! Lapack work space
   real(long)                :: work(8*N)
   integer                   :: ipvt(N), lwork, info
   lwork = 8*N

   if (N_monomer == 2) then ! do analytically
      det = A(1,1)*A(2,2) - A(1,2)*A(2,1)
      B(1,1) = A(1,1)
      B(1,2) = A(1,2)
      B(2,1) = A(2,1)
      B(2,2) = A(2,2)
      A(1,1) =  B(2,2)/det
      A(2,2) =  B(1,1)/det
      A(1,2) = -B(1,2)/det
      A(2,1) = -B(2,1)/det
   else ! Use LAPACK
      call dgetrf(N_monomer,N_monomer,A,N_monomer,ipvt,info)
      if (info/=0) stop "QR factorization failed in rpa_mod::invert"
      call dgetri(N_monomer,A,N_monomer,ipvt,work,lwork,info)
      if (info/=0) stop "Inversion failed in rpa_mod::invert"
   endif

   end subroutine invert
   !====================================================================

end module rpa_mod
