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
!****m scf/spinodal_mod
! MODULE
!    spinodal_mod
! PURPOSE
!    homogeneous spinodal sweep with Random Phase Approximation (RPA)
! COMMENT
! AUTHOR
! SOURCE
!-----------------------------------------------------------------------
module spinodal_mod
   use const_mod
   use io_mod
   use chemistry_mod
   use response_pd_mod
   implicit none

   private

   public :: rpa_homo_startup
   public :: rpa_homo
   public :: rpa_blend
   public :: triblock_rpa_homo
   !***

   ! Private variables
   integer, allocatable    :: ipvt(:)         ! Pivots for LU decomp 
   integer                 :: lwork           ! workspace dimension of inverse
   real(long), allocatable :: work(:)         ! workspace for matrix inverse
   real(long), allocatable :: foo(:,:)        ! temporary array
   real(long), allocatable :: qstar(:)        ! q**2 for correlation calc.

contains

   !---------------------------------------------------------------------
   !****p spinodal_mod/rpa_homo_startup
   ! SUBROUTINE
   !    rpa_homo_startup
   ! PURPOSE
   !   Allocate arrays needed by rpa_homo
   ! SOURCE
   !---------------------------------------------------------------------
   subroutine rpa_homo_startup
   !***

   ! local variables
   integer    :: info

   allocate(qstar(1))
   call input(qstar(1),'initial qstar^2')

   call init_response_pd(N_monomer,1)

   allocate(ipvt(N_monomer))
   allocate(work(N_monomer))
   allocate(foo(N_monomer,N_monomer))
   ! estimate the opitimal workspace dimension
   lwork = -1
   call dgetri(N_monomer,foo,N_monomer,ipvt,work,lwork,info)
   lwork = work(1)
   if (allocated(work)) deallocate(work)
   allocate(work(lwork))

   !call triblock_bimode

   end subroutine rpa_homo_startup
   !===================================================================

   !---------------------------------------------------------------------
   !****p spinodal_mod/rpa_homo
   ! SUBROUTINE
   !    rpa_homo(N)
   ! PURPOSE
   !    Find the spinodals of homogeneous phase
   !    Two searching strategies are implemented
   !    (1) search N at fixed compositions
   !    (2) search compositions at fixed N
   ! ARGUMENTS
   !    integer i_unit = input unit
   ! SOURCE
   !---------------------------------------------------------------------
   subroutine rpa_homo(i_unit)

   integer, intent(IN)    :: i_unit  ! input unit
   !***

   ! local variables
   real(long) :: block_0(N_blk_max,N_chain)
   real(long) :: d_block(N_blk_max,N_chain)
   real(long) :: phi_0(N_chain), d_phi(N_chain)
   real(long) :: chi_0(N_monomer,N_monomer)
   real(long) :: N_srch, scale_srch, qsq_srch
   real(long) :: res_tolerance
   logical    :: search_N

   real(long),parameter       :: qsq_tol=1.0d-5
   real(long),parameter       :: n_inc=1.0d-3
   real(long)                 :: oldres, newres
   character(len = 100)       :: comment_line
   integer                    :: Nitr, itrmax
   integer                    :: i, j, imax

   real(long),dimension(N_monomer,N_monomer)     :: R0, Ri
   real(long),dimension(N_monomer-1,N_monomer-1) :: gbar

   call input(imax,'total_steps')
   read(i_unit,*) comment_line
   call output(comment_line,f='N',j='L')
   do j = 1, N_chain
      call input(d_block(:,j),N_block(j),f='N')
   enddo
   read(i_unit,*) comment_line
   call output(comment_line,f='N',j='L')
   do j = 1, N_chain
      call input(d_phi(j),f='N')
   enddo
   call input(res_tolerance, 'res_tolerance')
   call input(search_N, 'search_N')

   itrmax = 2000
   qsq_srch = qstar(1)
   if (search_N) then
      chi_0  = chi
      N_srch = 1.0_long
      open(15, file="spinodal.dat", status='unknown')

      do i = 1, imax
         block_length = block_length + d_block
         phi_chain    = phi_chain + d_phi
   
         Nitr = 0
         chi  = N_srch * chi_0

         N_loop : do
           oldres = min_eigen(qsq_srch)
   
           if ( Nitr > itrmax ) then
              write(6,*) "exceeding max search in N_loop"
              exit N_loop
           else if ( abs( oldres ) < res_tolerance ) then
              exit N_loop
           end if
   
           chi = (N_srch + N_inc) * chi_0
           newres = min_eigen(qsq_srch)
           N_srch = N_srch - oldres * N_inc / (newres - oldres)
   
           chi = N_srch * chi_0
           Nitr = Nitr + 1
         end do N_loop

         write(15,'(2f10.5,2f18.5)') block_length(2,1), block_length(1,1), N_srch*chi_0(1,2), qsq_srch
      end do
      close(15)
   else
      block_0 = block_length
      phi_0   = phi_chain

      scale_srch = 0.0_long
      Nitr       = 0
      scale_loop : do
        oldres = min_eigen(qsq_srch)

        if ( Nitr > itrmax ) then
           write(6,*) "exceeding max search in N_loop"
           exit scale_loop
        else if ( abs( oldres ) < res_tolerance ) then
           exit scale_loop
        end if

        block_length = block_0 + (scale_srch + N_inc) * d_block
        phi_chain    = phi_0   + (scale_srch + N_inc) * d_phi
        newres       = min_eigen(qsq_srch)

        scale_srch   = scale_srch - oldres * N_inc / (newres - oldres)
        block_length = block_0 + scale_srch * d_block
        phi_chain    = phi_0   + scale_srch * d_phi

        Nitr = Nitr + 1
      end do scale_loop
      write(6,*) "-----------------------------------------------"
      write(6,*) "      N_iteration = ", Nitr
      write(6,*) "          residal = ", oldres
      write(6,*) "         qsq_srch = ", qsq_srch
      write(6,*) "           phi(1) = ", phi_chain(1)
      do j=1, N_block(1)
         write(6,*) "block_length(:,1) = ", block_length(j,1)
      enddo
      write(6,*) "-----------------------------------------------"
   end if

   contains

      !----------------------------------------------------------
      function f2(x)
      real(long)            :: f2, x
      real(long),parameter  :: eps2 = 1.0D-3
      f2 = ( f1(x+eps2) -  f1(x) ) / eps2
      end function f2
      !----------------------------------------------------------

      !----------------------------------------------------------
      function f1(x)
      real(long)            :: f1, x
      real(long),parameter  :: eps1 = 1.0D-3
      f1 = ( eigen(x+eps1) - eigen(x) )/ eps1
      end function f1
      !----------------------------------------------------------

      !----------------------------------------------------------
      function min_eigen(x)
      real(long) :: min_eigen,x
      integer    :: qitr
      qitr = 0
      q_loop : do
        if ( qitr > itrmax ) then
           write(6,*) "exceeding max search in min_eigen/q_loop"
           exit q_loop
        else if ( abs( f1(x) ) < qsq_tol ) then
           exit q_loop
        end if
        x = x - f1(x)/f2(x)
        qitr = qitr + 1
      end do q_loop
      min_eigen = eigen(x)
      end function min_eigen
      !----------------------------------------------------------

      !----------------------------------------------------------
      ! PURPOSE
      !   calculate the minimum eigenvalue of structure factor
      !   for triblocks
      !     R0   = 3x3, ideal Gaussian chain correlation function
      !     Ri   = 3x3, inverse of R0
      !     gbar = reduced Ri, by using incompressibility
      !----------------------------------------------------------
      function eigen(x)
      real(long) :: eigen, x
      real(long) :: b1, b2, c

      qstar(1) = x
      call make_correlation(1,qstar)
      R0 = corrlt(:,:,1) 
      Ri = inv( R0, N_monomer )

      if (N_monomer <= 1) then
         stop "At lease 2 monomer types are needed."
      else if (N_monomer == 2) then
         b1 = Ri(1,1) + Ri(2,2) - Ri(1,2) - Ri(2,1)
         eigen = b1 - 2.0_long * chi(1,2)
      else if (N_monomer == 3) then
         gbar(1,1) = Ri(1,1) - 2.0_long*Ri(1,2) + Ri(2,2) - 2.0_long*chi(1,2)
         gbar(2,2) = Ri(3,3) - 2.0_long*Ri(3,2) + Ri(2,2) - 2.0_long*chi(2,3)
         gbar(1,2) = Ri(1,3) + Ri(2,2) - Ri(1,2) - Ri(2,3) &
                     + chi(1,3) - chi(1,2) - chi(2,3)
         gbar(2,1) = gbar(1,2)

         b1 = gbar(1,1) +  gbar(2,2)
         b2 = dsqrt( (gbar(1,1)-gbar(2,2))**2 + 4.0_long*gbar(1,2)**2 )
         if( abs(gbar(1,1)-gbar(2,2)) < 1.0D-12 ) b2 = 2.0_long*gbar(1,2)

         eigen = (b1 - b2) / 2.0_long
      else
         stop "More than 3 monomer types are not implemented."
      end if
      end function eigen
      !----------------------------------------------------------
    
      !----------------------------------------------------------
      function inv(x,xdim)
      real(long) :: x(:,:)
      integer    :: xdim
      real(long) :: inv(xdim, xdim)
      integer    :: lda, info
      lda = xdim
      ! LU factorization
      call dgetrf(xdim,xdim,x,lda,ipvt(1:xdim),info)
      if(info/=0) stop "LU factorization failed."
      ! matrix inversion
      call dgetri(xdim,x,lda,ipvt(1:xdim),work,lwork,info)
      if(info/=0) stop "Matrix inversion failed."
      inv = x
      end function inv
      !----------------------------------------------------------

   end subroutine rpa_homo
   !===================================================================



   !---------------------------------------------------------------------
   !****p spinodal_mod/rpa_homo_old
   ! SUBROUTINE
   !   rpa_homo_old(i_unit)
   ! PURPOSE
   !   ....
   ! ARGUMENTS
   !   i_unit = # of q vectors
   ! SOURCE
   !---------------------------------------------------------------------
   subroutine rpa_homo_old(i_unit)
   implicit none

   integer, intent(IN)    :: i_unit  ! input unit #
   !***

   ! local variables
   real(long)            :: block0(N_blk_max,N_chain)
   real(long)            :: dblock(N_blk_max,N_chain)
   real(long)            :: eps, sigma
   real(long),dimension(N_monomer,N_monomer)      :: R0, R1, R
   real(long),dimension(N_monomer,N_monomer)      :: E, ones
   real(long),dimension(N_monomer-1,N_monomer-1)  :: Rcore, Ri, vl, vr
   real(long),dimension(N_monomer-1)              :: wr,wi
   real(long),dimension(N_monomer,N_monomer)      :: eigv
   character(len = 100)  :: comment_line

   real(long)  :: old, new
   real(long)  :: Hessen(2,2), gradold(2), gradnew(2)
   real(long)  :: qsq, lmd       ! q^2, lambda
   real(long)  :: dqsq, dlmd     ! change in q^2, lambda
   real(long)  :: qsqold, lmdold ! q^2, lambda
   real(long)  :: qsqnew, lmdnew ! q^2, lambda
   real(long)  :: Hesseninv(2,2), relax

   real(long)  :: maxerr
   integer     :: cdim, info  ! dimension of core response matrix
   integer     :: j, k, nsch, itr

   cdim   = N_monomer - 1
   eps    = 1.0D-4
   ones   = 1.0_long
   E      = 0.0_long
   do j = 1, N_monomer
      E(j, j) = 1.0_long
   enddo

   call input(nsch,'n_search')

   read(i_unit,*) comment_line
   call output(comment_line,f='N',j='L')
   do j = 1, N_chain
      call input(dblock(:,j),N_block(j),f='N')
   enddo

   call input(eps, 'eps_tolerance')

   block0 = block_length
   qsq    = qstar(1)
   lmd    = 0.0_long
   dqsq   = 1.0D-4      ! adjustable small parameters
   dlmd   = 1.0D-4      ! adjustable small parameters
   relax  = 1.0D-1

   maxerr = 1.0D8
   itr = 0
   call output('---------------',f='N',j='L')
   spinodal_loop : do
      itr = itr + 1

      old = res(qsq, lmd)
      if (old < eps) then
         write(6,*) "spinodal point found successful"
         exit spinodal_loop
      else if (old > maxerr) then
         write(6,*) "error too large, exit loop"
         exit spinodal_loop
      else if (itr > nsch) then
         write(6,*) "maximum searching # reached, exit loop"
         exit spinodal_loop
      endif

      qsqold  = qsq
      lmdold  = lmd
      gradold = num_grad(old,qsqold,lmdold)

      qsqnew  = qsqold + dqsq
      new     = res(qsqnew, lmdold)
      gradnew = num_grad(new,qsqnew,lmdold)
      Hessen(:,1) = (gradnew - gradold) / dqsq

      lmdnew  = lmdold + dlmd
      new     = res(qsqold, lmdnew)
      gradnew = num_grad(new,qsqold,lmdnew)
      Hessen(:,2) = (gradnew - gradold) / dlmd

      Hesseninv = inv(Hessen, 2)
      qsq = qsqold - relax*dot_product( Hesseninv(1,:), gradold )
      lmd = lmdold - relax*dot_product( Hesseninv(2,:), gradold )

   end do spinodal_loop
   old = res(qsq, lmd, eigv)
   call output('---------------',f='N',j='L')

   call output(itr, 'itr      =',o=6,f='L')
   call output(old, 'residual =',o=6,f='L')
   call output(qsq, 'qstar^2  =',o=6,f='L')
   write(6,*)
   block0 = block0 + lmd * dblock
   call output('block_length',f='N',j='L')
   do k = 1, N_chain
      call output(block0(:,k),N_block(k),f='N')
   enddo
   write(6,*)

   contains

      !----------------------------------------------------------
      ! numerical gradient calculation
      function num_grad(r0,qsq0,lmd0)
      implicit none
      real(long)  :: num_grad(2)
      real(long)  :: r0, qsq0, lmd0
      real(long)  :: r1, qsq1, lmd1

      qsq1 = qsq0 + dqsq
      r1   = res(qsq1, lmd0)
      num_grad(1) = (r1 - r0) / dqsq

      lmd1 = lmd0 + dlmd
      r1   = res(qsq0, lmd1)
      num_grad(2) = (r1 - r0) / dlmd
      end function num_grad
      !----------------------------------------------------------

      !----------------------------------------------------------
      ! residual calculation
      function res(qq,lmd,eigv)
      implicit none
      real(long)          :: res, qq, lmd
      real(long),optional :: eigv(:,:)
      integer             :: i

      qstar(1) = qq
      block_length = block0 + lmd * dblock
      call make_correlation(1,qstar)

      R0    = corrlt(:,:,1) 
      sigma = sum( R0 )
      R1    = E - matmul( ones, R0 ) / sigma
      R1    = matmul( R0, R1 )
      foo   = E + matmul( R1, chi )
      R     = inv( foo, N_monomer )
      R     = matmul( R, R1 )

      do i = 1, cdim
         Rcore(i,:) = R(i,1:cdim) - R(N_monomer,1:cdim)
      enddo
      Ri = inv( Rcore, cdim )
   
      if ( present(eigv) ) then
         call dgeev('N','V',cdim,Ri,cdim,wr,wi,vl,cdim,vr,cdim,work,lwork,info) 
         do i = 1, cdim
            write(6,*) "eigen ",i," = ", wr(i)
            write(6,*) vr(:,i)
            eigv(1:cdim,i) = vr(:,i)
            eigv(N_monomer,i) = wr(i) * dot_product(R(N_monomer,1:cdim),vr(:,i))
            eigv(1:cdim,i) = eigv(1:cdim,i) + eigv(N_monomer,i)
            write(6,*) eigv(:,i)
         enddo
      else
         call dgeev('N','N',cdim,Ri,cdim,wr,wi,vl,cdim,vr,cdim,work,lwork,info) 
      endif
      if(info/=0) stop "Eigenvalue search failed."
      res = minval(abs(wr))
      res = res * res
   
      end function res
      !----------------------------------------------------------

      !----------------------------------------------------------
      function inv(x,xdim)
      implicit none
      real(long) :: x(:,:)
      integer    :: xdim
      real(long) :: inv(xdim, xdim)
      integer    :: lda
      lda = xdim
      ! LU factorization
      call dgetrf(xdim,xdim,x,lda,ipvt(1:xdim),info)
      if(info/=0) stop "LU factorization failed."
      ! matrix inversion
      call dgetri(xdim,x,lda,ipvt(1:xdim),work,lwork,info)
      if(info/=0) stop "Matrix inversion failed."
      inv = x
      end function inv
      !----------------------------------------------------------
   end subroutine rpa_homo_old
   !===================================================================

   !---------------------------------------------------------------------
   !****p spinodal_mod/rpa_blend
   ! SUBROUTINE
   !    rpa_blend
   ! PURPOSE
   !    For blend of A/B/A-B where A-B might be polydisperse
   !    with fixed monomer volume fractions, calculate the
   !    spinodal wave vector and chain length
   ! COMMENTS
   !    First minimize w.r.t. to q --> qsq_star
   !    Then  minimize w.r.t. to N --> N_star
   ! SOURCE
   !---------------------------------------------------------------------
   subroutine rpa_blend(i_unit)
   implicit none

   integer,intent(IN)         :: i_unit
   !***

   ! local variables
   real(long)                 :: homo_chain_length(N_monomer)
   real(long)                 :: homo_chain_length0(N_monomer)
   real(long)                 :: homo_vol_fraction(N_monomer)
   real(long)                 :: homo_phi, d_homo_phi
   real(long)                 :: q_star
   real(long)                 :: N_star
   real(long)                 :: R0(N_monomer,N_monomer)
   real(long)                 :: R_homo(N_monomer,N_monomer)
   real(long)                 :: gbar
   real(long),parameter       :: N_inc=1.0D-3
   real(long),parameter       :: qsq_tol=1.0D-5
   real(long)                 :: lam_tol
   real(long)                 :: oldres,newres
   integer                    :: Nitr, itrmax
   integer                    :: i, imax
   character(len = 100)       :: comment_line

   read(i_unit,*) comment_line
   call output(comment_line,f='N',j='L')
   call input(homo_chain_length,N_monomer,f='N')
   homo_chain_length0 = homo_chain_length

   read(i_unit,*) comment_line
   call output(comment_line,f='N',j='L')
   call input(homo_vol_fraction,N_monomer,f='N')

   call input(homo_phi,'homo_phi')
   call input(d_homo_phi,'d_homo_phi')
   call input(imax,'n_search')

   call input(lam_tol, 'eps_tolerance')

   q_star = 1.0_long
   N_star = 1.0_long
   itrmax = 2000

   open(15, file="blendRPA.dat", status='unknown')
   do i = 0, imax
      if(homo_phi > 1.0_long) homo_phi = 1.0_long
      if(homo_phi < 0.0_long) homo_phi = 0.0_long

      oldres = min_eigen(q_star)
      write(15,*) homo_phi, q_star, oldres/2.0_long

      homo_phi = homo_phi + d_homo_phi
   end do
   close(15)

   contains

      !----------------------------------------------------------
      function f2(x)
      implicit none
      real(long)            :: f2, x
      real(long),parameter  :: eps2 = 1.0D-3
      f2 = ( f1(x+eps2) -  f1(x) ) / eps2
      end function f2
      !----------------------------------------------------------

      !----------------------------------------------------------
      function f1(x)
      implicit none
      real(long)            :: f1, x
      real(long),parameter  :: eps1 = 1.0D-3
      f1 = ( eigen(x+eps1) - eigen(x) )/ eps1
      end function f1
      !----------------------------------------------------------

      !----------------------------------------------------------
      function min_eigen(x)
      implicit none
      real(long) :: min_eigen, x
      integer    :: qitr
      qitr = 0
      q_loop : do
        if ( qitr > itrmax ) then
           write(6,*) "exceeding max search in min_eigen/q_loop"
           exit q_loop
        else if ( abs( f1(x) ) < qsq_tol ) then
           exit q_loop
        end if
        if (abs(x)<1.0D-3) then
           x = 0.0d-3
           exit q_loop
        else
           x = x - f1(x)/f2(x)
           x = abs(x)
        endif
        qitr = qitr + 1
      end do q_loop
      min_eigen = eigen(x)
      end function min_eigen
      !----------------------------------------------------------

      !----------------------------------------------------------
      ! PURPOSE
      !   calculate the eigenvalues of 2x2 structure matrix
      !   for triblocks
      !     R0   = 3x3, ideal Gaussian chain correlation function
      !     Ri   = 3x3, inverse of R0
      !     gbar = reduced Ri, by using incompressibility
      !----------------------------------------------------------
      function eigen(x)
      implicit none
      real(long) :: eigen, x
      real(long) :: b1, b2
      integer    :: j

      qstar(1) = x**2
      call make_correlation(1,qstar)
      R0 = (1.0d0 - homo_phi) * corrlt(:,:,1) 

      R_homo = 0.0d0
      do j = 1, N_monomer
         b1 = Kuhn(j)**2 * qstar(1) * homo_chain_length(j) / 6.0_long
         if ( b1 < 1.0d-8 ) then
             b2 = 1.0d0
         else
             b2 = 2.0d0 * ( b1 - 1.0d0 + exp(-b1) ) / b1 / b1
         endif
         R_homo(j,j) = homo_vol_fraction(j) * homo_chain_length(j) * b2
      end do
      R0 = R0 + homo_phi * R_homo

      b1 = R0(1,1) + 2.0_long * R0(1,2) + R0(2,2)
      b2 = R0(1,1) * R0(2,2) - R0(1,2) * R0(1,2)
      eigen = b1 / b2

      end function eigen
      !----------------------------------------------------------
   end subroutine rpa_blend
   !===================================================================



   !---------------------------------------------------------------------
   !****p spinodal_mod/triblock_rpa_homo
   ! SUBROUTINE
   !    triblock_rpa_homo
   ! PURPOSE
   !    Calculate critical the points where the eigenvalues of both AC 
   !    modulation and B modulation are 0
   ! SOURCE
   !---------------------------------------------------------------------
   subroutine triblock_rpa_homo(i_unit)
   implicit none

   integer,intent(IN)         :: i_unit
   !***

   ! local variables
   real(long)                 :: dblock(N_blk_max,N_chain)
   real(long)                 :: chi0(N_monomer,N_monomer)
   real(long),dimension(3,3)  :: R0, Ri
   real(long),dimension(2,2)  :: gbar
   real(long),parameter       :: qsq_tol=1.0D-5
   real(long)                 :: lam_pos, lam_neg, lam_tol
   real(long)                 :: qsq_pos, qsq_neg
   real(long)                 :: N_pos, N_neg
   real(long),parameter       :: N_inc=1.0D-3
   real(long)                 :: oldres,newres
   integer                    :: Nitr
   integer                    :: itr, itrmax
   integer                    :: i,j,imin,imax
   character(len = 100)       :: comment_line

   call input(imax,'n_search')
   read(i_unit,*) comment_line
   call output(comment_line,f='N',j='L')
   do j = 1, N_chain
      call input(dblock(:,j),N_block(j),f='N')
   enddo
   call input(lam_tol, 'eps_tolerance')

   qsq_pos = 1.0_long
   qsq_neg = 1.0_long
   N_pos   = 1.0_long
   N_neg   = 1.0_long
   itrmax  = 2000
   chi0    = chi

   open(15, file="trispindl.dat", status='unknown')
   open(16, file="qstar.dat", status='unknown')
   do i = 0, imax
      block_length = block_length + dblock

      Nitr = 0
      chi  = N_pos * chi0
      Npos_loop : do

        ! for asymmetric ABA triblock, use the negtive branch
        if ( chi(1,3) < 1.0D-12 .and. &
             abs(block_length(1,1)-block_length(3,1))>1.0D-12 ) &
             exit Npos_loop

        oldres = min_eigen(qsq_pos,1)

        if ( Nitr > itrmax ) then
           write(6,*) "exceeding max search in N_loop"
           exit Npos_loop
        else if ( abs( oldres ) < lam_tol ) then
           exit Npos_loop
        else if ( chi(1,3) >= 4.0_long * chi(1,2) ) then
           exit Npos_loop
        end if

        chi = (N_pos + N_inc) * chi0
        newres = min_eigen(qsq_pos,1)
        N_pos = N_pos - oldres * N_inc / (newres - oldres)

        chi = N_pos * chi0
        Nitr = Nitr + 1
      end do Npos_loop
      lam_pos = oldres

      Nitr = 0
      chi = N_neg * chi0
      Nneg_loop : do

        ! for symmetric ABA triblock, use the postive branch
        if ( chi(1,3) < 1.0D-12 .and. &
             abs(block_length(1,1)-block_length(3,1))<1.0D-12 ) &
             exit Nneg_loop

        oldres = min_eigen(qsq_neg,2)

        if ( Nitr > itrmax ) then
           write(6,*) "exceeding max search in N_loop"
           exit Nneg_loop
        else if ( abs( oldres ) < lam_tol ) then
           exit Nneg_loop
        end if

        chi = (N_neg + N_inc) * chi0
        newres = min_eigen(qsq_neg,2)
        N_neg = N_neg - oldres * N_inc / (newres - oldres)

        chi = N_neg * chi0
        Nitr = Nitr + 1
      end do Nneg_loop
      lam_neg = oldres

      if ( chi(1,3) > 1.0D-12 ) then
         write(15,*) block_length(2,1), N_pos*chi0(1,3), N_neg*chi0(1,3)
         write(16,*) block_length(2,1), qsq_pos, qsq_neg
      else
         write(15,*) block_length(2,1), N_pos*chi0(1,2), N_neg*chi0(1,2)
         write(16,*) block_length(2,1), qsq_pos, qsq_neg
      end if
   end do
   close(15)
   close(16)

!  block_length(:,1) = (/0.1840, 0.6320, 0.1840/)
!  block_length(:,1) = (/0.1886, 0.6228, 0.1886/)
!  call qsweep(42.48_long,1)
!  call chisweep(42.25_long,1)

!  block_length(:,1) = (/0.1242, 0.7516, 0.1242/)
!  call qsweep(37.41_long,2)
!  call chisweep(21.40_long,2)

   contains

      !----------------------------------------------------------
      subroutine chisweep(x,flag)
      implicit none
      real(long)           :: x
      integer              :: flag
      integer              :: i_chi
      real(long)           :: root
      if ( flag == 1 ) then
         open(13, file='1chiswp.dat', status='unknown')
      else
         open(13, file='2chiswp.dat', status='unknown')
      endif
      do i_chi = 10, 300
         chi = dble(i_chi) * 0.01_long * chi0
         root = eigen(x,flag)
         write(13,*) dble(i_chi)*0.01_long*chi0(1,3), root
      end do
      close(13)
      end subroutine chisweep
      !----------------------------------------------------------

      !----------------------------------------------------------
      subroutine qsweep(x,flag)
      implicit none
      real(long)           :: x
      integer              :: flag
      integer              :: i_qsq
      real(long)           :: root
      chi = (x/chi0(1,3))*chi0
      if ( flag == 1 ) then
         open(13, file='1qswp.dat', status='unknown')
      else
         open(13, file='2qswp.dat', status='unknown')
      endif
      do i_qsq = 10, 100
         root = eigen(dble(i_qsq),flag)
         write(13,*) i_qsq, root
      end do
      close(13)
      end subroutine qsweep
      !----------------------------------------------------------

      !----------------------------------------------------------
      function f2(x,flag)
      implicit none
      real(long)            :: f2, x
      integer               :: flag
      real(long),parameter  :: eps2 = 1.0D-3
      f2 = ( f1(x+eps2,flag) -  f1(x,flag) ) / eps2
      end function f2
      !----------------------------------------------------------

      !----------------------------------------------------------
      function f1(x,flag)
      implicit none
      real(long)            :: f1, x
      integer               :: flag
      real(long),parameter  :: eps1 = 1.0D-3
      f1 = ( eigen(x+eps1,flag) - eigen(x,flag) )/ eps1
      end function f1
      !----------------------------------------------------------

      !----------------------------------------------------------
      function min_eigen(x,flag)
      implicit none
      real(long) :: min_eigen,x
      integer    :: flag, qitr
      qitr = 0
      q_loop : do
        if ( qitr > itrmax ) then
           write(6,*) "exceeding max search in min_eigen/q_loop"
           exit q_loop
        else if ( abs( f1(x,flag) ) < qsq_tol ) then
           exit q_loop
        end if
        x = x - f1(x,flag)/f2(x,flag)
        qitr = qitr + 1
      end do q_loop
      min_eigen = eigen(x,flag)
      end function min_eigen
      !----------------------------------------------------------

      !----------------------------------------------------------
      ! PURPOSE
      !   calculate the eigenvalues of 2x2 structure matrix
      !   for triblocks
      !     R0   = 3x3, ideal Gaussian chain correlation function
      !     Ri   = 3x3, inverse of R0
      !     gbar = reduced Ri, by using incompressibility
      !----------------------------------------------------------
      function eigen(x,flag)
      implicit none
      real(long) :: eigen, x
      integer    :: flag
      real(long) :: b1, b2, c

      qstar(1) = x
      call make_correlation(1,qstar)
      R0 = corrlt(:,:,1) 
      Ri = inv( R0, 3 )
      gbar(1,1) = Ri(1,1) - 2*Ri(1,2) + Ri(2,2) - 2.0_long*chi(1,2)
      gbar(2,2) = Ri(3,3) - 2*Ri(3,2) + Ri(2,2) - 2.0_long*chi(2,3)
      gbar(1,2) = Ri(1,3) + Ri(2,2) - Ri(1,2) - Ri(2,3) &
                  + chi(1,3) - chi(1,2) - chi(2,3)
      gbar(2,1) = gbar(1,2)

      b1 = gbar(1,1) +  gbar(2,2)
      b2 = dsqrt( (gbar(1,1)-gbar(2,2))**2 + 4.0_long*gbar(1,2)**2 )
      ! For symmetric case, the two branches of eigenvalues will cross
      ! each other -- degenerate. To avoid the ambiguity of sign for
      ! b2, we fixed it when gbar(1,1) and gbar(2,2) become too close.
      if( abs(gbar(1,1)-gbar(2,2)) < 1.0D-12 ) b2 = 2.0_long*gbar(1,2)
      !if( abs(block_length(1,1)-block_length(3,1)) < 1.0D-12) &
      !        b2 = 2.0_long*gbar(1,2)

      if (flag == 1) then
         eigen = (b1 + b2) / 2.0_long
      else if (flag == 2) then
         eigen = (b1 - b2) / 2.0_long
      else
         stop "Wrong flag in triblock_rpa_homo/eigen ..."
      endif
      end function eigen
      !----------------------------------------------------------
    
      !----------------------------------------------------------
      function inv(x,xdim)
      implicit none
      real(long) :: x(:,:)
      integer    :: xdim
      real(long) :: inv(xdim, xdim)
      integer    :: lda, info
      lda = xdim
      ! LU factorization
      call dgetrf(xdim,xdim,x,lda,ipvt(1:xdim),info)
      if(info/=0) stop "LU factorization failed."
      ! matrix inversion
      call dgetri(xdim,x,lda,ipvt(1:xdim),work,lwork,info)
      if(info/=0) stop "Matrix inversion failed."
      inv = x
      end function inv
      !----------------------------------------------------------
   end subroutine triblock_rpa_homo
   !===================================================================



   !---------------------------------------------------------------------
   !****p spinodal_mod/triblock_bimode
   ! SUBROUTINE
   !    triblock_bimode
   ! PURPOSE
   !    Calculate critical the points where the eigenvalues of both AC 
   !    modulation and B modulation are 0
   ! SOURCE
   !---------------------------------------------------------------------
   subroutine triblock_bimode
   implicit none
   !***

   ! local variables
   real(long)                 :: dblock(N_blk_max,N_chain)
   real(long),dimension(3,3)  :: R0, Ri
   real(long),dimension(2,2)  :: gbar
   real(long)                 :: Bqsq, ACqsq, a, b
   real(long),parameter       :: tol = 1.0D-5
   integer                    :: itr, itrmax
   integer                    :: i,imin,imax

   qstar(1)    = 1.0_long
   dblock(:,1) = (/-0.005, 0.01, -0.005/)
!  dblock(:,1) = (/-0.05, 0.1, -0.05/)
   imin = int( 0.01/dblock(2,1) )
   imax = int( 0.98/dblock(2,1) )

   Bqsq  = 10.0_long
   ACqsq =  1.0_long

   write(6,*) "=========================================="
   itrmax = 2000
   open(15, file="bimod.dat", status='unknown')
   open(16, file="qstar.dat", status='unknown')
   do i = imin, imax
      block_length(:,1) = (/0.5, 0.0, 0.5/) + i * dblock(:,1)  ! (/0.267,0.466,0.267/)

      ! search AC modulation, min a(q) in Erukhimovich's notation
      itr = 0
      AC_loop : do
        if ( itr > itrmax ) then
           write(6,*) "exceeding max search in AC_loop"
           exit AC_loop
        else if ( abs( f1(ACqsq,1) ) < tol ) then
           exit AC_loop
        end if
        ACqsq  = ACqsq - f1(ACqsq,1)/f2(ACqsq,1)
        itr = itr + 1
      end do AC_loop

      ! search B modulation, min b(q) in Erukhimovich's notation
      itr = 0
      B_loop : do
        if ( itr > itrmax ) then
           write(6,*) "exceeding max search in B_loop"
           exit B_loop
        else if ( abs( f1(Bqsq,2) ) < tol ) then
           exit B_loop
        end if
        Bqsq  = Bqsq - f1(Bqsq,2)/f2(Bqsq,2)
        itr = itr + 1
      end do B_loop

      a = fAC(ACqsq)
      b = fB(Bqsq)
      !write(15,*) (a+b)/4.0_long, a, block_length(1,1)
      write(15,*) block_length(1,1),a/((a+b)/4.0_long)
      write(16,*) block_length(1,1), sqrt(ACqsq), sqrt(Bqsq)
   end do
   close(15)
   close(16)
   write(6,*) "=========================================="

   contains
      !----------------------------------------------------------
      function f2(x,flag)
      implicit none
      real(long)            :: f2, x
      integer               :: flag
      real(long),parameter  :: eps2 = 1.0D-3
      f2 = ( f1(x+eps2,flag) -  f1(x,flag) ) / eps2
      end function f2
      !----------------------------------------------------------

      !----------------------------------------------------------
      function f1(x,flag)
      implicit none
      real(long)            :: f1, x
      integer               :: flag
      real(long),parameter  :: eps1 = 1.0D-3
      if ( flag == 1 ) then
         f1 = ( fAC(x+eps1) - fAC(x) )/ eps1
      else if ( flag == 2 ) then
         f1 = ( fB(x+eps1) - fB(x)   )/ eps1
      else
         stop "wrong flag when calculating derivative"
      endif
      end function f1
      !----------------------------------------------------------

      !----------------------------------------------------------
      function fAC(x)
      implicit none
      real(long) :: fAC, x
      qstar(1) = x
      call make_correlation(1,qstar)
      R0 = corrlt(:,:,1) 
      Ri = inv( R0, 3 )
      gbar(1,1) = Ri(1,1) - 2*Ri(1,2) + Ri(2,2)
      gbar(2,2) = Ri(3,3) - 2*Ri(3,2) + Ri(2,2)
      gbar(1,2) = Ri(1,3) + Ri(2,2) - Ri(1,2) - Ri(2,3)
      gbar(2,1) = gbar(1,2)
      fAC = gbar(1,1) - gbar(1,2)
      end function fAC
      !----------------------------------------------------------

      !----------------------------------------------------------
      function fB(x)
      implicit none
      real(long) :: fB, x
      qstar(1) = x
      call make_correlation(1,qstar)
      R0 = corrlt(:,:,1) 
      Ri = inv( R0, 3 )
      gbar(1,1) = Ri(1,1) - 2*Ri(1,2) + Ri(2,2)
      gbar(2,2) = Ri(3,3) - 2*Ri(3,2) + Ri(2,2)
      gbar(1,2) = Ri(1,3) + Ri(2,2) - Ri(1,2) - Ri(2,3)
      gbar(2,1) = gbar(1,2)
      fB = gbar(1,1) + gbar(1,2)
      end function fB
      !----------------------------------------------------------
    
      !----------------------------------------------------------
      function inv(x,xdim)
      implicit none
      real(long) :: x(:,:)
      integer    :: xdim
      real(long) :: inv(xdim, xdim)
      integer    :: lda, info
      lda = xdim
      ! LU factorization
      call dgetrf(xdim,xdim,x,lda,ipvt(1:xdim),info)
      if(info/=0) stop "LU factorization failed."
      ! matrix inversion
      call dgetri(xdim,x,lda,ipvt(1:xdim),work,lwork,info)
      if(info/=0) stop "Matrix inversion failed."
      inv = x
      end function inv
      !----------------------------------------------------------
   end subroutine triblock_bimode
   !===================================================================

end module spinodal_mod
