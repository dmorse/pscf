
  program main
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Program reads in generators of a group and outputs the
! full group, in the same format.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  USE const_mod
  USE group_mod
  implicit none
 
  integer, parameter :: pdim = 2
  integer            :: i
  type(group_type)   :: g
  real(long)         :: R_basis(pdim,pdim), G_basis(pdim,pdim)
 
  dim = pdim

  !Set up basis (orthorombic example)
  R_basis = 0.0_long
  do i=1, dim
     R_basis(i,i) = 1.00_long
  enddo
  !call output(R_basis,6)
  !G_basis = 4.0*acos(0.0)*Transpose(Inverse(R_basis))
  G_basis  = 4.0*acos(0.0)
  !call output(G_basis,6)

  !Read generators from standard input
  call read_group(g,5)

  !Complete group and output result
  call make_group(g,R_basis,G_basis)
  call output(g,6)

end program main
