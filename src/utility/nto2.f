!----------------------------------------------------------------------
! PROGRAM
!    switch n monomer to 2 monomer
! PURPOSE
! SOURCE
!----------------------------------------------------------------------
program switch
   use const_mod
   implicit none

   integer      :: i, j
   integer      :: N, N_monomer
   integer      :: table(2)      !
   real*8       :: omega(:,:)    ! chemical potential field 
                                 ! omega(monomer,basis function)
   integer      :: indice(:)
   character*30 :: filename
   character*80 :: record

   ! Read driver input
   read(5,*) filename
   read(5,*) N_monomer
   do i=1, 2
      read(5,*) table(i)
   enddo

   ! Read Omega file
   open(unit=10 , file = trim(filename) )
   read(10,*) record
   write(6,*) trim(record)
   read(10,*) dim
   write(6,*) dim
   do i=1,9
      read(10,*) record
      write(6,*) trim(record)
   enddo
   read(10,*)
   write(6,*) 2
   read(10,*) record
   write(6,*) trim(record)
   read(10,*) N
   write(6,'(I15)') N
   allocate(omega(N_monomer,N))
   allocate(indice(dim+1))
   do i = 1, N
      read(10,*) omega(:,i), indice(1:dim+1)
      write(6,'(2ES20.12,4X,4I4)') omega(table(1),i),omega(table(2),i),indice(1:dim+1)
   enddo
   close(10)

end program switch
