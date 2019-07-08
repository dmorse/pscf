!----------------------------------------------------------------------
! PROGRAM
!    switch
! PURPOSE
! SOURCE
!----------------------------------------------------------------------
program switch
   use const_mod
   use string_mod
   implicit none

   integer      :: i, j
   integer      :: N, N_monomer
   integer      :: table(:)      !
   real*8       :: omega(:)      ! exchange chemical potential field 
                                 ! omega(monomer,basis function)
   real*8       :: swap(:)       ! temporary omega
   integer      :: indice(:)
   character*30 :: filename
   character*80 :: record

   ! Read driver input
   read(5,*) filename
   read(5,*) N_monomer
   allocate(table(N_monomer))
   do i=1, N_monomer
      read(5,*) table(i)
   enddo

   ! Read Omega file
   open(unit=10, file = trim(filename) )
   read(10,*) record
   write(6,*) trim(record)
   read(10,*) dim
   write(6,*) dim
   do i=1,11
      read(10,*) record
      write(6,*) trim(record)
   enddo
   read(10,*) N
   write(6,'(I15)') N
   allocate(omega(N_monomer))
   allocate(swap(N_monomer))
   allocate(indice(1+dim))
   do i = 1, N
      read(10,*) swap(:), indice(1:1+dim)
      do j=1, N_monomer
         omega(j) = swap(table(j))
      enddo
      write(6,'('//trim(int_string(N_monomer))//'ES20.12,4X,4I4)') omega(1:N_monomer), indice(1:1+dim)
   enddo
   close(10)

end program switch
