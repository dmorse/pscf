program homo_response_test
   use const_mod
   use io_mod
   use chemistry_mod
   use rpa_mod
   implicit none

   real(long) :: q, S(2,2)
   integer    :: i, j, k, i_chain, i_example, N_example, i_q, N_q

   ! Set defaults for parameter I/O module - see io_mod
   call set_echo(1)                ! echo inputs
   call set_com_style('A','A','A') ! comments on line above data
   call set_com_use('R')           ! replace comment in echoed output
   call set_io_units(i=5,o=6)      ! set standard in and echo units

   ! Loop over example chemistries
   read(5,*) N_example
   write(6,*) 
   do i_example = 1, N_example

      read(5,*)
      write(6,*)
      call input_chemistry(5,'F')  
      read(5,*)
      write(6,*)
      call input(N_q)
      do i_q=1, N_q
         read(5,*) q
         write(6,*) 'q = ', q
      
         do i_chain=1, N_chain
            write(6,*) 'Output S_block(i_chain,i,j,q)'
            do i=1, N_block(i_chain)
               do j=1, N_block(i_chain)
                  write(6,*) i, j, S_block(i_chain,i,j,q)
               enddo
            enddo
         
            write(6,*) 'Output S_chain'
            call S_chain(i_chain,q,S)
            do i=1, N_monomer
               do j=1, N_monomer
                  write(6,*) i, j, S(i,j)
               enddo
            enddo
         enddo
      
         write(6,*) 'Output S_ideal'
         call S_ideal(q,S,2)
         do i=1, N_monomer
            do j=1, N_monomer
               write(6,*) i, j, S(i,j)
            enddo
         enddo
      
         write(6,*) 'Output S_inc'
         call S_inc(q,S,2)
         do i=1, N_monomer
            do j=1, N_monomer
               write(6,*) i, j, S(i,j)
            enddo
         enddo
      
         write(6,*) 'Output S_rpa'
         call S_rpa(q,S,2)
         do i=1, N_monomer
            do j=1, N_monomer
               write(6,*) i, j, S(i,j)
            enddo
         enddo

      enddo

   enddo

end program homo_response_test
