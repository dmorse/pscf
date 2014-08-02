   program new_2dgroups
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! This program was used to generate the select case statement for
   ! for 2D plane groups in the space_groups routine of space_groups_mod.
   ! It is similar to the program new_group.f, except for addition of 
   ! a outer loop over the groups named in a driver file.
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   use const_mod
   use io_mod
   use group_mod
   use unit_cell_mod
   implicit none

   type(group_type) :: group
   character(60)    :: group_name     
   integer          :: i,j,k, i_group

   call set_echo(0)                 ! Echo off
   call set_com_style('A','A','A') 
   call set_com_use('R')         
   call set_io_units(i=5,o=6)  

   ! Complete group if required
   do i_group=1, 17

      read(5,*)
      call input(group_name,'group_name')    

      call input_unit_cell(5,'F')        
      call make_unit_cell 
      call make_G_basis(R_basis,G_basis)

      open(65,file = trim(group_name), status = 'old')
      call read_group(group,65)
      close(65)

      call make_group(group,R_basis,G_basis)
    
      ! Output group
      select case('F')
      case('S') ! Output using standard format
   
         call output_group(group,6)
   
      case('F') ! Output as fortran90
   
         write(6,*)
         write(6,"(6X, 'case(',A12,')')") group_name
         write(6,"(9X, 'g%order   =',i4)") group%order
         write(6,"(9X, 'do i=1, g%order')")
         write(6,"(12X,'g%s(i)%m = 0.0_long')")
         write(6,"(12X,'g%s(i)%v = 0.0_long')")
         write(6,"(9X, 'enddo')")
         do k=1, group%order
            write(6,*)
            do i=1, dim
               do j=1, dim
                  if (abs(group%s(k)%m(i,j)) > 1.0E-4) then
                     if (k < 10) then
                        write(6, &
                        "(9X,'g%s(',i1,')%m(',i1,',',i1,')   =',F6.2)")&
                        k,i,j,group%s(k)%m(i,j)
                     else
                        write(6,&
                        "(9X,'g%s(',i2 ,')%m(',i1,',',i1,')  =',F6.2)")&
                        k,i,j,group%s(k)%m(i,j)
                     endif
                  endif
               enddo
            enddo
            do j=1, dim
               if ( abs(group%s(k)%v(j)) > 1.0E-4 ) then
                  if (k < 10) then
                     write(6, &
                     "(9X,'g%s(',i1,')%v(',i1,')     =',F6.2)")&
                     k,j,group%s(k)%v(j)
                  else
                     write(6,&
                     "(9X,'g%s(',i2 ,')%v(',i1,')    =',F6.2)")&
                     k,j,group%s(k)%v(j)
                  endif
               endif
            enddo
         enddo
   
      end select

   end do

   end program new_2dgroups
