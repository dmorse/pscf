   program new_group
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   use const_mod
   use io_mod
   use group_mod
   use unit_cell_mod
   implicit none

   type(group_type) :: group
   character(60)    :: group_name     
   character*5      :: input_format, output_format
   integer          :: i,j,k
   logical          :: complete_group

   call set_echo(0)                 ! Echo off
   call set_com_style('A','A','A') 
   call set_com_use('R')         
   call set_io_units(i=5,o=6)  

   ! Read in group or group generators 
   call input(dim,'dim')
   call input(input_format,'input_format')
   select case(trim(input_format))
   case('S')
      call input(group_name,'group_name')    
      open(65,file = trim(group_name), status = 'old')
      call read_group(group,65)
      close(65)
   case default
      ! allow for other formats, not implemented
   end select 

   ! Complete group if required
   call input(complete_group,'complete_group') 
   if (complete_group) then
      call input_cell_param(5,'F')        
      call make_unit_cell 
      call make_G_basis(R_basis,G_basis)
      call make_group(group,R_basis,G_basis)
   endif
 
   ! Output group
   call input(output_format,'output_format')
   select case(trim(output_format))
   case('S') ! Output using standard format

      call output_group(group,6)

   case('F') ! Output as fortran90

      write(6,*)
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

   end program new_group
