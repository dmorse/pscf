   program main
!----------------------------------------------------
! Program to test group_mod
!----------------------------------------------------
   USE const_mod
   USE group_mod
   implicit none

   integer :: readdim 
   integer :: i,j,k,l
   real(long)                 :: r1, r2, r3, r4, r5
   integer,    dimension(3)   :: i1, i2, i3, i4, i5, i6 
   real(long), dimension(3)   :: v1, v2, v3, v4, v5, v6
   real(long), dimension(3,3) :: m1, m2, m3, m4, m5, m6, m7
   type(symmetry_type)        :: s1, s2, s3, s4, s5, s6, s7
   type(group_type)           :: g
   type(table_type)           :: table
 
   real(long), allocatable, dimension(:,:) :: R_basis, G_basis

   ! Choose dimension of group to be read from standard input
   ! If readdim = 0, this test is bypassed

   readdim = 3

   write(6,*) 
   write(6,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
   write(6,*)'%        Tests for dim =3                     %'
   write(6,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
   write(6,*) 
 
   dim = 3
   allocate(R_basis(3,3))
   allocate(G_basis(3,3))

   R_basis = 0.0_long
   G_basis = 0.0_long
   do i=1, dim
      R_basis(i,i) = 1.0_long
      G_basis(i,i) = 4.0*acos(0.0)
   enddo
 
   !----------------------------------------------
   ! Test vectors and matrix operations
   !----------------------------------------------

   r1      = 2.0
   r2      = 3.0
 
   v1(1) = 1.0
   v1(2) = 2.0
   v1(3) = 3.0
 
   v2(1) = 1.0
   v2(2) = 0.0
   v2(3) = 0.0
 
   m1      = 0.0_long
   m1(1,1) = 1.0
   m1(2,2) = 2.0
   m1(3,3) = 1.0
 
   m2      = 0.0_long
   m2(1,2) = 3.0
   m2(2,1) = 4.0
   m2(3,3) = 2.0
 
   write(6,*) "r1="
   call output(r1,6)
   write(6,*) "r2="
   call output(r2,6)
   write(6,*) "v1="
   call output(v1,6)
   write(6,*) "v2="
   call output(v2,6)
   write(6,*) "m1="
   call output(m1,6)
   write(6,*) "m2="
   call output(m2,6)
  
   r3 = v1.dot.v2
   write(6,*) "r3=v1*v2="
   call output(r3,6)
  
   v4 = m2.dot.v1
   write(6,*) "v4=m2*v1="
   call output(v4,6)
 
   m4 = m2.dot.m2
   write(6,*) "m4=m2*m2="
   call output(m4,6)
 
   m5 = Inverse(m2)
   write(6,*) "m5 = 1/m2="
   call output(m5,6)
  
   m6 = m5.dot.m2
   write(6,*) "m6 = m5*m2="
   call output(m6,6)
 
   m6 = m4.dot.Inverse(m4)
   write(6,*) "m6 = m4*Inverse(m4)"
   call output(m6,6)
 
   !-----------------------------------------------
   ! Test of operations on symmetries
   !-----------------------------------------------
 
   s1%basis = 'Bravais  '
   s1%m = m1
   s1%v = v1
   write(6,*) "s1="
   call output(s1,6)
 
   s2%basis = 'Bravais  '
   s2%m = m2
   s2%v = v2
   write(6,*) "s2="
   call output(s2,6)
 
   s3 = s1.dot.s2
   write(6,*) "s3 = s1*s2="
   call output(s3,6)
 
   s4 = s3.dot.Inverse(s3)
   write(6,*) "s4 =s3*Inverse(s3)"
   call output(s4,6)
 
   s4 = s3.dot.Inverse(s2)
   write(6,*) "s4 =s3*Inverse(s2)"
   call output(s4,6)
 
   s4 = Inverse(s1).dot.s3
   write(6,*) "s4 =Inverse(s1)*s3"
   call output(s4,6)
 
   !----------------------------------------
   !Test of shift_translation
   !----------------------------------------
 
   write(6,*) "Unshifted translation"
   call output(s3%v,6)
   call shift_translation(s3,v3,R_basis,G_basis)
   write(6,*) "Shifted Vector"
   call output(v3,6)
 
   !---------------------------------------
   !Tests of variants of equal function
   !---------------------------------------
  
   if (equal(m2,m1)) then 
       write(6,*)"m2.eq.m1"
   else
       write(6,*)"m2.neq.m1"
   endif
  
   if (equal(s4,s3,R_basis,G_basis)) then 
       write(6,*)"s4.eq.s3"
   else
       write(6,*)"s4.neq.s3"
   endif
  
   !write(6,*)
  
   If (equal(s3,s3,R_basis,G_basis)) then 
        write(6,*)"s3.eq.s3"
   else
       write(6,*)"s3.neq.s3"
   endif
   write(6,*)
  
   s4 = s3
   s4%m(1,2) = s4%m(1,2) + 0.001
   write(6,*) 'Code to define s4:'
   write(6,*) '   s4 = s3'
   write(6,*) '   s4%m(1,2) = s4%m(1,2) + 0.001'
   write(6,*) 'Result of test of equal function:'
   If (equal(s4,s3,R_basis,G_basis)) then 
       write(6,*)"s4.eq.s3"
   else
       write(6,*)"s4.neq.s3"
   endif
   write(6,*)
  
   s4 = s3
   s4%v(1) = s4%v(1) + 0.001
   write(6,*) 'Code to define s4:'
   write(6,*) '   s4 = s3'
   write(6,*) '   s4%v(1) = s4%v(1) + 0.001'
   If (equal(s4,s3,R_basis,G_basis)) then 
       write(6,*)"s4.eq.s3"
   else
       write(6,*)"s4.neq.s3"
   endif
   write(6,*)
 
   write(6,*) 'Test of .not.equal(s3,s3):'
   If (.not.equal(s3,s3,R_basis,G_basis)) then 
       write(6,*) "s3.neq.s3"
   else
       write(6,*)"s3.eq.s3"
   endif
   write(6,*)
  
   write(6,*) 'Test of .not.equal(s4,s3):'
   If (.not.equal(s4,s3,R_basis,G_basis)) then 
       write(6,*)"s4.neq.s3"
   else
       write(6,*)"s4.eq.s3"
   endif
   write(6,*)
  
   !-------------------------------------- 
   !Test of identity constructor
   !--------------------------------------

   g%order = 4
   g%s(1)  = make_identity('Bravais  ')
   write(6,*) "Constructed identity matrix="
   call output(g%s(1),6)
  
  
  !---------------------------------------
  !First test of group operations
  !---------------------------------------
  
   g%s(2)%basis = 'Bravais'
   g%s(2)%m = 0.0_long
   g%s(2)%v = 0.0_long
   g%s(2)%m(1,2) = 1.0_long
   g%s(2)%m(2,1) = 1.0_long
   g%s(2)%m(3,3) = 1.0_long
  
   g%s(3)%basis = 'Bravais'
   g%s(3)%m = 0.0_long
   g%s(3)%v = 0.0_long
   g%s(3)%m(1,1) =-1.0_long
   g%s(3)%m(2,2) = 1.0_long
   g%s(3)%m(3,3) = 1.0_long
  
   g%s(4)%basis = 'Bravais'
   g%s(4)%m = 0.0_long
   g%s(4)%v = 0.0_long
   g%s(4)%m(1,1) = -1.0_long
   g%s(4)%m(2,2) = -1.0_long
   g%s(4)%m(3,3) = -1.0_long
  
   g%s(5)%basis = 'Bravais'
   g%s(5)%m = 0.0_long
   g%s(5)%v = 0.0_long
   g%s(5)%m(1,2) = 1.0_long
   g%s(5)%m(2,3) = 1.0_long
   g%s(5)%m(3,1) = 1.0_long
  
   g%s(6)%basis = 'Bravais'
   g%s(6)%m = 0.0_long
   g%s(6)%v = 0.0_long
   g%s(6)%m(1,1) =  1.0_long
   g%s(6)%m(2,3) =  1.0_long
   g%s(6)%m(3,2) =  1.0_long
  
   write(6,*) "show proto group:"
   write(6,*)
   call output(g,6)
  
   write(6,*) "indices of inverses="
   do i=1, g%order
      j = index_of_inverse(g,i,R_basis,G_basis)
      write(6,FMT='(2I5)') i, j
      ! call output(Inverse(g%s(i)),6)
   enddo
   write(6,*)
  
   write(6,*) "show table"
   call make_table(g,table,R_basis,G_basis)
   call output(table,6)
  
   call make_group(g,R_basis,G_basis)
   write(6,*) "show new group"
   call output(g,6)
  
   write(6,*) "indices of inverses="
   do i=1, g%order
      write(6,FMT='(2I5)') i, index_of_inverse(g,i,R_basis,G_basis)
      ! call output(Inverse(g%s(i)),6)
   enddo
  
   !-----------------------------------------------------------!
   ! Read, complete, and write a group                         !
   !-----------------------------------------------------------!
  
   if (readdim == dim) then 

      call read_group(g,5)
      write(6,*) "show group read from file"
      write(6,*)
      call output(g,6)
   
      call make_group(g,R_basis,G_basis)
      write(6,*) "show new group"
      write(6,*)
      call output(g,6)
  
      write(6,*) "indices of inverses="
      do i=1, g%order
         write(6,FMT='(2I5)') i, index_of_inverse(g,i,R_basis,G_basis)
      enddo

   endif

   deallocate(R_basis) 
   deallocate(G_basis) 

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   write(6,*)
   write(6,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
   write(6,*) '%        Tests for dim = 2                    %'
   write(6,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
   write(6,*)
 
   dim = 2
   allocate(R_basis(2,2))
   allocate(G_basis(2,2))

   R_basis = 0.0_long
   G_basis = 0.0_long
   do i=1, dim
      R_basis(i,i) = 1.0_long
      G_basis(i,i) = 4.0*acos(0.0)
   enddo

 
   !----------------------------------------------
   ! Test vectors and matrix operations
   !----------------------------------------------
 
   r1      = 2.0
   r2      = 3.0
 
   v1(1) = 1.0
   v1(2) = 2.0
 
   v2(1) = 1.0
   v2(2) = 0.0
 
   m1      = 0.0_long
   m1(1,1) = 1.0
   m1(2,2) = 2.0
 
   m2      = 0.0_long
   m2(1,2) = 3.0
   m2(2,1) = 4.0
 
   write(6,*) "r1="
   call output(r1,6)
   write(6,*) "r2="
   call output(r2,6)
   write(6,*) "v1="
   call output(v1,6)
   write(6,*) "v2="
   call output(v2,6)
   write(6,*) "m1="
   call output(m1,6)
   write(6,*) "m2="
   call output(m2,6)
  
   r3 = v1.dot.v2
   write(6,*) "r3=v1*v2="
   call output(r3,6)
  
   v4 = m2.dot.v1
   write(6,*) "v4=m2*v1="
   call output(v4,6)
 
   m4 = m2.dot.m2
   write(6,*) "m4=m2*m2="
   call output(m4,6)
 
   m5 = Inverse(m2)
   write(6,*) "m5 = 1/m2="
   call output(m5,6)
  
   m6 = m5.dot.m2
   write(6,*) "m6 = m5*m2="
   call output(m6,6)
 
   m6 = m4.dot.Inverse(m4)
   write(6,*) "m6 = m4*Inverse(m4)"
   call output(m6,6)
 
   !-----------------------------------------------
   ! Test of operations on symmetries
   !-----------------------------------------------
 
   s1%basis = 'Bravais  '
   s1%m = m1
   s1%v = v1
   write(6,*) "s1="
   call output(s1,6)
 
   s2%basis = 'Bravais  '
   s2%m = m2
   s2%v = v2
   write(6,*) "s2="
   call output(s2,6)
 
   s3 = s1.dot.s2
   write(6,*) "s3 = s1 .dot. s2="
   call output(s3,6)
 
   s4 = s3.dot.Inverse(s3)
   write(6,*) "s4 =s3 .dot. Inverse(s3) ="
   call output(s4,6)
 
   s4 = s3.dot.Inverse(s2)
   write(6,*) "s4 =s3 .dot. Inverse(s2) ="
   call output(s4,6)
 
   s4 = Inverse(s1).dot.s3
   write(6,*) "s4 =Inverse(s1)*s3 ="
   call output(s4,6)
 
   !----------------------------------------
   ! Test of shift_translation
   !----------------------------------------

   write(6,*) "Unshifted translation = s3%v"
   call output(s3%v,6)
 
   call shift_translation(s3,v3,R_basis,G_basis)
   write(6,*) "Shifted Vector"
   call output(v3,6)
 
 
   !---------------------------------------
   ! Tests of variants of equal function
   !---------------------------------------
  
   if (equal(m2,m1)) then 
       write(6,*)"m2.eq.m1"
   else
       write(6,*)"m2.neq.m1"
   endif
  
   if (equal(s4,s3,R_basis,G_basis)) then 
       write(6,*)"s4.eq.s3"
   else
       write(6,*)"s4.neq.s3"
   endif
  
   !write(6,*)
  
   If (equal(s3,s3,R_basis,G_basis)) then 
        write(6,*)"s3.eq.s3"
   else
       write(6,*)"s3.neq.s3"
   endif
   write(6,*)
  
   s4 = s3
   s4%m(1,2) = s4%m(1,2) + 0.001
   write(6,*) 'Code to define s4:'
   write(6,*) '   s4 = s3'
   write(6,*) '   s4%m(1,2) = s4%m(1,2) + 0.001'
   write(6,*) 'Result of test of equal function:'
   If (equal(s4,s3,R_basis,G_basis)) then 
       write(6,*)"s4.eq.s3"
   else
       write(6,*)"s4.neq.s3"
   endif
   write(6,*)
  
   s4 = s3
   s4%v(1) = s4%v(1) + 0.001
   write(6,*) 'Code to define s4:'
   write(6,*) '   s4 = s3'
   write(6,*) '   s4%v(1) = s4%v(1) + 0.001'
   write(6,*) 'Result of test of equal function:'
   If (equal(s4,s3,R_basis,G_basis)) then 
       write(6,*)"s4.eq.s3"
   else
       write(6,*)"s4.neq.s3"
   endif
   write(6,*)
  
   write(6,*) 'Test of .not.equal(s3,s3):'
   If (.not.equal(s3,s3,R_basis,G_basis)) then 
       write(6,*) "s3.neq.s3"
   else
       write(6,*)"s3.eq.s3"
   endif
   write(6,*)
  
   write(6,*) 'Test of .not.equal(s3,s3):'
   If (.not.equal(s4,s3,R_basis,G_basis)) then 
       write(6,*)"s4.neq.s3"
   else
       write(6,*)"s4.eq.s3"
   endif
   write(6,*)
  
  !-------------------------------------- 
  !Tests of identity constructor
  !--------------------------------------
   g%order = 4
   g%s(1)  = make_identity('Bravais  ')
   write(6,*) "Constructed identity matrix="
   call output(g%s(1),6)
  
  
  !---------------------------------------
  !First test of group operations
  !---------------------------------------
  
   g%s(2)%basis = 'Bravais'
   g%s(2)%m = 0.0_long
   g%s(2)%v = 0.0_long
   g%s(2)%m(1,2) = 1.0_long
   g%s(2)%m(2,1) = 1.0_long
  
   g%s(3)%basis = 'Bravais'
   g%s(3)%m = 0.0_long
   g%s(3)%v = 0.0_long
   g%s(3)%m(1,1) =-1.0_long
   g%s(3)%m(2,2) = 1.0_long
  
   g%s(4)%basis = 'Bravais'
   g%s(4)%m = 0.0_long
   g%s(4)%v = 0.0_long
   g%s(4)%m(1,1) = -1.0_long
   g%s(4)%m(2,2) = -1.0_long
  
   write(6,*) "show proto group:"
   write(6,*) 
   call output(g,6)
  
   write(6,*) "indices of inverses="
   do i=1, g%order
      j = index_of_inverse(g,i,R_basis,G_basis)
      write(6,FMT='(2I5)') i, j
      ! call output(Inverse(g%s(i)),6)
   enddo
   write(6,*)
  
   write(6,*) "show table"
   call make_table(g,table,R_basis,G_basis)
   call output(table,6)
  
   call make_group(g,R_basis,G_basis)
   write(6,*) "show new group"
   call output(g,6)
  
   write(6,*) "indices of inverses="
   do i=1, g%order
      write(6,FMT='(2I5)') i, index_of_inverse(g,i,R_basis,G_basis)
      ! call output(Inverse(g%s(i)),6)
   enddo
   
   !-----------------------------------------------------------!
   ! Read, complete, and write groups                          !
   !-----------------------------------------------------------!
  
   if (readdim == dim) then 

      call read_group(g,5)
      write(6,*) "show group read from file"
      write(6,*)
      call output(g,6)
   
      call make_group(g,R_basis,G_basis)
      write(6,*) "show new group"
      write(6,*)
      call output(g,6)
  
      write(6,*) "indices of inverses="
      do i=1, g%order
         write(6,FMT='(2I5)') i, index_of_inverse(g,i,R_basis,G_basis)
      enddo

   endif

   deallocate(R_basis) 
   deallocate(G_basis) 

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   write(6,*)
   write(6,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
   write(6,*) '%        Tests for dim = 1                    %'
   write(6,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
   write(6,*)
 
   dim = 1
   allocate(R_basis(1,1))
   allocate(G_basis(1,1))

   R_basis = 0.0_long
   G_basis = 0.0_long
   do i=1, dim
      R_basis(i,i) = 1.0_long
      G_basis(i,i) = 4.0*acos(0.0)
   enddo

 
   !----------------------------------------------
   ! Test vectors and matrix operations
   !----------------------------------------------
 
   r1      = 2.0
   r2      = 0.5
 
   v1(1) = 1.0
 
   v2(1) = 2.5
 
   m1      = 3.0_long
 
   m2      = 4.0_long
 
   write(6,*) "r1="
   call output(r1,6)
   write(6,*) "r2="
   call output(r2,6)
   write(6,*) "v1="
   call output(v1,6)
   write(6,*) "v2="
   call output(v2,6)
   write(6,*) "m1="
   call output(m1,6)
   write(6,*) "m2="
   call output(m2,6)
  
   r3 = v1.dot.v2
   write(6,*) "r3=v1*v2="
   call output(r3,6)
  
   v4 = m2.dot.v1
   write(6,*) "v4=m2*v1="
   call output(v4,6)
 
   m4 = m2.dot.m2
   write(6,*) "m4=m2*m2="
   call output(m4,6)
 
   m5 = Inverse(m2)
   write(6,*) "m5 = 1/m2="
   call output(m5,6)
  
   m6 = m5.dot.m2
   write(6,*) "m6 = m5*m2="
   call output(m6,6)
 
   m6 = m4.dot.Inverse(m4)
   write(6,*) "m6 = m4*Inverse(m4)"
   call output(m6,6)
 
   !-----------------------------------------------
   ! Test of operations on symmetries
   !-----------------------------------------------
 
   s1%basis = 'Bravais  '
   s1%m = m1
   s1%v = v1
   write(6,*) "s1="
   call output(s1,6)
 
   s2%basis = 'Bravais  '
   s2%m = m2
   s2%v = v2
   write(6,*) "s2="
   call output(s2,6)
 
   s3 = s1.dot.s2
   write(6,*) "s3 = s1 .dot. s2="
   call output(s3,6)
 
   s4 = s3.dot.Inverse(s3)
   write(6,*) "s4 =s3 .dot. Inverse(s3) ="
   call output(s4,6)
 
   s4 = s3.dot.Inverse(s2)
   write(6,*) "s4 =s3 .dot. Inverse(s2) ="
   call output(s4,6)
 
   s4 = Inverse(s1).dot.s3
   write(6,*) "s4 =Inverse(s1)*s3 ="
   call output(s4,6)
 
   !----------------------------------------
   ! Test of shift_translation
   !----------------------------------------

   write(6,*) "Unshifted translation = s3%v"
   call output(s3%v,6)
 
   call shift_translation(s3,v3,R_basis,G_basis)
   write(6,*) "Shifted Vector"
   call output(v3,6)
 
 
   !---------------------------------------
   ! Tests of variants of equal function
   !---------------------------------------
  
   if (equal(m2,m1)) then 
       write(6,*)"m2.eq.m1"
   else
       write(6,*)"m2.neq.m1"
   endif
  
   if (equal(s4,s3,R_basis,G_basis)) then 
       write(6,*)"s4.eq.s3"
   else
       write(6,*)"s4.neq.s3"
   endif
  
   !write(6,*)
  
   If (equal(s3,s3,R_basis,G_basis)) then 
        write(6,*)"s3.eq.s3"
   else
       write(6,*)"s3.neq.s3"
   endif
   write(6,*)
  
   s4 = s3
   s4%m(1,1) = s4%m(1,1) + 0.001
   write(6,*) 'Code to define s4:'
   write(6,*) '   s4 = s3'
   write(6,*) '   s4%m(1,1) = s4%m(1,1) + 0.001'
   write(6,*) 'Result of test of equal function:'
   If (equal(s4,s3,R_basis,G_basis)) then 
       write(6,*)"s4.eq.s3"
   else
       write(6,*)"s4.neq.s3"
   endif
   write(6,*)
  
   s4 = s3
   s4%v(1) = s4%v(1) + 0.001
   write(6,*) 'Code to define s4:'
   write(6,*) '   s4 = s3'
   write(6,*) '   s4%v(1) = s4%v(1) + 0.001'
   write(6,*) 'Result of test of equal function:'
   If (equal(s4,s3,R_basis,G_basis)) then 
       write(6,*)"s4.eq.s3"
   else
       write(6,*)"s4.neq.s3"
   endif
   write(6,*)
  
   write(6,*) 'Test of .not.equal(s3,s3):'
   If (.not.equal(s3,s3,R_basis,G_basis)) then 
       write(6,*) "s3.neq.s3"
   else
       write(6,*)"s3.eq.s3"
   endif
   write(6,*)
  
   write(6,*) 'Test of .not.equal(s3,s3):'
   If (.not.equal(s4,s3,R_basis,G_basis)) then 
       write(6,*)"s4.neq.s3"
   else
       write(6,*)"s4.eq.s3"
   endif
   write(6,*)
  
  !-------------------------------------- 
  !Tests of identity constructor
  !--------------------------------------
   g%order = 2
   g%s(1)  = make_identity('Bravais  ')
   write(6,*) "Constructed identity matrix="
   call output(g%s(1),6)
  
  
  !---------------------------------------
  !First test of group operations
  !---------------------------------------
  
   g%s(2)%basis = 'Bravais'
   g%s(2)%m = 0.0_long
   g%s(2)%v = 0.0_long
   g%s(2)%m(1,1) = -1.0_long
  
   write(6,*) "show proto group:"
   write(6,*)
   call output(g,6)
  
   write(6,*) "indices of inverses="
   do i=1, g%order
      j = index_of_inverse(g,i,R_basis,G_basis)
      write(6,FMT='(2I5)') i, j
      ! call output(Inverse(g%s(i)),6)
   enddo
   write(6,*)
  
   write(6,*) "show table"
   write(6,*)
   call make_table(g,table,R_basis,G_basis)
   call output(table,6)
  
   call make_group(g,R_basis,G_basis)
   write(6,*) "show new group:"
   write(6,*)
   call output(g,6)
  
   write(6,*) "indices of inverses="
   do i=1, g%order
      write(6,FMT='(2I5)') i, index_of_inverse(g,i,R_basis,G_basis)
      ! call output(Inverse(g%s(i)),6)
   enddo
   
   !-----------------------------------------------------------!
   ! Read, complete, and write groups                          !
   !-----------------------------------------------------------!
   dim = 3
   print *, 'readdim = ', readdim
   print *, 'dim = ', dim
   if (readdim == dim) then 

      call read_group(g,5)
      write(6,*) "show group read from file"
      write(6,*)
      call output(g,6)
   
      call make_group(g,R_basis,G_basis)
      write(6,*) "show new group"
      write(6,*)
      call output(g,6)
  
      write(6,*) "indices of inverses="
      do i=1, g%order
         write(6,FMT='(2I5)') i, index_of_inverse(g,i,R_basis,G_basis)
      enddo

   endif

   end program main
