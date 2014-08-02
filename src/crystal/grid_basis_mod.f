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
!****m  scf/grid_basis_mod
! PURPOSE
!   grid -- to -- basis and reverse conversions
! SOURCE
!-----------------------------------------------------------------------
module grid_basis_mod
   use const_mod
   use grid_mod
   use basis_mod
   implicit none

   PRIVATE
   PUBLIC :: basis_to_kgrid
   PUBLIC :: kgrid_to_basis
   PUBLIC :: check_symmetry
   !***

contains

   !****p grid_basis_mod/basis_to_kgrid --------------------------------
   ! SUBROUTINE
   !   basis_to_kgrid(basis, kgrid, [karray_full], [no_parity])
   ! PURPOSE
   !   Converts representation of field as array "basis" of 
   !   coefficients of symmetrized basis functions into an FFT
   !   array kgrid containing corresponding coefficients for
   !   plane waves
   ! ARGUMENTS
   !   basis     - Coefficients in expansion of field in basis functions
   !   kgrid     - Fourier components of periodic field, on an FFT grid
   !   karray_full - (optional) If absent or false, use only k_x > 0, to
   !               represent real field. If present and true, use entire grid.
   !   no_parity - (optional) If absent or false, use real basis functions
   !               If present and true, use star basis functions.
   ! SOURCE
   !---------------------------------------------------------------------
   subroutine basis_to_kgrid(basis, kgrid, karray_full, no_parity)
   implicit none
   real(long),intent(IN)     :: basis(:)        ! basis(N_star)
   complex(long),intent(OUT) :: kgrid(0:,0:,0:) ! FFT k-space grid
   logical,optional          :: karray_full
   logical,optional          :: no_parity
   !***

   ! local variables
   integer      :: kfft(3)      ! (0:N(1)-1,...), natural order of fft
   integer      :: kbz(3)       ! First Brillouin zone k
   integer      :: i_wave,i_star
   integer      :: i1,i2,i3, i1max
   real(long)   :: a, b, sq2
   complex(long):: c, z
   sq2 = dsqrt(2.0_long)

   i1max = ngrid(1)/2
   if(present(karray_full)) then
      if (karray_full) i1max = ngrid(1)-1
   end if

   kgrid=0.0_long
   do i1=0,i1max
      kfft(1)=i1
      do i2=0,ngrid(2)-1
         kfft(2)=i2
         do i3=0,ngrid(3)-1
            kfft(3)=i3

            kbz    = G_to_bz(kfft)
            i_wave = which_wave(kbz(1),kbz(2),kbz(3))

            if( i_wave<1 .or. i_wave>N_wave ) then

                  z=(0.0D0,0.0D0)

            else

               c=coeff(i_wave)
               i_star=star_of_wave(i_wave)
               
               if ( present(no_parity) .and. ( no_parity ) ) then

                  z = basis(i_star) * c

               else

                  select case( star_invert(i_star) )
                  case(0)
                     z = basis(i_star) * c
                  case(1)
                     if(i_star==N_star)stop "star_invert(N_star)=1"
                     if(star_invert(i_star+1)/=-1)then
                        write(6,*) "star_invert(",i_star,  ") = 1"
                        write(6,*) "star_invert(",i_star+1,")/=-1"
                        stop 
                     endif
                     a = basis(i_star)
                     b = basis(i_star+1)
                     z = dcmplx(a,-b) * c / sq2
                  case(-1)
                     if(star_invert(i_star-1)/=1)then
                        write(6,*) "star_invert(",i_star,  ") = -1"
                        write(6,*) "star_invert(",i_star-1,")/=  1"
                        stop 
                     endif
                     a = basis(i_star-1)
                     b = basis(i_star)
                     z = dcmplx(a, b) * c / sq2
                  case default
                     stop "Illegal star_invert value in basis_to_kgrid."
                  end select

               end if

            end if

            kgrid(i1,i2,i3) = z
         end do
      end do
   end do
   end subroutine basis_to_kgrid
   !==============================================================


   !****p grid_basis_mod/kgrid_to_basis -------------------------------
   ! SUBROUTINE
   !   kgrid_to_basis(kgrid, basis, [no_parity])
   ! PURPOSE
   !   Converts Fourier representation of field into an array basis
   !   of coefficients of symmetrized basis functions
   ! ARGUMENTS
   !   kgrid     - Fourier components of periodic field, on an FFT grid
   !   basis     - Coefficients in expansion of field in basis functions
   !   no_parity - (optional) If absent or false, use real basis functions
   !               If present, use star basis functions.
   ! SOURCE
   !---------------------------------------------------------------------
   subroutine kgrid_to_basis(kgrid, basis, no_parity)
   implicit none
   complex(long),intent(IN)  :: kgrid(0:,0:,0:)
   real(long),intent(OUT)    :: basis(:)     ! basis(N_star)
   logical,optional          :: no_parity
   !***

   !local variables
   integer      :: kfft(3)      ! (0:N(1)-1,...), natural order of fft
   integer      :: kbz(3)       ! First Brillouin zone k
   integer      :: i_star,i_wave
   real(long)   :: sq2
   complex(long):: z, c

   sq2=dsqrt(2.0_long)
  
   basis=0.0_long
   do i_star=1,N_star
      kbz=0
      kbz(1:dim)=wave_of_star(:,i_star)
      i_wave=which_wave(kbz(1),kbz(2),kbz(3))
      c=coeff(i_wave)

      if ( present(no_parity) .and. ( no_parity ) ) then

         kfft=G_to_fft(kbz)
         z = kgrid( kfft(1), kfft(2), kfft(3) )

         basis(i_star) = dble(z/c)

      else

         kfft=G_to_fft(kbz)
         if (kfft(1) > ngrid(1)/2) then
            kfft = -kfft
            kfft = G_to_fft(kfft)
            z = conjg( kgrid(kfft(1),kfft(2),kfft(3)) )
         else
            z = kgrid( kfft(1), kfft(2), kfft(3) )
         endif

         select case( star_invert(i_star) )
         case(0)
            basis(i_star) = dble(z/c)
         case(1)
            basis(i_star) = dble(z/c) * sq2
         case(-1)
            basis(i_star) =aimag(z/c) * sq2
         case default
            stop "Illegal star_invert value in grid_to_basis."
         end select
       end if
   end do

   end subroutine kgrid_to_basis
   !==============================================================


   !****p grid_basis_mod/check_symmetry --------------------------
   ! SUBROUTINE 
   !   check_symmetry
   ! PURPOSE
   !   Need to explain purpose and arguments
   ! SOURCE
   !--------------------------------------------------------------   
   subroutine check_symmetry(kgrid,basis,prefix)
   implicit none
   complex(long), intent(IN)         :: kgrid(0:,0:,0:)
   real(long), intent(IN), optional  :: basis(:)
   character(*)                      :: prefix
   !***

   ! local variables
   integer  :: i,j,k
   integer  :: n_dead,n_live
   integer  :: G(3), Gbz(3)
   integer  :: i_wave,i_star

   n_dead=0
   n_live=0
   open(11,file=trim(prefix)//'.dead',status='unknown')
   do k = 0, ngrid(3)-1
      G(3)=k
      do j = 0, ngrid(2)-1
         G(2)=j
         do i = 0, ngrid(1)-1
            G(1)=i
            Gbz=G_to_bz(G) 
            i_wave=which_wave(Gbz(1),Gbz(2),Gbz(3))

            if ( i_wave<1 .OR. i_wave>N_wave ) then
               n_dead = n_dead + 1
               if( i <= ngrid(1)/2 )then
                  write(11,*) Gbz(1:dim),kgrid(i,j,k)
               endif
            else
               n_live = n_live + 1
            endif
         enddo
      enddo
   enddo
   close(11)
   write(6,*)
   if (n_live /= N_wave) write(6,*) "********",prefix," SERIOUS ********"
   write(6,*) 'n_live=',n_live
   write(6,*) 'N_wave=',N_wave
   write(6,*)

   open(12,file=trim(prefix)//'.live',status='unknown')
   G  = 0
   Gbz= 0
   do i = 1, N_star
      write(12,*) 'star=',i,'star_count=',star_count(i)
      do j = star_begin(i), star_end(i)
         Gbz(1:dim) = wave(1:dim, j)
         G = G_to_fft(Gbz)
         if (G(1) > ngrid(1)/2) then
            G = -G
            G = G_to_fft(G)
         end if
         if (present(basis)) then
            write(12,*) Gbz(1:dim),DBLE(kgrid(G(1),G(2),G(3))/coeff(j))-basis(i)
         else
            write(12,*) Gbz(1:dim),DBLE(kgrid(G(1),G(2),G(3))/coeff(j))
         endif
      end do
      write(12,*)
   end do
   close(12)

   end subroutine check_symmetry
   !==============================================================

end module grid_basis_mod
