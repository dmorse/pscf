! fortran_dialect=elf
!-----------------------------------------------------------------------
! PSCF - Polymer Self-Consistent Field Theory
! Copyright (2002-2016) Regents of the University of Minnesota
! contact: David Morse, morse012@umn.edu
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation. A copy of this license is included in
! the LICENSE file in the top-level PSCF directory.
!-----------------------------------------------------------------------
!****m* scf/basis_mod
! MODULE
!    basis_mod - create symmetry-adapted basis functions
!
! PURPOSE
!    Define data structures and procedures involving reciprocal
!    lattice vectors, stars, and symmetry-adapted basis functions
!
! AUTHOR
!
!    David Morse (2002) Original version, for specral algorithm
!    Chris Tyler (2002-2003) various corrections and additions
!    Jian Qin (2005-2007) modified for pseudo-spectral algorithm
!
! SOURCE
!-----------------------------------------------------------------------
module basis_mod
   use const_mod  ! Defines constant long, and variable dim=1, 2, or 3
   use io_mod,        only : output
   use string_mod,    only : int_string
   use version_mod,   only : version_type, output_version
   use unit_cell_mod, only : G_basis, make_G_basis, output_unit_cell
   use group_mod  ! Defines types and operations for crystal groups
   use space_groups_mod, only : space_groups
   use grid_mod   ! Utilities for fft grids
   implicit none

   private

   ! Public procedures
   public:: make_basis, make_dGsq ! needed by scf
   public:: make_waves, make_stars, release_basis
   public:: f_basis
   public:: reorder_basis
   public:: scattering_intensity
   public:: output_waves

   ! Public variables
   public:: group
   public:: N_star
   public:: N_wave, wave, star_of_wave, coeff
   public:: G_max, which_wave
   public:: star_begin, star_end, star_count
   public:: star_cancel, star_invert
   public:: wave_of_star
   public:: sort
   public:: valid_wave

   ! Declaration of public variables
   type(group_type) :: group      ! space group

   integer    :: N_wave           ! # of reciprocal lattice G vectors
   integer    :: N_star           ! # of stars
   integer    :: G_max(3)         ! maximum of values G(1), G(2), G(3)
   integer    :: wave(:,:)        ! (dim, N_wave) list of G vectors
                                  ! components in BZ, grouped by star
   real(long) :: Gsq(:)           ! (N_wave) list of values of |G|^{2}
   complex(long):: coeff(:)       ! (N_wave) coefficient of wave in star
   integer    :: star_of_wave(:)  ! (N_wave) index of star containing
                                  !          this wave
   integer    :: which_wave(:,:,:)! (-G_max(1):G_max(1), ..., ... )
                                  ! wave id listed by components in BZ
   integer    :: star_begin(:)    ! (N_star) index of first wave in star
   integer    :: star_end(:)      ! (N_star) index of last  wave in star
   integer    :: star_count(:)    ! (N_star) # of G vectors in star
   integer    :: wave_of_star(:,:)! (dim,N_star) components of one wave
                                  ! from the star, serves as a label
   logical    :: star_cancel(:)   ! (N_star) see wave_format
   integer    :: star_invert(:)   ! (N_star) see wave_format
   !***

   allocatable :: wave, Gsq, coeff, star_of_wave, which_wave
   allocatable :: star_begin, star_end, star_count, wave_of_star
   allocatable :: star_invert, star_cancel

   ! private parameters
   integer, parameter :: max_star = 48   ! Max. # of G vectors in star
   integer, parameter :: max_list = 512  ! Max # of G vecs in a list

   ! private variables
   logical  :: grid             ! true if PS method is used

   ! Documentation headers for public global variables
   !--------------------------------------------------------------------
   !****v basis_mod/N_wave
   ! VARIABLE
   ! integer N_wave = # of G vectors (reciprocal lattice vectors)
   !*** ----------------------------------------------------------------
   !****v basis_mod/N_star
   ! VARIABLE
   ! integer N_star = # of stars
   !*** ----------------------------------------------------------------
   !****v basis_mod/G_max
   ! VARIABLE
   ! integer G_max(3) = Max values of components of integer waves
   !*** ----------------------------------------------------------------
   !****v basis_mod/wave
   ! VARIABLE
   ! integer wave(:,:)       - allocatable, dimension(dim,N_wave)
   !         wave(1:dim, i)  = integer wavevectors # i, Bravais basis
   !*** ----------------------------------------------------------------
   !****v basis_mod/Gsq
   ! VARIABLE
   ! real(long) Gsq(:)       - allocatable, dimension(N_wave)
   !            Gsq(i)       = value of |G|^{2} for wavevector i
   !*** ----------------------------------------------------------------
   !****v basis_mod/coeff
   ! VARIABLE
   ! complex(long) coeff(:)  - allocatable, dimension(N_wave)
   !               coeff(i)  = complex coefficients of wavevector i
   !                            in symmetry adapted function for star
   !*** ----------------------------------------------------------------
   !****v basis_mod/star_of_wave
   ! VARIABLE
   ! integer star_of_wave(:) - allocatable, dimension(N_wave)
   !         star_of_wave(i) = index of star containing wave i
   !*** ----------------------------------------------------------------
   !****v basis_mod/which_wave
   ! VARIABLE
   ! integer which_wave(G1,G2,G3) = index of wavevector G=(G1,G2,G3)
   !
   ! If wave(i) = (G1,G2,G3), then which_wave(G1,G2,G3) = i
   ! allocatable, dimension( -G_max(1):G_max(1) , ... , ... )
   !*** ----------------------------------------------------------------
   !****v basis_mod/star_begin
   ! VARIABLE
   ! integer star_begin(:)   - allocatable, dimension(N_star)
   !         star_begin(i)   = index of first G in star i
   !*** ----------------------------------------------------------------
   !****v basis_mod/star_end
   ! VARIABLE
   ! integer star_end(:)     - allocatable, dimension(N_star)
   !         star_end(i)     =  index of last G in star i
   !*** ----------------------------------------------------------------
   !****v basis_mod/star_count
   ! VARIABLE
   ! integer star_count(:)   - allocatable, dimension(N_star)
   !         star_count(i)   = # of G vectors in star i
   !*** ----------------------------------------------------------------
   !****v basis_mod/star_invert
   ! VARIABLE
   ! integer star_invert(:)  - allocatable, dimension(N_star)
   !         star_invert(i)  = invert flag for star i, values (0,1,-1)
   !                           See discussion of wave_format
   !*** ----------------------------------------------------------------
   !****v basis_mod/wave_of_star
   ! VARIABLE
   ! integer wave_of_star(:,:) - dimension(dim, N_star)
   !         wave_of_star(:,i) = G vector from star i, serves as a label
   !                             See discussion of wave_format
   !*** ----------------------------------------------------------------
   !****v basis_mod/star_cancel
   ! VARIABLE
   ! logical  star_cancel(:) - allocatable, dimension(N_star)
   !          star_cancel(i) = true if star i is cancelled
   !*** ----------------------------------------------------------------


   !--------------------------------------------------------------------
   !****c* basis_mod/wave_format
   ! COMMENT
   !
   !     i) G vectors are listed in the array wave in non-decreasing
   !        order of |G|^2 (for some G_basis passed to make_waves),
   !        and (after output from make_stars) are grouped by stars,
   !        each of which forms contiguous block of wave vectors.
   !        Components of each wave are image within the first BZ,
   !        the aliased image that minimizes |G|^{2}.
   !
   !    ii) Cancelled waves and stars are included in lists of waves
   !        and stars only if logical variable keep_cancel is passed to
   !        make_waves with a value .true.  If keep_cancel is true, 
   !        then cancelled waves are asigned coeff= 0.0 and cancelled
   !        stars are assigned star_cancel = true. If keep_cancel is
   !        false, then cancelled waves and stars are excluded from
   !        these arrays, and star_cancel is false for all stars that
   !        are retained.
   !
   !   iii) Within each star, vectors are listed in "descending" order,
   !        if vectors in wave are read as dim=1, 2, or 3 digit numbers,
   !        with more significant digits on the left. For example the
   !        [111] star of a cubic structure is listed in the order:
   !
   !           1  1  1  (first)
   !           1  1 -1
   !           1 -1  1
   !           1 -1 -1
   !          -1  1  1
   !          -1  1 -1
   !          -1 -1  1
   !          -1 -1 -1  (last)
   !
   !    iv) If star number i is closed under inversion, so that the
   !        negative -G of each wave G in the star is also an element
   !        of that star, then we assign star_invert(i) = 0. For
   !        space groups that contain inversion symmetry (i.e., for
   !        centro-symmetric groups), all stars must be closed under
   !        inversion.
   !
   !     v) For non-centrosymmetric space groups, pairs of stars may be
   !        related by inversion, so that the -G of each vector in star
   !        number i is found not in star i but in a second star. Such
   !        pairs of stars that are related by inversion are always
   !        listed consecutively. The first star in each such a pair is
   !        assigned star_invert(i) = 1, and the second is assigned
   !        star_invert(i+1) = -1.
   !
   !    vi) The combination of conventions (iii) and (v) implies that
   !        the inverse of the first vector in a star with
   !        star_invert = 1 is the  last vector in the next star, which
   !        is the partner star with star_invert = -1. More generally,
   !        the inverse of the ith member of the first star in such a
   !        pair is the ith to last vector of the partner star.
   !
   !   vii) wave_of_star is the first wave in a star with star_invert=0
   !        or 1, and the last wave in a star with star_invert=-1
   !
   !  viii) The elements of the array coeff() are complex numbers that
   !        have the same absolute magnitude for waves in the same star,
   !        with absolute magnitude 1.0/sqrt(star_count) for a start
   !        containing star_count distinct wavevectors. Each star 
   !        generates an associated normalized basis function, referred
   !        to as a star function, given by a linear of complex plane 
   !        waves exp(i G.r) in which each such plane wave is multipled 
   !        by the corresponding element of array coeff.
   !
   !    ix) For any star that is closed under inversion (start_invert=0), 
   !        values of coeff associated with any vector G and its negation 
   !        -G are complex conjugates, so that the corresponding star
   !        function is a real function of position. For centro-symmetric 
   !        groups, all elements of coeff are positive and real, and all
   !        star functions are real, functions of position with even 
   !        parity.
   !
   !        For non-centrosymmetric groups and pairs of stars that are
   !        related by inverse, the coefficients of the first wave of
   !        the first star (star_invert = 1) and the last wave of the 
   !        second star (star_invert = -1) are chosen to be real and 
   !        positive. This is sufficient to guarantee that if G is a 
   !        vector in the first star and -G is a corresponding vector 
   !        in the second, the coefficients of G and -G are complex 
   !        conjugates. It also guarantees that the star functions
   !        generated by the first and second stars in such an 
   !        inversion-related pair are complex conjugates of one
   !        another.
   !
   !    x)  The N_star stars are used to construct N_star real basis 
   !        functions, denoted here by f_1(r), ... , f_{N_star}(r).
   !        These real basis functions are defined and indexed as 
   !        follows:
   !
   !        For each star that is closed under inversion, with star
   !        index j and star_invert(j) =- 0, we create a single real 
   !        star function phi_j(r) and use this as a basis function
   !
   !           f_j(r) = phi_j(r)  
   !
   !        that is assigned a basis function index equal to the star
   !        index j.
   !
   !        For each pairs of stars that are related by inversion, 
   !        with star indices j and j+1 and associated star functions
   !        phi_{j}(r) and phi_{j+1}(r), we construct one "cosine-like"
   !        basis function, with basis index j, given by the sum
   !
   !           f_j(r) = [ phi_{j}(r) + phi_{j+1}(r) ]/sqrt(2) ,
   !
   !        and a second "sine-like" basis function, with a basis 
   !        index j+1, given by
   !
   !           f_{j+1}(r) = -i[ phi_{j}(r) - phi_{j+1}(r) ]/sqrt(2)
   !
   !        where -i = sqrt(-1). 
   !
   !***
   !--------------------------------------------------------------------

contains

   !--------------------------------------------------------------------
   !****p basis_mod/release_basis
   !
   ! SUBROUTINE
   !      release_basis()
   ! PURPOSE
   !      Release the memory allocated by constructing basis information
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine release_basis
   !***

   integer     :: error      ! error index

   if (allocated(wave)) deallocate(wave,stat=error)
   if (error /= 0) stop "wave deallocation error!"

   if (allocated(Gsq)) deallocate(Gsq,stat=error)
   if (error /= 0) stop "Gsq deallocation error!"

   if (allocated(coeff)) deallocate(coeff,stat=error)
   if (error /= 0) stop "coeff deallocation error!"

   if (allocated(star_of_wave)) deallocate(star_of_wave,stat=error)
   if (error /= 0) stop "star_of_wave deallocation error!"

   if (allocated(which_wave)) deallocate(which_wave,stat=error)
   if (error /= 0) stop "which_wave deallocation error!"

   if (allocated(star_begin)) deallocate(star_begin,stat=error)
   if (error /= 0) stop "star_begin deallocation error!"

   if (allocated(star_end)) deallocate(star_end,stat=error)
   if (error /= 0) stop "star_end deallocation error!"

   if (allocated(star_count)) deallocate(star_count,stat=error)
   if (error /= 0) stop "star_count deallocation error!"

   if (allocated(wave_of_star)) deallocate(wave_of_star,stat=error)
   if (error /= 0) stop "wave_of_star deallocation error!"

   if (allocated(star_invert)) deallocate(star_invert,stat=error)
   if (error /= 0) stop "star_invert deallocation error!"

   if (allocated(star_cancel)) deallocate(star_cancel,stat=error)
   if (error /= 0) stop "star_cancel deallocation error!"

   end subroutine release_basis
   !====================================================================


   !--------------------------------------------------------------------
   !****p basis_mod/make_basis
   ! SUBROUTINE
   !
   !      make_basis(R_basis, G_basis, group_name, Gabs_max)
   !
   ! PURPOSE
   !
   ! Generates information about basis functions needed by scf code:
   !   i) Calls make_G_basis to generate reciprocal lattice basis vector
   !      array G_basis, using R_basis (Bravais basis vectors) as an
   !      input
   !  ii) Calls space_groups to create a group identified by the string
   !      identifier group_name
   !      function space_groups
   ! iii) Calls make_waves and make_stars to generate all waves
   !      and stars with |G| < Gabs_max, excluding cancelled stars
   !
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine make_basis(         &
                      R_basis,    &! real lattice basis vectors
                      G_basis,    &! reciprocal lattice basis vectors
                      group_name, &! file for generators of space group
                      N_grids,    &! maximum value of |G|
                      grid_flag   &! true if PS method used
                      )

   logical,      intent(IN)    :: grid_flag
   real(long),   intent(IN)    :: R_basis(:,:)   ! real lattice basis
   real(long),   intent(OUT)   :: G_basis(:,:)   ! Reciprocal basis
   character(*), intent(INOUT) :: group_name     ! file for group
   integer,      intent(IN)    :: N_grids(3)     ! # of grid points
   !***

   ! Local variables
   integer            :: i,j,k,l
   logical            :: keep_cancel
   character(len = 3) :: Nchar

   ! Set module variable grid equal to input variable grid_flag
   grid = grid_flag

   ! Make reciprocal lattice basis G_basis
   call make_G_basis(R_basis, G_basis)

   ! Select a space group by name, return as variable group
   ! First look for the identifier group_name in the hard-coded database,
   ! otherwise try to read elements from a file with name group_name

   call space_groups(group_name, group)

   ! Check that the group is complete, and complete it if necessary
   call make_group(group, R_basis, G_basis)

   ! Set policy: All cancelled stars and waves are discarded in make_waves
   keep_cancel = .false.

   ! Make all waves, discarding cancelled waves
   ! On return, arrays wave and Gsq are filled with all non-cancelled
   ! waves listed in non-descending order, starting with the zero wave
   call make_waves(         &
                G_basis,    &!  (dim,dim), basis for reciprocal lattice
                R_basis,    &!  (dim,dim), basis for Bravais lattice
                N_grids,    &!  (3), # of grid points
                keep_cancel &!  if true, keep cancelled waves
                  )

   ! Group waves into stars related by space group symmetries
   call make_stars

   if (grid) call make_ksq(G_basis)  ! true if PS method used

   end subroutine make_basis
   !====================================================================


   !--------------------------------------------------------------------
   !****ip basis_mod/make_waves
   ! SUBROUTINE
   !    make_waves(G_basis, R_basis, N_grids, keep_cancel)
   ! PURPOSE
   !    Subroutine fills arrays wave, Gsq with lists of all reciprocal
   !    lattice vectors in FFT grid, sort in ascending order of Gsq.
   ! ARGUMENTS
   !    Inputs
   !    real(long)  G_basis(:,:) - basis for reciprocal lattice
   !    real(long)  R_basis(:,:) - basis for Bravais lattice
   !    integer     N_grids(3)   - # of grid points in x, y, z
   !    logical     keep_cancel  - if true, keep cancelled waves
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine make_waves(   &
                G_basis,    &!  (dim,dim), basis for reciprocal lattice
                R_basis,    &!  (dim,dim), basis for Bravais lattice
                N_grids,    &!  (3), # of grid points
                keep_cancel &!  if true, keep cancelled waves
                        )
   real(long), intent(IN)  :: G_basis(:,:), R_basis(:,:)
   integer,    intent(IN)  :: N_grids(3)
   logical                 :: keep_cancel
   !***

   ! local variables
   real(long) :: Gabs_max
   real(long) :: Gsq_new, Gsq_max, G_vec(3), twopi, a_abs
   integer    :: i1, i2, i3, i_vec(3), i, j, k
   logical    :: cancel

   ! Temporary arrays to store index and values of wave, and Gsq
   integer,    allocatable, dimension(:)   :: index
   integer,    allocatable, dimension(:,:) :: wave_temp

   ! Compute Gabs_max = maximum of |G| on the FFT grid, shifted to BZ
   Gabs_max = max_Gabs(G_basis)  ! function grid_mod/max_Gabs

   ! Note: Computing Gabs_max from knowledge of the FFT grid presumes we
   ! are using the pseudo-spectral method, rather than a spectral method.
   ! Blocks of code intended only to support the spectral method have 
   ! thus become obsolete. These blocks are also always bypassed, because
   ! make_basis is now always called with grid_flag = .TRUE. 
   ! TODO: Clean up or remove code vestigial code for the spectral 
   ! method, then remove grid_flag and grid.

   ! Calculate array G_max(dim) : maximum wavevector indices in the BZ.
   ! The components of G_Max are used as dimensions of which_wave
   Gsq_max = Gabs_max*Gabs_max
   twopi   = 4.0_long*acos(0.0_long)
   G_max   = 0
   do i=1, dim
      a_abs = sqrt( R_basis(i,:) .dot. R_basis(i,:) )
      G_max(i) = int( Gabs_max * a_abs / twopi ) + 1
   enddo

   ! Calculate number of plane waves
   N_wave = 0

   if ( .not. grid ) then

      ! Use spectral method convention, keep all waves |G|^{2} < Gsq_max.
      ! This is an obsolete block that is bypassed because grid = .TRUE.

      do i1= -G_max(1), G_max(1)
         i_vec(1) = i1
         if ( dim == 1) then
            G_vec= i_vec .dot. G_basis
            Gsq_new = G_vec .dot. G_vec
            if ( Gsq_new <= Gsq_max ) then
               cancel = vector_cancel(i_vec)
               if ((.not.cancel).or.keep_cancel) then
                  N_wave = N_wave + 1
               endif
            endif
         else
            do i2= -G_max(2), G_max(2)
               i_vec(2) = i2
               if (dim == 2) then
                  G_vec= i_vec .dot. G_basis
                  Gsq_new = G_vec .dot. G_vec
                  if ( Gsq_new <= Gsq_max ) then
                     cancel = vector_cancel(i_vec)
                     if ((.not.cancel).or.keep_cancel) then
                        N_wave = N_wave + 1
                     endif
                  endif
               else
                  do i3= -G_max(3), G_max(3)
                     i_vec(3) = i3
                     G_vec    = i_vec .dot. G_basis
                     Gsq_new  = G_vec .dot. G_vec
                     ! write(6,*) i_vec(1),i_vec(2),i_vec(3), Gsq_new
                     if (Gsq_new <= Gsq_max) then
                        cancel = vector_cancel(i_vec)
                        if ((.not.cancel).or.keep_cancel) then
                           N_wave = N_wave + 1
                        endif
                     endif
                  enddo
               endif
            enddo
         endif
      enddo

   else   ! use shifted FFT grid to generate wave

      do i1 = 0, N_grids(1) - 1
         i_vec(1) = i1
         do i2 = 0, N_grids(2) - 1
            i_vec(2)=i2
            do i3 = 0, N_grids(3) - 1
               i_vec(3) = i3
               i_vec = G_to_bz(i_vec)
               cancel = vector_cancel(i_vec)
               if ((.not.cancel).or.keep_cancel) N_wave = N_wave+1
            enddo
         enddo
      enddo

   end if !  basis if-then

   ! Allocate module arrays
   if (allocated(wave)) deallocate(wave)
   allocate(wave(dim,N_wave), stat = j)
   if (j.ne.0) stop 'make_basis allocate wave error'

   if (allocated(Gsq)) deallocate(Gsq)
   allocate(Gsq(N_wave), stat = j)
   if (j.ne.0) stop 'make_basis allocate N_wave error'

   if (allocated(star_of_wave)) deallocate(star_of_wave)
   allocate(star_of_wave(N_wave), stat = j)
   if (j.ne.0) stop 'Make_basis allocate star_of_wave error'

   if (allocated(which_wave)) deallocate(which_wave)
   allocate(which_wave&
        (-G_max(1):G_max(1),-G_max(2):G_max(2),-G_max(3):G_max(3)),stat=j)
   if (j.ne.0) stop 'Make_basis allocate which_wave error'

   if (allocated(coeff)) deallocate(coeff)
   allocate(coeff(N_wave),stat=j)
   if (j.ne.0) stop 'Make_basis coef allocate error)'

   ! Allocate temporary arrays
   allocate(wave_temp(dim,N_wave))
   allocate(index(N_wave))

   ! Construct wave and Gsq arrays
   if ( .not. grid ) then

      ! Use spectral method convention.
      ! Note: This is an obsolete block that is always bypassed.

      j     = 0
      i_vec = 0
      do i1= -G_max(1), G_max(1)
         i_vec(1) = i1
         if (dim == 1) then
            G_vec= i_vec .dot. G_basis
            Gsq_new = G_vec .dot. G_vec
            if (Gsq_new <= Gsq_max) then
               cancel = vector_cancel(i_vec)
               if ((.not.cancel).or.keep_cancel) then
                  j = j + 1
                  Gsq(j) = Gsq_new
                  ! Set wave(:,j) = i_vec
                  call assign_ivec(wave(:,j),i_vec)
               endif
            endif
         else
            do i2= -G_max(2), G_max(2)
               i_vec(2) = i2
               if (dim == 2) then
                  G_vec= i_vec .dot. G_basis
                  Gsq_new = G_vec .dot. G_vec
                  if (Gsq_new <= Gsq_max) then
                     cancel = vector_cancel(i_vec)
                     if ((.not.cancel).or.keep_cancel) then
                        j = j + 1
                        Gsq(j) = Gsq_new
                        ! Set wave(:,j) = i_vec
                        call assign_ivec(wave(:,j),i_vec)
                     endif
                  endif
               else
                  do i3= -G_max(3), G_max(3)
                     i_vec(3) = i3
                     G_vec= i_vec .dot. G_basis
                     Gsq_new = G_vec .dot. G_vec
                     if (Gsq_new <= Gsq_max) then
                        cancel = vector_cancel(i_vec)
                        if ((.not.cancel).or.keep_cancel) then
                           j = j + 1
                           Gsq(j) = Gsq_new
                           ! Set wave(:,j) = i_vec
                           call assign_ivec(wave(:,j),i_vec)
                        endif
                     endif
                  enddo
               endif
            enddo
         endif
      enddo

   else  ! use shifted fft grid to generate wave

      j = 0
      do i1=0, N_grids(1) - 1
         i_vec(1)=i1
         do i2=0, N_grids(2) - 1
            i_vec(2)=i2
            do i3=0, N_grids(3) - 1
               i_vec(3) = i3
               i_vec = G_to_bz(i_vec)
               cancel = vector_cancel(i_vec)
               if ( (.not.cancel) .or. keep_cancel ) then
                  j = j + 1
                  Gsq(j) = norm(i_vec, G_basis)
                  ! Set wave(:,j) = i_vec
                  call assign_ivec(wave(:,j),i_vec)
               endif
            enddo
         enddo
      enddo

   end if  ! grid if--then

   ! Create ordered index
   do i=1, N_wave
      index(i)=i
   enddo

   ! Sort Gsq in ascending order, and reorder index accordingly
   call sort(N_wave, Gsq, index)

   ! Use reordered index to re-order wave array
   do i=1, N_wave
      wave_temp(:,i) = wave(:,index(i))
   enddo
   wave = wave_temp

   ! Deallocate temporary arrays
   deallocate(wave_temp)
   deallocate(index)

   end subroutine make_waves
   !====================================================================


   !--------------------------------------------------------------------
   !****ip basis_mod/make_stars
   ! SUBROUTINE
   !    make_stars
   ! PURPOSE
   !    Given output of make_waves, the subroutine divides the list of
   !    G vectors first into "lists" of G vectors with equal values of
   !    Gsq and then into "stars" of vectors that are related by the
   !    symmetry operations of the group.
   !
   ! On input:
   !   The arrays wave(dim,N_wave) and Gsq(N_wave) contain lists
   !   of all waves and corresponding values of |G|^{2}, ordered by
   !   nondecreasing values of Gsq, but not yet grouped into stars.
   !
   ! On output:
   !    1)  Arrays wave and Gsq are reordered so that members of each
   !        star are listed consecutively
   !    2)  star_of_wave(i) contains index of star for wave i
   !    3)  coeff(i) contains complex coefficient for wave i
   !    4)  wave(which_wave(G(1),G(2),G(3))) = G for any wave in list
   !    5)  star_begin(i) and star_end(i) contain wave indices for first
   !        and last wave in star i, and star_count(i) = # waves in star
   !    6)  star_invert(i) contains an integer flag giving information
   !        about effect of inversion upon star i (see below)
   !    7)  star_cancel(i) = .true. if star i is cancelled.
   !        (Relevant only if make_waves left cancelled waves in lists)
   !--------------------------------------------------------------------
   ! Stars related by inversion, and the value of star_invert:
   !
   !    1) Stars that are closed under inversion (i.e., that contain
   !       the inverse -G of every vector G in the star) have
   !
   !          star_invert(i) = 0
   !
   !       For centrosymmetric groups (i.e., groups that include
   !       inversion symmetry) all stars will be closed under inversion
   !
   !    2) Pairs of stars that are related by inversion are list
   !       consecutively, with:
   !
   !          star_invert(i)   =  1 for for the 1st star in the pair
   !          star_invert(i+1) = -1 for the 2nd
   !
   !       Such inversion-related pairs of stars can exist only for
   !       non-centrosymmetric groups
   !
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine make_stars
   !***

   ! Parameter
   real(long), parameter :: epsilon  = 1.0E-7_long

   ! Local variables
   integer    :: list(3, max_list), star(3, max_star)
   real(long) :: star_phase(max_star)
   real(long) :: Gsq_max, twopi
   integer    :: list_index(max_list), star_index(max_star)
   integer    :: G1(3), G2(3)
   integer    :: first_list, last_list, list_N, star_N
   integer    :: i, j, k, l, i_index, i_star, root, invert_flag
   integer    :: i1, i2
   complex(long) :: c_norm, c1, c2, d
   logical    :: new_list, cancel

   ! Arrays to store index and temporary values of wave, and Gsq
   integer     :: index(:)            ! (N_wave)
   integer     :: wave_temp(:,:)      ! (dim, N_wave)
   real(long)  :: Gsq_temp(:)         ! (N_wave)
   allocatable :: index, wave_temp, Gsq_temp
   integer     :: info

   allocate(index(N_wave), stat=info)
   if ( info /= 0 ) then
      write(6,*) "index(N_wave) allocation error: N_wave = ", N_wave
      stop
   end if

   twopi       = 4.0_long*acos(0.0_long)
   which_wave  = 0

   ! Check that the first vector in wave & Gsq is the zero vector
   do k = 1, dim
      if (wave(k, 1) /= 0) then
         write(6,*) "Error: wave(k, 1) neq 0, for k = ", k
         stop
      endif
   enddo
   if (Gsq(1) > epsilon) then
      write(6,*) "Error: Gsq(1) > 0.0"
      stop
   endif

   ! Set initial values for variables used in loop over wavevectors
   first_list   = 1          ! first list begins with index 1
   Gsq_max      = 0.0_long   ! first vector has |G| = 0
   i_index      = 0          ! new index of this wave, after reordering
   i_star       = 0          ! index of this star
   invert_flag  = 0          ! invert_flag value of the previous star

   !--------------------------------------------------------------------
   ! Loop over wavevectors:     
   ! Define: a "list" is collection of vectors of equal magnitude |G|.
   ! Identify lists by finding where value of |G|^{2} changes.        
   ! This loop creates arrays index, star_of_wave, which_wave, and coeff.
   ! The array index is used later to reorder wave and Gsq.
   !--------------------------------------------------------------------
   do i = 2, N_wave + 1

      !----------------------------------------------------------------
      ! Test for end of "list", due to increase in Gsq or end of array.
      ! If end of a list is found, set last_list and new_list == .true. 
      !----------------------------------------------------------------
      if (i > N_wave) then
         new_list = .false.
         if (first_list == N_wave) then
            new_list = .true.
            last_list = N_wave
         endif
      else
         if (Gsq(i) < 0) then
            write(6,*) 'Error in make_stars: Gsq < 0'
            stop
         endif
         if (Gsq(i) < Gsq_max - epsilon) then
            write(6,*) 'Error in make_stars: unordered Gsq array'
            stop
         endif
         new_list = .false.
         if (Gsq(i) > Gsq_max + epsilon) then
            Gsq_max = Gsq(i)
            last_list = i - 1
            new_list  = .true.
         else if (i == N_wave) then
            last_list = N_wave
            new_list  = .true.
         end if
      end if

      !--------------------------------------------------------------
      ! If new_list == true (i.e., if a new list was completed):
      !   last_list is index of last vector in the completed list
      !   first_list is index of first vector in the completed list
      ! If a change in |G|^2 was detected, then Gsq_max is the
      !   value of |G|^2 for the next list. 
      ! Note: 
      !   first_list will be reset to i after processing this list
      !--------------------------------------------------------------

      ! Begin processing the completed list, if any
      if (new_list) then

        ! Compute number list_N of vectors in list
        list_N = last_list - first_list + 1

        ! Check that list_N <= max_list
        if (list_N > max_list) then
           write(6,*) 'Error: list_N > max_list in make_stars'
           write(6,*) last_list, first_list, max_list, N_wave
           write(6,'(3i5)') list(:,first_list-1)
           write(6,'(3i5)') list(:,first_list:last_list)
           stop
        endif

        ! Initialize local arrays list and list_index
        list = 0
        do j = first_list, last_list
           k = j - first_list + 1
           call assign_ivec(list(:,k), wave(:,j))
           list_index(k) = j
        enddo

        !----------------------------------------------------------
        ! For index k = 1, ... , list_N
        !   list(:,k) contains a wave in the list
        !   list_index(k) is "old" index of that wave in wave, Gsq
        !----------------------------------------------------------

        ! Sort list and list_index by Miller indices (highest first)
        call G_sort(list_N, list, list_index)

        ! Output list_N, first and last waves (diagnostic)
        ! write(6,*) "  "
        ! write(6,*) "Processing list of ", list_N
        ! write(6,*) "list(1)      = ", list(:,1)
        ! write(6,*) "list(list_N) = ", list(:,list_N)
        ! write(6,*) "  "

        !----------------------------------------------------------
        ! Begin loop over stars within list:
        ! Pass list to routine get_star, which outputs a "star"
        !    constructed from vector # root in input list, and
        !    remaining "list" containing vectors from input list
        !    that are not in this star
        ! Add the resulting star to global arrays.
        ! Repeat if the remaining "list" is not empty.
        !----------------------------------------------------------
        root = 1
        do while (list_N > 0)
           call get_star(root, &
                      list, list_index, list_N, &
                      star, star_index, star_N, star_phase)
           ! Check that star is not empty
           if (star_N == 0) then
              write(6,*) ' star_N == 0 in make_star'
              stop
           endif

           ! Update i_star (index of new star)
           i_star = i_star  + 1

           ! Sort vectors in star and star_phase arrays
           call G_sort(star_N, star, star_index, reorder_real=star_phase)

           !------------------------------------------------------------
           ! At this point:
           !   array star holds an ordered star
           !   array star_phase holds phases in rads, relative to root
           !------------------------------------------------------------

           ! Sort vectors in remaining list
           call G_sort(list_N, list, list_index)

           ! Output contents of star & remaining list (diagnostic)
           ! write(6,*)
           ! write(6,*) 'i_star =', i_star
           ! write(6,*) 'star_N =', star_N
           ! write(6,*) 'star ='
           ! do j=1, star_N
           !   write(6,FMT='(2I5,3X,3I5)') j, star_index(j), star(:,j)
           ! enddo
           ! write(6,*)
           ! if (list_N > 0) then
           !   write(6,*) 'list_N =', list_N
           !   write(6,*) 'remaining list='
           !   do j=1, list_N
           !      write(6,FMT='(2I5,3X,3I5)') j, list_index(j), list(:,j)
           !   enddo
           ! endif

           ! Check for cancellation of first vector in star
           G1 = star(:,1)
           cancel = vector_cancel(G1)

           ! Check that entire star is either cancelled or not cancelled
           do j=1, star_N
              if (cancel .NEQV. vector_cancel(star(:,j))) then
                 write(6,*) 'Inconsistent results for cancellation:'
                 write(6,*) 'star(:,1), cancel=', G1, cancel
                 cancel = vector_cancel(star(:,j))
                 write(6,*) 'j, star(:,j), cancel=',star(:,j),cancel
                 stop
              endif
           enddo

           !------------------------------------------------------
           ! Compute normalization for un-cancelled stars:
           ! Coefficient of last wave in star is positive real if
           !    invert_flag == 1 for previous star, so that the
           !    current star is the 2nd star in an inversion pair
           ! Otherwise, star_invert = 0 or 1 for current star,
           !    and the coefficient of 1st wave is positive real
           ! In either case, abs(c_norm) = sqrt(star_N)
           !------------------------------------------------------
           if (.not.cancel) then
              if (invert_flag == 1) then
                  c_norm = exp((0.0_long,1.0_long)*star_phase(star_N))
               else
                  c_norm = exp((0.0_long,1.0_long)*star_phase(1))
              endif
              c_norm = c_norm * sqrt(real(star_N))
           endif

           !---------------------------------------------------------
           ! Add star to global arrays that are indexed by wave:
           !   Set block of values of index to old wave indices
           !   Add values to star_of_wave, which_wave, & coeff
           !---------------------------------------------------------
           do j=1, star_N

              ! Increment i_index (new global index of this wave)
              i_index = i_index + 1

              ! Set next element of index to old global index of wave
              index(i_index) = star_index(j)

              ! Set next element of star_of_wave to star index
              star_of_wave(i_index)  = i_star

              ! Set next element of which_wave
              G1 = star(:,j)
              which_wave(G1(1),G1(2),G1(3)) = i_index

              ! Set next element of coeff
              if (cancel) then 
                 coeff(i_index) = (0.0_long,0.0_long)
              else
                 coeff(i_index) = &
                     exp( (0.0_long,1.0_long)*star_phase(j) )/c_norm
              endif

           end do
           !---------------------------------------------------------
           ! Note: coeff is set postive real for the wave of star,
           ! i.e., for the first wave for closed star or first star 
           ! in pair related by inversion, or for the last wave in 
           ! the second star in a pair.
           !---------------------------------------------------------
        
           !---------------------------------------------------------
           ! Set invert_flag for this star & root id for next star:
           !
           ! Look for negation of the first vector of this star:
           !   - In this star or, if not in this star, them
           !   - In the remaining list
           !
           ! If star is closed under inversion (invert_flag = 0),
           !    or is second star of pair related by inversion
           !    (invert_flag = -1), set root = 1, so the first
           !    vector of the remaining list will be used as the
           !    root vector of the next star.
           !
           ! If star is not closed under inversion, and is the 
           !    1st star in pair related by inversion (invert_flag=1), 
           !    set root to the id within the remaining list of the 
           !    negation -G1 of the root vector G1 of the current 
           !    star, so that -G1 is used as the root of the next
           !    star.
           !---------------------------------------------------------
           G1 = star(:,1)
           G2 = -G1
           if (grid) G2 = G_to_bz(G2)

           if (in_list(G2, star, star_N) /= 0) then
              ! This star is closed
              if (invert_flag == 1) then
                 ! Expect second star in pair, shouldn't be closed'
                 write(6,*) &
                    'Error: Expect invert_flag = -1, but -G1 is in star'
                 stop
              else
                 invert_flag = 0
                 root        = 1
                 !write(6,*) 'Parity flag =0, -G1 in star'
              endif
           else if (invert_flag == 1) then
              ! This is the second open star in a pair
              invert_flag = -1
              root        =  1
              !write(6,*) 'Parity flag = -1'
           else
              ! This is the first open star in a pair
              invert_flag = 1
              ! Find the index of -G1 in the remaining list
              root = in_list(G2, list, list_N)
              if (root == 0) then
                 write(6,*) 'Error: -G1 in neither in star nor list'
                 stop
              endif
           endif
           !---------------------------------------------------------
           ! Upon exit of above section:
           !    invert_flag =  0 -> this star closed under inversion
           !    invert_flag =  1 -> this star is 1st in pair related
           !                        by inversion
           !    invert_flag = -1 -> this star is 2nd in pair related
           !                        by inversion
           !---------------------------------------------------------

        end do ! while (list_N > 0)
        ! End loop over stars within list of G's of equal |G|

        ! Set first member of next list = i
        first_list = i

        ! Special case - completed processing of all waves
        if (last_list == N_wave) then
           first_list = 0
        endif

      endif
      ! End processing of a complete list

   enddo
   ! End loop over all wavevectors

   ! Set N_star = total number of stars = final value of i_star
   N_star = i_star

   ! Check that all wavevectors are accounted for
   if (i_index /= N_wave) then
      write(6,*) 'Error: i_index /= N_wave after all stars made'
      write(6,*) 'i_index=', i_index
      write(6,*) 'N_wave =', N_wave
      stop
   endif

   ! Allocate work arrays wave_temp and Gsq_temp 
   allocate(wave_temp(dim,N_wave), stat=info)
   if ( info /= 0 ) then
      write(6,*) "wave_temp(dim,N_wave) allocation error: N_wave =", N_wave
      stop
   end if
   allocate(Gsq_temp(N_wave), stat=info)
   if ( info /= 0 ) then
      write(6,*) "Gsq_temp(N_wave) allocation error: N_wave =", N_wave
      stop
   end if

   ! Use index array to re-order wave and Gsq arrays
   wave_temp = wave
   Gsq_temp  = Gsq
   do i=1, N_wave
      j = index(i)
      wave(:,i) = wave_temp(:, j)
      Gsq(i)= Gsq_temp(j)
   enddo

   ! De-allocate temporary arrays used to re-order wave and Gsq
   if (allocated(index)) deallocate(index)
   if (allocated(wave_temp)) deallocate(wave_temp)
   if (allocated(Gsq_temp)) deallocate(Gsq_temp)

   ! Check consistency of wave and which_wave
   do i=1, N_wave
      call assign_ivec(G1, wave(:,i))
      j = which_wave(G1(1), G1(2), G1(3))
      if (i /= j) then
         write(6,*) &
              'Error: Inconsistent wave & which_wave in make_stars'
         write(6,*) "wave(i) = ", G1
         write(6,*) "i = ", i
         write(6,*) "which_wave(G) = ", j
         stop
      endif
   enddo

   !-------------------------------------------------------------------
   ! At this point:
   !   wave, Gsq, star_of_wave, coeff and which_wave are complete 
   !   wave, Gsq, star_of_wave, and coeff use the same indexing,
   !      in which waves in the same star are in a single block
   !   which_wave is consistent with wave
   !-------------------------------------------------------------------

   ! Allocate arrays that are indexed by star id

   if (allocated(star_begin)) deallocate(star_begin)
   allocate(star_begin(N_star))

   if ( allocated(star_end)) deallocate(star_end)
   allocate(star_end(N_star))

   if (allocated(star_count)) deallocate(star_count)
   allocate(star_count(N_star))

   if (allocated(star_invert)) deallocate(star_invert)
   allocate(star_invert(N_star))

   if (allocated(wave_of_star)) deallocate(wave_of_star)
   allocate(wave_of_star(dim, N_star))

   if (allocated(star_cancel)) deallocate(star_cancel)
   allocate(star_cancel(N_star))

   !----------------------------------------------------------------
   ! Loop over waves to construct star_begin, star_end, star_count
   ! by detecting changes in star_of_wave(i) with increasing i
   !----------------------------------------------------------------
   i_star = star_of_wave(1)
   star_begin(i_star) = 1
   do i=2, N_wave
      if (i_star /= star_of_wave(i)) then
         star_end(i_star)   = i-1
         star_count(i_star) = star_end(i_star) - star_begin(i_star) + 1
         i_star = star_of_wave(i)
         star_begin(i_star) = i
         if (i == N_wave) then
            star_end(i_star)   = i
            star_count(i_star) = 1
         endif
      else if (i == N_wave) then
         star_end(i_star)   = N_wave
         star_count(i_star) = star_end(i_star) - star_begin(i_star) + 1
      endif
   enddo

   !---------------------------------------------------------------
   ! Check consistency of star_begin, star_end, star_count
   !---------------------------------------------------------------
   i_star = 0
   do i=1, N_star
      i_star = i_star + star_count(i)
      if (star_count(i) /= star_end(i) - star_begin(i) + 1) then
         write(6,*) 'Error in make_stars: Invalid star_count, i=',i
         stop
      endif
      if (i > 1) then
         if (star_begin(i) /= star_end(i-1) + 1) then
            write(6,*) 'Error in make_stars: begin(i) /= end(i-1)+1'
            stop
         endif
      endif
   enddo
   if (i_star /= N_wave) then
      write(6,*) 'Error in make_stars: sum(star_count) /= N_wave'
      stop
   endif
   if (star_end(N_star) /= N_wave) then
      write(6,*) 'Error in make_stars: star_end(N_star) /= N_wave'
      stop
   endif

   !---------------------------------------------------------------
   ! Loop over stars to make star_cancel, star_invert, wave_of_star
   !---------------------------------------------------------------
   invert_flag = 0
   do i = 1, N_star   ! loop over stars

      ! Fill array star with corresponding block of array wave
      star = 0
      do j=1, star_count(i)
         k = star_begin(i) + j - 1
         call assign_ivec(star(:,j), wave(:,k))
      enddo

      ! Set G1 (root) to first wave in this star, set G2 = -G1
      G1 = star(:,1)
      G2 = -G1
      if (grid) G2 = G_to_bz(G2)

      i1 = star_begin(i)
      if (i1 /= which_wave(G1(1), G1(2), G1(3))) then
         write(6,*) 'Error: Inconsistent id i1, first wave in star'
      endif

      ! Assign star_cancel(i)
      star_cancel(i) = vector_cancel(G1)

      ! Search for negation of first wave in this star
      k = in_list(G2, star, star_count(i))
      if (k /= 0) then

         if (invert_flag == 1) then

            ! Expect second star in pair, shouldn't be closed'
            write(6,*) 'Error: invert_flag = 2, -G1 is in star'
            stop

         else ! closed star

            invert_flag    = 0
            star_invert(i) = 0
            wave_of_star(:,i) = wave(:,star_begin(i))

            ! Identify global index i2 of -G2 within wave
            i2 = i1 + k - 1
            if (i2 /= which_wave(G2(1), G2(2), G2(3))) then
               write(6,*) 'Error: Inconsistent id i2, last wave in star'
            endif

            ! If star is not cancelled, check and modify coefficients
            if (.not. star_cancel(i)) then

               c1 = coeff(i1)  ! Coeff of G1
               c2 = coeff(i2)  ! Coeff of G2
   
               ! Check that coefficient of 1st wave is real & positive
               if ( abs(aimag(coeff(i1))) > 1.0e-8 ) then
                  print *, &
                     'Error: 1st coeff in closed star is not real'
               else if ( real(coeff(i1)) < 0.0 ) then
                  print *, &
                     'Error: 1st coeff in closed star is negative'
               endif
   
               ! Check if G1 and G2 are complex conjugates
               ! If not, force them to become conjugates
               d = c1 - conjg(c2)
               if (abs(d) > 1.0e-8) then  
                  d = c2/abs(c2)
                  d = sqrt(d)
                  if (abs(aimag(d)) > 1.0E-5) then
                     if (aimag(d) < 0.0) then
                        d = -d
                     endif
                  endif
                  do j = star_begin(i), star_end(i)
                     coeff(j) = coeff(j)/d
                  enddo 
               endif
   
            endif

         endif

      else if (invert_flag == 1) then

         ! This is the second star in a pair related by inversion
         invert_flag    = -1
         star_invert(i) = -1
         wave_of_star(:,i) = wave(:,star_end(i))

      else

         ! If this star is not closed, and the previous star was not 
         ! the first star of a pair, this must be the first of a pair

         invert_flag    = 1
         star_invert(i) = 1
         wave_of_star(:,i) = wave(:,star_begin(i))

      endif

   enddo ! end loop over stars

   end subroutine make_stars
   !====================================================================


   !--------------------------------------------------------------------
   !****ip basis_mod/get_star
   ! SUBROUTINE
   !     get_star
   ! PURPOSE
   ! Takes an array "list" of list_N vectors, and generates a star by
   ! applying all elements of "group" to vector # root in input list.
   ! Note: Algorithm assumes that no duplicates exist in input list.
   !
   ! The list_index is an arbitrary list of integer labels that is
   ! divided and re-ordered the same way as the lists of wavevectors,
   ! so that the correspondence between a wave an integer label is
   ! maintained.
   !
   ! SYNPOSIS
   !    get_star(root, &
   !             list, list_index, list_N,  &
   !             star, star_index, star_N, star_phase)
   ! ARGUMENTS
   !
   !    Inputs:
   !    root       = index of vector to which symmetry elements of
   !                 the group are applied
   !
   !    Outputs:
   !    star       = the vectors of the resulting star
   !    star_index = indices of vectors of the resulting star
   !    star_N     = number of vectors in the resulting star
   !    list       = repacked list of set of vectors from the
   !                 input list that are not in the star
   !    list_index = indices of vectors of output list
   !    list_N     = # of vectors in output list
   !
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine get_star(root,  &
                       list, list_index, list_N,  &
                       star, star_index, star_N,star_phase)
   integer, intent(IN)       :: root
   integer, intent(INOUT)    :: list(3, max_list)
   integer, intent(INOUT)    :: list_index(max_list)
   integer, intent(INOUT)    :: list_N
   integer, intent(OUT)      :: star(3,max_star)
   integer, intent(OUT)      :: star_index(max_star)
   integer, intent(OUT)      :: star_N
   real(long), intent(OUT)   :: star_phase(max_star)
   !***

   ! Local variables
   real(long)  :: twopi
   integer     :: g1(3), g2(3), temp(3, max_list)
   integer     :: i, j, k, temp_index(max_list)
   logical     :: in_star(max_list)

   ! Debug variable
   integer :: check(max_star)
   ! End Debug

   ! Check precondition
   if ( list_N > max_list ) then
      write(6,*) 'Error: list_N > max_list in get_star'
      stop
   endif

   twopi   = 4.0_long*acos(0.0_long)
   star_N  = 0
   in_star = .false.
   star = 0

   ! g1 is the root wavevector of the new star
   g1 = list(:,root)

   do i=1, group%order

      ! Generate wavevector related by symmmetry element group%s(i)
      g2 = g1.dot.group%s(i)
      if (grid) g2 = G_to_bz(g2)

      ! Search for wavevector g2 in list. If not found, then j == 0.
      j = in_list(g2, list, list_N)

      if (j == 0) then ! Failure to find g2 in remaining list is fatal 

         write(6,*) 'Error: list does not contain entire star'
         write(6,*) 'g1 =',g1
         write(6,*) 'g2 =',g2
         write(6,*) in_list(g2,list,list_N), list_N
         write(6,'(3(i3))') list(:,:list_N)
         stop

      else ! Process symmetry related wavevector g2

         in_star(j) = .true.

         ! Search for g2 in the current star array
         k = in_list(g2, star, star_N)

         ! If g2 is not already in the star, add it
         if (k == 0) then

            ! Add to end of star and star_index arrays
            star_N = star_N + 1
            star(:, star_N) = g2
            star_index(star_N) = list_index(j)

            ! Compute phase and add to end of star_phase
            star_phase(star_N) = twopi*g1.dot.group%s(i)%v

            ! Debug Line
            check(star_N) = i
            ! End Debug

         endif

      endif
   enddo

   ! Make compressed list, list_index of waves that are not in this star
   temp = 0
   temp_index = 0
   j = 0
   do i=1, list_N
      if (.not.in_star(i)) then
         j = j + 1
         temp(:,j) = list(:,i)
         temp_index(j) = list_index(i)
      endif
   enddo
   list_N = j
   list = temp
   list_index = temp_index

   !!$ ! Debug Printing
   !!$	print *, 'get_star star'
   !!$	do i = 1,star_N
   !!$		print '(3i4,2f8.3,i4)',star(:,i),star_phase(i)/twopi, &
   !!$		   star(:,1) .dot. group%s(check(i))%v,check(i)
   !!$	enddo	

   end subroutine get_star
   !====================================================================


   !--------------------------------------------------------------------
   !****ip basis_mod/vector_cancel
   ! FUNCTION
   !    logical vector_cancel(G)
   ! ARGUMENTS
   !    integer G(3)
   ! RETURN
   !    True if cancellation occurs for vector G, False otherwise
   ! ALGORITHM
   !    Operate on G with all members of group.
   !    The function is true if any symmetry in the group both:
   !      1) Leaves G invariant (i.e., G.m = G), and
   !      2) Produces a nonzero phase, such that modulo(G.v, twopi) /= 0
   ! SOURCE
   !--------------------------------------------------------------------
   logical function vector_cancel(G)
   integer, intent(IN)             :: G(3)    ! reciprocal wavevector
   !***

   integer    :: i, Gp(3)
   real(long) :: phase, twopi

   vector_cancel = .false.
   do i=1, group%order
      Gp = G .dot. group%s(i)

      if (grid) Gp = G_to_bz(Gp)

      !write(6,FMT='(i5,3X,3I6,3X,3I6') i, G, Gp
      if (equal(Gp,G)) then
         phase = G .dot. group%s(i)%v
         phase = modulo(phase,1.0_long)
         if (equal(phase,1.0_long)) phase = phase - 1.0_long
         !write(6,FMT='("Phase=",F10.4)') phase
         if (.not.equal(phase,0.0_long)) then
            vector_cancel = .true.
            return
         endif
      endif
   enddo
   end function vector_cancel
   !====================================================================


   !--------------------------------------------------------------------
   !****ip basis_mod/in_list
   ! FUNCTION
   !    integer in_list(g,list,n)
   ! PURPOSE
   !    Searches for vector "g" in list of n vectors "list"
   ! RETURN
   !    in_list = i if g = list(i) for i <= n
   !    in_list = 0 if g is not in list(i) for i <= n
   !***
   !--------------------------------------------------------------------
   integer function in_list(g,list,n)
   integer   :: g(:),list(:,:),n
   integer   :: i

   ! Check array dimensions
   if (size(g) /= size(list,1)) then
      write(6,*) 'Incompatible array dimensions in in_list'
      stop
   endif

   in_list = 0
   if (n == 0) return
   do i=1, n
      if ( equal(g(:),list(:,i)) ) then
         in_list = i
         return
      endif
   enddo

   end function in_list
   !====================================================================


   !--------------------------------------------------------------------
   !****ip basis_mod/sort
   ! SUBROUTINE
   !     sort(n,a,index)
   ! ARGUMENTS
   !     integer    n         - # of elements
   !     real(long) a(:)      - array to be sorted
   !     integer    index(:)  - original indices of elements
   ! PURPOSE
   !    Sort an array a(n) by Shell variant of insertion.
   !    Algorithm: Based on Numerical Recipes in F77, pg. 323,
   !         Modified by D.M. to apply the permutation to an index
   !    On output,
   !       a is sorted in ascending order, a(1) <= a(2),...
   !       The same permuation has been applied to a(n) and index(n)
   !
   !       If on input, array index was ordered, with index(i) = i, 
   !       then on output index(i) contains original index of a(i),
   !       such that a_output(i) = a_input(index(i))
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine sort(n,a,index)
   integer    :: n, index(:)
   real(long) ::    a(:)
   !***

   integer    :: i, j, inc, k
   real(long) :: v

   ! Begin Shell method algorithm
   inc=1
 1 inc = 3*inc + 1
   if (inc.le.n) goto 1
 2 continue
      inc=inc/3
      do i=inc+1, n
         v=a(i)
         k=index(i)
         j=i
 3       if (a(j-inc).gt.v) then
            a(j)=a(j-inc)
            index(j)=index(j-inc)
            j=j-inc
            if (j.le.inc) goto 4
         goto 3
         endif
 4       a(j)=v
         index(j)=k
      enddo
   if (inc.gt.1) goto 2
   return
   end subroutine sort
   !====================================================================


   !--------------------------------------------------------------------
   !****ip basis_mod/G_sort
   !
   ! SUBROUTINE
   !    G_sort(n, a, index, reorder_real )
   !
   ! PURPOSE
   !
   !    Same as sort, except that a is an array of integer vectors
   !
   !    Sort an array a of integer vectors by Shell variant
   !    of insertion, using the module function G_lt to define
   !    a logical operator "less than" for two integer vectors
   !
   !    Source: Numerical recipes in F77, pg. 323,
   !            Modified by D.M. to produce index
   !            Modified by D.M. to order vectors
   !
   !    On output,
   !       a is sorted descending, a(:,1) >= a(:,2),...
   !       index(i) contains original index of a(:,i),
   !           i.e., a_output(i) = a_input(index(i))
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine G_sort(n, a, index, reorder_real)
   integer, intent(INOUT)    :: n, a(:,:), index(:)
   !***

   integer    :: i, j, inc, k, dim_a, v(size(a,1))
   real(long), intent(INOUT), optional :: reorder_real(:)
   real(long) :: reorder_index

   real(long), allocatable :: temp_real(:)
   integer, allocatable :: temp_int(:)
   integer, allocatable :: temp_wave(:,:)

   ! Check dimension of array a and array index
   dim_a = size(a,1)
   if ((dim_a < 0).or.(dim_a > 3)) then
       write(6,*) 'Invalid dim in G_sort'
       stop
   endif
   if (dim_a < dim) then
       write(6,*) 'Dimension of a < dim in G_sort'
       stop
   endif
   if (size(a,2).lt.n) then
      write(6,*) 'SIZE(a,2) < n in G_sort'
      stop
   endif
   if (size(index).lt.n) then
      write(6,*) 'SIZE(index) < n in G_sort'
      stop
   endif

   ! Begin Shell method algorithm
   inc=1
 1 inc = 3*inc + 1
   if (inc.le.n) goto 1
 2 continue
      inc=inc/3
      do i=inc+1, n
         v=a(:,i)
         k=index(i)
         if (PRESENT(reorder_real)) reorder_index = reorder_real(i)
         j=i
 3       if (G_lt(a(:,j-inc),v)) then
            a(:,j)=a(:,j-inc)
            index(j)=index(j-inc)
            if (PRESENT(reorder_real)) reorder_real(j)=reorder_real(j-inc)
            j=j-inc
            if (j.le.inc) goto 4
         goto 3
         endif
 4       a(:,j)=v
         index(j)=k
         if (PRESENT(reorder_real)) reorder_real(j) = reorder_index
      enddo
   if (inc.gt.1) goto 2

   return
   end subroutine G_sort
   !====================================================================


   !--------------------------------------------------------------------
   !****ip basis_mod/G_lt
   ! FUNCTION
   !    logical G_lt(g1,g2)
   ! PURPOSE
   !    Definition of "less than" for integer G vectors
   !    Treats elements of g1 and g2 as digits in a number
   ! RETURN
   !    true  if gl < g2
   !    false otherwise
   !***
   !--------------------------------------------------------------------
   logical function G_lt(g1,g2)
   integer, intent(IN) :: g1(3), g2(3)
   G_lt = .false.
   if (equal(g1,g2)) then
      G_lt = .false.
      return
   endif
   if (g1(1).lt.g2(1)) then
      G_lt = .true.
      return
   else if (g1(1).gt.g2(1)) then
      G_lt = .false.
      return
   else
      if (dim > 1) then
         if (g1(2).lt.g2(2)) then
            G_lt = .true.
            return
         else if (g1(2).gt.g2(2)) then
            G_lt = .false.
            return
         else
            if (dim == 3) then
               if (g1(3).lt.g2(3)) then
                  G_lt = .true.
                  return
               endif
            else
               write(6,*) 'Error in G_lt, dim=',dim
               stop
            endif
         endif
      endif
   endif
   end function G_lt
   !====================================================================


   !--------------------------------------------------------------------
   ! Utilities:
   !
   !  function f_basis:    returns value of basis function at point R
   !  routine assign_ivec: equates 3 and dim dimensional int vectors
   !  routine assign_rvec: equates 3 and dim dimensional real vectors
   !  function valid_wave: true if a wave in list
   !
   !--------------------------------------------------------------------


   !--------------------------------------------------------------------
   !****p basis_mod/f_basis
   ! FUNCTION
   !    real(long) f_basis(i_basis,G_basis,R)
   ! RETURN
   !    Value of basis function i_basis at Cartesian position R(3)
   !    a crystal with reciprocal basis vector array G_basis(3,3)
   ! SOURCE
   !--------------------------------------------------------------------
   real(long) function f_basis(  &
                  i_basis,       &! index for basis function
                  G_basis,       &! reciprocal basis vectors
                  R              &! Cartesian position
                  )
   ! Dummy arguments
   integer, intent(IN)     :: i_basis
   real(long), intent(IN)  :: G_basis(3,3), R(3)
   !***

   ! Parameter
   real(long) :: sqrt2

   ! Local variables
   integer      :: i_wave, i_star, iG(3)
   real(long)   :: G(3)
   complex(long):: phase

   sqrt2   = sqrt(2.0_long)
   f_basis = 0.0_long
   iG = 0
   iG(:dim) = wave_of_star(:,i_basis)
   i_wave = which_wave(iG(1),iG(2),iG(3))
   i_star  = star_of_wave(i_wave)

   do i_wave = star_begin(i_star), star_end(i_star)
      G = wave(:,i_wave).dot.G_basis(:,:)
      phase = coeff(i_wave)*exp( (0.0_long,1.0_long)*(R.dot.G))
      select case(star_invert(i_star))
      case(0)
         f_basis = f_basis + real(phase)
      case(+1)
         f_basis = f_basis + real(phase)/sqrt2
      case(-1)
         f_basis = f_basis - imag(phase)/sqrt2
      end select
   enddo
   end function f_basis
   !====================================================================


   !--------------------------------------------------------------------
   ! SUBROUTINE
   !    assign_ivec(v1,v2)
   ! PURPOSE
   !    Sets v1 = v2, for int vectors with variable dimensions
   !--------------------------------------------------------------------
   subroutine assign_ivec(v1,v2)
   integer, intent(IN)  :: v2(:)
   integer, intent(OUT) :: v1(:)
   integer :: dim1, dim2
   dim1 = size(v1)
   dim2 = size(v2)
   if (dim1 > dim2) then
      v1 = 0
      v1(1:dim2) = v2
   else if (dim1 < dim2) then
      v1 = v2(1:dim1)
   else
      v1 = v2
   endif
   end subroutine assign_ivec
   !====================================================================


   !--------------------------------------------------------------------
   ! SUBROUTINE
   !    assign_rvec(v1,v2)
   ! PURPOSE
   !    Sets v1 = v2, for real vectors with variable dimensions
   !--------------------------------------------------------------------
   subroutine assign_rvec(v1,v2)
   real(long), intent(IN)  :: v2(:)
   real(long), intent(OUT) :: v1(:)
   integer :: dim1, dim2
   dim1 = size(v1)
   dim2 = size(v2)
   if (dim1 > dim2) then
      v1         = 0.0_long
      v1(1:dim2) = v2
   else if (dim1 < dim2) then
      v1 = v2(1:dim1)
   else
      v1 = v2
   endif
   end subroutine assign_rvec
   !====================================================================


   logical function valid_wave(wave_k)
   !--------------------------------------------------------------------
   ! True if vector wave_k is listed in which_wave
   !--------------------------------------------------------------------
   integer, intent(IN) :: wave_k(3)
   integer :: k, stark
   valid_wave = .true.
   if ( any( abs(wave_k) > abs(G_max) ) ) then
      valid_wave = .false.
      return
   endif
   k = which_wave(wave_k(1),wave_k(2),wave_k(3))
   if (k == 0) then
      valid_wave = .false.
      return
   endif
   if ((k < 1).or.(k > N_wave)) then
      write(6,*)'Error in valid_wave: Invalid wave index =', k
      stop
   endif
   stark = star_of_wave(k)
   if ((stark < 1).or.(stark > N_star)) then
      write(6,*)'Error in valid_wave: Invalid star index'
      stop
   endif
   end function valid_wave
   !====================================================================

   !--------------------------------------------------------------------
   ! Output routines: output_waves, output_stars
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !****p basis_mod/output_waves
   ! SUBROUTINE
   !    output_waves(wave_port)
   ! ARGUMENTS
   !    wave_port = unit number for output file
   ! PURPOSE
   !    Output list of waves to file unit io_port
   !    Format: i, wave, star_of_wave, |coeff|, coeff/|coeff|, |G|
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine output_waves(wave_port, group_name)

   integer      :: wave_port
   character(*) :: group_name
   !***

   integer       :: i
   real(long)    :: G_vec(3), Gsq, G_abs, coeff_abs
   complex(long) :: coeff_phase
   character(80) :: fmt
   type(version_type) :: version


   version%major = 1
   version%minor = 0
   call output_version(version,wave_port)
   call output_unit_cell(wave_port,'F')
   call output(group_name,'group_name',o=wave_port,f='A')
   call output(N_wave,'N_wave',o=wave_port,f='A')
   do i=1, N_wave
       G_vec = wave(:,i) .dot. G_basis
       Gsq   = G_vec .dot. G_vec
       G_abs = dsqrt(Gsq)
       coeff_abs   = dble(coeff(i))**2 + aimag(coeff(i))**2
       coeff_abs   = dsqrt(coeff_abs)
       coeff_phase = coeff(i)/coeff_abs
       fmt = '(I8,2X,'//trim(int_string(dim))//'I5,I6,3F12.6,F15.6)'
       write(wave_port,FMT=trim(fmt)) i, wave(:,i), star_of_wave(i), &
                                      coeff_abs, coeff_phase, G_abs
   end do

   end subroutine output_waves
   !====================================================================


   !--------------------------------------------------------------------
   !****p basis_mod/output_stars
   ! SUBROUTINE
   !    output_stars(o_port,keep_cancel)
   ! ARGUMENTS
   !    o_port = io-port of the wave file
   ! PURPOSE
   !    Output list of stars to supplied file o_port
   !    Format: i, wave_of_star, star_count, star_invert
   !    If keep_cancel = .true., then also print cancel
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine output_stars(o_unit,keep_cancel)
   use io_mod,     only : output
   use string_mod, only : int_string

   integer, optional   :: o_unit
   logical, optional   :: keep_cancel
   !***

   integer             :: i, G(dim)
   character(50)       :: fmt
   call output(N_star,'N_star',o=o_unit,f='R')
   if (present(keep_cancel).and.keep_cancel) then
       write(o_unit,*) &
      '    i     wave_of_star      count    invert    cancel'
   else
      write(o_unit,*) &
      '    i     wave_of_star      count    invert'
   endif
   do i=1, N_star
      G = wave_of_star(:,i)
      if (present(keep_cancel).and.keep_cancel) then
         fmt = 'I6,2X,'//trim(int_string(dim))//'I5,2X,2I5,L5'
         write(o_unit,FMT='(I6,2X,3I5,4I10,L5)') &
               i, G, star_count(i), star_invert(i), star_cancel(i)
      else
         fmt = 'I6,2X,'//trim(int_string(dim))//'I5,2X,2I5'
         write(o_unit,fmt) &
               i, G, star_count(i), star_invert(i)
      endif
   enddo
   end subroutine output_stars
   !====================================================================


   !--------------------------------------------------------------------
   !****p basis_mod/reorder_basis
   ! SUBROUTINE
   !    reorder_basis(new_basis,old_wave,old_basis,swap,new_stars)
   ! PURPOSE
   !    Reorders basis functions in new_wave to match ordering in
   !    old_wave.
   ! COMMENT
   !    Assumes each wave in new_wave and old_wave corresponds to
   !    a single basis function, and that there are no repeats in
   !    old_wave or new_wave
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine reorder_basis(new_basis,old_wave,old_basis,swap,new_stars)
   integer, intent(IN)            :: old_wave(:,:)  ! (dim,N_old)
   real(long), intent(IN)         :: old_basis(:,:) ! (dim,N_new)
   real(long), intent(OUT)        :: new_basis(:,:) ! (dim,N_max)
   integer, intent(OUT), optional :: swap(:)        ! (N_max)
   logical, intent(OUT), optional :: new_stars(:)   ! (N_max)
   ! N_max = max(N_old,N_new)
   !***

   integer :: i, j, k
   integer :: N_new, N_old
   integer, allocatable, dimension(:) :: l_swap
   logical, allocatable, dimension(:) :: new
   ! local copy of swap

   integer, dimension(3) :: G
   ! true if new_basis(i) is filled

   if ( dim == 1 ) then
      print *, 'Reordering not available for 1-d systems'
      return
   elseif (dim == 1) then
      print *, 'Reordering not available for 2-d systems'
   endif

   N_new = size(new_basis,2)
   N_old = min(size(old_wave,2),size(old_basis,2))
   print *, N_old, N_new
   allocate( l_swap(max(N_old,N_new)), stat = j )
   if (j.ne.0) stop 'Error allocating l_swap in reorder_basis, basis_mod'
   allocate(new(N_new), stat = j)
   if (j.ne.0) stop 'Error allocating new in reorder_basis, basis_mod'

   ! Check for size arrays
   if ( size(new_basis,1) .ne. size(old_basis,1) ) then
      write(6,*) 'Error.  new_basis and old_basis of different sizes'
      write(6,*) 'shape(new_basis) = ', shape(new_basis)
      write(6,*) 'shape(old_basis) = ', shape(old_basis)
      stop
   endif

   print *, 'reorder N_old = ', N_old
   print *, 'reorder N_new = ', N_new

   l_swap = 0
   new_basis = 1.0e-6
   !!$  print *, G_max
   !!$  do i = 1,N_old
   !!$     G = old_wave(:,i)
   !!$     ! Find which wave in current waves corresponds to the old_wave
   !!$     if ( any(abs(G) > G_max)) then
   !!$        print *, G, i
   !!$        cycle
   !!$     endif
   !!$     k = star_of_wave( which_wave(G(1), g(2), g(3) ))
   !!$
   !!$     if ( (k<= N_new) .and. (k > 0) ) then
   !!$        new_basis(:,k) = old_basis(:,i)
   !!$        filled(k) = .true.
   !!$
   !!$        l_swap(k) = i
   !!$     endif
   !!$  enddo


   ! print *, 'old loop', dim
   do i = 1,N_old
      G(:) = 0
      G(:dim) = old_wave(:dim,i)
      if ( any(abs(G(:dim)) > G_max(:dim))) then
         cycle
      else
         k = star_of_wave( which_wave(G(1),G(2),G(3)))
         if ( (k<=N_new) .and. (k>0) ) then
            new_basis(:,k) = old_basis(:,i)
            l_swap(k) = i
            new(k) = .false.
         endif
      endif
   enddo


   !!$  ! Check that all of new_basis is filled
   !!$  if ( .not. all(filled)) then
   !!$     ! Check which is the first not filled
   !!$     do i = 1,min(N_new,N_old)
   !!$        if ( .not. filled(i)) then
   !!$           print '(a40,i4,2x,3i4,2x,3i4)', &
   !!$                'Error.  New_basis not filled at i = ', i, old_wave(:,i)
   !!$           G = old_wave(:,i)
   !!$           k = which_wave(G(1),g(2),g(3))
   !!$           print *, G, k, star_of_wave(k), wave_of_star(:,i)
   !!$           !              stop
   !!$        endif
   !!$     enddo
   !!$  endif

   if (PRESENT(swap)) then
      do i = 1,min(min(size(swap),N_old),N_new)
         swap(i) = l_swap(i)
      enddo
   endif
   if (PRESENT(new_stars)) then
      do i = 1,min( N_new, size(new_stars))
         new_stars(i) = new(i)
      enddo
   endif
   print *, 'return from reorder_basis'
   deallocate(l_swap)
   deallocate(new)
   end subroutine reorder_basis
   !====================================================================


   !--------------------------------------------------------------------
   !****p basis_mod/make_dGsq
   ! SUBROUTINE
   !    make_dGsq(dGsq,dGG)
   ! PURPOSE
   !    Calculate dGsq given dGG_basis
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine make_dGsq( dGsq, dGG)
   real(long), intent(OUT) :: dGsq(:)
   real(long), intent(IN) :: dGG(:,:)
   !***

   !Local Variabls
   integer :: star, j, k
   real(long), dimension(dim) :: Gvec

   dGsq = 0.0_long
   ! loop over all stars
   do star = 1,N_star
      Gvec = dble( wave_of_star(:dim,star) )
      do j = 1,dim
         do k = 1,dim
            dGsq(star) = dGsq(star) + &
                         Gvec(j) * Gvec(k) * dGG(j,k)
         enddo
      enddo
   enddo

   end subroutine make_dGsq
   !====================================================================


   !--------------------------------------------------------------------
   !****p basis_mod/scattering_intensity
   ! SUBROUTINE
   !   scattering_intensity( G_basis, field, contrast, q, scat )
   ! PURPOSE
   !   Calculate the scattering intensities give a set of fields
   !   (e.g. density fields of SCFT solution ) and some  scattering
   !   densities (e.g. electron densities)
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine scattering_intensity( G_basis, field, contrast, q, scat )
   implicit none
   real(long), intent(IN), dimension(:,:) :: G_basis
   real(long), intent(IN), dimension(:,:) :: field
   real(long), intent(IN), dimension(:) :: contrast
   real(long), intent(OUT), dimension(:) :: q
   real(long), intent(OUT), dimension(:) :: scat
   !***

   ! Local variables
   integer :: i, j, k
   integer :: N
   integer :: alpha, N_blocks, beta
   real(long), dimension(dim) :: vector

   N = min(N_star,size(q))
   N_blocks = size(field,1)

   ! Calculate q for each star
   do i = 1,N
      vector = wave_of_star(:,i) .dot. G_basis
      q(i) = sqrt( vector .dot. vector )
   enddo

   print *, 'contrast'
   print *, contrast

   ! Calculate the scattering intensity for each star
   scat = 0.0_long
   do i = 1,N
      do alpha = 1, N_blocks
         do beta = alpha+1,N_blocks
            scat(i) = ( field(i,alpha) * contrast(alpha) - &
                        field(i,beta) * contrast(beta) ) ** 2 &
                      + scat(i)
         enddo
      enddo
      ! correct for normalization conventions
      ! Basis functions are normalized such that <f_i | f_i > = 1
      ! Scattering is normalized so that <f_i | f_i > = starcount(i)

      scat(i) = scat(i) * star_count(i)
   enddo

   end subroutine scattering_intensity
   !====================================================================


   !--------------------------------------------------------------------
   ! Array access functions:
   !
   ! These functions access elements of arrays, and have names given
   ! by the corresponding array names prefixed by "fn_". The functions
   ! may used as replacements for the arrays themselves, to provide
   ! array access with automatic check of bounds of both the array
   ! indices and (when possible) the desired array element.
   !
   ! Syntactical subtleties:
   ! i)  fn_wave_of_star returns array of dimension 3, rather than dim
   ! ii) fn_which_wave requires an integer array with dimension 3,
   !     rather than 3 integer scalars, as its input.
   !
   !--------------------------------------------------------------------


   integer function fn_star_invert(i)
   !--------------------------------------------------------------------
   !  Returns array element star_invert(i), after checking bounds
   !--------------------------------------------------------------------
   integer, intent(IN) :: i
   if  ( (i >= 1) .and. (i <= N_star) ) then
      fn_star_invert = star_invert(i)
   else
      write(6,*)'Error in fn_star_invert: Invalid star index =', i
      stop
   endif
   if ( abs(fn_star_invert) > 1 ) then
      write(6,*)'Error: invalid fn_star_invert =', fn_star_invert
      stop
   endif
   end function fn_star_invert
   !====================================================================


   integer function fn_star_of_wave(i)
   !--------------------------------------------------------------------
   !  Returns array element star_of_wave(i), after checking bounds
   !--------------------------------------------------------------------
   integer, intent(IN) :: i
   if  ( (i >= 1) .and. (i <= N_wave) ) then
      fn_star_of_wave = star_of_wave(i)
   else
      write(6,*)'Error in fn_star_of_wave: Invalid wave index =', i
      stop
   endif
   if ( ( fn_star_of_wave < 1 ).or.( fn_star_of_wave > N_star )) then
      write(6,*)'Error: invalid fn_star_of_wave =', fn_star_of_wave
      stop
   endif
   end function fn_star_of_wave
   !====================================================================


   function fn_wave_of_star(i)
   !--------------------------------------------------------------------
   !  Returns array element wave_of_star(i), after checking bounds
   !
   !  Note: function has fixed dimension 3, whereas vectors in array
   !    wave have dimension dim, because array is allocatable, but
   !    F90 requires array functions to have a constant dimension.
   !--------------------------------------------------------------------
   integer              :: fn_wave_of_star(3)
   integer, intent(IN)  :: i
   fn_wave_of_star = 0.0_long
   if  ( (i >= 1) .and. (i <= N_star) ) then
      fn_wave_of_star(1:dim) = wave_of_star(:,i)
   else
      write(6,*)'Error in fn_wave_of_star: Invalid star index =', i
      stop
   endif
   if (.not.valid_wave(fn_wave_of_star)) then
      write(6,*) 'Error: Invalid fn_wave_of_star=', fn_wave_of_star
   endif
   end function fn_wave_of_star
   !====================================================================


   integer function fn_which_wave(vec)
   !--------------------------------------------------------------------
   !  Returns array element which_wave(vec(1),vec(2),vec(3)),
   !  after checking for validity of vec and bounds on result
   !--------------------------------------------------------------------
   integer, intent(IN) :: vec(3)
   if (valid_wave(vec)) then
      fn_which_wave = which_wave(vec(1),vec(2),vec(3))
   else
      write(6,*) 'Error in fn_which_wave: Invalid wave=',vec
      stop
   endif
   if ( (fn_which_wave < 1).or.(fn_which_wave > N_wave)) then
      write(6,*) 'Error: fn_which_wave =', fn_which_wave
   endif
   end function fn_which_wave
   !====================================================================


end module basis_mod
