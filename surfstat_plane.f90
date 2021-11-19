subroutine surfstat_plane
   ! This subroutine calculates surface states using
   ! iterative Green's function method  ---  J.Phys.F.Met.Phys.15(1985)851-858
   !
   ! History:
   !
   !         by Quan Sheng Wu on 4/20/2010
   !
   !            mpi version      4/21/2010

   ! revised in 2018 by Yuefeng Yin
   ! based on the original code of Quansheng Wu
   ! Extend the SS_NEGF calculation to 2D plane
   ! current version is not compatible with PHONON calculation. Use at your own risk.


   use wmpi
   use para
   implicit none



   integer :: ierr

   integer :: arclfile, arcrfile, arcbulkfile

   ! general loop index
   integer :: i,j,io

   ! kpoint loop index
   integer :: ikp

   ! final counter
   integer :: ikp_f


   real(dp) :: k(2),w


   real(dp) :: k1min_shape, k1max_shape, k2min_shape, k2max_shape

   real(dp) :: time_start, time_end

   real(dp), allocatable :: omega(:)
   real(dp), allocatable :: k12(:,:)
   real(dp), allocatable :: k12_shape(:,:)

   real(dp), allocatable :: dos_l(:,:)
   real(dp), allocatable :: dos_l_mpi(:,:)
   real(dp), allocatable :: dos_r(:,:)
   real(dp), allocatable :: dos_r_mpi(:,:)
   real(dp), allocatable :: dos_l_only(:,:)
   real(dp), allocatable :: dos_r_only(:,:)
   real(dp), allocatable :: dos_bulk(:,:)
   real(dp), allocatable :: dos_bulk_mpi(:,:)

   complex(dp), allocatable :: ones(:,:)
   complex(dp), allocatable :: GLL(:,:), GRR(:,:), GB(:,:)
   complex(dp), allocatable :: H00(:,:), H01(:,:)

   allocate( omega(omeganum)) !
   allocate( k12(2, nk1*nk2))
   allocate( k12_shape(2, nk1*nk2))
   allocate( dos_l(nk1*nk2,omeganum))
   allocate( dos_l_mpi(nk1*nk2,omeganum))
   allocate( dos_r(nk1*nk2,omeganum))
   allocate( dos_l_only(nk1*nk2,omeganum))
   allocate( dos_r_only(nk1*nk2,omeganum))
   allocate( dos_r_mpi(nk1*nk2,omeganum))
   allocate( dos_bulk(nk1*nk2,omeganum))
   allocate( dos_bulk_mpi(nk1*nk2,omeganum))
   allocate( GLL(ndim,ndim), GRR(ndim,ndim), GB(ndim, ndim))
   omega=0d0
   k12=0d0
   k12_shape=0d0
   dos_l=0d0
   dos_l_mpi=0d0
   dos_r=0d0
   dos_r_mpi=0d0
   dos_bulk=0d0
   dos_bulk_mpi=0d0

   ikp=0
   do i= 1, nk1
      do j= 1, nk2
         ikp=ikp+1
         k12(:, ikp)=K2D_start+ (i-1)*K2D_vec1/dble(nk1-1) &
                    + (j-1)*K2D_vec2/dble(nk2-1)
         k12_shape(:, ikp)= k12(1, ikp)* Ka2+ k12(2, ikp)* Kb2
      enddo
   enddo

   k1min_shape= minval(k12_shape(1, :))
   k2min_shape= minval(k12_shape(2, :))
   k1max_shape= maxval(k12_shape(1, :))
   k2max_shape= maxval(k12_shape(2, :))

   allocate(H00(Ndim, Ndim))
   allocate(H01(Ndim, Ndim))
   allocate(ones(Ndim, Ndim))
   GLL= 0d0
   GRR= 0d0
   H00= 0d0
   H01= 0d0
   ones= 0d0

   do i=1,Ndim
      ones(i,i)=1.0d0
   enddo


   eta=(omegamax- omegamin)/dble(omeganum)*1.5d0 !
   do i= 1, omeganum
      omega(i)=omegamin+(i-1)*(omegamax-omegamin)/dble(omeganum)
   enddo !

   !> deal with phonon system
   !> for phonon system, omega should be changed to omega^2
   if (index(Particle,'phonon')/=0) then
      omega= omega*omega
   endif

   time_start= 0d0
   time_end= 0d0
   do ikp= 1+cpuid, nk1*nk2, num_cpu
      if (cpuid==0.and. mod(ikp/num_cpu, 100)==0) &
         write(stdout, *) 'Arc, ik ', ikp, 'Nk',Nk1*Nk2, 'time left', &
         (nk1*nk2-ikp)/num_cpu*(time_end- time_start), ' s'
      call now(time_start)
      k(1)= k12(1, ikp)
      k(2)= k12(2, ikp)


      if (index(Particle,'phonon')/=0.and.LOTO_correction) then
         call ham_qlayer2qlayer_LOTO(k,H00,H01)
      else
         call ham_qlayer2qlayer(k,H00,H01)
      endif


      !> calculate surface green function
      ! there are two method to calculate surface green's function
      ! the method in 1985 is better, you can find the ref in the
      ! subroutine
      do  j=1,omeganum
        w=omega(j)
        call surfgreen_1985(w,GLL,GRR,GB,H00,H01,ones)
      ! call surfgreen_1984(w,GLL,GRR,H00,H01,ones)

      ! calculate spectral function
      do i= 1, NtopOrbitals
         io= TopOrbitals(i)
         dos_l(ikp,j)=dos_l(ikp,j)- aimag(GLL(io,io))
      enddo ! i
      do i= 1, NBottomOrbitals
         io= Ndim- Num_wann+ BottomOrbitals(i)
         dos_r(ikp,j)=dos_r(ikp,j)- aimag(GRR(io,io))
      enddo ! i
      do i= 1, Ndim
         dos_bulk(ikp,j)=dos_bulk(ikp,j)- aimag(GB (i,i))
      enddo ! i
    enddo ! new j
      call now(time_end)

   enddo ! ikp

   !> we don't have to do allreduce operation
#if defined (MPI)
   call mpi_reduce(dos_l, dos_l_mpi, size(dos_l),mpi_double_precision,&
                   mpi_sum, 0, mpi_comm_world, ierr)
   call mpi_reduce(dos_r, dos_r_mpi, size(dos_r),mpi_double_precision,&
                   mpi_sum, 0, mpi_comm_world, ierr)
   call mpi_reduce(dos_bulk, dos_bulk_mpi, size(dos_bulk), mpi_double_precision,&
                   mpi_sum, 0, mpi_comm_world, ierr)
#else
   dos_l_mpi= dos_l
   dos_r_mpi= dos_r
   dos_bulk_mpi= dos_bulk
#endif

   do ikp=1, Nk1*Nk2
     do j=1,omeganum
      dos_l_only(ikp,j)= dos_l_mpi(ikp,j)- dos_bulk_mpi(ikp,j)
      if (dos_l_only(ikp,j)<0) dos_l_only(ikp,j)=eps9
      dos_r_only(ikp,j)= dos_r_mpi(ikp,j)- dos_bulk_mpi(ikp,j)
      if (dos_r_only(ikp,j)<0) dos_r_only(ikp,j)=eps9
    enddo
   enddo

   outfileindex= outfileindex+ 1
   arclfile= outfileindex
   outfileindex= outfileindex+ 1
   arcrfile= outfileindex
   outfileindex= outfileindex+ 1
   arcbulkfile= outfileindex
   if (cpuid.eq.0)then
      open (unit=arclfile, file='arc.dat_l')
      open (unit=arcrfile, file='arc.dat_r')
      open (unit=arcbulkfile, file='arc.dat_bulk')
      ikp_f=0
      do ikp=1, nk1*nk2
        do j=1, omeganum
         ikp_f=ikp_f+1
         if (log(dos_l_only(ikp,j))>arc_threshold) write(arclfile, '(30f16.8)')k12_shape(:, ikp), omega(j),log(dos_l_mpi(ikp,j)), log(dos_l_only(ikp,j))
         if (mod(ikp_f, nk2*omeganum)==0) write(arclfile, *)' '
         if (log(dos_r_only(ikp,j))>arc_threshold) write(arcrfile, '(30f16.8)')k12_shape(:, ikp), omega(j),log(dos_r_mpi(ikp,j)), log(dos_r_only(ikp,j))
         if (mod(ikp_f, nk2*omeganum)==0) write(arcrfile, *)' '
         write(arcbulkfile, '(30f16.8)')k12_shape(:, ikp), omega(j), log(abs(dos_bulk_mpi(ikp,j)))
         if (mod(ikp_f, nk2*omeganum)==0) write(arcbulkfile, *)' '
       enddo
      enddo
      close(arclfile)
      close(arcrfile)
      close(arcbulkfile)

      write(stdout,*)'ndim',ndim
      write(stdout,*)'Nk1,Nk2,eta', Nk1, Nk2, eta
      write(stdout,*)'calculate density of state successfully'
   endif

!new plotting function surfstat_plane_l/r.gnu
!will also supply l only r only and bulk version
! for gnuplot
! based on ek_bulk2D.f90 and fermiarc
  outfileindex=outfileindex+1
  if (cpuid==0)then
   open(unit=outfileindex, file='surfstat_plane_l_only.gnu')
   write(outfileindex, '(a)')"set encoding iso_8859_1"
   write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
   write(outfileindex, '(a)')"#set output 'surfstat_plane_l_only.eps'"
   write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
      ' size 1920, 1680 font ",36"'
   write(outfileindex, '(a)')"set output 'surfstat_plane_l_only.png'"
   write(outfileindex, '(a)')'set palette defined (0  "white", 6 "red", 20 "black" )'
   write(outfileindex, '(a)')'unset key'
   write(outfileindex, '(a)')'set pm3d'
   write(outfileindex, '(a)')'set origin 0.2, 0'
   write(outfileindex, '(a)')'set size 0.8, 1'
   write(outfileindex, '(a)')'set border lw 3'
   write(outfileindex, '(a)')'#set xtics font ",24"'
   write(outfileindex, '(a)')'#set ytics font ",24"'
   write(outfileindex, '(a)')'set size ratio -1'
   write(outfileindex, '(a)')'set xtics'
   write(outfileindex, '(a)')'set ytics'
   write(outfileindex, '(a)')'set view 80,60'
   write(outfileindex, '(a)')'set xlabel "k_1"'
   write(outfileindex, '(a)')'set ylabel "k_2"'
   write(outfileindex, '(a)')'set zlabel "Energy (eV)" rotate by 90'
   write(outfileindex, '(a)')'unset colorbox'
   write(outfileindex, '(a)')'set autoscale fix'
   write(outfileindex, '(a)')'set pm3d interpolate 4,4'
   write(outfileindex, '(2a)')"splot 'arc.dat_l' u 1:2:3:(exp($5)) w pm3d"
!   write(outfileindex, '(2a)')"      'arc.dat_l' u 4:5:9 w pm3d"

   close(outfileindex)

endif ! cpuid

  outfileindex=outfileindex+1

  if (cpuid==0)then
   open(unit=outfileindex, file='surfstat_plane_r_only.gnu')
   write(outfileindex, '(a)')"set encoding iso_8859_1"
   write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
   write(outfileindex, '(a)')"#set output 'surfstat_plane_r_only.eps'"
   write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
      ' size 1920, 1680 font ",36"'
   write(outfileindex, '(a)')"set output 'surfstat_plane_r_only.png'"
   write(outfileindex, '(a)')'set palette defined (0  "white", 6 "red", 20 "black" )'
   write(outfileindex, '(a)')'unset key'
   write(outfileindex, '(a)')'set pm3d'
   write(outfileindex, '(a)')'set origin 0.2, 0'
   write(outfileindex, '(a)')'set size 0.8, 1'
   write(outfileindex, '(a)')'set border lw 3'
   write(outfileindex, '(a)')'#set xtics font ",24"'
   write(outfileindex, '(a)')'#set ytics font ",24"'
   write(outfileindex, '(a)')'set size ratio -1'
   write(outfileindex, '(a)')'set xtics'
   write(outfileindex, '(a)')'set ytics'
   write(outfileindex, '(a)')'set view 80,60'
   write(outfileindex, '(a)')'set xlabel "k_1"'
   write(outfileindex, '(a)')'set ylabel "k_2"'
   write(outfileindex, '(a)')'set zlabel "Energy (eV)" rotate by 90'
   write(outfileindex, '(a)')'unset colorbox'
   write(outfileindex, '(a)')'set autoscale fix'
   write(outfileindex, '(a)')'set pm3d interpolate 4,4'
   write(outfileindex, '(2a)')"splot 'arc.dat_r' u 1:2:3:(exp($5)) w pm3d"
!   write(outfileindex, '(2a)')"      'arc.dat_l' u 4:5:9 w pm3d"

   close(outfileindex)

endif ! cpuid

outfileindex=outfileindex+1

if (cpuid==0)then
 open(unit=outfileindex, file='surfstat_plane_r.gnu')
 write(outfileindex, '(a)')"set encoding iso_8859_1"
 write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
 write(outfileindex, '(a)')"#set output 'surfstat_plane_r.eps'"
 write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
    ' size 1920, 1680 font ",36"'
 write(outfileindex, '(a)')"set output 'surfstat_plane_r.png'"
 write(outfileindex, '(a)')'set palette defined ( -10 "#194eff", 0 "white", 10 "red" )'
 write(outfileindex, '(a)')'unset key'
 write(outfileindex, '(a)')'set pm3d'
 write(outfileindex, '(a)')'set origin 0.2, 0'
 write(outfileindex, '(a)')'set size 0.8, 1'
 write(outfileindex, '(a)')'set border lw 3'
 write(outfileindex, '(a)')'#set xtics font ",24"'
 write(outfileindex, '(a)')'#set ytics font ",24"'
 write(outfileindex, '(a)')'set size ratio -1'
 write(outfileindex, '(a)')'set xtics'
 write(outfileindex, '(a)')'set ytics'
 write(outfileindex, '(a)')'set view 80,60'
 write(outfileindex, '(a)')'set xlabel "k_1"'
 write(outfileindex, '(a)')'set ylabel "k_2"'
 write(outfileindex, '(a)')'set zlabel "Energy (eV)" rotate by 90'
 write(outfileindex, '(a)')'unset colorbox'
 write(outfileindex, '(a)')'set autoscale fix'
 write(outfileindex, '(a)')'set pm3d interpolate 4,4'
 write(outfileindex, '(2a)')"splot 'arc.dat_r' u 1:2:3:4 w pm3d"
!   write(outfileindex, '(2a)')"      'arc.dat_l' u 4:5:9 w pm3d"

 close(outfileindex)

endif ! cpuid

outfileindex=outfileindex+1

if (cpuid==0)then
 open(unit=outfileindex, file='surfstat_plane_l.gnu')
 write(outfileindex, '(a)')"set encoding iso_8859_1"
 write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
 write(outfileindex, '(a)')"#set output 'surfstat_plane_l.eps'"
 write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
    ' size 1920, 1680 font ",36"'
 write(outfileindex, '(a)')"set output 'surfstat_plane_l.png'"
 write(outfileindex, '(a)')'set palette defined ( -10 "#194eff", 0 "white", 10 "red" )'
 write(outfileindex, '(a)')'unset key'
 write(outfileindex, '(a)')'set pm3d'
 write(outfileindex, '(a)')'set origin 0.2, 0'
 write(outfileindex, '(a)')'set size 0.8, 1'
 write(outfileindex, '(a)')'set border lw 3'
 write(outfileindex, '(a)')'#set xtics font ",24"'
 write(outfileindex, '(a)')'#set ytics font ",24"'
 write(outfileindex, '(a)')'set size ratio -1'
 write(outfileindex, '(a)')'set xtics'
 write(outfileindex, '(a)')'set ytics'
 write(outfileindex, '(a)')'set view 80,60'
 write(outfileindex, '(a)')'set xlabel "k_1"'
 write(outfileindex, '(a)')'set ylabel "k_2"'
 write(outfileindex, '(a)')'set zlabel "Energy (eV)" rotate by 90'
 write(outfileindex, '(a)')'unset colorbox'
 write(outfileindex, '(a)')'set autoscale fix'
 write(outfileindex, '(a)')'set pm3d interpolate 4,4'
 write(outfileindex, '(2a)')"splot 'arc.dat_l' u 1:2:3:4 w pm3d"
!   write(outfileindex, '(2a)')"      'arc.dat_l' u 4:5:9 w pm3d"

 close(outfileindex)

endif ! cpuid


#if defined (MPI)
   call mpi_barrier(mpi_cmw, ierr)
#endif
   deallocate( omega)
   deallocate( k12)
   deallocate( k12_shape)
   deallocate( dos_l)
   deallocate( dos_l_mpi)
   deallocate( dos_r)
   deallocate( dos_l_only)
   deallocate( dos_r_only)
   deallocate( dos_r_mpi)
   deallocate( dos_bulk)
   deallocate( dos_bulk_mpi)
   deallocate( GLL, GRR, GB)
   deallocate(H00)
   deallocate(H01)
   deallocate(ones)


   return
end subroutine surfstat_plane