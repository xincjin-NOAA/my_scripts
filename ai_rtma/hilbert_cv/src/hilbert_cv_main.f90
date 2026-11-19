!=============================================================================
program hilbert_cv
!=============================================================================
! Standalone version of the GSI 2DVar (RTMA) Hilbert-curve cross-validation
! grouping. It reproduces the group assignment done in GSI by
! apply_hilbertcurve (src/gsi/support_2dvar.f90) for one ob type:
!
!   1. duplicates (same grid x,y within 1e-3) keep only the ob nearest the
!      analysis time;
!   2. grid coordinates are scaled to (0,1] by nlon, nlat;
!   3. if there are more than ngrpsfact0*ngroups obs, subroutine hilbert
!      (Purser's code, unchanged) assigns each ob to a validation group
!      1..ngroups, or 0 (not in any group, i.e. always assimilated).
!
! In GSI, the obs in ONE chosen group (ncgroup from convinfo, or the date-based
! choice when random_cvgrp) are withheld from the analysis (usage set to
! ncmiter or usagecv; their duplicates to usage_dup). Here the program only
! writes the group of every ob, so the caller can pick any group.
!
! Usage:   hilbert_cv <infile> <outfile> <ngroups> <nlon> <nlat> [ngrpsfact0]
! infile:  one ob per line, free format:
!            station_id kx lon_deg lat_deg x_grid y_grid time_offset_h
!          x_grid in [1,nlon], y_grid in [1,nlat] (GSI dlon, dlat).
!          GSI only cross-validates surface kx 180-188, 280-288, 192-195,
!          292-295 with usage < 6; filter the input the same way to match.
! outfile: the input columns plus "group": 1..ngroups, 0 = no group,
!          -1 = duplicate removed before grouping.
!=============================================================================
use kinds, only: r_kind, i_kind
implicit none

real(r_kind), parameter :: epsilon = 1.e-03_r_kind
character(len=512)  :: infile, outfile, arg
character(len=8), allocatable :: cstn(:)
integer(i_kind), allocatable  :: kx(:), flag(:), grp(:), ipoint(:), test_set(:)
real(r_kind), allocatable     :: alon(:), alat(:), dlon(:), dlat(:), tim(:), hx(:), hy(:)
integer(i_kind) :: ngroups, ngrpsfact0, nob, ncross, i, j, ios, n, nlon, nlat
character(len=256) :: line

if (command_argument_count() < 5) then
   print *, 'usage: hilbert_cv infile outfile ngroups nlon nlat [ngrpsfact0=10]'
   stop 1
endif
call get_command_argument(1, infile)
call get_command_argument(2, outfile)
call get_command_argument(3, arg); read(arg, *) ngroups
call get_command_argument(4, arg); read(arg, *) nlon
call get_command_argument(5, arg); read(arg, *) nlat
ngrpsfact0 = 10
if (command_argument_count() >= 6) then
   call get_command_argument(6, arg); read(arg, *) ngrpsfact0
endif

! count lines
open(10, file=trim(infile), status='old', action='read')
nob = 0
do
   read(10, '(a)', iostat=ios) line
   if (ios /= 0) exit
   if (len_trim(line) > 0) nob = nob + 1
enddo
rewind(10)
allocate(cstn(nob), kx(nob), alon(nob), alat(nob), dlon(nob), dlat(nob), tim(nob), &
         flag(nob), grp(nob), ipoint(nob), hx(nob), hy(nob))
do i = 1, nob
   read(10, *) cstn(i), kx(i), alon(i), alat(i), dlon(i), dlat(i), tim(i)
enddo
close(10)

! 1. duplicates: keep the ob nearest the analysis time (as in GSI)
flag = 1
do j = 1, nob
   do i = j+1, nob
      if (abs(dlon(i)-dlon(j)) < epsilon .and. abs(dlat(i)-dlat(j)) < epsilon) then
         if (abs(tim(i)) >= abs(tim(j))) then
            flag(i) = -1
         else
            flag(j) = -1
         endif
      endif
   enddo
enddo
ncross = 0
do i = 1, nob
   if (flag(i) > 0) then
      ncross = ncross + 1
      hx(ncross) = dlon(i) / nlon     ! 2. scale to (0,1]
      hy(ncross) = dlat(i) / nlat
      ipoint(ncross) = i
   endif
enddo

! 3. Hilbert-curve groups
grp = -1
where (flag > 0) grp = 0
if (ncross > ngrpsfact0*ngroups) then
   allocate(test_set(ncross))
   call hilbert(ngroups, ncross, hx(1:ncross), hy(1:ncross), test_set)
   do n = 1, ncross
      grp(ipoint(n)) = test_set(n)
   enddo
else
   print *, 'not enough obs for cross validation: ncross=', ncross, &
            ' <= ngrpsfact0*ngroups=', ngrpsfact0*ngroups
endif

open(20, file=trim(outfile), status='replace', action='write')
do i = 1, nob
   write(20, '(a8,1x,i4,4(1x,f12.5),1x,f8.3,1x,i3)') cstn(i), kx(i), alon(i), alat(i), &
        dlon(i), dlat(i), tim(i), grp(i)
enddo
close(20)

print '(a,i8,a,i8)', 'obs read:', nob, '   after duplicate removal:', ncross
do n = 0, ngroups
   print '(a,i3,a,i8)', 'group', n, ': ', count(grp == n)
enddo
print '(a,i8)', 'duplicates removed (group -1): ', count(grp == -1)
end program hilbert_cv
