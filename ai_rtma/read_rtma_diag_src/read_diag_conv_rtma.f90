PROGRAM read_diag_conv
!
!  This program is to show how to
!  read GSI diagnositic file for conventional data, which are
!  generated from subroutine:
!      setupps.f90
!      setupt.f90
!      setupq.f90
!      setuppw.f90
!      setupuv.f90
!      setupsst.f90
!      setupgps.f90
!
!  For example in setupt.f90:
!      the arrary contents disgnosis information is rdiagbuf.
!        cdiagbuf(ii)       ! station id
!        rdiagbuf(1,ii)     ! observation type
!        rdiagbuf(2,ii)     ! observation subtype
!        rdiagbuf(3,ii)     ! observation latitude (degrees)
!        rdiagbuf(4,ii)     ! observation longitude (degrees)
!        rdiagbuf(5,ii)     ! station elevation (meters)
!        rdiagbuf(6,ii)     ! observation pressure (hPa)
!        rdiagbuf(7,ii)     ! observation height (meters)
!        rdiagbuf(8,ii)     ! obs time (hours relative to analysis time)
!        rdiagbuf(9,ii)     ! input prepbufr qc or event mark
!        rdiagbuf(10,ii)    ! setup qc or event mark (currently qtflg only)
!        rdiagbuf(11,ii)    ! read_prepbufr data usage flag
!        rdiagbuf(12,ii)    ! analysis usage flag (1=use, -1=not used)
!        rdiagbuf(13,ii)    ! nonlinear qc relative weight
!        rdiagbuf(14,ii)    ! prepbufr inverse obs error (K**-1)
!        rdiagbuf(15,ii)    ! read_prepbufr inverse obs error (K**-1)
!        rdiagbuf(16,ii)    ! final inverse observation error (K**-1)
!        rdiagbuf(17,ii)    ! temperature observation (K)
!        rdiagbuf(18,ii)    ! obs-ges used in analysis (K)
!        rdiagbuf(19,ii)    ! obs-ges w/o bias correction (K) (future slot)
!
!  It is written out as:
!     write(7)'  t',nchar,nreal,ii,mype
!     write(7)cdiagbuf(1:ii),rdiagbuf(:,1:ii)
!

  use kinds, only: r_kind,r_single,i_kind
  use nc_diag_write_mod, only: nc_diag_init, nc_diag_header, nc_diag_metadata, &
       nc_diag_write, nc_diag_metadata_to_single

  implicit none

  real(r_kind) tiny_r_kind
!
! read in variables
!
  character(8),allocatable,dimension(:):: cdiagbuf,cprvstg,csprvstg
  character(8)::cprovider,csubprovider
  real(r_single),allocatable,dimension(:,:)::rdiagbuf
  integer(i_kind) nchar,nreal,ii,mype
  integer(i_kind) idate
  character(180) :: infilename        ! file from GSI running directory
  character(180) :: outfilename       ! file name saving results
  logical :: r15887_diagfile_fmt
!
! output variables
!
  character(len=3)  :: var
  character(8)      :: var8           ! var padded to 8 chars for ncdiag
  real :: rhgt,rlat,rlon,rprs,robs1,rdpt1,robs2,rdpt2,ruse,rerr
  real :: rdhr, ddiff,rusev,factw
  character(8) :: stationID
  integer :: itype,iuse,slm
  real, parameter :: missing = -9.99e9
!
!  misc.
!
  character ::  ch
  integer :: i,j,k,ios
  integer :: ic, iflg

  r15887_diagfile_fmt=.true.

!
  outfilename='diag_results.nc4'
  call getarg(1,infilename)
!
  OPEN (17,FILE=trim(infilename),STATUS='OLD',IOSTAT=ios,ACCESS='SEQUENTIAL',  &
             FORM='UNFORMATTED')
     if(ios > 0 ) then
       write(*,*) ' file is unavailabe: ', trim(infilename)
       stop 123
     endif
     read(17, ERR=999,iostat=ios) idate
     write(*,*) 'process date: ',idate

  call nc_diag_init(trim(outfilename))
  call nc_diag_header("date_time", idate)

   do
     read(17, ERR=999,iostat=ios,end=110) var, nchar,nreal,ii,mype
!     write(*,*) var, nchar,nreal,ii,mype
     if (ii > 0) then
          allocate(cdiagbuf(ii),rdiagbuf(nreal,ii),cprvstg(ii),csprvstg(ii))
          read(17,ERR=999,iostat=ios,end=110) cdiagbuf, rdiagbuf
          read(17)cprvstg,csprvstg
          do i=1,ii
             cprovider=cprvstg(i)
             csubprovider=csprvstg(i)
             itype=rdiagbuf(1,i)    ! observation type
             rlat=rdiagbuf(3,i)     ! observation latitude (degrees)
             rlon=rdiagbuf(4,i)     ! observation longitude (degrees)
             rprs=rdiagbuf(6,i)     ! observation pressure (hPa)
             rhgt=rdiagbuf(7,i)     ! observation height (meters)
             rdhr=rdiagbuf(8,i)     ! obs time (hours relative to analysis time)
             iuse=int(rdiagbuf(12,i))    ! analysis usage flag (1=use, -1=monitoring ) from setup
             rusev=rdiagbuf(11,i)    ! analysis usage flag ( value ) from
                                          !  read_prepbufr
             ddiff=rdiagbuf(18,i)   ! obs-ges used in analysis (K)
             rerr = 0
             if (rdiagbuf(16,i)  > 0) then   ! final inverse observation error (K**-1)
               rerr=1.0/rdiagbuf(16,i)
             end if
             robs1=rdiagbuf(17,i)    !  observation (K)
             rdpt1=rdiagbuf(18,i)    !  obs-ges used in analysis

             stationID = cdiagbuf(i)

!           Remove odd spaces in the station, provider, and subprovider names
             iflg = 0
             do ic=8,1,-1
              ch = stationID(ic:ic)
              if (ch > ' ' .and. ch <= 'z') then
                iflg = 1
              else
                 stationID(ic:ic) = ' '
              end if
              if (stationID(ic:ic) == ' '  .and. iflg == 1) then
                 stationID(ic:ic) = '_'
              endif
             enddo

             iflg = 0
             do ic=8,1,-1
              ch = cprovider(ic:ic)
              if (ch > ' ' .and. ch <= 'z') then
                iflg = 1
              else
                 cprovider(ic:ic) = ' '
              end if
              if (cprovider(ic:ic) == ' '  .and. iflg == 1) then
                 cprovider(ic:ic) = '_'
              endif
             enddo

             iflg = 0
             do ic=8,1,-1
              ch = csubprovider(ic:ic)
              if (ch > ' ' .and. ch <= 'z') then
                iflg = 1
              else
                 csubprovider(ic:ic) = ' '
              end if
              if (csubprovider(ic:ic) == ' '  .and. iflg == 1) then
                 csubprovider(ic:ic) = '_'
              endif
             enddo
!
!   When the data is q, unit convert kg/kg -> g/kg **/
             if (trim(adjustl(var)) == "q") then
                robs1 = robs1 * 1000.0
                rdpt1 = rdpt1 * 1000.0
                rerr = rerr * 1000.0
                !Carley bug fix to DTC version
                ddiff = ddiff * 1000.0
             end if
!   When the data is pw, replase the rprs to -999.0 **/
             if (trim(adjustl(var)) == "pw") rprs=-999.0
!
             if(robs1 > 1.0e8) then
               robs1=-99999.9
               ddiff=-99999.9
             endif
         if (cprovider(1:4)=='B7Hv' .or. cprovider(5:8)=='vH7B') then !this provider name comes with strange characters
             cprovider(1:4)='B7Hv'
             cprovider(5:8)='   '
         endif

         if (csubprovider(1:4)=='B7Hv' .or. csubprovider(5:8)=='vH7B') then
             csubprovider(1:4)='B7Hv'
             csubprovider(5:8)='   '
         endif

         if (itype==154 .and. trim(adjustl(var))=='tca') then
             stationID='GOESSKY'
             cprovider='GOESSKY'
             csubprovider='GOESSKY'
         end if

         if (trim(cprovider)=='') cprovider='EMPTY'
         if (trim(csubprovider)=='') csubprovider='EMPTY'

         ! If we have ceiling or total cloud amount obs less than 0, set to
         ! missing
         if ((trim(adjustl(var))=="tca" .or. trim(adjustl(var))=="cei" .or.   &
              trim(adjustl(var))=="vis" .or. trim(adjustl(var))=="gst") .and. &
              robs1 < 0.) then
           robs1=0.10000E+10
           ddiff=0.10000E+10
         end if

         ! Find the dominant surface type: 0 = water; 1 = land
         if (trim(adjustl(var))=="mitm" .or. trim(adjustl(var))=="mxtm" .or. trim(adjustl(var))=="hwv") then
             slm=rdiagbuf(20,i)
         else if (trim(adjustl(var))=="t" .or. trim(adjustl(var))=="ps" .or. trim(adjustl(var))=="gst" .or. trim(adjustl(var))=="vis" .or. trim(adjustl(var))=="cei" .or. trim(adjustl(var))=="wst" .or. trim(adjustl(var))=="tca") then
             slm=rdiagbuf(21,i)
         else if (trim(adjustl(var))=="q") then
             slm=rdiagbuf(22,i)
         else if (trim(adjustl(var))=="uv") then
             slm=rdiagbuf(26,i)
         else
             slm=-999
         end if

         ! Set wind factor and v-component fields; use missing for non-wind obs
         if (trim(adjustl(var)) .eq. "wst" .OR. trim(adjustl(var)) .eq. "gst") then
             factw  = rdiagbuf(20,i)
             robs2  = missing
             rdpt2  = missing
         else if (trim(adjustl(var)) .eq. "uv") then
             factw  = rdiagbuf(23,i)
             robs2  = rdiagbuf(20,i)
             rdpt2  = rdiagbuf(21,i)
         else
             factw  = missing
             robs2  = missing
             rdpt2  = missing
         end if

         ! Left-justify var into 8-char field for ncdiag
         var8 = adjustl(var)

         ! Write observation metadata to NetCDF diagnostic file
         call nc_diag_metadata("Variable",                      var8        )
         call nc_diag_metadata("Station_ID",                    stationID   )
         call nc_diag_metadata("Provider_Name",                 cprovider   )
         call nc_diag_metadata("Subprovider_Name",              csubprovider)
         call nc_diag_metadata("Observation_Type",              itype       )
         call nc_diag_metadata("Dominant_Sfc_Type",             slm         )
         call nc_diag_metadata("Time",                          rdhr        )
         call nc_diag_metadata("Latitude",                      rlat        )
         call nc_diag_metadata("Longitude",                     rlon        )
         call nc_diag_metadata("Pressure",                      rprs        )
         call nc_diag_metadata("Height",                        rhgt        )
         call nc_diag_metadata("Analysis_Use_Flag",             iuse        )
         call nc_diag_metadata("Observation",                   robs1       )
         call nc_diag_metadata("Obs_Minus_Forecast_adjusted",   ddiff       )
         call nc_diag_metadata("Observation_Error",             rerr        )
         call nc_diag_metadata("Prep_Use_Flag",                 rusev       )
         call nc_diag_metadata("Wind_Factor",                   factw       )
         call nc_diag_metadata("V_Observation",                 robs2       )
         call nc_diag_metadata("V_Obs_Minus_Forecast_adjusted", rdpt2       )

          enddo   ! i  end for one station

          deallocate(cdiagbuf,rdiagbuf,cprvstg,csprvstg)
     else
        read(17,end=110)
        if (r15887_diagfile_fmt) read(17,end=110)
     endif
     end do
110  continue

    close(17)
    call nc_diag_write

  STOP 0

999   PRINT *,'error read in diag file. IOSTAT=',ios
      stop

END PROGRAM read_diag_conv
