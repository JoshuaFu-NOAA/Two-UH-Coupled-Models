MODULE mo_grib

  ! Description:
  !
  ! This module contains the output routines based on GRIB edition 1.
  !
  ! For a detailed description of the contained subroutines look on their
  ! header.
  ! 
  ! Most important new features:
  ! - Use of a code table in section 1 which allows the use of 128 tables
  !   each containg 128 different variables.
  ! - Uses the sub center entry (255 for ECHAM). Center is still 98 for
  !   ECMWF
  ! - The century parameter is correctly used. Now, in 360 day mode, 
  !   25599 years can be simulated before an overflow in the date occures.
  ! - Uses correct dates in 365/366 day mode. Allows for clean forecast and
  !   nudging runs.
  ! - Either the intrinsic codegb5 C function can be used for encoding, or
  !   the ECMWF EMOS library. For this define for compile -DEMOS for CFLAGS
  !   and FFLAGS.
  ! - The output is now real standard !!!
  !
  ! Authors:
  !
  ! L. Kornblueh, MPI, December 1998
  ! L. Kornblueh, MPI, April 1998, added NWP forecast mode
  ! U. Schulzweida, MPI, April 2000, EMOS compatible
  !

  IMPLICIT NONE

  INTEGER :: kleng                      ! size of kgrib
  INTEGER :: klengo                     ! byte size of kgrib after codegb5
  INTEGER, ALLOCATABLE :: kgrib(:)      ! grib data set
             
  INTEGER :: ksec0(2), ksec1(43)
  INTEGER :: ksec2_gp(22), ksec2_sp(22)
  INTEGER :: ksec3(2), ksec4(42)

  REAL, ALLOCATABLE :: psec2(:), psec3(:), psec4(:)

  INTEGER :: kword    ! output size in words

  INTEGER :: klenp    ! size of psec4

  INTEGER, PARAMETER :: nint = 0
  INTEGER, PARAMETER :: nbit = BIT_SIZE(nint)
  INTEGER, PARAMETER :: iobyte = nbit/8

  INTEGER, PARAMETER :: local_table     = 128 !  local code table
  INTEGER, PARAMETER :: nudging_table   = 129 !  nudging code table
  INTEGER, PARAMETER :: g4x_table       = 130 !  G4X code table
  INTEGER, PARAMETER :: center_id       =  98 !  identification of centre
  INTEGER            :: model_id              !  model identification
  INTEGER, PARAMETER :: grid_type       = 255 !  grid definition
  INTEGER, PARAMETER :: nflag           = 128 !  flag(see code table 1)
  INTEGER            :: code_parameter        !  data field code 
                                              !  (see code table 2 )
  INTEGER            :: level_type            !  indicator of type of level 
                                              !  (see code table 3)
  INTEGER            :: level_p1              !  data level1 (see code table 3)
  INTEGER            :: level_p2              !  data level2 (see code table 3)

  ! reference time of data

  INTEGER            :: year                  !  year of century 
  INTEGER            :: month                 !  month 
  INTEGER            :: day                   !  day 
  INTEGER            :: hour                  !  hour
  INTEGER            :: minute                !  minute

  INTEGER            :: time_unit       =   0 ! unit of time range 
                                              ! (see code table 4)
  INTEGER            :: time_p1         =   0 ! time range 1
  INTEGER            :: time_p2         =   0 ! time range 2
  INTEGER            :: range_flag      =  10 ! time range flag
                                              ! (see code table 5)
  INTEGER            :: century               ! century
  INTEGER, PARAMETER :: subcenter_id    = 255 ! subcenter
  INTEGER, PARAMETER :: decimal_scale   =   0 ! decimal scale
  INTEGER, PARAMETER :: local_extension =   0 ! local extension flag

  INTEGER            :: forecast_hours  =   0 ! number of forecast hours in
                                              ! NWP mode

  ! Variables for block 2

  INTEGER            :: npvct         ! no. of vertical ccordinate parameters

  ! type gaussian latitudes ( data_type = 4)

  INTEGER, PARAMETER :: ngptype =   4 ! data representation type (code table 6)
  INTEGER            :: npalat        ! number of points along a latitude      
  INTEGER            :: npamer        ! number of points along a meridional    
  INTEGER            :: nglaor        ! latitude of origin in degree*1000      
  INTEGER, PARAMETER :: ngloor  =   0 ! longitude of origin in degree*1000     
  INTEGER, PARAMETER :: neinfl  = 128 ! resolution flag(code table 7)          
  INTEGER            :: nglaex        ! latitude  of extreme point degree*1000 
  INTEGER            :: ngloex        ! longitude of extreme point degree*1000 
  INTEGER            :: ngdi          ! latitude    increment in  degree*1000  
  INTEGER            :: nglpe         ! numder of latitude lines between a pole
                                      ! and the equator                        
  INTEGER, PARAMETER :: nscmfl  =   0 ! scanning mode(flags see table 8) 

  ! type  ( spherical harmonics data_type = 50)

  INTEGER, PARAMETER :: nsptype =  50 ! data representation type (code table 6)
  INTEGER            :: nfasmj        ! j - pentagonal resolution parameter
  INTEGER            :: nfasmk        ! k - pentagonal resolution parameter 
  INTEGER            :: nfasmm        ! m - pentagonal resolution parameter
  INTEGER, PARAMETER :: nrtsh   =   1 ! representation type (see table 9) 
  INTEGER, PARAMETER :: nrmsh   =   1 ! representation mode (see table 10)

CONTAINS

  SUBROUTINE griberror (name)

    USE mo_exception,  ONLY: finish

    CHARACTER (*) :: name

    CALL finish(name, 'GRIB I/O error - disk full?')

  END SUBROUTINE griberror

  SUBROUTINE set_output_time

    USE mo_constants,     ONLY: dayl
    USE mo_start_dataset, ONLY: nstep
    USE mo_control,       ONLY: ncbase, ntbase, dtime, lnwp
    USE mo_year,          ONLY: ic2ymd
    USE mo_exception,     ONLY: finish

    INTEGER, SAVE :: century0, year0, month0, day0, hour0, minute0

    !  Local scalars: 

    REAL :: zsecs
    INTEGER :: itbase, date

    LOGICAL, SAVE :: l_nwp_started = .FALSE.

    !  Intrinsic functions 
    INTRINSIC MOD


    IF (.NOT. l_nwp_started .OR. .NOT. lnwp ) THEN
       zsecs = (nstep+1)*dtime + ntbase
       date = ic2ymd(ncbase+INT(zsecs/dayl))
       itbase = MOD(zsecs,dayl)
       year = date/10000
       month = MOD(date/100,100)
       day = MOD(date,100)
       hour = itbase/3600
       minute = MOD(itbase/60,60)

       century = year/100+1
       year = MOD(year,100)

       IF (century > 255 .OR. century < 0) THEN
          CALL finish ('set_output_time', &
                       'This century is not supported in GRIB 1')
       END IF
    END IF

    IF (lnwp .AND. .NOT. l_nwp_started) THEN 
       century0 = century
       year0    = year
       month0   = month
       day0     = day
       hour0    = hour
       minute0  = minute-ANINT(dtime/60)   ! adjust for start of NWP forecast
    ENDIF

    IF (lnwp) THEN     ! NWP mode ..., not analysis
       l_nwp_started = .TRUE.
       range_flag = 0 
       time_unit =  1
       century = century0
       year    = year0
       month   = month0
       day     = day0
       hour    = hour0
       minute  = minute0
    ENDIF

    IF (l_nwp_started) THEN
       forecast_hours = ANINT((nstep+1)*dtime/3600)
    END IF

  END SUBROUTINE set_output_time

  SUBROUTINE close_grib_file

    USE mo_post,       ONLY: nunitdf

    INTEGER :: iret

    !  External subroutines 
    EXTERNAL pbclose


    CALL pbclose(nunitdf,iret)

  END SUBROUTINE close_grib_file

  SUBROUTINE open_grib_file

    USE mo_doctor,     ONLY: nout, nerr
    USE mo_post,       ONLY: nunitdf
    USE mo_filename,   ONLY: path_limit, standard_grib_file
    USE mo_control,    ONLY: lnwp
    USE mo_exception,  ONLY: finish

    INTEGER :: iret

    CHARACTER (3) :: fcst
    CHARACTER (path_limit+5) :: filename

    !  External subroutines 
    EXTERNAL pbopen


    IF (lnwp) THEN
       CALL set_output_time
       IF (forecast_hours > 744) THEN
          WRITE (nerr,*) 'NWP mode makes no sense for this time range.'
          WRITE (nerr,*) 'Please change to the climate mode.'
          CALL finish('open_grib_file','Run terminated.')
       END IF
       WRITE (fcst,'(i3.3)') forecast_hours
       filename = TRIM(standard_grib_file)//'+'//fcst
       WRITE (nout,*) 'NWP forecast (', forecast_hours, 'hr) GRIB output: ', &
                       TRIM(filename(INDEX(filename,'/',.TRUE.)+1:)) 
    ELSE
       filename = TRIM(standard_grib_file)
    END IF

    CALL pbopen(nunitdf,filename,'w',iret)
    IF (iret /= 0) THEN
       WRITE (nerr,*) 'Could not open file: ', filename
       CALL finish('open_grib_file','Run terminated.')
    END IF

  END SUBROUTINE open_grib_file

  SUBROUTINE init_grib

    USE mo_constants, ONLY: api
    USE mo_control,   ONLY: nvclev, nlon, ngl, nhgl, nn, nm, nk, lanalysis
#ifdef EMOS
    USE mo_control,   ONLY: vct
#endif
    USE mo_gaussgrid, ONLY: gmu
    USE mo_doctor,    ONLY: nvers
    USE mo_post,      ONLY: nlalo


    ! allocate grib work array

    nlalo = nlon*ngl
    kleng = nlalo + 1000  ! only for max 32 (64)bit packing
                          ! dependend on generic INTEGER size,
                          ! nlalo - packed field size
                          !  1000 - upper limit for header size
    IF (.NOT. ALLOCATED(kgrib)) ALLOCATE(kgrib(kleng))

    ! set all GRIB block entries to 0

    ksec0(:)    = 0
    ksec1(:)    = 0 
    ksec2_gp(:) = 0
    ksec2_sp(:) = 0
    ksec3(:)    = 0
    ksec4(:)    = 0

    ! set model version and reset century

    model_id = nvers
    century  = 0
    level_p1 = 0
    level_p2 = 0

    ! set fixed values in GRIB block 1, common to all GRIB sets

    ! to set standard output to analysis mode lnwp must be .false. and
    ! lanalysis .true.. Later with working adiabatic NMI this must be 
    ! changed to range_flag = 1. range_flag = 0 is default for output
    ! of the OI or 3DVAR analysis, this are given external. In case a 
    ! 4DVAR is running inside ECHAM this value must be adjusted to 
    ! range_flag = 0, as it is now.
    ! For usual climate mode runs range_flag is set to 10, which means
    ! The time given in the GRIB header in section 1 is the valid time
    ! and time_p1 and time_p2 do not mean anything. 

    IF (lanalysis) range_flag = 0    

    ksec1(2)  = center_id
    ksec1(3)  = model_id
    ksec1(4)  = grid_type        ! No WMO predefined
    ksec1(5)  = nflag
    ksec1(8)  = level_p1
    ksec1(9)  = level_p2
    ksec1(15) = time_unit
    ksec1(16) = time_p1
    ksec1(17) = time_p2
    ksec1(18) = range_flag
    ksec1(21) = century          ! preliminary value
    ksec1(22) = subcenter_id
    ksec1(23) = decimal_scale
    ksec1(24) = local_extension  ! Until now no extension, might be used
                                 ! for paleoclimate runs centuries

    npalat = nlon
    npamer = ngl
    nglaor = ASIN(gmu(1))*180000.0/api
    nglaex = -nglaor
    ngdi = 360000/npalat
    ngloex = 360000-ngdi
    nglpe = nhgl

    ksec2_gp(1)  = ngptype
    ksec2_gp(2)  = npalat
    ksec2_gp(3)  = npamer
    ksec2_gp(4)  = nglaor
    ksec2_gp(5)  = ngloor
    ksec2_gp(6)  = neinfl 
    ksec2_gp(7)  = nglaex
    ksec2_gp(8)  = ngloex
    ksec2_gp(9)  = ngdi
    ksec2_gp(10) = nglpe
    ksec2_gp(11) = nscmfl
    ksec2_gp(12) = 2*nvclev

    nfasmj = nn
    nfasmk = nk
    nfasmm = nm

    ksec2_sp(1)  = nsptype
    ksec2_sp(2)  = nfasmj
    ksec2_sp(3)  = nfasmk
    ksec2_sp(4)  = nfasmm
    ksec2_sp(5)  = nrtsh
    ksec2_sp(6)  = nrmsh
    ksec2_sp(12) = 2*nvclev

    ksec3(1:2) = 0  

    ksec4(1:42) = 0 

    ksec4(2) = 24   ! bits, predefined standard in ECHAM4 

#ifdef EMOS
    ALLOCATE (psec2(10+2*nvclev))
    ALLOCATE (psec3(0))    
    psec2(1:10) = 0.0
    psec2(11:10+2*nvclev) = vct(1:2*nvclev)

    ! ECMWF GRIB library does not except centuries < 20, 
    ! so checking has to be disabled (ly365 seems to be 
    ! correct with checking)

    CALL grsvck (0)

#endif

  END SUBROUTINE init_grib

  SUBROUTINE outsp

    ! Description:
    !
    ! Control postprocessing of spectral fields
    !
    ! Method:
    !
    ! Pack and write
    !
    ! *outsp* is called from subroutine *stepon*
    !
    ! Results:
    ! Packed spectral fields are written out
    !
    ! Externals:
    ! *codegb5* - pack to grib code
    !
    ! Authors:
    !
    ! U. Schlese, DKRZ, July 1994, original source
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! L. Kornblueh, MPI, November 1998, f90 rewrite, extended to GRIB 1
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_memory_sp,     ONLY: memory_info, sp, get_entry, get_info
    USE mo_control,       ONLY: nsp, lnwp, nlevp1, nvclev, vct, n2sp, nlev
    USE mo_post,          ONLY: lppt, lppp, lppd, lppvo, nunitdf
    USE mo_doctor,        ONLY: nout
    USE mo_mpi,           ONLY: p_pe, p_io
    USE mo_decomposition, ONLY: dcg => global_decomposition
    USE mo_transpose,     ONLY: gather_sp

    !  Local scalars: 
    INTEGER :: ierr, iret, jlev

    !  Local arrays: 
    REAL, POINTER :: ptr3d(:,:,:)
    REAL, POINTER :: z3d(:,:,:)
    REAL :: psec4(2,nsp)

    TYPE (memory_info) :: info

    !  External subroutines 
    EXTERNAL codegb5, pbwrite


    !  Executable statements 

    NULLIFY (z3d)

    ! Finish definition GRIB block 1

    ksec1(1) = local_table

    CALL set_output_time

    IF (lnwp) THEN
       ksec1(15) = time_unit
       ksec1(16) = forecast_hours
       ksec1(18) = range_flag
    END IF
    ksec1(10) = year
    ksec1(11) = month
    ksec1(12) = day
    ksec1(13) = hour
    ksec1(14) = minute
    ksec1(21) = century

    IF (p_pe == p_io) THEN
       IF (lnwp) THEN
          WRITE (nout,'(a,i2.2,a1,i2.2,2x,i2.2,a1,i2.2,a1,i4.4,a,i3.3,a)') &
               ' Output time: ', hour, ':', minute,                        &
               day, '.', month, '.', year+100*(century-1),                 &
               ', forecast at +', forecast_hours, ' hr'
       ELSE
          WRITE (nout,'(a,i2.2,a1,i2.2,2x,i2.2,a1,i2.2,a1,i4.4)') &
               ' Output time: ', hour, ':', minute,               &
               day, '.', month, '.', year+100*(century-1)
       END IF
    END IF

    IF (lppp .OR. lppt) THEN

       CALL get_entry (sp, 'STP', ptr3d)
       CALL get_info  (sp, 'STP', info)
       ALLOCATE (z3d(info%gdim_1, info%gdim_2, info%gdim_3))
       CALL gather_sp (z3d, ptr3d, dcg)

    END IF

    ! Pressure

    IF (lppp) THEN

       code_parameter = 152
       level_type = 1
       level_p1 = 0

       ksec1(6) = code_parameter
       ksec1(7) = level_type
       ksec1(8) = level_p1

       IF (p_pe == p_io) THEN
          psec4(:,:) = z3d(nlevp1,:,:)
#ifdef EMOS
          ksec4(1) = n2sp
          CALL gribex (ksec0, ksec1, ksec2_sp, psec2, ksec3, psec3, &
               ksec4, psec4, n2sp, kgrib, kleng, kword, 'C', ierr)
#else
          CALL codegb5 (psec4,n2sp,16,nbit,ksec1,ksec2_sp,vct,2*nvclev, &
               kgrib,klengo,kword,0,ierr)
#endif
          CALL pbwrite(nunitdf,kgrib,kword*iobyte,iret)

          IF (iret /= kword*iobyte) CALL griberror('outsp')

       END IF
    END IF

    ! Temperature

    IF (lppt) THEN

       code_parameter = 130
       level_type = 109

       ksec1(6) = code_parameter
       ksec1(7) = level_type

       DO jlev = 1, nlev

          level_p1 = jlev
          ksec1(8) = level_p1

          IF (p_pe == p_io) THEN
             psec4(:,:) = z3d(jlev,:,:)
#ifdef EMOS
             ksec4(1) = n2sp
             CALL gribex (ksec0, ksec1, ksec2_sp, psec2, ksec3, psec3, &
                  ksec4, psec4, n2sp, kgrib, kleng, kword, 'C', ierr)

#else
             CALL codegb5(psec4,n2sp,16,nbit,ksec1,ksec2_sp,vct,2*nvclev, &
                  kgrib,klengo,kword,0,ierr)
#endif
             CALL pbwrite(nunitdf,kgrib,kword*iobyte,iret)

             IF (iret /= kword*iobyte) CALL griberror('outsp')

          END IF
       END DO
    END IF

    IF (ASSOCIATED(z3d)) DEALLOCATE(z3d)

    ! Divergence

    IF (lppd) THEN

       CALL get_entry (sp, 'SD', ptr3d)
       CALL get_info (sp, 'SD', info)
       ALLOCATE (z3d(info%gdim_1, info%gdim_2, info%gdim_3))
       CALL gather_sp (z3d, ptr3d, dcg)

       code_parameter = 155
       level_type = 109

       ksec1(6) = code_parameter
       ksec1(7) = level_type

       DO jlev = 1, nlev

          level_p1 = jlev
          ksec1(8) = level_p1

          IF (p_pe == p_io) THEN
             psec4(:,:) = z3d(jlev,:,:)
#ifdef EMOS
             ksec4(1) = n2sp
             CALL gribex (ksec0, ksec1, ksec2_sp, psec2, ksec3, psec3, &
                  ksec4, psec4, n2sp, kgrib, kleng, kword, 'C', ierr)
#else
             CALL codegb5(psec4,n2sp,16,nbit,ksec1,ksec2_sp,vct,2*nvclev, &
                  kgrib,klengo,kword,0,ierr)
#endif
             CALL pbwrite(nunitdf,kgrib,kword*iobyte,iret)

             IF (iret /= kword*iobyte) CALL griberror('outsp')

          END IF
       END DO
    END IF

    IF (ASSOCIATED(z3d)) DEALLOCATE(z3d)

    ! Vorticity

    IF (lppvo) THEN

       CALL get_entry (sp, 'SVO', ptr3d)
       CALL get_info (sp, 'SVO', info)
       ALLOCATE (z3d(info%gdim_1, info%gdim_2, info%gdim_3))
       CALL gather_sp (z3d, ptr3d, dcg)

       code_parameter = 138
       level_type = 109

       ksec1(6) = code_parameter
       ksec1(7) = level_type

       DO jlev = 1, nlev

          level_p1 = jlev
          ksec1(8) = level_p1

          IF (p_pe == p_io) THEN
             psec4(:,:) = z3d(jlev,:,:)
#ifdef EMOS
             ksec4(1) = n2sp
             CALL gribex (ksec0, ksec1, ksec2_sp, psec2, ksec3, psec3, &
                  ksec4, psec4, n2sp, kgrib, kleng, kword, 'C', ierr)

#else
             CALL codegb5(psec4,n2sp,16,nbit,ksec1,ksec2_sp,vct,2*nvclev, &
                  kgrib,klengo,kword,0,ierr)
#endif
             CALL pbwrite(nunitdf,kgrib,kword*iobyte,iret)

             IF (iret /= kword*iobyte) CALL griberror('outsp')

          END IF
       END DO
    END IF

    IF (ASSOCIATED(z3d)) DEALLOCATE(z3d)

    RETURN

  END SUBROUTINE outsp

  SUBROUTINE outgpi

    ! Description:
    !
    ! Control postprocessing of gridpoint fields.
    !
    ! Method:
    !
    ! 2-dimensional fields are allocated
    ! all values are collected in these field sorted north to south
    ! at completion packing is done by codegb5.
    !
    ! *outgpi* is called from subroutine *scan1sl*
    !
    ! Results:
    ! The selected output fields are written out.
    !
    ! Externals:
    ! *codegb5*  - pack grid point field to grib code
    !
    ! Authors:
    !
    ! U. Schlese, DKRZ, July 1994, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_control,       ONLY: nlon, ngl, nptime, dtime, lnwp, nlevp1, nvclev, vct
    USE mo_post,          ONLY: ncdmin, ncdmax, laccu, yn, npplev, npbits, nlalo, &
                                nunitdf
    USE mo_memory_g3b,    ONLY: memory_info, g3b, get_entry, get_info
    USE mo_mpi,           ONLY: p_pe, p_io
    USE mo_decomposition, ONLY: dcg => global_decomposition
    USE mo_transpose,     ONLY: gather_gp
    USE mo_time_control,  ONLY: time_inc_sec, TIME_INC_DAYS
    USE mo_doctor,        ONLY: nout, nerr

    !  Local scalars: 
    REAL :: quot
    INTEGER :: ierr, ilev, ipbits, iret, jfield, jlev
    CHARACTER (8) :: yname

    !  Local arrays: 

    ! pointers for grid point space
    REAL, POINTER :: ptr3d(:,:,:)
    REAL, POINTER :: z3d(:,:,:)
    REAL :: psec4(nlon,ngl)

    TYPE (memory_info) :: info

    !  External subroutines 
    EXTERNAL codegb5, pbwrite

    !  Intrinsic functions 
    INTRINSIC MOD


    !  Executable statements 

    ! Definition grib blocks

!    quot = 1.0/(nptime*dtime)
    quot = 1.0/time_inc_sec(nptime, TIME_INC_DAYS)

    ! Define grib block 1

    ksec1(1) = local_table

    CALL set_output_time

    IF (lnwp) THEN
       ksec1(15) = time_unit
       ksec1(16) = forecast_hours
       ksec1(18) = range_flag
    END IF
    ksec1(10) = year
    ksec1(11) = month
    ksec1(12) = day
    ksec1(13) = hour
    ksec1(14) = minute
    ksec1(21) = century

    ! Loop over all fields

    DO jfield = ncdmin, ncdmax

       yname = yn(jfield)
       IF (yname(1:1)==' ') CYCLE
       IF (yname(1:3)=='G4X') CYCLE

       ipbits = npbits(jfield)

       code_parameter = jfield

       ksec1(6) = code_parameter 

       CALL get_entry (g3b, yname, ptr3d)
       CALL get_info (g3b, yname, info)
       ALLOCATE (z3d(info%gdim_1, info%gdim_2, info%gdim_3))
       CALL gather_gp (z3d, ptr3d, dcg)

       IF (p_pe == p_io) THEN

          DO jlev = 1, nlevp1
             ilev = npplev(jlev,jfield)

             IF (ilev == 0) EXIT

             ! Find first puzzle piece

             IF (ilev == -100) THEN
                IF (SIZE(z3d,3) == 1) THEN
                   psec4(1:nlon,:) =  z3d(1:nlon,:,1)
                ELSE IF (SIZE(z3d,2) == 1) THEN
                   psec4(1:nlon,:) =  z3d(1:nlon,1,:)
                END IF
             ELSE
                psec4(1:nlon,:) =  z3d(1:nlon,ilev,:)
             END IF
             IF (laccu(jfield)) psec4(:,:) = psec4(:,:)*quot

             ! Define variable parts of grib block 1

             ! Due to the limited layer information possibilities in
             ! GRIB 1 the soil depth is given as the height of the
             ! layer center (level_type 111).  

             IF (ilev == -100) THEN
                IF (code_parameter == 167 .OR. &  ! 2 m temperature
                     code_parameter == 168 .OR. & ! 2 m dew point
                     code_parameter == 201 .OR. & ! 2 m maximum temperature 
                     code_parameter == 202 .OR. & ! 2 m minimum temperature 
                     code_parameter == 216) THEN  ! 2 m maximum wind speed
                   level_type = 105
                   level_p1   = 2
                ELSE IF (code_parameter == 165 .OR. & ! 10 m u
                     code_parameter == 166 .OR. &     ! 10 m v
                     code_parameter == 171) THEN      ! 10 m wind speed
                   level_type = 105
                   level_p1   = 10
                ELSE IF (code_parameter == 207) THEN ! TD3 
                   level_type = 111  
                   level_p1   = 3
                ELSE IF (code_parameter == 208) THEN ! TD4
                   level_type = 111  
                   level_p1   = 19
                ELSE IF (code_parameter == 209) THEN ! TD5
                   level_type = 111  
                   level_p1   = 78
                ELSE IF (code_parameter == 170) THEN ! TD
                   level_type = 111  
                   level_p1   = 268
                ELSE IF (code_parameter == 183) THEN ! TDCL
                   level_type = 111  
                   level_p1   = 698
                ELSE
                   level_type = 1
                   level_p1   = 0
                ENDIF
             ELSE
                level_type = 109
                level_p1   = ilev
             END IF

             ksec1(7) = level_type
             ksec1(8) = level_p1

             ! Pack and write
#ifdef EMOS
             ksec4(1) = nlalo
             ksec4(2) = ipbits
             CALL gribex (ksec0, ksec1, ksec2_gp, psec2, ksec3, psec3, &
                  ksec4, psec4, nlalo, kgrib, kleng, kword, 'C', ierr)

#else
             CALL codegb5(psec4,nlalo,ipbits,nbit,ksec1,ksec2_gp,vct, &
                  2*nvclev, kgrib,klengo,kword,0,ierr)
#endif
             CALL pbwrite(nunitdf,kgrib,kword*iobyte,iret)

             IF (iret /= kword*iobyte) CALL griberror('outgpi')

          END DO
       END IF

       DEALLOCATE (z3d)

    END DO

  END SUBROUTINE outgpi

  SUBROUTINE outgpli

    ! Description:
    !
    ! Controls postprocessing of semi lagrangian gridpoint fields.
    !
    ! Method:
    !
    ! 2-dimensional fields are allocated
    ! all values are collected in these field sorted north to south
    ! at completion packing is done by codegb
    !
    ! *outgpli* is called from subroutine *scan1sl*
    !
    ! Results:
    ! The selected output fields are written out
    !
    ! Externals:
    ! *codegb5*  - pack grid point field to grib code
    !
    ! Authors:
    !
    ! U. Schlese, DKRZ, July 1994, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_control,       ONLY: nlon, ngl, nlp2, lnwp, nlev, nvclev, vct
    USE mo_post,          ONLY: lppq, lppx, nunitdf, nlalo
    USE mo_tracer,        ONLY: jptrac, jps, ntrac, lppxt, ntcode
    USE mo_memory_gl,     ONLY: memory_info, gl, get_entry, get_info
    USE mo_mpi,           ONLY: p_pe, p_io
    USE mo_decomposition, ONLY: dcg => global_decomposition
    USE mo_transpose,     ONLY: gather_gp

    !  Local scalars: 
    INTEGER :: ierr, iret, ivar, jfield, jlev, jt
    CHARACTER (8) :: yname

    !  Local arrays: 

    ! pointers for grid point space
    REAL, POINTER :: ptr3d(:,:,:), ptr4d(:,:,:,:)
    REAL, POINTER :: z3d(:,:,:), z4d(:,:,:,:)
    REAL :: psec4(nlon,ngl)

    TYPE (memory_info) :: info

    INTEGER :: icodem(jptrac+jps), illen(jptrac+jps)

    !  External subroutines 
    EXTERNAL codegb5, pbwrite


    !  Executable statements 

    NULLIFY (z4d)

    ! Create tables

    ivar = 0

    ! q-field

    IF (lppq) THEN
       ivar = ivar + 1
       icodem(ivar) = 133
       illen(ivar)  = nlp2
    END IF

    ! x-field

    IF (lppx) THEN
       ivar = ivar + 1
       icodem(ivar) = 153
       illen(ivar)  = nlp2
    END IF

    ! Additional tracers

    DO jt = 1, ntrac
       IF (lppxt(jt)) THEN
          ivar = ivar + 1
          icodem(ivar) = ntcode(jt)
          illen(ivar)  = nlon
       END IF
    END DO

    ! Define grib blocks

    ! Define grib block 1

    ksec1(1) = local_table

    level_type = 109
    ksec1(7)  = level_type

    CALL set_output_time

    IF (lnwp) THEN
       ksec1(15) = time_unit
       ksec1(16) = forecast_hours
       ksec1(18) = range_flag
    END IF
    ksec1(10) = year
    ksec1(11) = month
    ksec1(12) = day
    ksec1(13) = hour
    ksec1(14) = minute
    ksec1(21) = century

    ! Loop over all fields

    DO jfield = 1, ivar

       code_parameter = icodem(jfield)
       ksec1(6) = code_parameter

       IF (jfield == 1) yname = 'Q'
       IF (jfield == 2) yname = 'X'
       IF (jfield >= 3) yname = 'XT'

       IF (yname(1:2) == 'XT') THEN
          IF (jfield == 3) THEN
             CALL get_entry (gl, yname, ptr4d)
             CALL get_info (gl, yname, info)
             ALLOCATE (z4d(info%gdim_1, info%gdim_2, info%gdim_3, info%gdim_4))
             CALL gather_gp (z4d, ptr4d, dcg)
          END IF
          ALLOCATE (z3d(info%gdim_1, info%gdim_2, info%gdim_4))
          z3d(:,:,:) = z4d(:,:,jfield-jps,:)
       ELSE   
          CALL get_entry (gl, yname, ptr3d)
          CALL get_info (gl, yname, info)
          ALLOCATE (z3d(info%gdim_1, info%gdim_2, info%gdim_3))
          CALL gather_gp (z3d, ptr3d, dcg)
       END IF

       IF (p_pe == p_io) THEN

          DO jlev = 1, nlev

             level_p1  = jlev
             ksec1(8)  = level_p1

             psec4(1:nlon,:) =  z3d(1:nlon,jlev,:)

             ! Pack and write
#ifdef EMOS
             ksec4(1) = nlalo
             CALL gribex (ksec0, ksec1, ksec2_gp, psec2, ksec3, psec3, &
                  ksec4, psec4, nlalo, kgrib, kleng, kword, 'C', ierr)
#else
             CALL codegb5(psec4,nlalo,16,nbit,ksec1,ksec2_gp,vct,2*nvclev, &
                  kgrib,klengo,kword,0,ierr)
#endif
             CALL pbwrite(nunitdf,kgrib,kword*iobyte,iret)

             IF (iret /= kword*iobyte) CALL griberror('outgpli')

          END DO
       END IF
       DEALLOCATE (z3d)
    END DO

    IF (ASSOCIATED (z4d)) DEALLOCATE (z4d)

  END SUBROUTINE outgpli

  SUBROUTINE outgpx

    ! Description:
    !
    ! Control postprocessing of g4x-gridpoint fields.
    !
    ! Method:
    !
    ! 2-dimensional fields are allocated
    ! all values are collected in these field sorted north to south
    ! at completion packing is done by codegb
    !
    ! *outgpx* is called from subroutine *scan1sl*
    !
    ! Results:
    ! The selected output fields are written out
    !
    ! Externals:
    ! *codegb*  - pack grid point field to grib code
    !
    ! Authors:
    !
    ! U. Schlese, DKRZ, July 1994, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! 
    ! for more details see file AUTHORS
    !

!!$    USE mo_parameters
!!$    USE mo_start_dataset
!!$    USE mo_control
!!$    USE mo_post
!!$    USE mo_year
!!$
!!$    !  Local scalars: 
!!$    REAL :: quot
!!$    INTEGER :: iend, ierr, ig4base, ilev, &
!!$&              ipbits, iret, irow, istart, ix, jfield, jlev
!!$    CHARACTER (8) :: yname
!!$
!!$    !  Local arrays: 
!!$    ! pointers for grid point space
!!$    REAL, POINTER :: z1dim(:)
!!$    REAL :: psec4(nlon,ngl)
!!$    
!!$    !  External subroutines 
!!$    EXTERNAL codegb5, pbwrite
!!$    
!!$    !  Intrinsic functions 
!!$    INTRINSIC INT, MOD
!!$    
!!$
!!$    !  Executable statements   
!!$  
!!$    ! Define grib block 1
!!$
!!$    ksec1(1) = g4x_table
!!$
!!$    CALL set_output_time
!!$
!!$    IF (lnwp) THEN
!!$       ksec1(15) = time_unit
!!$       ksec1(16) = forecast_hours
!!$       ksec1(18) = range_flag
!!$    END IF 
!!$    ksec1(10) = year
!!$    ksec1(11) = month
!!$    ksec1(12) = day
!!$    ksec1(13) = hour
!!$    ksec1(14) = minute
!!$    ksec1(21) = century
!!$
!!$    ! Loop over all g4x-fields
!!$    
!!$    DO jfield = 1, 128
!!$       yname = yn(jfield)
!!$       IF (yname(1:1)==' ') CYCLE
!!$       IF (yname(1:3)/='G4X') CYCLE
!!$       ipbits = npbits(jfield)
!!$       
!!$       DO jlev = 1, nlevp1
!!$          ilev = npplev(jlev,jfield)
!!$          IF (ilev==0) EXIT
!!$          
!!$          ! Find first puzzle piece
!!$          
!!$          ig4base = m_offset_g4b
!!$          istart = ig4base
!!$          IF (ilev==0) THEN
!!$             istart = istart + (ilev-1)*nlon
!!$          END IF
!!$          iend = istart + nlon - 1
!!$          z1dim => mptglo(ng1a) %vp(istart:iend)
!!$          
!!$          ! Collect data
!!$          
!!$          DO irow = 1, ngl
!!$             IF (MOD(irow,2)==1) THEN
!!$                ix = (irow+1)/2
!!$             ELSE
!!$                ix = ngl + 1 - irow/2
!!$             END IF
!!$             
!!$             IF (laccu(jfield)) THEN
!!$                quot = 1.0/(n4ptime*dtime)
!!$                psec4(1:nlon,ix) = z1dim(1:nlon)*quot
!!$             ELSE
!!$                psec4(1:nlon,ix) = z1dim(1:nlon)
!!$             END IF
!!$             istart = istart + nlon
!!$             iend = iend + nlon - 1
!!$             z1dim => mptglo(ng1a) %vp(istart:iend)
!!$          END DO
!!$          
!!$          ! Pack and write
!!$          
!!$          code_parameter = jfield
!!$          IF (ilev==-100) THEN
!!$             level_type = 1
!!$             level_p1   = 0
!!$          ELSE
!!$             level_type = 109
!!$             level_p1   = ilev
!!$          END IF
!!$
!!$          ksec1(7) = level_type
!!$          ksec1(8) = level_p1
!!$
!!$#ifdef EMOS
!!$          ksec4(1) = nlalo
!!$          CALL gribex (ksec0, ksec1, ksec2_gp, psec2, ksec3, psec3, &
!!$&                      ksec4, psec4, nlalo, kgrib, kleng, kword, 'C', ierr)
!!$
!!$#else
!!$          CALL codegb5(psec4,nlalo,ipbits,nbit,ksec1,ksec2_gp,vct,2*nvclev,&
!!$&                      kgrib,klengo,kword,0,ierr)
!!$#endif
!!$          CALL pbwrite(nunitg4x,kgrib,kword*iobyte,iret)
!!$
!!$          IF (iret /= kword*iobyte) CALL griberror('outgpx')
!!$          
!!$       END DO
!!$    END DO

  END SUBROUTINE outgpx

END MODULE mo_grib
