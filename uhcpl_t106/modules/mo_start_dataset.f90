MODULE mo_start_dataset

  IMPLICIT NONE

  ! U. Schlese, DKRZ, December 1994, original source

  INTEGER :: nstart   !  time step for start/restart.
  INTEGER :: nstep    !  current time step.
  INTEGER :: ntwodt   !  2 * time step interval (milliseconds).
  INTEGER :: neps     !  1000 * time filter coefficient.
  INTEGER :: nf1a     !  `
  INTEGER :: nf3a     !   > indices for fourier buffers.
  INTEGER :: nf4a     !  '
  INTEGER :: ng1a     !  `
  INTEGER :: ng2a     !   > indices for grid point buffers.
  INTEGER :: ng3a     !  
  INTEGER :: ng4a     !  '
  INTEGER :: ngl1a    !  index for grid point slt buffer.
  INTEGER :: ndiahdf  !  logical unit for hdiff diagnostics.
  INTEGER :: nemi     !  logical unit for surface emission file

  INTEGER :: nice     !  logical unit for AMIP2 sea ice file
  INTEGER :: na2stat  !  logical unit for AMIP2 global statistics file
  INTEGER :: na2stre  !  logical unit for AMIP2 statistics rerun file

  INTEGER :: nini     !  logical unit for tracer initial file
  INTEGER :: nisp     !  logical unit  for initial spectral fields.
  INTEGER :: nigp     !  logical unit  for initial grid point fields.
  INTEGER :: nfl1     !  logical unit for optional file read at nstep=0
  INTEGER :: nfl2     !  logical unit for optional file read at nresum
  INTEGER :: nist     !  logical unit for surf.temp. file
  INTEGER :: nhf1     !  `
  INTEGER :: nhf3     !   > logical units for fourier history files.
  INTEGER :: nhf4     !  '
  INTEGER :: nhg1     !  `
  INTEGER :: nhg2     !   >logical units for grid point history files.
  INTEGER :: nhg3     !  
  INTEGER :: nhg4     !  '
  INTEGER :: nl1a     !  index for legendre coefficients buffer
  INTEGER :: nhgl1    !  logical unit for grid point slt work file
  INTEGER :: ngribs   !  logical unit for spectral grib file
  INTEGER :: ngribg   !  logical unit for gridpoint grib file
  INTEGER :: ngribx   !  logical unit for g4x grib file
  INTEGER :: njin     !  logical unit for "jobn" input file
  INTEGER :: njout    !  logical unit for "subjobn" output file
  INTEGER :: ndstart  !   start date of a run
  INTEGER :: ntstart  !   start time of a run
  INTEGER :: ntimeadj !  time offset adjustment in seconds

  LOGICAL :: lres     !  .TRUE. if forecast is restarted.
  LOGICAL :: ldebugio !  .TRUE. to debug IO
  LOGICAL :: ldebugmem!  .TRUE. to debug memory
  LOGICAL :: ly365    !  .TRUE. use 365 days
  LOGICAL :: lg3force   !  .TRUE. to force read from g3b if g3a not in restart
  LOGICAL :: lg3setzero !  .TRUE. to set G3 vars to 0. if not in restart

CONTAINS

  SUBROUTINE initsd

    ! Description:
    !
    ! Sets unit numbers and buffer indices.
    !
    ! Method:
    !
    ! Variables are assigned preset values which may be modified
    ! by namelist *sdsctl*.
    !
    ! Authors:
    !
    ! U. Schlese, DKRZ, December 1994, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! L. Kornblueh, MPI, June 1999, parallel version (MPI based)
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_mpi,     ONLY: p_parallel, p_parallel_io, p_io, p_bcast
    USE mo_doctor,  ONLY: nin

    INCLUDE 'sdsctl.inc'


    !  Executable statements 

    !-- 1. Pre-set values

    !-- 1.1 Basic constants

    nstart = 0
    nstep  = 0
    ntwodt = 0
    neps   = 0

    !-- 1.2 Logical units for files and buffer indices

    ! Subjob files "jobn" and "subjobn"

    njin  = 30
    njout = 39

    ! *Fourier work file and buffers.

    nf1a = 26
    nf3a = 43
    nf4a = 44

    ! Grid point work file and buffers.

    ng1a = 27
    ng2a = 62
    ng3a = 63
    ng4a = 64

    ! Gridpoint buffer index for slt-variables

    ngl1a = 28

    ! Spectral and grid point initial files.

    nisp = 23
    nigp = 24

    ! Legendre polynomial coefficients buffer index

    nl1a = 25

    ! Surface temperature annual cycle file

    nist = 20

    ! Sea ice file for AMIP2

    nice = 96

    ! File for surface emissions

    nemi = 15

    ! Optional files

    nini = 12
    nfl1 = 13
    nfl2 = 14

    ! Grib output files

    ngribs = 29
    ngribg = ngribs
    ngribx = 16

    ! History files.

    nhf1 = 31
    nhf3 = 33
    nhf4 = 34

    nhgl1 = 32

    nhg1 = 35
    nhg2 = 36
    nhg3 = 37
    nhg4 = 38

    ! File for diagnostics of horizontal diffusion

    ndiahdf = 11

    !-- 1.3 Initial values for logical flags

    lres = .TRUE.

    !-- 1.4 Optional variables for user defined purposes

    ldebugio  = .FALSE.
    ldebugmem = .FALSE.

    ! climate mode is 360 days per year (30 day months)
    ly365  = .FALSE.

    lg3force    = .FALSE.
    lg3setzero  = .FALSE.

    ndstart = 1
    ntstart = 0

    ! time offset adjustment in seconds
    ! ntimeadj = -1 : no adjustment
    ntimeadj = -1

    !-- 2. Read namelist sdsctl

    IF (p_parallel) THEN
       IF (p_parallel_io) THEN

!qbao====================================================c
!qbao  modify for kuroshio (apcc2) machine 11/29/2005

       open (nin,file='NAMELIST',&
                 & status='old')
!qbao====================================================c
       READ (nin,sdsctl)

       ENDIF
       CALL p_bcast (nstart, p_io)
       CALL p_bcast (nstep, p_io)
       CALL p_bcast (nf1a, p_io)
       CALL p_bcast (ng1a, p_io)
       CALL p_bcast (ngl1a, p_io)
       CALL p_bcast (nl1a, p_io)
       CALL p_bcast (nisp, p_io)
       CALL p_bcast (nigp, p_io)
       CALL p_bcast (ndiahdf, p_io)
       CALL p_bcast (nemi, p_io)
       CALL p_bcast (nist, p_io)
       CALL p_bcast (nhf1, p_io)
       CALL p_bcast (nhf3, p_io)
       CALL p_bcast (nhf4, p_io)
       CALL p_bcast (nhg1, p_io)
       CALL p_bcast (nhg2, p_io)
       CALL p_bcast (nhg3, p_io)
       CALL p_bcast (nhg4, p_io)
       CALL p_bcast (nhgl1, p_io)
       CALL p_bcast (nfl1, p_io)
       CALL p_bcast (nfl2, p_io)
       CALL p_bcast (nini, p_io)
       CALL p_bcast (ngribs, p_io)
       CALL p_bcast (ngribg, p_io)
       CALL p_bcast (ngribx, p_io)
       CALL p_bcast (njin, p_io)
       CALL p_bcast (njout, p_io)
       CALL p_bcast (ntimeadj, p_io)
       CALL p_bcast (lres, p_io)
       CALL p_bcast (ldebugio, p_io)
       CALL p_bcast (ldebugmem, p_io)
       CALL p_bcast (ly365, p_io)
       CALL p_bcast (lg3force, p_io)
       CALL p_bcast (lg3setzero, p_io)
       CALL p_bcast (ndstart, p_io)
       CALL p_bcast (ntstart, p_io)
    ELSE
       READ (nin,sdsctl)
    ENDIF

    IF (lg3setzero.AND..NOT.lg3force) lg3force = .TRUE.

  END SUBROUTINE initsd

END MODULE mo_start_dataset
