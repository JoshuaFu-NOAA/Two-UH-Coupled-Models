MODULE mo_parameters

  IMPLICIT NONE

  ! parameters controlling array sizes.

  INTEGER, PARAMETER :: jpm     = 106 ! max zonal wave number
  INTEGER, PARAMETER :: jpn     = 106 ! max meridional wave number for m=0
  INTEGER, PARAMETER :: jpk     = 106 ! max meridional wave number
  INTEGER, PARAMETER :: jpgl    = 160 ! number of gaussian latitudes
  INTEGER, PARAMETER :: jpnlon  = 320 ! max number of points on each latitude
  INTEGER, PARAMETER :: jpnlev  = 999 ! number of vertical levels.
                                      ! for ntrn only, f90 namelist restriction
                                      
  INTEGER, PARAMETER :: jpg3xf  = 99  ! maximum additional *g3 fields allowed
  INTEGER, PARAMETER :: jptrac  = 21  ! maximum number of prognostic tracers
  INTEGER, PARAMETER :: jps     = 2   ! basic slt variables without tracers
                                      ! i.e. humidity and cloud water
  INTEGER, PARAMETER :: jpmp1   = jpm + 1
  INTEGER, PARAMETER :: jphgl   = jpgl/2
  INTEGER, PARAMETER :: jpnlp2  = jpnlon + 2

END MODULE mo_parameters
