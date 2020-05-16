MODULE mo_constants

  IMPLICIT NONE

  ! ----------------------------------------------------------------
  !
  ! module *mo_constants* basic universal constants  and derived constants.
  !
  ! ----------------------------------------------------------------

  REAL, PARAMETER :: dayl = 86400.  ! length of the day (in seconds).

  REAL :: api     ! 2.*arcsin(1.).
  REAL :: a       ! radius of the earth.
  REAL :: omega   ! solid rotation velocity of the earth.
  REAL :: g       ! gravity acceleration.
  REAL :: cpd     ! specific heat at constant pressure (dry air).
  REAL :: cpv     !            idem               (water vapour).
  REAL :: rd      ! gas constant for dry air.
  REAL :: rv      !    idem      for water vapour.
  REAL :: rcpd    ! rcpd=1./cpd.
  REAL :: vtmpc1  ! vtmpc1=rv/rd-1.
  REAL :: vtmpc2  ! vtmpc2=cpv/cpd-1.
  REAL :: rhoh2o  ! density of liquid water.
  REAL :: alv     ! latent heat for vaporisation.
  REAL :: als     ! latent heat for sublimation.
  REAL :: alf     ! latent heat for fusion.
  REAL :: clw     ! specific heat for liquid water.
  REAL :: tmelt   ! temperature of fusion of ice.
  REAL :: solc    ! solar constant.
  REAL :: stbo    ! stephan boltzmann constant.
  REAL :: yearl   ! length of the year (in days).

  ! constants used for computation of saturation mixing ratio
  !   over liquid water(*c_les*) or ice(*c_ies*).

  REAL :: c1es    ! 610.78
  REAL :: c2es    ! 1es*rd/rv
  REAL :: c3les   ! 17.269
  REAL :: c3ies   ! 21.875
  REAL :: c4les   ! 35.86
  REAL :: c4ies   !  7.66
  REAL :: c5les   ! c3les*(tmelt-c4les)
  REAL :: c5ies   ! c3ies*(tmelt-c4ies)
  REAL :: c5alvcp ! c5les*alv/cpd
  REAL :: c5alscp ! c5ies*als/cpd
  REAL :: alvdcp  ! alv/cpd
  REAL :: alsdcp  ! als/cpd

CONTAINS

  SUBROUTINE inicon

    ! Description:
    ! Preset constants in mo_constants.
    !
    ! Method:
    !
    ! *inicon* is called from *setdyn*.
    !
    ! Authors:
    !
    ! M. Jarraud, ECMWF, December 1982, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! H.-S. Bauer, MPI, Jul 1998, changed
    ! A. Rhodin, MPI, Jan 1999, subroutine inicon put into module mo_constants
    !
    ! for more details see file AUTHORS
    !

    USE mo_start_dataset,   ONLY: ly365

    IMPLICIT NONE

    !  Intrinsic functions 
    INTRINSIC ASIN

    !  Executable statements 

    !-- 1. Preset constants

    api    = 2.*ASIN(1.)
    a      = 6371000.
    omega  = .7292E-4
    g      = 9.80665
    cpd    = 1005.46
    cpv    = 1869.46
    rd     = 287.05
    rv     = 461.51

    rcpd   = 1./cpd
    vtmpc1 = rv/rd - 1.
    vtmpc2 = cpv/cpd - 1.

    rhoh2o = 1000.
    alv    = 2.5008E6
    als    = 2.8345E6
    alf    = als - alv

    clw    = 4186.84
    tmelt  = 273.16

    solc   = 1365.
    stbo   = 5.67E-8

    IF (ly365) THEN
       yearl = 365.2422
    ELSE
       yearl = 360.
    ENDIF

    c1es    = 610.78
    c2es    = c1es*rd/rv
    c3les   = 17.269
    c3ies   = 21.875
    c4les   = 35.86
    c4ies   =  7.66
    c5les   = c3les*(tmelt-c4les)
    c5ies   = c3ies*(tmelt-c4ies)
    c5alvcp = c5les*alv/cpd
    c5alscp = c5ies*als/cpd
    alvdcp  = alv/cpd
    alsdcp  = als/cpd

  END SUBROUTINE inicon

END MODULE mo_constants
