MODULE mo_radiation

  IMPLICIT NONE

  !      -----------------------------------------------------------------
  !
  ! *    control variables for radiation: internal switch and indices
  !
  !      -----------------------------------------------------------------

  INTEGER :: nuaer
  INTEGER :: ntraer

  REAL    :: diff
  REAL    :: c10e
  REAL    :: zeelog    !   security threshold for absorber amount in laplace transform
  REAL    :: zepsc     !   security threshold for cloud cover
  REAL    :: zepsco    !   security threshold for ozone amount
  REAL    :: zepscq    !   security threshold for water vapor
  REAL    :: zepsct    !   security threshold for shortwave optical thickness
  REAL    :: zepscw    !   security threshold for cloud liquid water path

      DATA NUAER, NTRAER /29, 15/
      DATA ZEPSC, ZEPSCO, ZEPSCQ  /1.E-20, 1.E-10, 1.E-10/
      DATA ZEPSCT, ZEPSCW, ZEELOG /1.E-20, 1.E-20, 1.E-10/
      DATA DIFF, C10E /1.66, 0.4342945/

END MODULE mo_radiation
