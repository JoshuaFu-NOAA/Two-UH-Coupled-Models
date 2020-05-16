MODULE mo_physc1

  IMPLICIT NONE

  ! ----------------------------------------------------------------
  !
  ! module *mo_physc1* constants to communicate between the main program
  !                    and the radiation subroutines.
  !
  ! ----------------------------------------------------------------

  REAL :: cdisse       !   *solar constant normalised by its annual mean.
  REAL :: crae         !   *ratio:atmospheric height/radius of the earth.
  REAL :: czen1        !   *orbital parameters.
  REAL :: czen2
  REAL :: czen3
  REAL :: czen1m
  REAL :: czen2m
  REAL :: czen3m
  REAL :: cdissem
  REAL :: cosrad(128)  !   *cos of longitudes used in *solange*.
  REAL :: sinrad(128)  !   *sin of longitudes used in *solange*.

END MODULE mo_physc1
