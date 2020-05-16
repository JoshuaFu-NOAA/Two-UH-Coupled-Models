MODULE mo_vegetation

  IMPLICIT NONE

  ! _______________________________________________________________________
  !
  ! module *mo_vegetation* variables for computation of evapotranspiration
  !
  !
  !      Version 1  C.A Blondin  2/12/86  ECMWF
  !
  ! _______________________________________________________________________

  REAL :: cva     !   constant to define the stomatal resistance
  REAL :: cvb     !   constant to define the stomatal resistance
  REAL :: cvc     !   minimum stomatal resistance
  REAL :: cvbc    !   cvb*cvc
  REAL :: cvxpklt !   exp(cvklt)
  REAL :: cvxmklt !   exp(-cvklt)
  REAL :: cvkc    !   cvk*cvc
  REAL :: cvabc   !   (cva+cvbc)/cvc
  REAL :: cvroots !   percentage of roots in the 1st soil layer
  REAL :: cvrootd !   percentage of roots in the 2nd soil layer
  REAL :: cvrootc !   percentage of roots in the 3rd soil layer
  REAL :: cvrad   !   fraction of the net s.w radiation contributing to p.a.r
  REAL :: cvinter !   efficency of interception of precipitation as rain
  REAL :: cvk

END MODULE mo_vegetation
