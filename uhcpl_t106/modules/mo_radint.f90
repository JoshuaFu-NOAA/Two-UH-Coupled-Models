MODULE mo_radint

  IMPLICIT NONE

  !      -----------------------------------------------------------------
  !
  ! *    coefficients and security parameters within *radint*
  !
  !      -----------------------------------------------------------------

  REAL :: csdtsn    !
  REAL :: calbsea   !  open sea albedo
  REAL :: calbsno   !  thick snow albedo
  REAL :: ccardi    !  specific atmospheric content in co2
  REAL :: cemiss    !  surface longwave emissivity
  REAL :: zsnowal   !  snow depth (in equivalent water) for which snow is considered
                    !  as thick
  REAL :: zmetha
  REAL :: znitox
  REAL :: zcfc(16)
  REAL :: rch4
  REAL :: rn2o
  REAL :: rcfc11
  REAL :: rcfc12
  REAL :: zepalb    !  security to avoid zero albedos.
  REAL :: zepclc    !  security to avoid zero or one cloud covers
  REAL :: zeph2o    !  security to avoid water vapour content in a layer
                    !  to be more then the respective value at saturation.
  REAL :: zepsec    !  avoids 0/0 in the diagnostic of total cloud cover.

END MODULE mo_radint
