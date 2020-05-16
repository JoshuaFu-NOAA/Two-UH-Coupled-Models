MODULE mo_stat_zonal

  USE mo_parameters

  IMPLICIT NONE

  ! module *mo_stat_zonal* - zonal means for global statistics (dynamics)

  REAL, ALLOCATABLE :: delph(:)

  REAL :: voz(jpgl)
  REAL :: dz(jpgl)
  REAL :: qz(jpgl)
  REAL :: xz(jpgl)
  REAL :: tz(jpgl)
  REAL :: psz(jpgl)

  REAL :: gkez(jpgl)
  REAL :: gpez(jpgl)
  REAL :: glqz(jpgl)
  REAL :: gpsz(jpgl)
  REAL :: gtsz(jpgl)
  REAL :: gtdz(jpgl)
  REAL :: gwsz(jpgl)
  REAL :: gwdz(jpgl)
  REAL :: gsnz(jpgl)

END MODULE mo_stat_zonal
