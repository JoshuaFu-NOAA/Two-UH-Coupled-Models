MODULE mo_diagnostics_zonal

  USE mo_parameters

  IMPLICIT NONE

  !
  ! module *mo_diagnostics_zonal* - zonal means of quantities related
  !                                 to global physical diagnostics
  !

  REAL :: dadconz(jpgl)
  REAL :: dsrad0z(jpgl)
  REAL :: dtrad0z(jpgl)
  REAL :: dsradsz(jpgl)
  REAL :: dtradsz(jpgl)
  REAL :: dvdisz(jpgl)
  REAL :: dhfsz(jpgl)
  REAL :: devapz(jpgl)
  REAL :: dgwdisz(jpgl)
  REAL :: dcvgrz(jpgl)
  REAL :: dcverz(jpgl)
  REAL :: dcvgsz(jpgl)
  REAL :: dcvfrz(jpgl)
  REAL :: dcvmsz(jpgl)
  REAL :: dlsgrz(jpgl)
  REAL :: dlsgsz(jpgl)
  REAL :: dlsmsz(jpgl)
  REAL :: dlserz(jpgl)
  REAL :: dlsesz(jpgl)
  REAL :: dssradz (jpgl)
  REAL :: dstradz(jpgl)
  REAL :: dshflz(jpgl)
  REAL :: dsdtflz(jpgl)
  REAL :: dslsrz(jpgl)
  REAL :: dslssz(jpgl)
  REAL :: dscvrz(jpgl)
  REAL :: dscvsz(jpgl)
  REAL :: dsevwz(jpgl)
  REAL :: dseviz(jpgl)
  REAL :: dssnmtz(jpgl)
  REAL :: dsrosz(jpgl)

END MODULE mo_diagnostics_zonal
