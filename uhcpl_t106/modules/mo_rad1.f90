MODULE mo_rad1

  IMPLICIT NONE

  ! ----------------------------------------------------------------
  !
  ! module *mo_rad1* spectral distribution of aerosols and ozone.
  !                  (triangular *t10* truncationfor aerosols).
  !                  (triangular *t5* truncationfor ozone).
  !
  ! ----------------------------------------------------------------

  REAL :: caesc(66)  !   *caes_*   *refers to *sea aerosols.
  REAL :: caess(55)
  REAL :: caelc(66)  !   *cael_*   *refers to *land aerosols.
  REAL :: caels(55)
  REAL :: caeuc(66)  !   *caeu_*   *refers to *urban aerosols.
  REAL :: caeus(55)
  REAL :: caedc(66)  !   *caed_*   *refers to *desert aerosols.
  REAL :: caeds(55)
  REAL :: cozqc(21)  !   *coz__*   *refers to *ozone.
  REAL :: cozqs(15)
  REAL :: cozhc(21)
  REAL :: cozhs(15)
                     !   *c___c*   *refers to *cos component.
                     !   *c___s*   *refers to *sin component.

END MODULE mo_rad1
