MODULE mo_rad2

  USE mo_parameters

  IMPLICIT NONE

  ! ----------------------------------------------------------------
  !
  ! module *mo_rad2* constants internal to the radiation subroutines.
  !
  ! ----------------------------------------------------------------

  REAL :: caeops           !  *caeop_*  *constants used for aerosol computations.
  REAL :: caeopl
  REAL :: caeopu
  REAL :: caeopd
  REAL, ALLOCATABLE :: cvdaes(:)  !  *c___s*   *refers to *sea aerosols.
  REAL, ALLOCATABLE :: cvdael(:)  !  *c___l*   *refers to *land aerosols.
  REAL, ALLOCATABLE :: cvdaeu(:)  !  *c___u*   *refers to *urban aerosols.
  REAL, ALLOCATABLE :: cvdaed(:)  !  *c___d*   *refers to *desert aerosols.
                           !  *cvdae_*  *constants used for aerosol computations.
  REAL :: ctrbga           !  *c__bga*  *constants used for aerosol computations.
  REAL :: cvobga
  REAL :: cstbga
  REAL :: ctrpt            !  *ctrpt*   *constants used for aerosol computations.
  REAL :: caeadk(3)        !  *caeadk*  *constants used for aerosol computations.
  REAL :: caeadm           !  *caeadm*  *constants used for aerosol computations.

END MODULE mo_rad2
