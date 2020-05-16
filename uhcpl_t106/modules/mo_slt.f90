MODULE mo_slt

  IMPLICIT NONE

  !
  !  parameters common to many slt routines
  !
  INTEGER, PARAMETER :: ppdy   = 4      ! length of interpolation grid stencil
  LOGICAL, PARAMETER :: plimdr = .TRUE. ! flag to limit derivatives

END MODULE mo_slt
