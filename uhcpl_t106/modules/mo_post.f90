MODULE mo_post

  USE mo_kind, ONLY: cp
  USE mo_parameters

  IMPLICIT NONE

  !
  ! module *mo_post* - control variables for postprocessing
  !
  !      e. kirk     uni hamburg     2-mar-89
  !

  REAL    :: setval(256)

  INTEGER, ALLOCATABLE :: npplev(:,:)  !   selected levels
  INTEGER :: npbits(256)
  INTEGER :: nlevg3(256)
  INTEGER :: nxpbits(jpg3xf)
  INTEGER :: n4xpbits(jpg3xf)
  INTEGER (cp) :: nunitdf
  INTEGER (cp) :: nunitg4x
  INTEGER :: nunitsp              !   *unit number of output spectral fields
  INTEGER :: nunitun              !   *unit number of output uninterpolated fields
                                  !   - both units may be equal -
  INTEGER :: nlalo                !   nlon * ngl
  INTEGER :: nung4x
  INTEGER :: ncdmin
  INTEGER :: ncdmax

  LOGICAL :: laccu(256)
  LOGICAL :: lxaccu(jpg3xf)
  LOGICAL :: l4xaccu(jpg3xf)
  LOGICAL :: lppspe
  LOGICAL :: lppd
  LOGICAL :: lppvo
  LOGICAL :: lppt
  LOGICAL :: lppp
  LOGICAL :: lppq
  LOGICAL :: lppx

  CHARACTER (8) ::  yn(256)         !   array names

END MODULE mo_post
