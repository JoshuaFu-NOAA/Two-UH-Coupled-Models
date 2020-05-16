!+ initializes module mo_diff
!+ $Id: sudif.f90,v 1.7 1998/11/10 15:14:48 m214003 Exp $

SUBROUTINE sudif

  ! Description:
  !
  ! Initializes module *mo_diff* for horizontal diffusion subroutines
  !
  ! Authors:
  !
  ! M. Esch, MPI, September 1993, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_parameters
  USE mo_control
  USE mo_diff

  IMPLICIT NONE

  !  Local scalars: 
  INTEGER :: jn

  !  Local arrays: 
  INTEGER :: ihq106(nlev), ihq21(nlev), ihq30(nlev), ihq42(nlev)

  !  Executable statements 

!-- 1. Set parameter
  
  IF (lmidatm) THEN
    ihq106(:)  = 4

    ihq42(:)   = 5

    ihq30(:)   = 5

    ihq21(:)   = 2
  ELSE
    ihq106(1:2) = 1
    ihq106(3:4) = 2
    ihq106(5)   = 3
    ihq106(6:)  = 4

    ihq42(1:2)  = 1
    ihq42(3:4)  = 2
    ihq42(5)    = 3
    ihq42(6)    = 4
    ihq42(7:)   = 5

    ihq30(1:2)  = 1
    ihq30(3:4)  = 2
    ihq30(5)    = 3
    ihq30(6)    = 4
    ihq30(7:)   = 5

    ihq21(1)    = 1
    ihq21(2)    = 2
    ihq21(3)    = 3
    ihq21(4)    = 4
    ihq21(5)    = 6
    ihq21(6)    = 8
    ihq21(7:)   = 10
  ENDIF

  ncdif(1:nlev) = 0

!-- 2. Copy to iq

  DO jn = 1, nlev
    IF (nn==21) THEN
      iq(jn) = ihq21(jn)
    ELSE IF (nn==30) THEN
      iq(jn) = ihq30(jn)
    ELSE IF (nn==106) THEN
      iq(jn) = ihq106(jn)
    ELSE
      iq(jn) = ihq42(jn)
    END IF
  END DO

  RETURN
END SUBROUTINE sudif
