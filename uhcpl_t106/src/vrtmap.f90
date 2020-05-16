!+ map indices of an artificial evenly spaced (in log) vertical grid to
!  the indices of the log of the model vertical grid.
!+ $Id: vrtmap.f90,v 1.6 1999/07/19 09:38:49 m214030 Exp $

SUBROUTINE vrtmap(pkdim,pmap,sigln,dsigln,kdpmap)

  ! Description:
  !
  ! Map indices of an artificial evenly spaced (in log) vertical grid to
  ! the indices of the log of the model vertical grid
  !
  ! Method:
  !
  ! The resultant array of mapped indices will be used by "kdpfnd"
  ! to find the vertical location of any departure point relative
  ! to the model grid.
  !
  ! Authors:
  !
  ! Original version:  J. Olson
  ! Standardized:      J. Rosinski, June 1992
  ! Reviewed:          D. Williamson, P. Rasch, August 1992
  ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_exception, ONLY: finish
  USE mo_doctor,    ONLY: nout

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: pkdim, pmap

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: dsigln(pkdim)    ! intervals between model levels (log)
  REAL, INTENT (IN) :: sigln(pkdim)     ! model levels (log(eta))

  !  Array arguments with intent(Out):
  INTEGER, INTENT (OUT) :: kdpmap(pmap) ! array of mapped indices

  !  Local scalars: 
  REAL :: del                ! artificial grid interval
  REAL :: dp                 ! artificial departure point
  REAL :: eps                ! epsilon factor
  INTEGER :: imin, k, kk
  INTEGER ::  newmap         ! estimated value of "pmap"

  !  External functions 
  INTEGER, EXTERNAL :: ismin

  !  Intrinsic functions 
  INTRINSIC REAL


  !  Executable statements 

  eps = 1.E-05
  del = (sigln(pkdim)-sigln(1))/REAL(pmap)
  imin = ismin(pkdim-1,dsigln,1)
  IF (del+eps>=dsigln(imin)) THEN
    newmap = (sigln(pkdim)-sigln(1))/dsigln(imin) + 1
    WRITE(nout,'(A)') ' VRTMAP:  Not enough artificial grid intervals.'
    WRITE(nout,'(A,I20)') ' Currently, "pmap" is set to ', pmap
    WRITE(nout,'(A,I20)') ' Reset parameter "pmap" to at least ', newmap
    CALL finish('vrtmap','Run terminated.')
  END IF

  kdpmap(1) = 1
!CDIR LOOPCNT=20000
  DO kk = 2, pmap
    dp = sigln(1) + REAL(kk-1)*del
    DO k = 1, pkdim - 1
      IF (dp>sigln(k)+eps) THEN
        kdpmap(kk) = k
      END IF
    END DO
  END DO

  RETURN
END SUBROUTINE vrtmap
