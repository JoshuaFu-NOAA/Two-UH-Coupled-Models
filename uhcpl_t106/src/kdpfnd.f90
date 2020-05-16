!+ determine vertical departure point indices
!+ $Id: kdpfnd.f90,v 1.5 1999/08/27 17:00:24 m214030 Exp $

SUBROUTINE kdpfnd(pkdim,pmap,sig,sigdp,kdpmap,kdp)

  ! Description:
  !
  ! Determine vertical departure point indices.
  !
  ! Method:
  !
  ! Determine vertical departure point indices that point into a grid
  ! containing the full or half sigma levels.  Use an artificial evenly
  ! spaced vertical grid to map into the true model levels.
  !
  ! Indices are computed assuming the the sigdp values have
  ! been constrained so that sig(1) .le. sigdp(i,j) .lt. sig(pkdim).
  !
  ! Authors:
  !
  ! Standardized:      J. Rosinski, June 1992
  ! Reviewed:          D. Williamson, P. Rasch, August 1992
  ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_grid, ONLY: pgls

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: pkdim         ! dimension of "sig"
  INTEGER, INTENT (IN) :: pmap          ! dimension of "kdpmap"

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: sig(pkdim)       ! vertical grid coordinates
  REAL, INTENT (IN) ::  sigdp(pgls)     ! vertical coords. of departure points
  INTEGER, INTENT (IN) ::  kdpmap(pmap) ! array of model grid indices which
                                        ! are mapped into the artificial grid.

  !  Array arguments with intent(Out):
  INTEGER, INTENT (OUT) :: kdp(pgls)    ! vertical index for each dep. pt.

  !  Local scalars: 
  REAL :: rdel          ! Reciprocal of interval in artificial g
  REAL :: sig1ln        ! ln (sig(1))
  INTEGER :: i, ii      ! indices

  !  Intrinsic functions 
  INTRINSIC INT, LOG, MAX, MIN, REAL


  !  Executable statements 
  rdel = REAL(pmap)/(LOG(sig(pkdim))-LOG(sig(1)))
  sig1ln = LOG(sig(1))

  DO i = 1, pgls

    ! First guess of the departure point's location in the model grid

    ii = MAX(1,MIN(pmap,INT((LOG(sigdp(i))-sig1ln)*rdel+1.)))
    kdp(i) = kdpmap(ii)

    ! Determine if location is in next interval

    IF (sigdp(i)>=sig(kdp(i)+1)) THEN
      kdp(i) = kdp(i) + 1
    END IF
  END DO

  RETURN
END SUBROUTINE kdpfnd
