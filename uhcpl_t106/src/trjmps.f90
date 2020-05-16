!+ estimate mid-point interval of parcel trajectory
!+ $Id: trjmps.f90,v 1.5 1999/08/27 17:00:37 m214030 Exp $

SUBROUTINE trjmps(dt,upr,vpr,phimp,lampr,phipr)

  ! Description:
  !
  ! Estimate mid-point interval of parcel trajectory (global spherical
  ! coordinates).
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

  USE mo_grid, ONLY: pgls

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  REAL, INTENT (IN) :: dt

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: phimp(pgls), upr(pgls), vpr(pgls)

  !  Array arguments with intent(Out):
  REAL, INTENT (OUT) :: lampr(pgls), phipr(pgls)

  ! dt      Time interval that corresponds to the parcel trajectory.
  ! upr     u-coordinate of velocity corresponding to the most recent
  !         estimate of the trajectory mid-point.
  ! vpr     v-coordinate of velocity corresponding to the most recent
  !         estimate of the trajectory mid-point.
  ! phimp   Phi value of trajectory midpoint (most recent estimate).
  ! lampr   Longitude coordinate of trajectory mid-point relative to the
  !         arrival point.
  ! phipr   Latitude  coordinate of trajectory mid-point relative to the
  !         arrival point.

  !  Local scalars: 
  INTEGER :: i

  !  Intrinsic functions 
  INTRINSIC COS


  !  Executable statements 

  DO i = 1, pgls
    lampr(i) = -.5*dt*upr(i)/COS(phimp(i))
    phipr(i) = -.5*dt*vpr(i)
  END DO

  RETURN
END SUBROUTINE trjmps
