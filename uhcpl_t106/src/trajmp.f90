!+ estimate mid-point of parcel trajectory (geodesic coordinates) 
!+ $Id: trajmp.f90,v 1.5 1999/08/27 17:00:36 m214030 Exp $

SUBROUTINE trajmp(dt,upr,vpr,phipr,lampr)

  ! Description:
  !
  ! Estimate mid-point of parcel trajectory (geodesic coordinates) based
  ! upon horizontal wind field.
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
  REAL :: dt

  !  Array arguments with intent(In):
  REAL ::upr(pgls), vpr(pgls)

  !  Array arguments with intent(InOut):
  REAL, INTENT (INOUT) ::  phipr(pgls)

  !  Array arguments with intent(Out):
  REAL, INTENT (OUT) :: lampr(pgls)

  ! dt      Time interval that corresponds to the parcel trajectory.
  ! upr     u-coordinate of velocity corresponding to the most recent
  !         estimate of the trajectory mid-point (in geodesic system).
  ! vpr     v-coordinate of velocity corresponding to the most recent
  !         estimate of the trajectory mid-point (in geodesic system).
  ! phipr   Phi value at trajectory mid-point (geodesic coordinates).
  !         On entry this is the most recent estimate.
  ! lampr   Lambda value at trajectory mid-point (geodesic coordinates).

  !  Local scalars: 
  INTEGER :: i

  !  Intrinsic functions 
  INTRINSIC COS


  !  Executable statements 

  DO i = 1, pgls
    lampr(i) = -.5*dt*upr(i)/COS(phipr(i))
    phipr(i) = -.5*dt*vpr(i)
  END DO

  RETURN
END SUBROUTINE trajmp
