!+ map relative trajectory mid/departure point coordinates to global
!  latitude/longitude coordinates and test limits
!+ $Id: trjgl.f90,v 1.15 2002/12/19 09:47:51 m214003 Exp $

SUBROUTINE trjgl(jgc,finc,phicen,lam,phib,lampr,phipr,lamp,phip)

  ! Description:
  !
  ! Map relative trajectory mid/departure point coordinates to global
  ! latitude/longitude coordinates and test limits
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

  USE mo_exception,     ONLY: finish
  USE mo_doctor,        ONLY: nerr
  USE mo_grid,          ONLY: plon, plev, platd, pgls
  USE mo_decomposition, ONLY: dc=>local_decomposition, debug_seriell

  IMPLICIT NONE

  !  Scalar arguments with intent(in):
  REAL, INTENT (in) :: finc, phicen
  INTEGER, INTENT (in) :: jgc

  !  Array arguments with intent(in):
  REAL, INTENT (in) :: lam(plon), lampr(plon,plev), phib(platd), phipr(plon,plev)

  !  Array arguments with intent(Out):
  REAL, INTENT (out) :: lamp(plon,plev), phip(plon,plev)


  ! jgc     Index in model grid corresponding to latitude being forecast.
  ! finc    Time step factor (1. for midpoint, 2. for dep. point)
  ! phicen  Latitude value for current latitude being forecast.
  ! lam     Longitude values for the extended grid.
  ! phib    Latitude  values for the extended grid.
  ! lampr   Longitude coordinates (relative to the arrival point) of the
  !         trajectory mid-points of the parcels that correspond to the
  ! global grid points contained in the latitude slice being
  !         forecast.
  ! phipr   Latitude coordinates (relative to the arrival point) of the
  !         trajectory mid-points of the parcels that correspond to the
  !         global grid points contained in the latitude slice being
  !         forecast.
  ! lamp    Longitude coordinates of the trajectory mid-points of the
  !         parcels that correspond to the global grid points contained
  !         in the latitude slice being forecast.
  ! phip    Latitude  coordinates of the trajectory mid-points of the
  !         parcels that correspond to the global grid points contained
  !         in the latitude slice being forecast.

  !  Local scalars: 
  REAL :: pi, twopi
  INTEGER :: i, iimax, iimin, imax, imin, k, kkmax, kkmin

  !  External functions 
  INTEGER, EXTERNAL :: ismax, ismin

  !  Intrinsic functions 
  INTRINSIC ATAN


  !  Executable statements 

  pi = 4.*ATAN(1.)
  twopi = pi*2.

  DO k = 1, plev
    DO i = 1, plon
      lamp(i,k) = lam(i) + finc*lampr(i,k)
      phip(i,k) = phicen + finc*phipr(i,k)

      IF (debug_seriell .AND. dc%nprocb == 1) THEN
         IF (lamp(i,k) >= twopi) lamp(i,k) = lamp(i,k) - twopi
         IF (lamp(i,k) < 0.0) lamp(i,k) = lamp(i,k) + twopi
      END IF

    END DO
  END DO

  ! Test that the latitudinal extent of trajectory is NOT over the poles

  imax = ismax(pgls,phip,1)
  imin = ismin(pgls,phip,1)
  kkmax = (imax-1)/plon + 1
  iimax = imax - (kkmax-1)*plon
  kkmin = (imin-1)/plon + 1
  iimin = imin - (kkmin-1)*plon

  ! Since trjgl is called only if ABS(lat) < 70(degrees):
  ! Test that the latitudinal extent of trajectory is NOT over the poles
  !
  IF (phip(iimax,kkmax) >= pi/2.0)  THEN
    CALL errmsg
  ELSE IF (phip(iimin,kkmin) <= -pi/2.0) THEN
    CALL errmsg
  END IF

  RETURN

CONTAINS

  SUBROUTINE errmsg

    WRITE (nerr,'(/,A)') ' Model is blowing up ...'
    WRITE (nerr,'(/,A,I5,A,I5,A,I5,A)') &
&          ' Parcel associated with longitude ',iimax,', level ', kkmax, &
&          ' and latitude ',jgc,' is outside the model domain.'
    CALL finish('trjgl','Run terminated.')

  END SUBROUTINE errmsg

END SUBROUTINE trjgl
