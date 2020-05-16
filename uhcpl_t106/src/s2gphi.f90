!+ calculate transformed local geodesic latitude coordinate
!+ $Id: s2gphi.f90,v 1.5 1999/08/27 17:00:31 m214030 Exp $

SUBROUTINE s2gphi(lam,cosphi,sinphi,lamsc,phisc,phigc)

  ! Description:
  !
  ! Calculate transformed local geodesic latitude coordinates.
  !
  ! Method:
  !
  ! Calculate transformed local geodesic latitude coordinates for a set
  ! of points, each of which is associated with a grid point in a global
  ! latitude slice. Transformation is spherical to local geodesic.
  ! (Williamson and Rasch, 1991)
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

  USE mo_grid, ONLY: plon, plev

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  REAL, INTENT (IN) :: cosphi, sinphi

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: lam(plon), lamsc(plon,plev), phisc(plon,plev)

  !  Array arguments with intent(Out):
  REAL, INTENT (OUT) :: phigc(plon,plev) 

  ! lam    longitude coordinates of the global grid points in spherical
  !        system. The grid points in the global array are the reference
  !        points for the local geodesic systems.
  ! cosphi cosine of the latitude of the global latitude slice.
  ! sinphi sine of the latitude of the global latitude slice.
  ! lamsc  longitude coordinate of dep. points in spherical coordinates.
  ! phisc  latitude  coordinate of dep. points in spherical coordinates.
  ! phigc  latitude  coordinate of dep. points in local geodesic coords.

  !  Local scalars: 
  REAL :: clamsc, cphisc, sphisc
  INTEGER :: i, k

  !  Intrinsic functions 
  INTRINSIC ASIN, COS, SIN


  !  Executable statements 

  DO k = 1, plev
    DO i = 1, plon
      sphisc = SIN(phisc(i,k))
      cphisc = COS(phisc(i,k))
      clamsc = COS(lam(i)-lamsc(i,k))
      phigc(i,k) = ASIN(sphisc*cosphi-cphisc*sinphi*clamsc)
    END DO
  END DO

  RETURN
END SUBROUTINE s2gphi
