!+ transform velocity components at departure points
!+ $Id: s2gvel.f90,v 1.5 1999/08/27 17:00:32 m214030 Exp $

SUBROUTINE s2gvel(udp,vdp,lam,cosphi,sinphi,lamdp,phidp,upr,vpr)

  ! Description:
  !
  ! Transform velocity components at departure points.
  !
  ! Method:
  !
  ! Transform velocity components at departure points associated with a
  ! single latitude slice from spherical coordinates to local geodesic
  ! coordinates. (Williamson and Rasch, 1991)
  !
  ! Authors
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
  REAL, INTENT (IN) :: lam(plon), lamdp(plon,plev), phidp(plon,plev),  &
        &              udp(plon,plev), vdp(plon,plev)

  !  Array arguments with intent(Out):
  REAL, INTENT (OUT) :: upr(plon,plev), vpr(plon,plev)

  ! udp    u-component of departure point velocity in spherical coords.
  ! vdp    v-component of departure point velocity in spherical coords.
  ! lam    Longitude of arrival point position (model grid point) in
  !        spherical coordinates.
  ! cosphi Cos of latitude of arrival point positions (model grid pt).
  ! sinphi Sin of latitude of arrival point positions (model grid pt).
  ! lamdp  Longitude of departure point position in spherical
  !        coordinates.
  ! phidp  Latitude  of departure point position in spherical
  !        coordinates.
  ! upr    u-component of departure point velocity in geodesic coords.
  ! vpr    v-component of departure point velocity in geodesic coords.

  !  Local scalars: 
  REAL :: cdlam, clamp, cphid, cphip, dlam, sdlam, slamp, sphid, sphip
  INTEGER :: i, k

  !  Intrinsic functions 
  INTRINSIC ASIN, COS, SIN


  !  Executable statements 

  DO k = 1, plev
    DO i = 1, plon
      dlam = lam(i) - lamdp(i,k)
      sdlam = SIN(dlam)
      cdlam = COS(dlam)
      sphid = SIN(phidp(i,k))
      cphid = COS(phidp(i,k))
      sphip = sphid*cosphi - cphid*sinphi*cdlam
      cphip = COS(ASIN(sphip))
      slamp = -sdlam*cphid/cphip
      clamp = COS(ASIN(slamp))

      vpr(i,k) = (vdp(i,k)*(cphid*cosphi+sphid*sinphi*cdlam)-udp(i,k)*sinphi* &
&          sdlam)/cphip
      upr(i,k) = (udp(i,k)*cdlam+vdp(i,k)*sphid*sdlam+vpr(i,k)*slamp*sphip)/ &
&          clamp
    END DO
  END DO

  RETURN
END SUBROUTINE s2gvel
