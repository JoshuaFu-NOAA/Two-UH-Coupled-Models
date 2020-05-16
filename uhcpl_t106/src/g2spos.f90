!+ transform position coordinates for a set of points
!+ $Id: g2spos.f90,v 1.12 1999/09/27 13:26:02 k202057 Exp $

SUBROUTINE g2spos(ih,dttmp,lam,phib,phi,cosphi,sinphi,upr,vpr,lamgc,phigc,lamsc,phisc)

  ! Description:
  !
  ! Transform position coordinates for a set of points.
  !
  ! Method:
  !
  ! Transform position coordinates for a set of points, each of which is
  ! associated with a grid point in a global latitude slice, from local
  ! geodesic to spherical coordinates.
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

  USE mo_grid,          ONLY: plon, pgls, platd, plev, j1, plato2
  USE mo_decomposition, ONLY: dc=>local_decomposition, debug_seriell

  IMPLICIT NONE

  !  Local parameters: 
  ! 1. - eps, eps from machine precision
  REAL, PARAMETER :: fac = 1. - 1.E-12

  !  Scalar arguments with intent(In):
  REAL, INTENT (IN) :: cosphi, dttmp, phi, sinphi

integer, intent (in) :: ih

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: lam(plon), lamgc(pgls), phib(platd), phigc(pgls), &
&                      upr(pgls), vpr(pgls)

  !  Array arguments with intent(Out):
  REAL, INTENT (OUT) :: lamsc(plon,plev), phisc(pgls)

  ! dttmp  Time step over which midpoint/endpoint trajectory is
  !        calculated (seconds).
  ! lam    Longitude coordinates of the global grid points in spherical
  !        system.  The grid points in the global array are the reference
  ! points for the local geodesic systems.
  ! phib   Latitude values for the extended grid.
  ! phi    Latitude coordinate (in the global grid) of the current
  !        latitude slice.
  ! cosphi cos( phi )
  ! sinphi sin( phi )
  ! upr    zonal      velocity at departure point in local geodesic coord
  ! vpr    Meridional velocity at departure point in local geodesic coord
  ! lamgc  Longitude coordinate of points in geodesic coordinates.
  ! phigc  Latitude coordinate of points in geodesic coordinates.
  ! lamsc  Longitude coordinate of points in spherical coordinates.
  ! phisc  Latitude coordinate of points in spherical coordinates.

  !  Local scalars: 
  REAL :: clamgc        ! cos(lamgc)
  REAL :: coeff         ! tmp variable
  REAL :: cphigc        ! cos(phigc)
  REAL :: distmx        ! max distance
  REAL :: phipi2        ! tmp variable
  REAL :: pi            ! 4.*atan(1.)
  REAL :: pi2           ! pi/2
  REAL :: sgnphi        ! holds sign of phi
  REAL :: slam2         ! sin(lamgc)**2
  REAL :: sphigc        ! sin(phigc)
  REAL :: twopi         ! 2.*pi
  INTEGER :: i, ii, k
  INTEGER :: nval       ! number of values returned from whenfgt

  !  Local arrays: 
  REAL :: dist(pgls)    ! approx. distance traveled along traj.
  REAL :: dlam(pgls)    ! zonal extent of trajectory
  REAL :: slamgc(pgls)  ! sin(lamgc)
  INTEGER :: indx(pgls) ! index holder

  !  External subroutines 
  EXTERNAL whenfgt

  !  Intrinsic functions 
  INTRINSIC ABS, ASIN, ATAN, COS, SIGN, SIN, SQRT


  !  Executable statements 

  pi = 4.*ATAN(1.)
  twopi = pi*2.
  pi2 = pi/2.
  coeff = (1.1*dttmp)**2
  distmx = (SIGN(pi2,phi)-phi)**2/coeff
  sgnphi = SIGN(1.,phi)

  DO i = 1, pgls
    sphigc = SIN(phigc(i))
    cphigc = COS(phigc(i))
    slamgc(i) = SIN(lamgc(i))
    clamgc = COS(lamgc(i))
    phisc(i) = ASIN((sphigc*cosphi+cphigc*sinphi*clamgc)*fac)

!    IF (ABS(phisc(i))>=phib(j1+plat)*fac) phisc(i) = SIGN(phib(j1+plat), &
!&        phisc(i))*fac

!    IF (ABS(phisc(i))>=phib(j1+plato2)*fac) phisc(i) = SIGN(phib(j1+plato2), &
!&        phisc(i))*fac

!         if(phisc(i) .ge. phib(j1+plato2))   phisc(i) = phib(j1+plato2)*fac

    ! parallel version
    IF (dc%set_a ==1) THEN  ! If I am a pole PE
	IF (ih == 2) THEN
   	 IF (phisc(i) >= phib(j1+plato2)*fac) phisc(i) = phib(j1+plato2)*fac
	ELSE
   	 IF (phisc(i) <= phib(j1-1)*fac) phisc(i) = phib(j1-1)*fac
	END IF
    END IF

    dlam(i) = ASIN((slamgc(i)*cphigc/COS(phisc(i)))*fac)

    ! Compute estimated trajectory distance based upon winds alone

    dist(i) = upr(i)**2 + vpr(i)**2
  END DO

  ! Determine which trajectories may have crossed over pole

  CALL whenfgt(pgls,dist,1,distmx,indx,nval)

  ! Check that proper branch of arcsine is used for calculation of
  ! dlam for those trajectories which may have crossed over pole.

!DIR$ IVDEP
!OCL NOVREC
  DO ii = 1, nval
    i = indx(ii)
    slam2 = slamgc(i)**2
    phipi2 = ASIN((SQRT((slam2-1.)/(slam2-1./cosphi**2)))*fac)
    IF (sgnphi*phigc(i) > phipi2) THEN
      dlam(i) = SIGN(pi,lamgc(i)) - dlam(i)
    END IF
  END DO

  DO k = 1, plev
    DO i = 1, plon
      lamsc(i,k) = lam(i) + dlam((k-1)*plon+i)

      ! Restrict longitude to be in the range [0, twopi).
      ! But: Parallel code allows extension beyond twopi

      IF (debug_seriell .and. dc%nprocb == 1) THEN     
         IF (lamsc(i,k) >= twopi) lamsc(i,k) = lamsc(i,k) - twopi
         IF (lamsc(i,k) < 0.0)    lamsc(i,k) = lamsc(i,k) + twopi
      END IF

    END DO
  END DO

  RETURN
END SUBROUTINE g2spos
