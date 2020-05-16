!+ compute departure points for semi-Lagrangian transport on surface of sphere
!+ $Id: sphdep.f90,v 1.12 1999/11/15 08:46:54 m214003 Exp $

SUBROUTINE sphdep(ih,jcen,jgc,dt,ra,iterdp,locgeo,uvb,uxl,uxr,lam,phib,lbasiy, &
&      lammp,phimp,lamdp,phidp,idp,jdp,nxpt_a)

  ! Description:
  !
  ! Compute departure points for semi-Lagrangian transport on surface of sphere
  !
  ! Method:
  !
  ! Compute departure points for semi-Lagrangian transport on surface of
  ! sphere using midpoint quadrature.  Computations are done in:
  !
  !   1) "local geodesic"   coordinates for "locgeo" = .true.
  !   2) "global spherical" coordinates for "locgeo" = .false.
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

  USE mo_grid,          ONLY: plond, platd, plev, plon, i1
  USE mo_slt,           ONLY: ppdy
  USE mo_decomposition, ONLY: dc=>local_decomposition

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  REAL, INTENT (IN) :: dt, ra
  INTEGER, INTENT (IN) :: iterdp, ih, jcen, jgc
  LOGICAL, INTENT (IN) :: locgeo

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: lam(plond), lbasiy(4,2,platd), phib(platd), &
&      uvb(plond,plev,2,platd), uxl(plond,plev,2,platd), &
&      uxr(plond,plev,2,platd)
  INTEGER, INTENT (IN) :: nxpt_a(plev)

  !  Array arguments with intent(InOut):
  REAL, INTENT (INOUT) :: lammp(plon,plev), phimp(plon,plev)

  !  Array arguments with intent(Out):
  REAL, INTENT (OUT) :: lamdp(plon,plev), phidp(plon,plev)
  INTEGER, INTENT (OUT) :: idp(plon,plev), jdp(plon,plev)


  ! jcen    Index in extended grid corresponding to latitude being
  !         forecast.
  ! jgc     Index in model    grid corresponding to latitude being
  !         forecast.
  ! dt      Time interval that parameterizes the parcel trajectory.
  ! ra      Reciprocal of radius of earth.
  ! iterdp  Number of iterations used for departure point calculation.
  ! locgeo  Logical flag to indicate computation in "local geodesic" or
  !         "global spherical" space.
  ! ub      Longitudinal velocity components in spherical coordinates.
  ! uxl     uxl ( , , 1, ) = x-derivatives of u at the left  (west) edge
  !         of given interval
  !         uxl ( , , 2, ) = x-derivatives of v at the left  (west) edge
  !         of given interval
  ! uxr     uxr ( , , 1, ) = x-derivatives of u at the right (east) edge
  !         of given interval
  !         uxr ( , , 2, ) = x-derivatives of v at the right (east) edge
  !         of given interval
  ! lam     Longitude values for the extended grid.
  ! phib    Latitude  values for the extended grid.
  ! lbasiy  Weights for Lagrange cubic interpolation on the unequally
  !         spaced latitude grid.
  ! lammp   Longitude coordinates of the trajectory mid-points of the
  !         parcels that correspond to the global grid points contained
  !         in the latitude slice being forecast.  On entry lammp
  !         is an initial guess.
  ! phimp   Latitude coordinates of the trajectory mid-points of the
  !         parcels that correspond to the global grid points contained
  !         in the latitude slice being forecast.  On entry phimp
  !         is an initial guess.
  ! lamdp   Longitude coordinates of the departure points that correspond
  !         to the global grid points contained in the latitude slice
  !         being forecast.  lamdp is constrained so that
  !         0.0 .le. lamdp(i) .lt. 2*pi .
  ! phidp   Latitude coordinates of the departure points that correspond
  !         to the global grid points contained in the latitude slice
  !         being forecast.  If phidp is computed outside the latitudinal
  !         domain of the extended grid, then an abort will be called by
  !         subroutine "trjgl".
  ! idp     Longitude index of departure points.  This index points into
  !         the extended arrays, e.g.,
  !         lam (idp(i,k)) .le. lamdp(i,k) .lt. lam (idp(i,k)+1).
  ! jdp     Latitude  index of departure points.  This index points into
  !         the extended arrays, e.g.,
  !         phib(jdp(i,k)) .le. phidp(i,k) .lt. phib(jdp(i,k)+1).

  !  Local scalars: 
  REAL :: cphic       ! cos(phicen)
  REAL :: dlam        ! increment of grid in x-direction
  REAL :: dttmp       ! time step (seconds)
  REAL :: finc        ! time step factor
  REAL :: phicen      ! latitude coord of current lat slice
  REAL :: sphic       ! sin(phicen)
  REAL :: pi
  INTEGER :: i, iter, k

  !  Local arrays: 
  REAL :: fint(plon,plev,ppdy,2) ! u/v x-interpolants
  REAL :: lampr(plon,plev)       ! relative long coord of dep pt
  REAL :: phipr(plon,plev)       ! relative lat  coord of dep pt
  REAL :: upr(plon,plev)         ! u in local geodesic coords
  REAL :: uvmp(plon,plev,2)      ! u/v (spherical) interpltd to dep pt
  REAL :: vpr(plon,plev)         ! v in local geodesic coords

  !  External subroutines 
  EXTERNAL bandij, g2spos, hrvint, s2gphi, s2gvel, trajmp, trjgl, trjmps

  !  Intrinsic functions 
  INTRINSIC COS, SIN, ATAN, REAL


  !  Executable statements 

  !  dlam = lam(nxpt+2) - lam(nxpt+1)
  ! for parallel version
  pi = 4.0*ATAN(1.0)
  dlam = 2.0*pi/REAL(dc%nlon)

  phicen = phib(jcen)
  cphic = COS(phicen)
  sphic = SIN(phicen)

  ! Convert latitude coordinates of trajectory midpoints from spherical
  ! to local geodesic basis.

  IF (locgeo) CALL s2gphi(lam(i1),cphic,sphic,lammp,phimp,phipr)

  ! Loop over departure point iterates.

  DO iter = 1, iterdp

    ! Compute midpoint indices.

    CALL bandij(dlam,phib,lammp,phimp,idp,jdp,nxpt_a)

    ! Interpolate velocity fields to midpoint locations (tensor product
    ! Lagrange cubic interpolation).

    CALL hrvint(2,uvb,uxl,uxr,lam,lbasiy,lammp,phimp,idp,jdp,fint,uvmp)

    ! Put u/v on unit sphere

    DO k = 1, plev
      DO i = 1, plon
        uvmp(i,k,1) = uvmp(i,k,1)*ra
        uvmp(i,k,2) = uvmp(i,k,2)*ra
      END DO
    END DO

    ! For local geodesic:

    !   a) Convert velocity coordinates at trajectory midpoints from
    !      spherical coordinates to local geodesic coordinates,
    !   b) Estimate midpoint parcel trajectory,
    !   c) Convert back to spherical coordinates

    ! Else, for global spherical

    !      Estimate midpoint trajectory with no conversions

    IF (locgeo) THEN
      CALL s2gvel(uvmp(1,1,1),uvmp(1,1,2),lam(i1),cphic,sphic,lammp,phimp, &
&          upr,vpr)
      CALL trajmp(dt,upr,vpr,phipr,lampr)
      dttmp = 0.5*dt

      CALL g2spos(ih,dttmp,lam(i1),phib,phicen,cphic,sphic,upr,vpr, &
&                 lampr,phipr, lammp,phimp)

    ELSE

      CALL trjmps(dt,uvmp(1,1,1),uvmp(1,1,2),phimp,lampr,phipr)
      finc = 1.

      CALL trjgl(jgc,finc,phicen,lam(i1),phib,lampr,phipr,lammp,phimp)
    END IF

  END DO ! End of iter=1,iterdp loop

  ! Compute departure points in geodesic coordinates, and convert back
  ! to spherical coordinates.

  ! Else, compute departure points directly in spherical coordinates

  IF (locgeo) THEN
    DO k = 1, plev
      DO i = 1, plon
        lampr(i,k) = 2.*lampr(i,k)
        phipr(i,k) = 2.*phipr(i,k)
      END DO
    END DO
    dttmp = dt

    CALL g2spos(ih,dttmp,lam(i1),phib,phicen,cphic,sphic,upr,vpr, &
&               lampr,phipr,lamdp,phidp)
  ELSE
    finc = 2.
    CALL trjgl(jgc,finc,phicen,lam(i1),phib,lampr,phipr,lamdp,phidp)
  END IF

  ! Compute departure point indicies.

  CALL bandij(dlam,phib,lamdp,phidp,idp,jdp,nxpt_a)

  RETURN
END SUBROUTINE sphdep
