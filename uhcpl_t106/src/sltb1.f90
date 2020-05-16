!+ drive the slt algorithm on a given latitude slice
!+ $Id: sltb1.f90,v 1.11 1999/11/15 08:46:53 m214003 Exp $

#ifndef SLDIAG
SUBROUTINE sltb1(pmap,ih,jcen,jgc,dt,ra,iterdp,uvb,uxl,uxr,wb,fb,fxl,fxr,lam, &
&      phib,dphib,sig,sigh,dsig,dsigh,lbasdy,lbasdz,lbassd,lbasiy,kdpmpf, &
&      kdpmph,lammp,phimp,sigmp,fbout,nxpt_a)
#else
SUBROUTINE sltb1(pmap,ih,jcen,jgc,dt,ra,iterdp,uvb,uxl,uxr,wb,fb,fxl,fxr,lam, &
&      phib,dphib,sig,sigh,dsig,dsigh,lbasdy,lbasdz,lbassd,lbasiy,kdpmpf, &
&      kdpmph,lammp,phimp,sigmp,fbout,nxpt_a,diag)
#endif

  ! Description:
  !
  ! Drive the slt algorithm on a given latitude slice.
  !
  ! Method:
  !
  ! Drive the slt algorithm on a given latitude slice in the extended
  ! data arrays using information from the entire latitudinal extent
  ! of the arrays.
  !
  ! Authors:
  !
  ! Original version:  J. Olson
  ! Standardized:      J. Rosinski, June 1992
  ! Reviewed:          D. Williamson, P. Rasch, August 1992
  ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
  ! parallel version:  T. Diehl, DKRZ, July 1999
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_grid, ONLY: platd, plev, plevp1, plond, pcnst, plon, nxpt
  USE mo_slt,  ONLY: ppdy, plimdr

  IMPLICIT NONE

  !  Local parameters: 
  REAL, PARAMETER :: phigs = 1.221730  ! cut-off latitude: about 70 degree

  !  Scalar arguments with intent(In):
  REAL, INTENT (IN) :: dt, ra
  INTEGER, INTENT (IN) :: iterdp, ih, jcen, jgc, pmap

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: dphib(platd), dsig(plev), dsigh(plevp1), &
       &           fb(plond,plev,pcnst,platd), wb(plon,plevp1),  &
       &           fxl(plond,plev,pcnst,platd), fxr(plond,plev,pcnst,platd), &
       &           lam(plond), lbasdy(4,2,platd), lbasdz(4,2,plev), &
       &           lbasiy(4,2,platd), lbassd(4,2,plevp1), phib(platd),  &
       &           sig(plev), sigh(plevp1), uvb(plond,plev,2,platd),  &
       &           uxl(plond,plev,2,platd), uxr(plond,plev,2,platd)
  INTEGER, INTENT (IN) :: kdpmpf(pmap), kdpmph(pmap)

  !  Array arguments with intent(InOut):
  REAL, INTENT (INOUT) :: lammp(plon,plev), phimp(plon,plev), sigmp(plon,plev)

  !  Array arguments with intent(Out):
  REAL, INTENT (OUT) :: fbout(plond,plev,pcnst)
#ifdef SLDIAG
  REAL, INTENT (OUT) :: diag(plond,plev,pcnst)
#endif

  INTEGER :: nxpt_a(plev)

  ! pmap    Dimension of kdpmpX arrays
  ! jcen    Latitude index in extended grid corresponding to lat slice
  !         being forecasted.
  ! jgc     Latitude index in model    grid corresponding to lat slice
  !         being forecasted.
  ! dt      Time interval that parameterizes the parcel trajectory.
  ! ra      Reciprocal of radius of earth.
  ! iterdp  Number of iterations used for departure point calculation.
  ! ub      u-velocity component.
  ! uxl     uxl ( , , 1, ) = x-derivatives of u at the left  (west) edge
  !         of given interval
  ! uxl ( , , 2, ) = x-derivatives of v at the left  (west) edge
  !         of given interval
  ! uxr     uxr ( , , 1, ) = x-derivatives of u at the right (east) edge
  !         of given interval
  ! uxr ( , , 2, ) = x-derivatives of v at the right (east) edge
  !         of given interval
  ! wb      z-velocity component (eta-dot).
  ! fb      Scalar components to be transported.
  ! fxl     x-derivatives at the left  edge of each interval containing
  !         the departure point.
  ! fxr     x-derivatives at the right edge of each interval containing
  !         the departure point.
  ! lam     Longitude values for the extended grid.
  ! phib    Latitude  values for the extended grid.
  ! dphib   Interval between latitudes in the extended grid.
  ! sig     Hybrid eta values at the "full-index" levels.
  ! sigh    Half-index eta-levels including sigh(i,1) = eta(1/2) = 0.0
  !         and sigh(i,plev+1) = eta(plev+1/2) = 1.  Note that in general
  !         sigh(i,k) .lt. sig(i,k)  where sig(i,k) is the hybrid value
  !         at the k_th full-index level.
  ! dsig    Interval lengths in full-index hybrid level grid.
  ! dsigh   Interval lengths in half-index hybrid level grid.
  ! lbasdy  Weights for Lagrange cubic derivative estimates on the
  !         unequally spaced latitude grid.
  ! lbasdz  Weights for Lagrange cubic derivative estimates on the
  !         unequally spaced vertical grid (full levels).
  ! lbassd  Weights for Lagrange cubic derivative estimates on the
  !         unequally spaced vertical grid (half levels).
  ! lbasiy  Weights for Lagrange cubic interpolation on the unequally
  !         spaced latitude grid.
  ! kdpmpf  indices of artificial grid mapped into the full level grid
  ! kdpmph  indices of artificial grid mapped into the half level grid
  ! lammp   Longitude coordinates of the trajectory mid-points of the
  ! parcels that correspond to the global grid points contained
  !         in the latitude slice being forecasted.  On entry lammp
  !         is an initial guess.
  ! phimp   Latitude coordinates of the trajectory mid-points of the
  !         parcels that correspond to the global grid points contained
  !         in the latitude slice being forecasted.  On entry phimp
  !         is an initial guess.
  ! sigmp   Hybrid value at the trajectory midpoint for each gridpoint
  !         in a vertical slice from the global grid.  On entry sigmp is
  !         an initial guess.
  ! fbout   Extended array only one latitude of which, however, is filled
  !         with forecasted (transported) values.  This routine must be
  !         called multiple times to fill the entire array.  This is
  !         done to facilitate multi-tasking.
  ! diag    Time tendency due to horizontal advection for each scalar
  !         component at the global grid points that correspond to the
  !         latitude slice being forecasted.  As in "fbout", only one
  !         latitude slice at a time is filled.

  !  Local scalars: 
  INTEGER :: m
  INTEGER :: i, k
  INTEGER :: isafe, nxpt_max, iext ! for parallel version
  LOGICAL :: locgeo  ! flag indicating coordinate sys
  
  REAL :: fhr(plon,plev,pcnst)        ! horizontal interpolants
  REAL :: fhsb(plon,plev,pcnst)       ! derivative at bot of interval
  REAL :: fhst(plon,plev,pcnst)       ! derivative at top of interval
  REAL :: fint(plon,plev,ppdy,pcnst)  ! work space
  REAL :: fyb(plon,plev,pcnst)        ! work space
  REAL :: fyt(plon,plev,pcnst)        ! work space
  REAL :: lamdp(plon,plev)            ! zonal      departure pt. coord.
  REAL :: phidp(plon,plev)            ! meridional departure pt. coord.
  REAL :: sigdp(plon,plev)            ! vertical   departure pt. coord.
  REAL :: wsb(plon,plevp1)             ! w derivative at bot of interval
  REAL :: wst(plon,plevp1)             ! w derivative at top of interval

  REAL :: fdp(plon,plev,pcnst)

#ifdef SLDIAG
  REAL :: rdt                         ! 1./dt
#endif
  INTEGER :: idp(plon,plev)           ! zonal      dep point index
  INTEGER :: jdp(plon,plev)           ! meridional dep point index
  INTEGER :: kdp(plon,plev)           ! vertical   dep point index

  !  External subroutines 
  EXTERNAL cubzdr, herzin, hrintp, limdz, sphdep, vrtdep

  !  Intrinsic functions 
  INTRINSIC ABS


  !  Executable statements 

  ! Horizontal interpolation

  ! Compute departure points and corresponding indices.
  ! Poleward of latitude phigs (radians), perform the computation in
  ! local geodesic coordinates.
  ! Equatorward of latitude phigs, perform the computation in global
  ! spherical coordinates

  locgeo = ABS(phib(jcen)) >= phigs

  CALL sphdep(ih,jcen,jgc,dt,ra,iterdp,locgeo,uvb,uxl,uxr,lam,phib,lbasiy, &
              lammp,phimp,lamdp,phidp,idp,jdp,nxpt_a)

  ! Interpolate scalar fields to the departure points.

  CALL hrintp(pcnst,fb,fxl,fxr,lam,phib,dphib,lbasdy,lamdp,phidp,idp, &
&      jdp,plimdr,fint,fyb,fyt,fhr)

#ifdef SLDIAG

  ! Compute time tendency due to horizontal advection.

  rdt = 1./dt
  DO m = 1, pcnst
    DO k = 1, plev
      DO i = 1, plon
        diag(i,k,m) = (fhr(i,k,m)-fb(nxpt+i,k,m,jcen))*rdt
      END DO
    END DO
  END DO
#endif

  ! Vertical interpolation.
  ! Compute vertical derivatives of vertical wind

  CALL cubzdr(plon,plevp1,wb,lbassd,wst,wsb)

  ! Compute departure points and corresponding indices.

  CALL vrtdep(pmap,dt,iterdp,wb,wst,wsb,sig,sigh,dsigh,kdpmpf,kdpmph,sigmp, &
&      sigdp,kdp)

  ! Vertical derivatives of scalar fields.
  ! Loop over constituents.

  DO m = 1, pcnst
    CALL cubzdr(plon,plev,fhr(1,1,m),lbasdz,fhst(1,1,m),fhsb(1,1,m))
  END DO
  IF (plimdr) THEN
    CALL limdz(fhr,dsig,fhst,fhsb)
  END IF

  ! Vertical interpolation of scalar fields.

  CALL herzin(plev,pcnst,fhr,fhst,fhsb,sig,dsig,sigdp,kdp,fdp)


  ! Transfer transported values to extended array

  DO m = 1,pcnst
     DO k = 1,plev
        DO i = 1,plon
           fbout(nxpt+i,k,m) = fdp(i,k,m)
        END DO
     END DO
  END DO
  
  ! Locally update the nxpt_a array using idp
  ! The array values over the row of processors in longitude
  ! are synchronized in sltini

  isafe = 4
  DO k = 1,plev
     nxpt_max = 0
     DO i = 1,plon
        iext= i + nxpt
        nxpt_max = MAX(nxpt_max,ABS(iext - idp(i,k)) + isafe)
     END DO
     nxpt_a(k) = MIN(nxpt_max,nxpt)
  END DO
  
  RETURN
END SUBROUTINE sltb1
