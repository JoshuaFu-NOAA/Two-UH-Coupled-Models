!+ Initialize model and extended grid parameters
!+ $Id: grdini.f90,v 1.8 1999/09/24 08:10:12 k202057 Exp $

SUBROUTINE grdini(pmap,etamid,etaint,gravit,dlam,lam,phi,dphi,gw,sinlam, &
&      coslam,lbasdy,lbasdz,lbassd,lbasiy,detam,detai,kdpmpf,kdpmph,cwava)

  !
  ! Description:
  !
  ! Initialize model and extended grid parameters.
  ! Initialize weights for Lagrange cubic derivative estimates.
  ! Initialize weights for Lagrange cubic interpolant.
  !
  ! Method:
  !
  ! Authors:
  !
  ! Original version:  J. Olson
  ! Standardized:      J. Rosinski, June 1992
  ! Reviewed:          D. Williamson, P. Rasch, August 1992
  ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
  ! Parallel version:  T. Diehl, DKRZ, July 1999
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_decomposition, ONLY: lc => local_decomposition
  USE mo_grid, ONLY: plat, platd, plev, plevp1, plon, plond, dphibr

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  REAL, INTENT (IN) :: gravit
  INTEGER, INTENT (IN) :: pmap

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: etamid(plev), etaint(plevp1)

  !  Scalar arguments with intent(Out):
  REAL, INTENT (OUT) :: cwava, dlam

  !  Array arguments with intent(Out):
  REAL, INTENT (OUT) :: coslam(plon), detai(plevp1), detam(plev), &
&      dphi(platd,2), &
&      gw(plat), lam(plond), lbasdy(4,2,platd,2), lbasdz(4,2,plev), &
&      lbasiy(4,2,platd,2), lbassd(4,2,plevp1), phi(platd,2), sinlam(plon)
  INTEGER, INTENT (OUT) :: kdpmpf(pmap), kdpmph(pmap)

  ! pmap    Dimension of artificial evenly spaced vertical grid arrays
  ! etamid  Full-index hybrid-levels in vertical grid.
  ! etaint  Half-index hybrid-levels from sig(1/2) = etaint(1) = 0. to
  !         sig(plev+1/2) = etaint(plevp1) = 1.
  ! gravit  Gravitational constant.
  ! dlam    Length of increment in longitude grid.
  ! lam     Longitude values in the extended grid.
  ! phi     Latitude values in the extended grid.
  ! dphi    Interval between latitudes in the extended grid
  ! gw      Gauss weights for latitudes in the global grid.  (These sum
  !         to 2.0.)
  ! sinlam  Sine of longitudes in global grid (no extension points).
  ! coslam  Cosine of longitudes in global grid (no extension points).
  ! lbasdy  Weights for Lagrange cubic derivative estimates on the
  !         unequally spaced latitude grid
  ! lbasdz  Weights for Lagrange cubic derivative estimates on the
  !         unequally spaced vertical grid (corresponding to model
  !         full levels).
  ! lbassd  Weights for Lagrange cubic derivative estimates on the
  !         unequally spaced vertical grid (corresponding to model
  !         half levels).
  ! lbasiy  Weights for Lagrange cubic interpolation on the
  !         unequally spaced latitude grid
  ! detam   Increment between model mid-levels ("full" levels)
  ! detai   Increment between model interfaces ("half" levels).
  ! kdpmpf  Array of indicies of the model full levels which are mapped
  !         into an artificial evenly spaced vertical grid.  Used to aid
  !         in search for vertical position of departure point
  ! kdpmph  Array of indicies of the model half levels which are mapped
  !         into an artificial evenly spaced vertical grid.  Used to aid
  !         in search for vertical position of departure point
  ! cwava   Weight applied to global integrals 1./(plon*gravit)

  !  Local scalars: 
  INTEGER :: j, k, ih

  !  Local arrays: 
  REAL :: detailn(plevp1)    ! dlog(etaint)
  REAL :: detamln(plev)     ! dlog(etamid)
  REAL :: etailn(plevp1)     ! log(etaint)
  REAL :: etamln(plev)      ! log(etamid)

  !  External subroutines 
  EXTERNAL basdy, basdz, basiy, grdxy, vrtmap

  !  Intrinsic functions 
  INTRINSIC LOG


  !  Executable statements 

  ! Initialize extended horizontal grid coordinates.

  CALL grdxy(dlam,lam,phi,gw,sinlam,coslam)

  ! Basis functions for computing Lagrangian cubic derivatives
  ! on unequally spaced latitude and vertical grids.

  CALL basdy(phi(1,1),lbasdy(1,1,1,1))
  CALL basdy(phi(1,2),lbasdy(1,1,1,2))
  CALL basdz(plev,etamid,lbasdz)
  CALL basdz(plevp1,etaint,lbassd)

  ! Basis functions for computing weights for Lagrangian cubic
  ! interpolation on unequally spaced latitude grids.

  CALL basiy(phi(1,1),lbasiy(1,1,1,1))
  CALL basiy(phi(1,2),lbasiy(1,1,1,2))

  ! Compute interval lengths in latitudinal grid

  dphi(:,:)   = 0.
  dphibr      = 0.
  DO ih = 1,2
     DO j = 1, platd - 1
        dphi(j,ih) = phi(j+1,ih) - phi(j,ih)
        dphibr = MAX(ABS(dphi(j,ih)),dphibr)
     END DO
  END DO
  dphibr = 1./dphibr

  ! Compute interval lengths in vertical grids.

  DO k = 1, plev
    etamln(k) = LOG(etamid(k))
  END DO
  DO k = 1, plevp1
    etailn(k) = LOG(etaint(k))
  END DO
  DO k = 1, plev - 1
    detam(k) = etamid(k+1) - etamid(k)
    detamln(k) = etamln(k+1) - etamln(k)
  END DO
  DO k = 1, plev
    detai(k) = etaint(k+1) - etaint(k)
    detailn(k) = etailn(k+1) - etailn(k)
  END DO

  ! Build artificial evenly spaced vertical grid for use in determining
  ! vertical position of departure point.
  ! Build one grid for full model levels and one for half levels.

  CALL vrtmap(plev,pmap,etamln,detamln,kdpmpf)
  CALL vrtmap(plevp1,pmap,etailn,detailn,kdpmph)

  ! Compute moisture integration constant

  cwava = 1./(lc%nlon*gravit)

  RETURN
END SUBROUTINE grdini
