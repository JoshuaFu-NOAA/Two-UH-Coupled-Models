!+ define the "extended" grid used in the semi-Lagrangian transport scheme.
!+ $Id: grdxy.f90,v 1.7 1999/08/02 09:46:22 k202057 Exp $

SUBROUTINE grdxy(dlam,lam,phi,w,sinlam,coslam)

  ! Description:
  !
  ! Define the "extended" grid used in the semi-Lagrangian transport
  ! scheme.  The longitudes are equally spaced and the latitudes are
  ! Gaussian.  The global grid is extended to include "wraparound" points
  ! on all sides.
  !
  ! Method:
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

  USE mo_decomposition, ONLY: lc => local_decomposition
  USE mo_grid,          ONLY: istart, istop, jstart, jstop, &
&                             plat, platd, plon, plond, jintmx, nxpt, plato2
  USE mo_gaussgrid,     ONLY: gauaw ! module subroutine
  USE mo_control,       ONLY: ngl

  IMPLICIT NONE

  !  Scalar arguments with intent(Out):
  REAL, INTENT (OUT) :: dlam

  !  Array arguments with intent(Out):
  REAL, INTENT (OUT) :: coslam(plon), lam(plond), phi(platd,2),  &
&                       sinlam(plon), w(plat)

  ! dlam    Length of increment in longitude grid.
  ! lam     Longitude values in the extended grid.
  ! phi     Latitude values in the extended grid.
  ! w       Gauss weights for latitudes in the global grid.
  !         (These sum to 2.0 like the ones in CCM1.)
  ! sinlam  Sine of longitudes in global grid (no extension points).
  ! coslam  Cosine of longitudes in global grid (no extension points).

  !  Local scalars: 
  REAL    :: lam0, pi
  INTEGER :: i, ig, j
  INTEGER :: jg, jgn, myglat, ih, ioff, offset, ppid_current, plft_current
  INTEGER :: jstartg, jstopg

  !  Local arrays: 
  REAL :: wrk(platd)    ! work space

  ! Local global arrays
!  INTEGER, PARAMETER :: pglatd=lc%nlat + 2*jintmx + 2*nxpt
!  REAL :: phig(pglatd), wrkg(pglatd), wg(lc%nlat)
  INTEGER :: pglatd
  REAL :: phig(ngl+2*jintmx+2*nxpt)
  REAL :: wrkg(ngl+2*jintmx+2*nxpt)
  REAL :: wg(ngl+2*jintmx+2*nxpt)


  !  Intrinsic functions 
  INTRINSIC ASIN, ATAN, COS, REAL, SIN


  !  Executable statements 

  pglatd = ngl + 2*jintmx + 2*nxpt
  lam0 = 0.0
  pi = 4.*ATAN(1.)
  jstartg = jstart
  jstopg = jstart - 1 + lc%nlat

  ! Interval length in equally spaced longitude grid.

  dlam = 2.*pi/REAL(lc%nlon)

  ! Longitude values on extended grid.

!  ! compute offset for general case
!  offset = 0
!  ppid_current   = lc%pid
!  DO ioff=1,lc%pid_ln-1
!     plft_current         = gc(ppid_current)%pid_lft 
!     offset               = offset + gc(plft_current)%nglon
!     ppid_current         = plft_current
!  END DO

  DO i = 1, plond
     lam(i) = REAL(i - istart + lc%glons(1)-1)*dlam + lam0
!     lam(i) = REAL(i-istart + offset)*dlam + lam0
  END DO

  ! Compute Gauss latitudes and weights.  On return; phi contains the
  ! sine of the latitudes starting closest to the north pole and going
  ! toward the south; w contains the corresponding Gauss weights.

  CALL gauaw(phig,wg,lc%nlat)

  ! Reorder and compute latitude values.
  
  DO j = jstartg, jstopg
    wrkg(j) = ASIN(phig(jstopg-j+1))
  END DO
  phig(jstartg:jstopg)=wrkg(jstartg:jstopg)

  ! North and south poles.

  phig(jstartg-1) = -pi/2.0
  phig(jstopg+1) = pi/2.0

  ! Extend Gauss latitudes below south pole so that the spacing above
  ! the pole is symmetric, and phi is decreasing, i.e., phi < -pi/2

  IF (jstartg > 2) THEN
    DO j = 1, jstartg - 2
      phig(j) = -pi - phig(2*jstartg-2-j)
    END DO
  END IF

  ! Analogously for Northern Hemisphere

  IF (pglatd > jstopg+1) THEN
    DO j = jstopg + 2, pglatd
      phig(j) = pi - phig(2*jstopg+2-j)
    END DO
  END IF

  ! Pack into local phi and w arrays

  DO j=1,platd
     ! Southern Hemisphere
     jg = lc%glats(1) + j -1
     phi(j,1) = phig(jg)
     ! Northern hemisphere
     jgn = pglatd - lc%glats(1) - plato2 - 2*nxpt - 2*jintmx + j + 1
     phi(j,2) = phig(jgn)
  END DO
  DO j=1,plato2
     ! Southern Hemisphere
     jg = lc%glats(1) + j - 1
     w(j) = wg(jg)
     ! Northern Hemisphere
     jgn = lc%nlat - lc%glats(1) - plato2 + j + 1
     w(j+plato2) = wg(jgn)
  END DO

  ! Sine and cosine of longitude.

  ig = 0
  DO i = istart, istop
    ig = ig + 1
    sinlam(ig) = SIN(lam(i))
    coslam(ig) = COS(lam(i))
  END DO

  RETURN
END SUBROUTINE grdxy
