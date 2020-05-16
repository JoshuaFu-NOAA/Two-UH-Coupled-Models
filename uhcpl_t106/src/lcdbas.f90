!+ calculate weights
!+ $Id: lcdbas.f90,v 1.4 1998/10/28 12:30:43 m214003 Exp $

SUBROUTINE lcdbas(grd,dbas2,dbas3)

  ! Description:
  !
  ! Calculate weights.
  !
  ! Method:
  !
  ! Calculate weights used to evaluate derivative estimates at the
  ! inner grid points of a four point stencil based on Lagrange
  ! cubic polynomial through four unequally spaced points.
  !
  ! Authors:
  !
  ! Standardized:      J. Rosinski, June 1992
  ! Reviewed:          D. Williamson, P. Rasch, August 1992
  ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
  ! 
  ! for more details see file AUTHORS
  !

  IMPLICIT NONE

  !  Array arguments with intent(In):
  REAL, INTENT (IN) ::  grd(4)

  !  Array arguments with intent(Out):
  REAL, INTENT (OUT) :: dbas2(4), dbas3(4)

  ! grd    Coordinate values of four points in stencil.
  ! dbas2  Derivatives of the four basis functions at grid point 2.
  ! dbas3  Derivatives of the four basis functions at grid point 3.

  !  Local scalars: 
  REAL :: x1, x2, x3, x4                        ! - grid values
  REAL :: x1mx2,x1mx3,x1mx4,x2mx3,x2mx4,x3mx4   ! - differences of grid values


  !  Executable atatements 
  x1 = grd(1)
  x2 = grd(2)
  x3 = grd(3)
  x4 = grd(4)
  x1mx2 = x1 - x2
  x1mx3 = x1 - x3
  x1mx4 = x1 - x4
  x2mx3 = x2 - x3
  x2mx4 = x2 - x4
  x3mx4 = x3 - x4

  dbas2(1) = x2mx3*x2mx4/(x1mx2*x1mx3*x1mx4)
  dbas2(2) = -1./x1mx2 + 1./x2mx3 + 1./x2mx4
  dbas2(3) = -x1mx2*x2mx4/(x1mx3*x2mx3*x3mx4)
  dbas2(4) = x1mx2*x2mx3/(x1mx4*x2mx4*x3mx4)

  dbas3(1) = -x2mx3*x3mx4/(x1mx2*x1mx3*x1mx4)
  dbas3(2) = x1mx3*x3mx4/(x1mx2*x2mx3*x2mx4)
  dbas3(3) = -1./x1mx3 - 1./x2mx3 + 1./x3mx4
  dbas3(4) = -x1mx3*x2mx3/(x1mx4*x2mx4*x3mx4)

  RETURN
END SUBROUTINE lcdbas
