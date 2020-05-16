!+ evaluate the partial Lagrangian cubic basis functions
!+ $Id: lcbas.f90,v 1.4 1998/10/28 12:30:41 m214003 Exp $

SUBROUTINE lcbas(grd,bas1,bas2)

  ! Description:
  !
  ! Evaluate the partial Lagrangian cubic basis functions (denominator
  ! only ) for the grid points and gather grid values
  !
  ! Authors:
  !
  ! Standardized:      J. Rosinski, June 1992
  ! Reviewed:          D. Williamson, P. Rasch, August 1992
  ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
  ! 
  ! for more details see file AUTHORS

  IMPLICIT NONE

  !  Array arguments with intent(In):
  REAL, INTENT (IN) ::  grd(4)   ! Grid stencil

  !  Array arguments with intent(Out):
  REAL, INTENT (OUT) :: bas1(4)  ! Grid values on stencil
  REAL, INTENT (OUT) :: bas2(4)  ! Lagrangian basis functions

  !  Local scalars: 
  !  Grid value differences used in weights
  REAL :: x0mx1, x0mx2, x0mx3, x1mx2, x1mx3, x2mx3


  !  Executable statements 
  x0mx1 = grd(1) - grd(2)
  x0mx2 = grd(1) - grd(3)
  x0mx3 = grd(1) - grd(4)
  x1mx2 = grd(2) - grd(3)
  x1mx3 = grd(2) - grd(4)
  x2mx3 = grd(3) - grd(4)

  bas1(1) = grd(1)
  bas1(2) = grd(2)
  bas1(3) = grd(3)
  bas1(4) = grd(4)

  bas2(1) = 1./(x0mx1*x0mx2*x0mx3)
  bas2(2) = -1./(x0mx1*x1mx2*x1mx3)
  bas2(3) = 1./(x0mx2*x1mx2*x2mx3)
  bas2(4) = -1./(x0mx3*x1mx3*x2mx3)

  RETURN
END SUBROUTINE lcbas
