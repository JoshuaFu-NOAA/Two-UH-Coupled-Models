!+  Vertical derivative estimates for a vertical slice using Lagrangian
!   cubic formulas.
!+ $Id: cubzdr.f90,v 1.5 1999/08/27 17:00:19 m214030 Exp $

SUBROUTINE cubzdr(pidim,pkdim,f,lbasdz,dfz1,dfz2)

  ! Description:
  !
  ! Vertical derivative estimates for a vertical slice using Lagrangian
  ! cubic formulas.
  !
  ! Method:
  !
  ! Vertical derivative estimates for a vertical slice using Lagrangian
  ! cubic formulas.  Derivatives are set to zero at the top and bottom.
  ! At the "inner nodes" of the top and bottom intervals, a "one sided"
  ! estimate is used.
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

! USE mo_grid

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: pidim, pkdim

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: f(pidim,pkdim), lbasdz(4,2,pkdim)

  !  Array arguments with intent(Out):
  REAL, INTENT (OUT) :: dfz1(pidim,pkdim), dfz2(pidim,pkdim)

  ! pidim   Horizontal dimension of arrays. (pidim = plon)
  ! pkdim   Vertical dimension of arrays.
  ! f       Vertical slice of data for which derivative estimates are made
  ! lbasdz  Lagrangian cubic basis functions for evaluating the
  !         derivatives on the unequally spaced vertical grid.
  ! dfz1    dfz1 contains derivative estimates at the "top" edges of the
  !         intervals in the f array.
  ! dfz2    dfz2 contains derivative estimates at the "bottom" edges of
  !         the intervals in the f array.

  !  Local scalars: 
  INTEGER :: i, k   ! Indices


  !  Executable statements 

  DO k = 2, pkdim - 2
    DO i = 1, pidim

      ! Lagrangian derivative estimates (cubic) for the two center nodes in a
      ! four node stencil.

      dfz1(i,k) = lbasdz(1,1,k)*f(i,k-1) + lbasdz(2,1,k)*f(i,k) + &
&          lbasdz(3,1,k)*f(i,k+1) + lbasdz(4,1,k)*f(i,k+2)

      dfz2(i,k) = lbasdz(1,2,k)*f(i,k-1) + lbasdz(2,2,k)*f(i,k) + &
&          lbasdz(3,2,k)*f(i,k+1) + lbasdz(4,2,k)*f(i,k+2)
    END DO
  END DO

  ! Constrain derivatives to zero at top and bottom of vertical grid.
  ! At the interior nodes of the intervals at the top and bottom of the
  ! vertical grid, use the derivative estimate at that same node for the
  ! adjacent interval.  (This is a "one-sided" estimate for that node.)

  DO i = 1, pidim
    dfz1(i,1) = 0.0
    dfz2(i,1) = dfz1(i,2)
    dfz1(i,pkdim-1) = dfz2(i,pkdim-2)
    dfz2(i,pkdim-1) = 0.0
  END DO

  RETURN
END SUBROUTINE cubzdr
