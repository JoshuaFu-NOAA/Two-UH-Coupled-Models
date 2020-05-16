!+ Compute weights for the calculation of derivative estimates
!+ $Id: basdy.f90,v 1.4 1999/07/19 09:38:44 m214030 Exp $

SUBROUTINE basdy(phi,lbasdy)

  ! Description:
  !
  ! Compute weights for the calculation of derivative estimates
  !
  ! Method:
  !
  ! Compute weights for the calculation of derivative estimates at the two
  ! center points of the four point stencil for each interval in the
  ! unequally spaced latitude grid. Estimates are from differentiating
  ! a Lagrange cubic polynomial through the four point stencil.
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

  USE mo_grid, only: jfirst, jlast, platd

  IMPLICIT NONE

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: phi(platd)

  !  Array arguments with intent(Out):
  REAL, INTENT (OUT) :: lbasdy(4,2,platd)

  ! phi     Latitude values in the extended grid.
  ! lbasdy  Weights for derivative estimates based on Lagrange cubic
  !         polynomial on the unequally spaced latitude grid.
  !         If grid interval j (in extended grid) is surrounded by
  !         a 4 point stencil, then the derivative at the "bottom"
  !         of the interval uses the weights lbasdy(1,1,j),
  !         lbasdy(2,1,j), lbasdy(3,1,j), and lbasdy(4,1,j).
  !         The derivative at the "top" of the interval
  !         uses lbasdy(1,2,j), lbasdy(2,2,j), lbasdy(3,2,j),
  !         and lbasdy(4,2,j).

  !  Local scalars: 
  INTEGER :: jj

  !  External subroutines 
  EXTERNAL lcdbas


  !  Executable statements 

  DO jj = jfirst, jlast
    CALL lcdbas(phi(jj-1),lbasdy(1,1,jj),lbasdy(1,2,jj))
  END DO

  RETURN
END SUBROUTINE basdy
