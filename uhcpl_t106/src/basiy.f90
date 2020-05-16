!+ Compute weights used in Lagrange cubic polynomial interpolation
!+ $Id: basiy.f90,v 1.4 1999/07/19 09:38:44 m214030 Exp $

SUBROUTINE basiy(phi,lbasiy)

  ! Description:
  !
  ! Compute weights used in Lagrange cubic polynomial interpolation
  !
  ! Method:
  !
  ! Compute weights used in Lagrange cubic polynomial interpolation in
  ! the central interval of a four point stencil.  Done for each interval
  ! in the unequally spaced latitude grid.
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
  REAL, INTENT (IN) :: phi(platd)  ! grid values in extended gri

  !  Array arguments with intent(InOut):
  REAL, INTENT (OUT) :: lbasiy(4,2,platd)  ! Weights for Lagrange cubic interp

  !  Local scalars: 
  INTEGER :: jj

  !  External subroutines 
  EXTERNAL lcbas


  !  Executable statements 

  DO jj = jfirst, jlast
    CALL lcbas(phi(jj-1),lbasiy(1,1,jj),lbasiy(1,2,jj))
  END DO

  RETURN
END SUBROUTINE basiy
