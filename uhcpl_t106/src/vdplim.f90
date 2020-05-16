!+ restrict vertical departure points to be between the top and bottom
!  sigma levels of the "full-" or "half-" level grid.
!+ $Id: vdplim.f90,v 1.5 1999/08/27 17:00:37 m214030 Exp $

SUBROUTINE vdplim(pkdim,sig,sigdp)

  ! Description:
  !
  ! Restrict vertical departure points to be between the top and bottom
  ! sigma levels of the "full-" or "half-" level grid.
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

  USE mo_grid, ONLY: pgls

  IMPLICIT NONE

  !  Local parameters: 
  ! machine precision to restrict dep 
  ! point to be inside last grid point
  REAL, PARAMETER :: eps = 1.E-12

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: pkdim

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: sig(pkdim)

  !  Array arguments with intent(InOut):
  REAL, INTENT (INOUT) :: sigdp(pgls)

  ! pkdim   Vertical dimension of "sig"
  ! sig     Sigma values at the "full" or "half" model levels
  ! sigdp   Sigma value at the trajectory endpoint or midpoint for each
  !         gridpoint in a vertical slice from the global grid.  This
  !         routine restricts those departure points to within the
  !         model's vertical grid.

  !  Local scalars: 
  INTEGER :: i


  !  Executable statements 

  DO i = 1, pgls
    IF (sigdp(i)<sig(1)) THEN
      sigdp(i) = sig(1)
    END IF
    IF (sigdp(i)>=sig(pkdim)) THEN
      sigdp(i) = sig(pkdim) - eps
    END IF
  END DO

  RETURN
END SUBROUTINE vdplim
