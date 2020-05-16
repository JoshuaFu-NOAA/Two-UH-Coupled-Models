!+ Compute weights for the calculation of derivative estimates
!+ $Id: basdz.f90,v 1.3 1998/10/28 12:27:18 m214003 Exp $

SUBROUTINE basdz(pkdim,sig,lbasdz)

  ! Description:
  !
  ! Compute weights for the calculation of derivative estimates
  !
  ! Method:
  !
  ! Compute weights for the calculation of derivative estimates at two
  ! center points of the four point stencil for each interval in the
  ! unequally spaced vertical grid (as defined by the array sig).
  ! Estimates are from differentiating a Lagrange cubic polynomial
  ! through the four point stencil.
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
 
  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: pkdim

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: sig(pkdim)

  !  Array arguments with intent(Out):
  REAL, INTENT (OUT) :: lbasdz(4,2,pkdim)

  ! pkdim   Number of grid points in vertical grid.
  ! sig     Sigma values in the vertical grid.
  ! lbasdz  Weights for derivative estimates based on Lagrange cubic
  !         polynomial on the unequally spaced vertical grid.
  !         If grid interval j is surrounded by a 4 point stencil,
  !         then the derivative at the "top" of the interval (smaller
  !         sigma value) uses the weights lbasdz(1,1,j),lbasdz(2,1,j),
  !         lbasdz(3,1,j), and lbasdz(4,1,j).  The derivative at the
  !         "bottom" of the interval uses lbasdz(1,2,j), lbasdz(2,2,j),
  !         lbasdz(3,2,j), and lbasdz(4,2,j).  (Recall the vertical
  !         level indices increase from the top of the atmosphere
  !         towards the bottom.)

  !  Local scalars: 
  INTEGER :: kk

  !  External subroutines 
  EXTERNAL lcdbas


  !  Executable statements 

  DO kk = 2, pkdim - 2
    CALL lcdbas(sig(kk-1),lbasdz(1,1,kk),lbasdz(1,2,kk))
  END DO

  RETURN
END SUBROUTINE basdz
