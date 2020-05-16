!+ Compute Lagrangian cubic derivative estimates for data on an equally
!  spaced grid.
!+ $Id: cubxdr.f90,v 1.4 1998/10/28 12:28:05 m214003 Exp $

SUBROUTINE cubxdr(pidim,ibeg,len,dx,f,fxl,fxr)

  ! Description:
  !
  ! Compute Lagrangian cubic derivative estimates for data on an equally
  ! spaced grid.
  !
  ! Method:
  !
  ! Compute Lagrangian cubic derivative estimates for data on an equally
  ! spaced grid.  Suppose grid interval i is centered in a 4 point
  ! stencil consisting of grid points i-1, i, i+1, and i+2.  Then the
  ! derivative at the left edge of the interval (i.e., grid point i)
  ! is stored in fxl(i), and the derivative at the right edge of the
  ! interval (i.e., grid point i+1) is stored in fxr(i).  Note that
  ! fxl(i) is not necessarily equal to fxr(i-1) even though both of
  ! these values are estimates of the derivative at grid point i.
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
  REAL, INTENT (IN) :: dx
  INTEGER, INTENT (IN) :: ibeg, len, pidim

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: f(pidim)

  !  Array arguments with intent(Out):
  REAL, INTENT (OUT) :: fxl(pidim), fxr(pidim)

  ! pidim   Length of f, fxl, and fxr.
  ! ibeg    First interval of grid for which derivatives are computed.
  ! len     Number of grid intervals for which derivatives are computed.
  !         (There are pidim - 1 intervals between the pidim gridpoints
  !         represented in f, fxl, and fxr.)
  ! dx      Value of grid spacing.
  ! f       Values on equally spaced grid for which derivatives are
  !         computed.
  ! fxl     fxl(i) is the derivative at the left  edge of interval i.
  ! fxr     fxr(i) is the derivative at the right edge of interval i.

  !  Local scalars: 
  REAL :: rdx6        ! normalization weight
  INTEGER :: i, iend


  !  Executable statements 

  iend = ibeg + len - 1
  rdx6 = 1./(6.*dx)

  DO i = ibeg, iend
    fxl(i) = (-2.*f(i-1)-3.*f(i)+6.*f(i+1)-f(i+2))*rdx6
    fxr(i) = (f(i-1)-6.*f(i)+3.*f(i+1)+2.*f(i+2))*rdx6
  END DO

  RETURN
END SUBROUTINE cubxdr
