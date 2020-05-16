!+ apply SCM0 limiter to derivative estimates.
!+ $Id: scm0.f90,v 1.6 1999/08/27 17:00:32 m214030 Exp $

SUBROUTINE scm0(n,deli,df1,df2)

  ! Description:
  !
  ! Apply SCM0 limiter to derivative estimates.
  !
  ! Authors:
  !
  ! Original version:  J. Olson
  ! Standardized:      J. Rosinski, June 1992
  ! Reviewed:          D. Williamson, P. Rasch, August 1992
  ! Modified:          U. Schlese,  DKRZ - Hamburg,  May 1992
  ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
  ! 
  ! for more details see file AUTHORS
  !

!  USE mo_grid

  IMPLICIT NONE

  !  Local parameters: 
  ! machine precision parameter to prevent rounding
  ! negatives in interpolant
  REAL, PARAMETER :: eps = 1.E-12
  ! factor applied in limiter
  REAL, PARAMETER :: fac = 3.*(1.-eps)

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: n

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: deli(n) 

  !  Array arguments with intent(InOut):
  REAL, INTENT (INOUT) :: df1(n), df2(n)

  ! n      Dimension of input arrays.
  ! deli   deli(i) is the discrete derivative on interval i, i.e.,
  ! deli(i) = ( f(i+1) - f(i) )/( x(i+1) - x(i) ).
  ! df1    df1(i) is the limited derivative at the left  edge of interval
  ! df2    df2(i) is the limited derivative at the right edge of interval

  !  Local scalars: 
  REAL :: tmp1  ! derivative factor
  REAL :: tmp2  ! abs(tmp1)
  INTEGER :: i

  !  Intrinsic functions 
  INTRINSIC ABS, MERGE

  !  Executable statements 

  DO i = 1, n
    tmp1 = fac*deli(i)
    tmp2 = ABS(tmp1)

    df1(i) = MERGE(0.,df1(i),-(deli(i)*df1(i))>=0.)
    df2(i) = MERGE(0.,df2(i),-(deli(i)*df2(i))>=0.)
    df1(i) = MERGE(tmp1,df1(i),tmp2-ABS(df1(i))<0.)
    df2(i) = MERGE(tmp1,df2(i),tmp2-ABS(df2(i))<0.)

  END DO

  RETURN
END SUBROUTINE scm0
