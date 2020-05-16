!+
!+ $Id: extra.f90,v 1.11 1999/09/16 08:09:07 k202057 Exp $
FUNCTION ydate()

  IMPLICIT NONE

  !  Function return value 
  CHARACTER (8) :: ydate

  !  Local scalars: 
  CHARACTER (10) :: aclock, adate

  !  Intrinsic functions 
  INTRINSIC DATE_AND_TIME


  !  Executable statements 
  CALL DATE_AND_TIME(adate,aclock)
  ydate = adate(5:6) // '/' // adate(7:8) // '/' // adate(3:4)
  RETURN
END FUNCTION ydate

FUNCTION yclock()

  IMPLICIT NONE

  !  Function return value 
  CHARACTER (8) :: yclock

  !  Local scalars: 
  CHARACTER (10) :: aclock, adate

  !  Intrinsic functions 
  INTRINSIC DATE_AND_TIME


  !  Executable statements 
  CALL DATE_AND_TIME(adate,aclock)
  yclock = aclock(1:2) // ':' // aclock(3:4) // ':' // aclock(5:6)
  RETURN
END FUNCTION yclock

#if (! defined CRAY) || (defined SX) || (defined _CRAYMPP)
SUBROUTINE whenfgt(n,x,incr,target,index,nn)

  !
  ! Authors:
  !
  ! U. Schulzweida, MPI, September 1997
  !

  IMPLICIT NONE

  !  Scalar arguments 
  REAL :: target
  INTEGER :: incr, n, nn

  !  Array arguments 
  REAL :: x(*)
  INTEGER :: INDEX(*)

  !  Local scalars: 
  INTEGER :: ival, j


  !  Executable statements 

  ival = 0
  DO j = 1, n, incr
    IF (x(j)>target) THEN
      ival = ival + 1
      INDEX(ival) = j
    END IF
  END DO

  nn = ival

  RETURN
END SUBROUTINE whenfgt
#endif

#if (! defined CRAY) || (defined SX)
FUNCTION intmax(n,x,incr)

  !
  ! Authors:
  !
  ! U. Schulzweida, MPI, September 1997
  !

  IMPLICIT NONE

  !  Function Return value 
  INTEGER :: intmax

  !  Scalar arguments 
  INTEGER :: incr, n

  !  Array arguments 
  INTEGER :: x(*)

  !  Local scalars: 
  INTEGER :: ival, ivmax, j


  !  Executable statements 

  ival = 0
  ivmax = -1E9
  IF (incr>0) THEN
    DO j = 1, n, incr
      IF (x(j)>ivmax) THEN
        ivmax = x(j)
        ival = j
      END IF
    END DO
  ELSE IF (incr<0) THEN
    DO j = n, 1, incr
      IF (x(j)>ivmax) THEN
        ivmax = x(j)
        ival = j
      END IF
    END DO
  END IF

  intmax = ival

  RETURN
END FUNCTION intmax
#endif

#if (! defined CRAY) || (defined SX)
FUNCTION intmin(n,x,incr)

  !
  ! Authors:
  !
  ! U. Schulzweida, MPI, September 1997
  !

  IMPLICIT NONE

  !  Function return value 
  INTEGER :: intmin

  !  Scalar arguments 
  INTEGER :: incr, n

  !  Array arguments 
  INTEGER :: x(*)

  !  Local scalars: 
  INTEGER :: ival, ivmin, j


  !  Executable statements 

  ival = 0
  ivmin = 1E9
  IF (incr>0) THEN
    DO j = 1, n, incr
      IF (x(j)<ivmin) THEN
        ivmin = x(j)
        ival = j
      END IF
    END DO
  ELSE IF (incr<0) THEN
    DO j = n, 1, incr
      IF (x(j)>ivmin) THEN
        ivmin = x(j)
        ival = j
      END IF
    END DO
  END IF

  intmin = ival

  RETURN
END FUNCTION intmin
#endif

#if (! defined CRAY) || (defined SX)
FUNCTION ismax(n,x,incr)

  !
  ! Authors:
  !
  ! U. Schulzweida, MPI, September 1997
  !

  IMPLICIT NONE

  !  Function Return value 
  INTEGER :: ismax

  !  Scalar arguments 
  INTEGER :: incr, n

  !  Array arguments 
  REAL :: x(*)

  !  Local scalars: 
  REAL :: rvmax
  INTEGER :: ival, j


  !  Executable statements 

  ival = 0
  rvmax = -1.E20
  IF (incr>0) THEN
    DO j = 1, n, incr
      IF (x(j)>rvmax) THEN
        rvmax = x(j)
        ival = j
      END IF
    END DO
  ELSE IF (incr<0) THEN
    DO j = n, 1, incr
      IF (x(j)>rvmax) THEN
        rvmax = x(j)
        ival = j
      END IF
    END DO
  END IF

  ismax = ival

  RETURN
END FUNCTION ismax
#endif

#if (! defined CRAY) || (defined SX)
FUNCTION ismin(n,x,incr)

  !
  ! Authors:
  !
  ! U. Schulzweida, MPI, September 1997
  !

  IMPLICIT NONE

  !  Function return value 
  INTEGER :: ismin

  !  Scalar arguments 
  INTEGER :: incr, n

  !  Array arguments 
  REAL :: x(*)

  !  Local scalars: 
  REAL :: rvmin
  INTEGER :: ival, j


  !  Executable statements 

  ival = 0
  rvmin = 1.E10
  IF (incr>0) THEN
    DO j = 1, n, incr
      IF (x(j)<rvmin) THEN
        rvmin = x(j)
        ival = j
      END IF
    END DO
  ELSE IF (incr<0) THEN
    DO j = n, 1, incr
      IF (x(j)>rvmin) THEN
        rvmin = x(j)
        ival = j
      END IF
    END DO
  END IF

  ismin = ival

  RETURN
END FUNCTION ismin
#endif
