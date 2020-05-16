!+
!+ $Id: ichlen2.f90,v 1.3 1998/09/21 08:10:47 m214003 Exp $

FUNCTION ichlen2(ytext,klen)

  IMPLICIT NONE

  !  Function return value 
  INTEGER :: ichlen2

  !  Scalar arguments 
  INTEGER :: klen
  CHARACTER (*) :: ytext

  !  Local scalars: 
  INTEGER :: j


  !  Executable statements 

  DO j = klen, 2, -1
    IF (ytext(j:j)/=' ') EXIT
  END DO

  ichlen2 = j

  RETURN
END FUNCTION ichlen2
