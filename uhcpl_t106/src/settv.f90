!+ get cpu's
!+ $Id: settv.f90,v 1.4 1998/12/08 09:25:25 m214003 Exp $

SUBROUTINE settv(icpu)

  USE mo_doctor, ONLY: nout

  IMPLICIT NONE

  !  Scalar arguments 
  INTEGER :: icpu

  !  Local scalars: 
  INTEGER :: envalue

  !  External subroutines 
  EXTERNAL util_igetenv


  !  Executable statements 

  ! Get NCPU environment variable

#ifdef sun
  CALL util_igetenv('PARALLEL',1,envalue)
#else
  CALL util_igetenv('NCPUS',1,envalue)
#endif
  IF (envalue==1) THEN
#ifdef sun
    WRITE (nout,*) ' Environment Variable PARALLEL set to 1.'
#else
    WRITE (nout,*) ' Environment Variable NCPUS set to 1.'
#endif
    icpu = 1
  ELSE
    icpu = envalue
  END IF
  RETURN
END SUBROUTINE settv
