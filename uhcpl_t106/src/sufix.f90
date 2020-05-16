!+ set up mass fixer
!+ $Id: sufix.f90,v 1.8 2000/03/17 09:13:01 m214003 Exp $

SUBROUTINE sufix(kftype)

  ! Description:
  !
  ! Set up mass fixer
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, May 1995, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_exception,  ONLY: finish
  USE mo_parameters, ONLY: jps
  USE mo_grid,       ONLY: pcnst
  USE mo_tracer,     ONLY: lslt, nfix, nfixt, ntcode, ntrac
  USE mo_doctor,     ONLY: nout
  USE mo_mpi,        ONLY: p_pe, p_io

  IMPLICIT NONE

  !  Array arguments 
  INTEGER :: kftype(pcnst)

  !  Local scalars: 
  INTEGER :: ix, jt


  !  Executable statements 

  !-- 1. Determine fixer type for transported constituents

  ix = 0
  DO jt = 1, ntrac + jps

     IF (jt<=jps) THEN
        ix = ix + 1
        kftype(ix) = nfix(jt)
     ELSE
        IF (lslt(jt-jps)) THEN
           ix = ix + 1
           kftype(ix) = nfixt(jt-jps)
        END IF
     END IF

     IF (kftype(ix)<0 .OR. kftype(ix)>2) THEN
        WRITE (nout, *) ' Wrong fixer type for variable no. ', jt
        WRITE (nout, *) '   type= ', kftype(ix)
        CALL finish('sufix','Run terminated.')
     END IF
  END DO

  ! Print switches

  IF (p_pe == p_io) THEN

     WRITE (nout,'(/,4(a,/),a,i2,/,a,i2)') &
          ' -------------------------------------------------------', &
          ' SLT-setup:',                                              &
          '  Code   tracer-index   lslt   fixer-type',                &
          ' -------------------------------------------------------', &
          '  133         Q           T       ', nfix(1),              &
          '  153         X           T       ', nfix(2)
     DO jt = 1, ntrac
        WRITE (nout, '(i5,7x,i3,10x,l2,7x,i2)') &
             ntcode(jt), jt, lslt(jt), nfixt(jt)
     END DO

     WRITE (nout, '(a)') &
          ' -------------------------------------------------------'
  END IF

END SUBROUTINE sufix

