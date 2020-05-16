!+ submits a job at end of run
!+ $Id: subjob.f90,v 1.9 2000/04/07 11:42:26 m214003 Exp $

SUBROUTINE subjob(kjob)

  ! Description:
  !
  ! Submits a job at end of run.
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, April 1995, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_filename,       ONLY: nhy, nhm, nhd, yexp, dir_limit, ypath
  USE mo_control,        ONLY: lsub
  USE mo_start_dataset,  ONLY: njin, njout
  USE mo_doctor,         ONLY: nout

  IMPLICIT NONE

  !  Scalar arguments 
  INTEGER :: kjob

  !  Local scalars: 
  INTEGER :: ILEN, istat, iypt
  LOGICAL :: loex
  CHARACTER   (2) :: yoday, yojn, yomn
  CHARACTER   (4) :: yoifil, yoyr
  CHARACTER (132) :: yoline
  CHARACTER   (7) :: yoofil

  !  External functions 
  INTEGER, EXTERNAL :: ichlen2, util_system


  !  Executable statements 

  WRITE (yoyr, '(I4.4)') nhy
  WRITE (yomn, '(I2.2)') nhm
  WRITE (yoday,'(I2.2)') nhd
  WRITE (yojn, '(I1)')   kjob

  yoifil = 'job' // yojn
  yoofil = 'subjob' // yojn

  INQUIRE (file=yoifil, exist=loex)
  IF (loex) THEN
     OPEN (njin, file=yoifil)
  ELSE
     WRITE (nout,*) ' *** WARNING: no input-file ' // yoifil // ' found.', &
                    ' job ' // yoofil // ' not submitted!'
     RETURN
  END IF

  INQUIRE (file=yoofil, exist=loex)
  IF (loex) THEN
     OPEN (njout, file=yoofil, status='OLD')
  ELSE
     OPEN (njout, file=yoofil, status='NEW')
  END IF

  DO
     READ (njin,'(A)',END=20) yoline
     ILEN = ichlen2(yoline,132)
     WRITE (njout,'(A)') yoline(1:ILEN)
     IF (yoline(1:5) == '#YEAR') THEN
        WRITE (njout,'(A)') 'YEAR=' // yoyr
     ELSE IF (yoline(1:6) == '#MONTH') THEN
        WRITE (njout,'(A)') 'MONTH=' // yomn
     ELSE IF (yoline(1:4) == '#DAY') THEN
        WRITE (njout,'(A)') 'DAY=' // yoday
     ELSE IF (yoline(1:4) == '#EXP') THEN
        WRITE (njout,'(A)') 'EXP=' // yexp
     ELSE IF (yoline(1:6) == '#DPATH') THEN
        DO iypt = dir_limit, 1, -1
           IF (ypath(iypt:iypt) /= ' ') EXIT
        END DO
        WRITE (njout,'(A)') 'DPATH=' // ypath(1:iypt)
     END IF
  END DO

20 CONTINUE

  CLOSE (njin)
  CLOSE (njout)

  IF (lsub) THEN
     istat = util_system ('qsub '//yoofil)
     IF (istat==0) THEN
        WRITE (nout,*) ' Job ' // yoofil // ' submitted'
     ELSE
        WRITE (nout,*) ' Error submitting ' // yoofil //   & 
                       '. "util_system" returned ', istat
     END IF
  END IF

  RETURN
END SUBROUTINE subjob
