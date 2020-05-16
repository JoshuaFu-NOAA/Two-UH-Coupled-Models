!+ saves history files new version for echam3

SUBROUTINE savehis3(number)

  ! Description:
  !
  ! This routine saves history files for future use.
  !
  ! Method:
  !
  ! The history files are put together with "tar".
  ! The name of the tar file is: exp.number"re"yymm
  !
  ! *savehis* is called from *posts2*
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, September 1991, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_exception,     ONLY: finish
  USE mo_start_dataset, ONLY: nhf1
  USE mo_filename,      ONLY: nhm, nhy, yexp
  USE mo_doctor,        ONLY: nout
  USE mo_control,       ONLY: lamip2

  IMPLICIT NONE

  !  Scalar arguments 
  INTEGER :: number

  !  Local scalars: 
  INTEGER :: ilenc, ilenf1, ilenf2, istat, n, nunit
  CHARACTER (300) :: yocmnd
  CHARACTER (6) :: yodate
  CHARACTER (40) :: yofile1, yofile2
  CHARACTER (20) :: yofile3

  !  Local arrays: 
  CHARACTER (2) :: younit(8)

  !  External functions 
  INTEGER, EXTERNAL :: ichlen2, util_system

  !  Executable Statements 

  WRITE (yodate,'(I4.4,I2.2)') nhy, nhm

  yofile1 = yexp // '.'
  yofile2 = yexp(1:5) // 're' // yodate

  IF (lamip2) THEN
    yofile3 = yexp(:5)//'statrerun'
  ELSE
    yofile3 = ' '
  END IF

  ilenf1 = ichlen2(yofile1,40)
  ilenf2 = ichlen2(yofile2,40)

  DO n = 1, number
    nunit = nhf1 - 1 + n
    WRITE (younit(n),'(I2)') nunit
  END DO

  IF (number==7) THEN
    yocmnd = 'tar cvf ' // yofile2(:ilenf2) // ' ' // yofile1(:ilenf1) // &
             younit(1) // ' ' // yofile1(:ilenf1) // younit(2) // ' ' //  &
                  yofile1(:ilenf1) // younit(5) // ' ' //                 &
        yofile1(:ilenf1) // younit(6) // ' ' // yofile1(:ilenf1) // younit(7)
  ELSE
    yocmnd = 'tar cvf ' // yofile2(:ilenf2) // ' ' // yofile1(:ilenf1) // &
             younit(1) // ' ' // yofile1(:ilenf1) // younit(2) // ' ' //  &
             yofile1(:ilenf1) // younit(5) // ' ' //                      &
             yofile1(:ilenf1) // younit(6) // ' ' // yofile1(:ilenf1) //  &
             younit(7) // ' ' // yofile1(:ilenf1) // younit(8) // ' ' //  &
             yofile3
  END IF

  ilenc = ichlen2(yocmnd,300)

  istat = util_system(yocmnd(:ilenc))
  IF (istat/=0) THEN
    WRITE (nout,*) ' ERROR: tar of history files failed! Job aborted!'
    CALL finish('savehis3','Run terminated.')
  END IF

  WRITE (nout,*) number, ' history files saved as ' // yofile2(:ilenf2)

  RETURN
END SUBROUTINE savehis3
