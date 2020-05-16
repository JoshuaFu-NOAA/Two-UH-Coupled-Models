!+ program master main program.

PROGRAM master

  ! Description:
  !
  ! Call the control subroutine (*control*).
  !
  ! Externals:
  !
  ! *control*   called to control the run.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, February 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_doctor, ONLY: nout
  USE mo_mpi,    ONLY: p_start, p_stop, p_pe, p_io

  IMPLICIT NONE

  REAL :: zwtime

  !  External functions 
  REAL, EXTERNAL :: util_walltime

  !  External subroutines 
  EXTERNAL control


  !  Executable statements 

  ! Initialize wallclock timer

  zwtime = util_walltime()

  ! Start MPI

  CALL p_start

  ! Print version

  IF (p_pe == p_io) THEN
     WRITE (nout,'(78("-"),/)')

     WRITE (nout,'(a,/,a,/,a,/)')            &
          '  ECHAM         - Version  4.0 ', &
          '  Parallel, f90 - Revision 6.02', &
          '  Copyright by Max-Planck-Institute for Meteorology, 2000'

     WRITE (nout,'(/,78("-"),/)')
  END IF

  ! Call *control*

  CALL control

  ! Stop MPI  

  CALL p_stop

END PROGRAM master
