!+ control the restart of the forecast
!+ $Id: restart.f90,v 1.9 1999/09/08 16:40:51 m214030 Exp $

SUBROUTINE restart

  ! Description:
  !
  ! Control the setting up of work files, and the restart of
  ! the forecast.
  !
  ! Method:
  !
  ! Appropriate routines are called to perform the required
  ! tasks.
  !
  ! Externals:
  ! *iorestart*  - routine to perform required i/o.
  ! *nnsc1*      - routine to perform scan 1 of the forecast.
  ! *n-sc-*      - routine to perform other forecast scans.
  !
  ! Authors:
  !
  ! J. K. Gibson, ECMWF, December 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_exception
  USE mo_parameters
  USE mo_control
  USE mo_start_dataset
  USE mo_doctor, ONLY: nerr

  IMPLICIT NONE

  !  External subroutines 
  EXTERNAL inhysi , iorestart


  !  Executable statements 

!-- 1. Position files and re-set variables

  IF (nstep>0) THEN
    nstep = nstep + 1
  END IF
  IF (nstep > nstop) THEN
    WRITE (nerr,*) ' NSTOP is less than current timestep:'
    WRITE (nerr,*) ' NSTOP= ', nstop, '  NSTEP= ', nstep, '  Job aborted!'
    CALL finish('restart','Run terminated.')
  END IF

  CALL iorestart

  nresum = nstep

!-- 2. Set initial constants

  CALL inhysi

END SUBROUTINE restart
