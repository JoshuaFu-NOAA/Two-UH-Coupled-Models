!+ 2nd scan over the latitude lines controlling the inverse Legendre
!  transforms.
!+ $Id: scan2.f90,v 1.14 2000/03/24 14:34:29 m214030 Exp $

SUBROUTINE scan2

  ! Description:
  !
  ! 2nd scan over the latitude lines controlling the inverse Legendre
  !
  ! Method:
  !
  ! This subroutine scans over the latitude lines to perform
  ! inverse *legendre transforms(*lti*).
  !
  ! *scan2* is called from *stepon*
  !
  ! The spectral components are located in long term storage
  ! arrays. *posts2* allocate and release buffers,
  ! and perform other tasks associated with the i/o.
  !
  ! Externals:
  ! *lti*       inverse legendre transforms
  ! *posts2*    input/output f-buffer
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, January 1995, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_control,       ONLY: lwtime, lcolumn
  USE mo_start_dataset, ONLY: nstep
  USE mo_call_trans,    ONLY: spectral_to_legendre

  IMPLICIT NONE

  !  External subroutines 

  EXTERNAL lti, posts2

  !  Executable statements 

!-- Transpose: spectral space -> Legendre space

  CALL spectral_to_legendre

!-- Inverse *Legendre transforms

  CALL lti

!-- 3. Complete the scan

!-- 3.1 Release buffers, etc

  IF (lwtime .AND. nstep/=-1) THEN
    IF (.not.lcolumn) CALL posts2
  END IF

END SUBROUTINE scan2



