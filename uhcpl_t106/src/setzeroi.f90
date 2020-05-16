!+ resets accumulated fields
!+ $Id: setzeroi.f90,v 1.16 2000/01/05 14:21:34 m214003 Exp $

SUBROUTINE setzeroi

  ! Description:
  !
  ! Resets accumulated fields
  !
  ! Method:
  !
  ! All accumulated fiels are reset to default values
  ! this routine is called from the incore version of the model !
  !
  ! Locate every latitude row in the models wokspace.
  ! Then set all values to default.
  ! This is the incore version of setzero.
  !
  ! *setzeroi* is called from *nnsc1*.
  !
  ! Results:
  ! The accumulated fields are reset
  !
  ! Authors:
  !
  ! E. Kirk, MI, August 1990, original source
  ! U. Schlese, DKRZ, unknown, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_exception
  USE mo_parameters
  USE mo_start_dataset
  USE mo_control
  USE mo_io_tables
  USE mo_post
  USE mo_machine,        ONLY: prec
  USE mo_memory_g3b

  IMPLICIT NONE

  !  Local scalars: 
  REAL :: default
  INTEGER :: jfield, jx
  CHARACTER (8) :: yname

  !  Local arrays: 
  REAL, POINTER :: ptr3d(:,:,:)


  !  Executable statements 

  ! 1. Loop over all fields except g3x-fields

  DO jfield = 129, 256
    yname = yn(jfield)
    default = setval(jfield)
    IF (yname(1:1)==' ') CYCLE
    IF (ABS(default)<prec .AND. .NOT. laccu(jfield)) CYCLE

    CALL get_entry (g3b, yname, ptr3d)

    IF (nlevg3(jfield) == 1) THEN
       ptr3d(:,:,1) = default
    ELSE
       ptr3d(:,:,:) = default
    END IF
  END DO

  ! 2. Process extra g3-fields

  IF (ng3xp>0) THEN
    DO jx = 1, ng3xp
      WRITE(yname,'(a3,i2.2)' ) 'G3X',jx

      CALL get_entry (g3b, yname, ptr3d)

      IF (lxaccu(jx)) THEN
         ptr3d(:,:,:) = 0.
      END IF
    END DO
  END IF

  RETURN

END SUBROUTINE setzeroi
