!+ complete scan 1 input/output.
!+ $Id: posts1.f90,v 1.21 1999/07/20 08:48:33 m214030 Exp $

SUBROUTINE posts1

  ! Description:
  !
  ! *posts1* is used to release buffers and perform other non
  ! multi-tasked input/output related functions on leaving scan 1.
  !
  ! Method:
  !
  ! Any required final processing is done. in the case of work
  ! files on solid state storage device, this includes the output of
  ! grid point history files at write-up times. finally, the scan 1
  ! buffers are released.
  !
  ! Externals:
  ! *posts1* is called from *scan1sl*
  ! *posts1* calls: *IO_write_buffer*  - to write history files.
  !
  ! Reference:
  ! Further information is contained in a note on the i/o package
  !
  ! Authors:
  !
  ! J. K. Gibson, ECMWF, February 1982, original source
  ! U. Schlese, DKRZ, March 1994, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1999, netCDF version
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_parameters
  USE mo_io_tables
  USE mo_control
  USE mo_start_dataset
  USE mo_doctor, ONLY: nout
  USE mo_io
  USE mo_mpi

  IMPLICIT NONE


  !  Executable statements 

  ! Write restart grid point files

  CALL IO_write_buffer(nhg1)
  CALL IO_write_buffer(nhg2)
  CALL IO_write_buffer(nhg3)

  IF (p_pe == p_io) THEN
     WRITE (nout,'(a,i2,a,i2,a,i2)') &
          &       '  G1, G2, G3 buffers written to history files unit.', nhg1, &
          &       ', unit.', nhg2, ', unit.', nhg3

     IF (lg4x) THEN
        WRITE (nout,'(a,i2)') &
             &       '          G4 buffer  written to history file  unit.', nhg4
     END IF
  END IF

  ! Write restart slt file

  CALL IO_write_buffer(nhgl1)

  IF (p_pe == p_io) THEN
     WRITE (nout,'(a,i2)') &
          &       '          GL buffer  written to history file  unit.', nhgl1
  END IF

END SUBROUTINE posts1
