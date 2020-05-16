SUBROUTINE iorestart

  ! Description:
  !
  ! Reads netCDF history files for a resumed run.
  !
  ! Method:
  !
  ! *iorestart* positions data sets at the beginning of a rerun,
  ! writing data description records, and setting up necessary work
  ! files.
  !
  ! Information is written to the data description records of
  ! appropriate files, and work files are written if necessary.
  !
  !
  ! Authors:
  !
  ! L. Kornblueh, MPI, May 1999, f90 rewrite
  ! U. Schulzweida, MPI, May 1999, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_tracer,        ONLY: ntrac, nhtrac, xtini
  USE mo_doctor,        ONLY: nout
  USE mo_start_dataset, ONLY: nhg1, nhg2, nhg3, nhgl1, nhf1
  USE mo_io,            ONLY: IO_read_buffer
  USE mo_memory_gl,     ONLY: xt
  USE mo_memory_g1a,    ONLY: xtm1

  IMPLICIT NONE


  !  Executable statements 

  ! Restart from history files

  ! Read fourier buffer

  CALL IO_read_buffer(nhf1)

  ! Read slt buffer

  CALL IO_read_buffer(nhgl1)

  ! Read restart grid point files

  CALL IO_read_buffer(nhg1)
  CALL IO_read_buffer(nhg2)
  CALL IO_read_buffer(nhg3)

  ! initialize new optional tracer fields

  IF ( ntrac > nhtrac ) THEN
    CALL xtini(nhtrac+1, xt)
    xtm1(:,:,nhtrac+1:,:) = xt(:,:,nhtrac+1:,:)
  ENDIF

  WRITE (nout,'(/)')

END SUBROUTINE iorestart
