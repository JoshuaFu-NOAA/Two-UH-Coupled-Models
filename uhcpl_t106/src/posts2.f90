!+ complete scan 2 input/output.
!+ $Id: posts2.f90,v 1.24 1999/10/21 13:17:35 m214003 Exp $

SUBROUTINE posts2

  ! Description:
  !
  ! *posts2* is used to perform other non
  ! multi-tasked input/output related functions on leaving scan 2.
  !
  ! Method:
  !
  ! Any required final processing is done. In the case of work
  ! files on solid state storage device, this includes the output of
  ! grid point history files at write-up times.
  !
  ! Externals:
  ! *posts2* is called from *scan2*
  ! *posts2* calls: *IO_write_bufout*  - to write history files.
  !
  ! Reference:
  ! Further information is contained in a note on the i/o package
  !
  ! Authors:
  !
  ! J. K. Gibson, ECMWF, February 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1999, netCDF version
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_control,       ONLY: dtime, lg4x, ncbase, ntbase, nwlag
  USE mo_start_dataset, ONLY: nhf1, nstep, ly365
  USE mo_filename,      ONLY: nhm
  USE mo_constants,     ONLY: dayl
  USE mo_year,          ONLY: im2day, cd2dat
  USE mo_doctor,        ONLY: nout
  USE mo_io,            ONLY: IO_write_buffer
  USE mo_grib,          ONLY: close_grib_file
  USE mo_mpi,           ONLY: p_pe, p_io

  IMPLICIT NONE

  !  Local scalars: 
  INTEGER :: id, iday, im, iy, number

  !  External subroutines 
  EXTERNAL savehis3

  !  Intrinsic functions 
  INTRINSIC MOD


  !  Executable statements 


  ! Write restart fourier file

  CALL IO_write_buffer(nhf1)

  ! Close grib file

  IF (p_pe == p_io) THEN
     CALL close_grib_file
  END IF

  IF (p_pe == p_io) THEN
     WRITE (nout,'(a,i2)') &
          '           F buffer  written to history file  unit.', nhf1

     ! Save history files if required

     IF (MOD(nhm,nwlag)==0) THEN
        ! Check if last day in the month is reached
        iday = ncbase + (ntbase+dtime*nstep)/dayl + 0.01
        CALL cd2dat(iday,id,im,iy)
        IF (ly365) THEN
           IF (im2day(im,iy) == id) THEN
              IF (lg4x) THEN
                 number = 8
              ELSE
                 number = 7
              END IF
              CALL savehis3(number)
           END IF
        ELSE
           IF (MOD(id,30) == 0) THEN
              IF (lg4x) THEN
                 number = 8
              ELSE
                 number = 7
              END IF
              CALL savehis3(number)
           END IF
        END IF
     END IF
  END IF

END SUBROUTINE posts2
