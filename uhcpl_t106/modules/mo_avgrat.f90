MODULE mo_avgrat

  IMPLICIT NONE

  REAL, ALLOCATABLE :: avgrat(:,:,:)  ! (nlon,ngl,0:13)

CONTAINS

  SUBROUTINE readavgrat

    USE mo_doctor,        ONLY: nout
    USE mo_control,       ONLY: ngl, nlon
    USE mo_exception,     ONLY: finish
    USE mo_mpi,           ONLY: p_pe, p_io   
    USE mo_decomposition, ONLY: lc => local_decomposition, global_decomposition
    USE mo_transpose,     ONLY: scatter_gp
    USE mo_io

    REAL, ALLOCATABLE, TARGET :: zin(:,:,:)
    REAL, POINTER :: gl_avgrat(:,:,:)

    INTEGER    :: navgrat, i
    INTEGER    :: io_ngl, io_nlon


    navgrat = 91

    ! Allocate memory for avgrat per PE

    ALLOCATE (avgrat(lc%nglon, lc%nglat, 0:13))

    IF (p_pe == p_io) THEN

       WRITE (nout, '(/)')
       WRITE (nout,*) 'Reading avgrat from', ' unit ', navgrat

       ! Open file

       CALL IO_open_unit (navgrat, yearnc, IO_READ)
       IO_file_id = yearnc%nc_file_id

       ! Check resolution

       CALL IO_inq_dimid  (IO_file_id, 'lat', io_var_id)
       CALL IO_inq_dimlen (IO_file_id, io_var_id, io_ngl)
       CALL IO_inq_dimid  (IO_file_id, 'lon', io_var_id)
       CALL IO_inq_dimlen (IO_file_id, io_var_id, io_nlon)

       IF (io_nlon /= nlon .OR. io_ngl /= ngl) THEN
          WRITE(nerr,*) 'readavgrat: unexpected resolution ', io_nlon, io_ngl
          CALL finish ('readavgrat', 'unexpected resolution')
       END IF

       ! Allocate memory for avgrat global field

       ALLOCATE (zin(lc%nlon, lc%nlat, 0:13))

       ! Read data

       CALL IO_inq_varid (IO_file_id, 'VGRAT', io_var_id)
       CALL IO_get_var_double (IO_file_id, io_var_id, zin(:,:,1:12))

       ! Close file

       CALL IO_close(yearnc)

       zin(:,:,0)  = zin(:,:,12)
       zin(:,:,13) = zin(:,:,1)

    END IF

    NULLIFY (gl_avgrat)
    DO i = 0, 13
       IF (p_pe == p_io) gl_avgrat => zin(:,:,i:i)
       CALL scatter_gp (gl_avgrat, avgrat(:,:,i:i), global_decomposition)
    END DO

    IF (p_pe == p_io) THEN
       DEALLOCATE (zin)
    END IF

    RETURN
  END SUBROUTINE readavgrat

  SUBROUTINE clavgrat

    ! Description:
    !
    ! Passes climate vegetation ratio to atmosphere
    !
    ! Method:
    !
    ! This subroutine calculates the vegetation ratio for
    ! each time step and updates vgratm.
    !
    ! *clavgrat* is called from *gpc*.
    !
    ! Authors: 
    !
    ! U. Schulzweida, MPI, July 1999
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_memory_g3a,    ONLY: vgratm
    USE mo_control,       ONLY: nrow
    USE mo_rad_switches,  ONLY: nmonth
    USE mo_decomposition, ONLY: dc=>local_decomposition
    USE mo_timeint,       ONLY: wgt1, wgt2, nmw1, nmw2

    !  Local scalars: 
    INTEGER :: im, jrow, jn, nglon


    !  Executable Statements 

    jrow  = nrow(2)
    nglon = dc%nglon

    ! Update vegetation

    ! Annual cycle

    IF (nmonth == 0) THEN

       ! Interpolation in time

       DO jn = 1, nglon
          vgratm(jn,jrow) = wgt1*avgrat(jn,jrow,nmw1) + wgt2*avgrat(jn,jrow,nmw2)
       END DO

    ELSE

       ! Perpetual month

       im = nmonth
       DO jn = 1, nglon
          vgratm(jn, jrow) = avgrat(jn,jrow,im)
       END DO

    END IF

    RETURN
  END SUBROUTINE clavgrat

END MODULE mo_avgrat
