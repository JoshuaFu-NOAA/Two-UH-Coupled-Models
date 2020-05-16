MODULE mo_alai

  IMPLICIT NONE

  REAL, ALLOCATABLE :: alai(:,:,:)  ! (nlon,ngl,0:13)

CONTAINS

  SUBROUTINE readalai

    USE mo_doctor,        ONLY: nout
    USE mo_control,       ONLY: ngl, nlon
    USE mo_exception,     ONLY: finish
    USE mo_mpi,           ONLY: p_pe, p_io   
    USE mo_decomposition, ONLY: lc => local_decomposition, global_decomposition
    USE mo_transpose,     ONLY: scatter_gp
    USE mo_io

    REAL, ALLOCATABLE, TARGET :: zin(:,:,:)
    REAL, POINTER :: gl_alai(:,:,:)

    INTEGER    :: nalai, i
    INTEGER    :: io_ngl, io_nlon


    nalai = 90

    ! Allocate memory for alai per PE

    ALLOCATE (alai(lc%nglon, lc%nglat, 0:13))

    IF (p_pe == p_io) THEN

       WRITE (nout, '(/)')
       WRITE (nout,*) 'Reading alai from', ' unit ', nalai

       ! Open file

       CALL IO_open_unit (nalai, yearnc, IO_READ)
       IO_file_id = yearnc%nc_file_id

       ! Check resolution

       CALL IO_inq_dimid  (IO_file_id, 'lat', io_var_id)
       CALL IO_inq_dimlen (IO_file_id, io_var_id, io_ngl)
       CALL IO_inq_dimid  (IO_file_id, 'lon', io_var_id)
       CALL IO_inq_dimlen (IO_file_id, io_var_id, io_nlon)

       IF (io_nlon /= nlon .OR. io_ngl /= ngl) THEN
          WRITE(nerr,*) 'readalai: unexpected resolution ', io_nlon, io_ngl
          CALL finish ('readalai', 'unexpected resolution')
       END IF

       ! Allocate memory for alai global field

       ALLOCATE (zin(lc%nlon, lc%nlat, 0:13))

       ! Read data

       CALL IO_inq_varid (IO_file_id, 'VLT', io_var_id)
       CALL IO_get_var_double (IO_file_id, io_var_id, zin(:,:,1:12))

       ! Close file

       CALL IO_close(yearnc)

       zin(:,:,0)  = zin(:,:,12)
       zin(:,:,13) = zin(:,:,1)

       zin(:,:,:)  = MAX(zin(:,:,:), 0.001)

    END IF

    NULLIFY (gl_alai)
    DO i = 0, 13
       IF (p_pe == p_io) gl_alai => zin(:,:,i:i)
       CALL scatter_gp (gl_alai, alai(:,:,i:i), global_decomposition)
    END DO

    IF (p_pe == p_io) THEN
       DEALLOCATE (zin)
    END IF

    RETURN
  END SUBROUTINE readalai

  SUBROUTINE clalai

    ! Description:
    !
    ! Passes climate leaf area index to atmosphere
    !
    ! Method:
    !
    ! This subroutine calculates the leaf area index for
    ! each time step and updates vltm.
    !
    ! *clalai* is called from *gpc*.
    !
    ! Authors: 
    !
    ! U. Schulzweida, MPI, July 1999
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_memory_g3a,    ONLY: vltm
    USE mo_control,       ONLY: nrow
    USE mo_rad_switches,  ONLY: nmonth
    USE mo_decomposition, ONLY: dc=>local_decomposition
    USE mo_timeint,       ONLY: wgt1, wgt2, nmw1, nmw2

    !  Local scalars: 
    INTEGER :: im, jrow, jn, nglon


    !  Executable Statements 

    jrow  = nrow(2)
    nglon = dc%nglon

    ! Update leaf area index

    ! Annual cycle

    IF (nmonth == 0) THEN

       ! Interpolation in time

       DO jn = 1, nglon
          vltm(jn,jrow) = wgt1*alai(jn,jrow,nmw1) + wgt2*alai(jn,jrow,nmw2)
       END DO

    ELSE

       ! Perpetual month

       im = nmonth
       DO jn = 1, nglon
          vltm(jn,jrow) = alai(jn,jrow,im)
       END DO

    END IF

    RETURN
  END SUBROUTINE clalai

END MODULE mo_alai
