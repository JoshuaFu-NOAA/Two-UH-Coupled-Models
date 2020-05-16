MODULE mo_netCDF

  USE mo_doctor,        ONLY: nerr
  USE mo_exception,     ONLY: finish, message
  USE mo_start_dataset, ONLY: ldebugio

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  TYPE netCDF_file

     LOGICAL :: nc_opened                           ! open = .true. or closed = .FALSE.

     INTEGER :: nc_file_id                          ! netCDF file id 
     INTEGER :: nc_access_mode                      ! access mode for that file
     INTEGER :: ncdims(NF_MAX_VAR_DIMS) 

     CHARACTER (NF_MAX_NAME) :: nc_creation_program ! name of this program
     CHARACTER (NF_MAX_NAME) :: nc_creation_user    ! who has run this program
     CHARACTER (NF_MAX_NAME) :: nc_creation_date    ! date of creation of netCDF file
     CHARACTER (NF_MAX_NAME) :: nc_binary_source    ! binary data type of src (CRAY/IEEE)
     CHARACTER (NF_MAX_NAME) :: nc_file_type        ! initital or restart file ...
     CHARACTER (NF_MAX_NAME) :: nc_file_name        ! nc file name
     CHARACTER (NF_MAX_NAME) :: nc_title
  END TYPE netCDF_file

  INTEGER, PARAMETER :: max_dim_name = 10

  TYPE IO_dim
     INTEGER :: dim_id
     INTEGER :: dim_len
     CHARACTER (max_dim_name) :: dim_name
  END TYPE IO_dim

  INTEGER, PARAMETER :: max_dim_ids = 200
  INTEGER            :: IO_ndim_ids
  TYPE (IO_dim)      :: IO_dim_ids(max_dim_ids)

CONTAINS

  FUNCTION IO_get_varindx (dimname) RESULT(indx)

    CHARACTER (*), INTENT(in) :: dimname
    INTEGER :: indx


    DO indx = 1, max_dim_ids
       IF (IO_dim_ids(indx)%dim_name == dimname) RETURN
    END DO

    DO indx = 1, max_dim_ids
       IF (IO_dim_ids(indx)%dim_name(1:1) /= ' ') THEN
          WRITE(nerr,*) 'IO_get_varindx', indx, IO_dim_ids(indx)%dim_name
       END IF
    END DO
    WRITE (nerr,*) 'element ',dimname,' not available ...'
    CALL finish('IO_get_varindx','IO_dim_ids error')

  END FUNCTION IO_get_varindx

  SUBROUTINE IO_init_dims

    USE mo_control,    ONLY: ngl, nhgl, nlon, nlp2, nlev, nlevp1, nsp, nvclev, nmp1
    USE mo_io_tables,  ONLY: ng3xp, ng3xl
    USE mo_tracer,     ONLY: ntrac

    INTEGER :: jx
    CHARACTER (8) :: yname


    IO_dim_ids( 1)%dim_len  =  ngl
    IO_dim_ids( 1)%dim_name = "ngl"
    IO_dim_ids( 2)%dim_len  =  nhgl
    IO_dim_ids( 2)%dim_name = "nhgl"
    IO_dim_ids( 3)%dim_len  =  nlon
    IO_dim_ids( 3)%dim_name = "nlon"
    IO_dim_ids( 4)%dim_len  =  nlp2
    IO_dim_ids( 4)%dim_name = "nlp2"
    IO_dim_ids( 5)%dim_len  =  nlev
    IO_dim_ids( 5)%dim_name = "nlev"
    IO_dim_ids( 6)%dim_len  =  nlevp1
    IO_dim_ids( 6)%dim_name = "nlevp1"
    IO_dim_ids( 7)%dim_len  =  nsp
    IO_dim_ids( 7)%dim_name = "nsp"
    IO_dim_ids( 8)%dim_len  =  nvclev
    IO_dim_ids( 8)%dim_name = "nvclev"
    IO_dim_ids( 9)%dim_len  =  2
    IO_dim_ids( 9)%dim_name = "n2"
    IO_dim_ids(10)%dim_len  =  nmp1
    IO_dim_ids(10)%dim_name = "nmp1"
    IO_dim_ids(11)%dim_len  =  ntrac
    IO_dim_ids(11)%dim_name = "nhtrac"
    IO_dim_ids(12)%dim_len  =  53
    IO_dim_ids(12)%dim_name = "nswitches"

    IO_ndim_ids = 12

    IF (ng3xp>0) THEN
       DO jx = 1, ng3xp
          WRITE(yname,'(a5,i2.2)' ) 'ng3xl',jx
          IO_dim_ids(IO_ndim_ids+jx)%dim_len   =  ng3xl(jx)
          IO_dim_ids(IO_ndim_ids+jx)%dim_name  =  yname
       END DO
    END IF

  END SUBROUTINE IO_init_dims

  SUBROUTINE IO_DEF_DIM (ncid, name, len, dimid)

    INTEGER :: ncid, len, dimid, status, indx
    CHARACTER*(*) name


    status = NF_DEF_DIM (ncid, name, len, dimid)
    IF (status /= NF_NOERR) THEN
       WRITE(nerr,*) 'IO_DEF_DIM :', ncid, name, len
       CALL message ('IO_DEF_DIM', NF_STRERROR(status))
       CALL finish  ('IO_DEF_DIM', 'Run terminated.')
    END IF

    DO indx = 1, max_dim_ids
       IF (name == IO_dim_ids(indx)%dim_name) THEN
          IO_dim_ids(indx)%dim_id = dimid
          RETURN
       END IF
    END DO

    WRITE (nerr,*) 'element ',name,' not available ...'
    CALL finish('IO_DEF_DIM','IO_dim_ids error')

  END SUBROUTINE IO_DEF_DIM

  SUBROUTINE IO_INQ_DIMID (ncid, name, dimid)

    INTEGER :: ncid, dimid, status
    CHARACTER*(*) name


    status = NF_INQ_DIMID (ncid, name, dimid)
    IF (status /= NF_NOERR) THEN
       WRITE(nerr,*) 'IO_INQ_DIMID :', name
       CALL message ('IO_INQ_DIMID', NF_STRERROR(status))
       CALL finish  ('IO_INQ_DIMID', 'Run terminated.')
    END IF

  END SUBROUTINE IO_INQ_DIMID

  SUBROUTINE IO_INQ_DIMLEN (ncid, dimid, len)

    INTEGER :: ncid, dimid, len, status


    status = NF_INQ_DIMLEN (ncid, dimid, len)
    IF (status /= NF_NOERR) THEN
       WRITE(nerr,*) 'IO_INQ_DIMLEN :', ncid, dimid
       CALL message ('IO_INQ_DIMLEN', NF_STRERROR(status))
       CALL finish  ('IO_INQ_DIMLEN', 'Run terminated.')
    END IF

  END SUBROUTINE IO_INQ_DIMLEN

  SUBROUTINE IO_INQ_VARID (ncid, name, varid)

    INTEGER :: ncid, varid, status
    CHARACTER*(*) name


    status = NF_INQ_VARID (ncid, name, varid)

    IF (ldebugio) THEN
       WRITE(nerr,*) 'IO_INQ_VARID : ','     Id=',ncid,' varid=',varid,'   ',name
    END IF

    IF (status /= NF_NOERR) THEN
       WRITE(nerr,*) 'IO_INQ_VARID :', ncid, name
       CALL message ('IO_INQ_VARID', NF_STRERROR(status))
       CALL finish  ('IO_INQ_VARID', 'Run terminated.')
    END IF

  END SUBROUTINE IO_INQ_VARID

  SUBROUTINE IO_DEF_VAR (ncid, name, xtype, nvdims, vdims, varid)

    INTEGER :: ncid, varid, xtype, nvdims
    CHARACTER*(*) name
    INTEGER :: vdims(*)
    INTEGER :: status


    status = NF_DEF_VAR (ncid, name, xtype, nvdims, vdims, varid)
    IF (status /= NF_NOERR) THEN
       WRITE(nerr,*) 'IO_DEF_VAR :', ncid, name, xtype, nvdims,vdims(1:nvdims)
       CALL message ('IO_DEF_VAR', NF_STRERROR(status))
       CALL finish  ('IO_DEF_VAR', 'Run terminated.')
    END IF

  END SUBROUTINE IO_DEF_VAR

  SUBROUTINE IO_GET_ATT_TEXT (ncid, varid, name, text)

    INTEGER :: ncid, varid
    CHARACTER*(*) name, text
    INTEGER :: status


    status = NF_GET_ATT_TEXT (ncid, varid, name, text)
    IF (status /= NF_NOERR) THEN
       WRITE(nerr,*) 'IO_GET_ATT_TEXT :', ncid, varid, name
       CALL message ('IO_GET_ATT_TEXT', NF_STRERROR(status))
       CALL finish  ('IO_GET_ATT_TEXT', 'Run terminated.')
    END IF

  END SUBROUTINE IO_GET_ATT_TEXT

  SUBROUTINE IO_PUT_ATT_TEXT (ncid, varid, name, text)

    INTEGER :: ncid, varid
    CHARACTER*(*) name, text
    INTEGER :: status
    INTEGER :: len


    DO len = 80-1, 1, -1
       IF (text(len:len) /= ' ') EXIT
    END DO

    status = NF_PUT_ATT_TEXT (ncid, varid, name, len, text)
    IF (status /= NF_NOERR) THEN
       WRITE(nerr,*) 'IO_PUT_ATT_TEXT :', ncid, varid, name, len, text
       CALL message ('IO_PUT_ATT_TEXT', NF_STRERROR(status))
       CALL finish  ('IO_PUT_ATT_TEXT', 'Run terminated.')
    END IF

  END SUBROUTINE IO_PUT_ATT_TEXT

  SUBROUTINE IO_GET_ATT_INT (ncid, varid, name, ival)

    INTEGER :: ncid, varid
    CHARACTER*(*) name
    INTEGER :: ival
    INTEGER :: status


    status = NF_GET_ATT_INT (ncid, varid, name, ival)
    IF (status /= NF_NOERR) THEN
       WRITE(nerr,*) 'IO_GET_ATT_INT :', ncid, varid, name
       CALL message ('IO_GET_ATT_INT', NF_STRERROR(status))
       CALL finish  ('IO_GET_ATT_INT', 'Run terminated.')
    END IF

  END SUBROUTINE IO_GET_ATT_INT

  SUBROUTINE IO_PUT_ATT_INT (ncid, varid, name, ival)

    INTEGER :: ncid, varid
    CHARACTER*(*) name
    INTEGER :: ival
    INTEGER :: xtype, len
    INTEGER :: status


    len = 1
    xtype = NF_INT
    status = NF_PUT_ATT_INT (ncid, varid, name, xtype, len, ival)
    IF (status /= NF_NOERR) THEN
       WRITE(nerr,*) 'IO_PUT_ATT_INT :', ncid, varid, name
       CALL message ('IO_PUT_ATT_INT', NF_STRERROR(status))
       CALL finish  ('IO_PUT_ATT_INT', 'Run terminated.')
    END IF

  END SUBROUTINE IO_PUT_ATT_INT

  SUBROUTINE IO_GET_ATT_DOUBLE (ncid, varid, name, rval)

    INTEGER :: ncid, varid
    CHARACTER*(*) name
    REAL :: rval
    INTEGER :: status


    status = NF_GET_ATT_DOUBLE (ncid, varid, name, rval)
    IF (status /= NF_NOERR) THEN
       WRITE(nerr,*) 'IO_GET_ATT_DOUBLE :', ncid, varid,  name
       CALL message ('IO_GET_ATT_DOUBLE', NF_STRERROR(status))
       CALL finish  ('IO_GET_ATT_DOUBLE', 'Run terminated.')
    END IF

  END SUBROUTINE IO_GET_ATT_DOUBLE

  SUBROUTINE IO_PUT_ATT_DOUBLE (ncid, varid, name, rval)

    INTEGER :: ncid, varid
    CHARACTER*(*) name
    REAL :: rval
    INTEGER :: xtype, len
    INTEGER :: status


    len = 1
    xtype = NF_DOUBLE

    status = NF_PUT_ATT_DOUBLE (ncid, varid, name, xtype, len, rval)
    IF (status /= NF_NOERR) THEN
       WRITE(nerr,*) 'IO_PUT_ATT_DOUBLE :', ncid, varid, name
       CALL message ('IO_PUT_ATT_DOUBLE', NF_STRERROR(status))
       CALL finish  ('IO_PUT_ATT_DOUBLE', 'Run terminated.')
    END IF

  END SUBROUTINE IO_PUT_ATT_DOUBLE

  SUBROUTINE IO_GET_VAR_DOUBLE (ncid, varid, dvals)

    INTEGER :: ncid, varid
    REAL    :: dvals(*)
    INTEGER :: status


    status = NF_GET_VAR_DOUBLE (ncid, varid, dvals)

    IF (ldebugio) THEN
       WRITE(nerr,*) 'IO_GET_VAR_DOUBLE :',' Id=',ncid,' varid=',varid
    END IF

    IF (status /= NF_NOERR) THEN
       WRITE(nerr,*) 'IO_GET_VAR_DOUBLE :', ncid, varid
       CALL message ('IO_GET_VAR_DOUBLE', NF_STRERROR(status))
       CALL finish  ('IO_GET_VAR_DOUBLE', 'Run terminated.')
    END IF

  END SUBROUTINE IO_GET_VAR_DOUBLE

  SUBROUTINE IO_GET_VARA_DOUBLE (ncid, varid, start, count, dvals)

    INTEGER :: ncid, varid
    INTEGER :: start(*), COUNT(*)
    REAL :: dvals(*)
    INTEGER :: status

    status = NF_GET_VARA_DOUBLE (ncid, varid, start, count, dvals)
    IF (status /= NF_NOERR) THEN
       WRITE(nerr,*) 'IO_GET_VARA_DOUBLE :', ncid, varid, start(1:4), COUNT(1:4)
       CALL message ('IO_GET_VARA_DOUBLE', NF_STRERROR(status))
       CALL finish  ('IO_GET_VARA_DOUBLE', 'Run terminated.')
    END IF

  END SUBROUTINE IO_GET_VARA_DOUBLE

  SUBROUTINE IO_PUT_VARA_DOUBLE (ncid, varid, start, count, dvals)

    INTEGER :: ncid, varid
    INTEGER :: start(*), COUNT(*)
    REAL :: dvals(*)
    INTEGER :: status
    CHARACTER (80) :: varname


    status = NF_PUT_VARA_DOUBLE (ncid, varid, start, count, dvals)
    IF (status /= NF_NOERR) THEN
       status = NF_INQ_VARNAME (ncid, varid, varname)
       WRITE(nerr,*) 'IO_PUT_VARA_DOUBLE :', ncid, varid, varname
       CALL message ('IO_PUT_VARA_DOUBLE', NF_STRERROR(status))
       CALL finish  ('IO_PUT_VARA_DOUBLE', 'Run terminated.')
    END IF

  END SUBROUTINE IO_PUT_VARA_DOUBLE

  SUBROUTINE IO_PUT_VAR_DOUBLE (ncid, varid, dvals)

    INTEGER :: ncid, varid
    REAL :: dvals(*)
    INTEGER :: status
    CHARACTER (80) :: varname


    IF (ldebugio) THEN
       status = NF_INQ_VARNAME (ncid, varid, varname)
       WRITE(nerr,*) 'IO_PUT_VAR_DOUBLE :', ncid, varid, varname
       IF (status /= NF_NOERR) THEN
          WRITE(nerr,*) 'IO_PUT_VAR_DOUBLE :', ncid, varid
          CALL message ('IO_PUT_VAR_DOUBLE', NF_STRERROR(status))
          CALL finish  ('IO_PUT_VAR_DOUBLE', 'Run terminated.')
       END IF
    END IF

    status = NF_PUT_VAR_DOUBLE (ncid, varid, dvals)
    IF (status /= NF_NOERR) THEN
       WRITE(nerr,*) 'IO_PUT_VAR_DOUBLE :', ncid, varid
       CALL message ('IO_PUT_VAR_DOUBLE', NF_STRERROR(status))
       CALL finish  ('IO_PUT_VAR_DOUBLE', 'Run terminated.')
    END IF

  END SUBROUTINE IO_PUT_VAR_DOUBLE

  SUBROUTINE IO_ENDDEF (ncid)

    INTEGER :: ncid, status


    status = NF_ENDDEF (ncid)
    IF (status /= NF_NOERR) THEN
       CALL message ('IO_ENDDEF', NF_STRERROR(status))
       CALL finish  ('IO_ENDDEF', 'Run terminated.')
    END IF

  END SUBROUTINE IO_ENDDEF

END MODULE mo_netCDF
