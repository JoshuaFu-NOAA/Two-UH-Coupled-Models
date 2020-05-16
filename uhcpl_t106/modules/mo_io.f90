MODULE mo_io

  USE mo_netCDF
  USE mo_parameters,    ONLY: jpg3xf
  USE mo_io_tables,     ONLY: ng3xl, ng3xp
  USE mo_exception,     ONLY: message
  USE mo_start_dataset, ONLY: ldebugio
  USE mo_mpi,           ONLY: p_io, p_pe, p_bcast
  USE mo_doctor,        ONLY: nerr, nout

  IMPLICIT NONE

  TYPE (netCDF_file), SAVE :: yearnc
  TYPE (netCDF_file), SAVE :: ini_ozon
  TYPE (netCDF_file), SAVE :: ini_field
  TYPE (netCDF_file), SAVE :: sstnc0, sstnc1, sstnc2
  TYPE (netCDF_file), SAVE :: sicnc0, sicnc1, sicnc2
  TYPE (netCDF_file), SAVE :: header
  TYPE (netCDF_file), SAVE :: ini_surf
  TYPE (netCDF_file), SAVE :: ini_spec
  TYPE (netCDF_file), SAVE :: restart(31:38)

  REAL, ALLOCATABLE, SAVE :: vlon(:), vlat(:)   !fu++

  ! IO variables

  INTEGER, PARAMETER :: IO_READ = 1, IO_WRITE = 2

  INTEGER :: IO_file_id, IO_var_id, IO_dim_id
  INTEGER :: IO_dims(4)
  INTEGER :: IO_timestep = -1

  ! time variables

  INTEGER :: forecast_date, forecast_time, verification_date, verification_time

  INTEGER :: no_ens = 0   ! no of ensemble

CONTAINS

  SUBROUTINE IO_close(info)

    TYPE (netCDF_file), INTENT(INOUT)  :: info
    INTEGER :: status


    IF (p_pe == p_io) THEN

       IF (ldebugio) THEN
          WRITE(nerr,'(A,A12,A,I3)') ' IO_close : ',info%nc_file_name, &
                                     ' Id=',info%nc_file_id
       END IF

       status = NF_CLOSE(info%nc_file_id)

       IF (status == NF_NOERR) THEN
          info%nc_opened = .FALSE.
       ELSE
          CALL message ('IO_close', NF_STRERROR(status));
          CALL finish  ('IO_close', 'Run terminated.')
       END IF

    END IF

  END SUBROUTINE IO_close

  SUBROUTINE IO_open(filename, info, mode)

    INTEGER,            INTENT(IN)     :: mode
    CHARACTER (*),      INTENT(IN)     :: filename
    TYPE (netCDF_file), INTENT(INOUT)  :: info

    INTEGER :: status, ncmode


    IF (p_pe == p_io) THEN
       IF (info%nc_opened) THEN
          WRITE(nerr,*) 'IO_open : file ',info%nc_file_name,' already open'
          CALL finish  ('IO_open', 'Run terminated.')
       END IF

       IF (mode == IO_READ) THEN
          ncmode = NF_NOWRITE
       ELSE IF (mode == IO_WRITE) THEN
          ncmode = NF_WRITE
       ELSE
          CALL message ('IO_open', 'unexpected mode');
          CALL finish  ('IO_open', 'Run terminated.')
       END IF

       info%nc_file_name = filename

       IF (mode == IO_READ) THEN
          status = NF_OPEN (filename, ncmode, info%nc_file_id)
       ELSE
          status = NF_CREATE (filename, NF_CLOBBER, info%nc_file_id)
       END IF

       IF (ldebugio) THEN
          WRITE(nerr,*) 'IO_open  : ', filename, ' mode=', mode, &
                        ' ncmode=', ncmode, ' Id=', info%nc_file_id
       END IF

       IF (status == NF_NOERR) THEN
          info%nc_opened = .TRUE.
       ELSE
          CALL message ('IO_open', filename);
          CALL message ('IO_open', NF_STRERROR(status));
          CALL finish  ('IO_open', 'Run terminated.')
       END IF

    END IF

  END SUBROUTINE IO_open

  SUBROUTINE IO_open_unit(unit, info, mode)

    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(IN) :: mode
    TYPE (netCDF_file), INTENT(INOUT)  :: info
    CHARACTER (13) :: filename


    IF (p_pe == p_io) THEN
       IF (info%nc_opened) THEN
          WRITE(nerr,*) 'IO_open_unit: unit ',unit,' allready assigned to ', &
                         info%nc_file_name
          CALL finish  ('IO_open_unit', 'Run terminated.')
       END IF

       IF (unit < 10 .OR. unit > 99) THEN
          WRITE(nerr,*) 'IO_open_unit: unit ',unit,' out of range'
          CALL finish  ('IO_open_unit', 'Run terminated.')
       END IF
    END IF

!ik++
!fu++ change for large ensembles
    if ( no_ens > 0) then
       WRITE (filename,'(A5,I3.3,a1,i2.2)') 'unit.',no_ens,'.',unit 
    else
       WRITE (filename,'(A5,I2.2)') 'unit.',unit 
    end if
!ik++

    CALL IO_open(filename, info, mode)

  END SUBROUTINE IO_open_unit

  SUBROUTINE IO_write_header(io_unit, io_list)

    USE mo_doctor,        ONLY: ylabel1, ylabel2, ylabel3, ylabel4, &
                                ylabel5, ylabel6, ylabel7, ylabel8
    USE mo_control,       ONLY: nvclev, twodt, ntbase, ncbase, vct, dtime, nm, nn, nk
    USE mo_year,          ONLY: ic2ymd, isec2hms
    USE mo_start_dataset, ONLY: nhg3, nstep
    USE mo_filename,      ONLY: yomdn
    USE mo_linked_list,   ONLY: list_element

    INTEGER, INTENT(IN) :: io_unit
    INTEGER :: ifcday, ifcsec
    INTEGER :: idim, ndim
    REAL :: zfcsec
    CHARACTER (10) :: io_name

    INTEGER :: io_vlat_id, io_vlon_id, io_vcta_id, io_vctb_id 
    INTEGER :: io_nswitches_id
    REAL, ALLOCATABLE :: restart_diag(:)
    TYPE (list_element), POINTER :: io_list
    INTEGER :: nswitches


    IF (p_pe == p_io) THEN

       header = restart(io_unit)
       IO_file_id = header%nc_file_id

       DO idim = 1, IO_ndim_ids
          CALL IO_DEF_DIM (IO_file_id, IO_dim_ids(idim)%dim_name, &
                           IO_dim_ids(idim)%dim_len,  &
                           IO_dim_ids(idim)%dim_id)
       END DO

       IF (ng3xp>0 .AND. io_unit==nhg3) THEN
          DO idim = IO_ndim_ids+1, IO_ndim_ids+ng3xp
             CALL IO_DEF_DIM (IO_file_id, IO_dim_ids(idim)%dim_name, &
                              IO_dim_ids(idim)%dim_len,  &
                              IO_dim_ids(idim)%dim_id)
          END DO
       END IF


       WRITE(header%nc_file_type,'(A,I2,A)') 'Restart history file (unit:',io_unit,')'
       header%nc_binary_source    = 'IEEE'
       header%nc_creation_program = yomdn
       header%nc_creation_user    = ylabel7(2:)
       header%nc_creation_date    = ylabel6(2:)

       !    CALL IO_PUT_ATT_TEXT (IO_file_id, NF_GLOBAL, 'title',  header%nc_title)
       CALL IO_PUT_ATT_TEXT (IO_file_id, NF_GLOBAL, 'file_type',   header%nc_file_type)
       CALL IO_PUT_ATT_TEXT (IO_file_id, NF_GLOBAL, 'source_type', header%nc_binary_source)
       CALL IO_PUT_ATT_TEXT (IO_file_id, NF_GLOBAL, 'history',     header%nc_creation_program)
       CALL IO_PUT_ATT_TEXT (IO_file_id, NF_GLOBAL, 'user',        header%nc_creation_user)
       CALL IO_PUT_ATT_TEXT (IO_file_id, NF_GLOBAL, 'created',     header%nc_creation_date)

       CALL IO_PUT_ATT_TEXT (IO_file_id, NF_GLOBAL, 'label_1', ylabel1(2:))
       CALL IO_PUT_ATT_TEXT (IO_file_id, NF_GLOBAL, 'label_2', ylabel2(2:))
       CALL IO_PUT_ATT_TEXT (IO_file_id, NF_GLOBAL, 'label_3', ylabel3(2:))
       CALL IO_PUT_ATT_TEXT (IO_file_id, NF_GLOBAL, 'label_4', ylabel4(2:))
       CALL IO_PUT_ATT_TEXT (IO_file_id, NF_GLOBAL, 'label_5', ylabel5(2:))
       CALL IO_PUT_ATT_TEXT (IO_file_id, NF_GLOBAL, 'label_6', ylabel6(2:))
       CALL IO_PUT_ATT_TEXT (IO_file_id, NF_GLOBAL, 'label_7', ylabel7(2:))
       CALL IO_PUT_ATT_TEXT (IO_file_id, NF_GLOBAL, 'label_8', ylabel8(2:))

       ! put data reference times

       forecast_time = isec2hms(ntbase)
       forecast_date = ic2ymd(ncbase)

       zfcsec = (nstep+1)*twodt*0.5+ntbase
       ifcday = INT(zfcsec/(24*3600))
       ifcsec = INT(MOD(zfcsec,REAL(24*3600)))

       verification_time = isec2hms(ifcsec)
       verification_date = ic2ymd(ncbase+ifcday)

       CALL IO_PUT_ATT_INT (IO_file_id, NF_GLOBAL, 'fdate', forecast_date)
       CALL IO_PUT_ATT_INT (IO_file_id, NF_GLOBAL, 'ftime', forecast_time)
       CALL IO_PUT_ATT_INT (IO_file_id, NF_GLOBAL, 'vdate', verification_date)
       CALL IO_PUT_ATT_INT (IO_file_id, NF_GLOBAL, 'vtime', verification_time)

       ! put spherical truncations ...

       CALL IO_PUT_ATT_INT (IO_file_id, NF_GLOBAL, 'spherical_truncation_n', nn)
       CALL IO_PUT_ATT_INT (IO_file_id, NF_GLOBAL, 'spherical_truncation_m', nm)
       CALL IO_PUT_ATT_INT (IO_file_id, NF_GLOBAL, 'spherical_truncation_k', nk)

       ! put nstep

       CALL IO_PUT_ATT_INT (IO_file_id, NF_GLOBAL, 'nstep', nstep)

       ! put timestep

       CALL IO_PUT_ATT_DOUBLE (IO_file_id, NF_GLOBAL, 'timestep', dtime)

       ! tracer

       !    IF (nhtrac == 0) THEN
       !       CALL IO_PUT_ATT_TEXT (IO_file_id, NF_GLOBAL, 'tracer_definition', 'No tracer available')
       !    ELSE
       !       CALL IO_PUT_ATT_TEXT (IO_file_id, NF_GLOBAL, 'tracer_definition', 'Tracer available')
       !    END IF
       IF (ng3xp>0 .AND. io_unit==nhg3) THEN
          CALL IO_PUT_ATT_INT (IO_file_id, NF_GLOBAL, 'ng3xp', ng3xp)
       END IF

       IO_dims(1) = IO_dim_ids( IO_get_varindx("ngl") )%dim_id
       CALL IO_DEF_VAR (IO_file_id, "lat", NF_DOUBLE, 1, IO_dims, io_vlat_id)
       IO_dims(1) = IO_dim_ids( IO_get_varindx("nlon") )%dim_id
       CALL IO_DEF_VAR (IO_file_id, "lon", NF_DOUBLE, 1, IO_dims, io_vlon_id)
       IO_dims(1) = IO_dim_ids( IO_get_varindx("nvclev") )%dim_id
       CALL IO_DEF_VAR (IO_file_id, "vct_a", NF_DOUBLE, 1, IO_dims, io_vcta_id)
       CALL IO_DEF_VAR (IO_file_id, "vct_b", NF_DOUBLE, 1, IO_dims, io_vctb_id)
       IO_dims(1) = IO_dim_ids( IO_get_varindx("nswitches") )%dim_id
       CALL IO_DEF_VAR (IO_file_id, "switches", NF_DOUBLE, 1, IO_dims, io_nswitches_id)


       DO WHILE (ASSOCIATED(io_list))
          IF (io_list%field%info%restart) THEN
             ndim = io_list%field%info%ndim
             io_name = io_list%field%info%name
             CALL IO_DEF_VAR (IO_file_id, io_name, NF_DOUBLE, ndim,   &
                              io_list%field%info%IO_var_indx, &
                              io_list%field%info%IO_var_id)
          END IF
          io_list => io_list%next_list_element
       END DO

       CALL IO_ENDDEF(IO_file_id)

       CALL IO_PUT_VAR_DOUBLE (IO_file_id, io_vlat_id, vlat)
       CALL IO_PUT_VAR_DOUBLE (IO_file_id, io_vlon_id, vlon)
       CALL IO_PUT_VAR_DOUBLE (IO_file_id, io_vcta_id, vct(1:nvclev))
       CALL IO_PUT_VAR_DOUBLE (IO_file_id, io_vctb_id, vct(nvclev+1:2*nvclev))

       nswitches = 53
       ALLOCATE (restart_diag(nswitches))

!!!!  setting of "restart_diag" removed  
       restart_diag(:) = 0
!!!!  remove "restart_diag" later

       CALL IO_PUT_VAR_DOUBLE (IO_file_id, io_nswitches_id, restart_diag)
       DEALLOCATE (restart_diag)

    END IF

  END SUBROUTINE IO_write_header

  SUBROUTINE IO_read_header()

    USE mo_doctor,    ONLY: ylabel1, ylabel2, ylabel3, ylabel4, &
                            ylabel5, ylabel6, ylabel7, ylabel8


    IF (p_pe == p_io) THEN
       IO_file_id = header%nc_file_id
       CALL IO_GET_ATT_TEXT (IO_file_id, NF_GLOBAL, 'file_type', &
                             header%nc_file_type)

       IF ( header%nc_file_type(1:3) /= "Ini" .AND. &
            header%nc_file_type(1:3) /= "Res" ) THEN
          CALL message ('IO_read_header', header%nc_file_type);
          CALL message ('IO_read_header', 'No ECHAM initial or restart file.');
          CALL finish  ('IO_read_header', 'Run terminated.')
       END IF

       !    CALL IO_GET_ATT_TEXT (IO_file_id, NF_GLOBAL, 'title', header%nc_title)
       CALL IO_GET_ATT_TEXT (IO_file_id, NF_GLOBAL, 'file_type', header%nc_file_type)
       CALL IO_GET_ATT_TEXT (IO_file_id, NF_GLOBAL, 'source_type', header%nc_binary_source)
       CALL IO_GET_ATT_TEXT (IO_file_id, NF_GLOBAL, 'history', header%nc_creation_program)
       CALL IO_GET_ATT_TEXT (IO_file_id, NF_GLOBAL, 'user', header%nc_creation_user)
       CALL IO_GET_ATT_TEXT (IO_file_id, NF_GLOBAL, 'created', header%nc_creation_date)

       ylabel1(:) = ' '; ylabel2(:) = ' '; ylabel3(:) = ' '; ylabel4(:) = ' ';
       ylabel5(:) = ' '; ylabel6(:) = ' '; ylabel7(:) = ' '; ylabel8(:) = ' ';

       CALL IO_GET_ATT_TEXT (IO_file_id, NF_GLOBAL, 'label_1', ylabel1)
       CALL IO_GET_ATT_TEXT (IO_file_id, NF_GLOBAL, 'label_2', ylabel2)
       CALL IO_GET_ATT_TEXT (IO_file_id, NF_GLOBAL, 'label_3', ylabel3)
       CALL IO_GET_ATT_TEXT (IO_file_id, NF_GLOBAL, 'label_4', ylabel4)
       CALL IO_GET_ATT_TEXT (IO_file_id, NF_GLOBAL, 'label_5', ylabel5)
       CALL IO_GET_ATT_TEXT (IO_file_id, NF_GLOBAL, 'label_6', ylabel6)
       CALL IO_GET_ATT_TEXT (IO_file_id, NF_GLOBAL, 'label_7', ylabel7)
       CALL IO_GET_ATT_TEXT (IO_file_id, NF_GLOBAL, 'label_8', ylabel8)

       IF (ldebugio) THEN
          WRITE (nerr, *)
          IF (header%nc_file_type(1:12) == "Initial file") THEN
             WRITE (nerr, '(6(1x,a,/))') ylabel1, ylabel2, ylabel3, ylabel4, &
                                         ylabel5, ylabel6
          ELSE
             WRITE (nerr, '(8(1x,a,/))') ylabel1, ylabel2, ylabel3, ylabel4, &
                                         ylabel5, ylabel6, ylabel7, ylabel8 
          END IF
          WRITE (nerr, *)
       END IF
    END IF

    CALL p_bcast (ylabel1, p_io)
    CALL p_bcast (ylabel2, p_io)
    CALL p_bcast (ylabel3, p_io)
    CALL p_bcast (ylabel4, p_io)
    CALL p_bcast (ylabel5, p_io)
    CALL p_bcast (ylabel6, p_io)
    CALL p_bcast (ylabel7, p_io)
    CALL p_bcast (ylabel8, p_io)

  END SUBROUTINE IO_read_header

  SUBROUTINE IO_init

    USE mo_control,       ONLY: ngl, nhgl, nlon, nlp2, nlev, nlevp1, nsp, nvclev, nmp1, &
                                ntbase, ncbase, vct, nm, nn, nk, ntimst, ntdata, ncdata, &
                                n2sp, n4mp1, n2mp1, nnp1, nkp1, maxrow
    USE mo_start_dataset, ONLY: lres, nhf1, nisp, nstep, ly365, ndstart, ntstart
    USE mo_year,          ONLY: iymd2c, cd2dat, ihms2sec
!   USE m_alloc_mods,     ONLY: alloc_mods ! module subroutine
    USE mo_tracer,        ONLY: nhtrac

    REAL, ALLOCATABLE :: restart_diag(:)

    INTEGER :: idv
    INTEGER :: nswitches
    INTEGER :: icd, id, im, iy, ih, imin

    CHARACTER (3) :: month(12) = (/ 'jan', 'feb', 'mar', 'apr', 'may', 'jun', &
                                    'jul', 'aug', 'sep', 'oct', 'nov', 'dec' /)


    yearnc%nc_opened = .FALSE.
    sstnc0%nc_opened = .FALSE.
    sstnc1%nc_opened = .FALSE.
    sstnc2%nc_opened = .FALSE.
    sicnc0%nc_opened = .FALSE.
    sicnc1%nc_opened = .FALSE.
    sicnc2%nc_opened = .FALSE.
    header%nc_opened = .FALSE.
    ini_surf%nc_opened = .FALSE.
    ini_spec%nc_opened = .FALSE.
    ini_ozon%nc_opened = .FALSE.
    ini_field%nc_opened = .FALSE.
    restart(31:38)%nc_opened = .FALSE.

    IF (lres) THEN
       idv = nhf1
    ELSE
       idv = nisp
    END IF

    CALL IO_open_unit(idv, header, IO_READ)

    CALL IO_read_header()

    IF (p_pe == p_io) THEN
       IO_file_id = header%nc_file_id

       ! get data reference times

       CALL IO_GET_ATT_INT (IO_file_id, NF_GLOBAL, 'fdate', forecast_date)
       CALL IO_GET_ATT_INT (IO_file_id, NF_GLOBAL, 'ftime', forecast_time)
       CALL IO_GET_ATT_INT (IO_file_id, NF_GLOBAL, 'vdate', verification_date)
       CALL IO_GET_ATT_INT (IO_file_id, NF_GLOBAL, 'vtime', verification_time)
    END IF
    CALL p_bcast (forecast_date, p_io)
    CALL p_bcast (forecast_time, p_io)
    CALL p_bcast (verification_date, p_io)
    CALL p_bcast (verification_time, p_io)


    IF (ndstart == 1) THEN
      CALL reset_year
    ELSE IF (ndstart == 0) THEN
      CONTINUE
    ELSEIF(ndstart >= 101) THEN
      !  Set initial date and time to namelist value
      forecast_date = ndstart    !  YYYYMMDD
      forecast_time = ntstart    !    HHMMSS
      WRITE(nout,'(A,I8.8)') 'Initial date set to ', ndstart
      WRITE(nout,'(A,I8.6)') 'Initial time set to ', ntstart
    ELSE
      WRITE(nout,*) 'Improper initial date', ndstart
      CALL finish('IO_init', 'Improper initial date')
    ENDIF
    
    icd = iymd2c(forecast_date)
    CALL cd2dat (icd, id, im, iy)

    ih   = forecast_time/10000 
    imin = (forecast_time -ih*10000)/100 

    IF (p_pe == p_io) THEN

       WRITE(nout,'(a,i2.2,a1,i2.2,a,i2,1x,a3,i5,a,i15)') &
                  ' Initial data is at    ', ih, ':', imin, &
                  ' on ', id, month(im), iy, ' - century/julian day: ', icd

       ! get spherical truncations ...

       CALL IO_GET_ATT_INT (IO_file_id, NF_GLOBAL, 'spherical_truncation_n', nn)
       CALL IO_GET_ATT_INT (IO_file_id, NF_GLOBAL, 'spherical_truncation_m', nm)
       CALL IO_GET_ATT_INT (IO_file_id, NF_GLOBAL, 'spherical_truncation_k', nk)

       ! get nstep

       IF (lres) CALL IO_GET_ATT_INT (IO_file_id, NF_GLOBAL, 'nstep', nstep)

       ! get timestep

       IF (lres) CALL IO_GET_ATT_INT (IO_file_id, NF_GLOBAL, 'timestep', IO_timestep)

       ! inquire for dimensions and get values

       CALL IO_INQ_DIMID  (IO_file_id, 'ngl', IO_dim_id)
       CALL IO_INQ_DIMLEN (IO_file_id, IO_dim_id, ngl)
       !    IF (lres) CALL IO_INQ_DIMID  (IO_file_id, 'nhgl', IO_dim_id)
       !    IF (lres) CALL IO_INQ_DIMLEN (IO_file_id, IO_dim_id, nhgl)
       CALL IO_INQ_DIMID  (IO_file_id, 'nlon', IO_dim_id)
       CALL IO_INQ_DIMLEN (IO_file_id, IO_dim_id, nlon)
       !    CALL IO_INQ_DIMID  (IO_file_id, 'nlp2', IO_dim_id)
       !    CALL IO_INQ_DIMLEN (IO_file_id, IO_dim_id, nlp2)
       CALL IO_INQ_DIMID  (IO_file_id, 'nlev', IO_dim_id)
       CALL IO_INQ_DIMLEN (IO_file_id, IO_dim_id, nlev)
       !    CALL IO_INQ_DIMID  (IO_file_id, 'nlevp1', IO_dim_id)
       !    CALL IO_INQ_DIMLEN (IO_file_id, IO_dim_id, nlevp1)
       CALL IO_INQ_DIMID  (IO_file_id, 'nsp', IO_dim_id)
       CALL IO_INQ_DIMLEN (IO_file_id, IO_dim_id, nsp)
       CALL IO_INQ_DIMID  (IO_file_id, 'nvclev', IO_dim_id)
       CALL IO_INQ_DIMLEN (IO_file_id, IO_dim_id, nvclev)
       !    CALL IO_INQ_DIMID  (IO_file_id, 'n2', IO_dim_id)
       !    CALL IO_INQ_DIMLEN (IO_file_id, IO_dim_id, n2)
       IF (lres) CALL IO_INQ_DIMID  (IO_file_id, 'nhtrac', IO_dim_id)
       IF (lres) CALL IO_INQ_DIMLEN (IO_file_id, IO_dim_id, nhtrac)

    END IF
    CALL p_bcast (nn, p_io)
    CALL p_bcast (nm, p_io)
    CALL p_bcast (nk, p_io)
    CALL p_bcast (nstep, p_io)
    CALL p_bcast (IO_timestep, p_io)
    CALL p_bcast (ngl, p_io)
    CALL p_bcast (nlon, p_io)
    CALL p_bcast (nlev, p_io)
    CALL p_bcast (nsp, p_io)
    CALL p_bcast (nvclev, p_io)  
    CALL p_bcast (nhtrac, p_io)

    ! derive dependend dimensions 

    maxrow = ngl
    nkp1   = nk+1
    nmp1   = nm+1
    nnp1   = nn+1
    n2mp1  = nmp1+nmp1
    n4mp1  = n2mp1+n2mp1
    nlevp1 = nlev+1
    nhgl   = ngl/2
    nlp2   = nlon+2
    n2sp   = nsp+nsp

    ncbase = iymd2c(forecast_date)
    ntbase = ihms2sec(forecast_time)
    ncdata = iymd2c(verification_date)
    ntdata = ihms2sec(verification_time)
    ntimst = 1

    ! read lon

    ALLOCATE (vlon(nlon))
    IF (p_pe == p_io) THEN
       CALL IO_INQ_VARID (IO_file_id, 'lon', IO_var_id)
       CALL IO_GET_VAR_DOUBLE (IO_file_id, IO_var_id, vlon)
    END IF
    CALL p_bcast (vlon, p_io)

    ! read lat

    ALLOCATE (vlat(ngl))
    IF (p_pe == p_io) THEN
       CALL IO_INQ_VARID (IO_file_id, 'lat', IO_var_id)
       CALL IO_GET_VAR_DOUBLE (IO_file_id, IO_var_id, vlat)
    END IF
    CALL p_bcast (vlat, p_io)

    ! read vct 

    ALLOCATE (vct(nvclev*2))
    IF (p_pe == p_io) THEN
       CALL IO_INQ_VARID (IO_file_id, 'vct_a', IO_var_id)
       CALL IO_GET_VAR_DOUBLE (IO_file_id, IO_var_id, vct(1:nvclev))

       CALL IO_INQ_VARID (IO_file_id, 'vct_b', IO_var_id)   
       CALL IO_GET_VAR_DOUBLE (IO_file_id, IO_var_id, vct(nvclev+1:2*nvclev))
    END IF
    CALL p_bcast (vct, p_io)

    IF (lres) THEN

       IF (p_pe == p_io) THEN
          CALL IO_INQ_DIMID  (IO_file_id, 'nswitches', IO_dim_id)
          CALL IO_INQ_DIMLEN (IO_file_id, IO_dim_id, nswitches)
       END IF
       CALL p_bcast (nswitches, p_io)

       ALLOCATE (restart_diag(nswitches))
       IF (p_pe == p_io) THEN
          CALL IO_INQ_VARID  (IO_file_id, 'switches', IO_var_id)
          CALL IO_GET_VAR_DOUBLE (IO_file_id, IO_var_id, restart_diag)
       END IF
       CALL p_bcast (restart_diag, p_io)

!!!!  setting of physical budgets removed  
!!!!  remove "restart_diag" later

       icd   = iymd2c(verification_date)
       CALL cd2dat(icd,id,im,iy)

       ih   = verification_time/10000
       imin = (verification_time -ih*10000)/100 

       WRITE(nout,'(a,i2.2,a1,i2.2,a,i2,1x,a3,i5,a,i15)') &
                  ' Experiment resumed at ', ih, ':', imin, &
                  ' on ', id, month(im), iy, ' - century/julian day: ', icd

       DEALLOCATE (restart_diag)
    END IF

    CALL IO_close(header)

    ! Allocate arrays which depend on nlev

!    CALL alloc_mods  postponed to inictl

    CALL IO_init_dims

  END SUBROUTINE IO_init

  SUBROUTINE reset_year

    IF (header%nc_file_type(1:12) == "Initial file") THEN
       forecast_date     = forecast_date-10000*(forecast_date/10000-1)
       verification_date = verification_date-10000*(verification_date/10000-1)
    ELSE
       forecast_date = forecast_date-10000*(forecast_date/10000-1)   
    END IF

  END SUBROUTINE reset_year

  INTEGER FUNCTION IO_xl (nxl)

    USE mo_start_dataset, ONLY: nhg3

    INTEGER, INTENT(IN) :: nxl
    INTEGER :: io_unit
    CHARACTER (8) :: yname


    io_unit = nhg3

    WRITE(yname,'(a5,i2.2)' ) 'ng3xl', nxl

    IF (p_pe == p_io) THEN
       CALL IO_open_unit(io_unit, restart(io_unit), IO_READ)
       IO_file_id = restart(io_unit)%nc_file_id
       CALL IO_INQ_DIMID  (IO_file_id, yname, IO_dim_id)
       CALL IO_INQ_DIMLEN (IO_file_id, IO_dim_id, IO_xl)
       CALL IO_close(restart(io_unit))
    END IF
    CALL p_bcast (IO_xl, p_io)

  END FUNCTION IO_xl

  INTEGER FUNCTION IO_ng3xp ()

    USE mo_start_dataset, ONLY: nhg3

    INTEGER :: io_unit


    io_unit = nhg3
    IF (p_pe == p_io) THEN
       CALL IO_open_unit(io_unit, restart(io_unit), IO_READ)
       CALL IO_GET_ATT_INT (restart(io_unit)%nc_file_id, NF_GLOBAL, 'ng3xp', IO_ng3xp)
       CALL IO_close(restart(io_unit))
    END IF
    CALL p_bcast (IO_ng3xp, p_io)

  END FUNCTION IO_ng3xp

  REAL FUNCTION IO_dt ()

    IF (IO_timestep == -1) CALL finish('IO_dt','timestep was not read')
    IO_dt = IO_timestep

  END FUNCTION IO_dt

  SUBROUTINE IO_write_buffer(io_unit)

    USE mo_start_dataset, ONLY: nhg1, nhg2, nhg3, nhgl1, nhf1
    USE mo_memory_f,      ONLY: f
    USE mo_memory_gl,     ONLY: gl
    USE mo_memory_g1a,    ONLY: g1a
    USE mo_memory_g2a,    ONLY: g2a
    USE mo_memory_g3a,    ONLY: g3a
    USE mo_linked_list,   ONLY: list_element
    USE mo_decomposition, ONLY: dcg => global_decomposition
    USE mo_transpose,     ONLY: gather_gp, gather_sa

    INTEGER, INTENT(IN) :: io_unit
    REAL, POINTER :: zout(:,:,:,:), zptr(:,:,:,:)
    TYPE (list_element), POINTER :: io_list


    IF (io_unit == nhf1) THEN
       io_list => f%first_list_element
    ELSE IF (io_unit == nhgl1) THEN
       io_list => gl%first_list_element
    ELSE IF (io_unit == nhg1) THEN
       io_list => g1a%first_list_element
    ELSE IF (io_unit == nhg2) THEN
       io_list => g2a%first_list_element
    ELSE IF (io_unit == nhg3) THEN
       io_list => g3a%first_list_element
    ELSE
       WRITE(nerr, *) 'IO_write_buffer: io_unit=',io_unit
       CALL finish ('IO_write_buffer', 'io_unit unexpected')
    END IF

    ! Open restart file
    
    CALL IO_open_unit(io_unit, restart(io_unit), IO_WRITE)

    ! Write data description record
    
    CALL IO_write_header(io_unit, io_list)

    ! Necessary, because IO_write_header changes io_list

    IF (io_unit == nhf1) THEN
       io_list => f%first_list_element
    ELSE IF (io_unit == nhgl1) THEN
       io_list => gl%first_list_element
    ELSE IF (io_unit == nhg1) THEN
       io_list => g1a%first_list_element
    ELSE IF (io_unit == nhg2) THEN
       io_list => g2a%first_list_element
    ELSE IF (io_unit == nhg3) THEN
       io_list => g3a%first_list_element
    END IF

    ! Write buffer

    IF (p_pe == p_io) THEN
       IO_file_id = restart(io_unit)%nc_file_id
    END IF

    DO WHILE (ASSOCIATED(io_list))
       IF (io_list%field%info%restart) THEN
          ALLOCATE(zout(io_list%field%info%gdim_1, &
                        io_list%field%info%gdim_2, &
                        io_list%field%info%gdim_3, & 
                        io_list%field%info%gdim_4))
          zptr => io_list%field%ptr(:,:,:,:)
          IF (io_unit == nhf1) THEN
             CALL gather_sa (zout, zptr, dcg)
          ELSE   
             CALL gather_gp (zout, zptr, dcg)
          END IF   
          IF (p_pe == p_io) THEN
             IO_var_id = io_list%field%info%IO_var_id
             CALL IO_PUT_VAR_DOUBLE (IO_file_id, IO_var_id, zout)
          END IF
          DEALLOCATE(zout)
       END IF

       io_list => io_list%next_list_element
    END DO
    
    ! Close restart file

    CALL IO_close(restart(io_unit))

  END SUBROUTINE IO_write_buffer

  SUBROUTINE IO_read_buffer (io_unit)

    USE mo_start_dataset, ONLY: nhf1, nhg1, nhg2, nhg3, nhgl1, lg3setzero, lg3force
    USE mo_memory_f,      ONLY: f
    USE mo_memory_gl,     ONLY: gl
    USE mo_memory_g1a,    ONLY: g1a
    USE mo_memory_g2a,    ONLY: g2a
    USE mo_memory_g3a,    ONLY: g3a
    USE mo_linked_list,   ONLY: list_element
    USE mo_mpi,           ONLY: p_pe, p_io
    USE mo_decomposition, ONLY: dcg => global_decomposition
    USE mo_transpose,     ONLY: scatter_gp, scatter_sa
    USE mo_tracer,        ONLY: ntrac, nhtrac

    INTEGER, INTENT(IN) :: io_unit
    REAL, POINTER :: zin(:,:,:,:), zptr(:,:,:,:)
    REAL, POINTER :: zvarp(:,:)
    TYPE (list_element), POINTER :: io_list
    LOGICAL :: lhtrac
    INTEGER :: start(4), COUNT(4)
    INTEGER :: ntracdim
    INTEGER :: itrac
    INTEGER       :: i, status
    LOGICAL       :: lsetzero
    CHARACTER (8) :: yname
    LOGICAL       :: lrvarp


    lsetzero = .FALSE.
    lrvarp   = .FALSE.

    IF (io_unit == nhf1) THEN
       io_list => f%first_list_element
    ELSE IF (io_unit == nhgl1) THEN
       io_list => gl%first_list_element
    ELSE IF (io_unit == nhg1) THEN
       io_list => g1a%first_list_element
    ELSE IF (io_unit == nhg2) THEN
       io_list => g2a%first_list_element
    ELSE IF (io_unit == nhg3) THEN
       io_list => g3a%first_list_element
    ELSE
       WRITE(nerr, *) 'IO_read_buffer: io_unit=',io_unit
       CALL finish ('IO_read_buffer', 'io_unit unexpected')
    END IF

    ! Open restart file

    CALL IO_open_unit(io_unit, restart(io_unit), IO_READ)

    ! Read buffer

    IF (p_pe == p_io) THEN
       IO_file_id = restart(io_unit)%nc_file_id
    END IF

    DO WHILE (ASSOCIATED(io_list))
       IF (io_list%field%info%restart) THEN
          ALLOCATE(zin(io_list%field%info%gdim_1, &
                       io_list%field%info%gdim_2, &
                       io_list%field%info%gdim_3, &
                       io_list%field%info%gdim_4))
          IF (p_pe == p_io) THEN
             IF (io_unit == nhg3 .AND. lg3force) THEN
                status = NF_INQ_VARID (IO_file_id, io_list%field%info%name, IO_var_id)
                IF (status /= NF_NOERR) THEN
                   yname = io_list%field%info%name
                   DO i = 1, 8
                      IF (yname(i:i) == ' ') EXIT
                   END DO
                   yname(i-1:i-1) = ''
                   status = NF_INQ_VARID (IO_file_id, yname, IO_var_id)
                   IF (status /= NF_NOERR) THEN
                      IF (lg3setzero) THEN
                         lsetzero = .TRUE.
                         WRITE(nerr,*) 'IO_read_buffer : ', io_list%field%info%name(1:8), &
                              ' not found ... set to zero !!!'
                      ELSE
                         ! Only for error handling
                         CALL IO_INQ_VARID (IO_file_id, io_list%field%info%name, IO_var_id)
                      END IF
                   ELSE
                      WRITE(nerr,*) 'IO_read_buffer : ', io_list%field%info%name(1:8), &
                           ' not found ... read ', yname, ' instead !!!'
                   END IF
                END IF
             ELSE
                status = NF_INQ_VARID (IO_file_id, io_list%field%info%name, IO_var_id)
                IF (io_unit == nhg3 .AND. status /= NF_NOERR) THEN
                   yname = io_list%field%info%name
                   lrvarp = .FALSE.
                   IF (yname(1:4)=='EWOV' .OR. yname(1:4)=='NSOV' .OR.  &
                       yname(1:4)=='NWOV' .OR. yname(1:4)=='NEOV') THEN
                      ! If unpacked directional orographic variance not found
                      ! try to read packed directional orographic variance
                      ! and unpack it (not on CRAY's !)
                      status = NF_INQ_VARID (IO_file_id, 'VARPM', IO_var_id)
                      IF (status /= NF_NOERR) THEN
                         ! Only for error handling
                         CALL IO_INQ_VARID (IO_file_id, 'VARPM', IO_var_id)
                      ELSE
#ifdef CRAY
                         CALL finish  ('IO_read_buffer', 'This restart file gives wrong results on CRAY')
#endif
                         WRITE(nout,*) 'IO_read_buffer: use unpacked ',yname(1:4),' from VARP'
                         lrvarp = .TRUE.
                      END IF
                   END IF
                ELSE
                   ! Only for error handling
                   CALL IO_INQ_VARID (IO_file_id, io_list%field%info%name, IO_var_id)
                END IF
             END IF

             lhtrac = .FALSE.
             IF (ntrac > 0 .AND. ntrac /= nhtrac) THEN
                IF (io_list%field%info%name == 'XT' .OR. &
                    io_list%field%info%name == 'XTM1') lhtrac = .TRUE.
             END IF
             IF (lhtrac) THEN
                WRITE(*,*) 'lhtrac: ',io_list%field%info%name
                IF (ntrac > nhtrac) THEN
                   ntracdim = nhtrac
                ELSE
                   ntracdim = io_list%field%info%gdim_3
                END IF
                COUNT(:) = (/ io_list%field%info%gdim_1, &
                              io_list%field%info%gdim_2, &
                              1,                  &
                              io_list%field%info%gdim_4 /)
                DO itrac = 1, ntracdim
                   start(:) = (/ 1, 1, itrac, 1 /)
                   CALL IO_GET_VARA_DOUBLE (IO_file_id, IO_var_id, start, count, zin(:,:,itrac,:))
                ENDDO
             ELSE
                IF (io_unit == nhg3 .AND. lsetzero) THEN
                   zin = 0.
                   lsetzero = .FALSE.
                ELSE
                   IF (lrvarp) THEN
                      lrvarp = .FALSE.
                      ALLOCATE(zvarp(io_list%field%info%gdim_1,io_list%field%info%gdim_2))
                      CALL IO_GET_VAR_DOUBLE (IO_file_id, IO_var_id, zvarp)

                      CALL unpack_varp(io_list%field%info%gdim_1*io_list%field%info%gdim_2, &
                                  io_list%field%info%name, zvarp, zin)

                      DEALLOCATE (zvarp)
                   ELSE
                      CALL IO_GET_VAR_DOUBLE (IO_file_id, IO_var_id, zin)
                   END IF
                END IF
             END IF
          END IF
          zptr => io_list%field%ptr(:,:,:,:)
          IF (io_unit == nhf1) THEN
             CALL scatter_sa (zin, zptr, dcg)
          ELSE   
             CALL scatter_gp (zin, zptr, dcg)
          END IF
          DEALLOCATE(zin)
       END IF
       io_list => io_list%next_list_element
    END DO

    ! Close restart file

    CALL IO_close(restart(io_unit))

    RETURN

  CONTAINS

    SUBROUTINE unpack_varp(ndim, cname, zvarp, zin)

      USE mo_exception,       ONLY: finish

      CHARACTER (len=4) :: cname
      INTEGER :: ndim, nrec
      REAL    :: zvarp(ndim), zin(ndim)
      REAL    :: ztmp(ndim,4)

      EXTERNAL util_gwunpk

      CALL util_gwunpk(ztmp(1,3),ztmp(1,1),ztmp(1,2),ztmp(1,4),zvarp,ndim)

      IF      (cname == 'EWOV') THEN
         nrec = 3
      ELSE IF (cname == 'NSOV') THEN
         nrec = 1
      ELSE IF (cname == 'NWOV') THEN
         nrec = 4
      ELSE IF (cname == 'NEOV') THEN
         nrec = 2
      ELSE
         CALL finish ('unpack_varp', cname)
      END IF

      zin(:) = ztmp(:,nrec)

    END SUBROUTINE unpack_varp

  END SUBROUTINE IO_read_buffer

END MODULE mo_io
