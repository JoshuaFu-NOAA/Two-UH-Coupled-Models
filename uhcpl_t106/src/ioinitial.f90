SUBROUTINE ioinitial

  ! Description:
  !
  ! Read initial data
  !
  ! Method:
  !
  ! An abstracttion layer is used to access netCDF data files and
  ! retrieve data.   
  ! 
  ! Further information is contained in a note on the IO package
  !
  ! Authors:
  !
  ! U. Schulzweida, MPI, May 1999, original version
  !
  !

  USE mo_io,            ONLY: ini_spec, ini_surf, io_file_id, io_read,    &
                              io_var_id, io_get_var_double, io_open_unit, &
                              io_inq_varid, io_close
  USE mo_netCDF,        ONLY: NF_INQ_VARID, NF_NOERR
  USE mo_mpi,           ONLY: p_io, p_pe
  USE mo_doctor,        ONLY: nerr, nout
  USE mo_control,       ONLY: nlevp1, lvctch, lamip2
  USE mo_start_dataset, ONLY: nigp, nisp, ldebugio
  USE mo_memory_sp,     ONLY: sd, sp, stp, su0, svo, memory_info, get_info
  USE mo_memory_gl,     ONLY: gl, q, x, xt
  USE mo_memory_g3a,    ONLY: g3a, get_entry
  USE mo_tracer,        ONLY: ntrac, xtini
  USE mo_hyb,           ONLY: apsurf
  USE mo_decomposition, ONLY: global_decomposition
  USE mo_transpose,     ONLY: scatter_sp, scatter_gp
#ifdef CRAY
  USE mo_exception,     ONLY: finish
#endif

  IMPLICIT NONE

  !  Local scalars: 

  ! number of codes read from surface initialfile
  INTEGER, PARAMETER :: nrec_surf = 18
  CHARACTER (8) :: cname, csurf(nrec_surf)

  INTEGER :: nsvoid, nsdid, nstpid, nqid
  INTEGER :: irec, i
  INTEGER :: status

  REAL, POINTER :: zin(:,:,:), zsu0(:,:), zptr(:,:,:)
  REAL, POINTER :: zvarp(:,:,:)

  TYPE (memory_info) :: info

  !  Intrinsic functions 
  INTRINSIC LOG


  !  Executable Statements

  ! Initial file information already read in initialise (CALL IO_init)

  ! 1. Process spectral files

  ! skip if column model runs with changed hybrid levels

  IF (.NOT.lvctch) THEN

    IF (p_pe == p_io) THEN

      CALL IO_open_unit(nisp, ini_spec, IO_READ)

      IO_file_id = ini_spec%nc_file_id

      CALL IO_INQ_VARID (IO_file_id, 'SVO', nsvoid)
      CALL get_info (sp, 'SVO', info)
      ALLOCATE (zin(info%gdim_1, info%gdim_2, info%gdim_3))
      CALL IO_GET_VAR_DOUBLE (IO_file_id, nsvoid, zin)

    END IF

    CALL scatter_sp (zin, svo, global_decomposition)

    ! 1.1 derive su0 from svo, saves one transpose operation and is 
    !     fast enough

    IF (p_pe == p_io) THEN

      CALL get_info (sp, 'SU0', info)
      ALLOCATE (zsu0(info%gdim_1, info%gdim_2))
      CALL init_su0 (zin, zsu0)      

    END IF

    CALL scatter_sp (zsu0, su0, global_decomposition)

    ! finish setup of svo and calculation of su0  

    IF (p_pe == p_io) DEALLOCATE (zin, zsu0)

    IF (p_pe == p_io) THEN

      CALL IO_INQ_VARID (IO_file_id, 'SD', nsdid)
      CALL get_info (sp, 'SD', info)
      ALLOCATE (zin(info%gdim_1, info%gdim_2, info%gdim_3))
      CALL IO_GET_VAR_DOUBLE (IO_file_id, nsdid, zin)

    END IF

    CALL scatter_sp (zin, sd, global_decomposition)

    IF (p_pe == p_io) DEALLOCATE (zin)

    IF (p_pe == p_io) THEN

      CALL IO_INQ_VARID (IO_file_id, 'STP', nstpid)
      CALL get_info (sp, 'STP', info)
      ALLOCATE (zin(info%gdim_1, info%gdim_2, info%gdim_3))
      CALL IO_GET_VAR_DOUBLE (IO_file_id, nstpid, zin)

      ! Set global mean of surface pressure for initial field: STP(nlevp1,1,1)
     
      IF (lamip2) THEN
        zin(nlevp1,1,1) = LOG(apsurf-286.)
      ELSE
        zin(nlevp1,1,1) = LOG(apsurf)
      END IF
     
    END IF

    CALL scatter_sp (zin, stp, global_decomposition)

    IF (p_pe == p_io) DEALLOCATE (zin)

    ! 2. Read grid point data, advected by SL

    IF (p_pe == p_io) THEN

      CALL IO_INQ_VARID (IO_file_id, 'Q', nqid)
      CALL get_info (gl, 'Q', info)
      ALLOCATE (zin(info%gdim_1, info%gdim_2, info%gdim_3))
      CALL IO_GET_VAR_DOUBLE (IO_file_id, nqid, zin)

    END IF

    IF (p_pe == p_io) CALL IO_close(ini_spec)

    CALL scatter_gp (zin, q, global_decomposition)

    IF (p_pe == p_io) DEALLOCATE (zin)

  ELSE          ! column model runs with changed hybrid levels
    q   = 0.
  ENDIF

  ! 2.1 liquid water set initial to zero

  x(:,:,:) = 0.

  ! 2.2 initialize optional tracer fields

  IF (ntrac > 0) CALL xtini(1, xt)

  ! 3. Prepare grid point surface fields

  IF (p_pe == p_io) THEN

    CALL IO_open_unit(nigp, ini_surf, IO_READ)
    IO_file_id = ini_surf%nc_file_id

  END IF

  ! Codes read from surface initialfile (unit:24)

  csurf( 1) = 'GEOSP'    ! Surface geopotential
  csurf( 2) = 'TS'       ! Surface temperature
  csurf( 3) = 'WS'       ! Surface soil wetness
  csurf( 4) = 'WL'       ! Skin reservoir content
  csurf( 5) = 'SN'       ! Snow depth
  csurf( 6) = 'SLM'      ! Land sea mask
  csurf( 7) = 'ALB'      ! Albedo
  csurf( 8) = 'AZ0'      ! Surface roughness
  csurf( 9) = 'EWOV'     !  E-W  orographic variance
  csurf(10) = 'NSOV'     !  N-S  orographic variance
  csurf(11) = 'NWOV'     ! NW-SE orographic variance
  csurf(12) = 'NEOV'     ! NE-SW orographic variance
  csurf(13) = 'VAROR'    ! Orographic variance (for surface runoff)
  csurf(14) = 'FOREST'   ! Vegetation type
  csurf(15) = 'VGRAT'    ! Vegetation ratio
  csurf(16) = 'VLT'      ! Leaf area index
  csurf(17) = 'WSMX'     ! Field capacity of soil
  csurf(18) = 'FAO'      ! Fao data set

  DO irec = 1, nrec_surf

    cname = csurf(irec)
    DO i=1, 8
      IF (cname(i:i).EQ.'') THEN
        cname(i:i) = 'M'
        EXIT
      END IF
    END DO

    IF (p_pe == p_io) THEN

      IF (ldebugio) WRITE(nerr,*) 'IO_initial : read ', csurf(irec)

      status = NF_INQ_VARID (IO_file_id, csurf(irec), IO_var_id)
      IF (status /= NF_NOERR) THEN
        IF (csurf(irec)(1:4)=='EWOV' .OR. csurf(irec)(1:4)=='NSOV' .OR.  &
            csurf(irec)(1:4)=='NWOV' .OR. csurf(irec)(1:4)=='NEOV') THEN
          ! If unpacked directional orographic variance not found
          ! try to read packed directional orographic variance
          ! and unpack it (not on CRAY's !)
          status = NF_INQ_VARID (IO_file_id, 'VARP', IO_var_id)
          IF (status /= NF_NOERR) THEN
            ! Only for error handling
            CALL IO_INQ_VARID (IO_file_id, csurf(irec), IO_var_id)
          ELSE
#ifdef CRAY
            CALL finish  ('ioinitial', 'This surface initial data gives wrong results on CRAY')
#endif
            WRITE(nout,*) 'ioinitial: use unpacked ',csurf(irec),' from VARP'
            CALL get_info (g3a, cname, info)
            ALLOCATE (zin(info%gdim_1, info%gdim_2, info%gdim_3))

            ALLOCATE (zvarp(info%gdim_1, info%gdim_2, info%gdim_3))

            CALL IO_GET_VAR_DOUBLE (IO_file_id, IO_var_id, zvarp)

            CALL unpack_varp(info%gdim_1*info%gdim_2*info%gdim_3, csurf(irec), zvarp, zin)

            DEALLOCATE (zvarp)

          END IF
        ELSE
          ! Only for error handling
          CALL IO_INQ_VARID (IO_file_id, csurf(irec), IO_var_id)
        END IF
      ELSE
        CALL get_info (g3a, cname, info)
        ALLOCATE (zin(info%gdim_1, info%gdim_2, info%gdim_3))

        CALL IO_GET_VAR_DOUBLE (IO_file_id, IO_var_id, zin)
      END IF

    END IF

    CALL get_entry (g3a, cname, zptr)

    CALL scatter_gp (zin, zptr, global_decomposition)

    IF (p_pe == p_io)  DEALLOCATE (zin)

  END DO

  IF (p_pe == p_io) CALL IO_close (ini_surf)

  CALL init_g3()

  WRITE (nout,'(/)')

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

  SUBROUTINE init_su0 (psvo, psu0)

    ! Description:
    !
    ! Compute initial spectral components for the zonal mean wind used
    ! in the linearization of the vorticity and humidity equations.
    !
    ! Method:
    !
    ! This subroutine computes initial spectral components
    ! of the mean zonal wind used in the semi-implicit treatment of
    ! vorticity and humidity equations from the vorticity zonal
    ! spectral components.
    !
    ! *inisu0* is called from *ioinitial*.
    !
    ! Authors:
    !
    ! M. Jarraud, ECMWF, February 1983, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_control,    ONLY: nlev, nn, nnp1
    USE mo_constants,  ONLY: a

    IMPLICIT NONE

    REAL, POINTER :: psvo(:,:,:), psu0(:,:)

    !  Local scalars: 
    REAL :: zeps1, zeps2, zn
    INTEGER :: jlev, jn

    !  Intrinsic functions 
    INTRINSIC SQRT


    !  Executable statements 

    !-- 1. Set up *su0*

    DO jlev = 1, nlev

      zeps2 = a/SQRT(3.)
      psu0(jlev,1) = zeps2*psvo(jlev,1,2)

      DO jn = 2, nn
        zeps1 = zeps2
        zn = 4.0*jn*jn-1.0
        zeps2 = a/SQRT(zn)
        psu0(jlev,jn) = -zeps1*psvo(jlev,1,jn-1)+zeps2*psvo(jlev,1,jn+1)
      END DO

      zeps1 = zeps2

      psu0(jlev,nnp1) = -zeps1*psvo(jlev,1,nn)
    END DO

  END SUBROUTINE init_su0

  SUBROUTINE init_g3 ()

    !
    ! init_g3a - initialize parameterisation scheme data.
    !
    ! J. K. Gibson, ECMWF, April 1983
    !
    ! Purpose: to prepare the *g3a* work buffer from an netCDF initial file.
    !
    ! Method: Initial values are set for appropriate variables.
    !

    USE mo_io_tables,     ONLY: ng3xp
    USE mo_post,          ONLY: lxaccu
    USE mo_physc2,        ONLY: ctfreez
    USE mo_memory_g3a,    ONLY: aprlm, aprcm, aprsm, sradsm, tradsm, srad0m, trad0m, &
                                vdism, ustrm, vstrm, ahfsm, evapm, ahflm, wind10m,   &
                                teffm, ustrgwm, vstrgwm, vdisgwm, temp2m, dew2m,     &
                                u10m, v10m, tsurfm, runoffm, srad0um, tradsum,       &
                                ustar3m, sradsum, t2maxm, t2minm, wimaxm, topmaxm,   &
                                snmelm, runtocm, tslinm, dsnacm, aclcvm, aclcovm,    &
                                qvim, alwcvim, auxil1m, auxil2m, tsnm, td3m, td4m,   &
                                td5m, tdm, sicedm, wlm1m, snm1m, emterm, trsolm,     &
                                aclcm, aclcacm, tkem, tkem1m, seaicem, tsm1m, tdm1m, &
                                td3m1m, td4m1m, td5m1m, tdclm1m, tsnm1m, g3m, glacm, &
                                faom, rgcgnm, sodifm, wsm, wsm1m, slmm, wsmxm, snm,  &
                                tsm, tsminm, tsmaxm, tdclm, wlm
    USE mo_doctor,        ONLY: nerr
    USE mo_decomposition, ONLY: lc => local_decomposition

    IMPLICIT NONE

    !  Local scalars: 

    INTEGER :: jl, jg, jx


    !  Executable Statements

    ! Initialize *g3a* variables not read

    aprlm(:,:)     = 0.0
    aprcm(:,:)     = 0.0
    aprsm(:,:)     = 0.0
    sradsm(:,:)    = 0.0
    tradsm(:,:)    = 0.0
    srad0m(:,:)    = 0.0
    trad0m(:,:)    = 0.0
    vdism(:,:)     = 0.0
    ustrm(:,:)     = 0.0
    vstrm(:,:)     = 0.0
    ahfsm(:,:)     = 0.0
    evapm(:,:)     = 0.0
    ahflm(:,:)     = 0.0
    wind10m(:,:)   = 0.0
    teffm(:,:)     = 0.0
    ustrgwm(:,:)   = 0.0
    vstrgwm(:,:)   = 0.0
    vdisgwm(:,:)   = 0.0
    temp2m(:,:)    = 0.0
    dew2m(:,:)     = 0.0
    u10m(:,:)      = 0.0
    v10m(:,:)      = 0.0
    tsurfm(:,:)    = 0.0
    runoffm(:,:)   = 0.0
    srad0um(:,:)   = 0.0
    tradsum(:,:)   = 0.0
    ustar3m(:,:)   = 0.0
    sradsum(:,:)   = 0.0
    t2maxm(:,:)    = 0.0
    t2minm(:,:)    = 999.0
    tsmaxm(:,:)    = 0.0
    tsminm(:,:)    = 999.0
    wimaxm(:,:)    = 0.0
    topmaxm(:,:)   = 99999.0
    snmelm(:,:)    = 0.0
    runtocm(:,:)   = 0.0
    tslinm(:,:)    = 0.0
    dsnacm(:,:)    = 0.0
    aclcvm(:,:)    = 0.0
    aclcovm(:,:)   = 0.0
    qvim(:,:)      = 0.0
    alwcvim(:,:)   = 0.0
    auxil1m(:,:)   = tsm(:,:)
    auxil2m(:,:)   = tsm(:,:)
    tsnm(:,:)      = tsm(:,:)
    td3m(:,:)      = tsm(:,:)
    td4m(:,:)      = tdm(:,:)
    td5m(:,:)      = tdclm(:,:)
    tdm(:,:)       = tdclm(:,:)
    sicedm(:,:)    = 0.0
    wlm1m(:,:)     = wlm(:,:)
    snm1m(:,:)     = snm(:,:)

    emterm(:,:,:)  = 0.0
    trsolm(:,:,:)  = 0.0

    aclcm(:,:,:)   = 0.0
    aclcacm(:,:,:) = 0.0
    tkem(:,:,:)    = 1.0e-4
    tkem1m(:,:,:)  = tkem(:,:,:)

    ! Set soil temperatures

    DO jg = 1, lc%nglat
      DO jl = 1, lc%nglon
        IF(tsm(jl,jg) < ctfreez .AND. slmm(jl,jg) < 0.5) THEN
          seaicem(jl,jg) = 1.0
        ELSE
          seaicem(jl,jg) = 0.0
        ENDIF
      END DO
    END DO

    tsm1m(:,:)   = tsm(:,:)
    tdm1m(:,:)   = tdm(:,:)
    td3m1m(:,:)  = td3m(:,:)
    td4m1m(:,:)  = td4m(:,:)
    td5m1m(:,:)  = td5m(:,:)
    tdclm1m(:,:) = tdclm(:,:)
    tsnm1m(:,:)  = tsnm(:,:)

    ! Extra g3-fields

    IF(ng3xp.GT.0) THEN
      DO jx = 1,ng3xp
        IF(lxaccu(jx)) THEN
          g3m(jx)%x(:,:,:) = 0.0
        ENDIF
      ENDDO
    ENDIF

    ! Extra g4-fields

!!$  IF(ng4xp.GT.0) THEN
!!$     i1 = 1
!!$     DO jx = 1,ng4xp
!!$        i2 = nlon*ng4xl(jx)+i1-1
!!$        IF(l4xaccu(jx)) THEN
!!$           DO jl = i1,i2
!!$              !mem              g4x(jl,1,1) = 0.
!!$           ENDDO
!!$        ENDIF
!!$        i1 = i2+1
!!$     ENDDO
!!$  ENDIF

    ! Minor corrections to initial snow field. Corrections for 
    ! snow accumulation found in T21 and T42 in long term integrations
    ! seems to allow to define some places as glacier instead of simple
    ! snow coverage.
    ! Therefore the input fields from ECMWF are patched. snm > 9 defines
    ! a grid point as glacier.
    ! T21: snm(51:61,2) and snm1m((51:61,2)
    ! T42: WHERE slmm(101:120,3) > 0.5 => snm(101:120,3)
    !      snm(101:102:4) 
    !      snm(25,57) 
    !       
    ! NOTE: T30 must be set as well but is not defined yet.

    ! Glacier mask

    DO jg = 1, lc%nglat
      DO jl = 1, lc%nglon
        IF(snm(jl,jg) > 9.0 .AND. slmm(jl,jg) > 0.5) THEN
          glacm(jl,jg) = 1.
        ELSE
          glacm(jl,jg) = 0.
        ENDIF
      ENDDO
    ENDDO

    ! Setting of array of variable soil characteristics to be use
    ! in *surf*
    ! Input: FAO soils interpolated from 0.5 degree resolution
    !        to model resolution (simple average).
    ! and setting of array of variable available water storage capacity
    ! to be used in *surf*
    ! Input: Patterson data interpolated from 0.5 degree resolution
    !        to model resolution (simple average).

    DO jg = 1, lc%nglat
      DO jl = 1, lc%nglon
        IF (NINT(faom(jl,jg)) == 1) THEN
          rgcgnm(jl,jg) = 1.93e+06
          sodifm(jl,jg) = 8.7e-7
        ELSE IF (NINT(faom(jl,jg)) == 2) THEN
          rgcgnm(jl,jg) = 2.10e+06
          sodifm(jl,jg) = 8.0e-7
        ELSE IF (NINT(faom(jl,jg)) == 3) THEN
          rgcgnm(jl,jg) = 2.25e+06
          sodifm(jl,jg) = 7.4e-7
        ELSE IF (NINT(faom(jl,jg)) == 4) THEN
          rgcgnm(jl,jg) = 2.36e+06
          sodifm(jl,jg) = 7.1e-7
        ELSE IF (NINT(faom(jl,jg)) == 5) THEN
          rgcgnm(jl,jg) = 2.48e+06
          sodifm(jl,jg) = 6.7e-7
        ELSE
          IF (NINT(faom(jl,jg)) == 0) THEN
            rgcgnm(jl,jg) = 2.25e+06
            sodifm(jl,jg) = 7.4e-7
          ELSE
            WRITE (nerr,*) 'faom(',jl,',',jg,') = ',faom(jl,jg)
          END IF
        END IF
        wsm(jl,jg)   = MIN(wsm(jl,jg), wsmxm(jl,jg))
        wsm1m(jl,jg) = wsm(jl,jg)
      END DO
    END DO

  END SUBROUTINE init_g3

END SUBROUTINE ioinitial
