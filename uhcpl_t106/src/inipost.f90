SUBROUTINE inipost

  ! Description:
  !
  ! Preset constants in mo_control.
  !
  ! Method:
  !
  ! Preset module *mo_post*.
  !
  ! *inipost* is called from subroutine *initialize*
  !
  ! all variables of module *mo_post* are defined
  !
  ! Authors:
  !
  ! E. Kirk, MI, March 1989, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! L. Kornblueh, MPI, April 1998, added NWP forecast mode
  ! L. Kornblueh, MPI, June 1999, parallel version (MPI based)
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_exception
  USE mo_parameters
  USE mo_doctor
  USE mo_post
  USE mo_mpi,    ONLY: p_pe, p_io, p_parallel, p_parallel_io, p_bcast
  USE mo_control
  USE mo_start_dataset
  USE mo_hdiff
  USE mo_constants
  USE mo_gaussgrid
  USE mo_filename
  USE mo_year
  USE mo_grib, ONLY: init_grib, open_grib_file

  IMPLICIT NONE

  !  Local scalars:
  INTEGER :: i, icode, id, iday, ilast, im, inx, iret, iy, iypt, j, jfield, &
             jlev, name_length
  LOGICAL :: lexists
  CHARACTER (1) :: ycom
  CHARACTER (8) :: ynam, ystat
  CHARACTER (path_limit) :: ypph,ypps,yres  ! The latter two are added for AMIP2

  !  Local arrays: 
  INTEGER :: ilev(nlevp1)

  !  External subroutines 
  EXTERNAL pbopen

  INCLUDE 'postctl.inc'

  !  Executable statements 

  !-- 1. Define variables

  ! Logical units for grib output files

  CALL init_grib

  nunitsp = ngribs
  nunitun = ngribg
  nung4x = ngribx

  DO j = 1, 256
    yn(j) = ' '
    nlevg3(j) = 1
    setval(j) = 0.
    DO jlev = 1, nlevp1
      npplev(jlev,j) = 0
    END DO
  END DO

  lppspe = .TRUE.
  lppt   = .TRUE.
  lppd   = .TRUE.
  lppvo  = .TRUE.
  lppq   = .TRUE.
  lppp   = .TRUE.
  lppx   = .TRUE.

  DO j = 1, 256
    laccu(j) = .FALSE.
    npbits(j) = 16
  END DO

  IF (lamip2) THEN
    ! ----> increase number of pack bits for sensitive water budget variables
    npbits(142) = 24   ! aprl
    npbits(143) = 24   ! aprc
    npbits(144) = 24   ! aprs
    npbits(182) = 24   ! evap
    npbits(160) = 24   ! runoff
    npbits(221) = 24   ! dsnac
  END IF
  npbits(162) = 8   ! aclc
  npbits(164) = 8   ! aclcov
  npbits(172) = 8   ! slm
  npbits(174) = 8   ! alb
  npbits(175) = 8   ! albedo
  npbits(210) = 8   ! seaice
  npbits(212) = 8   ! forest
  npbits(163) = 8   ! aclcv
  npbits(223) = 8   ! aclcac

  laccu(142) = .TRUE.    ! aprl
  laccu(143) = .TRUE.    ! aprc
  laccu(144) = .TRUE.    ! aprs
  laccu(145) = .TRUE.    ! vdis
  laccu(146) = .TRUE.    ! ahfs
  laccu(147) = .TRUE.    ! ahfl

  laccu(159) = .TRUE.    ! ustar3
  laccu(160) = .TRUE.    ! runoff
  laccu(161) = .TRUE.    ! drain

  laccu(164) = .TRUE.    ! aclcov
  IF (.NOT. lamip2) THEN
    laccu(165) = .TRUE.  ! u10
    laccu(166) = .TRUE.  ! v10
    laccu(167) = .TRUE.  ! temp2
    laccu(168) = .TRUE.  ! dew2
  END IF
  laccu(169) = .TRUE.    ! tsurf
  laccu(171) = .TRUE.    ! wind10 !!

  laccu(176) = .TRUE.    ! srads
  laccu(177) = .TRUE.    ! trads
  laccu(178) = .TRUE.    ! srad0
  laccu(179) = .TRUE.    ! trad0
  laccu(180) = .TRUE.    ! ustr
  laccu(181) = .TRUE.    ! vstr
  laccu(182) = .TRUE.    ! evap

  laccu(185) = .TRUE.    ! srafs
  laccu(186) = .TRUE.    ! trafs
  laccu(187) = .TRUE.    ! sraf0
  laccu(188) = .TRUE.    ! traf0
  laccu(189) = .TRUE.    ! sclfs
  laccu(190) = .TRUE.    ! tclfs
  laccu(191) = .TRUE.    ! sclf0
  laccu(192) = .TRUE.    ! tclf0

  laccu(195) = .TRUE.    ! ustrgw
  laccu(196) = .TRUE.    ! vstrgw
  laccu(197) = .TRUE.    ! vdisgw

  laccu(203) = .TRUE.    ! srad0u
  laccu(204) = .TRUE.    ! sradsu
  laccu(205) = .TRUE.    ! tradsu

  laccu(213) = .TRUE.    ! teff

  laccu(218) = .TRUE.    ! snmel
  laccu(219) = .TRUE.    ! runtoc
  laccu(220) = .TRUE.    ! tslin
  laccu(221) = .TRUE.    ! dsnac
  laccu(223) = .TRUE.    ! aclcac
  laccu(230) = .TRUE.    ! qvi
  laccu(231) = .TRUE.    ! alwcvi

  yn(128+1) = 'GEOSP'         ! Surface geopotential
  yn(128+6) = 'APS'           ! Surface pressure
  yn(128+11) = 'TS'           ! Surface temperature
  yn(128+12) = 'WS'           ! Surface soil wetness
  yn(128+13) = 'SN'           ! Snow depth
  yn(128+14) = 'APRL'         ! Large scale rain
  yn(128+15) = 'APRC'         ! Convective rain
  yn(128+16) = 'APRS'         ! Snow fall
  yn(128+17) = 'VDIS'         ! Boundary layer dissipation
  yn(128+18) = 'AHFS'         ! Surface sensible heat flux
  yn(128+19) = 'AHFL'         ! Surface latent heat flux
  yn(128+36) = 'ACLCOV'       ! Cloud cover
  yn(128+37) = 'U10'          ! 10m u
  yn(128+38) = 'V10'          ! 10m v
  yn(128+39) = 'TEMP2'        ! 2m Temperature
  yn(128+40) = 'DEW2'         ! 2m dew point temperature
  yn(128+42) = 'TD'           ! Deep soil temperature
  yn(128+43) = 'WIND10'       ! 10M WINDSPEED
  yn(128+44) = 'SLM'          ! Land sea mask
  yn(128+45) = 'AZ0'          ! Surface roughness
  yn(128+46) = 'ALB'          ! Albedo
  yn(128+48) = 'SRADS'        ! Surface solar radiation
  yn(128+49) = 'TRADS'        ! Surface thermal radiation
  yn(128+50) = 'SRAD0'        ! Top solar radiation
  yn(128+51) = 'TRAD0'        ! Top thermal radiation
  yn(128+52) = 'USTR'         ! u-stress
  yn(128+53) = 'VSTR'         ! v-stress
  yn(128+54) = 'EVAP'         ! Evaporation
  yn(128+55) = 'TDCL'         ! Climatological deep soil temperature
  yn(159) = 'USTAR3'   ! ustar**3
  yn(160) = 'RUNOFF'   ! Surface runoff
  yn(161) = 'DRAIN'    ! Drainage
  yn(162) = 'ACLC'     ! Cloud cover
  yn(163) = 'ACLCV'    ! Total cloud cover (not acc.)

  yn(169) = 'TSURF'    ! Surface temperature

  yn(175) = 'ALBEDO'   ! Surface albedo

  yn(185) = 'SRAFS'    ! Surface solar radiation (clear sky)
  yn(186) = 'TRAFS'    ! Surface thermal radiation (clear sky)
  yn(187) = 'SRAF0'    ! Top solar radiation (clear sky)
  yn(188) = 'TRAF0'    ! Top thermal radiation (clear sky)
  yn(189) = 'SCLFS'    ! SW surf. cloud forcing
  yn(190) = 'TCLFS'    ! LW surf. cloud forcing
  yn(191) = 'SCLF0'    ! SW top cloud forcing
  yn(192) = 'TCLF0'    ! LW top cloud forcing
  yn(193) = 'WL'       ! Skin reservoir content
  yn(194) = 'WLM1'     ! Skin reservoir content at t-dt
  yn(195) = 'USTRGW'   ! u-gravity wave stress
  yn(196) = 'VSTRGW'   ! v-gravity wave stress
  yn(197) = 'VDISGW'   ! Gravity wave dissipation
  yn(198) = 'VGRAT'    ! Vegetation ratio
  yn(199) = 'VAROR'    ! Orographic variance (for surface runoff)
  yn(200) = 'VLT'      ! Leaf area index
  yn(201) = 'T2MAX'    ! Maximum 2m temperature
  yn(202) = 'T2MIN'    ! Minimum 2m temperature
  yn(203) = 'SRAD0U'   ! Top upward solar radiation
  yn(204) = 'SRADSU'   ! Surface upward solar radiation
  yn(205) = 'TRADSU'   ! Surface upward thermal radiation
  yn(206) = 'TSN'      ! Surface snow temperature
  yn(207) = 'TD3'      ! Deep soil temperature
  yn(208) = 'TD4'      !        "
  yn(209) = 'TD5'      !        "
  yn(210) = 'SEAICE'   ! Sea ice cover
  yn(211) = 'SICED'    ! Sea ice depth
  yn(212) = 'FOREST'   ! Vegetation type
  yn(213) = 'TEFF'     ! (Effective) skin temperature of sea ice
  yn(214) = 'TSMAX'    ! Maximum surface temperature
  yn(215) = 'TSMIN'    ! Minimum surface temperature
  yn(216) = 'WIMAX'    ! Maximum 10m wind speed
  yn(217) = 'TOPMAX'   ! Maximum convective cloud tops
  yn(218) = 'SNMEL'    ! Snow melt
  yn(219) = 'RUNTOC'   ! Surface runoff into ocean
  yn(220) = 'TSLIN'    ! Lin. surface temperature component + residuum           
  yn(221) = 'DSNAC'    ! Snow depth change
  yn(223) = 'ACLCAC'   ! Cloud cover (acc.)
  yn(224) = 'TKE'      ! Turbulent kinetic energy
  yn(225) = 'TKEM1'    ! Turbulent kinetic energy (t-1)
  yn(226) = 'FAO'      ! Fao data set
  yn(227) = 'RGCGN'    ! Heat capacity of soil
  yn(228) = 'SODIF'    ! Soil diffusivity
  yn(229) = 'WSMX'     ! Field capacity of soil
  yn(230) = 'QVI'      ! Vertical integrated specific humidity
  yn(231) = 'ALWCVI'   ! Vertical integrated liquid water content
  yn(232) = 'GLAC'     ! Glacier mask
  yn(233) = 'RUNLND'   ! Surface runoff not running into ocean

  nlevg3(162) = nlev   ! aclc
  nlevg3(223) = nlev   ! aclcac

  setval(201) = -99.   ! t2max
  setval(202) = 999.   ! t2min
  setval(214) = -99.   ! tsmax
  setval(215) = 999.   ! tsmin
  setval(216) = -99.   ! wimax
  setval(217) = 99999. ! topmax

  !-- 2. Read namelist postctl

  IF (p_parallel) THEN
    IF (p_parallel_io) THEN
      READ (nin,postctl)
    ENDIF
    CALL p_bcast (nunitdf, p_io)
    CALL p_bcast (nunitsp, p_io)
    CALL p_bcast (nunitun, p_io)
    CALL p_bcast (nung4x, p_io)
    CALL p_bcast (lppspe, p_io)
    CALL p_bcast (lppd, p_io)
    CALL p_bcast (lppvo, p_io)
    CALL p_bcast (lppt, p_io)
    CALL p_bcast (lppp, p_io)
    CALL p_bcast (lppq, p_io)
    CALL p_bcast (lppx, p_io)   
  ELSE
    READ (nin,postctl)
  ENDIF

  IF ( .NOT. lppspe) THEN
    lppd = .FALSE.
    lppvo = .FALSE.
    lppt = .FALSE.
    lppp = .FALSE.
    lppq = .FALSE.
    lppx = .FALSE.
  END IF

  !-- 3. Write namelist postctl

#ifdef DEBUG
  WRITE (nerr,postctl)
#endif

  !-- 4. Write selected fields (code,name,levels)

  iday = ncbase + (ntbase+dtime*(nstep+1))/dayl + 0.01
  CALL cd2dat(iday,id,im,iy)
  nhm = im
  nhy = iy
  nhd = id

  ! reset start time to initial time for output filename

  IF (lnwp .AND. lres) THEN
    iday = ncbase + (ntbase+dtime)/dayl + 0.01
    CALL cd2dat(iday,id,im,iy)
  END IF

  IF (p_parallel) THEN
    IF (p_parallel_io) THEN
      READ (nin,'(A)') ypath
    END IF
    CALL p_bcast (ypath, p_io)
  ELSE
    READ (nin,'(A)') ypath
  ENDIF

  iypt = 0
  DO i = 1, dir_limit
    IF (ypath(i:i)/=' ') THEN
      iypt = iypt + 1
      standard_grib_file(iypt:iypt) = ypath(i:i)
    END IF
  END DO
  ypph = standard_grib_file

  IF (p_parallel) THEN
    IF (p_parallel_io) THEN
      READ (nin,'(A9)') yexp
    END IF
    CALL p_bcast (yexp, p_io)
  ELSE
    READ (nin,'(A9)') yexp
  ENDIF

  ypph(iypt+1:iypt+10) = yexp(1:5) // 'HDIF'

  DO i = 1, 9
    iypt = iypt + 1
    IF (yexp(i:i)==' ') THEN
      standard_grib_file(iypt:iypt) = '_'
    ELSE
      standard_grib_file(iypt:iypt) = yexp(i:i)
    END IF
  END DO

  IF (lnwp) THEN
    WRITE (standard_grib_file(iypt+1:path_limit),'("_",I4.4,I2.2,I2.2)') &
         iy, im, id
  ELSE
    WRITE (standard_grib_file(iypt+1:path_limit),'("_",I4.4,I2.2,".",I2.2)') &
         nhy, nhm, nhd
  END IF
  extended_grib_file = standard_grib_file
  WRITE (extended_grib_file(iypt+1:path_limit),'("X",I4.4,I2.2,".",I2.2)') &
       nhy, nhm, nhd
  name_length = iypt+10

  IF (.NOT. lnwp .AND. p_pe == p_io) THEN
    CALL open_grib_file
  END IF

  ! Open optional file for g4x variables

  IF (lg4x .AND. nung4x/=nunitun) THEN
    CALL pbopen(nunitg4x,extended_grib_file,'w',iret)
    IF (iret/=0) THEN
      WRITE (nerr,*) 'Could not open file: ', extended_grib_file
      CALL finish('inipost','Run terminated.')
    END IF
  END IF

  IF (p_pe == p_io) THEN
    IF (.NOT. lnwp) THEN
      WRITE (nout, *) 'Experiment: ', yexp, ' - data file : ', &
           TRIM(standard_grib_file(INDEX(standard_grib_file,'/',.TRUE.)+1:))
    ELSE
      WRITE (nout, *) 'Experiment: ', yexp, ' - data file basename: ', &
           TRIM(standard_grib_file(INDEX(standard_grib_file,'/',.TRUE.)+1:))
    END IF
    WRITE (nout, *) 'Data path: ', &
         ADJUSTL(standard_grib_file(1:INDEX(standard_grib_file,'/',.TRUE.)))

    WRITE (nout,'(/)')
    IF (lg4x .AND. nung4x/=nunitun) THEN
      WRITE (nout, *) ' G4X file : ', extended_grib_file(1:name_length)
    END IF
  END IF

  ! Open optional file for statistics of horizontal diffusion

  IF (ldiahdf) THEN
    WRITE (ypph(iypt+1:path_limit),'("_",I4.4,I2.2,".",I2.2)') nhy, nhm, nhd
    IF (ypph==standard_grib_file) THEN
      WRITE (nerr, *) 'Data file and hdif-file have the same name!'
      WRITE (nerr, *) 'Choose other character section for ' // &
                      ' experiment name!'
      CALL finish('inipost','Run terminated.')
    END IF
    INQUIRE (file=ypph,exist=lexists)
    IF (lexists) THEN
      ystat = 'OLD'
    ELSE
      ystat = 'NEW'
    END IF
    OPEN (ndiahdf,file=ypph,status=ystat,form='formatted')
    IF (p_pe == p_io) THEN
      WRITE (nout, *) 'Statistics for horizontal diffusion on file ', &
           ypph(1:name_length)
    END IF
  END IF

  DO
    IF (p_parallel) THEN
      IF (p_parallel_io) THEN
        READ (nin,'(A1,I3,1X,A8,1X,99I3)',END=20) ycom, icode, ynam, &
             (ilev(i),i=1,nlevp1)
      END IF
      CALL p_bcast (ycom, p_io)
      CALL p_bcast (icode, p_io)
      CALL p_bcast (ynam, p_io)
      CALL p_bcast (ilev, p_io)
    ELSE
      READ (nin,'(A1,I3,1X,A8,1X,99I3)',END=20) ycom, icode, ynam, &
           (ilev(i),i=1,nlevp1)
    ENDIF

    IF (ycom/=' ') CYCLE
    IF (icode>255 .OR. ynam=='END') EXIT
    IF (icode<0) THEN
      CALL finish('inipost','Run terminated.')
    END IF
    IF (ilev(1)<0) THEN
      ilast = -ilev(1)
      IF (ilast>nlevp1) ilast = nlevp1
      DO i = 1, ilast
        ilev(i) = i
      END DO
    END IF

    IF (ynam(1:3)=='G3X') THEN
      IF (icode>=129) THEN
        WRITE (nout, *) 'Warning'
        WRITE (nout, *) 'Code ', icode, ' of extra field ' // ynam // &
                        ' conflicts with a preassigned code of a standard variable'
        icode = icode - 128
        WRITE (nout, *) 'Code of ' // ynam // ' changed to ', icode
      END IF
      READ (ynam(4:5),'(I2.2)') inx
      laccu(icode) = lxaccu(inx)
      IF (nxpbits(inx)==16 .OR. nxpbits(inx)==8 .OR. nxpbits(inx)==24 .OR. &
                        nxpbits(inx)==32) THEN
        npbits(icode) = nxpbits(inx)
      ELSE
        WRITE (nerr, *) ' Error: nxpbits of extra field ' // ynam // &
                        ' must have values  8, 16, 24 or 32 ! '
        CALL finish('inipost','Run terminated.')
      END IF
    END IF

    ! Read g4x codes

    IF (ynam(1:3)=='G4X') THEN
      IF (icode>=129) THEN
        WRITE (nout, *) 'Warning!'
        WRITE (nout, *) 'Code ', icode, ' of extra field ' // ynam // &
                        ' conflicts with a preassigned code of a standard variable'
        icode = icode - 128
        WRITE (nout, *) 'Code of ' // ynam // ' changed to ', icode
      END IF
      READ (ynam(4:5),'(I2.2)') inx
      laccu(icode) = l4xaccu(inx)
      IF (n4xpbits(inx)==16 .OR. n4xpbits(inx)==8 .OR. nxpbits(inx)==24 .OR. &
                       nxpbits(inx)==32) THEN
        npbits(icode) = n4xpbits(inx)
      ELSE
        WRITE (nerr, *) ' Error: n4xpbits of extra field ' // ynam // &
                        ' must have values  8, 16, 24 or 32 ! '
        CALL finish('inipost','Run terminated.')
      END IF
    END IF

    IF (yn(icode)==' ' .OR. yn(icode)==ynam) THEN
      yn(icode) = ynam
      IF (ilev(1)==0) ilev(1) = -100
      DO jlev = 1, nlevp1
        IF ((ilev(jlev)/=-100) .AND. (ilev(jlev)<0 .OR. ilev(jlev)>nlevp1)) THEN
          WRITE (nerr, *) ' Illegal level ', ilev(jlev), ' for ' // ynam
          CALL finish('inipost','Run terminated.')
        END IF
        npplev(jlev,icode) = ilev(jlev)
      END DO
    ELSE
      WRITE (nerr, *) ' Code - name conflict in inipost'
      WRITE (nerr, *) ' Code           : ', icode
      WRITE (nerr, *) ' Expected  name : ', yn(icode)
      WRITE (nerr, *) ' Requested name : ', ynam
      CALL finish('inipost','Run terminated.')
    END IF

    IF (icode == 999) EXIT

  END DO

20 CONTINUE

  ! Force postprocessing of surface geopotential
  npplev(1,129) = -100

  ncdmin = 0
  ncdmax = 129
  DO jfield = 1, 234
    IF (npplev(1,jfield)/=0) THEN
      IF (yn(jfield)(1:3)/='G4X') THEN
        IF (ncdmin==0) ncdmin = jfield
        ncdmax = jfield
      END IF
    END IF
  END DO

  IF (p_pe == p_io) THEN
    WRITE (nout,'(a)') &
         ' Postprocessing grid point fields'
    WRITE (nout,'(a)') &
         ' -------------------------------------------------------'

    DO jfield = 1, 234
      IF (npplev(1,jfield)/=0) THEN
        WRITE (nout,'(1x,i3,1x,a,i4,99i3)') jfield, yn(jfield), &
             (npplev(i,jfield),i=1,nlevp1)
      END IF
    END DO
  END IF

END SUBROUTINE inipost
