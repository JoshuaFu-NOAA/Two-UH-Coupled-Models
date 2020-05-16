MODULE mo_couple
  ! interface between ocean and atmosphere
  !
  ! I.Kirchner, IfM, FU-Berlin, Sep/2005


  ! transfer from atmospher to ocean
  !
  !   wind speed 10 meter
  !   surface wind stress (u,v)
  !   heat fluxes (latent, sensible, solar, thermal)
  !
  ! transfer from ocean to atmosphere
  !
  !   sea surface temperatur in tropics
  !
   
  USE mo_exception,     ONLY: message, finish
  USE mo_mpi,           ONLY: p_pe, p_io              !fu++
  USE mo_doctor,        ONLY: nout
  USE mo_constants,     ONLY: dayl
  USE mo_kind,          ONLY: dp
  USE mo_nudging,       ONLY: sstn, nudg_stop          !fu++

  USE mo_io,            ONLY: vlon, vlat
  USE mo_memory_g3a,    ONLY: slmm, tsm
  USE mo_landsea,       ONLY: bzone, sstatm, amask

  USE mo_year,          ONLY: cd2dat, ic2ymd
  USE mo_transpose,     ONLY: gather_gp, scatter_gp
  USE mo_decomposition, ONLY: dc=>local_decomposition, dcg=>global_decomposition
  USE mo_control,       ONLY: dtime, lctime, lwtime, nctime,  &
                           ngl, nlon, lcouple, imtt, ncbase, ntbase, lonudg

  USE mo_start_dataset, ONLY: nstep, lres
  USE mo_filename,      ONLY: path_limit, standard_grib_file
  USE mo_time_control,  ONLY: time_inc_sec, TIME_INC_DAYS

  IMPLICIT NONE

  PRIVATE

  REAL(dp), POINTER, PUBLIC, SAVE :: ssto(:,:,:)
  REAL(dp), POINTER, PUBLIC, SAVE :: bzo(:,:,:)

  INTEGER, PARAMETER :: OCE_RESTART_WRITE = 1
  INTEGER, PARAMETER :: OCE_RESTART_READ  = 2

  PUBLIC :: oce_init          ! called in control
  PUBLIC :: oce_run           ! called in stepon

  REAL, POINTER :: slmfu(:,:)
  REAL, POINTER :: tsmfu(:,:)
  REAL, POINTER :: gl_sst(:,:,:)
  REAL, POINTER :: gl_bzone(:,:,:)

  REAL    :: runday, dta, dto, day1, tseg
  INTEGER :: numcpl, numatm, numocn
  COMMON /fudate/runday,dta,dto,day1,tseg,numcpl,numatm,numocn

  CHARACTER(path_limit+4) :: ocean_filename

  CHARACTER(256) :: mess_text = ''

  REAL,TARGET,ALLOCATABLE :: sstfu(:,:,:)
  REAL,TARGET,ALLOCATABLE :: bzfu(:,:,:)

  REAL, ALLOCATABLE, SAVE, PUBLIC, TARGET :: flatfx(:,:)
  REAL, ALLOCATABLE, SAVE, PUBLIC, TARGET :: fsenfx(:,:)
  REAL, ALLOCATABLE, SAVE, PUBLIC, TARGET :: fsolfx(:,:)
  REAL, ALLOCATABLE, SAVE, PUBLIC, TARGET :: fradfx(:,:)
  REAL, ALLOCATABLE, SAVE, PUBLIC, TARGET :: ftaux(:,:)
  REAL, ALLOCATABLE, SAVE, PUBLIC, TARGET :: ftauy(:,:)
  REAL, ALLOCATABLE, SAVE, PUBLIC, TARGET :: fwind(:,:)
  !fu++
  REAL, ALLOCATABLE, SAVE, PUBLIC, TARGET :: fu_sst(:,:)
  REAL, ALLOCATABLE, SAVE, PUBLIC, TARGET :: sstnfu(:,:)
  ! local fields for oce_collect

  ! pick up in VDIFF (m_vdiff.f90)
  REAL, ALLOCATABLE, PUBLIC :: wind10_local(:,:)
  REAL, ALLOCATABLE, PUBLIC ::   vstr_local(:,:)
  REAL, ALLOCATABLE, PUBLIC ::   ustr_local(:,:)
  REAL, ALLOCATABLE, PUBLIC ::   ahfs_local(:,:)
  REAL, ALLOCATABLE, PUBLIC ::   ahfl_local(:,:)

  ! pick up in RADHEAT (m_radheat.f90)
  REAL, ALLOCATABLE, PUBLIC ::  trads_local(:,:)
  REAL, ALLOCATABLE, PUBLIC ::  srads_local(:,:)

  REAL, PUBLIC :: accu_time = 0.0      ! accumulation intervall in seconds
  INTEGER      :: ioerr

!fu++  LOGICAL, PUBLIC :: ioflag = .FALSE.  ! do i/o on that PE

  ! fields for restart
  REAL, ALLOCATABLE, PUBLIC :: restart_u1(:,:,:)
  REAL, ALLOCATABLE, PUBLIC :: restart_u2(:,:,:)
  REAL, ALLOCATABLE, PUBLIC :: restart_v1(:,:,:)
  REAL, ALLOCATABLE, PUBLIC :: restart_v2(:,:,:)
  REAL, ALLOCATABLE, PUBLIC :: restart_h1(:,:,:)
  REAL, ALLOCATABLE, PUBLIC :: restart_h2(:,:,:)
  REAL, ALLOCATABLE, PUBLIC :: restart_tmp(:,:,:)
  INTEGER, PUBLIC           :: restart_itt
  !fu++
  REAL, ALLOCATABLE, PUBLIC :: restart_sst(:,:)

  INTEGER, PARAMETER :: oce_xdim = 720
  INTEGER, PARAMETER :: oce_ydim = 121

  INTEGER, PARAMETER :: restart_unit = 51
  INTEGER, PARAMETER, PUBLIC :: output_unit  = 117

CONTAINS

  !================================================================================
  ! initialisation of ocean

  SUBROUTINE oce_init

    INTEGER :: iday, id1, im1, iy

!fu++    ioflag = (p_parallel .AND. p_parallel_io) .OR. (.NOT.p_parallel)

    IF (p_pe == p_io) THEN
       ALLOCATE (slmfu(nlon,ngl))
       ALLOCATE (tsmfu(nlon,ngl))

       ALLOCATE (sstfu(nlon,ngl,0:0))
       ALLOCATE (bzfu(nlon,ngl,0:0))
    END IF

    CALL gather_gp (slmfu,slmm,dcg)
    CALL gather_gp (tsmfu,tsm,dcg)

    ALLOCATE (ssto(dc%nglon,dc%nglat,0:0))
    ALLOCATE (bzo(dc%nglon,dc%nglat,0:0))

    ! allocate ocean transfer data
    ALLOCATE (wind10_local (dc%nglon, dc%nglat))
    ALLOCATE (  vstr_local (dc%nglon, dc%nglat))
    ALLOCATE (  ustr_local (dc%nglon, dc%nglat))
    ALLOCATE (  ahfs_local (dc%nglon, dc%nglat))
    ALLOCATE (  ahfl_local (dc%nglon, dc%nglat))
    ALLOCATE ( trads_local (dc%nglon, dc%nglat))
    ALLOCATE ( srads_local (dc%nglon, dc%nglat))
    ALLOCATE (   sstnfu    (dc%nglon, dc%nglat))

    wind10_local(:,:) = 0.0
    vstr_local  (:,:) = 0.0
    ustr_local  (:,:) = 0.0
    ahfs_local  (:,:) = 0.0
    ahfl_local  (:,:) = 0.0
    trads_local (:,:) = 0.0
    srads_local (:,:) = 0.0

    ! initialise ocean

    IF (p_pe == p_io) THEN

       WRITE(nout,*) 'lon=',vlon(1:nlon)
       WRITE(nout,*) 'lat=',vlat(1:ngl)
       WRITE(nout,*) 'lat,slm=',ngl/2,slmfu(1:nlon,ngl/2)


       !==> establish common block "/fudate/"=================>
       runday = 360.                             !one-year
       dta    = dtime                            !second ATM step length
       dto    = 14400.                           !second OCE step length
       day1   = 24.*3600.                        !1 day
!ik    tseg   = 1.*day1                          !coupling interval
       tseg   =  time_inc_sec(nctime, TIME_INC_DAYS)

       !==> total segments for ocean & atmosphere
       numcpl = INT(runday*day1/tseg)
       numatm = INT(tseg/dta)
       numocn = INT(tseg/dto)

       WRITE(mess_text,*) 'c****************************************c'
       CALL message('oce_init',mess_text)
       WRITE(mess_text,*) 'coupling segment interval',tseg
       CALL message('oce_init',mess_text)
       WRITE(mess_text,*) 'ATM steps/segment=', numatm
       CALL message('oce_init',mess_text)
       WRITE(mess_text,*) 'OCN steps/segment=', numocn
       CALL message('oce_init',mess_text)
       WRITE(mess_text,*) 'c****************************************c'
       CALL message('oce_init',mess_text)

       !========================================================>
       !for deposit ocean outputs

       iday = ncbase + (ntbase+dtime*(nstep+1))/dayl + 0.01 !input data
       CALL cd2dat(iday,id1,im1,iy)

       WRITE(mess_text,*) 'initial run? imtt=(1,1000),yr=',imtt,iy
       CALL message('oce_init',mess_text)

       ocean_filename = TRIM(standard_grib_file)//'_oce'

       OPEN(output_unit,file=ocean_filename,form='unformatted',&
            access='sequential')

       WRITE(mess_text,*) 'open ocean outputfile ',TRIM(ocean_filename)
       CALL message('oce_init',mess_text)

       OPEN(restart_unit,form='unformatted', access='sequential')
       CALL message('oce_init','open ocean restart file')

       !========================================================> 
       !>>get the coupled land-sea masks for ocean and atmosphere

       CALL amask(slmfu,vlat,vlon)
       !>>set the initial condition for intermediate ocean model

       ! allocate restart fields
       ALLOCATE(restart_u1 (oce_xdim,oce_ydim,2))
       ALLOCATE(restart_u2 (oce_xdim,oce_ydim,2))
       ALLOCATE(restart_v1 (oce_xdim,oce_ydim,2))
       ALLOCATE(restart_v2 (oce_xdim,oce_ydim,2))
       ALLOCATE(restart_h1 (oce_xdim,oce_ydim,2))
       ALLOCATE(restart_h2 (oce_xdim,oce_ydim,2))
       ALLOCATE(restart_tmp(oce_xdim,oce_ydim,2))
       ALLOCATE(restart_sst(nlon,ngl))   !fu++

    END IF

    ! read restart data of ocean

    IF ( lres .OR. (imtt > 0) ) THEN
       CALL oce_restart(OCE_RESTART_READ)

       IF (p_pe == p_io) THEN
          WRITE(mess_text,*) ' OCE RESTART ITT ',restart_itt
          CALL message('oce_init',mess_text)
          WRITE(mess_text,*) ' OCE RESTART U1  ',MINVAL(restart_u1), MAXVAL(restart_u1)
          CALL message('oce_init',mess_text)
          WRITE(mess_text,*) ' OCE RESTART U2  ',MINVAL(restart_u2), MAXVAL(restart_u2)
          CALL message('oce_init',mess_text)
          WRITE(mess_text,*) ' OCE RESTART V1  ',MINVAL(restart_v1), MAXVAL(restart_v1)
          CALL message('oce_init',mess_text)
          WRITE(mess_text,*) ' OCE RESTART V2  ',MINVAL(restart_v2), MAXVAL(restart_v2)
          CALL message('oce_init',mess_text)
          WRITE(mess_text,*) ' OCE RESTART H1  ',MINVAL(restart_h1), MAXVAL(restart_h1)
          CALL message('oce_init',mess_text)
          WRITE(mess_text,*) ' OCE RESTART H2  ',MINVAL(restart_h2), MAXVAL(restart_h2)
          CALL message('oce_init',mess_text)
          WRITE(mess_text,*) ' OCE RESTART TMP ',MINVAL(restart_tmp),MAXVAL(restart_tmp)
          CALL message('oce_init',mess_text)
          WRITE(mess_text,*) ' OCE RESTART SST ',MINVAL(restart_sst),MAXVAL(restart_sst) !fu++
          CALL message('oce_init',mess_text)  !fu++

       END IF

    ELSE
       CALL message('oce_init','initialize ocean')

    END IF

    IF (p_pe == p_io) THEN

       CALL setupo(imtt)
       !>>get the initial SST for the AGCM

       ALLOCATE(sstatm(1:nlon,1:ngl))

       !fu++ CALL sstint(nlon,ngl,sstatm)
       !using restart SST to start the rerun

       sstatm = restart_sst

       !define pointers
       sstfu(:,:,0) = sstatm(:,:)
       bzfu(:,:,0)  = bzone(:,:)

       gl_sst   => sstfu(:,:,0:0)
       gl_bzone => bzfu(:,:,0:0) 

    END IF


    !scatter sstatm and bzone
    CALL scatter_gp (gl_sst,  ssto(:,:,0:0), dcg)
    CALL scatter_gp (gl_bzone, bzo(:,:,0:0), dcg)

    ! variables for data exchange, global fields

    IF (p_pe == p_io) THEN
       IF(.NOT.ALLOCATED(flatfx)) ALLOCATE( flatfx (nlon,ngl))
       IF(.NOT.ALLOCATED(fsenfx)) ALLOCATE( fsenfx (nlon,ngl))
       IF(.NOT.ALLOCATED(fsolfx)) ALLOCATE( fsolfx (nlon,ngl))
       IF(.NOT.ALLOCATED(fradfx)) ALLOCATE( fradfx (nlon,ngl))
       IF(.NOT.ALLOCATED(ftaux))  ALLOCATE( ftaux  (nlon,ngl))
       IF(.NOT.ALLOCATED(ftauy))  ALLOCATE( ftauy  (nlon,ngl))
       IF(.NOT.ALLOCATED(fwind))  ALLOCATE( fwind  (nlon,ngl))
     !fu++
       IF(.NOT.ALLOCATED(fu_sst))  ALLOCATE( fu_sst  (nlon,ngl))
    END IF

  END SUBROUTINE oce_init


  !================================================================================
  ! write and read ocean restart file

  SUBROUTINE oce_restart(flag)

    INTEGER, INTENT(in) :: flag

    IF (p_pe == p_io) THEN

       SELECT CASE(flag)

       CASE(OCE_RESTART_READ)
          CALL message('oce_restart','read ocean restart data')

          READ(restart_unit) &
               restart_itt, restart_u2, restart_v2, restart_h2, &
               restart_tmp, restart_u1, restart_v1, restart_h1, restart_sst
          READ(restart_unit,iostat=ioerr) &
               accu_time, wind10_local, vstr_local, ustr_local, &
               ahfs_local, ahfl_local, trads_local, srads_local
          IF (ioerr /= 0) &
               CALL message('oce_restart','missing second record in restartfile')

          WRITE(mess_text,*) 'read ocean step ',restart_itt,' from historyfile'
          CALL message('oce_restart',mess_text)

       CASE(OCE_RESTART_WRITE)
          CALL message('oce_restart','write ocean restart data')

          WRITE(mess_text,*) ' OCE RESTART ITT ',restart_itt 
          CALL message('oce_restart',mess_text)
          WRITE(mess_text,*) ' OCE RESTART U1  ',MINVAL(restart_u1), MAXVAL(restart_u1)
          CALL message('oce_restart',mess_text)
          WRITE(mess_text,*) ' OCE RESTART U2  ',MINVAL(restart_u2), MAXVAL(restart_u2)
          CALL message('oce_restart',mess_text)
          WRITE(mess_text,*) ' OCE RESTART V1  ',MINVAL(restart_v1), MAXVAL(restart_v1)
          CALL message('oce_restart',mess_text)
          WRITE(mess_text,*) ' OCE RESTART V2  ',MINVAL(restart_v2), MAXVAL(restart_v2)
          CALL message('oce_restart',mess_text)
          WRITE(mess_text,*) ' OCE RESTART H1  ',MINVAL(restart_h1), MAXVAL(restart_h1)
          CALL message('oce_restart',mess_text)
          WRITE(mess_text,*) ' OCE RESTART H2  ',MINVAL(restart_h2), MAXVAL(restart_h2)
          CALL message('oce_restart',mess_text)
          WRITE(mess_text,*) ' OCE RESTART TMP ',MINVAL(restart_tmp),MAXVAL(restart_tmp)
          CALL message('oce_restart',mess_text)
          !fu++ 
	  WRITE(mess_text,*) ' OCE RESTART SST ',MINVAL(restart_sst),MAXVAL(restart_sst)
          CALL message('oce_restart',mess_text)

          REWIND(restart_unit)
          WRITE(restart_unit) &
               restart_itt, restart_u2, restart_v2, restart_h2, &
               restart_tmp, restart_u1, restart_v1, restart_h1, restart_sst  !fu++
          WRITE(restart_unit) &
               accu_time, wind10_local, vstr_local, ustr_local, &
               ahfs_local, ahfl_local, trads_local, srads_local

          WRITE(mess_text,*) 'save ocean step ',restart_itt,' to historyfile'
          CALL message('oce_restart',mess_text)

       END SELECT

    END IF

  END SUBROUTINE oce_restart


  !================================================================================
  ! collect data for the ocean

  SUBROUTINE oce_collect

    USE mo_control,       ONLY: nctime
    USE mo_transpose,     ONLY: gather_gp

    REAL          :: quot
    REAL, POINTER :: collect_field(:,:)

    ! gather input arrays

    ! AHFS surface sensible heat flux [W/m*m]
    collect_field => fsenfx(:,:)
    CALL gather_gp(collect_field,ahfs_local,dcg)

    ! AHFL surface latent heat flux [W/m*m]
    collect_field => flatfx(:,:)
    CALL gather_gp(collect_field,ahfl_local,dcg)

    ! WIND10 10 meter wind speed [m/s]
    collect_field => fwind(:,:)
    CALL gather_gp(collect_field,wind10_local,dcg)

    ! SRADS solar surface flux [W/m*m]
    collect_field => fsolfx(:,:)
    CALL gather_gp(collect_field,srads_local,dcg)

    ! TRADS thermal surface flux [W/m*m]
    collect_field => fradfx(:,:)
    CALL gather_gp(collect_field,trads_local,dcg)

    ! USTR surface u-stress [Pa] -> 0.1 Pa
    collect_field => ftaux(:,:)
    CALL gather_gp(collect_field,ustr_local,dcg)

    ! VSTR surface v-stress [Pa] -> 0.1 Pa
    collect_field => ftauy(:,:)
    CALL gather_gp(collect_field,vstr_local,dcg)

    !fu++ NCEP Nudging SST (K) 
     sstn => sstnfu(:,:)
    collect_field => fu_sst(:,:)
    CALL gather_gp(collect_field,sstnfu,dcg)

    IF (p_pe == p_io) THEN

       IF ( accu_time > 0.0 ) THEN
          quot = 1.0/accu_time
       ELSE
          CALL finish('oce_collect','wrong setting of accumulation time')
       END IF

       WRITE(mess_text,*) ' accumulation factor [1/sec] ',quot
       CALL message('oce_collect',mess_text)

       fsenfx(:,:) = fsenfx(:,:)*quot
       flatfx(:,:) = flatfx(:,:)*quot
       fwind(:,:)  = fwind(:,:) *quot
       fsolfx(:,:) = fsolfx(:,:)*quot
       fradfx(:,:) = fradfx(:,:)*quot
       ftaux(:,:)  = ftaux(:,:) *quot*10.
       ftauy(:,:)  = ftauy(:,:) *quot*10.
       fu_sst(:,:) = fu_sst(:,:)-273.16

    END IF

    ! zero exchange fields

    wind10_local(:,:) = 0.0
      vstr_local(:,:) = 0.0
      ustr_local(:,:) = 0.0
      ahfs_local(:,:) = 0.0
      ahfl_local(:,:) = 0.0
     trads_local(:,:) = 0.0
     srads_local(:,:) = 0.0

     accu_time = 0.0

    RETURN

  END SUBROUTINE oce_collect



  !================================================================================
  ! run the ocean model

  SUBROUTINE oce_run

    REAL :: fsflux(nlon,ngl)
    INTEGER :: iday, sheadd

    IF (lctime) THEN
       

       !collect input data
       CALL oce_collect

       IF (p_pe == p_io) THEN

          !current ocean nudging 
          iday = ncbase + (ntbase+dtime*(nstep+1))/dayl + 0.01 !input data
          sheadd =ic2ymd(iday)
          if(lonudg) then
          if(nudg_stop < sheadd) lonudg = .FALSE.
          end if
          write(nout,*) "nudging NCEP SST?", lonudg

          fsflux = flatfx + fsenfx + fsolfx + fradfx

          WRITE(mess_text,*) 'atmospheric SST at (140W,EQ):', tsmfu(79,32)-273.16
          CALL message('oce_run',mess_text)
          WRITE(mess_text,*) 'surface heat flux at (140W,EQ):', fsflux(79,32)
          CALL message('oce_run',mess_text)


          WRITE(nout,*) ' solar Min,Max ',MINVAL(fsolfx),MAXVAL(fsolfx)
          WRITE(nout,*) ' taux Min,Max ',MINVAL(ftaux),MAXVAL(ftaux)
          WRITE(nout,*) ' tauy Min,Max ',MINVAL(ftauy),MAXVAL(ftauy)
          WRITE(nout,*) ' fu_sst Min,Max ',MINVAL(fu_sst),MAXVAL(fu_sst)

          !run IOM and return SST for AGCM
          !inputs:lwtime,ftaux,ftauy,fsflux,fsolfx,fwind,nlon,ngl
          !output:sstatm

          CALL ocean(sstatm,ftaux,ftauy,fsflux,fsolfx,fwind,nlon,ngl,fu_sst)

          !define pointers
          sstfu(:,:,0) = sstatm(:,:)
          bzfu(:,:,0)  =  bzone(:,:)

          gl_sst   => sstfu(:,:,0:0)
          gl_bzone =>  bzfu(:,:,0:0) 

       END IF

       !scatter sstatm and bzone
       CALL scatter_gp (gl_sst,  ssto(:,:,0:0), dcg)
       CALL scatter_gp (gl_bzone, bzo(:,:,0:0), dcg)

    END IF

    ! write restart data of ocean

    IF (lwtime) CALL oce_restart(OCE_RESTART_WRITE)

  END SUBROUTINE oce_run

END MODULE mo_couple
