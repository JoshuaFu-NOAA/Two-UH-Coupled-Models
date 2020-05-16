MODULE mo_column
!
!+ $Id: mo_column.f90,v 1.25 2000/08/05 16:55:00 m214030 Exp $
!
! This module gathers all routines and variables requred to run
! ECHAM in column model mode.
!
! Working principle: 
! All terms which cannot be calculated by the column model are stored 
! on a file by the 3D model and read in later by the column model.
!
! Authors:
!
! A. Rhodin, MPI, November 1999, original source
!
  !
  ! modules used
  !
#ifdef NAG
  USE f90_unix,         ONLY: flush
#endif
  USE mo_decomposition, ONLY: gdc => global_decomposition, &
                              ldc => local_decomposition
  USE mo_transpose,     ONLY: tag_gather_gp
  USE mo_mpi,           ONLY: p_pe, p_io,                  & 
                              p_bcast, p_send, p_recv,     &
                              p_parallel, p_parallel_io
  USE mo_doctor,        ONLY: nout, nerr, nin
  USE mo_exception,     ONLY: finish
  USE mo_start_dataset, ONLY: nstep, lres
  USE mo_control,       ONLY: nresum, vct, nvclev, nlevp1, dtime
  !
  ! Some more modules are accessed by module subroutines to store/restore
  ! module variables:
  !
  ! mo_scan_buffer, mo_memory_g1a, mo_memory_g2a, mo_sc1, mo_gaussgrid
  !
  IMPLICIT NONE
  !
  ! public entities
  !
  PRIVATE
  PUBLIC :: inicolumn    ! subroutine: read column model namelist
  PUBLIC :: setcolumn    ! subroutine: set module variables within this module
  PUBLIC :: resetcolumn  ! subroutine: deinitialize column model
  PUBLIC :: get_col_ffti ! subroutine: get information for column model
  PUBLIC :: get_col_tran ! subroutine: get information for column model
  PUBLIC :: cal_col_expl ! subroutine: explicit integration for column model
  PUBLIC :: write_column ! subroutine: write field to plot file

  PUBLIC :: lat_1d       ! global latitude  index of column
  PUBLIC :: lon_1d       ! global longitude index of column
  PUBLIC :: lat_1db      ! global latitude  index (edges of area)
  PUBLIC :: lon_1db      ! global longitude index (edges of area)
  !
  ! module variables (namelist columnctl)
  ! 
  INTEGER,PARAMETER:: mc        = 100! max.number of columns in namelist
  INTEGER,PARAMETER:: ml        = 100! max.number of levels (nudging coeffs.)

  CHARACTER(len=8) :: comode    = '' ! mode
  INTEGER          :: lat_1d(mc)= -1 ! global latitude index of column
  INTEGER          :: lon_1d(mc)= -1 ! global longitude index of column
  INTEGER          :: lat_1db(2)= -1 ! global latitude index of column (area)
  INTEGER          :: lon_1db(2)= -1 ! global longitude index of column (area)
  INTEGER          :: nforce    =  1 ! force until timestep nforce
                                     ! flags for column model integration
  LOGICAL          :: lrewind=.FALSE.! rewind before reading each timestep
  INTEGER          :: int_t     = 1  ! temperature
  INTEGER          :: int_ps    = 1  ! log surface pressure
  INTEGER          :: int_uv    = 1  ! wind fields
                                     ! flags for column model forcing
  INTEGER          :: for_t     = 1  ! temperature
  INTEGER          :: for_ps    = 1  ! log surface pressure
  INTEGER          :: for_uv    = 1  ! wind fields
  INTEGER          :: for_q     = 0  ! specific humidity
  INTEGER          :: for_x     = 0  ! cloud water content
  INTEGER          :: for_xt    = 0  ! tracer
  INTEGER          :: for_def   = 1  ! default value
  INTEGER          :: for_tr    = 1  ! qe,xe,xte tendencies (due to transport) 
  INTEGER          :: for_trp   = 1  ! surf. pressure modifications(transport)
  INTEGER          :: for_te    = 0  ! additional temperature tendency term
  INTEGER          :: for_d     = 1  ! divergence
  INTEGER          :: for_vo    = 1  ! vorticity
  INTEGER          :: for_dt    = 1  ! temperature gradient
  INTEGER          :: for_dp    = 1  ! log surface pressure gradient
  INTEGER          :: for_zm    = 0  ! zonal means calculated in spectral space
  REAL             :: ug_1d     = 0. ! geostrophic wind component
  REAL             :: vg_1d     = 0. ! geostrophic wind component
  INTEGER          :: k_vg      = 0  ! estimate geostr. wind above this level

  REAL             :: nudguv(ml)= 0. ! nudging weigth for u,v
  REAL             :: nudgt (ml)= 0. ! nudging weigth for temperature
  REAL             :: nudgp     = 0. ! nudging weigth for log sfc pressure
  REAL             :: nudgq (ml)= 0. ! nudging weigth for water vapour
  REAL             :: nudgx (ml)= 0. ! nudging weigth for liq. water + ice
  REAL             :: nudgxt(ml)= 0. ! nudging weigth for tracer
  !
  ! other module variables
  ! 
  LOGICAL          :: lnudge=.FALSE.! nudging applied
  INTEGER          :: p_col    = -1 ! PE running the column model
  INTEGER          :: p_ful    = -1 ! full model PE holding the column
  INTEGER          :: nlat          ! number of latitudes
  INTEGER          :: nlon          ! number of longitudes
  INTEGER          :: ncol          ! number of columns
  INTEGER          :: nlcol         ! number of columns on this pe
  INTEGER          :: iread         ! index of column in input files
  INTEGER,ALLOCATABLE :: lats(:)
  INTEGER,ALLOCATABLE :: lons(:)
  INTEGER,ALLOCATABLE :: pes (:)
  INTEGER,ALLOCATABLE :: llats(:)
  INTEGER,ALLOCATABLE :: llons(:)
  INTEGER          :: unit_f1d = 51 ! fortran unit for forcing file
  INTEGER          :: unit_r1d = 52 ! fortran unit for residuum file
  INTEGER          :: unit_b1d = 53 ! fortran unit for results before correct. 
  INTEGER          :: unit_a1d = 54 ! fortran unit for results after correct.
  LOGICAL          :: fapp  =.FALSE.! forcing applied in this step
  REAL,ALLOCATABLE :: ug(:,:,:)     ! geostrophic wind component profile
  REAL,ALLOCATABLE :: vg(:,:,:)     ! geostrophic wind component profile
  !
  ! interfaces
  !
  INTERFACE write_column            ! write field x to plot file:
    MODULE PROCEDURE write_column3  ! ( x(lon,lev,lat) ,name)
    MODULE PROCEDURE write_column2  ! ( x(lon,    lat) ,name)
    MODULE PROCEDURE write_column3r ! ( x(lon,lev)     ,name ,jlat)
    MODULE PROCEDURE write_column2r ! ( x(lon)         ,name ,jlat)
  END INTERFACE

CONTAINS
!==============================================================================
! Initialization routines
!==============================================================================
  SUBROUTINE inicolumn (lcolumn, lvctch, nlev, nlevp1, nvclev, vct)
  LOGICAL ,INTENT (out)   :: lcolumn ! .true. if a column model is run
  LOGICAL ,INTENT (out)   :: lvctch  ! .true. if vct table is changed
  INTEGER ,INTENT (inout) :: nlev    ! levels (changed if lvctch == .true.)
  INTEGER ,INTENT (inout) :: nlevp1  ! levels (changed if lvctch == .true.)
  INTEGER ,INTENT (inout) :: nvclev  ! levels (changed if lvctch == .true.)
  REAL    ,POINTER        :: vct (:) ! hybrid level coefficients     ( '' )
    REAL ,POINTER :: x (:)
    !
    ! Initialize module variables from namelist.
    !
    INCLUDE 'columnctl.inc'
    !
    ! default values
    !
    lat_1d  = -1
    lon_1d  = -1
    lat_1db = -1
    lon_1db = -1
    nforce  =  1
    comode  = ''
    lrewind = .FALSE.
    lvctch  = .FALSE.
    int_t   =  1
    int_ps  =  1
    int_uv  =  1
    for_t   =  1
    for_ps  =  1
    for_uv  =  1
    for_q   =  0
    for_x   =  0
    for_xt  =  0
    for_def =  1
    for_d   = -1
    for_vo  = -1
    for_dt  = -1
    for_dp  = -1
    for_tr  = -1
    for_trp = -1
    for_zm  = -1
    for_te  =  0
    ug_1d   =  0.
    vg_1d   =  0.
    k_vg    =  0
    nudguv  =  0.
    nudgt   =  0.
    nudgp   =  0.
    nudgq   =  0.
    nudgx   =  0.
    nudgxt  =  0.
    !
    ! read namelist
    !
    IF (p_parallel) THEN
      IF (p_parallel_io) THEN
        READ (nin,columnctl)
      ENDIF
      CALL p_bcast (lat_1d,  p_io)
      CALL p_bcast (lon_1d,  p_io)
      CALL p_bcast (lat_1db, p_io)
      CALL p_bcast (lon_1db, p_io)
      CALL p_bcast (comode,  p_io)
      CALL p_bcast (nforce,  p_io)
      CALL p_bcast (lrewind, p_io)
      CALL p_bcast (lvctch,  p_io)
      CALL p_bcast (int_t,   p_io)
      CALL p_bcast (int_ps,  p_io)
      CALL p_bcast (int_uv,  p_io)
      CALL p_bcast (for_t,   p_io)
      CALL p_bcast (for_ps,  p_io)
      CALL p_bcast (for_uv,  p_io)
      CALL p_bcast (for_q,   p_io)
      CALL p_bcast (for_x,   p_io)
      CALL p_bcast (for_xt,  p_io)
      CALL p_bcast (for_def, p_io)
      CALL p_bcast (for_d,   p_io)
      CALL p_bcast (for_vo,  p_io)
      CALL p_bcast (for_dt,  p_io)
      CALL p_bcast (for_dp,  p_io)
      CALL p_bcast (for_tr,  p_io)
      CALL p_bcast (for_trp, p_io)
      CALL p_bcast (for_zm,  p_io)
      CALL p_bcast (for_te,  p_io)
      CALL p_bcast (ug_1d,   p_io)
      CALL p_bcast (vg_1d,   p_io)
      CALL p_bcast (k_vg    ,p_io)
      CALL p_bcast (nudguv  ,p_io)
      CALL p_bcast (nudgt   ,p_io)
      CALL p_bcast (nudgp   ,p_io)
      CALL p_bcast (nudgq   ,p_io)
      CALL p_bcast (nudgx   ,p_io)
      CALL p_bcast (nudgxt  ,p_io)
    ELSE
      READ (nin,columnctl)
    ENDIF
    !
    ! set flags for individual variables
    !
    IF (for_d  ==-1) for_d   = for_def
    IF (for_vo ==-1) for_vo  = for_def
    IF (for_dt ==-1) for_dt  = for_def
    IF (for_dp ==-1) for_dp  = for_def
    IF (for_tr ==-1) for_tr  = for_def
    IF (for_trp==-1) for_trp = for_def
    IF (for_zm ==-1) for_zm  = 0
    lnudge = ANY(nudguv/=0.) .OR. ANY(nudgt/=0.) .OR.     nudgp /=0. &
        .OR. ANY(nudgq /=0.) .OR. ANY(nudgx/=0.) .OR. ANY(nudgxt/=0.)
    !
    ! calculate ncol,nlat,nlon
    !
    IF (ANY(lat_1d/=-1).OR.ANY(lon_1d/=-1)) THEN
      lat_1db    = -1
      lon_1db    = -1
      nlat       = COUNT(lat_1d/=-1)
      nlon       = COUNT(lon_1d/=-1)
      IF(nlat/=nlon)CALL finish('inicolumn','error in lat_1db/lon_1db')
      ncol=nlat
    ELSE
      IF (lat_1db(2)==-1) lat_1db(2) = lat_1db(1)
      nlat = MAX(0,lat_1db(2)-lat_1db(1)+1)
      IF (ANY(lat_1db==-1)) nlat = 0
      IF (lon_1db(2)==-1) lon_1db(2) = lon_1db(1)
      nlon = MAX(0,lon_1db(2)-lon_1db(1)+1)
      IF (ANY(lon_1db==-1)) nlon = 0
      ncol  = nlat * nlon
    ENDIF
    !
    ! check value of comode
    !
    SELECT CASE (comode)
    CASE ('')
    CASE ('traject')
      IF(ncol==0) CALL finish('inicolumn','location of columns not specified') 
    CASE ('resid','force','add','free')
      IF(ncol/=1) CALL finish('inicolumn',&
                              'exactly one column location must be specified')
    CASE default
      CALL finish('inicolumn','comode='//comode)
    END SELECT
    !
    ! set return arguments, change vertical levels
    !
    lcolumn = (comode /= '' .AND. comode /= 'traject')
    lvctch  = lvctch .AND. lcolumn
    IF(lvctch) THEN
      NULLIFY (x)
      IF (p_pe==p_io) THEN
        OPEN (unit_f1d,file='forcing' ,form='unformatted')
        CALL read_colx (x, 'AK', unit_f1d)
        nvclev = SIZE(x)
        DEALLOCATE (vct)
        ALLOCATE (vct (2*nvclev))
        vct(1:nvclev) = x
        CALL read_colx (x, 'BK', unit_f1d)
        vct(nvclev+1:nvclev+nvclev) = x
        DEALLOCATE (x)
      ENDIF
      CALL p_bcast (nvclev, p_io)
      nlev   = nvclev-1
      nlevp1 = nvclev
      IF (p_pe/=p_io) THEN
        DEALLOCATE (vct)
        ALLOCATE (vct (2*nvclev))
      ENDIF
      CALL p_bcast (vct, p_io)            
      IF (p_pe==p_io) THEN
        WRITE(nout,*)
        WRITE(nout,'(a)') REPEAT('-',72)
        WRITE(nout,'(a)')' Column model has changed hybrid levels:'
        WRITE(nout,'(a)')' this option does not work with restart files.'
        WRITE(nout,'(a)')' this option does not work with nudging.'
        WRITE(nout,*)
        WRITE(nout,'(a,i4)')  '      nlev      ='  ,nlev
        WRITE(nout,'(a,i4)')  '      nlevp1    ='  ,nlevp1
        WRITE(nout,'(a,i4)')  '      nvclev    ='  ,nvclev
        WRITE(nout,*)
        IF(lres) &
          CALL finish('inicolumn','lvctch=.true. conflicts with lres=.true.')
      ENDIF
    ENDIF
  END SUBROUTINE inicolumn
!------------------------------------------------------------------------------
  SUBROUTINE setcolumn
  !
  ! Set module variables within this module.
  ! This routine must be called after the parallel decomposition has
  ! been set (after init_decomposition)
  !
    INTEGER :: p, i, j, k, l
    INTEGER :: lon, lat
    REAL    :: tlat(ncol)
    REAL    :: tlon(ncol)
    REAL    :: tpes(ncol)

    p_ful    = -1
    p_col    = -1
    IF(ldc% pe==p_io) THEN
      ALLOCATE (lats  (ncol))  ;lats  = -1
      ALLOCATE (lons  (ncol))  ;lons  = -1
      ALLOCATE (pes   (ncol))  ;pes   = -1
    ELSE
      ALLOCATE (lats  (0))
      ALLOCATE (lons  (0))
      ALLOCATE (pes   (0))
    ENDIF
    iread = 1
    !
    ! loop over pe's
    !
    DO p=1,SIZE(gdc)
      !
      ! set p_col, llat_1d = 1, llon_1d = 1
      !
      IF(gdc(p)% col_1d)THEN
        IF(ncol>1) CALL finish('setcolumn','more than one gp in column model')
        p_col = gdc(p)% pe
        IF(ldc% pe == gdc(p)% pe) THEN
          nlcol    = 1
          ALLOCATE (llats (nlcol)) ;llats = 1
          ALLOCATE (llons (nlcol)) ;llons = 1
          IF(ldc% pe==p_io) THEN
            lats = lat_1d(1)
            lons = lon_1d(1)
          ENDIF
        ENDIF
      ENDIF
      !
      ! set p_ful, llat_1d, llon_1d
      !
      IF(.NOT.gdc(p)% col_1d) THEN
        l = 0
        i=0; j=0; k=0
        DO
         IF (lat_1db(1)/=-1) THEN
           k=MAX(lon_1db(1),k)
           j=MAX(lat_1db(1)-1,j)
           j=j+1
           IF(j>lat_1db(2)) THEN
             j=lat_1db(1)
             k=k+1
             IF(k>lon_1db(2)) EXIT
           ENDIF
           lon=k
           lat=j
         ELSE
           k=k+1
           IF(k>ncol) EXIT
           IF(lon_1d(k)==-1) CYCLE
           lat=lat_1d(k)
           lon=lon_1d(k)
         ENDIF
         i=i+1
         IF(gdc(p)%glats(1)<=lat .AND. gdc(p)% glate(1)>=lat .AND. &
            gdc(p)%glons(1)<=lon .AND. gdc(p)% glone(1)>=lon) THEN
           p_ful   = gdc(p)% pe
           l   = l+1
           tlat(l) = lat - gdc(p)% glats(1) + 1
           tlon(l) = lon - gdc(p)% glons(1) + 1
           tpes(l) = gdc(p)% pe
         ENDIF
         IF(gdc(p)%glats(2)<=lat .AND. gdc(p)% glate(2)>=lat .AND. &
            gdc(p)%glons(2)<=          gdc(p)% glone(2)      .AND. &
            gdc(p)%glons(2)<=lon .AND. gdc(p)% glone(2)>=lon) THEN
           p_ful   = gdc(p)% pe
           l   = l+1
           tlat(l) = lat - gdc(p)% glats(2) + 1
           tlon(l) = lon - gdc(p)% glons(2) + 1
           tpes(l) = gdc(p)% pe
         ENDIF
         IF(gdc(p)%glats(2)<=lat .AND. gdc(p)% glate(2)>=lat .AND. &
            gdc(p)%glons(2)>           gdc(p)% glone(2)      .AND. &
           (gdc(p)%glons(2)<=lon .OR.  gdc(p)% glone(2)>=lon))THEN
           p_ful   = gdc(p)% pe
           l   = l+1
           tlat(l) = lat - gdc(p)% glats(2) + 1
           tlon(l) = lon - gdc(p)% glons(2) + 1
           IF (tlon(l)<1) tlon(l) = tlon(l) + gdc(p)% nlon
           tpes(l) = gdc(p)% pe
         ENDIF
         IF(ldc% pe==p_io) THEN
           lats(i)  = lat
           lons(i)  = lon
           pes (i)  = gdc(p)% pe
         ENDIF
        END DO
        IF (ldc% pe==gdc(p)% pe) THEN
          nlcol = l
          ALLOCATE (llats (nlcol)) ;llats = tlat(:nlcol)
          ALLOCATE (llons (nlcol)) ;llons = tlon(:nlcol)
        ENDIF
      ENDIF
    END DO
    !
    ! allocate, check parameters
    !
    IF (.NOT.ALLOCATED(llats)) ALLOCATE (llats (0))
    IF (.NOT.ALLOCATED(llons)) ALLOCATE (llons (0))
    IF (ldc% col_1d) THEN
      ALLOCATE (ug (ldc%nglon, ldc%nlev, ldc%nglat)); ug = ug_1d
      ALLOCATE (vg (ldc%nglon, ldc%nlev, ldc%nglat)); vg = vg_1d
      k_vg  = MIN (k_vg , ldc%nlev) 
    ENDIF
    !
    ! printout
    !
    IF (p_col == -1 .AND. p_ful == -1) RETURN
    IF (p_pe == p_io) THEN
      WRITE(nout,*)
      WRITE(nout,'(a)') REPEAT('-',72)
      WRITE(nout,'(a)')' Column model initialised:'
      WRITE(nout,*)
      WRITE(nout,'(a,1x,a)')'      comode    =', comode
      WRITE(nout,'(a,i4)')  '      p_col     =', p_col
      WRITE(nout,'(a,i4)')  '      p_ful     =', p_ful
      WRITE(nout,*)
      WRITE(nout,'(a,i4)')  '      ncol      =', ncol
      WRITE(nout,'(a,i4,7i7/(14x,8i7))') &
                            '      jlat      =', lats
      WRITE(nout,'(a,i4,7i7/(14x,8i7))') &
                            '      ilon      =', lons
      WRITE(nout,'(a,8f7.2/(17x,8f7.2))') &
                            '      longitude :',360./ldc%nlon*(lons-1)
      WRITE(nout,'(a,8f7.2/(17x,8f7.2))') &
                            '      approx.lat:',90.-180.*(lats-0.5)/ldc%nlat
      WRITE(nout,*)
      WRITE(nout,'(a,l4)')  '      lrewind   ='  ,lrewind
      WRITE(nout,*)
      WRITE(nout,'(a,f7.2)')'      ug_1d     =', ug_1d
      WRITE(nout,'(a,f7.2)')'      vg_1d     =', vg_1d
      WRITE(nout,'(a,i4)')  '      k_vg      =', k_vg 
      WRITE(nout,*)
      WRITE(nout,'(a,i4)')  '      int_t     =', int_t
      WRITE(nout,'(a,i4)')  '      int_ps    =', int_ps
      WRITE(nout,'(a,i4)')  '      int_uv    =', int_uv
      WRITE(nout,*)
      WRITE(nout,'(a,i4)')  '      for_t     =', for_t
      WRITE(nout,'(a,i4)')  '      for_ps    =', for_ps
      WRITE(nout,'(a,i4)')  '      for_uv    =', for_uv
      WRITE(nout,'(a,i4)')  '      for_q     =', for_q
      WRITE(nout,'(a,i4)')  '      for_x     =', for_x
      WRITE(nout,'(a,i4)')  '      for_xt    =', for_xt
      WRITE(nout,*)
      WRITE(nout,'(a,i4)')  '      for_tr    =', for_tr
      WRITE(nout,'(a,i4)')  '      for_trp   =', for_trp
      WRITE(nout,'(a,i4)')  '      for_d     =', for_d
      WRITE(nout,'(a,i4)')  '      for_vo    =', for_vo
      WRITE(nout,'(a,i4)')  '      for_dt    =', for_dt
      WRITE(nout,'(a,i4)')  '      for_dp    =', for_dp
      WRITE(nout,'(a,i4)')  '      for_zm    =', for_zm
      WRITE(nout,'(a,i4)')  '      for_te    =', for_te
      WRITE(nout,*)
      WRITE(nout,'(a,l4)')  '      lnudge    =', lnudge
      DO i=ml,1,-1
      IF(nudguv(i)==0. .AND. nudgt (i)==0 .AND. nudgq(i)==0 .AND. &
         nudgx (i)==0. .AND. nudgxt(i)==0 .AND.       i > 1       ) CYCLE
      WRITE(nout,'(a,(10f5.2)/(17x,10f5.2))')'      nudguv    =', nudguv(1:i)
      WRITE(nout,'(a,(10f5.2)/(17x,10f5.2))')'      nudgt     =', nudgt (1:i)
      WRITE(nout,'(a,   f5.2              )')'      nudgp     =', nudgp
      WRITE(nout,'(a,(10f5.2)/(17x,10f5.2))')'      nudgq     =', nudgq (1:i)
      WRITE(nout,'(a,(10f5.2)/(17x,10f5.2))')'      nudgx     =', nudgx (1:i)
      WRITE(nout,'(a,(10f5.2)/(17x,10f5.2))')'      nudgxt    =', nudgxt(1:i)
      EXIT
      END DO
      WRITE(nout,'(a)') REPEAT('-',72)
      WRITE(nout,*)
      !
      ! open files
      !
      IF (p_col/=-1 .OR. p_ful/=-1) &
        OPEN (unit_f1d,file='forcing' ,form='unformatted')
      IF (p_col/=-1) THEN
        OPEN (unit_r1d,file='residui' ,form='unformatted')
        OPEN (unit_b1d,file='result1' ,form='unformatted')
        OPEN (unit_a1d,file='result2' ,form='unformatted')
      ENDIF
    ENDIF
  END SUBROUTINE setcolumn
!------------------------------------------------------------------------------
  SUBROUTINE resetcolumn
    !
    ! clean up module variables
    !
    comode = ' '
    lat_1d =  0
    lon_1d =  0
    p_ful  = -1
    p_col  = -1
    ncol   =  0
    nlcol  =  0
    !
    ! close output files
    !
    IF (p_pe == p_io) THEN
      CLOSE (unit_f1d)
      CLOSE (unit_r1d)
      CLOSE (unit_b1d)
      CLOSE (unit_a1d)
    ENDIF
    !
    ! deallocate
    !
    IF (ldc% col_1d) DEALLOCATE (ug, vg)
    IF (ALLOCATED (lats) ) DEALLOCATE (lats)
    IF (ALLOCATED (lons) ) DEALLOCATE (lons)
    IF (ALLOCATED (llats)) DEALLOCATE (llats)
    IF (ALLOCATED (llons)) DEALLOCATE (llons)
    IF (ALLOCATED (pes)  ) DEALLOCATE (pes)    
  END SUBROUTINE resetcolumn
!==============================================================================
! Routines to be called during grid point calculations (scan1sl)
!==============================================================================
  SUBROUTINE get_col_tran (p_dest)
  USE mo_scan_buffer, ONLY: alps_scb, qe_scb, xe_scb, xte_scb, te_scb
  USE mo_memory_g1a,  ONLY: alpsm1
  INTEGER ,INTENT(in) ,OPTIONAL :: p_dest ! destination PE
  !
  ! Prepare or forcing fields for column model.
  ! To be called from scan1sl after the semilagrangian transport
  ! (after CALL slt2)
  !
    INTEGER               :: dest, i
    CHARACTER(len=6)      :: name = 'XTEnnn'
    REAL     ,ALLOCATABLE :: te (:,:,:)
    dest = p_col; IF(PRESENT(p_dest)) dest = p_dest
    !
    ! fetch from other PEs
    !
    IF (comode/=' ') THEN
      IF(for_tr >0) CALL process_col3 (qe_scb,   'QE'      ,for_tr ==1)
      IF(for_tr >0) CALL process_col3 (xe_scb,   'XE'      ,for_tr ==1)
      IF(for_trp>0) CALL process_col2 (alps_scb, 'ALPS_TR' ,for_trp==1)
      IF(for_trp>0) CALL process_col2 (alpsm1,   'ALPSM1'  ,for_trp==1)
      IF(for_tr >0) THEN
        DO i=1,SIZE(xte_scb,3)
          WRITE (name(4:6),'(i3.3)') i
          CALL process_col3 (xte_scb(:,:,i,:), name ,for_tr==1)
        END DO
      ENDIF
      !
      ! allow for an additional temperature tendency term
      !
      IF(for_te>0) THEN
        ALLOCATE (te (SIZE(te_scb,1),SIZE(te_scb,2),SIZE(te_scb,3)))
        te = 0.
        CALL process_col3 (te, 'TE', for_te==1)
        te_scb = te_scb + te
        DEALLOCATE (te)
      ENDIF
    ENDIF
  END SUBROUTINE get_col_tran
!------------------------------------------------------------------------------
  SUBROUTINE get_col_ffti (p_dest)
  USE mo_scan_buffer, ONLY: alps_scb, d_scb, dalpsl_scb, dalpsm_scb,  &
                            dtl_scb, dtm_scb, du0_scb, t_scb, u0_scb, &
                            u_scb, ul_scb, v_scb, vo_scb
  USE mo_memory_gl, ONLY : q, x, xt
  INTEGER ,INTENT(in) ,OPTIONAL :: p_dest ! destination PE
  !
  ! Prepare or get forcing fields for column model.
  ! To be called from scan1sl after the inverse fft 
  ! (after CALL fourier_to_gridpoint).
  !    
    CHARACTER(len=5) :: name = 'XTnnn'
    INTEGER          :: i
    INTEGER          :: dest
    dest = p_col; IF(PRESENT(p_dest)) dest = p_dest
    !
    ! fetch from other PEs
    !
    IF (comode/=' ') THEN
      !
      ! flush files
      !
      IF (p_pe == p_io) THEN
        IF (p_ful/=-1)        CALL flush (unit_f1d)
        IF (p_col/=-1) THEN
          IF(comode=='resid') CALL flush (unit_r1d)
          CALL flush (unit_b1d)
          CALL flush (unit_a1d)
        ENDIF
      ENDIF
      !
      ! rewind forcing, residui files
      !
      IF (lrewind .AND. p_pe == p_io .AND. p_col/=-1 .AND. &
         (comode=='force' .OR. comode=='add')) THEN
        REWIND (unit_f1d)
        REWIND (unit_r1d)
      ENDIF
      !
      ! write nstep, latitude/longitude index
      !
      CALL process_nstep_lonlat
      !
      ! prognostic variables in column model
      !
      IF (for_ps >0) CALL process_col2 ( alps_scb, 'ALPS'   ,for_ps==1 ,nudgp )
      IF (for_t  >0) CALL process_col3 (    t_scb, 'T'      ,for_t ==1 ,nudgt )
      IF (for_uv >0) CALL process_col3 (    u_scb, 'U'      ,for_uv==1 ,nudguv)
      IF (for_uv >0) CALL process_col3 (    v_scb, 'V'      ,for_uv==1 ,nudguv)
      IF (for_q  >0) CALL process_col3 (  q      , 'Q'      ,for_q ==1 ,nudgq )
      IF (for_x  >0) CALL process_col3 (  x      , 'X'      ,for_x ==1 ,nudgx )
      IF (for_xt >0) THEN
        DO i=1,SIZE(xt,3)
          WRITE (name(3:5),'(i3.3)') i
          CALL process_col3 (xt(:,:,i,:), name ,for_xt==1, nudgxt)
        END DO
      ENDIF
      !
      ! geostrophic wind
      !
      IF (k_vg  > 0 .AND. fapp) THEN
        ug (:,:k_vg    ,:) =      u_scb (:,:k_vg ,:)
        vg (:,:k_vg    ,:) =      v_scb (:,:k_vg ,:)
        ug (:, k_vg +1:,:) = SPREAD (ug (:, k_vg ,:), 2, ldc%nlev-k_vg )
        vg (:, k_vg +1:,:) = SPREAD (vg (:, k_vg ,:), 2, ldc%nlev-k_vg )
      ENDIF
      !
      ! diagnostic variables in column model
      !
      IF (for_dt>0) CALL process_col3 (   dtm_scb, 'DTM'    ,for_dt==1)
      IF (for_dt>0) CALL process_col3 (   dtl_scb, 'DTL'    ,for_dt==1)
      IF (for_vo>0) CALL process_col3 (    vo_scb, 'VO'     ,for_vo==1)
      IF (for_d >0) CALL process_col3 (     d_scb, 'D'      ,for_d ==1)
      IF (for_dp>0) CALL process_col2 (dalpsl_scb, 'DALPSL' ,for_dp==1)
      IF (for_dp>0) CALL process_col2 (dalpsm_scb, 'DALPSM' ,for_dp==1)
      !
      ! zonal means
      !
      IF (for_zm>0) CALL process_colzm(  u0_scb, 'U0'     ,for_zm==1)
      IF (for_zm>0) CALL process_colzm(  ul_scb, 'UL'     ,for_zm==1)
      IF (for_zm>0) CALL process_colzm( du0_scb, 'DU0'    ,for_zm==1)
    ENDIF
  END SUBROUTINE get_col_ffti
!------------------------------------------------------------------------------
  SUBROUTINE cal_col_expl (ztodt, jlat)
  !
  ! Perform explicit time integration in a column model run 
  ! for variables treated implicitly in the full model.
  ! To be called from scan1sl before the call to 'si1'
  !
  USE mo_scan_buffer, ONLY: u_scb, v_scb, t_scb, alps_scb,   &! new values
                            d_scb, vo_scb, dtm_scb, dtl_scb, &
                            dalpsl_scb, dalpsm_scb
  USE mo_memory_g2a,  ONLY: um1, vm1                          ! old values
  USE mo_memory_g1a,  ONLY: alpsm1, tm1                       ! old values
  USE mo_sc1,         ONLY: vol, vom, te, alpse               ! tendencies
  USE mo_gaussgrid,   ONLY: coriol                            ! Coriolis param.
  REAL    ,INTENT(in) :: ztodt                                ! time step
  INTEGER ,INTENT(in) :: jlat                                 ! latit. index
    INTEGER :: jglat, irow ! globel indices: north -> south, ping pong
    INTEGER :: ngl         ! number of gaussian latitudes
    REAL    :: zcorio      ! coriolis parameter
    IF(ldc% col_1d) THEN
      !
      ! get coriolis parameter
      !
      ngl   = ldc% nlat                      ! number of gaussian latitudes
      jglat = ldc% glat(jlat)                ! global continuous north -> south
      irow  = MIN(2*jglat-1,2*(ngl+1-jglat)) ! global ping pong index
      zcorio = coriol(irow)
      !
      ! explicit integration for prognostic variables
      !
      IF (int_uv==1) u_scb    (:,:,jlat) = um1   (:,:,jlat) + ztodt *(vom  &
                                    - zcorio * vg(:,:,jlat))
      IF (int_uv==1) v_scb    (:,:,jlat) = vm1   (:,:,jlat) + ztodt *(vol  &
                                    + zcorio * ug(:,:,jlat))
      IF (int_t ==1) t_scb    (:,:,jlat) = tm1   (:,:,jlat) + ztodt * te
      IF (int_ps==1) alps_scb (:,  jlat) = alpsm1(:,  jlat) + ztodt * alpse
      !
      ! zero prescribed terms if flag == 0
      !
      IF (for_d ==0) d_scb      (:,:,jlat) = 0.
      IF (for_vo==0) vo_scb     (:,:,jlat) = 0.
      IF (for_dt==0) dtm_scb    (:,:,jlat) = 0.
      IF (for_dt==0) dtl_scb    (:,:,jlat) = 0.
      IF (for_dp==0) dalpsl_scb (:,  jlat) = 0.
      IF (for_dp==0) dalpsm_scb (:,  jlat) = 0.
    ENDIF
  END SUBROUTINE cal_col_expl
!==============================================================================
! Routines to be called for diagnostics
!==============================================================================
  SUBROUTINE write_column3 (x, name)
  REAL             ,INTENT(inout) :: x (:,:,:)
  CHARACTER(len=*) ,INTENT(in)    :: name
      CALL process_col3 (x, name, .FALSE.)
  END SUBROUTINE write_column3
!------------------------------------------------------------------------------
  SUBROUTINE write_column2 (x, name)
  REAL             ,INTENT(inout) :: x (:,:)
  CHARACTER(len=*) ,INTENT(in)    :: name
      CALL process_col2 (x, name, .FALSE.)
  END SUBROUTINE write_column2
!------------------------------------------------------------------------------
  SUBROUTINE write_column3r (x, name, jlat)
  REAL             ,INTENT(in) :: x (:,:)
  CHARACTER(len=*) ,INTENT(in) :: name
  INTEGER          ,INTENT(in) :: jlat
    REAL :: tmp (SIZE(x,1),SIZE(x,2),ldc%nglat)
    IF (ncol==1) THEN
      IF (ANY(jlat==llats).OR.ANY(jlat==lats)) THEN
        tmp = 0.
        tmp (:,:,jlat) = x
        CALL process_col3 (tmp, name, .FALSE.)
      ENDIF
    ENDIF
  END SUBROUTINE write_column3r
!------------------------------------------------------------------------------
  SUBROUTINE write_column2r (x, name, jlat)
  REAL             ,INTENT(in) :: x (:)
  CHARACTER(len=*) ,INTENT(in) :: name
  INTEGER          ,INTENT(in) :: jlat
    REAL :: tmp (SIZE(x),ldc%nglat)
    IF (ncol==1) THEN
      IF (ANY(jlat==llats).OR.ANY(jlat==lats)) THEN
        tmp = 0.
        tmp (:,jlat) = x
        CALL process_col2 (tmp, name, .FALSE.)
      ENDIF
    ENDIF
  END SUBROUTINE write_column2r
!==============================================================================
! Private routines within this module
!==============================================================================
!------------------------------------------------------------------------------
  SUBROUTINE gather_col3 (co, gp, p_dest)
  REAL    ,INTENT(in)           :: gp (:,:,:) ! local field (lon ,lev,lat)
  REAL    ,POINTER              :: co (:,:,:) ! column      (ncol,lev, 1)
  INTEGER ,INTENT(in)           :: p_dest     ! destination PE
  !
  ! get a column from the full model (distributed domain)
  !
    REAL    :: tmp (ncol,SIZE(gp,2))
    INTEGER :: i
    IF (p_pe == p_dest) THEN
      IF(SIZE(co,1)/=ncol.OR.SIZE(co,3)/=1) THEN
        WRITE(nerr,*)'gather_col3: p_pe, p_dest, p_ful=',p_pe, p_dest, p_ful
        CALL finish ('gather_col3: shape(co)/=(/ncol,:,1/)')
      ENDIF
      IF (p_ful /= p_dest) THEN
        CALL p_recv (co, p_ful, tag_gather_gp)
      ELSE
        DO i=1,nlcol
          co(i,:,1) = gp(llons(i),:,llats(i))
        END DO
      ENDIF 
    ELSE IF (p_pe == p_ful) THEN
      DO i=1,nlcol
        tmp(i,:) = gp(llons(i),:,llats(i))
      END DO
      CALL p_send (tmp, p_dest, tag_gather_gp)
    ENDIF
  END SUBROUTINE gather_col3
!------------------------------------------------------------------------------
  SUBROUTINE gather_col2 (co, gp, p_dest)
  REAL    ,INTENT(in)           :: gp (:,:) ! local field (lon ,lat)
  REAL    ,POINTER              :: co (:,:) ! column      (ncol, 1 )
  INTEGER ,INTENT(in)           :: p_dest   ! destination PE
  !
  ! get surface field from the full model (distributed domain)
  !
    REAL    :: tmp (ncol)
    INTEGER :: i
    IF (p_pe == p_dest) THEN
      IF(SIZE(co,1)/=ncol.OR.SIZE(co,2)/=1) &
        CALL finish ('gather_col2: shape(co)/=(/ncol,1/)')
      IF (p_ful /= p_dest) THEN
        CALL p_recv (co, p_ful, tag_gather_gp)
      ELSE
        DO i=1,nlcol
          co(i,1) = gp(llons(i),llats(i))
        END DO
      ENDIF 
    ELSE IF (p_pe == p_ful) THEN
      DO i=1,nlcol
        tmp(i) = gp(llons(i),llats(i))
      END DO
      CALL p_send (tmp, p_dest, tag_gather_gp)
    ENDIF
  END SUBROUTINE gather_col2
!------------------------------------------------------------------------------
  SUBROUTINE gather_colzm (co, gp, p_dest)
  REAL    ,INTENT(in)           :: gp (:,:) ! local field (lev ,lat)
  REAL    ,POINTER              :: co (:,:) ! column      (ncol,lev)
  INTEGER ,INTENT(in)           :: p_dest   ! destination PE
  !
  ! get zonal mean from the full model (distributed domain)
  !
    REAL    :: tmp (ncol,SIZE(gp,1))
    INTEGER :: i
    IF (p_pe == p_dest) THEN
      IF(SIZE(co,1)/=ncol) CALL finish ('gather_col2: size(co,1)/=ncol')
      IF (p_ful /= p_dest) THEN
        CALL p_recv (co, p_ful, tag_gather_gp)
      ELSE
        DO i=1,nlcol
          co(i,:) = gp(:,llats(i))
        END DO
      ENDIF
    ELSE IF (p_pe == p_ful) THEN
        DO i=1,nlcol
          tmp(i,:) = gp(:,llats(i))
        END DO
      CALL p_send (tmp, p_dest, tag_gather_gp)
    ENDIF
  END SUBROUTINE gather_colzm
!------------------------------------------------------------------------------
  SUBROUTINE process_nstep_lonlat
    !
    ! write time step, longitude/latitude index
    !
    REAL ,ALLOCATABLE :: ijlonlat (:,:)
    REAL              :: rstep
    CHARACTER(len=16) :: nam16
    IF (p_pe==p_io) THEN
      !
      ! write time step
      !
      rstep = nstep
      nam16 = 'NSTEP'
      IF(comode=='traject'.OR.comode=='test') THEN
        WRITE (unit_f1d) nam16, 1, 1, 0, 0
        WRITE (unit_f1d) rstep
      ENDIF
      IF (p_col/=-1) THEN
        IF(comode=='resid') THEN
          WRITE (unit_r1d) nam16, 1, 1, 0, 0
          WRITE (unit_r1d) rstep
        ENDIF
        WRITE (unit_a1d) nam16, 1, 1, 0, 0
        WRITE (unit_a1d) rstep
        WRITE (unit_b1d) nam16, 1, 1, 0, 0
        WRITE (unit_b1d) rstep
      ENDIF
      !
      ! write longitude/latitude index, cvt table, timestep
      !
      IF (nstep==nresum) THEN
        ALLOCATE (ijlonlat (ncol,2))
        ijlonlat (:,1) = lons
        ijlonlat (:,2) = lats
        SELECT CASE (comode)
        CASE ('')
        CASE ('traject','test')
          nam16 = 'ILON_JLAT' 
          WRITE (unit_f1d) nam16, 2, SHAPE(ijlonlat), 0
          WRITE (unit_f1d) ijlonlat
          nam16 = 'DTIME'
          WRITE (unit_f1d) nam16, 2, 1, 1, 0, 0
          WRITE (unit_f1d) dtime
          nam16 = 'AK'
          WRITE (unit_f1d) nam16, 2, 1, nlevp1, 0, 0
          WRITE (unit_f1d) vct(1:nlevp1)
          nam16 = 'BK'
          WRITE (unit_f1d) nam16, 2, 1, nlevp1, 0, 0
          WRITE (unit_f1d) vct(nvclev+1:nvclev+nlevp1)
        CASE default
        
        END SELECT
      ENDIF
    ENDIF
  END SUBROUTINE process_nstep_lonlat
!------------------------------------------------------------------------------
  SUBROUTINE process_col3 (gp, name, lforce, wnudge)
  REAL              ,INTENT(inout) :: gp (:,:,:)
  CHARACTER (len=*) ,INTENT(in)    :: name
  LOGICAL           ,INTENT(in)    :: lforce
  REAL    ,OPTIONAL ,INTENT(in)    :: wnudge(:)
    REAL, POINTER      :: f (:,:,:)
    REAL, POINTER      :: r (:,:,:)
    REAL, POINTER      :: m (:,:,:)
    REAL, POINTER      :: b (:,:,:)
    INTEGER            :: k
    CHARACTER (len=16) :: nam16
    nam16 = name
    !
    ! allocate temporaries
    !
    NULLIFY (f,r,m)
    IF (p_pe==p_io) THEN
      ALLOCATE (f(ncol,SIZE(gp,2),1))
      ALLOCATE (r(ncol,SIZE(gp,2),1))
      ALLOCATE (m(ncol,SIZE(gp,2),1))
      ALLOCATE (b(ncol,SIZE(gp,2),1))
      f = 0.
      r = 0.
    ENDIF
    !
    ! get forcing if full model is running
    !
    IF(p_ful>-1.AND.(p_pe==p_io.OR.p_pe==p_ful)) CALL gather_col3 (f,gp,p_io)
    !
    ! read forcing
    !
    IF (lforce) THEN
      IF(p_pe==p_io) THEN
        SELECT CASE (comode)
        CASE ('force','resid')
          CALL read_col3 (f, name, unit_f1d)
        CASE ('add')
          CALL read_col3 (f, name, unit_f1d)
          CALL read_col3 (r, name, unit_r1d)
        CASE ('free')
          IF (nstep < nresum + nforce) CALL read_col3 (f, name, unit_f1d)
        END SELECT
      ENDIF
    ENDIF
    !
    ! get column model result if column model is running
    !
    IF (p_col>-1.AND.(p_pe==p_io.OR.p_pe==p_col)) THEN
      IF(p_pe==p_io) THEN
        IF(p_io==p_col) THEN
          m = gp
        ELSE
          CALL p_recv (m ,p_col, tag_gather_gp)
        ENDIF
        b = m
      ENDIF
      IF(p_pe==p_col .AND. p_io/=p_col) CALL p_send (gp, p_io, tag_gather_gp)
    ENDIF
    !
    ! process forcing
    !
    fapp = .FALSE.
    IF (lforce) THEN
      SELECT CASE (comode)
      CASE ('resid')
        IF(p_pe==p_io) r = f - m
        IF(p_pe==p_io) m = f
        fapp = .TRUE.
      CASE ('free')
      CASE ('add')
        IF(p_pe==p_io) THEN
          m = m + r
          IF(lnudge .AND. PRESENT(wnudge)) THEN
            DO k=1,SIZE(m,2)
              m(:,k,:) = wnudge(k) * f(:,k,:) + (1.-wnudge(k)) * m(:,k,:)
            END DO
          END IF
        END IF
      CASE default
        IF(p_pe==p_io) m = f
        fapp = .TRUE.
      END SELECT
      IF (nstep < nresum + nforce .AND. p_pe==p_io) m = f
      IF (nstep < nresum + nforce) fapp = .TRUE.
    ENDIF
    !
    ! write forcing, residui and results
    !
    IF(p_pe==p_io) THEN
      SELECT CASE (comode)
      CASE ('traject','test')
        WRITE (unit_f1d) nam16, 3, SHAPE(f)
!WRITE (0,*) nam16, 3, SHAPE(f)
        WRITE (unit_f1d) f
      END SELECT
      SELECT CASE (comode)
      CASE ('resid','test')
        WRITE (unit_r1d) nam16, 3, SHAPE(r)
        WRITE (unit_r1d) r
      END SELECT
      IF (p_col>-1) THEN
        WRITE (unit_b1d) nam16, 3, SHAPE(r)
        WRITE (unit_b1d) b
        WRITE (unit_a1d) nam16, 3, SHAPE(r)
        WRITE (unit_a1d) m
      ENDIF
    ENDIF
    !
    ! apply forcing
    !
    IF (lforce) THEN
      SELECT CASE (comode)
      CASE ('add','traject','resid','force','free')
        IF(p_pe==p_io) THEN
          IF(p_io==p_col) THEN
            gp = m
          ELSE
            CALL p_send (m ,p_col, tag_gather_gp)
          ENDIF
        ENDIF
        IF(p_pe==p_col .AND. p_io/=p_col) CALL p_recv (gp, p_io, tag_gather_gp)
      END SELECT
    ENDIF
    !
    ! deallocate temporaries
    !
    IF (p_pe==p_io)  DEALLOCATE (f,r,m,b)
  END SUBROUTINE process_col3
!------------------------------------------------------------------------------
  SUBROUTINE process_col2 (gp, name, lforce, wnudge)
  REAL              ,INTENT(inout) :: gp (:,:)
  CHARACTER (len=*) ,INTENT(in)    :: name
  LOGICAL           ,INTENT(in)    :: lforce
  REAL    ,OPTIONAL ,INTENT(in)    :: wnudge
    REAL, POINTER      :: f (:,:)
    REAL, POINTER      :: r (:,:)
    REAL, POINTER      :: m (:,:)
    REAL, POINTER      :: b (:,:)
    CHARACTER (len=16) :: nam16
    nam16 = name
    !
    ! allocate temporaries
    !
    NULLIFY (f,r,m)
    IF (p_pe==p_io) THEN
      ALLOCATE (f(ncol,1))
      ALLOCATE (r(ncol,1))
      ALLOCATE (m(ncol,1))
      ALLOCATE (b(ncol,1))
      f = 0.
      r = 0.
    ENDIF
    !
    ! get forcing if full model is running
    !
    IF(p_ful>-1.AND.(p_pe==p_io.OR.p_pe==p_ful)) CALL gather_col2 (f,gp,p_io)
    !
    ! read forcing
    !
    IF (lforce) THEN
      IF(p_pe==p_io) THEN
        SELECT CASE (comode)
        CASE ('force','resid')
          CALL read_col2 (f, name, unit_f1d)
        CASE ('add')
          CALL read_col2 (f, name, unit_f1d)
          CALL read_col2 (r, name, unit_r1d)
        CASE ('free')
          IF (nstep < nresum + nforce) CALL read_col2 (f, name, unit_f1d)
        END SELECT
      ENDIF
    ENDIF
    !
    ! get column model result if column model is running
    !
    IF (p_col>-1.AND.(p_pe==p_io.OR.p_pe==p_col)) THEN
      IF(p_pe==p_io) THEN
        IF(p_io==p_col) THEN
          m = gp
        ELSE
          CALL p_recv (m ,p_col, tag_gather_gp)
        ENDIF
        b = m
      ENDIF
      IF(p_pe==p_col .AND. p_io/=p_col) CALL p_send (gp, p_io, tag_gather_gp)
    ENDIF
    !
    ! process forcing
    !
    IF (lforce) THEN
      IF(p_pe==p_io) THEN
        SELECT CASE (comode)
        CASE ('resid')
          r = f - m
          m = f
        CASE ('free')
        CASE ('add')
          m = m + r
          IF(lnudge .AND. PRESENT(wnudge)) THEN
            m(:,:) = wnudge * f(:,:) + (1.-wnudge) * m(:,:)
          END IF
        CASE default
          m = f
        END SELECT
        IF (nstep < nresum + nforce) m = f
      ENDIF
    ENDIF
    !
    ! write forcing
    !
    IF(p_pe==p_io) THEN
      SELECT CASE (comode)
      CASE ('traject','test')
        WRITE (unit_f1d) nam16, 2, SHAPE(f), 0
        WRITE (unit_f1d) f
      END SELECT
      SELECT CASE (comode)
      CASE ('resid','test')
        WRITE (unit_r1d) nam16, 2, SHAPE(r), 0
        WRITE (unit_r1d) r
      END SELECT
      IF (p_col>-1) THEN
        WRITE (unit_a1d) nam16, 2, SHAPE(r), 0
        WRITE (unit_a1d) m
        WRITE (unit_b1d) nam16, 2, SHAPE(r), 0
        WRITE (unit_b1d) b
      ENDIF
    ENDIF
    !
    ! apply forcing
    !
    IF (lforce) THEN
      SELECT CASE (comode)
      CASE ('add','traject','resid','force','free')
        IF(p_pe==p_io) THEN
          IF(p_io==p_col) THEN
            gp = m
          ELSE
            CALL p_send (m ,p_col, tag_gather_gp)
          ENDIF
        ENDIF
        IF(p_pe==p_col .AND. p_io/=p_col) CALL p_recv (gp, p_io, tag_gather_gp)
      END SELECT
    ENDIF
    !
    ! deallocate temporaries
    !
    IF (p_pe==p_io)  DEALLOCATE (f,r,m,b)
  END SUBROUTINE process_col2
!------------------------------------------------------------------------------
  SUBROUTINE process_colzm (gp, name, lforce)
  REAL              ,INTENT(inout) :: gp (:,:)
  CHARACTER (len=*) ,INTENT(in)    :: name
  LOGICAL           ,INTENT(in)    :: lforce
    REAL, POINTER      :: f (:,:)
    REAL, POINTER      :: r (:,:)
    REAL, POINTER      :: m (:,:)
    REAL, POINTER      :: b (:,:)
    CHARACTER (len=16) :: nam16
    nam16 = name
    !
    ! allocate temporaries
    !
    NULLIFY (f,r,m)
    IF (p_pe==p_io) THEN
      ALLOCATE (f(ncol,SIZE(gp,1)))
      ALLOCATE (r(ncol,SIZE(gp,1)))
      ALLOCATE (m(ncol,SIZE(gp,1)))
      ALLOCATE (b(ncol,SIZE(gp,1)))
      f = 0.
      r = 0.
    ENDIF
    !
    ! get forcing if full model is running
    !
    IF(p_ful>-1.AND.(p_pe==p_io.OR.p_pe==p_ful)) CALL gather_colzm(f,gp,p_io)
    !
    ! read forcing
    !
    IF (lforce) THEN
      IF(p_pe==p_io) THEN
        SELECT CASE (comode)
        CASE ('force','resid')
          CALL read_col2 (f, name, unit_f1d)
        CASE ('add')
          CALL read_col2 (f, name, unit_f1d)
          CALL read_col2 (r, name, unit_r1d)
        CASE ('free')
          IF (nstep < nresum + nforce) CALL read_col2 (f, name, unit_f1d)
        END SELECT
      ENDIF
    ENDIF
    !
    ! get column model result if column model is running
    !
    IF (p_col>-1.AND.(p_pe==p_io.OR.p_pe==p_col)) THEN
      IF(p_pe==p_io) THEN
        IF(p_io==p_col) THEN
          m = gp
        ELSE
          CALL p_recv (m ,p_col, tag_gather_gp)
        ENDIF
        b = m
      ENDIF
      IF(p_pe==p_col .AND. p_io/=p_col) CALL p_send (gp, p_io, tag_gather_gp)
    ENDIF
    !
    ! process forcing
    !
    IF (lforce) THEN
      IF(p_pe==p_io) THEN
        SELECT CASE (comode)
        CASE ('resid')
          r = f - m
          m = f
        CASE ('free')
        CASE ('add')
          m = m + r
        CASE default
          m = f
        END SELECT
        IF (nstep < nresum + nforce) m = f
      ENDIF
    ENDIF
    !
    ! write forcing
    !
    IF(p_pe==p_io) THEN
      SELECT CASE (comode)
      CASE ('traject','test')
        WRITE (unit_f1d) nam16, 2, SHAPE(f), 0
        WRITE (unit_f1d) f
      END SELECT
      SELECT CASE (comode)
      CASE ('resid','test')
        WRITE (unit_r1d) nam16, 2, SHAPE(r), 0
        WRITE (unit_r1d) r
      END SELECT
      IF (p_col>-1) THEN
        WRITE (unit_a1d) nam16, 2, SHAPE(r), 0
        WRITE (unit_a1d) m
        WRITE (unit_b1d) nam16, 2, SHAPE(r), 0
        WRITE (unit_b1d) b
      ENDIF
    ENDIF
    !
    ! apply forcing
    !
    IF (lforce) THEN
      SELECT CASE (comode)
      CASE ('add','traject','resid','force','free')
        IF(p_pe==p_io) THEN
          IF(p_io==p_col) THEN
            gp = m
          ELSE
            CALL p_send (m ,p_col, tag_gather_gp)
          ENDIF
        ENDIF
        IF(p_pe==p_col .AND. p_io/=p_col) CALL p_recv (gp, p_io, tag_gather_gp)
      END SELECT
    ENDIF
    !
    ! deallocate temporaries
    !
    IF (p_pe==p_io)  DEALLOCATE (f,r,m,b)
  END SUBROUTINE process_colzm
!------------------------------------------------------------------------------
  SUBROUTINE read_col3 (x, name, iunit)
  REAL             ,INTENT(out) :: x(:,:,:)
  CHARACTER(len=*) ,INTENT(in)  :: name 
  INTEGER          ,INTENT(in)  :: iunit
    CHARACTER(len=16) :: nam16
    INTEGER           :: rank, shap(3), ios
    REAL ,ALLOCATABLE :: z(:,:,:)
    INTEGER           :: ireadold
!write(0,*)'read_col3:',name
l1: DO
      READ (iunit, iostat=ios) nam16, rank, shap
      IF (ios/=0) EXIT
      IF (nam16=='ILON_JLAT') THEN
        ALLOCATE (z(shap(1),shap(2),1))
        READ (iunit) z
        ireadold = iread
        DO iread = 1,shap(1)
!         WRITE(nout,*)'read_col2: ilon,jlat=',z(iread,:,1)
          IF (z(iread,1,1)==lon_1d(1).AND.z(iread,2,1)==lat_1d(1)) THEN
            IF (ireadold /= iread)                           &
              WRITE(nout,*)'read_col2: iread set to ',iread, &
              ',ilon,jlat=',lon_1d(1),lat_1d(1)
            DEALLOCATE (z)
            CYCLE l1
          ENDIF
        END DO
        WRITE(nerr,*)'read_col2: ilon,jlat=',lon_1d(1),lat_1d(1)
        CALL finish('read_col2','lat/lon-index not found in forcing file')
      ENDIF
      ALLOCATE (z(shap(1),shap(2),shap(3)))
      READ (iunit) z
      IF (nam16==name) THEN
        IF (rank /= 3 .OR. ANY (shap(2:3)/=SHAPE(x(1,:,:)))) &
          CALL finish ('read_col3', 'invalid rank/shape of field '//name)
        IF (shap(1)>1) THEN
          x = z (iread:iread,:,:)
        ELSE
          x = z
        ENDIF
        DEALLOCATE (z)
        RETURN
      ENDIF
      DEALLOCATE (z)
    END DO l1
    CALL finish ('read_col3', 'field not found: '//name)
  END SUBROUTINE read_col3
!------------------------------------------------------------------------------
  SUBROUTINE read_col2 (x, name, iunit)
  REAL             ,INTENT(out) :: x(:,:)
  CHARACTER(len=*) ,INTENT(in)  :: name
  INTEGER          ,INTENT(in)  :: iunit
    CHARACTER(len=16) :: nam16
    INTEGER           :: rank, shap(3), ios
    REAL ,ALLOCATABLE :: z(:,:)
!write(0,*)'read_col3:',name
l1: DO
      READ (iunit, iostat=ios) nam16, rank, shap
      IF (ios/=0) EXIT
      ALLOCATE (z(shap(1),shap(2)))
      READ (iunit) z
      IF (nam16=='ILON_JLAT') THEN
        DO iread = 1,shap(1)
!         WRITE(nout,*)'read_col2: ilon,jlat=',z(iread,:)
          IF (z(iread,1)==lon_1d(1).AND.z(iread,2)==lat_1d(1)) THEN
            WRITE(nout,*)'read_col2: iread set to ',iread,&
                                      ',ilon,jlat=',lon_1d(1),lat_1d(1)
            DEALLOCATE (z)
            CYCLE l1
          ENDIF
        END DO
        WRITE(nerr,*)'read_col2: ilon,jlat=',lon_1d(1),lat_1d(1)
        CALL finish('read_col2','lat/lon-index not found in forcing file')
      ENDIF
      IF (nam16==name) THEN
        IF (rank /= 2 .OR. shap(2)/=SIZE(x,1)) &
          CALL finish ('read_col2', 'invalid rank/shape of field '//name)
        IF (shap(1)>1) THEN
          x = z (iread:iread,:)
        ELSE
          x = z
        ENDIF
        DEALLOCATE (z)
        RETURN
      ENDIF
      DEALLOCATE (z)
    END DO l1
    CALL finish ('read_col2', 'field not found: '//name)
  END SUBROUTINE read_col2
!------------------------------------------------------------------------------
  SUBROUTINE read_colx (x, name, iunit)
  REAL             ,POINTER    :: x(:)
  CHARACTER(len=*) ,INTENT(in) :: name
  INTEGER          ,INTENT(in) :: iunit
    CHARACTER(len=16) :: nam16
    INTEGER           :: rank, shap(3), ios
l1: DO
      IF(ASSOCIATED(x)) DEALLOCATE (x)
      READ (iunit, iostat=ios) nam16, rank, shap
      IF (ios/=0) EXIT
      IF (nam16==name) THEN
        ALLOCATE (x(PRODUCT(shap(1:rank))))
        READ (iunit) x
        RETURN
      ELSE
        READ (iunit)
      ENDIF
    END DO l1
    CALL finish ('read_colx', 'field not found: '//name)      
  END SUBROUTINE read_colx
!------------------------------------------------------------------------------
END MODULE mo_column
