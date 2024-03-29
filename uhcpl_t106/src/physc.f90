!OCL NOALIAS

SUBROUTINE physc

  ! Description:
  !
  ! This subroutine controls the calls to the various
  ! physical subroutines.
  !
  ! Method:
  !
  ! *physc* is called from *gpc*.
  !
  ! Externals:
  ! *geopot*   computes full level geopotentials.
  ! *pres*     computes half level pressures.
  ! *presf*    computes full level pressures.
  ! *vdiff*    vertical exchange of u,v,t,q by turbulence.
  ! *cond*     large scale water phase changes.
  ! *radint*   controls radiation computations.
  ! *radheat*  adds radiation tendencies.
  ! *cucall*   controls mass-flux scheme
  ! *surf*     compute new surface values.
  ! *gwdrag *  gravity wave drag scheme
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, January 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, August 1998, tendency diagnostics
  ! A. Rhodin,   MPI, January 1999, subroutines mofed into modules
  ! U.Niemeier,  MPI, July 1999, change nlon, nlp2 nrow etc for parallelization
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_exception,         ONLY: finish
  USE mo_sc1,               ONLY: alnpr, alpha, alpse, qe, te, vol, vom,    &
                                  vervel, xe, xte
  USE mo_memory_g1a,        ONLY: xm1, tm1, qm1, alpsm1, xtm1
  USE mo_memory_g2a,        ONLY: vm1, um1
  USE mo_memory_g3a,        ONLY: aclcacm, aclcm,aclcovm, aclcvm, ahflm,    &
                                  ahfsm, albedom, albm, alwcvim,            &
                                  aprfluxm, aprlm, auxil1m, auxil2m,        &
                                  az0m,  dew2m,  dsnacm, emtefm,  emterm,   &
                                  evapm,forestm, geospm, glacm, qvim,       &
                                  rgcgnm, runoffm, sclf0m, sclfsm, slmm,    &
                                  snm, snm1m, snmelm, sodifm, srad0m,       &
                                  srad0um, sradsm, sradsum, sraf0m, srafsm, &
                                  t2maxm,  t2minm, tclf0m,  tclfsm,   td3m, &
                                  td3m1m, td4m,  td4m1m,  td5m, td5m1m,     &
                                  tdclm, tdclm1m, tdm, tdm1m, teffm, temp2m,&
                                  tkem, tkem1m, topmaxm, trad0m,            &
                                  tradsm, tradsum,  traf0m, trafsm, trsofm, &
                                  trsolm, tslinm, tsm, tsm1m, tsmaxm,       &
                                  tsminm, tsnm, tsnm1m, tsurfm, u10m,       &
                                  ustar3m, ustrm,  v10m, vdism, vltm,       &
                                  vstrm, wimaxm, wind10m,  wlm,             &
                                  wlm1m,  wsm, wsm1m, wsmxm, drainm, g3m
  USE mo_memory_g3b,        ONLY: aclc, aclcac, aclcov, aclcv,  ahfl,       &
                                  ahfs, alb, albedo, alwcvi, aprc, aprflux, &
                                  aprl, aprs, auxil1,  auxil2, az0,         &
                                  dew2,  dsnac, emtef, emter, evap,         &
                                  qvi, rgcgn, runoff, sclf0, sclfs,         &
                                  seaice, siced, sn, snm1, snmel,           &
                                  sodif, srad0, srad0u, srads, sradsu,      &
                                  sraf0, srafs, t2max, t2min, tclf0, tclfs, &
                                  td, td3, td3m1, td4, td4m1, td5, td5m1,   &
                                  tdcl, tdclm1, tdm1, teff,  temp2, tke,    &
                                  tkem1, topmax, trad0,  trads, tradsu,     &
                                  traf0, trafs, trsof, trsol, ts, tslin,    &
                                  tsm1, tsmax, tsmin, tsn, tsnm1, tsurf,    &
                                  u10, ustar3, ustr, ustrgw, v10, varor,    &
                                  ewov, nsov, nwov, neov,                   &
                                  vdis, vdisgw, vgrat, vlt, vstr,           &
                                  vstrgw, wimax, wind10, wl, wlm1, ws,      &
                                  wsm1, wsmx, drain, g3
  USE mo_tmp_buffer,        ONLY: ahfli, aphm1, aphp1, apm1, app1, cvs, cvw, &
                                  dhfq, dhfqs, dhfqw, dhft, geom1, loglac,   &
                                  loland, qhfl, rsfc, rsfl, srfl, ssfc, ssfl,&
                                  thfl, wlmx, xhfl, xtec
  USE mo_control,           ONLY: lpci, lmidatm, ltdiag, nlev, nlevp1,     &
                                  nresum, nrow, twodt, dtime,    &
                                  ngl, lamip2, lstratiform
  USE mo_gaussgrid,         ONLY: aslm, budw, sqcst
  USE mo_hyb,               ONLY: delb, nlevm1
  USE mo_diagnostics,       ONLY: ldiap
  USE mo_rad_switches,      ONLY: nradfr
  USE mo_param_switches,    ONLY: lcond, lconv
  USE mo_constants,         ONLY: cpd, g, rhoh2o, vtmpc1, vtmpc2
  USE mo_start_dataset,     ONLY: nstart, nstep
  USE mo_tracer,            ONLY: ntrac, xtsink, xttropo
  USE mo_diagnostics_zonal, ONLY: dadconz, dcverz, dcvfrz, dcvgrz, dcvgsz,   &
                                  dcvmsz, dlserz, dlsesz, dlsgrz, dlsgsz,    &
                                  dlsmsz
  USE mo_soil_impl
  USE mo_midatm,            ONLY: cccgwd
  USE mo_field,             ONLY: field2
  USE mo_decomposition,     ONLY: dc => local_decomposition
  USE mo_global_op,         ONLY: sum_zonal

  USE mo_diag_tendency,     ONLY: pdiga
  !  USE m_vdiff,           ONLY: vdiff   ! module subroutine
  !  USE m_surf,            ONLY: surf    ! module subroutine
  !  USE m_cond,            ONLY: cond    ! module subroutine
  !  USE m_cond5,           ONLY: cond5   ! module subroutine
  USE mo_skintem,           ONLY: skintem ! module subroutine
  !  USE m_gwdrag,          ONLY: gwdrag  ! module subroutine
  !  USE m_cucall,          ONLY: cucall  ! module subroutine
  !  USE m_radheat,         ONLY: radheat ! module subroutine
!fu++
  USE mo_stratiform,        ONLY: press_meso
 
  IMPLICIT NONE

  !  Local array bounds
  INTEGER :: nglon, nglpx, nglat, nlof
  !  Local scalars: 
  REAL :: zbud, zcst, zdiagt, zdiagw, zrcst, zsdiss, zsevap, zsheat, zsmelt, &
          zsrain, zsum, ztmst, ztwodt, zdifiz, zrici
  INTEGER :: irow, jrow, jlev, jlon, kfdia, kidia, ktdia, n2lp2, jl, jk

  !  Local arrays: 
  REAL :: zdadc(dc%nglpx), zdpsdt(dc%nglpx), zgeo(dc%nglpx), &
          ztvm1(dc%nglpx,nlev), gorsvar (dc%nglpx)
  INTEGER :: ilab(dc%nglpx,nlev), itropo(dc%nglon), itype(dc%nglpx)

  !  Implicit soil scheme 

  REAL :: zdifav(dc%nglon)     ! ground temperature diffusivity  [m**2/s]
  REAL :: zcapav(dc%nglon)     ! ground volumetric heat capacity [J/m**3/K]
  REAL :: zqsnew(dc%nglon)     ! new surf. sat. spec. humidity
  REAL :: zqsnfl(dc%nglon)     ! sublimation
  REAL :: zqlfl(dc%nglon)      ! interception loss
  REAL :: zqvfl(dc%nglon)      ! transpiration
  REAL :: zqgfl(dc%nglon)      ! bare soil evaporation

  ! Implicit soil scheme

  ! Pointers of variables which are needed in restart files 

  REAL, POINTER :: rcapc(:), rgrhx(:), rcgrd(:,:), rdgrd(:,:),     &
                   rcapcm(:), rgrhxm(:), rcgrdm(:,:), rdgrdm(:,:)

  REAL, POINTER :: g3xc5(:),g3xc5m(:)

  !  External subroutines 
  EXTERNAL geopot, pres, presf, radint, statp, radheat, surf, gwdrag

  !  Intrinsic functions 
  INTRINSIC EXP, MOD, SUM


  !  Executable statements

  !  Local array bounds
  irow  = nrow(1)        ! global ping pong latitude index
  jrow  = nrow(2)        ! local continuous latitude index
  nglon = dc%nglon       ! local number of longitudes
  nglpx = dc%nglpx       ! local number of longitudes allocated
  nglat = dc%nglat       ! local number of latitudes
  nlof  = dc%glon(jrow)  ! local longitude offset to global fields

  !-- 1. Allocate storage

  ALLOCATE (thfl(nglpx))
  ALLOCATE (qhfl(nglpx))
  ALLOCATE (dhft(nglpx))
  ALLOCATE (dhfq(nglpx))
  ALLOCATE (dhfqw(nglpx))
  ALLOCATE (dhfqs(nglpx))
  ALLOCATE (cvs(nglpx))
  ALLOCATE (cvw(nglpx))
  ALLOCATE (wlmx(nglpx))
  ALLOCATE (geom1(nglpx,nlev))
  ALLOCATE (aphm1(nglpx,nlevp1))
  ALLOCATE (apm1(nglpx,nlev))
  ALLOCATE (aphp1(nglpx,nlevp1))
  ALLOCATE (app1(nglpx,nlev))
  ALLOCATE (loland(nglpx))
  ALLOCATE (loglac(nglpx))

  !-- 2. Compute some fields needed by the physical routines

  IF (lpci) THEN
     g3xc5 =>  g3(1)%x(:,1,jrow)
     g3xc5m => g3m(1)%x(:,1,jrow)
  END IF

  IF (lamip2) THEN
     ! --------------------- IMPLICIT SOIL SCHEME ------------------------

     ! Assign pointers of variables which are needed in the restart files

     rcgrd => g3(1)%x(:,1:ngrndmx,jrow)  ! cgrnd coefficient in SOIL2ECH
     rdgrd => g3(2)%x(:,1:ngrndmx,jrow)  ! dgrnd coefficient in SOIL2ECH
     rcapc => g3(3)%x(:,1,jrow)          ! ground heat capacity
     rgrhx => g3(4)%x(:,1,jrow)          ! ground heat flux

     rcgrdm => g3m(1)%x(:,1:ngrndmx,jrow)
     rdgrdm => g3m(2)%x(:,1:ngrndmx,jrow)
     rcapcm => g3m(3)%x(:,1,jrow)
     rgrhxm => g3m(4)%x(:,1,jrow)


     ! Physical constants
     zdifiz = 12.E-07      ! temperature diffusivity of ice  [m**2/s]
     zrici  =  2.09E+06    ! volumetric heat capacity of ice [J/m**3/K]

     ! Computational constants
     ztmst  = twodt
     IF(nstep == nstart) ztmst = 0.5*twodt
     zdiagt = 0.5*twodt
     zdiagw = zdiagt/rhoh2o
  END IF

  !-- 2. Compute some fields needed by the physical routines

  !-- 2.1 Compute virtual temperature at t-dt and set *zgeo* to 0

  DO jlev = 1, nlev
     DO jlon = 1, nglon
        ztvm1(jlon,jlev) = tm1(jlon,jlev,jrow) *    &
                          (1.+vtmpc1*qm1(jlon,jlev,jrow)-xm1(jlon,jlev,jrow))
     END DO
  END DO

  zgeo(:) = 0.

  !-- 2.2 Compute (phi-phis) at t-dt using ln(p) at t

  CALL geopot(geom1,ztvm1,alnpr,alpha,zgeo,nglpx,nglon)

  !-- 2.3 Compute pressure at full and half levels at t-dt

  DO jlon = 1, nglon
     aphm1(jlon,nlevp1) = EXP(alpsm1(jlon,jrow))
  END DO

  CALL pres(aphm1,nglpx,aphm1(1,nlevp1),nglon)

  CALL presf(apm1,nglpx,aphm1,nglon)

  IF (lstratiform) THEN
    do jlon=1,nglon
    do jlev=1,nlev
    press_meso(jlon,jlev,jrow) = apm1(jlon,jlev)
    end do
    end do
  END IF

  !-- 2.4 Compute real winds and wind tendencies

  zrcst = 1./sqcst(irow)
  DO jlev = 1, nlev
     DO jlon = 1, nglon
        um1(jlon,jlev,jrow) = zrcst*um1(jlon,jlev,jrow)
        vm1(jlon,jlev,jrow) = zrcst*vm1(jlon,jlev,jrow)
        vol(jlon,jlev) = zrcst*vol(jlon,jlev)
        vom(jlon,jlev) = zrcst*vom(jlon,jlev)
        IF (lpci) THEN
           xte(jlon,jlev,1)=0.
           xe(jlon,jlev)=0.
        END IF
     END DO
  END DO

  !-- 2.5 Estimate adiabatic conversion of potential energy

  zbud = budw(irow)*(.5*twodt/g)
  zdadc(1:nglon) = 0.

  DO jlon = 1, nglon
     zdpsdt(jlon) = aphm1(jlon,nlevp1)*alpse(jlon)
     zdadc(jlon) = zdadc(jlon) + geospm(jlon,jrow)*zdpsdt(jlon)
  END DO

  DO jlev = 1, nlev
     DO jlon = 1, nglon
        zdadc(jlon) = zdadc(jlon) + (1. + vtmpc2*qm1(jlon,jlev,jrow))*cpd*    &
                     (te(jlon,jlev)*(aphm1(jlon,jlev+1) - aphm1(jlon,jlev)) + &
                      tm1(jlon,jlev,jrow)*delb(jlev)*zdpsdt(jlon))
     END DO
  END DO

  zsum = zbud*SUM(zdadc(1:nglon))

  dadconz(irow) = -zsum

  !-- 2.6 Compute logical mask for land and glacier

  DO jlon = 1, nglon
     loland(jlon) = slmm(jlon,jrow) > 0.5
     loglac(jlon) = glacm(jlon,jrow) > 0.5
  END DO

  !-- 2.7 Prepare grid point diagnostics computations

  zbud   = budw(irow)
  zdiagt = 0.5*twodt
  zdiagw = zdiagt/rhoh2o

  !-- 2.8 Set loop values for physical parameterizations

  kidia = 1
  kfdia = nglon
  ktdia = 1

  ! 2.9 Determine tropopause height and mass budgets for tracers

  IF (ntrac > 0) THEN

     CALL xttropo ( nglon, nglpx, nlev, irow, ntrac,         &
                    tm1(:,:,jrow),    apm1,   aphm1,  geom1, &
                    xtm1(:,:,:,jrow), itropo)

  END IF

  IF (lamip2) THEN

     ! Implicit soil scheme

     ! 2.10 Initialisation of soil temperature scheme

     IF (nstep == nresum) THEN

        IF (.NOT. ALLOCATED(rnatur))   ALLOCATE (rnatur(nglon,ngl))
        IF (.NOT. ALLOCATED(cgrnd))    ALLOCATE (cgrnd(nglon,ngl,ngrndmx))
        IF (.NOT. ALLOCATED(dgrnd))    ALLOCATE (dgrnd(nglon,ngl,ngrndmx))
        IF (.NOT. ALLOCATED(pcapcal))  ALLOCATE (pcapcal(nglon,ngl))
        IF (.NOT. ALLOCATED(pfluxgrd)) ALLOCATE (pfluxgrd(nglon,ngl))
        IF (.NOT. ALLOCATED(ptn))      ALLOCATE (ptn(nglon,ngl,ngrndmx))

        DO jl = 1,nglon
           rnatur(jl,jrow) = 0.0
           IF ( slmm(jl,jrow) > 0.5) rnatur(jl,jrow) = 1.0
           IF (glacm(jl,jrow) > 0.5) rnatur(jl,jrow) = 2.0
        END DO

        ! 2.11 Restart: recover fields of the soil temperature scheme of prev. time step

        IF (nstep /= nstart) THEN
           DO jl = 1,nglon
              ! Land surface or glacier
              IF (rnatur(jl,jrow) > 0.0) THEN
                 pcapcal (jl,jrow) = rcapcm(jl)
                 pfluxgrd(jl,jrow) = rgrhxm(jl)
                 DO jk = 1,ngrndmx
                    cgrnd(jl,jrow,jk) = rcgrdm(jl,jk)
                    dgrnd(jl,jrow,jk) = rdgrdm(jl,jk)
                 END DO
              ELSE
                 pcapcal (jl,jrow) = 0.
                 pfluxgrd(jl,jrow) = 0.
                 DO jk = 1,ngrndmx
                    cgrnd(jl,jrow,jk) = 0.
                    dgrnd(jl,jrow,jk) = 0.
                 END DO
              END IF
           END DO
        END IF

     ENDIF
     IF (nstep == nstart) THEN

        ! 2.12 Prepare ground temperature diffusivity and vol. heat capacity

        DO jl = 1,nglon
           zdifav(jl) = MERGE(zdifiz,sodifm(jl,jrow),loglac(jl))
           zcapav(jl) = MERGE(zrici,rgcgnm(jl,jrow),loglac(jl))
        END DO

        ! 2.13 Initialize the soil temperature profile

        CALL inisoil2ech(jrow, nglon, dtime, tsm1m(1:nglon,jrow), td4m1m(1:nglon,jrow), &
             td5m1m(1:nglon,jrow), tdm1m(1:nglon,jrow), tdclm1m(1:nglon,jrow), &
             zdifav, zcapav, snm1m(1:nglon,jrow))

     END IF
  END IF

  !-- 3. Radiation parameterisation

  !-- 3.1 Compute radiation tendencies

  IF (MOD(nstep,nradfr)==0) THEN

     CALL radint

  END IF

  !-- 3.2 Add radiation tendencies every time step

  ALLOCATE (srfl(nglpx))

  ! prepare first field
  IF (ltdiag) pdiga(1:nglon,:,15,jrow) = pdiga(1:nglon,:,15,jrow) - te(1:nglon,:)

  CALL radheat(kidia,kfdia,nglon,nglpx,nlof,ktdia,nlev,nlevp1,aphm1,apm1,                  &
              emtefm(1,1,jrow),emterm(1,1,jrow),qm1(1,1,jrow),tm1(1,1,jrow),               &
              trsofm(1,1,jrow),trsolm(1,1,jrow),aclcvm(1,jrow),albedom(1,jrow),            &
              albm(1,jrow),sclf0m(1,jrow),sclfsm(1,jrow),srad0m(1,jrow),                   &
              srad0um(1,jrow),sradsm(1,jrow),sradsum(1,jrow),sraf0m(1,jrow),               &
              srafsm(1,jrow),tclf0m(1,jrow),tclfsm(1,jrow),trad0m(1,jrow),                 &
              tradsm(1,jrow),tradsum(1,jrow),traf0m(1,jrow),trafsm(1,jrow),tsm1m(1,jrow),  &
              emtef(1,1,jrow),emter(1,1,jrow),trsof(1,1,jrow),trsol(1,1,jrow),             &
              aclcv(1,jrow),albedo(1,jrow),alb(1,jrow),sclf0(1,jrow),sclfs(1,jrow),        &
              srad0(1,jrow),srad0u(1,jrow),srads(1,jrow),sradsu(1,jrow),                   &
              sraf0(1,jrow),srafs(1,jrow),srfl,tclf0(1,jrow),tclfs(1,jrow), trad0(1,jrow), &
              trads(1,jrow),tradsu(1,jrow),traf0(1,jrow),trafs(1,jrow),te(1,1))

  IF (ltdiag) THEN
     ! store RADHEAT increment
     pdiga(1:nglon,:,15,jrow) = pdiga(1:nglon,:,15,jrow) + te (1:nglon,:)
     ! prepare next fields
     pdiga(1:nglon,:, 3,jrow) = pdiga(1:nglon,:, 3,jrow) - vom(1:nglon,:) / zrcst
     pdiga(1:nglon,:, 8,jrow) = pdiga(1:nglon,:, 8,jrow) - vol(1:nglon,:) / zrcst
     pdiga(1:nglon,:,16,jrow) = pdiga(1:nglon,:,16,jrow) - te (1:nglon,:)
  ENDIF

  !-- 4. Vertical exchange of u,v,t,q by turbulence

  !-- 4.1 Locate and allocate some space

  ALLOCATE (rsfc(nglpx))
  ALLOCATE (ssfc(nglpx))
  ALLOCATE (rsfl(nglpx))
  ALLOCATE (ssfl(nglpx))
  ALLOCATE (xhfl(nglpx))

  ! Compute pressure at full and half levels at t+dt.

  ztwodt = twodt
  IF (nstep==nstart) ztwodt = ztwodt*.5
  DO jlon = 1, nglon
     aphp1(jlon,nlevp1) = EXP(alpsm1(jlon,jrow)+ztwodt*alpse(jlon))
  END DO

  CALL pres(aphp1,nglpx,aphp1(1,nlevp1),nglon)

  CALL presf(app1,nglpx,aphp1,nglon)

  CALL vdiff(kidia,kfdia,nglon,nglpx,ktdia,nlev,nlevm1,nlevp1,ntrac,                     &
             xtm1(:,:,:,jrow),aclcm(:,:,jrow),aphm1,apm1,geom1,qm1(:,:,jrow),            &
             tkem(:,:,jrow),tkem1m(:,:,jrow),tm1(:,:,jrow),um1(:,:,jrow),                &
             vm1(:,:,jrow),xm1(:,:,jrow),emterm(:,:,jrow),ztvm1,loland,ahflm(:,jrow),ahfsm(:,jrow), &
             az0m(:,jrow),dew2m(:,jrow),evapm(:,jrow),forestm(:,jrow),seaice(:,jrow),    &
             snm1m(:,jrow),srfl,temp2m(:,jrow),tsm1m(:,jrow),t2maxm(:,jrow),             &
             t2minm(:,jrow),ustar3m(:,jrow),ustrm(:,jrow),u10m(:,jrow),vdism(:,jrow),    &
             vstrm(:,jrow),v10m(:,jrow),wimaxm(:,jrow),wind10m(:,jrow),wlm1m(:,jrow),    &
             wsm1m(:,jrow),wsmxm(:,jrow),vltm(:,jrow),tke(:,:,jrow),tkem1(:,:,jrow),     &
             itropo,ahfl(:,jrow),ahfs(:,jrow),az0(:,jrow),                               &
             cvs,cvw,dew2(:,jrow),dhfqs,dhfqw,dhft,evap(:,jrow),ts(:,jrow),qhfl,rsfl,    &
             temp2(:,jrow),thfl,t2max(:,jrow),t2min(:,jrow),ustar3(:,jrow),ustr(:,jrow), &
             u10(:,jrow),vdis(:,jrow),vstr(:,jrow),v10(:,jrow),wimax(:,jrow),            &
             wlmx,wind10(:,jrow),xhfl,zqsnew,zqsnfl,zqlfl,zqvfl,zqgfl,xte,vol(1,1),vom(1,1),   &
             qe(1,1),te(1,1),xe(1,1),vgrat(:,jrow))

  IF (ltdiag) THEN
     ! store VDIFF increment
     pdiga(1:nglon,:, 3,jrow)  = pdiga(1:nglon,:, 3,jrow) + vom(1:nglon,:)/zrcst
     pdiga(1:nglon,:, 8,jrow)  = pdiga(1:nglon,:, 8,jrow) + vol(1:nglon,:)/zrcst
     pdiga(1:nglon,:,16,jrow)  = pdiga(1:nglon,:,16,jrow) + te (1:nglon,:)
     ! prepare next fields
     pdiga(1:nglon,:, 4,jrow)  = pdiga(1:nglon,:, 4,jrow) - vom(1:nglon,:)/zrcst
     pdiga(1:nglon,:, 9,jrow)  = pdiga(1:nglon,:, 9,jrow) - vol(1:nglon,:)/zrcst
     pdiga(1:nglon,:,17,jrow)  = pdiga(1:nglon,:,17,jrow) - te (1:nglon,:)
  ENDIF

  IF (ntrac > 0) THEN
     ztmst = twodt
     IF (nstep == nstart) ztmst = 0.5*twodt
     zdiagt = 0.5*twodt

     CALL xtsink ( nglon, nlev, ztmst, xtm1(:,:,:,jrow), xte )
  END IF

  ! Gravity wave drag parameterisation

  IF (lmidatm) THEN
     IF (nstep .EQ. nresum) THEN
        aprflux (:,jrow) = 0. 
        aprfluxm(:,jrow) = 0. 
     ENDIF

     DO jlon=1,nglon
        gorsvar(jlon) = field2(jlon,1,jrow)
     END DO

     CALL cccgwd ( nglon,  ktdia,    nlev, nlevm1,    nlevp1,         nglpx,         &
                   aphm1,            apm1,            geom1,                         &
                   tm1(:,:,jrow),    um1(:,:,jrow),   vm1(:,:,jrow),                 &
                   ewov(:,jrow),     nsov(:,jrow),    nwov(:,jrow),   neov(:,jrow),  &
                   ustrgw(:,jrow),   vstrgw(:,jrow),  vdisgw(:,jrow),                &
                   te,               vol,             vom,                           &
                   aprfluxm(:,jrow), gorsvar )

  ELSE

     CALL gwdrag ( nglon,  ktdia,    nlev, nlevm1,    nlevp1,         nglpx,         &
                   aphm1,            apm1,            geom1,                         &
                   tm1(:,:,jrow),    um1(:,:,jrow),   vm1(:,:,jrow),                 &
                   ewov(:,jrow),     nsov(:,jrow),    nwov(:,jrow),   neov(:,jrow),  &
                   ustrgw(:,jrow),   vstrgw(:,jrow),  vdisgw(:,jrow),                &
                   te,               vol,             vom )

  END IF

  IF (ltdiag) THEN
     ! store GWDRAG increment
     pdiga(1:nglon,:, 4,jrow)  = pdiga(1:nglon,:, 4,jrow) + vom(1:nglon,:)/zrcst
     pdiga(1:nglon,:, 9,jrow)  = pdiga(1:nglon,:, 9,jrow) + vol(1:nglon,:)/zrcst
     pdiga(1:nglon,:,17,jrow)  = pdiga(1:nglon,:,17,jrow) + te (1:nglon,:)
     ! prepare next fields
     pdiga(1:nglon,:, 5,jrow)  = pdiga(1:nglon,:, 5,jrow) - vom(1:nglon,:)/zrcst
     pdiga(1:nglon,:,10,jrow)  = pdiga(1:nglon,:,10,jrow) - vol(1:nglon,:)/zrcst
     pdiga(1:nglon,:,18,jrow)  = pdiga(1:nglon,:,18,jrow) - te (1:nglon,:)
  ENDIF

  !-- 5. Convection parameterisation

  !-- 5.1 Allocate some space

  ALLOCATE (xtec(nglpx,nlev))

  !-- 5.2 Compute *t* and *q* tendencies by moist convection
  !       and shallow convection.

  itype(:) = 0

  !-- 5.2.1 Initialize arrays for convective precipitation
  !         and copy arrays for convective cloud parameters
  !         -----------------------------------------------

  xtec(:,:) = 0.

  DO jlon = 1, nglon
     rsfc(jlon) = 0.
     ssfc(jlon) = 0.
  END DO

  !-- 5.2.2 Call subroutine cucall for cumulus parameterization

  IF (lconv) THEN

     zsdiss = 0.
     zsrain = 0.
     zsevap = 0.
     zsheat = 0.
     zsmelt = 0.

     n2lp2 = nglpx

     CALL cucall(nglpx,n2lp2,nglon,nlev,nlevp1,nlevm1,ilab,ntrac,xtm1(:,:,:,jrow), &
                 xte,tm1(:,:,jrow),qm1(:,:,jrow),um1(:,:,jrow),                    &
                 vm1(:,:,jrow),xm1(:,:,jrow),te(1,1),qe(1,1),                      &
                 vom(1,1),vol(1,1),xe(1,1),vervel(1,1),qhfl,                       &
                 xtec,app1,aphp1,geom1,rsfc,ssfc,aprc(:,jrow),aprs(:,jrow),        &
                 itype,loland,topmax(:,jrow),topmaxm(:,jrow),                      &
                 zsrain,zsevap,zsheat,zsdiss,zsmelt)

     !-- 5.2.3 Do zonal budget for convection

     dcvgrz(irow) = zdiagw*zbud*zsrain
     dcverz(irow) = zdiagw*zbud*zsevap
     dcvgsz(irow) = 0.
     dcvmsz(irow) = zdiagt*zbud*zsmelt
     dcvfrz(irow) = -zdiagt*zbud*zsdiss

  ELSE
     ! Necessary computations if massflux is by-passed

     dcvgrz(irow) = 0.
     dcverz(irow) = 0.
     dcvgsz(irow) = 0.
     dcvmsz(irow) = 0.

     dcvfrz(irow) = 0.

     ilab(1:nglon,1:nlev) = 0

  END IF

  IF (ltdiag) THEN
     ! store CUCALL increment (massflux)
     pdiga(1:nglon,:, 5,jrow)  = pdiga(1:nglon,:, 5,jrow) + vom(1:nglon,:)/zrcst
     pdiga(1:nglon,:,10,jrow)  = pdiga(1:nglon,:,10,jrow) + vol(1:nglon,:)/zrcst
     pdiga(1:nglon,:,18,jrow)  = pdiga(1:nglon,:,18,jrow) + te (1:nglon,:)
     ! prepare next fields
     pdiga(1:nglon,:,19,jrow)  = pdiga(1:nglon,:,19,jrow) - te (1:nglon,:)
  ENDIF

  !-- 6. Large scale condensation

  IF (lcond) THEN
     IF (lpci) THEN
        IF (ntrac > 0) THEN
           CALL cond5(kidia,kfdia,nglon,nglpx,ktdia,nlev,nlevp1,ilab,aclcacm,aphp1,         &
                      apm1,app1,qm1(:,:,jrow),tm1(:,:,jrow),xm1(:,:,jrow),xtm1(:,:,1,jrow), &
                      xtec,loland,itype,aclcovm,aprlm,                                      &
                      qvim, alwcvim, g3xc5m, aclc, aclcac, aclcov,aprl,qvi, alwcvi, g3xc5,  &
                      ssfl,qe,te,xe,xte(:,:,1),aprs(:,jrow),rsfl)
        ELSE
           CALL finish('physc','ntrac is not set in namelist TRACTL')
        END IF
     ELSE
        CALL cond(kidia,kfdia,nglon,nglpx,ktdia,nlev,nlevp1,ilab,aclcacm(:,:,jrow),      &
                  aphm1,aphp1,apm1,app1,geom1,qm1(:,:,jrow),tm1(:,:,jrow),xm1(:,:,jrow), &
                  xtec,loland,itype,aclcovm(:,jrow),alwcvim(:,jrow),aprlm(:,jrow),       &
                  qvim(:,jrow),aclc(:,:,jrow),aclcac(:,:,jrow),aclcov(:,jrow),           &
                  alwcvi(:,jrow),aprl(:,jrow),qvi(:,jrow),ssfl,te(1,1),                  &
                  qe(1,1),xe(1,1),aprs(:,jrow),rsfl)
     ENDIF

  ELSE

     ! Necessary computations if *cond* is by-passed.

     DO jlon = kidia, kfdia
        ssfl(jlon) = 0.
     END DO

     DO jlev = ktdia, nlev
        DO jlon = kidia, kfdia
           aclc(jlon,jlev,jrow) = 0.
        END DO
     END DO

     dlsgrz(irow) = 0.
     dlsgsz(irow) = 0.
     dlsmsz(irow) = 0.
     dlserz(irow) = 0.
     dlsesz(irow) = 0.
  END IF

  ! store COND increment
  IF (ltdiag) pdiga(1:nglon,:, 19,jrow)  = pdiga(1:nglon,:, 19,jrow) + te(1:nglon,:)

  IF (lmidatm) THEN
     DO jlon=1,nglon
        aprflux(jlon,jrow)=rsfl(jlon)+ssfl(jlon)+rsfc(jlon)+ssfc(jlon)
     END DO
  END IF

  ! Return workspace.

  DEALLOCATE (xtec)

  !-- 7. Computation of new surface values

  !-- 7.1 Compute number of land points on each latitude line

  ! now being set in ioinitial ... what to do for restart ...

  IF (nstep==nresum) THEN
     aslm(irow) = sum_zonal(slmm(1:nglon,jrow),jrow)
  END IF

  !-- 7.2 Call *surf*

  IF (ldiap) CALL statp

  CALL surf     (                 &
       nglon, nglpx, nlevp1,             &
       wsmx(1,jrow) , wsmxm(1,jrow),                   &
       rgcgn(1,jrow), rgcgnm(1,jrow),                  &
       sodif(1,jrow), sodifm(1,jrow),                  &
       vlt(1,jrow), vltm(1,jrow),                      &
       loland, loglac, slmm(1,jrow) ,          &
       ts(1,jrow), tsm(1,jrow), tsm1(1,jrow), tsm1m(1,jrow),  & ! -- surf temp and moist fields  --
       ws(1,jrow), wsm(1,jrow), wsm1(1,jrow), wsm1m(1,jrow),  &
       sn(1,jrow), snm(1,jrow), snm1(1,jrow), snm1m(1,jrow),  &
       wl(1,jrow), wlm(1,jrow), wlm1(1,jrow), wlm1m(1,jrow),  &
       td(1,jrow), tdm(1,jrow), tdm1(1,jrow), tdm1m(1,jrow),  &
       drain(1,jrow), drainm(1,jrow),  &
       srfl, thfl, qhfl, xhfl,         & ! -- surface fluxes  --
       rsfc, ssfc, rsfl, ssfl,         &
       ahfli,                          &
       dhft, dhfqw, dhfqs,             & ! -- flux derivatives --
       evap(1,jrow), evapm(1,jrow),            & ! -- evaporation ---
       vgrat(1,jrow), cvs, wlmx,               & ! -- vegetation  --
       emterm(1,1,jrow),                       & ! -- emissivities --
       dsnac(1,jrow), dsnacm(1,jrow),          & ! -- remaining elements in 'surf' --
       runoff(1,jrow), runoffm(1,jrow),                &
       snmel(1,jrow), snmelm(1,jrow),                  &
       tdcl(1,jrow) , tdclm(1,jrow) , tdclm1(1,jrow), tdclm1m(1,jrow), &
       td3(1,jrow)  , td3m(1,jrow)  , td3m1(1,jrow) , td3m1m(1,jrow) , &
       td4(1,jrow)  , td4m(1,jrow)  , td4m1(1,jrow) , td4m1m(1,jrow) , &
       td5(1,jrow)  , td5m(1,jrow)  , td5m1(1,jrow) , td5m1m(1,jrow) , &
       tslin(1,jrow), tslinm(1,jrow),                  &
       tsmax(1,jrow), tsmaxm(1,jrow),                  &
       tsmin(1,jrow), tsminm(1,jrow),                  &
       tsn(1,jrow)  , tsnm(1,jrow)  , tsnm1(1,jrow) , tsnm1m(1,jrow) , &
       tsurf(1,jrow), tsurfm(1,jrow), varor(1,jrow)            )

  ! Unlock space (previously done in surf)

  DEALLOCATE (qhfl)
  DEALLOCATE (dhfqw)
  DEALLOCATE (dhfqs)
  DEALLOCATE (cvs)
  DEALLOCATE (cvw)
  DEALLOCATE (wlmx)
  DEALLOCATE (rsfl)
  DEALLOCATE (ssfl)

  !-- 7.2.1 Call for skintem

  CALL skintem (nglpx, nlevp1, nglon, albedom(:,jrow), emterm(:,:,jrow), &
       seaice(:,jrow), siced(:,jrow), slmm(:,jrow),                      &
       teff(:,jrow), teffm(:,jrow), tsm1m(:,jrow), dhft, srfl, thfl,     &
       tskin=auxil1(:,jrow), tskinm=auxil1m(:,jrow),                     &
       tskin1=auxil2(:,jrow), tskin1m=auxil2m(:,jrow))

  ! Unlock space (previously done in skintem)

  DEALLOCATE (thfl)
  DEALLOCATE (dhft)
  DEALLOCATE (srfl)

  !-- 7.3 Unlock space

  DEALLOCATE (dhfq)
  DEALLOCATE (rsfc)
  DEALLOCATE (ssfc)
  DEALLOCATE (ahfli)
  DEALLOCATE (xhfl)
  DEALLOCATE (geom1)
  DEALLOCATE (aphm1)
  DEALLOCATE (apm1)
  DEALLOCATE (aphp1)
  DEALLOCATE (app1)

  !-- 8. Restore winds and wind tendencies

  zcst = sqcst(irow)
  DO jlev = 1, nlev
     DO jlon = 1, nglon
        um1(jlon,jlev,jrow) = zcst*um1(jlon,jlev,jrow)
        vm1(jlon,jlev,jrow) = zcst*vm1(jlon,jlev,jrow)
        vol(jlon,jlev) = zcst*vol(jlon,jlev)
        vom(jlon,jlev) = zcst*vom(jlon,jlev)
     END DO
  END DO

  ! 9. -- Prepare soil temperature variables for restart

  IF (lamip2) THEN
     DO jl = 1,nglon
        rcapc(jl) = pcapcal(jl,jrow)
        rgrhx(jl) = pfluxgrd(jl,jrow)
        DO jk = 1,ngrndmx
           rcgrd(jl,jk) = cgrnd(jl,jrow,jk)
           rdgrd(jl,jk) = dgrnd(jl,jrow,jk)
        END DO
     END DO
  END IF

  !-- 10. Release space

  DEALLOCATE (loland)
  DEALLOCATE (loglac)

END SUBROUTINE physc
