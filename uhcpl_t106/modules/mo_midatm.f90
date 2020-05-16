!#define G3X

MODULE mo_midatm

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: shigh     ! (nkp1,nlev)                         used in hdiff
  PUBLIC :: damhih    !                                             hdiff, setdyn
  PUBLIC :: noz       ! number of vertical levels of ozone distrib. radint
  PUBLIC :: scloz     ! (noz)                                       radint
  PUBLIC :: scloza    ! (nlon,nlev)                                 radint
  PUBLIC :: spe       ! (noz)                                       radint
  PUBLIC :: ozone     ! (noz,ngl,0:13)                              radint
  PUBLIC :: spdrag    ! upper sponge layer coefficient (sec)-1      setdyn
  PUBLIC :: enspodi   ! factor from one level to next level         setdyn
  PUBLIC :: nlvspd1   ! last (uppermost) layer of upper sponge      setdyn
  PUBLIC :: nlvspd2   ! first (lowest) layer of upper sponge        setdyn
  PUBLIC :: cccgwd    ! subroutine - gravity wave param.            physc
  PUBLIC :: uspnge    ! subroutine - relaxation on zonal waves      stepon
  PUBLIC :: readozone ! subroutine - reads climat. ozone distrib.   control

  REAL, ALLOCATABLE  :: ozone(:,:,:)
  ! noz, number of vertical levels of ozone distribution 
  ! nioz=unit number associated WITH ozone file
  INTEGER :: noz
  INTEGER, PARAMETER :: nioz = 21

  REAL               :: salp(66), sfaes(21), sfael(21)
  REAL, ALLOCATABLE  :: scloz(:)
  REAL, ALLOCATABLE  :: spe(:)
  REAL, ALLOCATABLE  :: scloza(:,:)
  REAL, ALLOCATABLE  :: shigh(:,:)

  REAL :: damhih
  REAL :: spdrag, enspodi
  ! enspodi  factor by which   upper sponge layer 
  !          coefficient is increased from one        
  !          level to next level above. 
  ! spdrag   upper sponge layer coefficient (sec)-1

  INTEGER :: nlvspd1, nlvspd2
  ! nlvspd1  last (uppermost) layer of upper sponge  
  ! nlvspd2  first (lowest) layer of upper sponge  

CONTAINS

  SUBROUTINE readozone

    !
    !  readozone - reads zonal climatological ozone distribution
    !
    !  U. Schulzweida, MPI, Oct 1999
    !

    USE mo_control,       ONLY: ngl, nkp1, nlev
    USE mo_doctor,        ONLY: nout, nerr
    USE mo_mpi,           ONLY: p_pe, p_io, p_bcast
    USE mo_decomposition, ONLY: dc => local_decomposition
    USE mo_exception,     ONLY: finish
    USE mo_io,            ONLY: IO_open_unit, IO_close, IO_READ, &
                                IO_var_id, IO_file_id, ini_ozon
    USE mo_netCDF,        ONLY: IO_inq_dimid, IO_inq_dimlen,     &
                                IO_inq_varid, IO_get_var_double, &
                                IO_get_vara_double

    INTEGER  :: jk, io_nlon, io_ngl
    INTEGER  :: start(4), COUNT(4)


    ! ALLOCATE memory 
    IF ( .NOT. ALLOCATED(scloza)) ALLOCATE(scloza(dc%nglon,nlev))
    IF ( .NOT. ALLOCATED(shigh))  ALLOCATE(shigh(nkp1,nlev))

    IF (p_pe==p_io) THEN

       ! Read OZONE file

       WRITE(nout,'(/,A,I2)') ' Read OZONE data from unit ', nioz

       CALL IO_open_unit (nioz, ini_ozon, IO_READ)
       IO_file_id = ini_ozon%nc_file_id

       ! Check resolution
       CALL IO_inq_dimid  (IO_file_id, 'lat', io_var_id)
       CALL IO_inq_dimlen (IO_file_id, io_var_id, io_ngl)
       CALL IO_inq_dimid  (IO_file_id, 'lon', io_var_id)
       CALL IO_inq_dimlen (IO_file_id, io_var_id, io_nlon)
       CALL IO_inq_dimid  (IO_file_id, 'level', io_var_id)
       CALL IO_inq_dimlen (IO_file_id, io_var_id, noz)

       IF (io_nlon /= 1 .OR. io_ngl /= ngl) THEN
          WRITE(nerr,*) 'readozone: unexpected resolution ', io_nlon, io_ngl
          CALL finish ('readozone','unexpected resolution')
       END IF
    END IF

    CALL p_bcast (noz, p_io)

    ! ALLOCATE memory 
    IF ( .NOT. ALLOCATED(ozone))  ALLOCATE(ozone(noz,ngl,0:13))
    IF ( .NOT. ALLOCATED(spe))    ALLOCATE(spe(noz))
    IF ( .NOT. ALLOCATED(scloz))  ALLOCATE(scloz(noz))

    IF (p_pe == p_io) THEN

       CALL IO_inq_varid      (IO_file_id, 'level', io_var_id)
       CALL IO_get_var_double (IO_file_id, io_var_id, spe(:))

       CALL IO_inq_varid (IO_file_id, 'OZON', io_var_id)
       DO jk = 1, noz
          COUNT(:) = (/ 1, ngl,  1, 12 /)
          start(:) = (/ 1,   1, jk,  1 /)
          CALL IO_get_vara_double (IO_file_id,io_var_id,start,count,ozone(jk,:,1:12))
       END DO

       CALL IO_close (ini_ozon)

       WRITE(nout,'(I6,A)')  noz, ' level: '
       WRITE(nout,'(10F8.0)')  spe(1:noz)

       ozone(:,:,0)  = ozone(:,:,12)
       ozone(:,:,13) = ozone(:,:,1)

    ENDIF

    CALL p_bcast (spe, p_io)
    CALL p_bcast (ozone, p_io)

  END SUBROUTINE readozone

  SUBROUTINE uspnge
    !
    !****   *uspnge* upper sponge for divergence and vorticity 
    !
    !   E. Manzini, MPI-HH, 1995, original version
    !   T. Diehl, DKRZ, July 1999, parallel version
    !
    !   purpose
    !   ------
    !   
    !   sponge layer for vorticity and divergence:
    !   upper layer(s) linear relaxation on zonal waves 
    !
    !**   interface
    !     ---------
    !   *uspnge* is called from *stepon*
    !

    USE mo_decomposition, ONLY: lc => local_decomposition
    USE mo_control,       ONLY: nlev, twodt
    USE mo_start_dataset, ONLY: nstart, nstep
    USE mo_memory_sp,     ONLY: sd, stp, svo

    REAL    :: zlf(nlev), zs(nlev,2)
    REAL    :: ztwodt, zspdrag, znul
    INTEGER :: jk, is, jlev, jr, snsp
    INTEGER :: mymsp(lc%snsp)

    snsp = lc%snsp
    mymsp = lc%mymsp

    ! 1.  compute linear relaxation
    !     ------- ------ ----------
    ztwodt = twodt
    IF(nstep.EQ.nstart) ztwodt = 0.5*twodt
    zspdrag = spdrag

    znul = 0.0
    DO jk = 1,nlev
       zlf(jk)=znul
    END DO

    zlf(nlvspd2) = zspdrag

    DO jk = nlvspd2-1,nlvspd1,-1
       zspdrag = zspdrag*enspodi
       zlf(jk) = zspdrag
    END DO

    DO jk = 1,nlev
       zs(jk,1) = 1./(1.+zlf(jk)*ztwodt)
       zs(jk,2) = 1./(1.+zlf(jk)*ztwodt)
    END DO

    ! 3.  modify divergence and vorticity         
    !     ------ ---------- --- ---------         

    DO is = 1, snsp
       IF (mymsp(is) /= 0) THEN
          !          DO jlev = 1, nlev                                               
          sd (:,:,is)  = sd(:,:,is)*zs(:,:)
          svo(:,:,is) = svo(:,:,is)*zs(:,:)                   
          !          END DO

          DO jr = 1, 2
             DO jlev = 1, nlev
                stp(jlev,jr,is) = stp(jlev,jr,is)*zs(jlev,1)
             END DO
          END DO

       END IF
    END DO

  END SUBROUTINE uspnge

  SUBROUTINE cccgwd ( klon, ktdia, klev, klevm1, klevp1,   klp2,  &
                      paphm1,      papm1,        pgeom1,          &
                      ptm1,        pum1,         pvm1,            &
                      pewov,       pnsov,        pnwov,    pneov, &
                      pustrgw,     pvstrgw,      pvdisgw,         &
                      ptte,        pvol,         pvom,            &
                      paprfluxm,   porsvar )

    !
    !**** *cccgwd - Replaces the ECMWF gravity wave parameterisation
    !               with the ccc/mam scheme.            
    !
    !     Subject.
    !     --------
    !
    !          This routine computes the physical tendencies of the
    !          prognostic variables u,v  and t due to vertical transports by
    !          subgridscale  gravity waves
    !
    !
    !**   Interface.
    !     ----------
    !
    !          *cccgwd* is called from *physc*.
    !
    !     Method.
    !     -------
    !         The scheme consists of two parts,the calculation of gravity
    !         wave stress and the stress profile in the vertical.
    !         the stress is computed using a low-level wind,static stability
    !         and an orographic variance.four components of variance are
    !         available,the choice determined by the wind direction.
    !         a wave richardson number is computed at each level and by
    !         requiring that its value is never less than a critical one
    !         a value of stress is determined at each model level.
    !
    !         The gravity wave stress comprises two parts depending on
    !         the critical froude number of the low level flow , and an
    !         orographic anisotropy function (a measure of two-dim) is
    !         used for the supercritical component.
    !
    !     Reference.
    !     ----------
    !         See model documentation
    !
    !     Author.
    !     -------
    !
    !      N. Mcfarlane,   DKRZ, May 1995
    !      E. Manzini      MPI,  May 1997,   modified
    !      U. Schulzweida, MPI,  March 2000, use unpack orographic variances
    !      
    !
    !      Based on a combination of the orographic scheme by N.Mcfarlane 1987
    !      and the hines scheme as coded by C. Mclandress 1995.                       
    !

    USE mo_control,           ONLY: nlat, nresum, nrow, twodt
    USE mo_constants,         ONLY: api, cpd, g, rd, rhoh2o, dayl
    USE mo_start_dataset,     ONLY: nstart, nstep
    USE mo_param_switches,    ONLY: lgwdrag
    USE mo_gaussgrid,         ONLY: twomu
#ifdef G3X
    USE mo_memory_g3a,        ONLY: g3m
    USE mo_memory_g3b,        ONLY: g3
#endif

    INTEGER ,INTENT(IN) :: klevm1, ktdia, klev, klevp1, klon, klp2

    !  Array arguments with intent(In):
    ! Input 1D
    REAL, TARGET, INTENT (IN) :: pewov(klp2) !  E-W  orographic variance
    REAL, TARGET, INTENT (IN) :: pnsov(klp2) !  N-S  orographic variance
    REAL, TARGET, INTENT (IN) :: pnwov(klp2) ! NW-SE orographic variance
    REAL, TARGET, INTENT (IN) :: pneov(klp2) ! NE-SW orographic variance
    REAL ,INTENT (IN) :: porsvar   (klp2)
    REAL ,INTENT (IN) :: paprfluxm (klp2)
    ! Input 2D
    REAL, INTENT (IN) :: paphm1(klp2,klevp1) ! half level pressure (t-dt)
    REAL, INTENT (IN) :: papm1(klp2,klev)    ! full level pressure (t-dt)
    REAL, INTENT (IN) :: pgeom1(klp2,klev)   ! geopotential above surface (t-dt)
    REAL, INTENT (IN) :: ptm1(klp2,klev)     ! temperature (t-dt)
    REAL, INTENT (IN) :: pum1(klp2,klev)     ! zonal wind (t-dt)
    REAL, INTENT (IN) :: pvm1(klp2,klev)     ! meridional wind (t-dt)

    !  Array arguments with intent(InOut):
    ! Input 1D
    REAL, INTENT (INOUT) :: pustrgw(klp2)    ! u-gravity wave stress (accumulated, new value)
    REAL, INTENT (INOUT) :: pvstrgw(klp2)    ! v-gravity wave stress (accumulated, new value)
    REAL, INTENT (INOUT) :: pvdisgw(klp2)    ! dissipation by gravity wave drag (accumulated, new value)
    ! Input 2D
    REAL, INTENT (INOUT) :: ptte(klp2,klev)  ! tendency of temperature
    REAL, INTENT (INOUT) :: pvol(klp2,klev)  ! tendency of meridional wind
    REAL, INTENT (INOUT) :: pvom(klp2,klev)  ! tendency of zonal wind

    !  Temporary arrays for ccc/mam gwd scheme

    !  Vertical positioning arrays.                                       

    REAL sgj(klon,klev),     shj(klon,klev),  &        
         shxkj(klon,klev),   dsgj(klon,klev),    th(klon,klev)

    !  Work arrays.                                
    REAL bvfreq(klon,klev),  veln(klon,klev),   &           
         ub(klon),           vb(klon),          &           
         vmodb(klon),        &           
         depfac(klon,klev),  hitesq(klon),      &           
         ampbsq(klon),       denfac(klon),      &            
         rmswind(klon), density(klon,klev),     &          
         utendgw(klon,klev), vtendgw(klon,klev), env(klon),        &            
         alt(klon,klev),     &             
         pressg(klon), visc_mol(klon,klev), anisof(klon),  &
         zpr(klon), theta(klon)

    REAL :: zvar(klon,4)

    !  Local scalars: 
    INTEGER       :: ilat, ilevh, ilevm1, ilevm2, irow, isector, isnorm, jk, jl, lref,    &
                     lrefp, nazmth, jrow
    REAL          :: rgocp, zargt1, zargt2, zcons1, zcons2, zcons4, zcons5, zdiagt, zlat, &
                     zpcons, ztheta, ztmst, zzpos, coslat
#ifdef G3X
    REAL, POINTER :: g3x01(:,:), g3x02(:,:), g3xm01(:,:), g3xm02(:,:) 
#endif

    ! Computational constants.

    ilevm2 = klev-2
    ilevm1 = klev-1
    ilevh  = klev/2
    ztmst  = twodt
    IF (nstep.EQ.nstart) ztmst=0.5*twodt
    zdiagt = 0.5*twodt

    zcons1 = 1./rd
    zcons2 = g**2/cpd
    zcons4 = 1./(g*ztmst)
    zcons5 = 1.5*api
    nazmth = 8
    zpcons = (1000.*dayl)/rhoh2o

    irow = nrow(1)
    ilat = nlat(1)

    jrow = nrow(2)

#ifdef G3X
    g3x01  =>  g3(1)%x(:,1:klev,jrow)
    g3xm01 => g3m(1)%x(:,1:klev,jrow)
    g3x02  =>  g3(2)%x(:,1:klev,jrow)
    g3xm02 => g3m(2)%x(:,1:klev,jrow)

    !  0. set to zero g3x's

    IF (nstep .EQ. nresum) THEN
       g3x01(:,:)  = 0.
       g3xm02(:,:) = 0.
       g3x01(:,:)  = 0.
       g3xm02(:,:) = 0.
    ENDIF
    g3x01  = 0.
    g3x02  = 0.
#endif

    IF (lgwdrag) THEN

       !-- 1. Initialize work array with directional orographic variances

       zvar(:,3) = pewov(1:klon)
       zvar(:,1) = pnsov(1:klon)
       zvar(:,2) = pneov(1:klon)
       zvar(:,4) = pnwov(1:klon)

       !  define constants and arrays needed for the ccc/mam gwd scheme

       !  Constants:
       rgocp = rd/cpd
       lrefp = klevm1
       lref  = lrefp-1

       !  Arrays
       DO jk = ktdia, klev
          DO jl = 1, klon
             shj(jl,jk)   =  papm1(jl,jk) / paphm1(jl,klevp1)
             sgj(jl,jk)   =  papm1(jl,jk) / paphm1(jl,klevp1)
             dsgj(jl,jk)  = (paphm1(jl,jk+1)-paphm1(jl,jk)) / paphm1(jl,klevp1)
             shxkj(jl,jk) = (papm1(jl,jk) / paphm1(jl,klevp1))**rgocp 
             th(jl,jk)    =  ptm1(jl,jk)
          END DO
       END DO

       DO jl = 1, klon
          pressg(jl) = paphm1(jl,klevp1)
       END DO

       DO jl = 1, klon
          env(jl) = porsvar(jl)*1.e+4  
       END DO

       DO jl = 1, klon
          zargt1 = zvar(jl,4) - zvar(jl,2)
          zargt2 = zvar(jl,3) - zvar(jl,1)
          zzpos = 0.25*api
          IF (zargt2.EQ.0.) THEN
             IF (zargt1.GT.0.) ztheta =  zzpos
             IF (zargt1.LT.0.) ztheta = -zzpos
             IF (zargt1.EQ.0.) ztheta =  0.
          ELSE
             ztheta = 0.5*ATAN2(zargt1,zargt2)
          ENDIF
          isector = MOD(INT((ztheta/api + 0.625)*4.),4) + 1
          theta(jl) = ztheta
          IF(isector.LE.2) THEN
             isnorm = isector+2
          ELSE
             isnorm = isector-2
          ENDIF
          anisof(jl) = zvar(jl,isector) / (4.e+4+zvar(jl,isnorm)) 
       END DO

       DO jl = 1,klon
          zpr(jl) = zpcons*paprfluxm(jl)
       END DO

       zlat   = ASIN(0.5*twomu(irow))
       coslat = COS(zlat)    

       CALL gwdorexv (pum1,pvm1,th,ptm1,pressg,env,sgj,shj,dsgj,     &
                      shxkj,utendgw,vtendgw,                         &
                      rd,rgocp,zdiagt,klev,lrefp,1,klon, klp2,       &  
                      nazmth,.FALSE.,.TRUE.,.TRUE.,                  &
                      bvfreq,veln,depfac,ub,vb,rmswind,density,      &       
                      vmodb,hitesq,                                  &
                      ampbsq,denfac,visc_mol,alt, anisof, theta,     &
                      zpr,nstep,nresum,coslat)

       !   update tendencies fdue to gwd. the arrays utendgw, vtendgw are the
       !   tendencies (m/s) due to gravity wave drag (orographic + hines).    
       DO jk=ktdia, klev
          DO jl=1, klon
             pvom(jl,jk)  = pvom(jl,jk) + utendgw(jl,jk)
             pvol(jl,jk)  = pvol(jl,jk) + vtendgw(jl,jk)
#ifdef G3X
             g3x01(jl,jk) = g3xm01(jl,jk) + zdiagt*utendgw(jl,jk)
             g3x02(jl,jk) = g3xm02(jl,jk) + zdiagt*vtendgw(jl,jk)
#endif
          END DO
       END DO

    END IF

  END SUBROUTINE cccgwd

  SUBROUTINE gwdorexv (u,v,th,tsg,pressg,env,sgj,shj,dsgj,         &
                       shxkj,utendgw,vtendgw,                      &
                       rgas,rgocp,delt,ilev,lrefp,il1,il2,ilg,     &
                       nazmth,gspong,envelop,extro,                &
                       bvfreq,veln,depfac,ub,vb,rmswind,density,   &
                       vmodb,hitesq,                               &
                       ampbsq,denfac,visc_mol,alt, anisof, theta,  &
                       zpr,nstep,nresum,coslat)

    !
    !     * aug. 14/95 - c. mclandress.
    !     * sep.    95   n. mcfarlane.
    !
    !     * this routine calculates the horizontal wind tendencies
    !     * due to mcfarlane's orographic gw drag scheme, hines'
    !     * doppler spread scheme for "extrowaves" and adds on
    !     * roof drag. it is based on the routine gwdflx8.
    !
    !     * lrefp is the index of the model level below the reference level

    !     * i/o arrays passed from main.
    !     * (pressg = surface pressure)
    !
    !
#ifdef G3X
    USE mo_memory_g3a, ONLY: g3m
    USE mo_memory_g3b, ONLY: g3
    USE mo_control,    ONLY: nrow
#endif
    USE mo_exception,  ONLY: finish

    INTEGER ,INTENT(in) :: ilev
    INTEGER ,INTENT(in) :: il1
    INTEGER ,INTENT(in) :: il2
    INTEGER ,INTENT(in) :: ilg
    INTEGER ,INTENT(in) :: nazmth

    REAL u(ilg,ilev),       v(ilg,ilev),       th(il2,ilev),   &
         tsg(ilg,ilev),        &
         utendgw(il2,ilev), vtendgw(il2,ilev), env(il2),   &
         urow(il2,ilev),    vrow(il2,ilev),    pressg(il2),   &
         uhs(il2,ilev),     vhs(il2,ilev), zpr(il2)

    !     * vertical positioning arrays.

    REAL sgj(il2,ilev), shj(il2,ilev), shxkj(il2,ilev), dsgj(il2,ilev)

    !     * logical switches to control roof drag, envelop gw drag and
    !     * hines' doppler spreading extrowave gw drag.
    !     * lozpr is true for zpr enhancement

    LOGICAL lodrag(il2), lozpr, lorms(il2)
    LOGICAL gspong, envelop, extro

    !     * work arrays.

    REAL m_alpha(il2,ilev,nazmth),     v_alpha(il2,ilev,nazmth),   &
         sigma_alpha(il2,ilev,nazmth), sigsqh_alpha(il2,ilev,nazmth),   &
         drag_u(il2,ilev),   drag_v(il2,ilev),  flux_u(il2,ilev),   &
         flux_v(il2,ilev),   heat(il2,ilev),    diffco(il2,ilev),   &
         bvfreq(il2,ilev),   density(il2,ilev), sigma_t(il2,ilev),   &
         visc_mol(il2,ilev), alt(il2,ilev),     veln(il2,ilev),    &
         depfac(il2,ilev), sigsqmcw(il2,ilev,nazmth),    &
         sigmatm(il2,ilev),    &
         ak_alpha(il2,nazmth),       k_alpha(il2,nazmth),   &
         mmin_alpha(il2,nazmth),     i_alpha(il2,nazmth),   &
         ub(il2),       vb(il2),     vmodb(il2),   hitesq(il2),     & 
         ampbsq(il2),   denfac(il2), deldfac(il2),   &
         rmswind(il2), bvfbot(il2), densbot(il2),   &
         cosaz(il2), sinaz(il2), anisof(il2), phasbs(il2),    &
         ubef(il2), vbef(il2), weight(il2), phase(il2), frnb(il2),   &
         theta(il2), dttdsf(il2)

    !     * thes are the input parameters for hines routine and
    !     * are specified in routine hines_setup. since this is called
    !     * only at first call to this routine these variables must be saved
    !     * for use at subsequent calls. this can be avoided by calling
    !     * hines_setup in main program and passing the parameters as
    !     * subroutine arguements.
    !

    REAL    :: rmscon
    INTEGER :: nmessg, iprint, ilrms

    INTEGER :: naz,icutoff,nsmax,iheatcal
    REAL    :: slope,f1,f2,f3,f5,f6,kstar,alt_cutoff,smco

    INTEGER :: i
    INTEGER :: iaz
    INTEGER :: ierror
    INTEGER :: iil2
    INTEGER :: ilevm
    INTEGER :: ilevp
    INTEGER :: iplev
    INTEGER :: j
    INTEGER :: l
    INTEGER :: len
    INTEGER :: lengw
    INTEGER :: levbot
    INTEGER :: lref
    INTEGER :: lrefp
    INTEGER :: nresum
    INTEGER :: nstep

    REAL    :: apibt
    REAL    :: bvf
    REAL    :: bvfb
    REAL    :: coslat
    REAL    :: cpart
    REAL    :: delfrsq
    REAL    :: delt
    REAL    :: delt2
    REAL    :: denom
    REAL    :: deux
    REAL    :: dfac
    REAL    :: dflxm
    REAL    :: dflxp
    REAL    :: dlnonfac
    REAL    :: dzed
    REAL    :: eta
    REAL    :: fcrit
    REAL    :: grav
    REAL    :: hmin
    REAL    :: hscal
    REAL    :: hsq
    REAL    :: hsqmax
    REAL    :: hsqprop
    REAL    :: pcons
    REAL    :: pcrit
    REAL    :: ratio
    REAL    :: rgas
    REAL    :: rgocp
    REAL    :: sigb
    REAL    :: taufac
    REAL    :: tendfac
    REAL    :: un
    REAL    :: v0
    REAL    :: veltan
    REAL    :: vmin
    REAL    :: wind
    REAL    :: zero

#ifdef G3X
    INTEGER       :: jrow
    REAL, POINTER :: g3x03(:,:), g3x04(:,:), g3xm03(:,:), g3xm04(:,:)
#endif

    !
    !     * constants values defined in data statement are :

    !     * vmin     = miminum wind in the direction of reference level
    !     *            wind before we consider breaking to have occured.
    !     * dmpscal  = damping time for gw drag in seconds.
    !     * taufac   = 1/(length scale).
    !     * hmin     = miminum envelope height required to produce gw drag.
    !     * v0       = value of wind that approximates zero.

    DATA    vmin  /    5.0 /, v0       / 1.e-10 /,   &
            taufac/  5.e-6 /, hmin     /   40000. /,   &
            zero  /    0.0 /, un       /    1.0 /, deux  /    2.0 /,   &
            grav  /9.80616 /, apibt / 1.5708 /,   &
            cpart /    0.7 /, fcrit    / 1. /

    !     * hines extrowave gwd constants defined in data statement are:

    !     * rmscon = root mean square gravity wave wind at lowest level (m/s).
    !     * nmessg  = unit number for printed messages.
    !     * iprint  = 1 to do print out some hines arrays.
    !     * ifl     = first call flag to hines_setup ("save" it)
    !     * pcrit = critical value of zpr (mm/d)
    !     * iplev = level of application of prcit
    !     * pcons = factor of zpr enhancement
    !
    DATA pcrit / 5. /,  pcons / 4.75 /

    DATA    rmscon  / 1. /, iprint   /  0  /, nmessg  /   6   /

    iplev = lrefp-1

    lozpr = .FALSE.

#ifdef G3X
    jrow=nrow(2)
    g3x03  => g3(3)%x(:,1:ilev,jrow)
    g3xm03 => g3m(3)%x(:,1:ilev,jrow)
    g3x04  => g3(4)%x(:,1:ilev,jrow)
    g3xm04 => g3m(4)%x(:,1:ilev,jrow)

    !-----------------------------------------------------------------------
    !
    !      * set to zero g3x's   
    !
    IF (nstep .EQ. nresum) THEN
       g3x03(:,:)  = 0.
       g3xm03(:,:) = 0.
       g3x04(:,:)  = 0.
       g3xm04(:,:) = 0.
    ENDIF
    g3x03(:,:) = 0.
    g3x04(:,:) = 0.
#endif

    !
    !     * set error flag

    ierror = 0

    !     * specify various parameters for hines routine at very first call.
    !     * (note that array k_alpha is specified so make sure that
    !     * it is not overwritten later on).
    !
    CALL hines_setup (naz,slope,f1,f2,f3,f5,f6,kstar,           &
                      icutoff,alt_cutoff,smco,nsmax,iheatcal,   &
                      k_alpha,ierror,nmessg,il2,nazmth,coslat)
    IF (ierror.NE.0) THEN

       !     * if error detected then abort.
       WRITE (nmessg,'(/a)')   ' execution aborted in gwdorexv'
       WRITE (nmessg,'(a,i4)') '     error flag =',ierror
       CALL finish('gwdorexv','Run terminated.')
    END IF
    !
    !     * start gwd calculations.

    delt2 = deux*delt
    ilevp = ilev + 1
    ilevm = ilev - 1
    lref  = lrefp-1
    len = il2 - il1 +1

    !     * vmod is the reference level wind and ( ub, vb)
    !     * are it's unit vector components.

!OCL VCT(MASK)
    DO i=il1,il2
       ub(i)    = MAX(u(i,lref),v0)
       vb(i)    = MAX(v(i,lref),v0)
       vmodb(i) = SQRT(ub(i)**2 + vb(i)**2)
       ub(i)    = ub(i)/vmodb(i)
       vb(i)    = vb(i)/vmodb(i)
       weight(i) = -1.
       cosaz(i) = 0.84339
       sinaz(i) = SQRT(1.-cosaz(i)**2)
       IF(anisof(i).GT.2.) THEN
          cosaz(i) = ub(i)*COS(theta(i))+vb(i)*SIN(theta(i))
          sinaz(i) = ub(i)*SIN(theta(i))-vb(i)*COS(theta(i))
       ENDIF
    END DO

    DO j=1,nazmth
       DO l=1,ilev
          DO i=il1,il2
             sigsqmcw(i,l,j) = 0.
          END DO
       END DO
    END DO

    !     * drag references the points where orographic gw calculations
    !     * will be done, that is (a- if over land, b- if bottom wind  vmin,
    !     * c- if we ask for it and c- if enveloppe height = hmin )
    !     * drag=1. for points that satisfy the above, else drag=0.

    lengw = 0
    !
    DO i=il1,il2
       IF ( vmodb(i).GT.vmin .AND. env(i).GT.hmin) THEN
          lengw   = lengw + 1
          lodrag(i) = .TRUE.
       ELSE
          lodrag(i) = .FALSE.
       ENDIF
    END DO

    !     * initialize necessary arrays.
    !
    DO l=1,ilev
       DO i=il1,il2
          utendgw(i,l) = zero
          vtendgw(i,l) = zero
          urow(i,l)    = u(i,l)
          vrow(i,l)    = v(i,l)
          uhs(i,l) = zero
          vhs(i,l) = zero
          depfac(i,l) = zero
       END DO
    END DO
    !
    !     * if using hines scheme then calculate b v frequency at all points
    !     * and smooth bvfreq.

    DO l=2,ilev
       DO i=il1,il2
          dttdsf(i)=(th(i,l)/shxkj(i,l)-th(i,l-1)/shxkj(i,l-1))/(shj(i,l)-shj(i,l-1))
          dttdsf(i)=MIN(dttdsf(i), -5./sgj(i,l))
       END DO
       DO i=il1,il2
          dttdsf(i)=dttdsf(i)*(sgj(i,l)**rgocp)
       END DO
       DO i=il1,il2
          !          bvfreq(i,l)=SQRT(-dttdsf(i)*sgj(i,l)*(sgj(i,l)**rgocp)/rgas)*grav/tsg(i,l)
          bvfreq(i,l)=SQRT(-dttdsf(i)*sgj(i,l)/rgas)*grav/tsg(i,l)
       END DO
    END DO

    DO l=1,ilev
       DO i=il1,il2
          IF (l.EQ.1) THEN
             bvfreq(i,l) = bvfreq(i,l+1)
          ENDIF
          IF (l.GT.1) THEN
             ratio=5.*LOG(sgj(i,l)/sgj(i,l-1))
             bvfreq(i,l) = (bvfreq(i,l-1) + ratio*bvfreq(i,l))/(1.+ratio)
          ENDIF
       END DO
    END DO
    !
    ! * if no orographic gwd then skip this section (300 continue)
    !
    IF (envelop .AND. lengw.GT.0) THEN
       !
       !
       !     * define square of wind magnitude relative to reference level.
       !
       !
       DO iaz = 1,2
          DO l=1,ilev
             DO i=il1,il2
                IF (lodrag(i)) THEN
                   ubef(i) = ub(i)*cosaz(i)+weight(i)*vb(i)*sinaz(i)
                   vbef(i) = vb(i)*cosaz(i)-weight(i)*ub(i)*sinaz(i)
                   veltan=u(i,l)*ubef(i)+v(i,l)*vbef(i)
                   veln(i,l)=MAX(veltan,v0)
                ENDIF
                bvfreq(i,l) = MAX(bvfreq(i,l), 0.001)
             END DO
          END DO

          !     * calculate effective square launching
          !     * height, reference b v frequency, etc.
          !     * env(i) is the sub-grid scale variance field.

          DO i=il1,il2
             phasbs(i) = 1000.
             phase(i) = 0.
             frnb(i) = 0.
             IF (lodrag(i)) THEN
                sigb=sgj(i,lref)
                bvfb = bvfreq(i,lref)
                hsqmax= cpart*fcrit*(veln(i,lref)/bvfb)**2
                IF(anisof(i).LE. 2.) hsqmax=0.5*hsqmax
                frnb(i) = SQRT(2.*env(i)/hsqmax)
                hsq=MIN(2.*env(i),hsqmax)
                IF(anisof(i).LE. 2.) hsq=MIN(4.*env(i),hsqmax)
                hscal=rgas*tsg(i,lref)/grav
                dfac=bvfb*sigb*veln(i,lref)/hscal
                ampbsq(i)=dfac
                hitesq(i)=hsq
                dzed = hscal*dsgj(i,lref)/sigb
                delfrsq = (4.*env(i)-hsq)/hsq
                dlnonfac=4.0*frnb(i)*(apibt)**2
                deldfac(i) = dlnonfac/(1. +  frnb(i)**2)
             ENDIF
          END DO
          !
          !     * calculate terms from the bottom-up.

          DO l=lref,1,-1
             DO i=il1,il2
                IF (lodrag(i)) THEN
                   wind=veln(i,l)
                   bvf=bvfreq(i,l)
                   hscal=rgas*tsg(i,l)/grav
                   hsqmax=cpart*fcrit*(wind/bvf)**2
                   IF(anisof(i).LE. 2.) hsqmax=0.5*hsqmax
                   dzed=hscal*dsgj(i,l)/sgj(i,l)
                   phase(i)=phase(i)+0.5*dzed*bvf/MAX(wind,un)
                   phase(i)=phase(i)+0.5*dzed*bvf/MAX(wind,un)
                   IF(veln(i,l).LT.un)   hsqmax=zero
                   dfac=bvf*sgj(i,l)*wind/hscal
                   ratio=ampbsq(i)/dfac
                   hsqprop = ratio*hitesq(i)
                   hsq=MIN(hsqprop,hsqmax)
                   phasbs(i) = MIN( phasbs(i), phase(i) )
                   hitesq(i)=hsq
                   ampbsq(i)=dfac
                   depfac(i,l)  =taufac*dfac*hsq
                ENDIF
             END DO
          END DO

          !     * calculate gw tendencies (top and bottom layers first).
          !     * bottom layer keeps initialized value of zero.
          !     * apply gw drag on (urow,vrow) work arrays to be passed to vrtdfs.

          DO i=il1,il2
             IF (lodrag(i)) THEN
                DO l=lrefp,ilev
                   depfac(i,l)=depfac(i,lref)
                END DO
                wind=veln(i,1)
                eta=zero
                dflxp=depfac(i,2)-depfac(i,1)
                IF(dflxp.GT.zero) eta=un
                dfac=6.*delt*depfac(i,1)*eta/wind
                denom=2.*dsgj(i,1)+dfac
                tendfac=dflxp/denom
                denfac(i)=dfac*tendfac
                utendgw(i,1)=utendgw(i,1)-tendfac*ubef(i)
                vtendgw(i,1)=vtendgw(i,1)-tendfac*vbef(i)
             ENDIF
          END DO

          DO l=2,lref
             DO i=il1,il2
                IF (lodrag(i)) THEN
                   wind=veln(i,l)
                   eta=zero
                   dflxp=depfac(i,l+1)-depfac(i,l)
                   dflxm=depfac(i,l)-depfac(i,l-1)
                   IF(dflxp.GT.zero) eta=un
                   dfac=6.*delt*depfac(i,l)*eta/wind
                   denom=2.*dsgj(i,l)+dfac
                   tendfac=(dflxp+dflxm+denfac(i))/denom
                   denfac(i)=dfac*tendfac
                   utendgw(i,l) =utendgw(i,l) -tendfac*ubef(i)
                   vtendgw(i,l) =vtendgw(i,l) -tendfac*vbef(i)
                ENDIF
             END DO
          END DO
          !
          DO i=il1,il2
             IF(anisof(i).LE. 2.)  weight(i)= weight(i) + 2.
          END DO
          !
       END DO
       !
       DO l=1,ilev                                                       
          DO i=il1,il2 
             urow(i,l) = urow(i,l) + delt2*utendgw(i,l)
             vrow(i,l) = vrow(i,l) + delt2*vtendgw(i,l)      
#ifdef G3X
             g3x03(i,l) = g3xm03(i,l)+delt*utendgw(i,l) 
             g3x04(i,l) = g3xm04(i,l)+delt*vtendgw(i,l)       
#endif
          END DO
       END DO
       !
    ENDIF
    !
    ! * end orographic gwd calculation 
    ! 
    !     * calculate gw drag due to hines' extrowaves
    !
    !     * set molecular viscosity to a very small value.
    !     * if the model top is greater than 100 km then the actual
    !     * viscosity coefficient could be specified here.

    DO l=1,ilev
       DO i=il1,il2
          visc_mol(i,l) = 1.5e-5
          drag_u(i,l) = 0.
          drag_v(i,l) = 0.
          flux_u(i,l) = 0.
          flux_v(i,l) = 0.
          heat(i,l) = 0.
          diffco(i,l) = 0.
       END DO
    END DO

    !     * altitude and density at bottom.

    DO i=il1,il2
       hscal = rgas * tsg(i,ilev) / grav
       density(i,ilev) = sgj(i,ilev) * pressg(i) / (grav*hscal)
       alt(i,ilev) = 0.
    END DO

    !     * altitude and density at remaining levels.

    DO l=ilevm,1,-1
       DO i=il1,il2
          hscal = rgas * tsg(i,l) / grav
          alt(i,l) = alt(i,l+1) + hscal * dsgj(i,l) / sgj(i,l)
          density(i,l) = sgj(i,l) * pressg(i) / (grav*hscal)
       END DO
    END DO
    !
    !     * initialize switches for hines gwd calculation
    !
    ilrms = 0

    DO i=il1,il2 
       lorms(i) = .FALSE.
    END DO

    !     * defile bottom launch level
    !
    levbot = iplev
    !
    !     * background wind minus value at bottom launch level.
    !
    DO l=1,levbot
       DO i=il1,il2 
          uhs(i,l) = u(i,l) - u(i,levbot)
          vhs(i,l) = v(i,l) - v(i,levbot)
       END DO
    END DO
    !
    !     * specify root mean square wind at bottom launch level.
    !
    DO i=il1,il2 
       rmswind(i) = rmscon
    END DO

    IF (lozpr) THEN
       DO i=il1,il2 
          IF (zpr(i) .GT. pcrit) THEN
             rmswind(i) = rmscon + ( (zpr(i)-pcrit)/zpr(i) )*pcons
          ENDIF
       END DO
    ENDIF

    DO i=il1,il2 
       IF (rmswind(i) .GT. 0.0) THEN
          ilrms = ilrms+1
          lorms(i) = .TRUE.
       ENDIF
    END DO

    iil2=il2
    !
    !     * calculate gwd (note that diffusion coefficient and
    !     * heating rate only calculated if iheatcal = 1).
    !
    IF (extro .AND. ilrms.GT.0)       THEN                    
       !
       CALL hines_extro0 (drag_u,drag_v,heat,diffco,flux_u,flux_v,    & 
                          uhs,vhs,bvfreq,density,visc_mol,alt,        & 
                          rmswind,k_alpha,m_alpha,v_alpha,            & 
                          sigma_alpha,sigsqh_alpha,ak_alpha,          & 
                          mmin_alpha,i_alpha,sigma_t,densbot,bvfbot,  & 
                          1,iheatcal,icutoff,iprint,nsmax,            & 
                          smco,alt_cutoff,kstar,slope,                & 
                          f1,f2,f3,f5,f6,naz,sigsqmcw,sigmatm,        & 
                          il1,il2,1,levbot,iil2,ilev,nazmth,          & 
                          lorms)

       !     * add on hines' gwd tendencies to orographic tendencies and
       !     * apply hines' gw drag on (urow,vrow) work arrays.

       DO l=1,ilev
          DO i=il1,il2
             utendgw(i,l) = utendgw(i,l) + drag_u(i,l)
             vtendgw(i,l) = vtendgw(i,l) + drag_v(i,l)
          END DO
       END DO
       !
       !   * plot out cut off wavenumber 
       !


       !     * end of hines calculations.
       !
    ENDIF

    !     * apply roof drag.
    !
    !     * finished.

    RETURN
    !-----------------------------------------------------------------------
  END SUBROUTINE gwdorexv

  SUBROUTINE hines_extro0 (drag_u,drag_v,heat,diffco,flux_u,flux_v,  &
                           vel_u,vel_v,bvfreq,density,visc_mol,alt,  &
                           rmswind,k_alpha,m_alpha,v_alpha,          &
                           sigma_alpha,sigsqh_alpha,ak_alpha,        &
                           mmin_alpha,i_alpha,sigma_t,densb,bvfb,    &
                           iorder,iheatcal,icutoff,iprint,nsmax,     &
                           smco,alt_cutoff,kstar,slope,              &
                           f1,f2,f3,f5,f6,naz,sigsqmcw,sigmatm,      &
                           il1,il2,lev1,lev2,nlons,nlevs,nazmth,     &
                           lorms)
    !
    !  main routine for hines' "extrowave" gravity wave parameterization based
    !  on hines' doppler spread theory. this routine calculates zonal
    !  and meridional components of gravity wave drag, heating rates
    !  and diffusion coefficient on a longitude by altitude grid.
    !  no "mythical" lower boundary region calculation is made so it
    !  is assumed that lowest level winds are weak (i.e, approximately zero).
    !
    !  aug. 13/95 - c. mclandress
    !  sept. /95  - n.mcfarlane
    !
    !  modifications:
    !
    !  output arguements:
    !
    !     * drag_u = zonal component of gravity wave drag (m/s^2).
    !     * drag_v = meridional component of gravity wave drag (m/s^2).
    !     * heat   = gravity wave heating (k/sec).
    !     * diffco = diffusion coefficient (m^2/sec)
    !     * flux_u = zonal component of vertical momentum flux (pascals)
    !     * flux_v = meridional component of vertical momentum flux (pascals)
    !
    !  input arguements:
    !
    !     * vel_u      = background zonal wind component (m/s).
    !     * vel_v      = background meridional wind component (m/s).
    !     * bvfreq     = background brunt vassala frequency (radians/sec).
    !     * density    = background density (kg/m^3) 
    !     * visc_mol   = molecular viscosity (m^2/s)
    !     * alt        = altitude of momentum, density, buoyancy levels (m)
    !     *              (note: levels ordered so that alt(i,1) > alt(i,2), etc.)
    !     * rmswind   = root mean square gravity wave wind at lowest level (m/s).
    !     * k_alpha    = horizontal wavenumber of each azimuth (1/m).
    !     * iorder	   = 1 means vertical levels are indexed from top down 
    !     *              (i.e., highest level indexed 1 and lowest level nlevs);
    !     *           .ne. 1 highest level is index nlevs.
    !     * iheatcal   = 1 to calculate heating rates and diffusion coefficient.
    !     * iprint     = 1 to print out various arrays.
    !     * icutoff    = 1 to exponentially damp gwd, heating and diffusion 
    !     *              arrays above alt_cutoff; otherwise arrays not modified.
    !     * alt_cutoff = altitude in meters above which exponential decay applied.
    !     * smco       = smoothing factor used to smooth cutoff vertical 
    !     *              wavenumbers and total rms winds in vertical direction
    !     *              before calculating drag or heating
    !     *              (smco >= 1 ==> 1:smco:1 stencil used).
    !     * nsmax      = number of times smoother applied ( >= 1),
    !     *            = 0 means no smoothing performed.
    !     * kstar      = typical gravity wave horizontal wavenumber (1/m).
    !     * slope      = slope of incident vertical wavenumber spectrum
    !     *              (slope must equal 1., 1.5 or 2.).
    !     * f1 to f6   = hines's fudge factors (f4 not needed since used for
    !     *              vertical flux of vertical momentum).
    !     * naz        = actual number of horizontal azimuths used.
    !     * il1        = first longitudinal index to use (il1 >= 1).
    !     * il2        = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * lev1       = index of first level for drag calculation.
    !     * lev2       = index of last level for drag calculation 
    !     *              (i.e., lev1 < lev2 <= nlevs).
    !     * nlons      = number of longitudes.
    !     * nlevs      = number of vertical levels.
    !     * nazmth     = azimuthal array dimension (nazmth >= naz).
    ! 
    !  work arrays.
    !
    !     * m_alpha      = cutoff vertical wavenumber (1/m).
    !     * v_alpha      = wind component at each azimuth (m/s) and if iheatcal=1
    !     *                holds vertical derivative of cutoff wavenumber.
    !     * sigma_alpha  = total rms wind in each azimuth (m/s).
    !     * sigsqh_alpha = portion of wind variance from waves having wave
    !     *                normals in the alpha azimuth (m/s).
    !     * sigma_t      = total rms horizontal wind (m/s).
    !     * ak_alpha     = spectral amplitude factor at each azimuth 
    !     *                (i.e.,{ajkj}) in m^4/s^2.
    !     * i_alpha      = hines' integral.
    !     * mmin_alpha   = minimum value of cutoff wavenumber.
    !     * densb        = background density at bottom level.
    !     * bvfb         = buoyancy frequency at bottom level and
    !     *                work array for icutoff = 1.
    !
    !     * lorms       = .true. for drag computation 
    !

    USE mo_doctor, ONLY: nout

    INTEGER :: naz, nlons, nlevs, nazmth, il1, il2, lev1, lev2
    INTEGER :: icutoff, nsmax, iorder, iheatcal, iprint
    REAL    :: kstar, f1, f2, f3, f5, f6, slope, alt_cutoff, smco
    REAL    :: drag_u(nlons,nlevs),   drag_v(nlons,nlevs) 
    REAL    :: heat(nlons,nlevs),     diffco(nlons,nlevs)
    REAL    :: flux_u(nlons,nlevs),   flux_v(nlons,nlevs)
    REAL    :: vel_u(nlons,nlevs),    vel_v(nlons,nlevs)
    REAL    :: bvfreq(nlons,nlevs),   density(nlons,nlevs)
    REAL    :: visc_mol(nlons,nlevs), alt(nlons,nlevs)
    REAL    :: rmswind(nlons),       bvfb(nlons),   densb(nlons)
    REAL    :: sigma_t(nlons,nlevs), sigsqmcw(nlons,nlevs,nazmth)
    REAL    :: sigma_alpha(nlons,nlevs,nazmth), sigmatm(nlons,nlevs)
    REAL    :: sigsqh_alpha(nlons,nlevs,nazmth)
    REAL    :: m_alpha(nlons,nlevs,nazmth), v_alpha(nlons,nlevs,nazmth)
    REAL    :: ak_alpha(nlons,nazmth),      k_alpha(nlons,nazmth)
    REAL    :: mmin_alpha(nlons,nazmth) ,   i_alpha(nlons,nazmth)
    REAL    :: smoothr1(nlons,nlevs), smoothr2(nlons,nlevs)

    LOGICAL :: lorms(nlons)
    !
    !  internal variables.
    !
    INTEGER :: levbot, levtop, i, n, l, lev1p, lev2m
    INTEGER :: ilprt1, ilprt2
    !----------------------------------------------------------------------- 
    !
    lev1p = lev1 + 1
    lev2m = lev2 - 1
    !
    !  index of lowest altitude level (bottom of drag calculation).
    !
    levbot = lev2
    levtop = lev1
    IF (iorder.NE.1)  THEN
       WRITE (nout,'(a)') '  error: iorder not one! '
    END IF
    !
    !  buoyancy and density at bottom level.
    !
    DO i = il1,il2
       bvfb(i)  = bvfreq(i,levbot)
       densb(i) = density(i,levbot)
    END DO
    !
    !  initialize some variables
    !
    DO n = 1,naz
       DO l=lev1,lev2
          DO i=il1,il2
             m_alpha(i,l,n) = 0.0
          END DO
       END DO
    END DO
    DO l=lev1,lev2
       DO i=il1,il2
          sigma_t(i,l) = 0.0
       END DO
    END DO
    DO n = 1,naz
       DO i=il1,il2
          i_alpha(i,n) = 0.0
       END DO
    END DO
    !
    !  compute azimuthal wind components from zonal and meridional winds.
    !
    CALL hines_wind ( v_alpha,             & 
                      vel_u, vel_v, naz,   &
                      il1, il2, lev1, lev2, nlons, nlevs, nazmth )
    !
    !  calculate cutoff vertical wavenumber and velocity variances.
    !
    CALL hines_wavnum ( m_alpha, sigma_alpha, sigsqh_alpha, sigma_t,   &
                        ak_alpha, v_alpha, visc_mol, density, densb,   &
                        bvfreq, bvfb, rmswind, i_alpha, mmin_alpha,    &
                        kstar, slope, f1, f2, f3, naz, levbot,         &
                        levtop,il1,il2,nlons,nlevs,nazmth, sigsqmcw,   &
                        sigmatm,lorms)
    !
    !  smooth cutoff wavenumbers and total rms velocity in the vertical 
    !  direction nsmax times, using flux_u as temporary work array.
    !   
    IF (nsmax.GT.0)  THEN
       DO n = 1,naz
          DO l=lev1,lev2
             DO i=il1,il2
                smoothr1(i,l) = m_alpha(i,l,n)
             END DO
          END DO
          CALL vert_smooth (smoothr1,smoothr2, smco, nsmax,  &
                            il1, il2, lev1, lev2, nlons, nlevs )
          DO l=lev1,lev2
             DO i=il1,il2
                m_alpha(i,l,n) = smoothr1(i,l)
             END DO
          END DO
       END DO
       CALL vert_smooth ( sigma_t, smoothr2, smco, nsmax, &
                          il1, il2, lev1, lev2, nlons, nlevs )
    END IF
    !
    !  calculate zonal and meridional components of the
    !  momentum flux and drag.
    !
    CALL hines_flux ( flux_u, flux_v, drag_u, drag_v,             &
                      alt, density, densb, m_alpha,               &
                      ak_alpha, k_alpha, slope, naz,              &
                      il1, il2, lev1, lev2, nlons, nlevs, nazmth, &
                      lorms)
    !
    !  cutoff drag above alt_cutoff, using bvfb as temporary work array.
    !
    IF (icutoff.EQ.1)  THEN	
       CALL hines_exp ( drag_u, bvfb, alt, alt_cutoff, iorder, &
                        il1, il2, lev1, lev2, nlons, nlevs )
       CALL hines_exp ( drag_v, bvfb, alt, alt_cutoff, iorder, &
                        il1, il2, lev1, lev2, nlons, nlevs )
    END IF
    !
    !  print out various arrays for diagnostic purposes.
    !
    IF (iprint.EQ.1)  THEN
       ilprt1 = 15
       ilprt2 = 16
       CALL hines_print ( flux_u, flux_v, drag_u, drag_v, alt,      &
                          sigma_t, sigma_alpha, v_alpha, m_alpha,   &
                          1, 1, 6, ilprt1, ilprt2, lev1, lev2,      &
                          naz, nlons, nlevs, nazmth)
    END IF
    !
    !  if not calculating heating rate and diffusion coefficient then finished.
    !
    IF (iheatcal.NE.1)  RETURN
    !
    !  calculate vertical derivative of cutoff wavenumber (store
    !  in array v_alpha) using centered differences at interior gridpoints
    !  and one-sided differences at first and last levels.
    ! 
    DO n = 1,naz
       DO l = lev1p,lev2m
          DO i = il1,il2
             v_alpha(i,l,n) = ( m_alpha(i,l+1,n) - m_alpha(i,l-1,n) )  &
                            / ( alt(i,l+1) - alt(i,l-1) )
          END DO
       END DO
       DO i = il1,il2
          v_alpha(i,lev1,n) = ( m_alpha(i,lev1p,n) - m_alpha(i,lev1,n) )   &
                            / ( alt(i,lev1p) - alt(i,lev1) )
       END DO
       DO i = il1,il2
          v_alpha(i,lev2,n) = ( m_alpha(i,lev2,n) - m_alpha(i,lev2m,n) )   &
                            / ( alt(i,lev2) - alt(i,lev2m) )
       END DO
    END DO
    !
    !  heating rate and diffusion coefficient.
    !
    CALL hines_heat ( heat, diffco,                                 &
                      m_alpha, v_alpha, ak_alpha, k_alpha,          &
                      bvfreq, density, densb, sigma_t, visc_mol,    &
                      kstar, slope, f2, f3, f5, f6, naz,            &
                      il1, il2, lev1, lev2, nlons, nlevs, nazmth)
    !
    !  finished.
    !
    RETURN
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_extro0

  SUBROUTINE hines_wavnum (m_alpha,sigma_alpha,sigsqh_alpha,sigma_t,   &
                           ak_alpha,v_alpha,visc_mol,density,densb,    &
                           bvfreq,bvfb,rms_wind,i_alpha,mmin_alpha,    &
                           kstar,slope,f1,f2,f3,naz,levbot,levtop,     &
                           il1,il2,nlons,nlevs,nazmth,sigsqmcw,        &
                           sigmatm,lorms)
    !
    !  this routine calculates the cutoff vertical wavenumber and velocity
    !  variances on a longitude by altitude grid for the hines' doppler 
    !  spread gravity wave drag parameterization scheme.
    !  note: (1) only values of four or eight can be used for # azimuths (naz).
    !        (2) only values of 1.0, 1.5 or 2.0 can be used for slope (slope). 
    !
    !  aug. 10/95 - c. mclandress
    !
    !  output arguements:
    !
    !     * m_alpha      = cutoff wavenumber at each azimuth (1/m).
    !     * sigma_alpha  = total rms wind in each azimuth (m/s).
    !     * sigsqh_alpha = portion of wind variance from waves having wave
    !     *                normals in the alpha azimuth (m/s).
    !     * sigma_t      = total rms horizontal wind (m/s).
    !     * ak_alpha     = spectral amplitude factor at each azimuth 
    !     *                (i.e.,{ajkj}) in m^4/s^2.
    !
    !  input arguements:
    !
    !     * v_alpha  = wind component at each azimuth (m/s). 
    !     * visc_mol = molecular viscosity (m^2/s)
    !     * density  = background density (kg/m^3).
    !     * densb    = background density at model bottom (kg/m^3).
    !     * bvfreq   = background brunt vassala frequency (radians/sec).
    !     * bvfb     = background brunt vassala frequency at model bottom.
    !     * rms_wind = root mean square gravity wave wind at lowest level (m/s).
    !     * kstar    = typical gravity wave horizontal wavenumber (1/m).
    !     * slope    = slope of incident vertical wavenumber spectrum
    !     *            (slope = 1., 1.5 or 2.).
    !     * f1,f2,f3 = hines's fudge factors.
    !     * naz      = actual number of horizontal azimuths used (4 or 8).
    !     * levbot   = index of lowest vertical level.
    !     * levtop   = index of highest vertical level 
    !     *            (note: if levtop < levbot then level index 
    !     *             increases from top down).
    !     * il1      = first longitudinal index to use (il1 >= 1).
    !     * il2      = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * nlons    = number of longitudes.
    !     * nlevs    = number of vertical levels.
    !     * nazmth   = azimuthal array dimension (nazmth >= naz).
    !
    !     * lorms       = .true. for drag computation 
    !
    !  input work arrays:
    !
    !     * i_alpha    = hines' integral at a single level.
    !     * mmin_alpha = minimum value of cutoff wavenumber.
    !

    USE mo_doctor, ONLY: nout

    INTEGER :: naz, levbot, levtop, il1, il2, nlons, nlevs, nazmth
    REAL    :: slope, kstar, f1, f2, f3
    REAL    :: m_alpha(nlons,nlevs,nazmth)
    REAL    :: sigma_alpha(nlons,nlevs,nazmth)
    REAL    :: sigalpmc(nlons,nlevs,nazmth)
    REAL    :: sigsqh_alpha(nlons,nlevs,nazmth)
    REAL    :: sigsqmcw(nlons,nlevs,nazmth)
    REAL    :: sigma_t(nlons,nlevs)
    REAL    :: sigmatm(nlons,nlevs)
    REAL    :: ak_alpha(nlons,nazmth)
    REAL    :: v_alpha(nlons,nlevs,nazmth)
    REAL    :: visc_mol(nlons,nlevs)
    REAL    :: f2mod(nlons,nlevs)
    REAL    :: density(nlons,nlevs),  densb(nlons)
    REAL    :: bvfreq(nlons,nlevs),   bvfb(nlons),  rms_wind(nlons)
    REAL    :: i_alpha(nlons,nazmth), mmin_alpha(nlons,nazmth)

    LOGICAL :: lorms(nlons)
    !
    ! internal variables.
    !
    INTEGER :: i, l, n, lstart, lend, lincr, lbelow
    REAL    :: m_sub_m_turb, m_sub_m_mol, m_trial
    REAL    :: visc, visc_min, azfac, sp1, f2mfac
    REAL    :: n_over_m(1000), sigfac(1000)
    DATA  visc_min / 1.e-10 / 
    !-----------------------------------------------------------------------     
    !
    sp1 = slope + 1.
    !
    !  indices of levels to process.
    !
    IF (levbot.GT.levtop)  THEN
       lstart = levbot - 1     
       lend   = levtop         
       lincr  = -1
    ELSE
       WRITE (nout,'(a)') '   error: iorder not one! '
    END IF
    !
    !  use horizontal isotropy to calculate azimuthal variances at bottom level.
    !
    azfac = 1. / float(naz)
    DO n = 1,naz
       DO i = il1,il2
          sigsqh_alpha(i,levbot,n) = azfac * rms_wind(i)**2
       END DO
    END DO
    !
    !  velocity variances at bottom level.
    !
    CALL hines_sigma ( sigma_t, sigma_alpha,          &
                       sigsqh_alpha, naz, levbot,     &
                       il1, il2, nlons, nlevs, nazmth)

    CALL hines_sigma ( sigmatm, sigalpmc,         &
                       sigsqmcw, naz, levbot,     &
                       il1, il2, nlons, nlevs, nazmth)
    !
    !  calculate cutoff wavenumber and spectral amplitude factor 
    !  at bottom level where it is assumed that background winds vanish
    !  and also initialize minimum value of cutoff wavnumber.
    !
    DO n = 1,naz
       DO i = il1,il2
          IF (lorms(i)) THEN
             m_alpha(i,levbot,n) =  bvfb(i) /    & 
                                    ( f1 * sigma_alpha(i,levbot,n)    & 
                                    + f2 * sigma_t(i,levbot) )
             ak_alpha(i,n)   = sigsqh_alpha(i,levbot,n)    & 
                                 / ( m_alpha(i,levbot,n)**sp1 / sp1 )
             mmin_alpha(i,n) = m_alpha(i,levbot,n)
          ENDIF
       END DO
    END DO
    !
    !  calculate quantities from the bottom upwards, 
    !  starting one level above bottom.
    !
    DO l = lstart,lend,lincr
       !
       !  level beneath present level.
       !
       lbelow = l - lincr 
       !
       !  calculate n/m_m where m_m is maximum permissible value of the vertical
       !  wavenumber (i.e., m > m_m are obliterated) and n is buoyancy frequency.
       !  m_m is taken as the smaller of the instability-induced 
       !  wavenumber (m_sub_m_turb) and that imposed by molecular viscosity
       !  (m_sub_m_mol). since variance at this level is not yet known
       !  use value at level below.
       !

!OCL VCT(MASK)
       DO i = il1,il2
          IF (lorms(i)) THEN

             f2mfac=sigmatm(i,lbelow)**2
             f2mod(i,lbelow) =1.+ 2.*f2mfac / ( f2mfac+sigma_t(i,lbelow)**2 )

             visc = amax1 ( visc_mol(i,l), visc_min )
             m_sub_m_turb = bvfreq(i,l) / ( f2 *f2mod(i,lbelow)*sigma_t(i,lbelow))
             m_sub_m_mol = (bvfreq(i,l)*kstar/visc)**0.33333333/f3
             IF (m_sub_m_turb .LT. m_sub_m_mol)  THEN
                n_over_m(i) = f2 *f2mod(i,lbelow)*sigma_t(i,lbelow)
             ELSE
                n_over_m(i) = bvfreq(i,l) / m_sub_m_mol 
             END IF
          ENDIF
       END DO
       !
       !  calculate cutoff wavenumber at this level.
       !
       DO n = 1,naz
          DO i = il1,il2
             IF (lorms(i)) THEN
                !
                !  calculate trial value (since variance at this level is not yet known
                !  use value at level below). if trial value is negative or if it exceeds 
                !  minimum value (not permitted) then set it to minimum value. 
                !                                                                      
                m_trial = bvfb(i) / ( f1 * ( sigma_alpha(i,lbelow,n)+   & 
                          sigalpmc(i,lbelow,n)) + n_over_m(i) + v_alpha(i,l,n) )
                IF (m_trial.LE.0. .OR. m_trial.GT.mmin_alpha(i,n))  THEN
                   m_trial = mmin_alpha(i,n)
                END IF
                m_alpha(i,l,n) = m_trial
                !
                !  reset minimum value of cutoff wavenumber if necessary.
                !
                IF (m_alpha(i,l,n) .LT. mmin_alpha(i,n))  THEN
                   mmin_alpha(i,n) = m_alpha(i,l,n)
                END IF

             ENDIF
          END DO
       END DO
       !
       !  calculate the hines integral at this level.
       !
       CALL hines_intgrl ( i_alpha,                              &
                           v_alpha, m_alpha, bvfb, slope, naz,   &
                           l, il1, il2, nlons, nlevs, nazmth,    &
                           lorms )

       !
       !  calculate the velocity variances at this level.
       !
       DO i = il1,il2
          sigfac(i) = densb(i) / density(i,l) * bvfreq(i,l) / bvfb(i) 
       END DO
       DO n = 1,naz
          DO i = il1,il2
             sigsqh_alpha(i,l,n) = sigfac(i) * ak_alpha(i,n) * i_alpha(i,n)
          END DO
       END DO
       CALL hines_sigma ( sigma_t, sigma_alpha, sigsqh_alpha, naz, l, &
                          il1, il2, nlons, nlevs, nazmth )

       CALL hines_sigma ( sigmatm, sigalpmc, sigsqmcw, naz, l,   &
                          il1, il2, nlons, nlevs, nazmth )
       !
       !  end of level loop.
       !
    END DO
    !
    RETURN
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_wavnum

  SUBROUTINE hines_wind (v_alpha,vel_u,vel_v,  &
                         naz,il1,il2,lev1,lev2,nlons,nlevs,nazmth)
    !
    !  this routine calculates the azimuthal horizontal background wind components 
    !  on a longitude by altitude grid for the case of 4 or 8 azimuths for
    !  the hines' doppler spread gwd parameterization scheme.
    !
    !  aug. 7/95 - c. mclandress
    !
    !  output arguement:
    !
    !     * v_alpha   = background wind component at each azimuth (m/s). 
    !     *             (note: first azimuth is in eastward direction
    !     *              and rotate in counterclockwise direction.)
    !
    !  input arguements:
    !
    !     * vel_u     = background zonal wind component (m/s).
    !     * vel_v     = background meridional wind component (m/s).
    !     * naz       = actual number of horizontal azimuths used (must be 4 or 8).
    !     * il1       = first longitudinal index to use (il1 >= 1).
    !     * il2       = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * lev1      = first altitude level to use (lev1 >=1). 
    !     * lev2      = last altitude level to use (lev1 < lev2 <= nlevs).
    !     * nlons     = number of longitudes.
    !     * nlevs     = number of vertical levels.
    !     * nazmth    = azimuthal array dimension (nazmth >= naz).
    !
    !  constants in data statements.
    !
    !     * cos45 = cosine of 45 degrees. 		
    !     * umin  = minimum allowable value for zonal or meridional 
    !     *         wind component (m/s).
    !
    !  subroutine arguements.
    !
    INTEGER  naz, il1, il2, lev1, lev2
    INTEGER  nlons, nlevs, nazmth
    REAL  v_alpha(nlons,nlevs,nazmth)
    REAL  vel_u(nlons,nlevs), vel_v(nlons,nlevs)
    !
    !  internal variables.
    !
    INTEGER  i, l
    REAL u, v, cos45, umin

    DATA  cos45 / 0.7071068 /
    DATA  umin / 0.001 /
    !-----------------------------------------------------------------------     
    !
    !  case with 4 azimuths.
    !
    IF (naz.EQ.4)  THEN
       DO l = lev1,lev2
          DO i = il1,il2
             u = vel_u(i,l)
             v = vel_v(i,l)
             IF (ABS(u) .LT. umin)  u = umin 
             IF (ABS(v) .LT. umin)  v = umin 
             v_alpha(i,l,1) = u 
             v_alpha(i,l,2) = v
             v_alpha(i,l,3) = - u
             v_alpha(i,l,4) = - v
          END DO
       END DO
    END IF
    !
    !  case with 8 azimuths.
    !
    IF (naz.EQ.8)  THEN
       DO l = lev1,lev2
          DO i = il1,il2
             u = vel_u(i,l)
             v = vel_v(i,l)
             IF (ABS(u) .LT. umin)  u = umin  
             IF (ABS(v) .LT. umin)  v = umin  
             v_alpha(i,l,1) = u 
             v_alpha(i,l,2) = cos45 * ( v + u )
             v_alpha(i,l,3) = v
             v_alpha(i,l,4) = cos45 * ( v - u )
             v_alpha(i,l,5) = - u
             v_alpha(i,l,6) = - v_alpha(i,l,2)
             v_alpha(i,l,7) = - v
             v_alpha(i,l,8) = - v_alpha(i,l,4)
          END DO
       END DO
    END IF
    !
    RETURN
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_wind

  SUBROUTINE hines_flux (flux_u,flux_v,drag_u,drag_v,alt,density,   &
                         densb,m_alpha,ak_alpha,k_alpha,slope,      &
                         naz,il1,il2,lev1,lev2,nlons,nlevs,nazmth,  &
                         lorms)
    !
    !  calculate zonal and meridional components of the vertical flux 
    !  of horizontal momentum and corresponding wave drag (force per unit mass)
    !  on a longitude by altitude grid for the hines' doppler spread 
    !  gwd parameterization scheme.
    !  note: only 4 or 8 azimuths can be used.
    !
    !  aug. 6/95 - c. mclandress
    !
    !  output arguements:
    !
    !     * flux_u = zonal component of vertical momentum flux (pascals)
    !     * flux_v = meridional component of vertical momentum flux (pascals)
    !     * drag_u = zonal component of drag (m/s^2).
    !     * drag_v = meridional component of drag (m/s^2).
    !
    !  input arguements:
    !
    !     * alt       = altitudes (m).
    !     * density   = background density (kg/m^3).
    !     * densb     = background density at bottom level (kg/m^3).
    !     * m_alpha   = cutoff vertical wavenumber (1/m).
    !     * ak_alpha  = spectral amplitude factor (i.e., {ajkj} in m^4/s^2).
    !     * k_alpha   = horizontal wavenumber (1/m).
    !     * slope     = slope of incident vertical wavenumber spectrum.
    !     * naz       = actual number of horizontal azimuths used (must be 4 or 8).
    !     * il1       = first longitudinal index to use (il1 >= 1).
    !     * il2       = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * lev1      = first altitude level to use (lev1 >=1). 
    !     * lev2      = last altitude level to use (lev1 < lev2 <= nlevs).
    !     * nlons     = number of longitudes.
    !     * nlevs     = number of vertical levels.
    !     * nazmth    = azimuthal array dimension (nazmth >= naz).
    !
    !     * lorms       = .true. for drag computation 
    !
    !  constant in data statement.
    !
    !     * cos45 = cosine of 45 degrees. 		
    !
    !  subroutine arguements.
    !
    INTEGER  naz, il1, il2, lev1, lev2, lev2p
    INTEGER  nlons, nlevs, nazmth
    REAL  slope
    REAL  flux_u(nlons,nlevs), flux_v(nlons,nlevs)
    REAL  drag_u(nlons,nlevs), drag_v(nlons,nlevs)
    REAL  alt(nlons,nlevs),    density(nlons,nlevs), densb(nlons)
    REAL  m_alpha(nlons,nlevs,nazmth)
    REAL  ak_alpha(nlons,nazmth), k_alpha(nlons,nazmth)

    LOGICAL lorms(nlons)
    !
    !  internal variables.
    !
    INTEGER  i, l, lev1p, lev2m
    REAL  cos45, prod2, prod4, prod6, prod8, dendz, dendz2
    DATA  cos45 / 0.7071068 /   
    !-----------------------------------------------------------------------
    !
    lev1p = lev1 + 1
    lev2m = lev2 - 1
    lev2p = lev2 + 1
    !
    !  sum over azimuths for case where slope = 1.
    !
    IF (slope.EQ.1.)  THEN
       !
       !  case with 4 azimuths.
       !
       IF (naz.EQ.4)  THEN
          DO l = lev1,lev2
             DO i = il1,il2
                flux_u(i,l) = ak_alpha(i,1)*k_alpha(i,1)*m_alpha(i,l,1)  &
                            - ak_alpha(i,3)*k_alpha(i,3)*m_alpha(i,l,3)
                flux_v(i,l) = ak_alpha(i,2)*k_alpha(i,2)*m_alpha(i,l,2)  &
                            - ak_alpha(i,4)*k_alpha(i,4)*m_alpha(i,l,4)
             END DO
          END DO
       END IF
       !
       !  case with 8 azimuths.
       !
       IF (naz.EQ.8)  THEN
          DO l = lev1,lev2
             DO i = il1,il2
                prod2 = ak_alpha(i,2)*k_alpha(i,2)*m_alpha(i,l,2)
                prod4 = ak_alpha(i,4)*k_alpha(i,4)*m_alpha(i,l,4)
                prod6 = ak_alpha(i,6)*k_alpha(i,6)*m_alpha(i,l,6)
                prod8 = ak_alpha(i,8)*k_alpha(i,8)*m_alpha(i,l,8)
                flux_u(i,l) =  ak_alpha(i,1)*k_alpha(i,1)*m_alpha(i,l,1)   &
                            - ak_alpha(i,5)*k_alpha(i,5)*m_alpha(i,l,5)    &
                            + cos45 * ( prod2 - prod4 - prod6 + prod8 )
                flux_v(i,l) =  ak_alpha(i,3)*k_alpha(i,3)*m_alpha(i,l,3)   &
                            - ak_alpha(i,7)*k_alpha(i,7)*m_alpha(i,l,7)    &
                            + cos45 * ( prod2 + prod4 - prod6 - prod8 )
             END DO
          END DO
       END IF

    END IF
    !
    !  sum over azimuths for case where slope not equal to 1.
    !
    IF (slope.NE.1.)  THEN
       !
       !  case with 4 azimuths.
       !
       IF (naz.EQ.4)  THEN
          DO l = lev1,lev2
             DO i = il1,il2
                flux_u(i,l) = ak_alpha(i,1)*k_alpha(i,1)*m_alpha(i,l,1)**slope  &
                            - ak_alpha(i,3)*k_alpha(i,3)*m_alpha(i,l,3)**slope
                flux_v(i,l) = ak_alpha(i,2)*k_alpha(i,2)*m_alpha(i,l,2)**slope  &
                            - ak_alpha(i,4)*k_alpha(i,4)*m_alpha(i,l,4)**slope
             END DO
          END DO
       END IF
       !
       !  case with 8 azimuths.
       !
       IF (naz.EQ.8)  THEN
          DO l = lev1,lev2
             DO i = il1,il2
                prod2 = ak_alpha(i,2)*k_alpha(i,2)*m_alpha(i,l,2)**slope
                prod4 = ak_alpha(i,4)*k_alpha(i,4)*m_alpha(i,l,4)**slope
                prod6 = ak_alpha(i,6)*k_alpha(i,6)*m_alpha(i,l,6)**slope
                prod8 = ak_alpha(i,8)*k_alpha(i,8)*m_alpha(i,l,8)**slope
                flux_u(i,l) = ak_alpha(i,1)*k_alpha(i,1)*m_alpha(i,l,1)**slope   &
                            - ak_alpha(i,5)*k_alpha(i,5)*m_alpha(i,l,5)**slope   &
                            + cos45 * ( prod2 - prod4 - prod6 + prod8 )
                flux_v(i,l) = ak_alpha(i,3)*k_alpha(i,3)*m_alpha(i,l,3)**slope   &
                            - ak_alpha(i,7)*k_alpha(i,7)*m_alpha(i,l,7)**slope   &
                            + cos45 * ( prod2 + prod4 - prod6 - prod8 )
             END DO
          END DO
       END IF

    END IF
    !
    !  calculate flux from sum.
    !
    DO l = lev1,lev2
       DO i = il1,il2
          flux_u(i,l) = flux_u(i,l) * densb(i) / slope
          flux_v(i,l) = flux_v(i,l) * densb(i) / slope
       END DO
    END DO
    !
    !  calculate drag at intermediate levels using centered differences 
    !      
    DO l = lev1p,lev2m
       DO i = il1,il2
          IF (lorms(i)) THEN
             !ccc       dendz2 = density(i,l) * ( alt(i,l+1) - alt(i,l-1) )
             dendz2 = density(i,l) * ( alt(i,l-1) - alt(i,l) ) 
             !ccc       drag_u(i,l) = - ( flux_u(i,l+1) - flux_u(i,l-1) ) / dendz2
             drag_u(i,l) = - ( flux_u(i,l-1) - flux_u(i,l) ) / dendz2
             !ccc       drag_v(i,l) = - ( flux_v(i,l+1) - flux_v(i,l-1) ) / dendz2
             drag_v(i,l) = - ( flux_v(i,l-1) - flux_v(i,l) ) / dendz2       
          ENDIF
       END DO
    END DO
    !
    !  drag at first and last levels using one-side differences.
    ! 
    DO i = il1,il2
       IF (lorms(i)) THEN
          dendz = density(i,lev1) * ( alt(i,lev1) - alt(i,lev1p) ) 
          drag_u(i,lev1) =  flux_u(i,lev1)  / dendz
          drag_v(i,lev1) =  flux_v(i,lev1)  / dendz
       ENDIF
    END DO
    DO i = il1,il2
       IF (lorms(i)) THEN
          dendz = density(i,lev2) * ( alt(i,lev2m) - alt(i,lev2) )
          drag_u(i,lev2) = - ( flux_u(i,lev2m) - flux_u(i,lev2) ) / dendz
          drag_v(i,lev2) = - ( flux_v(i,lev2m) - flux_v(i,lev2) ) / dendz
       ENDIF
    END DO
    IF (nlevs .GT. lev2) THEN
       DO i = il1,il2
          IF (lorms(i)) THEN
             dendz = density(i,lev2p) * ( alt(i,lev2) - alt(i,lev2p) )
             drag_u(i,lev2p) = -  flux_u(i,lev2)  / dendz
             drag_v(i,lev2p) = - flux_v(i,lev2)  / dendz
          ENDIF
       END DO
    ENDIF

    RETURN
  END SUBROUTINE hines_flux

  SUBROUTINE hines_heat (heat,diffco,m_alpha,dmdz_alpha,             &
                         ak_alpha,k_alpha,bvfreq,density,densb,      &
                         sigma_t,visc_mol,kstar,slope,f2,f3,f5,f6,   &
                         naz,il1,il2,lev1,lev2,nlons,nlevs,nazmth)
    !
    !  this routine calculates the gravity wave induced heating and 
    !  diffusion coefficient on a longitude by altitude grid for  
    !  the hines' doppler spread gravity wave drag parameterization scheme.
    !
    !  aug. 6/95 - c. mclandress
    !
    !  output arguements:
    !
    !     * heat   = gravity wave heating (k/sec).
    !     * diffco = diffusion coefficient (m^2/sec)
    !
    !  input arguements:
    !
    !     * m_alpha     = cutoff vertical wavenumber (1/m).
    !     * dmdz_alpha  = vertical derivative of cutoff wavenumber.
    !     * ak_alpha    = spectral amplitude factor of each azimuth 
    !                     (i.e., {ajkj} in m^4/s^2).
    !     * k_alpha     = horizontal wavenumber of each azimuth (1/m).
    !     * bvfreq      = background brunt vassala frequency (rad/sec).
    !     * density     = background density (kg/m^3).
    !     * densb       = background density at bottom level (kg/m^3).
    !     * sigma_t     = total rms horizontal wind (m/s).
    !     * visc_mol    = molecular viscosity (m^2/s).
    !     * kstar       = typical gravity wave horizontal wavenumber (1/m).
    !     * slope       = slope of incident vertical wavenumber spectrum.
    !     * f2,f3,f5,f6 = hines's fudge factors.
    !     * naz         = actual number of horizontal azimuths used.
    !     * il1         = first longitudinal index to use (il1 >= 1).
    !     * il2         = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * lev1        = first altitude level to use (lev1 >=1). 
    !     * lev2        = last altitude level to use (lev1 < lev2 <= nlevs).
    !     * nlons       = number of longitudes.
    !     * nlevs       = number of vertical levels.
    !     * nazmth      = azimuthal array dimension (nazmth >= naz).
    !
    INTEGER  naz, il1, il2, lev1, lev2, nlons, nlevs, nazmth
    REAL  kstar, slope, f2, f3, f5, f6
    REAL  heat(nlons,nlevs), diffco(nlons,nlevs)
    REAL  m_alpha(nlons,nlevs,nazmth), dmdz_alpha(nlons,nlevs,nazmth)
    REAL  ak_alpha(nlons,nazmth), k_alpha(nlons,nazmth)
    REAL  bvfreq(nlons,nlevs), density(nlons,nlevs),  densb(nlons) 
    REAL  sigma_t(nlons,nlevs), visc_mol(nlons,nlevs)
    !
    ! internal variables.
    !
    INTEGER  i, l, n
    REAL  m_sub_m_turb, m_sub_m_mol, m_sub_m, heatng
    REAL  visc, visc_min, cpgas, zero, sm1
    DATA  zero / 0. /
    !
    ! specific heat at constant pressure
    !
    DATA  cpgas / 1004. / 
    !             
    ! minimum permissible viscosity
    !
    DATA  visc_min / 1.e-10 /       
    !-----------------------------------------------------------------------     
    !
    !  initialize heating array.
    !
    DO l = 1,nlevs
       DO i = 1,nlons
          heat(i,l) = zero
       END DO
    END DO
    !
    !  perform sum over azimuths for case where slope = 1.
    !
    IF (slope.EQ.1.)  THEN
       DO n = 1,naz
          DO l = lev1,lev2
             DO i = il1,il2
                heat(i,l) = heat(i,l) + ak_alpha(i,n) * k_alpha(i,n)   &
                                      * dmdz_alpha(i,l,n) 
             END DO
          END DO
       END DO
    END IF
    !
    !  perform sum over azimuths for case where slope not 1.
    !
    IF (slope.NE.1.)  THEN
       sm1 = slope - 1.
       DO n = 1,naz
          DO l = lev1,lev2
             DO i = il1,il2
                heat(i,l) = heat(i,l) + ak_alpha(i,n) * k_alpha(i,n)   &
                                    * m_alpha(i,l,n)**sm1 * dmdz_alpha(i,l,n) 
             END DO
          END DO
       END DO
    END IF
    !
    !  heating and diffusion.
    !
    DO l = lev1,lev2
       DO i = il1,il2
          !
          !  maximum permissible value of cutoff wavenumber is the smaller 
          !  of the instability-induced wavenumber (m_sub_m_turb) and 
          !  that imposed by molecular viscosity (m_sub_m_mol).
          !
          visc    = amax1 ( visc_mol(i,l), visc_min )
          m_sub_m_turb = bvfreq(i,l) / ( f2 * sigma_t(i,l) )
          m_sub_m_mol  = (bvfreq(i,l)*kstar/visc)**0.33333333/f3
          m_sub_m      = amin1 ( m_sub_m_turb, m_sub_m_mol )

          heatng = - heat(i,l) * f5 * bvfreq(i,l) / m_sub_m   &
                               * densb(i) / density(i,l) 
          diffco(i,l) = f6 * heatng**0.33333333 / m_sub_m**1.33333333
          heat(i,l)   = heatng / cpgas

       END DO
    END DO

    RETURN
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_heat

  SUBROUTINE hines_sigma (sigma_t,sigma_alpha,sigsqh_alpha,  &
                          naz,lev,il1,il2,nlons,nlevs,nazmth)
    !
    !  this routine calculates the total rms and azimuthal rms horizontal 
    !  velocities at a given level on a longitude by altitude grid for 
    !  the hines' doppler spread gwd parameterization scheme.
    !  note: only four or eight azimuths can be used.
    !
    !  aug. 7/95 - c. mclandress
    !
    !  output arguements:
    !
    !     * sigma_t      = total rms horizontal wind (m/s).
    !     * sigma_alpha  = total rms wind in each azimuth (m/s).
    !
    !  input arguements:
    !
    !     * sigsqh_alpha = portion of wind variance from waves having wave
    !     *                normals in the alpha azimuth (m/s).
    !     * naz       = actual number of horizontal azimuths used (must be 4 or 8).
    !     * lev       = altitude level to process.
    !     * il1       = first longitudinal index to use (il1 >= 1).
    !     * il2       = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * nlons     = number of longitudes.
    !     * nlevs     = number of vertical levels.
    !     * nazmth    = azimuthal array dimension (nazmth >= naz).
    !
    !  subroutine arguements.
    !
    INTEGER  lev, naz, il1, il2
    INTEGER  nlons, nlevs, nazmth
    REAL  sigma_t(nlons,nlevs)
    REAL  sigma_alpha(nlons,nlevs,nazmth)
    REAL  sigsqh_alpha(nlons,nlevs,nazmth)
    !
    !  internal variables.
    !
    INTEGER  i, n
    REAL  sum_even, sum_odd 
    !-----------------------------------------------------------------------     
    !
    !  calculate azimuthal rms velocity for the 4 azimuth case.
    !
    IF (naz.EQ.4)  THEN
       DO i = il1,il2
          sigma_alpha(i,lev,1) = SQRT(sigsqh_alpha(i,lev,1)+sigsqh_alpha(i,lev,3))
          sigma_alpha(i,lev,2) = SQRT(sigsqh_alpha(i,lev,2)+sigsqh_alpha(i,lev,4))
          sigma_alpha(i,lev,3) = sigma_alpha(i,lev,1)
          sigma_alpha(i,lev,4) = sigma_alpha(i,lev,2)
       END DO
    END IF
    !
    !  calculate azimuthal rms velocity for the 8 azimuth case.
    !
    IF (naz.EQ.8)  THEN
       DO i = il1,il2
          sum_odd  = ( sigsqh_alpha(i,lev,1) + sigsqh_alpha(i,lev,3)   &
                   + sigsqh_alpha(i,lev,5) + sigsqh_alpha(i,lev,7) ) / 2.
          sum_even = ( sigsqh_alpha(i,lev,2) + sigsqh_alpha(i,lev,4)   &
                   + sigsqh_alpha(i,lev,6) + sigsqh_alpha(i,lev,8) ) / 2.
          sigma_alpha(i,lev,1) = SQRT( sigsqh_alpha(i,lev,1)   &
                               + sigsqh_alpha(i,lev,5) + sum_even )
          sigma_alpha(i,lev,2) = SQRT( sigsqh_alpha(i,lev,2)   &
                               + sigsqh_alpha(i,lev,6) + sum_odd )
          sigma_alpha(i,lev,3) = SQRT( sigsqh_alpha(i,lev,3)   &
                               + sigsqh_alpha(i,lev,7) + sum_even )
          sigma_alpha(i,lev,4) = SQRT( sigsqh_alpha(i,lev,4)   &
                               + sigsqh_alpha(i,lev,8) + sum_odd )
          sigma_alpha(i,lev,5) = sigma_alpha(i,lev,1)
          sigma_alpha(i,lev,6) = sigma_alpha(i,lev,2)
          sigma_alpha(i,lev,7) = sigma_alpha(i,lev,3)
          sigma_alpha(i,lev,8) = sigma_alpha(i,lev,4)
       END DO
    END IF
    !
    !  calculate total rms velocity.
    !
    DO i = il1,il2
       sigma_t(i,lev) = 0.
    END DO
    DO n = 1,naz
       DO i = il1,il2
          sigma_t(i,lev) = sigma_t(i,lev) + sigsqh_alpha(i,lev,n)
       END DO
    END DO
    DO i = il1,il2
       sigma_t(i,lev) = SQRT( sigma_t(i,lev) )
    END DO

    RETURN
    !-----------------------------------------------------------------------     
  END SUBROUTINE hines_sigma

  SUBROUTINE hines_intgrl (i_alpha,v_alpha,m_alpha,bvfb,slope,  &
                           naz,lev,il1,il2,nlons,nlevs,nazmth,lorms)
    !
    !  this routine calculates the vertical wavenumber integral
    !  for a single vertical level at each azimuth on a longitude grid
    !  for the hines' doppler spread gwd parameterization scheme.
    !  note: (1) only spectral slopes of 1, 1.5 or 2 are permitted.
    !        (2) the integral is written in terms of the product qm
    !            which by construction is always less than 1. series
    !            solutions are used for small |qm| and analytical solutions
    !            for remaining values.
    !
    !  aug. 8/95 - c. mclandress
    !
    !  output arguement:
    !
    !     * i_alpha = hines' integral.
    !
    !  input arguements:
    !
    !     * v_alpha = azimuthal wind component (m/s). 
    !     * m_alpha = azimuthal cutoff vertical wavenumber (1/m).
    !     * bvfb    = background brunt vassala frequency at model bottom.
    !     * slope   = slope of initial vertical wavenumber spectrum 
    !     *           (must use slope = 1., 1.5 or 2.)
    !     * naz     = actual number of horizontal azimuths used.
    !     * lev     = altitude level to process.
    !     * il1     = first longitudinal index to use (il1 >= 1).
    !     * il2     = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * nlons   = number of longitudes.
    !     * nlevs   = number of vertical levels.
    !     * nazmth  = azimuthal array dimension (nazmth >= naz).
    !
    !     * lorms       = .true. for drag computation 
    !
    !  constants in data statements:
    !
    !     * qmin = minimum value of q_alpha (avoids indeterminant form of integral)
    !     * qm_min = minimum value of q_alpha * m_alpha (used to avoid numerical
    !     *          problems).
    !

    USE mo_doctor, ONLY: nout

    INTEGER  lev, naz, il1, il2, nlons, nlevs, nazmth
    REAL  i_alpha(nlons,nazmth)
    REAL  v_alpha(nlons,nlevs,nazmth)
    REAL  m_alpha(nlons,nlevs,nazmth)
    REAL  bvfb(nlons), slope

    LOGICAL lorms(nlons)
    LOGICAL lerror(nlons)
    !
    !  internal variables.
    !
    INTEGER  i, n
    REAL  q_alpha, qm, sqrtqm, q_min, qm_min

    DATA  q_min / 1.0 /, qm_min / 0.01 /
    !-----------------------------------------------------------------------     
    !
    !  for integer value slope = 1.
    !
    IF (slope .EQ. 1.)  THEN

       DO n = 1,naz
!OCL VCT(MASK)
          DO i = il1,il2
             IF (lorms(i)) THEN

                q_alpha = v_alpha(i,lev,n) / bvfb(i)
                qm = q_alpha * m_alpha(i,lev,n)
                !
                !  if |qm| is small then use first 4 terms series of taylor series
                !  expansion of integral in order to avoid indeterminate form of integral,
                !  otherwise use analytical form of integral.
                !
                IF (ABS(q_alpha).LT.q_min .OR. ABS(qm).LT.qm_min)  THEN  
                   IF (q_alpha .EQ. 0.)  THEN
                      i_alpha(i,n) = m_alpha(i,lev,n)**2 / 2.
                   ELSE
                      i_alpha(i,n) = (qm**2/2. + qm**3/3. + qm**4/4. + qm**5/5.)/q_alpha**2
                   END IF
                ELSE
                   i_alpha(i,n) = - ( LOG(1.-qm) + qm ) / q_alpha**2
                END IF
             ENDIF
          END DO
       END DO
    END IF
    !
    !  for integer value slope = 2.
    !
    IF (slope .EQ. 2.)  THEN

       DO n = 1,naz
!OCL VCT(MASK)
          DO i = il1,il2
             IF (lorms(i)) THEN
                q_alpha = v_alpha(i,lev,n) / bvfb(i)
                qm = q_alpha * m_alpha(i,lev,n)
                !
                !  if |qm| is small then use first 4 terms series of taylor series
                !  expansion of integral in order to avoid indeterminate form of integral,
                !  otherwise use analytical form of integral.
                !
                IF (ABS(q_alpha).LT.q_min .OR. ABS(qm).LT.qm_min)  THEN  
                   IF (q_alpha .EQ. 0.)  THEN
                      i_alpha(i,n) = m_alpha(i,lev,n)**3 / 3.
                   ELSE
                      i_alpha(i,n) = (qm**3/3. + qm**4/4. + qm**5/5.+ qm**6/6.)/q_alpha**3
                   END IF
                ELSE
                   i_alpha(i,n) = - ( LOG(1.-qm) + qm + qm**2/2.) / q_alpha**3
                END IF
             ENDIF
          END DO
       END DO
    END IF
    !
    !  for real value slope = 1.5
    !
    IF (slope .EQ. 1.5)  THEN
       DO n = 1,naz
!OCL VCT(MASK)
          DO i = il1,il2
             IF (lorms(i)) THEN
                q_alpha = v_alpha(i,lev,n) / bvfb(i)
                qm = q_alpha * m_alpha(i,lev,n)       
                !
                !  if |qm| is small then use first 4 terms series of taylor series
                !  expansion of integral in order to avoid indeterminate form of integral,
                !  otherwise use analytical form of integral.
                !
                IF (ABS(q_alpha).LT.q_min .OR. ABS(qm).LT.qm_min)  THEN  
                   IF (q_alpha .EQ. 0.)  THEN
                      i_alpha(i,n) = m_alpha(i,lev,n)**2.5 / 2.5
                   ELSE
                      i_alpha(i,n) = ( qm/2.5 + qm**2/3.5      &
                                   + qm**3/4.5 + qm**4/5.5 )   &
                                   * m_alpha(i,lev,n)**1.5 / q_alpha
                   END IF
                ELSE
                   qm     = ABS(qm)
                   sqrtqm = SQRT(qm)
                   IF (q_alpha .GE. 0.)  THEN
                      i_alpha(i,n) = ( LOG( (1.+sqrtqm)/(1.-sqrtqm) )  &
                                   - 2.*sqrtqm*(1.+qm/3.) ) / q_alpha**2.5
                   ELSE
                      i_alpha(i,n) = 2. * ( ATAN(sqrtqm) + sqrtqm*(qm/3.-1.) ) &
                                   / ABS(q_alpha)**2.5
                   END IF
                END IF
             ENDIF
          END DO
       END DO
    END IF
    !
    !  if integral is negative (which in principal should not happen) then
    !  print a message and some info since execution will abort when calculating
    !  the variances.
    !
    DO n = 1, naz
      lerror(:) = .FALSE.
      DO i = il1, il2
        IF (i_alpha(i,n) < 0.)  THEN
          lerror(i) = .TRUE.
        END IF
      END DO

      IF (ANY(lerror)) THEN
        WRITE (nout,*) 
        WRITE (nout,*) '******************************'
        WRITE (nout,*) 'hines integral i_alpha < 0 '
        WRITE (nout,*) '  longitude i=',i
        WRITE (nout,*) '  azimuth   n=',n
        WRITE (nout,*) '  level   lev=',lev
        WRITE (nout,*) '  i_alpha =',i_alpha(i,n)
        WRITE (nout,*) '  v_alpha =',v_alpha(i,lev,n)
        WRITE (nout,*) '  m_alpha =',m_alpha(i,lev,n)
        WRITE (nout,*) '  q_alpha =',v_alpha(i,lev,n)/bvfb(i)
        WRITE (nout,*) '  qm      =',v_alpha(i,lev,n)/bvfb(i)*m_alpha(i,lev,n)
        WRITE (nout,*) '******************************'
      END IF
    END DO

    RETURN
  END SUBROUTINE hines_intgrl

  SUBROUTINE hines_setup (naz,slope,f1,f2,f3,f5,f6,kstar,            &
                          icutoff,alt_cutoff,smco,nsmax,iheatcal,    &
                          k_alpha,ierror,nmessg,nlons,nazmth,coslat)
    !
    !  this routine specifies various parameters needed for the
    !  the hines' doppler spread gravity wave drag parameterization scheme.
    !
    !  aug. 8/95 - c. mclandress
    !
    !  output arguements:
    !
    !     * naz        = actual number of horizontal azimuths used
    !     *              (code set up presently for only naz = 4 or 8).
    !     * slope      = slope of incident vertical wavenumber spectrum
    !     *              (code set up presently for slope 1., 1.5 or 2.).
    !     * f1         = "fudge factor" used in calculation of trial value of
    !     *              azimuthal cutoff wavenumber m_alpha (1.2 <= f1 <= 1.9).
    !     * f2         = "fudge factor" used in calculation of maximum
    !     *              permissible instabiliy-induced cutoff wavenumber 
    !     *              m_sub_m_turb (0.1 <= f2 <= 1.4).
    !     * f3         = "fudge factor" used in calculation of maximum 
    !     *              permissible molecular viscosity-induced cutoff wavenumber 
    !     *              m_sub_m_mol (0.1 <= f2 <= 1.4).
    !     * f5         = "fudge factor" used in calculation of heating rate
    !     *              (1 <= f5 <= 3).
    !     * f6         = "fudge factor" used in calculation of turbulent 
    !     *              diffusivity coefficient.
    !     * kstar      = typical gravity wave horizontal wavenumber (1/m)
    !     *              used in calculation of m_sub_m_turb.
    !     * icutoff    = 1 to exponentially damp off gwd, heating and diffusion 
    !     *              arrays above alt_cutoff; otherwise arrays not modified.
    !     * alt_cutoff = altitude in meters above which exponential decay applied.
    !     * smco       = smoother used to smooth cutoff vertical wavenumbers
    !     *              and total rms winds before calculating drag or heating.
    !     *              (==> a 1:smco:1 stencil used; smco >= 1.).
    !     * nsmax      = number of times smoother applied ( >= 1),
    !     *            = 0 means no smoothing performed.
    !     * iheatcal   = 1 to calculate heating rates and diffusion coefficient.
    !     *            = 0 means only drag and flux calculated.
    !     * k_alpha    = horizontal wavenumber of each azimuth (1/m) which
    !     *              is set here to kstar.
    !     * ierror     = error flag.
    !     *            = 0 no errors.
    !     *            = 10 ==> naz > nazmth
    !     *            = 20 ==> invalid number of azimuths (naz must be 4 or 8).
    !     *            = 30 ==> invalid slope (slope must be 1., 1.5 or 2.).
    !     *            = 40 ==> invalid smoother (smco must be >= 1.)
    !
    !  input arguements:
    !
    !     * nmessg  = output unit number where messages to be printed.
    !     * nlons   = number of longitudes.
    !     * nazmth  = azimuthal array dimension (nazmth >= naz).
    !
    INTEGER  naz, nlons, nazmth, iheatcal, icutoff
    INTEGER  nmessg, nsmax, ierror
    REAL  kstar, slope, f1, f2, f3, f5, f6, alt_cutoff, smco, coslat
    REAL  k_alpha(nlons,nazmth)
    REAL  ksmin, ksmax
    !
    ! internal variables.
    !
    INTEGER  i, n
    !-----------------------------------------------------------------------     
    !
    !  specify constants.
    !
    naz   = 8
    slope = 1.
    f1    = 1.5 
    f2    = 0.3 
    f3    = 1.0 
    f5    = 3.0 
    f6    = 1.0       
    ksmin = 1.e-5       
    ksmax = 1.e-4       
    kstar = ksmin/( coslat+(ksmin/ksmax) )      
    icutoff    = 1   
    alt_cutoff = 105.e3
    smco       = 2.0 
    !   smco       = 1.0 
    nsmax      = 5
    !   nsmax      = 2
    iheatcal   = 0 
    !
    !  print information to output file.
    !
    !      WRITE (nmessg,6000)
    ! 6000 FORMAT (/' subroutine hines_setup:')
    !      WRITE (nmessg,*)  '  slope = ', slope
    !      WRITE (nmessg,*)  '  naz = ', naz
    !      WRITE (nmessg,*)  '  f1,f2,f3  = ', f1, f2, f3
    !      WRITE (nmessg,*)  '  f5,f6     = ', f5, f6
    !      WRITE (nmessg,*)  '  kstar     = ', kstar
    !     >           ,'  coslat     = ', coslat
    !      IF (icutoff .EQ. 1)  THEN
    !        WRITE (nmessg,*) '  drag exponentially damped above ',
    !     &                       alt_cutoff/1.e3
    !     END IF
    !      IF (nsmax.LT.1 )  THEN
    !        WRITE (nmessg,*) '  no smoothing of cutoff wavenumbers, etc'
    !      ELSE
    !        WRITE (nmessg,*) '  cutoff wavenumbers and sig_t smoothed:'
    !        WRITE (nmessg,*) '    smco  =', smco
    !        WRITE (nmessg,*) '    nsmax =', nsmax
    !     END IF
    !
    !  check that things are setup correctly and log error if not
    !
    ierror = 0
    IF (naz .GT. nazmth)                                  ierror = 10
    IF (naz.NE.4 .AND. naz.NE.8)                          ierror = 20
    IF (slope.NE.1. .AND. slope.NE.1.5 .AND. slope.NE.2.) ierror = 30
    IF (smco .LT. 1.)                                     ierror = 40
    !
    !  use single value for azimuthal-dependent horizontal wavenumber.
    !
    DO n = 1,naz
       DO i = 1,nlons
          k_alpha(i,n) = kstar
       END DO
    END DO

    RETURN
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_setup

  SUBROUTINE hines_print (flux_u,flux_v,drag_u,drag_v,alt,sigma_t,    &
                          sigma_alpha,v_alpha,m_alpha,                &
                          iu_print,iv_print,nmessg,                   &
                          ilprt1,ilprt2,levprt1,levprt2,              &
                          naz,nlons,nlevs,nazmth)
    !
    !  print out altitude profiles of various quantities from
    !  hines' doppler spread gravity wave drag parameterization scheme.
    !  (note: only for naz = 4 or 8). 
    !
    !  aug. 8/95 - c. mclandress
    !
    !  input arguements:
    !
    !     * iu_print = 1 to print out values in east-west direction.
    !     * iv_print = 1 to print out values in north-south direction.
    !     * nmessg   = unit number for printed output.
    !     * ilprt1   = first longitudinal index to print.
    !     * ilprt2   = last longitudinal index to print.
    !     * levprt1  = first altitude level to print.
    !     * levprt2  = last altitude level to print.
    !
    INTEGER  naz, ilprt1, ilprt2, levprt1, levprt2
    INTEGER  nlons, nlevs, nazmth
    INTEGER  iu_print, iv_print, nmessg
    REAL  flux_u(nlons,nlevs), flux_v(nlons,nlevs)
    REAL  drag_u(nlons,nlevs), drag_v(nlons,nlevs)
    REAL  alt(nlons,nlevs), sigma_t(nlons,nlevs)
    REAL  sigma_alpha(nlons,nlevs,nazmth)
    REAL  v_alpha(nlons,nlevs,nazmth), m_alpha(nlons,nlevs,nazmth)
    !
    !  internal variables.
    !
    INTEGER  n_east, n_west, n_north, n_south
    INTEGER  i, l
    !-----------------------------------------------------------------------
    !
    !  azimuthal indices of cardinal directions.
    !
    n_east = 1
    IF (naz.EQ.4)  THEN
       n_west  = 3       
       n_north = 2
       n_south = 4       
    ELSE IF (naz.EQ.8)  THEN
       n_west  = 5       
       n_north = 3
       n_south = 7       
    END IF
    !
    !  print out values for range of longitudes.
    !
    DO i = ilprt1,ilprt2
       !
       !  print east-west wind, sigmas, cutoff wavenumbers, flux and drag.
       !
       IF (iu_print.EQ.1)  THEN
          WRITE (nmessg,*) 
          WRITE (nmessg,'(a,i3)') 'hines gw (east-west) at longitude i =',i
          WRITE (nmessg,6005) 
6005      FORMAT (15x,' u ',2x,'sig_e',2x,'sig_t',3x,'m_e',  &
                            4x,'m_w',4x,'fluxu',5x,'gwdu')
          DO l = levprt1,levprt2
             WRITE (nmessg,6701) alt(i,l)/1.e3, v_alpha(i,l,n_east),   &
                                 sigma_alpha(i,l,n_east), sigma_t(i,l),  &
                                 m_alpha(i,l,n_east)*1.e3,   &
                                 m_alpha(i,l,n_west)*1.e3,  &
                                 flux_u(i,l)*1.e5, drag_u(i,l)*24.*3600.
          END DO
6701      FORMAT (' z=',f7.2,1x,3f7.1,2f7.3,f9.4,f9.3)
       END IF
       !
       !  print north-south winds, sigmas, cutoff wavenumbers, flux and drag.
       !
       IF (iv_print.EQ.1)  THEN
          WRITE(nmessg,*) 
          WRITE(nmessg,'(a,i3)') 'hines gw (north-south) at longitude i =',i
          WRITE(nmessg,6006) 
6006      FORMAT (15x,' v ',2x,'sig_n',2x,'sig_t',3x,'m_n',   &
                            4x,'m_s',4x,'fluxv',5x,'gwdv')
          DO l = levprt1,levprt2
             WRITE (nmessg,6701) alt(i,l)/1.e3, v_alpha(i,l,n_north),    &
                                 sigma_alpha(i,l,n_north), sigma_t(i,l),   &
                                 m_alpha(i,l,n_north)*1.e3,    &
                                 m_alpha(i,l,n_south)*1.e3,   &
                                 flux_v(i,l)*1.e5, drag_v(i,l)*24.*3600.
          END DO
       END IF
       !
    END DO
    !
    RETURN
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_print

  SUBROUTINE hines_exp (darr,data_zmax,alt,alt_exp,iorder,  &
                        il1,il2,lev1,lev2,nlons,nlevs)
    !
    !  this routine exponentially damps a longitude by altitude array 
    !  of darr above a specified altitude.
    !
    !  aug. 13/95 - c. mclandress
    !
    !  output arguements:
    !
    !     * darr = modified data array.
    !
    !  input arguements:
    !
    !     * darr    = original data array.
    !     * alt     = altitudes.
    !     * alt_exp = altitude above which exponential decay applied.
    !     * iorder	= 1 means vertical levels are indexed from top down 
    !     *           (i.e., highest level indexed 1 and lowest level nlevs);
    !     *           .ne. 1 highest level is index nlevs.
    !     * il1     = first longitudinal index to use (il1 >= 1).
    !     * il2     = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * lev1    = first altitude level to use (lev1 >=1). 
    !     * lev2    = last altitude level to use (lev1 < lev2 <= nlevs).
    !     * nlons   = number of longitudes.
    !     * nlevs   = number of vertical
    !
    !  input work arrays:
    !
    !     * data_zmax = data values just above altitude alt_exp.
    !
    INTEGER  iorder, il1, il2, lev1, lev2, nlons, nlevs
    REAL  alt_exp
    REAL  darr(nlons,nlevs), data_zmax(nlons), alt(nlons,nlevs)
    !
    ! internal variables.
    !
    INTEGER  levbot, levtop, lincr, i, l
    REAL  hscale
    DATA  hscale / 5.e3 /
    !-----------------------------------------------------------------------     
    !
    !  index of lowest altitude level (bottom of drag calculation).
    !
    levbot = lev2
    levtop = lev1
    lincr  = 1
    IF (iorder.NE.1)  THEN
       levbot = lev1
       levtop = lev2
       lincr  = -1
    END IF
    !
    !  data values at first level above alt_exp.
    !
    DO i = il1,il2
       DO l = levtop,levbot,lincr
          IF (alt(i,l) .GE. alt_exp)  THEN
             data_zmax(i) = darr(i,l) 
          END IF
       END DO
    END DO
    !
    !  exponentially damp field above alt_exp to model top at l=1.
    !
    DO l = 1,lev2 
       DO i = il1,il2
          IF (alt(i,l) .GE. alt_exp)  THEN
             darr(i,l) = data_zmax(i) * EXP( (alt_exp-alt(i,l))/hscale )
          END IF
       END DO
    END DO
    !
    RETURN
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_exp

  SUBROUTINE vert_smooth (darr,work,coeff,nsmooth,  &
                          il1,il2,lev1,lev2,nlons,nlevs)
    !
    !  smooth a longitude by altitude array in the vertical over a
    !  specified number of levels using a three point smoother. 
    !
    !  note: input array darr is modified on output!
    !
    !  aug. 3/95 - c. mclandress
    !
    !  output arguement:
    !
    !     * darr    = smoothed array (on output).
    !
    !  input arguements:
    !
    !     * darr    = unsmoothed array of data (on input).
    !     * work    = work array of same dimension as darr.
    !     * coeff   = smoothing coefficient for a 1:coeff:1 stencil.
    !     *           (e.g., coeff = 2 will result in a smoother which
    !     *           weights the level l gridpoint by two and the two 
    !     *           adjecent levels (l+1 and l-1) by one).
    !     * nsmooth = number of times to smooth in vertical.
    !     *           (e.g., nsmooth=1 means smoothed only once, 
    !     *           nsmooth=2 means smoothing repeated twice, etc.)
    !     * il1     = first longitudinal index to use (il1 >= 1).
    !     * il2     = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * lev1    = first altitude level to use (lev1 >=1). 
    !     * lev2    = last altitude level to use (lev1 < lev2 <= nlevs).
    !     * nlons   = number of longitudes.
    !     * nlevs   = number of vertical levels.
    !
    !  subroutine arguements.
    !
    INTEGER  nsmooth, il1, il2, lev1, lev2, nlons, nlevs
    REAL  coeff
    REAL  darr(nlons,nlevs), work(nlons,nlevs)
    !
    !  internal variables.
    !
    INTEGER  i, l, ns, lev1p, lev2m
    REAL  sum_wts
    !-----------------------------------------------------------------------     
    !
    !  calculate sum of weights.
    !
    sum_wts = coeff + 2.
    !
    lev1p = lev1 + 1
    lev2m = lev2 - 1
    !
    !  smooth nsmooth times
    !
    DO ns = 1,nsmooth
       !
       !  copy darr into work array.
       !
       DO l = lev1,lev2
          DO i = il1,il2
             work(i,l) = darr(i,l)
          END DO
       END DO
       !
       !  smooth array work in vertical direction and put into darr.
       !
       DO l = lev1p,lev2m
          DO i = il1,il2
             darr(i,l) = (work(i,l+1)+coeff*work(i,l)+work(i,l-1) ) / sum_wts 
          END DO
       END DO
    END DO

    RETURN
    !-----------------------------------------------------------------------
  END SUBROUTINE vert_smooth

END MODULE mo_midatm
