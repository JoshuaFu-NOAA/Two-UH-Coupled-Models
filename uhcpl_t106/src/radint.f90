!+ organises the radiation full computations.

!OCL NOALIAS

!#define DEBUG

SUBROUTINE radint

  ! Description:
  !
  ! Organises the radiation full computations.
  !
  ! Method:
  !
  ! This routine organises the input/output for the black-box
  ! radiation computations performed in *radlsw* every time there is a
  ! full radiation time step. Input are prognostic model variables at
  ! time step t-1, surface values of short-wave albedo and long-wave
  ! emissivity and climatological values for aerosols and ozone (time
  ! of the year dependent). Output are flux transmissivities and
  ! emissivities at all the half levels of the grid (respectively
  ! ratio solar flux/solar input and ratio thermal flux/local
  ! black-body flux). This output will be used in *radheat* at all
  ! time steps until the next full radiation time step.
  !
  ! A call to subroutine *solang* gives fields of solar zenith
  ! angles and relative day length (results depending on the switch on
  ! or off) of the diurnal cycle. The consistency of these values with
  ! the use they will later have in *radheat* was ensured by giving to
  ! *solang* an input corresponding to the middle of the period during
  ! which the results of *radint* will be valid.
  ! *Legendre and *Fourier transforms allow to go from the t5 and
  ! t10 spectral definitions of ozone and aerosols' climatologies to
  ! the models' grid. The computation of the *Legendre polynomials is
  ! done in subroutine *legtri* that gets as input the sine of
  ! latitude.
  ! Temperatures and relative humidities are inter-/extrapolated
  ! from the prognostic levels to the levels' boundaries with a method
  ! (non linear in p) consistent with the one used in *radheat* for
  ! the same purpose.
  ! The actual handling of the input/output for *radlsw* is rather
  ! straightforward if one knows the choices made for the vertical
  ! distributions of the radiative agents: clouds, aerosols, gases
  ! and temperatures.
  !
  ! *radint* is called from *physc*.
  ! The communications with *radheat* have been described above
  ! and those with *radlsw* are via dummy list.
  !
  ! Externals:
  ! *solang*, *legtri* and *radlsw*
  !
  ! Reference:
  ! See radiation's part of the model's documentation for details
  ! about the mathematics of this routine.
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, July 1993, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! A. Rhodin, MPI, December 1998, lookup tables removed 
  ! U. Niemeier, MPI, July 1999, chance nlon, nlp2 nrow etc for parallelization
  !
  ! for more details see file AUTHORS
  !

  USE mo_memory_gl,      ONLY: x
  USE mo_memory_g1a,     ONLY: qm1, tm1
  USE mo_memory_g3a,     ONLY: aclcm, aclcvm, albedom, albm, emtefm, emterm, &
                               forestm, snm1m, trsofm, trsolm, tsm1m
  USE mo_memory_g3b,     ONLY: seaice
  USE mo_tmp_buffer,     ONLY: aphm1, apm1, dia1, diag, loglac, loland
  USE mo_control,        ONLY: dtime, lmidatm, maxrow, ncbase, nlev, nlevp1, &
                               nrow, ntbase, ntimst, lamip2
  USE mo_doctor,         ONLY: nout
  USE mo_gaussgrid,      ONLY: budw, coslon, sinlon, sqcst, twomu
  USE mo_rad_switches,   ONLY: laer, lcfc, lgadsrh, lrad, lsolc, nmonth, &
                               nradfr, nradpfr, nradpla
  USE mo_physc1,         ONLY: cdissem, crae, czen1m, czen2m, czen3m
  USE mo_constants,      ONLY: api, cpd, dayl, g, solc, stbo, tmelt, vtmpc1, yearl
  USE mo_start_dataset,  ONLY: nstep
  USE mo_rad1,           ONLY: caedc, caeds, caelc, caels, caesc, caess, &
                               caeuc, caeus, cozhc, cozhs, cozqc, cozqs
  USE mo_rad2,           ONLY: caeopd, caeopl, caeops, caeopu, cstbga, &
                               ctrbga, ctrpt, cvdaed, cvdael, cvdaes, &
                               cvdaeu, cvobga
  USE mo_radint,         ONLY: cemiss, calbsea, calbsno, zepalb, zepclc, &
                               zeph2o, zepsec, zsnowal
  USE mo_aerosols,       ONLY: newaer
#ifndef NOLOOKUP
  USE mo_convect_tables, ONLY: tlucua, jptlucu1, jptlucu2, &
                               lookuperror, lookupoverflow
#else
  USE mo_constants,      ONLY: c2es, c3ies, c3les, c4ies, c4les
#endif
  USE mo_midatm,         ONLY: noz, scloz, scloza, spe, ozone
  USE mo_year,           ONLY: cd2dat, cd2dy
!  USE m_solang,          ONLY: solang ! module procedure
  USE mo_decomposition,  ONLY: dc => local_decomposition
  ! Added for amip2
  USE mo_timeint,        ONLY: wgt1, wgt2, nmw1, nmw2
  USE mo_physc2
  USE mo_ozone,          ONLY: jpozlev


  IMPLICIT NONE

  !  Local array bounds
  INTEGER :: nglon, nglpx, nglat, nlof

  !  Local scalars:
  REAL :: caeros, zaeqdn, zaeqdo, zaeqln, zaeqlo, zaeqsn, zaeqso, zaequn,         &
          zaequo, zaetr, zalbmax, zalbmin, zalbmn0, zalbmn1, zalbmx0, zalbmx1,    &
          zalbpr, zalbsn, zalte, zbud, zcos1, zcos10, zcos2, zcos3, zcos4, zcos5, &
          zcos6, zcos7, zcos8, zcos9, zcphn3, zcpho3, zcrae, zdalb, zdegday,      &
          zdials, zdift1, zdift2, zdift3, zdift4, zdift5, zlatd, zqofn, zqofo,    &
          zrh, zsct, zsdpn3, zsdpo3, zsin, zsin1, zsin10, zsin2, zsin3, zsin4,    &
          zsin5, zsin6, zsin7, zsin8, zsin9, ztalb, ztim1, ztim2, ztim3,          &
          ztsnmelt, zzdi, zzozh

  ! Added for AMIP2 
  REAL :: zconvoz, ztsalb

  INTEGER :: imm, imnc, imns, irow, jrow, jaer, jk, jkk, jkl, jl, jmm, &
             jn, jnn, kaer, kcfc, kewaer, kmode, nobase
#ifndef NOLOOKUP
  INTEGER :: it
#else
  LOGICAL :: lo
#endif
  LOGICAL :: lodia, lodiap

  !  Local arrays: 
  REAL :: zaed(dc%nglon)  , zael(dc%nglon),      zaes(dc%nglon),      zaetrn(dc%nglon),   &
          zaetro(dc%nglon), zaeu(dc%nglon),      zalp(66),             zalso(dc%nglon),   &
          zamu0(dc%nglpx),  zclc(dc%nglon,nlev), zclwa(dc%nglon,nlev),                    &
          zdia(dc%nglon),   zdia1(7),             zdiaf(dc%nglon,7),   zdiag(nlevp1,10),  &
          zdiat(nlevp1),     zdp(dc%nglon,nlev),  zfaed(21),            zfael(21),        &
          zfaes(21),         zfaeu(21),            zfls(dc%nglon,nlevp1),                 &
          zflsc(dc%nglon,2),zflt(dc%nglon,nlevp1), zfltc(dc%nglon,2), zfozh(11),          &
          zfozq(11),         zhti(dc%nglon,nlevp1), zmu0(dc%nglon),    zozh(dc%nglon),    &
          zozq(dc%nglon),   zq(dc%nglon,nlev),   zqof(dc%nglon,nlev), zqs(dc%nglon,nlev), &
          zrdayl(dc%nglpx), zsaer(dc%nglon,nlev,5+newaer)
  INTEGER :: iaerh(dc%nglon,nlev)


  !  Added for AMIP2 
  REAL ::  zdpe(jpozlev)

  ! lmidatm
  REAL :: zmomid, zmonthl, zyt, zmosec, zw1, zw2, zw3, zkp1, zk 
  INTEGER :: imonth, iday, id, im, iy, im1, im2, im3, k, krow

  !  External subroutines 
  EXTERNAL legtri, radlsw

  !  Intrinsic functions 
  INTRINSIC ASIN, ABS, MAX, MIN, MOD, SQRT, SUM


  !  Executable statements 

#ifndef NOLOOKUP
  lookupoverflow = .FALSE.
#endif

  ! Data statements
  ! *calbsea* and *calbsno* are albedo values for
  ! open sea and thick snow respectively. *zsnowal*
  ! is a value of snow depth (in equivalent water) for which the
  ! snow starts beeing considered as thick. *zcardi* is the specific
  ! atmospheric content in co2.

  DATA zalbmn0, zalbmx0, zalbmn1, zalbmx1/.4, .8, .3, .4/

  ! Pressure thicknesses of the 2d Fortuin ozone distribution

  DATA zdpe / &
          40.,    35.,    75.,   100.,   150.,    200.,   250., &
         650.,  1000.,  1500.,  2000.,  2500.,   4000.,  5000., &
        7500., 15000., 20000., 30000., 10000. /
  ! ------------------------------------------------------------------

  ! Security parameter.

  ! *zepsec* avoids 0/0 in the diagnostic of total cloud cover.
  ! *zepclc* is a security to avoid zero or one cloud covers and
  ! *zeph2o* is a security to avoid water vapour content in a layer
  ! *         to be more then the respective value at saturation.
  ! *zepalb* is a security to avoid zero albedos.

  ! Computational constants.

  ztalb = tmelt - 10.
  ztsalb = ctfreez - 10.  ! Added for amip2

  zsct = cdissem*solc
  zdegday = g*dayl/cpd
  zcrae = crae*(crae+2.)
  ! Number of additional aerosols
  kewaer = newaer

  !  Local array bounds
  irow = nrow(1)        ! global ping pong latitude index
  jrow = nrow(2)        ! local continuous latitude index
  krow = dc%glat(jrow)  ! global continuous latitude index
  nglon = dc%nglon      ! local number of longitudes
  nglpx = dc%nglpx      ! local number of longitudes allocated
  nglat = dc%nglat      ! local number of latitudes
  nlof  = dc%glon(jrow) ! longitude offset to global field

  ! Conversion factor: volume mixing ratio -> mass mixing ratio of ozone
  IF (lamip2) zconvoz=((3.*16.)/28.84)*1.E-6

  zbud = budw(irow)
  lodiap = .FALSE.

  IF (lrad) THEN

!-- 1. Preparation of the radiation diagnostics

    lodia = .FALSE.
    IF (nradpfr/=0) THEN
      lodia = MOD(nstep,nradpfr) == 0
    END IF
    IF (nradpla/=0 .AND. lodia) THEN
      lodiap = MOD(ABS(maxrow+1-2*((irow+1)/2)),nradpla) == 0
    END IF
    IF (lodia) THEN
      DO jn = 1, 8
        DO jk = 1, nlevp1
          zdiag(jk,jn) = 0.
        END DO
      END DO
      DO jn = 1, 7
        zdia1(jn) = 0.
      END DO
      zdials = 1./nglon
    END IF

!-- 2. Solar angle and ozone/aerosol parameters computations.

!-- 2.1 Introduce the latitude dependency

    zsin = .5*twomu(irow)

    ! Use mean value between two full radiation time-steps for orbital
    ! parameters.

    ztim1 = czen1m*zsin
    ztim2 = -czen2m*sqcst(irow)
    ztim3 = czen3m*sqcst(irow)

!-- 2.2 Call to *solang* for zenith angle and daylength

    CALL solang(nglon,nglpx,nlof,ztim1,ztim2,ztim3,zamu0,zrdayl)

    DO jl = 1, nglon
      zmu0(jl) = crae/(SQRT(zamu0(jl)**2+zcrae)-zamu0(jl))
    END DO

    ! No legendre transformation for ozone for the middle atmosphere
    ! and AMIP2 version
    IF (.NOT.lmidatm .AND. .NOT.lamip2) THEN    

!-- 2.3 Call to legtri

      CALL legtri(zsin,6,zalp)

!-- 2.4 *Legendre transform for ozone

      DO jmm = 1, 11
        zfozq(jmm) = 0.
        zfozh(jmm) = 0.
      END DO
      imm = 0
      imnc = 0
      imns = 0
      DO jmm = 1, 6
        imm = imm + 1
        DO jnn = jmm, 6
          imnc = imnc + 1
          zfozq(imm) = zfozq(imm) + zalp(imnc)*cozqc(imnc)
          zfozh(imm) = zfozh(imm) + zalp(imnc)*cozhc(imnc)
        END DO
        IF (jmm/=1) THEN
          imm = imm + 1
          DO jnn = jmm, 6
            imns = imns + 1
            zfozq(imm) = zfozq(imm) + zalp(imns+6)*cozqs(imns)
            zfozh(imm) = zfozh(imm) + zalp(imns+6)*cozhs(imns)
          END DO
        END IF
      END DO
    END IF

    ! 2.5 Call to legtri

    CALL legtri(zsin,11,zalp)

    ! 2.6 *Legendre transform for aerosols

    caeros = 1.E-15
    DO jmm = 1, 21
      zfaes(jmm) = 0.
      zfael(jmm) = 0.
      zfaeu(jmm) = 0.
      zfaed(jmm) = 0.
    END DO
    imm = 0
    imnc = 0
    imns = 0
    DO jmm = 1, 11
      imm = imm + 1
      DO jnn = jmm, 11
        imnc = imnc + 1
        zfaes(imm) = zfaes(imm) + zalp(imnc)*caesc(imnc)
        zfael(imm) = zfael(imm) + zalp(imnc)*caelc(imnc)
        zfaeu(imm) = zfaeu(imm) + zalp(imnc)*caeuc(imnc)
        zfaed(imm) = zfaed(imm) + zalp(imnc)*caedc(imnc)
      END DO
      IF (jmm/=1) THEN
        imm = imm + 1
        DO jnn = jmm, 11
          imns = imns + 1
          zfaes(imm) = zfaes(imm) + zalp(imns+11)*caess(imns)
          zfael(imm) = zfael(imm) + zalp(imns+11)*caels(imns)
          zfaeu(imm) = zfaeu(imm) + zalp(imns+11)*caeus(imns)
          zfaed(imm) = zfaed(imm) + zalp(imns+11)*caeds(imns)
        END DO
      END IF
    END DO

!-- 3. Prepare input for radiation

!-- 3.1 Humidity

    DO jk = 1, nlev
      DO jl = 1, nglon
#ifndef NOLOOKUP
        it = INT(tm1(jl,jk,jrow)*1000.)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
#ifdef DEBUG
        IF (it<jptlucu1 .OR. it>jptlucu2) THEN
          WRITE(*,*) 'radint : 1 it=',it,jl,jk
        END IF
#endif
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqs(jl,jk) = tlucua(it)/apm1(jl,jk)
#else
        lo = tm1(jl,jk,jrow) > tmelt
        zqs(jl,jk) = c2es*EXP(MERGE(c3les,c3ies,lo)*(tm1(jl,jk,jrow)-tmelt) &
                / (tm1(jl,jk,jrow)-MERGE(c4les,c4ies,lo)))/apm1(jl,jk)
#endif
        zqs(jl,jk) = MIN(zqs(jl,jk),0.5)
        zqs(jl,jk) = zqs(jl,jk)/(1.-vtmpc1*zqs(jl,jk))
        zqs(jl,jk) = MAX(2.*zeph2o,zqs(jl,jk))
        zq(jl,jk)  = MAX(qm1(jl,jk,jrow),zeph2o)
      END DO
    END DO
#ifndef NOLOOKUP
    IF (lookupoverflow) CALL lookuperror ('radint      ')
#endif

!-- 3.2 Define relative humidity classes for *gads* aerosol
!       optical parameters.

    IF (newaer>0) THEN
      IF (lgadsrh) THEN
        ! Variable r.h. classes:
        DO jk = 1, nlev
          DO jl = 1, nglon
            zrh = zq(jl,jk)/zqs(jl,jk)
            IF (zrh<0.25) THEN
              iaerh(jl,jk) = 1
            ELSE IF (zrh<0.60) THEN
              iaerh(jl,jk) = 2
            ELSE IF (zrh<0.75) THEN
              iaerh(jl,jk) = 3
            ELSE IF (zrh<0.85) THEN
              iaerh(jl,jk) = 4
            ELSE IF (zrh<0.925) THEN
              iaerh(jl,jk) = 5
            ELSE IF (zrh<0.965) THEN
              iaerh(jl,jk) = 6
            ELSE IF (zrh<0.985) THEN
              iaerh(jl,jk) = 7
            ELSE
              iaerh(jl,jk) = 8
            END IF
          END DO
        END DO
      ELSE
        ! Fixed 80% r.h. class (default):
        DO jk = 1, nlev
          DO jl = 1, nglon
            iaerh(jl,jk) = 4
          END DO
        END DO
      END IF
    END IF

!-- 3.3 Half-level temperatures

    DO jk = 2, nlev
      DO jl = 1, nglon
        zhti(jl,jk) = (tm1(jl,jk-1,jrow)*apm1(jl,jk-1)*(apm1(jl,jk)-aphm1(jl,jk)) &
            +tm1(jl,jk,jrow)*apm1(jl,jk)*(aphm1(jl,jk)-apm1(jl,jk-1)))            &
            *(1./(aphm1(jl,jk)*(apm1(jl,jk)-apm1(jl,jk-1))))
      END DO
    END DO

    DO jl = 1, nglon
      zhti(jl,nlevp1) = tsm1m(jl,jrow)
      zhti(jl,1) = tm1(jl,1,jrow) - apm1(jl,1)*(tm1(jl,1,jrow)-zhti(jl,2))/ &
                  (apm1(jl,1)-aphm1(jl,2))
    END DO

!-- 3.4 Clouds

    DO jk = 1, nlev
      DO jl = 1, nglon
        zclc(jl,jk) = aclcm(jl,jk,jrow)
        zclc(jl,jk) = MIN(MAX(zclc(jl,jk),zepclc),1.-zepclc)
        zclwa(jl,jk) = MAX(x(jl,jk,jrow),0.)
      END DO
    END DO

!DIR$ IVDEP
!OCL NOVREC
    DO jl = 1, nglon
      aclcvm(jl,jrow) = 1. - aclcm(jl,1,jrow)
    END DO
    DO jk = 2, nlev
!DIR$ IVDEP
!OCL NOVREC
      DO jl = 1, nglon
        aclcvm(jl,jrow) = aclcvm(jl,jrow)*(1.-MAX(aclcm(jl,jk,jrow),aclcm(jl,jk-1,jrow)))/ &
                         (1.-MIN(aclcm(jl,jk-1,jrow),1.-zepsec))
      END DO
    END DO
    DO jl = 1, nglon
      aclcvm(jl,jrow) = 1. - aclcvm(jl,jrow)
    END DO

!-- 3.5 Input diagnostics for temperature and clouds

    IF (lodia) THEN
      zzdi = SUM(tsm1m(1:nglon,jrow))
      zdiag(1,9) = zbud*zzdi
      zdiag(1,1) = zzdi
      DO jk = 1, nlev
        jkk = nlev + 2 - jk
        zzdi = SUM(tm1(1:nglon,jk,jrow))
        zdiag(jkk,9) = zbud*zzdi
        zdiag(jkk,1) = zzdi
      END DO
      zzdi = SUM(aclcvm(1:nglon,jrow))
      zdiag(1,10) = zbud*zzdi
      zdiag(1,3) = zdiag(1,3) + zzdi
      DO jk = 1, nlev
        jkk = nlev + 2 - jk
        zzdi = SUM(zclc(1:nglon,jk))
        zdiag(jkk,10) = zbud*zzdi
        zdiag(jkk,3) = zdiag(jkk,3) + zzdi
      END DO

!DIR$ IVDEP
!OCL NOVREC
      DO jk = nlevp1, 1, -1
        diag(jk,1) = diag(jk,1) + zdiag(jk,9)
        diag(jk,2) = diag(jk,2) + zdiag(jk,10)
      END DO

    END IF

!-- 3.6 Albedo

    zalte = 1.0 - cemiss
!DIR$ IVDEP
!OCL NOVREC
    DO jl = 1, nglon
      zalso(jl) = albm(jl,jrow)
      IF (loland(jl)) THEN
        ztsnmelt = tmelt
        IF (loglac(jl)) THEN
          zalbmin = 0.6
          zalbmax = 0.8
        ELSE
          zalbmin = (1.-forestm(jl,jrow))*zalbmn0 + forestm(jl,jrow)*zalbmn1
          zalbmax = (1.-forestm(jl,jrow))*zalbmx0 + forestm(jl,jrow)*zalbmx1
        END IF
        IF (tsm1m(jl,jrow)>=ztsnmelt) THEN
          zalbsn = zalbmin
        ELSE IF (tsm1m(jl,jrow)<ztalb) THEN
          zalbsn = zalbmax
        ELSE
          zdalb  = (zalbmax-zalbmin)/(ztsnmelt-ztalb)
          zalbsn = zalbmin + zdalb*(ztsnmelt-tsm1m(jl,jrow))
        END IF
        IF (lamip2) THEN
          zalso(jl) = zalso(jl) + (zalbsn-zalso(jl))*(1.-EXP(-snm1m(jl,jrow)*cqsncr))
        ELSE
          zalso(jl) = zalso(jl) + (zalbsn-zalso(jl))*(snm1m(jl,jrow)/(snm1m(jl,jrow)+zsnowal))
        END IF
      ELSE IF (( .NOT. loland(jl)) .AND. seaice(jl,jrow)>0.5) THEN        
        IF (lamip2) THEN
          ztsnmelt = ctfreez
        ELSE
          ztsnmelt = tmelt
        END IF
        zalbmin = 0.5
        zalbmax = 0.75
        IF (tsm1m(jl,jrow)>=ztsnmelt) THEN
          zalso(jl) = zalbmin
        ELSE IF (lamip2 .AND. tsm1m(jl,jrow)<ztsalb) THEN
          zalso(jl) = zalbmax
        ELSE IF (.NOT. lamip2 .AND. tsm1m(jl,jrow)<ztalb) THEN
          zalso(jl) = zalbmax
        ELSE
          IF (lamip2) THEN
            zdalb = (zalbmax-zalbmin)/(ztsnmelt-ztsalb)
          ELSE
            zdalb = (zalbmax-zalbmin)/(ztsnmelt-ztalb)
          END IF
          zalso(jl) = zalbmin + zdalb*(ztsnmelt-tsm1m(jl,jrow))
        END IF
      ELSE
        zalso(jl) = calbsea
      END IF
      zalso(jl) = MIN(calbsno,MAX(zepalb,zalso(jl)))
      albedom(jl,jrow) = zalso(jl)

    END DO

    IF (lmidatm) THEN

!  3.9 ozone (from data file).

!  3.9.1 interpolation in time from monthly mean fields

      ! set parameters
      imonth=12
      zmomid=15.*dayl
      zmonthl=30.*dayl     

      ! ------ annual cycle switched on (nmonth=0).

      IF(nmonth.EQ.0) THEN

        CALL cd2dy(ncbase,nobase,iy)

        iday=ncbase+(ntbase+dtime*(nstep+1))/dayl+0.01
        CALL cd2dat(iday,id,im,iy)

        im1 = im-1
        im2 = im
        im3 = im+1

        zyt = MOD(((nobase-1.)+(ntbase*ntimst+nstep*dtime  &
              +0.5*nradfr*dtime)/dayl)/yearl,1.)
        zmosec = MOD(zyt*yearl*dayl,30.*dayl)

        zw1 = MAX(-(zmosec-zmomid),0.)/zmonthl
        zw3 = MAX(  zmosec-zmomid ,0.)/zmonthl
        zw2 = 1.-(zw1+zw3)

        DO jk=1,noz
          scloz(jk) = zw1*ozone(jk,krow,im1) + zw2*ozone(jk,krow,im2)  &
                    + zw3*ozone(jk,krow,im3)
        ENDDO

      ELSE
      ! ------ annual cycle switched off (nmonth.gt.0)

        im = nmonth
        DO jk = 1,noz
          scloz(jk) = ozone(jk,krow,im)
        ENDDO

      ENDIF

! 3.9.2 interpolation in pressure from fixed noz-grid

      DO k = 1,noz-1
        DO jk = 1,nlev 
          DO jl = 1,nglon 
            IF ( apm1(jl,jk).GE.spe(k) .AND. apm1(jl,jk) .LT. spe(k+1) ) THEN
              zkp1 = spe(k+1) - apm1(jl,jk)
              zk = apm1(jl,jk) - spe(k)
              scloza(jl,jk) = (zk*scloz(k+1)+zkp1*scloz(k))/(zk+zkp1)
            ENDIF 
          ENDDO
        ENDDO
      ENDDO

      DO jk = 1,nlev 
        DO jl = 1,nglon 
          IF ( apm1(jl,jk).LT.spe(1) ) scloza(jl,jk) = scloz(1)
          IF ( apm1(jl,jk).GE.spe(noz) ) scloza(jl,jk) = scloz(noz)
        ENDDO
      ENDDO

! 3.9.3 from mass mixing ratio to (kg/kg)*pa 
                           
      DO jk=1,nlev
        DO jl=1,nglon 
          zqof(jl,jk) = scloza(jl,jk)*(aphm1(jl,jk+1)-aphm1(jl,jk))
        ENDDO
      ENDDO

    ELSE IF (lamip2) THEN
      !  ozone for 19 level model
      !  3.9 ozone (from data file).
      !
      !  3.9.1 interpolation in time from monthly mean fields
      !
      ! ------ annual cycle switched on (nmonth=0).

      IF(nmonth == 0) THEN

!CDIR$ SHORTLOOP
        DO jk = 1,noz
          scloz(jk) = wgt1*ozone(jk,krow,nmw1)+wgt2*ozone(jk,krow,nmw2)
        END DO
      ELSE
        ! ------ annual cycle switched off (nmonth.gt.0)
        im = nmonth
        DO jk = 1,noz
          scloz(jk) = ozone(jk,krow,im)
        ENDDO
      ENDIF

      ! 3.9.2 Interpolation in pressure from fixed noz-grid

      !       Interpolation for the first two model layers (p-levels)
      !       which contain the first ten ozone layers using zdpe as 
      !       weighting function.
      !
      !       The weighting function differs from zdpe at the interface
      !         between model level 1 and 2.
      !
      !         (====use only for 19 level model version ====)
      !
      !       vertical structure of the uppermost 2 model layers
      !       and corresponding 10 ozone layers [hpa]:
      !
      !
      !      K           OZONE               DP               MODEL

      ! --------------------------------------------------------------------
      !
      !      0    -----------------          0    ----------------------
      !
      !  1   0.3   - - - - - - - -     0.4
      !
      !      0.4  -----------------
      !
      !  2   0.5  - - - - - - - - -    0.35
      !
      !      0.75 -----------------
      !
      !  3   1    - - - - - - - - -    0.75
      !
      !      1.5 ------------------
      !
      !  4   2   - - - - - - - - -     1
      !
      !      2.5 ------------------
      !
      !  5   3   - - - - - - - - -     1.5
      !
      !      4  -------------------
      !
      !  6   5  - - - - - - - - - -    2
      !
      !      6  -------------------
      !
      !  7   7  - - - - - - - - - -    2.5
      !
      !      8.5 -------------------
      !
      !  8  10  - - - - - - - - - -    6.5   10   - - - - - - - - - - - - - - 
      !
      !     15  ---------------------
      !
      !  9  20  - - - - - - - - - - - 10     20  -----------------------------
      !
      !     25  ---------------------
      !
      ! 10  30  - - - - - - - - - - - 15     30  - - - - - - - - - - - - - - -
      !
      !     40 ----------------------        40  -----------------------------
      !
      !
      !
      !  Interpolation for the first model layer. I.e. integration of 
      !  all ozone layers down to first model layer boundary
      !
      DO jl = 1,nglon
        scloza(jl,1) = 0.
      ENDDO

      DO k = 1,8
        DO jl = 1,nglon
          scloza(jl,1) = scloza(jl,1)+scloz(k)*zdpe(k)
        ENDDO
      ENDDO

      !  Above layer boundary use 5 hpa weight from ozone layer 9

      DO jl = 1,nglon
        scloza(jl,1) = (scloza(jl,1)+scloz(9)*500.)/2000.
      ENDDO

      !  Interpolation for the second model layer
      !  5 hpa from ozone layer 9 and 15 hpa from ozone layer 10

      DO jl = 1,nglon
        scloza(jl,2) = (scloz(9)*500.+scloz(10)*1500.)/2000.
      ENDDO

      ! Interpolation for the remaining layers

      DO k = 10,noz-1
        DO jk = 3,nlev 
          DO jl = 1,nglon 
            IF ( apm1(jl,jk) >= spe(k) .AND. apm1(jl,jk) < spe(k+1) ) THEN
              zkp1 = spe(k+1) - apm1(jl,jk)
              zk = apm1(jl,jk) - spe(k)
              scloza(jl,jk) = (zk*scloz(k+1)+zkp1*scloz(k))/(zk+zkp1)
            ENDIF 
          ENDDO
        ENDDO
      ENDDO

      DO jk = 3,nlev 
        DO jl = 1,nglon 
          !!! IF (apm1(jl,jk)<spe(1)) scloza(jl,jk)=zcloz(1)
          IF (apm1(jl,jk)>=spe(noz)) scloza(jl,jk)=scloz(noz)
        ENDDO
      ENDDO

     ! 3.9.3 from volume mixing ratio to mass mixing ratio * dp
                           
      DO jk = 1,nlev
        DO jl = 1,nglon 
          zqof(jl,jk) = scloza(jl,jk)*zconvoz*(aphm1(jl,jk+1)-aphm1(jl,jk))
        ENDDO
      ENDDO
    ENDIF

!-- 4. Fourier transform for ozone and aerosols

    DO jl = 1, nglon
      zcos1 = coslon(jl+nlof)
      zsin1 = sinlon(jl+nlof)
      zcos2 = zcos1*zcos1 - zsin1*zsin1
      zsin2 = zsin1*zcos1 + zcos1*zsin1
      zcos3 = zcos2*zcos1 - zsin2*zsin1
      zsin3 = zsin2*zcos1 + zcos2*zsin1
      zcos4 = zcos3*zcos1 - zsin3*zsin1
      zsin4 = zsin3*zcos1 + zcos3*zsin1
      zcos5 = zcos4*zcos1 - zsin4*zsin1
      zsin5 = zsin4*zcos1 + zcos4*zsin1

      !  -----> Fourier transformation for ozone

      IF (.NOT.lmidatm .AND. .NOT.lamip2) THEN
        zozq(jl) = zfozq(1) + 2.*(zfozq(2)*zcos1+zfozq(3)*zsin1+zfozq(4)*zcos2+ &
                   zfozq(5)*zsin2+zfozq(6)*zcos3+zfozq(7)*zsin3+zfozq(8)*zcos4+ &
                   zfozq(9)*zsin4+zfozq(10)*zcos5+zfozq(11)*zsin5)
        zzozh = zfozh(1) + 2.*(zfozh(2)*zcos1+zfozh(3)*zsin1+zfozh(4)*zcos2+ &
                zfozh(5)*zsin2+zfozh(6)*zcos3+zfozh(7)*zsin3+zfozh(8)*zcos4+ &
                zfozh(9)*zsin4+zfozh(10)*zcos5+zfozh(11)*zsin5)
        zozh(jl) = SQRT(zzozh)**3
      ENDIF

      zcos6 = zcos5*zcos1 - zsin5*zsin1
      zsin6 = zsin5*zcos1 + zcos5*zsin1
      zcos7 = zcos6*zcos1 - zsin6*zsin1
      zsin7 = zsin6*zcos1 + zcos6*zsin1
      zcos8 = zcos7*zcos1 - zsin7*zsin1
      zsin8 = zsin7*zcos1 + zcos7*zsin1
      zcos9 = zcos8*zcos1 - zsin8*zsin1
      zsin9 = zsin8*zcos1 + zcos8*zsin1
      zcos10 = zcos9*zcos1 - zsin9*zsin1
      zsin10 = zsin9*zcos1 + zcos9*zsin1

      zaes(jl) = zfaes(1) + 2.*(zfaes(2)*zcos1+zfaes(3)*zsin1+zfaes(4)*zcos2+      &
                 zfaes(5)*zsin2+zfaes(6)*zcos3+zfaes(7)*zsin3+zfaes(8)*zcos4+      &
                 zfaes(9)*zsin4+zfaes(10)*zcos5+zfaes(11)*zsin5+zfaes(12)*zcos6+   &
                 zfaes(13)*zsin6+zfaes(14)*zcos7+zfaes(15)*zsin7+zfaes(16)*zcos8+  &
                 zfaes(17)*zsin8+zfaes(18)*zcos9+zfaes(19)*zsin9+zfaes(20)*zcos10+ &
                 zfaes(21)*zsin10)
      zael(jl) = zfael(1) + 2.*(zfael(2)*zcos1+zfael(3)*zsin1+zfael(4)*zcos2+      &
                 zfael(5)*zsin2+zfael(6)*zcos3+zfael(7)*zsin3+zfael(8)*zcos4+      &
                 zfael(9)*zsin4+zfael(10)*zcos5+zfael(11)*zsin5+zfael(12)*zcos6+   &
                 zfael(13)*zsin6+zfael(14)*zcos7+zfael(15)*zsin7+zfael(16)*zcos8+  &
                 zfael(17)*zsin8+zfael(18)*zcos9+zfael(19)*zsin9+zfael(20)*zcos10+ &
                 zfael(21)*zsin10)
      zaeu(jl) = zfaeu(1) + 2.*(zfaeu(2)*zcos1+zfaeu(3)*zsin1+zfaeu(4)*zcos2+      &
                 zfaeu(5)*zsin2+zfaeu(6)*zcos3+zfaeu(7)*zsin3+zfaeu(8)*zcos4+      &
                 zfaeu(9)*zsin4+zfaeu(10)*zcos5+zfaeu(11)*zsin5+zfaeu(12)*zcos6+ &
                 zfaeu(13)*zsin6+zfaeu(14)*zcos7+zfaeu(15)*zsin7+zfaeu(16)*zcos8+ &
                 zfaeu(17)*zsin8+zfaeu(18)*zcos9+zfaeu(19)*zsin9+zfaeu(20)*zcos10+ &
                 zfaeu(21)*zsin10)
      zaed(jl) = zfaed(1) + 2.*(zfaed(2)*zcos1+zfaed(3)*zsin1+zfaed(4)*zcos2+      &
                 zfaed(5)*zsin2+zfaed(6)*zcos3+zfaed(7)*zsin3+zfaed(8)*zcos4+      &
                 zfaed(9)*zsin4+zfaed(10)*zcos5+zfaed(11)*zsin5+zfaed(12)*zcos6+   &
                 zfaed(13)*zsin6+zfaed(14)*zcos7+zfaed(15)*zsin7+zfaed(16)*zcos8+  &
                 zfaed(17)*zsin8+zfaed(18)*zcos9+zfaed(19)*zsin9+zfaed(20)*zcos10+ &
                 zfaed(21)*zsin10)

    END DO

!-- 5.1 Input: co2, o3 and aerosols
!       (aerosol is stored upside down)

    DO jl = 1, nglon
      zaetro(jl) = 1.
    END DO
    DO jk = 1, nlev
      jkl = nlev + 1 - jk
      DO jl = 1, nglon
        IF (.NOT.lmidatm .AND. .NOT.lamip2) THEN
          zcpho3 = aphm1(jl,jk)**3
          zcphn3 = aphm1(jl,jk+1)**3
          zsdpo3 = SQRT(zcpho3)
          zsdpn3 = SQRT(zcphn3)
        ENDIF
        zaeqso = caeops*zaes(jl)*cvdaes(jk)
        zaeqsn = caeops*zaes(jl)*cvdaes(jk+1)
        zaeqlo = caeopl*zael(jl)*cvdael(jk)
        zaeqln = caeopl*zael(jl)*cvdael(jk+1)
        zaequo = caeopu*zaeu(jl)*cvdaeu(jk)
        zaequn = caeopu*zaeu(jl)*cvdaeu(jk+1)
        zaeqdo = caeopd*zaed(jl)*cvdaed(jk)
        zaeqdn = caeopd*zaed(jl)*cvdaed(jk+1)
        zaetrn(jl) = zaetro(jl)*(MIN(1.,zhti(jl,jk)/zhti(jl,jk+1)))**ctrpt
        zaetr = SQRT(zaetrn(jl)*zaetro(jl))
        IF (.NOT.lmidatm .AND. .NOT.lamip2) THEN
          zqofo = zozq(jl)*zsdpo3/(zsdpo3+zozh(jl))
          zqofn = zozq(jl)*zsdpn3/(zsdpn3+zozh(jl))
          zqof(jl,jk) = zqofn - zqofo
        ENDIF
        zdp(jl,jk) = aphm1(jl,jk+1) - aphm1(jl,jk)
        zsaer(jl,jkl,1) = (1.-zaetr)*(ctrbga*zdp(jl,jk)+zaeqln-zaeqlo+zaeqdn- &
                          zaeqdo)
        zsaer(jl,jkl,2) = (1.-zaetr)*(zaeqsn-zaeqso)
        zsaer(jl,jkl,3) = (1.-zaetr)*(zaequn-zaequo)
        zsaer(jl,jkl,4) = zaetr*cvobga*zdp(jl,jk)
        zsaer(jl,jkl,5) = zaetr*cstbga*zdp(jl,jk)
        zaetro(jl) = zaetrn(jl)
      END DO
    END DO

    ! Aerosol switch

    IF (laer) THEN
      kaer = 1
    ELSE
      kaer = 0
    END IF

    ! Kaer is the multiplication factor for aerosol optical depths

    ! kaer=0:  Set tanre aerosols to zero  (zsaer( , ,15) )
    ! kaer=1:  Use zsaer values except for volcanic type (index 4)

    IF (kaer==1) THEN
      DO jk = 1, nlev
        DO jl = 1, nglon
          zsaer(jl,jk,1) = MAX(zsaer(jl,jk,1),caeros)
          zsaer(jl,jk,2) = MAX(zsaer(jl,jk,2),caeros)
          zsaer(jl,jk,3) = MAX(zsaer(jl,jk,3),caeros)
          zsaer(jl,jk,4) = caeros
          zsaer(jl,jk,5) = MAX(zsaer(jl,jk,5),caeros)
        END DO
      END DO
    ELSE
      DO jaer = 1, 5
        DO jk = 1, nlev
          DO jl = 1, nglon
            zsaer(jl,jk,jaer) = 0.
          END DO
        END DO
      END DO
    END IF

    ! Set additional aerosols (optional code example)

    ! (The aerosol *zsaer* is stored upside down)

    ! DO ja=1,newaer
    ! DO jk=1,nlev
    ! jkl=nlev+1-jk
    ! DO jl=1,nglon
    ! zsaer(jl,jkl,5+ja)=xtm1(jl,jk,ja)  ! mass mixing ratio
    ! END DO
    ! END DO
    ! END DO

    ! Switch for solar clear sky diagnostic

    IF (lsolc) THEN
      kmode = 1
    ELSE
      kmode = 0
    END IF

    ! Cfc switch

    IF (lcfc) THEN
      kcfc = 1
    ELSE
      kcfc = 0
    END IF

!-- 5.2 Call to *radlsw*

    CALL radlsw(nglon,nlev,nglpx,kewaer,iaerh,kmode,kaer,kcfc,zsct,loland,         &
               loglac,zsaer,zalso,zclc,zmu0,zqof,aphm1(1,nlevp1),zdp,zclwa,zq,zqs, &
               aphm1,tm1(:,:,jrow),zhti,zfls,zflt,zflsc,zfltc)


!-- 5.3 Storage of the output

    DO jk = 1, nlevp1
!DIR$ IVDEP
!OCL NOVREC
      DO jl = 1, nglon
        trsolm(jl,jk,jrow) = zfls(jl,jk)/(zsct*zmu0(jl))
        IF (lmidatm) THEN
          emterm(jl,jk,jrow) = zflt(jl,jk)
        ELSE
          emterm(jl,jk,jrow) = zflt(jl,jk)/(stbo*zhti(jl,jk)**4)
        ENDIF
      END DO
    END DO

    ! Clear sky fluxes

!DIR$ IVDEP
!OCL NOVREC
    DO jl = 1, nglon
      IF (lmidatm) THEN
        emtefm(jl,1,jrow) = zfltc(jl,1)
        emtefm(jl,2,jrow) = zfltc(jl,2)
      ELSE
        emtefm(jl,1,jrow) = zfltc(jl,1)/(stbo*zhti(jl,1)**4)
        emtefm(jl,2,jrow) = zfltc(jl,2)/(stbo*zhti(jl,nlevp1)**4)
      ENDIF
      trsofm(jl,1,jrow) = zflsc(jl,1)/(zsct*zmu0(jl))
      trsofm(jl,2,jrow) = zflsc(jl,2)/(zsct*zmu0(jl))
    END DO

!-- 6. Output diagnostics for heating rates and fluxes

    IF (lodia) THEN
      DO jl = 1, nglon
        zdiaf(jl,1) = zrdayl(jl)*zsct*zmu0(jl)
        zdiaf(jl,2) = zdiaf(jl,1)*(1.-trsolm(jl,1,jrow))
        IF (lmidatm) THEN
          zdiaf(jl,3) = -emterm(jl,1,jrow)
        ELSE
          zdiaf(jl,3) = -emterm(jl,1,jrow)*stbo*zhti(jl,1)**4
        ENDIF
        zdiaf(jl,4) = zrdayl(jl)*zmu0(jl)*zsct*trsolm(jl,nlevp1,jrow)/ &
                     (1.-zalso(jl))
        zdiaf(jl,5) = zdiaf(jl,4) - zrdayl(jl)*zmu0(jl)*zsct*trsolm(jl,nlevp1,jrow)
        IF (lmidatm) THEN
          zdiaf(jl,6) = stbo*zhti(jl,nlevp1)**4+emterm(jl,nlevp1,jrow)*zalte/cemiss
          zdiaf(jl,7) = zdiaf(jl,6) + emterm(jl,nlevp1,jrow)
        ELSE
          zdiaf(jl,6) = stbo*zhti(jl,nlevp1)**4*(1.+emterm(jl,nlevp1,jrow)*zalte/cemiss)
          zdiaf(jl,7) = zdiaf(jl,6) + stbo*emterm(jl,nlevp1,jrow)*zhti(jl,nlev+1)**4
        ENDIF
      END DO
      DO jn = 1, 7
        zzdi = SUM(zdiaf(1:nglon,jn))
        zdia1(jn) = zbud*zzdi
      END DO
      DO jl = 1, nglon
        zdia(jl) = zdegday*zrdayl(jl)*zmu0(jl)*zsct* &
            (trsolm(jl,1,jrow)-trsolm(jl,nlevp1,jrow))/(aphm1(jl,nlevp1)-aphm1(jl,1))
      END DO
      zzdi = SUM(zdia(1:nglon))
      zdiag(1,9) = zbud*zzdi
      zdiag(1,5) = zdiag(1,5) + zzdi
      IF (lmidatm) THEN
        DO jl = 1, nglon
          zdia(jl) = zdegday*(emterm(jl,1,jrow)-emterm(jl,nlevp1,jrow)) &
                    /(aphm1(jl,nlevp1)-aphm1(jl,1))
        END DO
      ELSE
        DO jl = 1, nglon
          zdia(jl) = zdegday*stbo*(emterm(jl,1,jrow)*zhti(jl,1)**4-emterm(jl,nlevp1,jrow) &
                    *zhti(jl,nlevp1)**4)/(aphm1(jl,nlevp1)-aphm1(jl,1))
        END DO
      ENDIF
      zzdi = SUM(zdia(1:nglon))
      zdiag(1,10) = zbud*zzdi
      zdiag(1,7) = zdiag(1,7) + zzdi
      DO jk = 1, nlev
        jkk = nlev + 2 - jk
        DO jl = 1, nglon
          zdia(jl) = zdegday*zrdayl(jl)*zmu0(jl)*zsct* &
                  (trsolm(jl,jk,jrow)-trsolm(jl,jk+1,jrow))/(aphm1(jl,jk+1)-aphm1(jl,jk))
        END DO
        zzdi = SUM(zdia(1:nglon))
        zdiag(jkk,9) = zbud*zzdi
        zdiag(jkk,5) = zdiag(jkk,5) + zzdi
        IF (lmidatm) THEN
          DO jl = 1, nglon
            zdia(jl)=zdegday*(emterm(jl,jk,jrow)-emterm(jl,jk+1,jrow)) &
                    /(aphm1(jl,jk+1)-aphm1(jl,jk))
          END DO
        ELSE
          DO jl = 1, nglon
            zdia(jl) = zdegday*stbo*(emterm(jl,jk,jrow)*zhti(jl,jk)**4- &
                 emterm(jl,jk+1,jrow)*zhti(jl,jk+1)**4)/(aphm1(jl,jk+1)-aphm1(jl,jk))
          END DO
        ENDIF
        zzdi = SUM(zdia(1:nglon))
        zdiag(jkk,10) = zbud*zzdi
        zdiag(jkk,7) = zdiag(jkk,7) + zzdi
      END DO

      DO jn = 1, 7
        dia1(jn) = dia1(jn) + zdia1(jn)
      END DO

!DIR$ IVDEP
!OCL NOVREC
      DO jk = nlevp1, 1, -1
        diag(jk,3) = diag(jk,3) + zdiag(jk,9)
        diag(jk,4) = diag(jk,4) + zdiag(jk,10)
      END DO

    END IF

!-- 7. Necessary computations if subroutine is by-passed

  ELSE
    DO jk = 1, nlevp1
      DO jl = 1, nglon
        trsolm(jl,jk,jrow) = 0.
        emterm(jl,jk,jrow) = 0.
      END DO
    END DO
    DO jk = 1, 2
      DO jl = 1, nglon
        emtefm(jl,jk,jrow) = 0.
        trsofm(jl,jk,jrow) = 0.
      END DO
    END DO
  END IF

!-- 8. Print diagnostics if necessary and return workspace

  IF (lodiap) THEN
    zlatd = 180.*ASIN(.5*twomu(irow))/api
    WRITE (nout, '(/,a,f6.1,/,a,/)')                        &
        ' Radiation zonally averaged for latitude:', zlatd, &
        '   (Fluxes - T=Top, B=Bottom, D=Down, U=Up, S=Solar, L=Thermal, N=Net)'
    DO jn = 1, 7, 2
      DO jk = 1, nlevp1
        zdiag(jk,jn) = zdiag(jk,jn)*zdials
      END DO
    END DO
    DO jn = 1, 7
      zdia1(jn) = zdia1(jn)*zdials/zbud
    END DO
    DO jk = 1, nlevp1
      zdiat(jk) = zdiag(jk,5) + zdiag(jk,7)
    END DO
    zdift1 = zdia1(1) - zdia1(2)
    zdift2 = zdia1(4) - zdia1(5)
    zdift3 = zdia1(6) - zdia1(7)
    zdift4 = zdift1 - zdia1(3)
    zdift5 = zdift2 - zdift3
    IF (zdia1(1)>0.01) THEN
      zalbpr = zdia1(2)/zdia1(1)*100.
    ELSE
      zalbpr = 0.
    END IF
    WRITE (nout, '(a)') ' Temperature  [K]'
    WRITE (nout, '(10f6.1)') zdiag(1:nlev,1)
    WRITE (nout, '(a)') ' Cloud Cover     '
    WRITE (nout, '(10f6.1)') zdiag(1:nlev,2)
    WRITE (nout, '(a)') ' Shortwave heating rates [K/day]  '
    WRITE (nout, '(10f6.1)') zdiag(1:nlev,3)
    WRITE (nout, '(a)') ' Longwave heating rates [K/day]  '
    WRITE (nout, '(10f6.1)') zdiag(1:nlev,4)
    WRITE (nout, '(a)') ' Net heating rates [K/day]  ' 
    WRITE (nout, '(10f6.1)') zdiat(1:nlev)

    WRITE (nout, '(a)') ' Fluxes:'
    WRITE (nout, '(3(a,f6.1))') &
        '   TSD = ', zdia1(1), ' TSU = ', zdia1(2), ' TLU = ', zdia1(3)
    WRITE (nout, '(4(a,f6.1))') &
        '   BSD = ', zdia1(4), ' BSU = ', zdia1(5), &
        ' BLU = ', zdia1(6), ' BLD = ', zdia1(7)
    WRITE (nout, '(5(a,f6.1))') '   TS  = ', zdift1, &
        ' BS  = ',zdift2 ,' BL  = ', zdift3, ' TND = ', zdift4,' BND = ', zdift5
    WRITE (nout, '(a,f5.1,a)') '   Albedo = ', zalbpr, ' %'
    WRITE (nout, '(2/)')

  END IF

  RETURN
END SUBROUTINE radint
