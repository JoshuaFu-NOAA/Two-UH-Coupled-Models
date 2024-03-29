!+ computes large-scale water phase changes and cloud cover.

!#define DEBUG

!MODULE m_cond
!
!CONTAINS

SUBROUTINE cond(kidia,kfdia,klon,klp2,ktdia,klev,klevp1,klab,paclcacm,paphm1,         &
                paphp1,papm1,papp1,pgeom1,pqm1,ptm1,pxm1,pxtec,laland,ktype,paclcovm, &
                palwcvim,paprlm,pqvim,paclc,paclcac,paclcov,palwcvi,paprl,pqvi,pssfl, &
                ptte,pqte,pxte,paprs,prsfl)

  ! Description:
  !
  ! Computes large-scale water phase changes and cloud cover.
  !
  ! Method:
  !
  ! This routine computes the physical tendencies of the three
  ! prognostic variables t, q, and x (cloud water) due to water phase
  ! changes (condensation, evaporation of falling precipitation in
  ! unsaturated layers and melting or freezing of the falling water)
  ! and precipitaion formation (coalescense, sedimentation).
  !
  ! Rain, snowfall, surface fluxes, cloud water and cloud cover
  ! later to be used for soil processes and radiation are stored.
  !
  ! References:
  !
  ! 1. Large scale phase changes' part of the model's documentation.
  ! 2. Roeckner, E. and U. Schlese (1985), ECMWF-workshop on
  !    "cloud cover parameterization in numerical models", pp87-108.
  ! 3. Roeckner et al. (1991) ECMWF/WRCP workshop on
  !    "clouds, radiative transfer and the hydrological cycle",
  !    199-222, ECMWF, Reading, u.k.
  ! 4. Sundqvist, H. (1978), qjrms, 104, 677-690.
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, November 1992, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! A. Rhodin, MPI, December 1998, lookup tables removed
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_constants,         ONLY: als, alv, c2es, c4ies, c4les, c5ies, &
                                  c5les, cpd, g, rd, rhoh2o, tmelt,    &
                                  vtmpc1, vtmpc2
  USE mo_diagnostics_zonal, ONLY: dlserz, dlsesz, dlsgrz, dlsgsz, dlsmsz
  USE mo_control,           ONLY: lmidatm, nn, nrow, twodt, lamip2, ldsst
  USE mo_gaussgrid,         ONLY: budw
  USE mo_start_dataset,     ONLY: nstart, nstep
  USE mo_machine,           ONLY: prec
#ifndef NOLOOKUP
  USE mo_convect_tables,    ONLY: tlucua, jptlucu1, jptlucu2, &
                                  lookuperror, lookupoverflow
#else
  USE mo_constants,         ONLY: c3les, c3ies
#endif
  USE mo_exception,         ONLY: finish
!fu++
  USE mo_dsst,              ONLY: prlarg_diurnal

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: kfdia, kidia, klev, klevp1, klon, klp2, ktdia

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: paclcacm(klp2,klev), paclcovm(klp2),   &
                       palwcvim(klp2), paphm1(klp2,klevp1),   &
                       paphp1(klp2,klevp1), papm1(klp2,klev), &
                       papp1(klp2,klev), paprlm(klp2),        &
                       pgeom1(klp2,klev), pqm1(klp2,klev),    &
                       pqvim(klp2), ptm1(klp2,klev),          &
                       pxm1(klp2,klev), pxtec(klp2,klev)
  INTEGER, INTENT (IN) :: klab(klp2,klev), ktype(klp2)
  LOGICAL, INTENT (IN) :: laland(klp2)

  ! klab     : convection flag (0: no convection, 1:  
  ! paclcacm : cloud cover, accumulated (old value)
  ! paphm1   : pressure at half levels (t-dt)
  ! paphp1   : pressure at half levels (t+dt)
  ! papm1    : pressure at full levels (t-dt)
  ! papp1    : pressure at full levels (t+dt)
  ! pgeom1   : geopotential at full levels (t-dt)
  ! pqm1     : specific humidity (t-dt)
  ! ptm1     : temperature (t-dt)
  ! pxm1     : cloud water (t-dt)
  ! pxtec    : tendency of detrained convective cloud water
  ! ktype    : type of convection
  ! laland   : logical mask for land
  ! paclcovm : total cloud cover, accumulated (old value)
  ! palwcvim : vertically integrated cloud water, accumulated (old value)
  ! paprlm   : stratiform precipitaion, accumulated (old value)
  ! pqvim    : vertically integrated spec. humidity, accumulated (old value)

  !  Array arguments with intent(InOut):
  REAL, INTENT (INOUT) :: ptte(klp2,klev), pqte(klp2,klev), pxte(klp2,klev),  &
                          paprs(klp2), prsfl(klp2)

  ! ptte     : tendency of temperature
  ! pqte     : tendency of specific humidity
  ! pxte     : tendency of cloud water
  ! paprs    : snow fall, accumulated
  ! prsfl    : surface rain flux

  !  Array arguments with intent(Out):
  REAL, INTENT (OUT) :: paclc(klp2,klev), paclcac(klp2,klev),       &
                        paclcov(klp2), palwcvi(klp2), paprl(klp2),  &
                        pqvi(klp2), pssfl(klp2)
  ! paclc    : cloud cover
  ! paclcac  : cloud cover, accumulated
  ! paclcov  : total cloud cover
  ! palwcvi  : vertically integrated cloud water, accumulated
  ! paprl    : stratiform precipitation, accumulated
  ! pqvi     : vertically integrated spec. humidity, accumulated
  ! pssfl    : surface snow flux

  !  Local scalars: 
  REAL :: rdcpd, z2es, zaci, zacl, zalpha, zalwcold, zbud, zc1, zcaa, zcab,     &
          zclcm1, zclear, zcloia, zclp1, zcnull, zcons, zcons1, zcons2, zconsi, &
          zcor, zcpn, zcvm4, zcvm5, zcvml, zdaci, zdacl, zdelp, zdiagt, zdiagw, &
          zdnopr, zdqdt, zdt, zdtdt, zdthdp, zdxcor, zdxdt, zepclc, zepflm,     &
          zepqp1, zepsec, zepzwp, zfrac, zl, zmeltp2, zpcc, zqpa, zqq, zqr,     &
          zqsm1, zqtmst, zrcp, zrflo, zrho, zrhsc, zrs, zrtc, zrtl, zs, zsatsc, &
          zsflo, zsnmlt, ztc, ztmsc, ztmst, zvtfac, zwi, zwma, zwmax, zwmvn,    &
          zwpa, zwps, zwpvn, zwpvni, zwpvnl, zwrc, zwrl, zwrs, zwrsc, zxetmst,  &
          zzdr, zzep, zzepcons, zzfac, zzqt, zzrh, zztc

  !------------------------------------------------------------------------------
  !  Local scalars new for AMIP2 and the modifications of Lenderink et al., 1998   
  !
  REAL :: zepclw, zanvcc, zqconm1, zlvcdqsdt, zdqsdt, zclcaux, zqldt,           &
          zqldtstar, zqvdt, zdtdtstar,zdqsat, zevape, zcondc, ztp1tmp, zqp1tmp, &
          zqsp1tmp, zoversat, zdconratcor, zdxc
  !------------------------------------------------------------------------------


  INTEGER :: irow, jb, jk, jl, klev2, klev2p1, klevm1, nexc, nexl, klev_m, jrow
#ifndef NOLOOKUP
  INTEGER :: it
#endif
  LOGICAL :: lo, lo1

  !  Local arrays: 
  REAL :: zalwc(klon,klev), zbas(klon,klev), zclcp1(klon,klev),                   &
          zconrat(klon), zdac(klon), zdp(klon), zdthmin(klon), zep(klon),         &
          zepr(klon), zeps(klon), zfln(klon), zgeoh(klon,klevp1),                 &
          zlsdcp(klon,klev), zlvdcp(klon,klev), zqcon(klon,klev),                 &
          zqp1(klon,klev), zqsp1(klon,klev), zrain(klon), zrfl(klon),             &
          zrfln(klon), zrhc(klon,klev), zsat(klon,klev), zsfl(klon), zsfln(klon), &
          zsnmt(klon), zsnow(klon), ztheta(klon,klev), ztop(klon,klev),           &
          ztp1(klon,klev)
  INTEGER :: invb(klon)

  !  Intrinsic functions 
#ifdef ECLIB
  REAL :: EXPHF
!DIR$ VFUNCTION EXPHF
#define EXP(x)  EXPHF(x)
#else
  INTRINSIC EXP
#endif
  INTRINSIC MAX, MERGE, MIN, SQRT, SUM


  !  Executable statements


#ifndef NOLOOKUP
  lookupoverflow = .FALSE.
#endif

  ! Security parameters.

  ! zepflm is a minimum flux to avoid dividing by zero in the ice
  ! proportion calculations.

  zepflm = 1.E-24
  zepsec = 1.E-12
  zepqp1 = 0.
  zepzwp = 1.E-20
  zepclc = 1.E-2
  zdnopr = 1.E4
  zepclw = 1.E-7   ! new for AMIP2

  ! *    Computational constants.

  ztmst = twodt
  IF (nstep==nstart) ztmst = 0.5*twodt
  zqtmst = 1./ztmst
  zdiagt = 0.5*twodt
  zdiagw = zdiagt/rhoh2o

  zcons1  = cpd*vtmpc2
  zcons2  = 1./(ztmst*g)
  zmeltp2 = tmelt + 2.
  z2es    = c2es
  zdelp   = 1.E4
  zvtfac  = 3.29

  ! Attention: zsatsc>zrhsc

  zrs = 0.99

  IF (lamip2) THEN
    ztmsc    = 0.85
    zcnull   = 1.0E-4*ztmst
    zzepcons = 0.001
    zrhsc    = 0.85
    zsatsc   = 0.85
    zrtc     = 0.7
  ELSE
    ztmsc    = 0.60
    zcnull   = 2.0E-4*ztmst
    zzepcons = 0.015
    zrhsc    = 0.4
    zsatsc   = 0.8
    zrtc     = 0.6
  END IF

  zrtl = 0.6
  nexc = 4
  nexl = 4
  IF (nn == 106) THEN
    ztmsc    = 0.65
    zcnull   = 2.5E-4*ztmst
    zzepcons = 0.01
    zrhsc    = 0.6
    zsatsc   = 0.9
  END IF
  IF (nn == 30) THEN
    IF (lmidatm) THEN
      ztmsc    = 0.75
      zzepcons = 0.01
      zcnull   = 2.5e-4*ztmst
    ELSE
      ztmsc  = 0.56
      zcnull = 1.5E-4*ztmst
    END IF
  END IF
  IF (nn == 21) THEN
    IF (lmidatm) THEN
      ztmsc    = 0.65
      zzepcons = 0.01
    ELSE
      ztmsc    = 0.50
      zcnull   = 1.E-4*ztmst
      zzepcons = 0.02
      zrhsc    = 0.4
      zsatsc   = 0.8
    END IF
  END IF

  zconsi = zvtfac*g*ztmst*ztmsc/zdelp
  zcaa   = 0.0059
  zcab   = 0.003102

  IF (lamip2) THEN
    zc1  = ztmst
  ELSE
    zc1  = 2.*ztmst
  END IF

  zwrsc  = 0.5E-3
  zwrl   = 0.5E-3

  IF (lamip2) THEN
    zwrs = 0.5e-3
  ELSE
    zwrs = 0.3E-3
  END IF

  zanvcc = zepclc   ! new for AMIP2

  klevm1  = klev - 1
  klev2   = klev/2
  klev2p1 = klev2 + 1
  rdcpd   = rd/cpd
  zcloia  = 1.0E+02

  irow = nrow(1)
  jrow = nrow(2)

!-- 1. Top boundary conditions and quantities needed for
!      --- -------- ---------- --- ---------- ------ ---
!      condensation and precipitation calculations.
!      ------------ --- ------------- ------------

!-- 1.1 Set to zero precipitation fluxes at the top

  DO jl = kidia, kfdia
    zrfl(jl)    = 0.
    zsfl(jl)    = 0.
    zrain(jl)   = 0.
    zsnow(jl)   = 0.
    zepr(jl)    = 0.
    zeps(jl)    = 0.
    zsnmt(jl)   = 0.
    zrfln(jl)   = 0.
    zsfln(jl)   = 0.
    zdthmin(jl) = 0.0
    invb(jl)    = 1
  END DO
  zs = 0.

!-- 1.2 Calculate potential temperatures

  DO jk = klev, klev2, -1
    DO jl = kidia, kfdia
      ztheta(jl,jk) = ptm1(jl,jk)*(1.0E5/papm1(jl,jk))**rdcpd
    END DO
  END DO

!-- 1.3 Check for presence of low-level inversion
!       ( sea points only )

  IF (lamip2) THEN
    klev_m = 15
  ELSE
    klev_m = klev2p1
  END IF
  DO jk = klev, klev_m, -1
    DO jl = kidia, kfdia
      IF ((.NOT. laland(jl)) .AND. ktype(jl) == 0) THEN
        zdthdp = (ztheta(jl,jk)-ztheta(jl,jk-1))*zcloia/ &
                 (papm1(jl,jk)-papm1(jl,jk-1))
        lo = zdthdp < zdthmin(jl)
        zdthmin(jl) = MERGE(zdthdp,zdthmin(jl),lo)
        invb(jl) = MERGE(jk,invb(jl),lo)
      END IF
    END DO
  END DO

  DO jk = ktdia, klev
    DO jl = kidia, kfdia

!-- 1.4 t, q and qs provisional values at t+dt
!       effectives l for vaporisation and sublimation

      ztp1(jl,jk) = ptm1(jl,jk) + ztmst*ptte(jl,jk)
      zqp1(jl,jk) = pqm1(jl,jk) + ztmst*pqte(jl,jk)
      zqp1(jl,jk) = MAX(zqp1(jl,jk),zepqp1)

      zrcp = 1./(cpd+zcons1*zqp1(jl,jk))
      zlvdcp(jl,jk) = alv*zrcp
      zlsdcp(jl,jk) = als*zrcp

      zztc = tmelt - ztp1(jl,jk)
      lo = zztc<0.
      zcvm4 = MERGE(c4les,c4ies,lo)
      zcvm5 = MERGE(c5les*zlvdcp(jl,jk),c5ies*zlsdcp(jl,jk),lo)
      zzqt = 1./(ztp1(jl,jk)-zcvm4)
#ifndef NOLOOKUP
      it = ztp1(jl,jk)*1000.
      IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
#ifdef DEBUG
      IF (it<jptlucu1 .OR. it>jptlucu2) THEN
        WRITE(*,*) 'cond   : 1 it=',it,jl,jk
      END IF
#endif
      it = MAX(MIN(it,jptlucu2),jptlucu1)
      zqsp1(jl,jk) = tlucua(it)/papp1(jl,jk)
#else
      zqsp1(jl,jk) = c2es*EXP(MERGE(c3les,c3ies,lo)*(-zztc) &
                   / (ztp1(jl,jk)-zcvm4))/papp1(jl,jk)
#endif
      zqsp1(jl,jk) = MIN(zqsp1(jl,jk),0.5)
      zcor = 1./(1.-vtmpc1*zqsp1(jl,jk))
      zqsp1(jl,jk) = zqsp1(jl,jk)*zcor
      zqcon(jl,jk) = 1./(1.+zcvm5*zqsp1(jl,jk)*zcor*zzqt**2)

!-- 1.5 Threshold relative humidity

      jb = invb(jl)

      IF (lamip2) THEN
        lo = (jb==16 .OR. jb==17)
        IF (ktype(jl) == 0) THEN
          lo1 = (jk >= jb .AND. jk <= jb+1)
          zrhc(jl,jk) = zrtl+(zrs-zrtl)*EXP(1.-(paphp1(jl,klevp1)/ &
                        papp1(jl,jk))**nexl)
          zsat(jl,jk) = 1.
          IF (lo .AND. lo1) THEN
            zsat(jl,jk) = zsatsc*zrhc(jl,jk)
            zrhc(jl,jk) = zrhsc*zrhc(jl,jk)
          ENDIF
        ELSE
          zrhc(jl,jk) = zrtc + (zrs-zrtc)*EXP(1.-(paphp1(jl,klevp1)/ &
                        papp1(jl,jk))**nexc)
          zsat(jl,jk) = 1.
        END IF
      ELSE
        lo = (jb>klev-6 .AND. jb<klevm1)
        IF (klab(jl,jk)==2) THEN
          zrhc(jl,jk) = zrtc + (zrs-zrtc)*EXP(1.-(paphp1(jl,klevp1)/ &
                        papp1(jl,jk))**nexc)
          zsat(jl,jk) = 1.
        ELSE
          IF (lmidatm) THEN
            IF (jk>=jb+1 .AND. lo) THEN
              zrhc(jl,jk) = zrhsc
              zsat(jl,jk) = zsatsc
            ELSE
              zrhc(jl,jk) = zrtl + (zrs-zrtl)*EXP(1.-(paphp1(jl,klevp1)/ &
                            papp1(jl,jk))**nexl)
              zsat(jl,jk) = 1.
            END IF
          ELSE
            IF (jk<=jb+1 .AND. lo) THEN
              zrhc(jl,jk) = zrhsc
              zsat(jl,jk) = zsatsc
            ELSE
              zrhc(jl,jk) = zrtl + (zrs-zrtl)*EXP(1.-(paphp1(jl,klevp1)/ &
                            papp1(jl,jk))**nexl)
              zsat(jl,jk) = 1.
            END IF
          END IF
        END IF
      END IF

!-- 1.6 Cloud cover at t+dt

      ! --------------------- CHANGED FOR AMIP2 ----------------------------------
      !
      IF (lamip2) THEN
        !
        ! 1.6a  cloud cover at t+dt based on t-dt variables
        !       -----> Lenderink et al, 1998
        !
#ifndef NOLOOKUP
        it = ptm1(jl,jk)*1000.
#ifdef DEBUG
        IF (it<jptlucu1 .OR. it>jptlucu2) THEN
          WRITE(*,*) 'condgl : 1 it=',it,jl,jk
        END IF
#endif
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqsm1 = tlucua(it)/papm1(jl,jk)
#else
        WRITE(*,*) 'This is a special version of ECHAM4 (AMIP2)'
        CALL finish ('cond:','-DNOLOOKUP must not defined at compile time')
#endif
        zqsm1 = MIN(zqsm1,0.5)
        zcor  = 1./(1.-vtmpc1*zqsm1)
        zqsm1 = zqsm1*zcor
        zqr   = zqsm1*zsat(jl,jk)*zrhc(jl,jk)
        zclcp1(jl,jk) = (pqm1(jl,jk)-zqr)/(zqsm1*zsat(jl,jk)-zqr)
        zclcp1(jl,jk) = MAX(zclcp1(jl,jk),0.)
        zclcp1(jl,jk) = MIN(zclcp1(jl,jk),1.)
        zclcp1(jl,jk) = 1. - SQRT(1.-zclcp1(jl,jk))
      ELSE
        ! 
        ! 1.6b Cloud cover at t+dt (standard scheme)
        !
        zqr = zqsp1(jl,jk)*zsat(jl,jk)*zrhc(jl,jk)
        zclcp1(jl,jk) = (zqp1(jl,jk)-zqr)/(zqsp1(jl,jk)*zsat(jl,jk)-zqr)
        zclcp1(jl,jk) = MAX(zclcp1(jl,jk),0.)
        zclcp1(jl,jk) = MIN(zclcp1(jl,jk),1.)
        zclcp1(jl,jk) = 1. - SQRT(1.-zclcp1(jl,jk))
      END IF
    END DO
  END DO
#ifndef NOLOOKUP
  IF (lookupoverflow) CALL lookuperror ('cond (1)    ')
#endif

!-- 1.8 Cloud tops/bases and thickness

  ! 1.8.1 Geopotential at half levels

  DO jk = ktdia, klev - 1
    DO jl = kidia, kfdia
      zgeoh(jl,jk+1) = (pgeom1(jl,jk+1)+pgeom1(jl,jk))*0.5
    END DO
  END DO

  DO jl = kidia, kfdia
    zgeoh(jl,1) = 2.*pgeom1(jl,1) - zgeoh(jl,2)
    zgeoh(jl,klevp1) = 0.
  END DO

  ! 1.8.2  Cloud tops

  DO jl = kidia, kfdia
    IF (zclcp1(jl,1)>=zepclc) THEN
      ztop(jl,1) = zgeoh(jl,1)
    ELSE
      ztop(jl,1) = -1.
    END IF
  END DO

  DO jk = ktdia + 1, klev
    DO jl = kidia, kfdia
      IF (zclcp1(jl,jk)>=zepclc .AND. ztop(jl,jk-1)<0.) THEN
        ztop(jl,jk) = zgeoh(jl,jk)
      ELSE IF (zclcp1(jl,jk)>=zepclc) THEN
        ztop(jl,jk) = ztop(jl,jk-1)
      ELSE
        ztop(jl,jk) = -1.
      END IF
    END DO
  END DO

  ! 1.8.3 Cloud bases

  DO jl = kidia, kfdia
    IF (zclcp1(jl,klev)>=zepclc) THEN
      zbas(jl,klev) = zgeoh(jl,klevp1)
    ELSE
      zbas(jl,klev) = -1.
    END IF
  END DO

  DO jk = klev - 1, ktdia, -1
    DO jl = kidia, kfdia
      IF (zclcp1(jl,jk)>=zepclc .AND. zbas(jl,jk+1)<0.) THEN
        zbas(jl,jk) = zgeoh(jl,jk+1)
      ELSE IF (zclcp1(jl,jk)>=zepclc) THEN
        zbas(jl,jk) = zbas(jl,jk+1)
      ELSE
        zbas(jl,jk) = -1.
      END IF
    END DO
  END DO

  DO jk = ktdia, klev

!-- 2. Snow melt and saturation specific humidities at t-dt

    IF (lamip2) THEN
!DIR$ IVDEP
!OCL NOVREC,NOFMADD
      DO jl = kidia, kfdia

!-- 2.1 Melting of incoming snow

        zdp(jl) = paphp1(jl,jk+1) - paphp1(jl,jk)
        zcons = zcons2*zdp(jl)/(zlsdcp(jl,jk)-zlvdcp(jl,jk))
        zsnmlt = MIN(zsfl(jl),zcons*MAX(0.,(ztp1(jl,jk)-zmeltp2)))
        zsnmt(jl) = zsnmt(jl) + zsnmlt
        zrfln(jl) = zrfl(jl) + zsnmlt
        zsfln(jl) = zsfl(jl) - zsnmlt
        zdt = -zsnmlt/zcons
        ztp1(jl,jk) = ztp1(jl,jk) + zdt
        ptte(jl,jk) = ptte(jl,jk) + zdt*zqtmst

!-- 2.3 Saturation specific humidity for t-dt

#ifndef NOLOOKUP
        it = ptm1(jl,jk)*1000.
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
#ifdef DEBUG
        IF (it<jptlucu1 .OR. it>jptlucu2) THEN
          WRITE(*,*) 'cond   : 2 it=',it,jl,jk
        END IF
#endif
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqsm1 = tlucua(it)/papm1(jl,jk)
#else
        lo = ptm1(jl,jk) > tmelt
        zqsm1 = c2es*EXP(MERGE(c3les,c3ies,lo)*(ptm1(jl,jk)-tmelt) &
              / (ptm1(jl,jk)-MERGE(c4les,c4ies,lo)))/papm1(jl,jk)
#endif
        zqsm1 = MIN(zqsm1,0.5)
        zcor  = 1./(1.-vtmpc1*zqsm1)
        zqsm1 = zqsm1*zcor

        !------------------ CHANGES FOR AMIP2 --------------------------------
        ! ------------> Lenderink et al., 1998 (KNMI)
        !
        !zl =  l_{v/s} / c_p
        !zlvcdqsdt = l dq_sat / c_p dt
        !zdqsdt = dq_sat / dt
        !
        zztc = tmelt-ptm1(jl,jk)
        IF (zztc.GE.0.) THEN
          zl    = als / (cpd+zcons1*pqm1(jl,jk))
          zcvm4 = c4ies
          zcvm5 = c5ies*zl
        ELSE
          zl    = alv / (cpd+zcons1*pqm1(jl,jk))
          zcvm4 = c4les
          zcvm5 = c5les*zl
        END IF
        zzqt      = 1./(ptm1(jl,jk)-zcvm4)
        zqconm1   = 1./(1.+zcvm5*zqsm1*zcor*zzqt**2)
        zlvcdqsdt = (1./zqconm1-1.)
        zdqsdt    = zlvcdqsdt / zl
        !-------------------------------------------------------------

!-- 3. Provisional values of cloud cover, humidity and cloud water.

!-- 3.1 Cloud cover at t-dt

        zqr = zqsm1*zsat(jl,jk)*zrhc(jl,jk)
        zclcm1 = (pqm1(jl,jk)-zqr)/(zqsm1*zsat(jl,jk)-zqr)
        zclcm1 = MAX(zclcm1,0.)
        zclcm1 = MIN(zclcm1,1.)
        zclcm1 = 1. - SQRT(1.-zclcm1)

!-- 3.2 Humidity and liquid water in the cloudy and cloud-free
!       part for t-1.

        IF (zclcm1 > zepclc) THEN
          zqpa = MAX(zqsm1,pqm1(jl,jk))
        ELSE
          zqpa = pqm1(jl,jk)
        END IF

        IF (zclcp1(jl,jk) > zepclc .AND. pxm1(jl,jk) > 0.) THEN
          zwpa = pxm1(jl,jk)/zclcp1(jl,jk)
          zwma = 0.
        ELSE
          zwpa = pxm1(jl,jk)
          zwma = zwpa
        END IF


!-- 4. Condensation/evaporation rate

        !  CHANGES FOR AMIP2 

        !
        ! ----> Lenderink et al., 1998
        !
        zclcaux = zclcp1(jl,jk)
        zwpa    = pxm1(jl,jk)/MAX(zclcaux,zepclc)
        zqldt   = pxte(jl,jk)*ztmst
        !------ ADDED FOR AMIP2 (E. Roeckner, 29.10.98) -------
        IF(zclcaux <= zanvcc) THEN
          zqldt = zqldt+pxtec(jl,jk)*ztmst
        END IF
        IF (zqldt > 0) THEN
          zqldtstar = zqldt
          zwpa      = zwpa+zqldt
        ELSE
          zqldtstar = 0.
          zwpa      = zwpa+zqldt/MAX(zclcaux,zepclc)
        ENDIF
        zwpa      = MAX(zwpa, 0.)  ! ADDED FOR AMIP2 (E. Roeckner, 29.10.98)

        zqvdt = pqte(jl,jk)*ztmst
        zdtdt = ptte(jl,jk)*ztmst

        zdtdtstar = zdtdt+zl*(zclcaux*zqvdt-(1-zclcaux)*zqldtstar)
        zdqsat    = zdtdtstar*zdqsdt/(1.+zclcaux*zlvcdqsdt)

        zevape      = zqldtstar
        zcondc      = MAX(-zwpa,zqvdt-zdqsat)
        zconrat(jl) = zclcaux*zcondc+(1.-zclcaux)*(-zevape)
        zwpvn       = zwpa+zcondc

!-- 5. Cloud physics and precipitation fluxes

!-- 5.1 Autoconversion and accretion of cloud droplets in warm
!       clouds and sedimentation of ice crystals in cold clouds

        ! ---------------------- CHANGES FOR AMIP2 -----------------------
        ! --------------> Lenderink et al., 1998
        !
 
        ztp1tmp = ztp1(jl,jk)+zl*zconrat(jl)
        zqp1tmp = zqp1(jl,jk)-zconrat(jl)
        zzqt = 1./(ztp1tmp-zcvm4)
#ifndef NOLOOKUP
        it = ztp1tmp*1000
#ifdef DEBUG
        IF (it<jptlucu1 .OR. it>jptlucu2) THEN
          WRITE(*,*) 'condgl : 2 it=',it,jl,jk
        END IF
#endif
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqsp1tmp = tlucua(it)/papp1(jl,jk)
#else
        CALL finish ('cond:','-DNOLOOKUP must not defined at compile time')
#endif
        zqsp1tmp = MIN(zqsp1tmp,0.5)
        zcor = 1./(1.-vtmpc1*zqsp1tmp)
        zqsp1tmp = zqsp1tmp*zcor

        zoversat = 1.e-3*zqsp1tmp       ! added for numerical reasons
        IF (zqp1tmp > zqsp1tmp+zoversat) THEN
          zdconratcor = zqp1tmp-(zqsp1tmp+zoversat)
          zdconratcor = zdconratcor*zqcon(jl,jk)
          zconrat(jl) = zconrat(jl)+zdconratcor
        END IF

        IF(zconrat(jl) > 0.) THEN
          zconrat(jl) = MIN(zconrat(jl),zqp1tmp)
        END IF

        ztp1(jl,jk) = ztp1(jl,jk) + zl*zconrat(jl)
        zqp1(jl,jk) = zqp1(jl,jk) - zconrat(jl)

        zrcp = 1./(cpd+zcons1*zqp1(jl,jk))
        zlvdcp(jl,jk) = alv*zrcp
        zlsdcp(jl,jk) = als*zrcp

        ztc = ztp1(jl,jk)-tmelt
        IF (ztc > 0.) THEN
          zl=zlvdcp(jl,jk)
        ELSE
          zl=zlsdcp(jl,jk)
        END IF
 
        zztc = tmelt-ztp1(jl,jk)
        IF (zztc >= 0.) THEN
          zcvm4 = c4ies
          zcvm5 = c5ies*zlsdcp(jl,jk)
        ELSE
          zcvm4 = c4les
          zcvm5 = c5les*zlvdcp(jl,jk)
        END IF
        zzqt = 1./(ztp1(jl,jk)-zcvm4)
#ifndef NOLOOKUP
        it = ztp1(jl,jk)*1000.
#ifdef DEBUG
        IF (it<jptlucu1 .OR. it>jptlucu2) THEN
          WRITE(*,*) 'condgl : 3 it=',it,jl,jk
        END IF
#endif
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqsp1(jl,jk) = tlucua(it)/papp1(jl,jk)
#else
        CALL finish ('cond:','-DNOLOOKUP must not defined at compile time')
#endif
        zqsp1(jl,jk) = MIN(zqsp1(jl,jk),0.5)
        zcor = 1./(1.-vtmpc1*zqsp1(jl,jk))
        zqsp1(jl,jk) = zqsp1(jl,jk)*zcor
        zqcon(jl,jk) = 1./(1.+zcvm5*zqsp1(jl,jk)*zcor*zzqt**2)

        IF(zclcp1(jl,jk) > zanvcc) zwpvn=zwpvn+pxtec(jl,jk)*ztmst/zclcp1(jl,jk)

        IF (zwpvn*zclcp1(jl,jk)>zepzwp) THEN
          zwrc = MERGE(zwrl,zwrs,laland(jl))
          lo = ztop(jl,jk) - zbas(jl,jk) < zdnopr
          zwrc = MERGE(zwrsc,zwrc,lo)
          ztc = ztp1(jl,jk) - tmelt
          zrho = papp1(jl,jk)/(rd*ztp1(jl,jk))
          zzfac = MERGE(1.,0.,ztc<0.)
          zalpha = (1.-zzfac) + zzfac*(zcaa+(1.-zcaa)*EXP(-zcab*ztc**2))
          zwpvni = (1.-zalpha)*zwpvn
          zaci = zconsi*(zwpvni*zrho)**1.16
          zwpvnl = zalpha*zwpvn
          zqq = (zwpvnl/zwrc)**2
          zqq = MAX(zqq,1.E-6)
          zqq = MIN(zqq,100.)
          zacl = zcnull*zwpvnl*(1.-EXP(-zqq))
          zpcc = zwpvnl*(zrfln(jl)+zsfln(jl))*zc1
          zpcc = MAX(zpcc,0.)
          zdacl = MIN((zacl+zpcc),zwpvnl)

          zdaci = MIN(zaci,zwpvni)

          zdac(jl) = (zdacl+zdaci)*zclcp1(jl,jk)
        ELSE
          zdac(jl) = zwpvn*zclcp1(jl,jk)
        END IF

        zep(jl) = 0.

!-- 5.2 New precipitation fluxes after consideration of
!       coalescence and sedimentation processes.

        zzdr = MAX(0.,zcons2*zdp(jl)*zdac(jl))
        IF (ztp1(jl,jk) > tmelt) THEN
          zrfln(jl) = zrfln(jl) + zzdr
          zrain(jl) = zrain(jl) + zzdr
        ELSE
          zsfln(jl) = zsfln(jl) + zzdr
          zsnow(jl) = zsnow(jl) + zzdr
        END IF
        zfln(jl) = zrfln(jl) + zsfln(jl)
      END DO
    ELSE
!DIR$ IVDEP
!OCL NOVREC,NOFMADD
      DO jl = kidia, kfdia

!-- 2.1 Melting of incoming snow

        zdp(jl) = paphp1(jl,jk+1) - paphp1(jl,jk)
        zcons = zcons2*zdp(jl)/(zlsdcp(jl,jk)-zlvdcp(jl,jk))
        zsnmlt = MIN(zsfl(jl),zcons*MAX(0.,(ztp1(jl,jk)-zmeltp2)))
        zsnmt(jl) = zsnmt(jl) + zsnmlt
        zrfln(jl) = zrfl(jl) + zsnmlt
        zsfln(jl) = zsfl(jl) - zsnmlt
        zdt = -zsnmlt/zcons
        ztp1(jl,jk) = ztp1(jl,jk) + zdt
        ptte(jl,jk) = ptte(jl,jk) + zdt*zqtmst

!-- 2.3 Saturation specific humidity for t-dt

#ifndef NOLOOKUP
        it = ptm1(jl,jk)*1000.
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
#ifdef DEBUG
        IF (it<jptlucu1 .OR. it>jptlucu2) THEN
          WRITE(*,*) 'cond   : 2 it=',it,jl,jk
        END IF
#endif
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqsm1 = tlucua(it)/papm1(jl,jk)
#else
        lo = ptm1(jl,jk) > tmelt
        zqsm1 = c2es*EXP(MERGE(c3les,c3ies,lo)*(ptm1(jl,jk)-tmelt) &
              / (ptm1(jl,jk)-MERGE(c4les,c4ies,lo)))/papm1(jl,jk)
#endif
        zqsm1 = MIN(zqsm1,0.5)
        zcor  = 1./(1.-vtmpc1*zqsm1)
        zqsm1 = zqsm1*zcor

!-- 3. Provisional values of cloud cover, humidity and cloud water.

!-- 3.1 Cloud cover at t-dt

        zqr = zqsm1*zsat(jl,jk)*zrhc(jl,jk)
        zclcm1 = (pqm1(jl,jk)-zqr)/(zqsm1*zsat(jl,jk)-zqr)
        zclcm1 = MAX(zclcm1,0.)
        zclcm1 = MIN(zclcm1,1.)
        zclcm1 = 1. - SQRT(1.-zclcm1)

!-- 3.2 Humidity and liquid water in the cloudy and cloud-free
!       part for t-1.

        IF (zclcm1 > 1.E-10) THEN
          zqpa = MAX(zqsm1,pqm1(jl,jk))
        ELSE
          zqpa = pqm1(jl,jk)
        END IF

        IF (zclcp1(jl,jk) > 1.E-10 .AND. pxm1(jl,jk) > 0.) THEN
          zwpa = pxm1(jl,jk)/zclcp1(jl,jk)
          zwma = 0.
        ELSE
          zwpa = pxm1(jl,jk)
          zwma = zwpa
        END IF


!-- 4. Condensation/evaporation rate

        zcpn    = (zqpa+pqte(jl,jk)*ztmst-zqsp1(jl,jk))*zqcon(jl,jk)
        zxetmst = pxte(jl,jk)*ztmst
        zwpvn = zwpa + zxetmst
        zwmvn = zwma + zxetmst
        lo    = zwpvn >= 0.
        lo1   = lo .AND. -zcpn > zwpvn
        zcpn  = MERGE(-zwpvn,zcpn,lo1)
        zwi   = -zwpvn
        zwmax = MAX(zcpn,0.)
        zwps  = zwpvn + zcpn
        zwpvn = MERGE(zwps,zwmax,lo)
        zcpn  = MERGE(zcpn,zwpvn+zwi,lo)
        zconrat(jl) = zclcp1(jl,jk)*zcpn + (1.-zclcp1(jl,jk))*(-zwmvn)

!-- 5. Cloud physics and precipitation fluxes

!-- 5.1 Autoconversion and accretion of cloud droplets in warm
!       clouds and sedimentation of ice crystals in cold clouds

        ztc = ztp1(jl,jk) - tmelt
        zl = MERGE(zlvdcp(jl,jk),zlsdcp(jl,jk),-ztc<0.)
        IF (zconrat(jl) > 0.) THEN
          zconrat(jl) = MIN(zconrat(jl),zqp1(jl,jk))
        END IF
        ztp1(jl,jk) = ztp1(jl,jk) + zl*zconrat(jl)
        zqp1(jl,jk) = zqp1(jl,jk) - zconrat(jl)
        zwpvn = zwpvn + pxtec(jl,jk)*ztmst
        IF (zwpvn*zclcp1(jl,jk)>zepzwp) THEN
          zwrc = MERGE(zwrl,zwrs,laland(jl))
          lo = ztop(jl,jk) - zbas(jl,jk) < zdnopr
          zwrc = MERGE(zwrsc,zwrc,lo)
          ztc = ztp1(jl,jk) - tmelt
          zrho = papp1(jl,jk)/(rd*ztp1(jl,jk))
          zzfac = MERGE(1.,0.,ztc<0.)
          zalpha = (1.-zzfac) + zzfac*(zcaa+(1.-zcaa)*EXP(-zcab*ztc**2))
          zwpvni = (1.-zalpha)*zwpvn
          zaci = zconsi*(zwpvni*zrho)**1.16
          zwpvnl = zalpha*zwpvn
          zqq = (zwpvnl/zwrc)**2
          zqq = MAX(zqq,1.E-6)
          zqq = MIN(zqq,100.)
          zacl = zcnull*zwpvnl*(1.-EXP(-zqq))
          zpcc = zwpvnl*(zrfln(jl)+zsfln(jl))*zc1
          zpcc = MAX(zpcc,0.)
          zdacl = MIN((zacl+zpcc),zwpvnl)
          zdaci = MIN(zwpvni,MAX(0.1*zwpvni,zaci))
          zdac(jl) = (zdacl+zdaci)*zclcp1(jl,jk)
        ELSE
          zdac(jl) = zwpvn*zclcp1(jl,jk)
        END IF

        zep(jl) = 0.

!-- 5.2 New precipitation fluxes after consideration of
!       coalescence and sedimentation processes.

        zzdr = MAX(0.,zcons2*zdp(jl)*zdac(jl))
        IF (ztp1(jl,jk) > tmelt) THEN
          zrfln(jl) = zrfln(jl) + zzdr
          zrain(jl) = zrain(jl) + zzdr
        ELSE
          zsfln(jl) = zsfln(jl) + zzdr
          zsnow(jl) = zsnow(jl) + zzdr
        END IF
        zfln(jl) = zrfln(jl) + zsfln(jl)
      END DO
    END IF

#ifndef NOLOOKUP
    IF (lookupoverflow) CALL lookuperror ('cond (2)    ')
#endif

!-- 5.3 At top layer or if there is no precipitation avoid evaporation.

    IF (ABS(zs)<prec .AND. jk>1) THEN
      zs = SUM(zfln(1:klon))
    END IF

    IF (jk>1 .AND. ABS(zs)>0.) THEN

!-- 6. Evaporation of precipitation

      DO jl = kidia, kfdia
        zclear = 1. - zclcp1(jl,jk)
        zzrh = zrhc(jl,jk)*zsat(jl,jk) + zclcp1(jl,jk)*(1.-zrhc(jl,jk)*zsat(jl,jk))
        zzep = MIN(0.,(MIN(zqp1(jl,jk),zqsp1(jl,jk)*zzrh)-zqsp1(jl,jk)))
        zzep = zzepcons*zzep
        zzep = MAX(-(zrfl(jl)+zsfl(jl)),zzep*zqcon(jl,jk)*zcons2*zdp(jl))
        zfln(jl) = zfln(jl) + zzep*zclear
        zep(jl) = zzep/(zcons2*zdp(jl))*zclear
        zfrac = zfln(jl)/MAX(zrfln(jl)+zsfln(jl),zepflm)
        zrflo = zrfln(jl)
        zsflo = zsfln(jl)
        zrfln(jl) = zrfln(jl)*zfrac
        zsfln(jl) = zsfln(jl)*zfrac
        zepr(jl) = zepr(jl) + (zrflo-zrfln(jl))
        zeps(jl) = zeps(jl) + (zsflo-zsfln(jl))

      END DO

!-- 7. Incrementation of t,q and x tendencies and fluxes' swap.

!-- 7.1 Resume computations for the top layer

    END IF

!-- 7.2 Modification of the t,q and x tendencies,
!       cloud cover for diagnostics and radiation.

    IF (lamip2) THEN
      DO jl = kidia, kfdia
        zdqdt = -(zconrat(jl)+zep(jl))*zqtmst
        zztc = ztp1(jl,jk) - tmelt
        lo = -zztc<0.
        zcvml = MERGE(zlvdcp(jl,jk),zlsdcp(jl,jk),lo)
        zdtdt = -zcvml*zdqdt
        zdxdt = (zconrat(jl)-zdac(jl))*zqtmst
        pxte(jl,jk) = pxte(jl,jk) + pxtec(jl,jk) + zdxdt

!-- AMIP2
        pqte(jl,jk) = pqte(jl,jk) + zdqdt
        ptte(jl,jk) = ptte(jl,jk) + zdtdt
        zqp1(jl,jk) = pqm1(jl,jk) + ztmst*pqte(jl,jk)
        ztp1(jl,jk) = ptm1(jl,jk) + ztmst*ptte(jl,jk)
#ifndef NOLOOKUP
        it = ztp1(jl,jk)*1000.
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
#ifdef DEBUG
        IF (it<jptlucu1 .OR. it>jptlucu2) THEN
          WRITE(*,*) 'cond   : 3 it=',it,jl,jk
        END IF
#endif
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqsp1(jl,jk) = tlucua(it)/papp1(jl,jk)
#else
        zqsp1(jl,jk) = c2es*EXP(MERGE(c3les,c3ies,lo)*(zztc) &
                     / (ztp1(jl,jk)-MERGE(c4les,c4ies,lo)))/papp1(jl,jk)
#endif
        zqsp1(jl,jk) = MIN(zqsp1(jl,jk),0.5)
        zcor = 1./(1.-vtmpc1*zqsp1(jl,jk))
        zqsp1(jl,jk) = zqsp1(jl,jk)*zcor

!-- AMIP2
        zoversat = zqsp1(jl,jk)*1.e-3
        IF (zqp1(jl,jk) > zqsp1(jl,jk)+zoversat) THEN
          zdxc = zqp1(jl,jk)-(zqsp1(jl,jk)+zoversat)
          zdxc = zdxc*zqcon(jl,jk)/ztmst
          zqp1(jl,jk) = zqp1(jl,jk)-zdxc*ztmst
        ELSE
          zdxc = 0.
        ENDIF

        zqr = zqsp1(jl,jk)*zsat(jl,jk)*zrhc(jl,jk)
        zclp1 = (zqp1(jl,jk)-zqr)/(zqsp1(jl,jk)*zsat(jl,jk)-zqr)
        zclp1 = MAX(zclp1,0.)
        zclp1 = MIN(zclp1,1.)
        zclp1 = 1. - SQRT(1.-zclp1)

!-- AMIP2
        zalwc(jl,jk) = pxm1(jl,jk) + ztmst*(pxte(jl,jk)+zdxc)

        zalwcold = zalwc(jl,jk)

!-- AMIP2
        lo = (zalwc(jl,jk) < zepclw .OR. zclp1 < zepclc)

        zalwc(jl,jk) = MERGE(0.,zalwc(jl,jk),lo)
        zdxcor = (zalwc(jl,jk)-zalwcold)/ztmst

!-- AMIP2
        pxte(jl,jk) = pxte(jl,jk) + zdxcor + zdxc
        pqte(jl,jk) = pqte(jl,jk) - zdxcor - zdxc
        ptte(jl,jk) = ptte(jl,jk) + (zdxcor+zdxc)*zcvml

        paclc(jl,jk)   = MERGE(0.,zclp1,lo)
        paclcac(jl,jk) = paclcacm(jl,jk) + paclc(jl,jk)*zdiagt
      END DO
    ELSE
      DO jl = kidia, kfdia
        zdqdt = -(zconrat(jl)+zep(jl))*zqtmst
        zztc = ztp1(jl,jk) - tmelt
        lo = -zztc<0.
        zcvml = MERGE(zlvdcp(jl,jk),zlsdcp(jl,jk),lo)
        zdtdt = -zcvml*zdqdt
        zdxdt = (zconrat(jl)-zdac(jl))*zqtmst
        pxte(jl,jk) = pxte(jl,jk) + pxtec(jl,jk) + zdxdt

#ifndef NOLOOKUP
        it = ztp1(jl,jk)*1000.
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
#ifdef DEBUG
        IF (it<jptlucu1 .OR. it>jptlucu2) THEN
          WRITE(*,*) 'cond   : 3 it=',it,jl,jk
        END IF
#endif
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqsp1(jl,jk) = tlucua(it)/papp1(jl,jk)
#else
        zqsp1(jl,jk) = c2es*EXP(MERGE(c3les,c3ies,lo)*(zztc) &
                     / (ztp1(jl,jk)-MERGE(c4les,c4ies,lo)))/papp1(jl,jk)
#endif
        zqsp1(jl,jk) = MIN(zqsp1(jl,jk),0.5)
        zcor = 1./(1.-vtmpc1*zqsp1(jl,jk))
        zqsp1(jl,jk) = zqsp1(jl,jk)*zcor

        zqr = zqsp1(jl,jk)*zsat(jl,jk)*zrhc(jl,jk)
        zclp1 = (zqp1(jl,jk)-zqr)/(zqsp1(jl,jk)*zsat(jl,jk)-zqr)
        zclp1 = MAX(zclp1,0.)
        zclp1 = MIN(zclp1,1.)
        zclp1 = 1. - SQRT(1.-zclp1)

        zalwc(jl,jk) = pxm1(jl,jk) + ztmst*pxte(jl,jk)
        zalwcold = zalwc(jl,jk)
        lo = (zalwc(jl,jk) < zepzwp .OR. zclp1 < zepclc)
        zalwc(jl,jk) = MERGE(0.,zalwc(jl,jk),lo)
        zdxcor = (zalwc(jl,jk)-zalwcold)/ztmst
        pxte(jl,jk) = pxte(jl,jk) + zdxcor
        pqte(jl,jk) = pqte(jl,jk) - zdxcor + zdqdt
        ptte(jl,jk) = ptte(jl,jk) + zcvml*zdxcor + zdtdt
        paclc(jl,jk)   = MERGE(0.,zclp1,lo)
        paclcac(jl,jk) = paclcacm(jl,jk) + paclc(jl,jk)*zdiagt
      END DO
    END IF
#ifndef NOLOOKUP
    IF (lookupoverflow) CALL lookuperror ('cond (3)    ')
#endif

!-- 7.3 Swap of the fluxes and end of the vertical loop

    DO jl = kidia, kfdia
      zrfl(jl) = zrfln(jl)
      zsfl(jl) = zsfln(jl)
    END DO
  END DO

!-- 7.4 Surface fluxes

  DO jl = kidia, kfdia
    prsfl(jl) = prsfl(jl) + zrfl(jl)
    pssfl(jl) = zsfl(jl)
    paprl(jl) = paprlm(jl) + zdiagw*(prsfl(jl)+pssfl(jl))
!fu++ SST diurnal cycle
    IF (ldsst) THEN
    prlarg_diurnal(jl,jrow) = prlarg_diurnal(jl,jrow)+ &
                zdiagw*(prsfl(jl)+pssfl(jl))  
    END IF

    paprs(jl) = paprs(jl) + zdiagw*pssfl(jl)
  END DO

!-- 7.5 Accumulated total cloudcover

  DO jl = kidia, kfdia
    paclcov(jl) = 1. - paclc(jl,1)
  END DO
  DO jk = 2, klev
    DO jl = kidia, kfdia
      paclcov(jl) = paclcov(jl)*(1.-MAX(paclc(jl,jk),paclc(jl,jk-1))) &
                  / (1.-MIN(paclc(jl,jk-1),1.-zepsec))
    END DO
  END DO
  DO jl = kidia, kfdia
    paclcov(jl) = 1. - paclcov(jl)
    paclcov(jl) = paclcovm(jl) + zdiagt*paclcov(jl)
  END DO

!-- 7.6 Vertically integrated humidity and cloud water

  DO jl = kidia, kfdia
    pqvi(jl) = 0.
    palwcvi(jl) = 0.
  END DO

  DO jk = ktdia, klev
    DO jl = kidia, kfdia
      pqvi(jl) = pqvi(jl) + pqm1(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl,jk))
      palwcvi(jl) = palwcvi(jl) + zalwc(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
    END DO
  END DO

  DO jl = kidia, kfdia
    pqvi(jl) = pqvim(jl) + zdiagt*pqvi(jl)/g
    palwcvi(jl) = palwcvim(jl) + zdiagt*palwcvi(jl)/g
  END DO

  zbud = budw(irow)
  dlsgrz(irow) = zdiagw*zbud*SUM(zrain(1:klon))
  dlsgsz(irow) = zdiagw*zbud*SUM(zsnow(1:klon))
  dlsmsz(irow) = zdiagw*zbud*SUM(zsnmt(1:klon))
  dlserz(irow) = zdiagw*zbud*SUM(zepr(1:klon))
  dlsesz(irow) = zdiagw*zbud*SUM(zeps(1:klon))

  RETURN
END SUBROUTINE cond

!END MODULE m_cond
