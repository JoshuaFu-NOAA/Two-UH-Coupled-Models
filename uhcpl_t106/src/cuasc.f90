!+ does the calculations for cloud ascents for cumulus parameterization

SUBROUTINE cuasc(klp2,k2lp2,klon,klev,klevp1,klevm1,ptenh,pqenh,pxenh,puen,                  &
                 pven,ktrac,pxtenh,pxten,pxtu,pmfuxt,pten,pqen,pqsen,pgeo,pgeoh,papp1,       &
                 paphp1,pverv,klwmin,ldcum,ldland,ktype,klab,ptu,pqu,plu,puu,pvu,pmfu,pmfub, &
                 pentr,pmfus,pmfuq,pmful,plude,pdmfup,khmin,phhatt,phcbase,pqsenh,kcbot,     &
                 kctop,kctop0,kcum)

  ! ---------------> ldland included for AMIP2 !!

  ! Description:
  !
  ! Does the calculations for cloud ascents for cumulus parameterization.
  !
  ! Method:
  !
  ! Produce cloud ascents for cu-parametrization
  ! (vertical profiles of t,q,l,u and v and corresponding
  ! fluxes as well as precipitation rates)
  !
  ! This routine is called from *cumastr*.
  !
  ! Lift surface air dry-adiabatically to cloud base and then calculate
  ! moist ascent for entraining/detraining plume.
  ! Entrainment and detrainment rates differ for shallow and deep cumulus
  ! convection.
  ! In case there is no penetrative or shallow convection check for
  ! possibility of mid level convection (cloud base values calculated in
  ! *cubasmc*)
  !
  ! Externals:
  !    *cuadjtq*  adjust t and q due to condensation in ascent
  !    *cuentr*   calculate entrainment/detrainment rates
  !    *cubasmc*  calculate cloud base values for midlevel convection
  !
  ! Authors:
  !
  ! M. Tiedtke, ECMWF, in July 1986 and December 1989, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !
  USE mo_start_dataset, ONLY: nstep,   & ! current time step
                              nstart     ! time step for start/restart
  USE mo_control,       ONLY: twodt,   & ! 2.*dtime
                              lamip2       
  USE mo_constants,     ONLY: g,       & ! gravity acceleration
                              cpd,     & ! specific heat at constant pressure
                              rd,      & ! gas constant for dry air
                              rv,      & !   idem      for water vapour
                              alv,     & ! latent heat for vaporisation
                              vtmpc1,  & ! vtmpc1=rv/rd-1
                              rcpd       ! rcpd=1./cpd
  USE mo_cumulus_flux,  ONLY: lmfdudv, & !
                              lmfmid,  & !
                              cmfcmin, & !
                              cprcon,  & !
                              cmfctop
  USE mo_machine,       ONLY: prec

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: k2lp2, klev, klevm1, klevp1, klon, klp2, ktrac

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: paphp1(klp2,klevp1), papp1(klp2,klev), pgeo(klp2,klev), &
                       pgeoh(klp2,klev), phcbase(klp2), phhatt(klp2,klev),     &
                       pqen(klp2,klev), pqenh(klp2,klev), pqsen(klp2,klev),    &
                       pqsenh(klp2,klev), pten(klp2,klev), ptenh(klp2,klev),   &
                       puen(k2lp2,klev), pven(k2lp2,klev), pverv(klp2,klev),   &
                       pxenh(klp2,klev), pxten(klon,klev,ktrac),               &
                       pxtenh(klon,klev,ktrac)
  INTEGER, INTENT (IN) :: khmin(klp2), klwmin(klp2)
  LOGICAL, INTENT (IN) :: ldland(klp2)  ! Included for AMIP2

  !  Scalar arguments with intent(InOut):
  INTEGER, INTENT (INOUT) :: kcum

  !  Array arguments with intent(InOut):
  REAL, INTENT (INOUT) :: pdmfup(klp2,klev), pentr(klp2), plu(klp2,klev),       &
                          plude(klp2,klev), pmfu(klp2,klev), pmfub(klp2),       &
                          pmful(klp2,klev), pmfuq(klp2,klev), pmfus(klp2,klev), &
                          pmfuxt(klon,klev,ktrac), pqu(klp2,klev),              &
                          ptu(klp2,klev), puu(klp2,klev), pvu(klp2,klev),       &
                          pxtu(klon,klev,ktrac)
  INTEGER, INTENT (INOUT) :: kcbot(klp2), kctop(klp2), kctop0(klp2), &
                             klab(klp2,klev), ktype(klp2)
  LOGICAL, INTENT (INOUT) :: ldcum(klp2)

  !  Local scalars: 
  REAL :: zbuo, zbuoyz, zcons2, zdmfdu, zdmfeu, zdnoprc, zdprho, zdrodz, zdt, &
          zdz, zfac, zga, zlnew, zmfmax, zmftest, zmfulk, zmfuqk, zmfusk,     &
          zmfuxtk, zmse, znevn, zodmax, zprcon, zqcod, zqeen, zqude, zscde,   &
          zscod, zseen, ztmst, zxteen, zxtude, zz, zzdmf
  INTEGER :: icall, ik, ikb, ikt, is, jk, jl, jt

  !  Local arrays: 
  REAL :: zbuoy(klp2), zdmfde(klp2), zdmfen(klp2), zmfuu(klp2), zmfuv(klp2), &
          zodetr(klp2,klev), zoentr(klp2,klev), zpbase(klp2), zph(klp2),     &
          zqold(klp2)
  LOGICAL :: loflag(klp2)

  !  External subroutines 
  EXTERNAL cuadjtq, cubasmc, cuentr

  !  Intrinsic functions 
  INTRINSIC LOG, MAX, MERGE, MIN


  !  Executable statements 

!-- 1. Specify parameters

  ztmst = twodt
  IF (nstep==nstart) ztmst = 0.5*twodt
  zcons2 = 1./(g*ztmst)

!-- 2. Set default values

  DO jl = 1, klon
    zmfuu(jl) = 0.
    zmfuv(jl) = 0.
    IF ( .NOT. ldcum(jl)) ktype(jl) = 0
  END DO
  DO jk = 1, klev
    DO jl = 1, klon
      plu(jl,jk) = 0.
      pmfu(jl,jk) = 0.
      pmfus(jl,jk) = 0.
      pmfuq(jl,jk) = 0.
      pmful(jl,jk) = 0.
      plude(jl,jk) = 0.
      pdmfup(jl,jk) = 0.
      IF ( .NOT. ldcum(jl) .OR. ktype(jl)==3) klab(jl,jk) = 0
      IF ( .NOT. ldcum(jl) .AND. paphp1(jl,jk)<4.E4) kctop0(jl) = jk
    END DO
    DO jt = 1, ktrac
      DO jl = 1, klon
        pmfuxt(jl,jk,jt) = 0.
      END DO

    END DO
  END DO

  IF (.NOT. lamip2) zdnoprc = 1.5E4

  DO jk = 1, klev
    DO jl = 1, klon
      zoentr(jl,jk) = 0.
      zodetr(jl,jk) = 0.
    END DO
  END DO

!-- 3. Initialize values at lifting level

  DO jl = 1, klon
    kctop(jl) = klevm1
    IF ( .NOT. ldcum(jl)) THEN
      kcbot(jl) = klevm1
      pmfub(jl) = 0.
      pqu(jl,klev) = 0.
    END IF
    pmfu(jl,klev) = pmfub(jl)
    pmfus(jl,klev) = pmfub(jl)*(cpd*ptu(jl,klev)+pgeoh(jl,klev))
    pmfuq(jl,klev) = pmfub(jl)*pqu(jl,klev)
    IF (lmfdudv) THEN
      zmfuu(jl) = pmfub(jl)*puu(jl,klev)
      zmfuv(jl) = pmfub(jl)*pvu(jl,klev)
    END IF
  END DO

  DO jt = 1, ktrac
    DO jl = 1, klon
      IF ( .NOT. ldcum(jl)) THEN
        pxtu(jl,klev,jt) = 0.
      END IF
      pmfuxt(jl,klev,jt) = pmfub(jl)*pxtu(jl,klev,jt)
    END DO
  END DO

  DO jl = 1, klon
    ldcum(jl) = .FALSE.
  END DO

!-- 3.1 Find organized entrainment at cloud base

  DO jl = 1, klon
    IF (ktype(jl)==1) THEN
      ikb = kcbot(jl)
      zbuoy(jl) = g*(ptu(jl,ikb)-ptenh(jl,ikb))/ptenh(jl,ikb) + &
                  g*0.608*(pqu(jl,ikb)-pqenh(jl,ikb))
      IF (zbuoy(jl)>0.) THEN
        zdz = (pgeo(jl,ikb-1)-pgeo(jl,ikb))/g
        zdrodz = -LOG(pten(jl,ikb-1)/pten(jl,ikb))/zdz - g/(rd*ptenh(jl,ikb))
        ! nb zoentr is here a fractional value
        zoentr(jl,ikb-1) = zbuoy(jl)*0.5/(1.+zbuoy(jl)*zdz) + zdrodz
        IF (lamip2) THEN
          zoentr(jl,ikb-1) = MIN(zoentr(jl,ikb-1),3.E-4)
        ELSE
          zoentr(jl,ikb-1) = MIN(zoentr(jl,ikb-1),1.E-3)
        END IF
        zoentr(jl,ikb-1) = MAX(zoentr(jl,ikb-1),0.)
      END IF
    END IF
  END DO

!-- 4.  Do ascent: Subcloud layer (klab=1) ,clouds (klab=2)
!       by doing first dry-adiabatic ascent and then
!       by adjusting t,q and l accordingly in *cuadjtq*,
!       then check for buoyancy and set flags accordingly

  DO jk = klevm1, 2, -1

    ! Specify cloud base values for midlevel convection
    ! in *cubasmc* in case there is not already convection

    ik = jk
    IF (lmfmid .AND. ik<klevm1 .AND. ik>klev-10) THEN
      CALL cubasmc(klp2,k2lp2,klon,klev,ik,pten,pqen,pqsen,puen,pven,ktrac,        &
                   pxten,pxtu,pmfuxt,pverv,pgeo,pgeoh,ldcum,ktype,klab,pmfu,pmfub, &
                   pentr,kcbot,ptu,pqu,plu,puu,pvu,pmfus,pmfuq,pmful,pdmfup,zmfuu, &
                   zmfuv)
    END IF

    is = 0
    DO jl = 1, klon
      is = is + klab(jl,jk+1)
      IF (klab(jl,jk+1) == 0) klab(jl,jk) = 0
      loflag(jl) = klab(jl,jk+1) > 0
      zph(jl) = paphp1(jl,jk)
      IF (ktype(jl)==3 .AND. jk==kcbot(jl)) THEN
        zmfmax = (paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
        IF (pmfub(jl)>zmfmax) THEN
!fu++     zfac = pmfub(jl)/zmfmax
          zfac = zmfmax/pmfub(jl)    !original Tiedtke
          pmfu(jl,jk+1) = pmfu(jl,jk+1)*zfac
          pmfus(jl,jk+1) = pmfus(jl,jk+1)*zfac
          pmfuq(jl,jk+1) = pmfuq(jl,jk+1)*zfac
          zmfuu(jl) = zmfuu(jl)*zfac
          zmfuv(jl) = zmfuv(jl)*zfac
        END IF
      END IF
    END DO
    DO jt = 1, ktrac
      DO jl = 1, klon
        IF (ktype(jl)==3 .AND. jk==kcbot(jl)) THEN
          zmfmax = (paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
          IF (pmfub(jl)>zmfmax) THEN
!fu++       zfac = pmfub(jl)/zmfmax
            zfac = zmfmax/pmfub(jl)    !original Tiedtke
            pmfuxt(jl,jk+1,jt) = pmfuxt(jl,jk+1,jt)*zfac
          END IF
        END IF
      END DO
    END DO

    ! Reset pmfub if necessary

    DO jl = 1, klon
      IF (ktype(jl)==3 .AND. jk==kcbot(jl)) THEN
        zmfmax = (paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
        pmfub(jl) = MIN(pmfub(jl),zmfmax)
      END IF
    END DO

    IF (is==0) CYCLE

!   Specify turbulent entrainment and detrainments
!   rates plus organized detrainment rates in *cuentr*

    ik = jk
    CALL cuentr(klp2,klon,klev,klevp1,ik,ptenh,paphp1,papp1,klwmin,ldcum,       &
                ktype,kcbot,kctop0,zpbase,pmfu,pentr,zodetr,khmin,pgeoh,zdmfen, &
                zdmfde)

    ! Do adiabatic ascent for entraining/detraining plume
    ! the cloud ensemble entrains environmental values
    ! in turbulent detrainment cloud ensemble values are detrained
    ! in organized detrainment the dry static energy and
    ! moisture that are neutral compared to the
    ! environmental air are detrained

    DO jl = 1, klon
      IF (loflag(jl)) THEN
        IF (jk<kcbot(jl)) THEN
          zmftest = pmfu(jl,jk+1) + zdmfen(jl) - zdmfde(jl)
          zmfmax = MIN(zmftest,(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2)
          zdmfen(jl) = MAX(zdmfen(jl)-MAX(zmftest-zmfmax,0.),0.)
        END IF
        zdmfde(jl) = MIN(zdmfde(jl),0.75*pmfu(jl,jk+1))
        pmfu(jl,jk) = pmfu(jl,jk+1) + zdmfen(jl) - zdmfde(jl)
        IF (jk<kcbot(jl)) THEN
          zdprho = (pgeoh(jl,jk)-pgeoh(jl,jk+1))/g
          zoentr(jl,jk) = zoentr(jl,jk)*zdprho*pmfu(jl,jk+1)
          zmftest = pmfu(jl,jk) + zoentr(jl,jk) - zodetr(jl,jk)
          zmfmax = MIN(zmftest,(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2)
          zoentr(jl,jk) = MAX(zoentr(jl,jk)-MAX(zmftest-zmfmax,0.),0.)
        END IF
        IF (ktype(jl)==1 .AND. jk<kcbot(jl) .AND. jk<=khmin(jl)) THEN
          ! limit organized detrainment to not allowing for too
          ! deep clouds
          zmse = cpd*ptu(jl,jk+1) + alv*pqu(jl,jk+1) + pgeoh(jl,jk+1)
          ikt = kctop0(jl)
          znevn = (pgeoh(jl,ikt)-pgeoh(jl,jk+1))*(zmse-phhatt(jl,jk+1))/g
          IF (znevn<=0.) znevn = 1.
          zdprho = (pgeoh(jl,jk)-pgeoh(jl,jk+1))/g
          zodmax = ((phcbase(jl)-zmse)/znevn)*zdprho*pmfu(jl,jk+1)
          zodmax = MAX(zodmax,0.)
          zodetr(jl,jk) = MIN(zodetr(jl,jk),zodmax)
        END IF
        zodetr(jl,jk) = MIN(zodetr(jl,jk),0.75*pmfu(jl,jk))
        pmfu(jl,jk) = pmfu(jl,jk) + zoentr(jl,jk) - zodetr(jl,jk)
        zqeen = pqenh(jl,jk+1)*zdmfen(jl)
        zqeen = zqeen + pqenh(jl,jk+1)*zoentr(jl,jk)
        zseen = (cpd*ptenh(jl,jk+1)+pgeoh(jl,jk+1))*zdmfen(jl)
        zseen = zseen + (cpd*ptenh(jl,jk+1)+pgeoh(jl,jk+1))*zoentr(jl,jk)
        zscde = (cpd*ptu(jl,jk+1)+pgeoh(jl,jk+1))*zdmfde(jl)
        ! find moist static energy that give nonbuoyant air
        zga = alv*pqsenh(jl,jk+1)/(rv*(ptenh(jl,jk+1)**2))
        zdt = (plu(jl,jk+1)-0.608*(pqsenh(jl,jk+1)-pqenh(jl,jk+1))) &
            /(1./ptenh(jl,jk+1)+0.608*zga)
        zscod = cpd*ptenh(jl,jk+1) + pgeoh(jl,jk+1) + cpd*zdt
        zscde = zscde + zodetr(jl,jk)*zscod
        zqude = pqu(jl,jk+1)*zdmfde(jl)
        zqcod = pqsenh(jl,jk+1) + zga*zdt
        zqude = zqude + zodetr(jl,jk)*zqcod
        plude(jl,jk) = plu(jl,jk+1)*zdmfde(jl)
        plude(jl,jk) = plude(jl,jk) + plu(jl,jk+1)*zodetr(jl,jk)
        zmfusk = pmfus(jl,jk+1) + zseen - zscde
        zmfuqk = pmfuq(jl,jk+1) + zqeen - zqude
        zmfulk = pmful(jl,jk+1) - plude(jl,jk)
        plu(jl,jk) = zmfulk*(1./MAX(cmfcmin,pmfu(jl,jk)))
        pqu(jl,jk) = zmfuqk*(1./MAX(cmfcmin,pmfu(jl,jk)))
        ptu(jl,jk) = (zmfusk*(1./MAX(cmfcmin,pmfu(jl,jk)))-pgeoh(jl,jk))*rcpd
        ptu(jl,jk) = MAX(100.,ptu(jl,jk))
        ptu(jl,jk) = MIN(400.,ptu(jl,jk))
        zqold(jl) = pqu(jl,jk)
      END IF
    END DO

    DO jt = 1, ktrac
      DO jl = 1, klon
        IF (loflag(jl)) THEN
          zxteen = pxtenh(jl,jk+1,jt)*(zdmfen(jl)+zoentr(jl,jk))
          zxtude = pxtu(jl,jk+1,jt)*(zdmfde(jl)+zodetr(jl,jk))
          zmfuxtk = pmfuxt(jl,jk+1,jt) + zxteen - zxtude
          pxtu(jl,jk,jt) = zmfuxtk*(1./MAX(cmfcmin,pmfu(jl,jk)))
        END IF
      END DO
    END DO

    ! Do corrections for moist ascent
    ! by adjusting t,q and l in *cuadjtq*

    ik = jk
    icall = 1
    CALL cuadjtq(klp2,klon,klev,ik,zph,ptu,pqu,loflag,icall)

    IF (lamip2) THEN
!DIR$ IVDEP
!OCL NOVREC
      DO jl = 1, klon
        IF (loflag(jl)) THEN
          IF (ABS(pqu(jl,jk)-zqold(jl))>0.) THEN
            klab(jl,jk) = 2
            plu(jl,jk) = plu(jl,jk) + zqold(jl) - pqu(jl,jk)
            zbuo = ptu(jl,jk)*(1.+vtmpc1*pqu(jl,jk)-plu(jl,jk)) - &
                   ptenh(jl,jk)*(1.+vtmpc1*pqenh(jl,jk)-pxenh(jl,jk))
            IF (klab(jl,jk+1)==1) zbuo = zbuo + 0.5
            IF (zbuo>0. .AND. pmfu(jl,jk)>=0.01*pmfub(jl) &
                        .AND. jk>=kctop0(jl)) THEN
              kctop(jl) = jk
              ldcum(jl) = .TRUE.

!-- AMIP2
              zdnoprc = MERGE(3.E4,1.5E4,ldland(jl))  ! <--- Hier eine Ursache fuer 
                                                      !      den deutlich besseren
                                                      !      Monsun !!!!!

              zprcon = MERGE(0.,cprcon,zpbase(jl)-paphp1(jl,jk)<zdnoprc)
              zlnew = plu(jl,jk)/(1.+zprcon*(pgeoh(jl,jk)-pgeoh(jl,jk+1)))
              pdmfup(jl,jk) = MAX(0.,(plu(jl,jk)-zlnew)*pmfu(jl,jk))
              plu(jl,jk) = zlnew
            ELSE
              klab(jl,jk) = 0
              pmfu(jl,jk) = 0.
            END IF
          END IF
        END IF
      END DO
    ELSE
!DIR$ IVDEP
!OCL NOVREC
      DO jl = 1, klon
        IF (loflag(jl)) THEN
          IF (ABS(pqu(jl,jk)-zqold(jl))>0.) THEN
            klab(jl,jk) = 2
            plu(jl,jk) = plu(jl,jk) + zqold(jl) - pqu(jl,jk)
            zbuo = ptu(jl,jk)*(1.+vtmpc1*pqu(jl,jk)-plu(jl,jk)) - &
                   ptenh(jl,jk)*(1.+vtmpc1*pqenh(jl,jk)-pxenh(jl,jk))
            IF (klab(jl,jk+1) == 1) zbuo = zbuo + 0.5
            IF (zbuo>0. .AND. pmfu(jl,jk)>=0.01*pmfub(jl) &
                        .AND. jk>=kctop0(jl)) THEN
              kctop(jl) = jk
              ldcum(jl) = .TRUE.

              zprcon = MERGE(0.,cprcon,zpbase(jl)-paphp1(jl,jk)<zdnoprc)
              zlnew  = plu(jl,jk)/(1.+zprcon*(pgeoh(jl,jk)-pgeoh(jl,jk+1)))
              pdmfup(jl,jk) = MAX(0.,(plu(jl,jk)-zlnew)*pmfu(jl,jk))
              plu(jl,jk) = zlnew
            ELSE
              klab(jl,jk) = 0
              pmfu(jl,jk) = 0.
            END IF
          END IF
        END IF
      END DO
    END IF

    DO jl = 1, klon
      IF (loflag(jl)) THEN
        pmful(jl,jk) = plu(jl,jk)*pmfu(jl,jk)
        pmfus(jl,jk) = (cpd*ptu(jl,jk)+pgeoh(jl,jk))*pmfu(jl,jk)
        pmfuq(jl,jk) = pqu(jl,jk)*pmfu(jl,jk)
      END IF
    END DO
    DO jt = 1, ktrac
      DO jl = 1, klon
        IF (loflag(jl)) THEN
          pmfuxt(jl,jk,jt) = pxtu(jl,jk,jt)*pmfu(jl,jk)
        END IF
      END DO
    END DO

    IF (lmfdudv) THEN
      DO jl = 1, klon
        zdmfen(jl) = zdmfen(jl) + zoentr(jl,jk)
        zdmfde(jl) = zdmfde(jl) + zodetr(jl,jk)
      END DO
      DO jl = 1, klon
        IF (loflag(jl)) THEN
          IF (ktype(jl)==1 .OR. ktype(jl)==3) THEN
            zz = MERGE(3.,2.,ABS(zdmfen(jl))<prec)
          ELSE
            zz = MERGE(1.,0.,ABS(zdmfen(jl))<prec)
          END IF
          zdmfeu = zdmfen(jl) + zz*zdmfde(jl)
          zdmfdu = zdmfde(jl) + zz*zdmfde(jl)
          zdmfdu = MIN(zdmfdu,0.75*pmfu(jl,jk+1))
          zmfuu(jl) = zmfuu(jl) + zdmfeu*puen(jl,jk) - zdmfdu*puu(jl,jk+1)
          zmfuv(jl) = zmfuv(jl) + zdmfeu*pven(jl,jk) - zdmfdu*pvu(jl,jk+1)
          IF (pmfu(jl,jk)>0.) THEN
            puu(jl,jk) = zmfuu(jl)*(1./pmfu(jl,jk))
            pvu(jl,jk) = zmfuv(jl)*(1./pmfu(jl,jk))
          END IF
        END IF
      END DO
    END IF

    ! Compute organized entrainment
    ! for use at next level

    DO jl = 1, klon
      IF (loflag(jl) .AND. ktype(jl)==1) THEN
        zbuoyz = g*(ptu(jl,jk)-ptenh(jl,jk))/ptenh(jl,jk) + &
                 g*0.608*(pqu(jl,jk)-pqenh(jl,jk)) - g*plu(jl,jk)
        zbuoyz = MAX(zbuoyz,0.0)
        zdz = (pgeo(jl,jk-1)-pgeo(jl,jk))/g
        zdrodz = -LOG(pten(jl,jk-1)/pten(jl,jk))/zdz - g/(rd*ptenh(jl,jk))
        zbuoy(jl) = zbuoy(jl) + zbuoyz*zdz
        zoentr(jl,jk-1) = zbuoyz*0.5/(1.+zbuoy(jl)) + zdrodz
        IF (lamip2) THEN
          zoentr(jl,jk-1) = MIN(zoentr(jl,jk-1),3.E-4)
        ELSE
          zoentr(jl,jk-1) = MIN(zoentr(jl,jk-1),1.E-3)
        END IF
        zoentr(jl,jk-1) = MAX(zoentr(jl,jk-1),0.)

      END IF
    END DO

  END DO

!-- 5.  Determine convective fluxes above non-buoyancy level
!       (Note: Cloud variables like t,q and l are not
!              affected by detrainment and are already known
!              from previous calculations above)

  DO jl = 1, klon
    IF (kctop(jl)==klevm1) ldcum(jl) = .FALSE.
    kcbot(jl) = MAX(kcbot(jl),kctop(jl))
  END DO
  is = 0
  DO jl = 1, klon
    is = is + MERGE(1,0,ldcum(jl))
  END DO
  kcum = is

  IF (is==0) RETURN

!DIR$ IVDEP
!OCL NOVREC
  DO jl = 1, klon
    IF (ldcum(jl)) THEN
      jk = kctop(jl) - 1
      zzdmf = cmfctop
      zdmfde(jl) = (1.-zzdmf)*pmfu(jl,jk+1)
      plude(jl,jk) = zdmfde(jl)*plu(jl,jk+1)
      pmfu(jl,jk) = pmfu(jl,jk+1) - zdmfde(jl)
      pdmfup(jl,jk) = 0.
      pmfus(jl,jk) = (cpd*ptu(jl,jk)+pgeoh(jl,jk))*pmfu(jl,jk)
      pmfuq(jl,jk) = pqu(jl,jk)*pmfu(jl,jk)
      pmful(jl,jk) = plu(jl,jk)*pmfu(jl,jk)
      plude(jl,jk-1) = pmful(jl,jk)
    END IF
  END DO
  DO jt = 1, ktrac
    DO jl = 1, klon
      IF (ldcum(jl)) THEN
        jk = kctop(jl) - 1
        pmfuxt(jl,jk,jt) = pxtu(jl,jk,jt)*pmfu(jl,jk)
      END IF
    END DO
  END DO

  IF (lmfdudv) THEN
!DIR$ IVDEP
!OCL NOVREC
    DO jl = 1, klon
      IF (ldcum(jl)) THEN
        jk = kctop(jl) - 1
        puu(jl,jk) = puu(jl,jk+1)
        pvu(jl,jk) = pvu(jl,jk+1)
      END IF
    END DO
  END IF

  RETURN
END SUBROUTINE cuasc
