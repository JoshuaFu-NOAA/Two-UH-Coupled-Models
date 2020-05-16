!+ radiation interface

SUBROUTINE radlsw(kdlon,kflev,klp2,kewaer,kaerh,kmode,kaer,kcfc,zsct,laland,   &
                 laglac,paer,palbs,pclfr,pmu0,pozon,psurp,pdp,pqw,pq,pqs,paph, &
                 pt,pth,pfls,pflt,pflsc,pfltc)

  ! Description:
  !
  ! Controls radiation computations
  !
  ! Method:
  !
  ! Explicit arguments :
  ! pmu0   : (kdlon)               ; solar angle
  ! palbs  : (kdlon)               ; albedo in the two intervals
  !                                  .25-.68 and .68-4.
  ! ******************************************************** only 1]
  ! pcco2  :                       ; concentration in co2 (pa/pa)
  ! pozon  : (kdlon ,kflev)        ; concentration in ozone (pa/pa)
  ! pth    : (kdlon ,kflev+1)      ; half level temperature
  ! paph   : (klp2  ,kflev+1)      ; half level pressure
  ! psurp  : (klp2)                ; surface pressure
  ! pq     : (kdlon ,kflev)        ; specific humidity pa/pa
  ! pqw    : (kdlon ,kflev)        ; liquid water kg/kg
  ! pqs    : (kdlon ,kflev)        ; saturation specific humidity (kg/kg)
  ! pclfr  : (kdlon ,kflev)        ; cloud fractional cover
  ! paer   : (kdlon,kflev,5+kewaer); aerosol optical thickness (1,,5)
  !                                  tanre et al., 1984
  !                                  aerosol mass mixing ratio (kg/kg)
  !                                  (6,.,5+kewaer) computed in echam4
  ! ==== outputs ===
  ! pflt   : (kdlon ,kflev+1)      ; long wave fluxes
  ! pfls   : (kdlon ,kflev+1)      ; short wave fluxes
  ! pfltc  : (kdlon,2)             ; clear sky lw fluxes (top, bottom)
  ! pflsc  : (kdlon,2)             ; clear sky sw fluxes (top, bottom)
  !                                  only if lsolc=.true.
  !
  ! Authors:
  !
  ! J.-J. Morcrette, ECMWF, February 1988, original source
  ! U. Schlese, DKRZ, July 1993, changed
  ! R. Van Dorland, KNMI, May 1995, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !
  USE mo_control,      ONLY:  nrow,  &! latitude
                              dtime, &! timestep in seconds
                              nlev,  &! number of vertical levels
                              lamip2
  USE mo_constants
  USE mo_cfc
  USE mo_radiation
  USE mo_radint
  USE mo_memory_g3a
  USE mo_memory_g3b

  IMPLICIT NONE

  !  Scalar arguments 
  REAL :: zsct
  INTEGER :: kaer, kcfc, kdlon, kewaer, kflev, klp2, kmode

  !  Array arguments 
  REAL :: paer(kdlon,kflev,5+kewaer), palbs(kdlon), paph(klp2,kflev+1),     &
          pclfr(kdlon,kflev), pdp(kdlon,kflev), pfls(kdlon,kflev+1),        &
          pflsc(kdlon,2), pflt(kdlon,kflev+1), pfltc(kdlon,2), pmu0(kdlon), &
          pozon(kdlon,kflev), pq(kdlon,kflev), pqs(kdlon,kflev),            &
          pqw(kdlon,kflev), psurp(klp2), pt(klp2,kflev), pth(kdlon,kflev+1)
  INTEGER :: kaerh(kdlon,kflev)
  LOGICAL :: laglac(klp2), laland(klp2)

  !  Local scalars: 
  REAL :: tlab, zasic, zasymx1, zasymx2, zcaa, zcab, zcco2, zcdnc, zeffir,        &
          zefflr, zfcfc, zfpi, zg1i, zg1l, zg2i, zg2l, zhey, zhpbl, zicewp, zip1, &
          zip2, zip3, zip4, zliqwp, zlogi, zlogp2, zlogp3, zlp1, zlp2, zlp3,      &
          zlwc, zlwgkg, zmaci, zmacl, zmeanr, zo1i, zo1l, zo2i, zo2l, zomgmx1,    &
          zomgmx2, zpmbm, zprat, zref, zrex, zrh2o, zri0, ztaumx1, ztaumx2, ztc,  &
          ztmelt, zto1i, zto1l, zto2i, zto2l
  INTEGER :: ikl, jcfc, jk, jkl, jklp1, jl, jm, nexp
  LOGICAL :: lo2

  !  Local arrays: 
  REAL :: zalbsu(kdlon,2), zcfcabs(4), zcg(kdlon,2,kflev),                        &
          zcldlw(kdlon,kflev), zcldsw(kdlon,kflev), zfice(kdlon),                 &
          zflux(kdlon,2,kflev+1), zfluxc(kdlon,2,kflev+1), zflwp(kdlon),          &
          zfsdwn(kdlon,kflev+1), zfsup(kdlon,kflev+1), zkap(kdlon), zn1(kdlon),   &
          zn2(kdlon), zomega(kdlon,2,kflev+1), zpmb(kdlon,kflev+1), zpsol(kdlon), &
          zradip(kdlon), zradlp(kdlon), ztau(kdlon,2,kflev), ztave(kdlon,kflev),  &
          ztl(kdlon,kflev+1)

  !  External subroutines 
  EXTERNAL lw, sw

  !  Intrinsic functions 
  INTRINSIC EXP, LOG10, MAX, MERGE, MIN


  !  Executable statements 

  ! Assignment of the pointers for the additionsl G3X field

!DIR$ NOBOUNDS ZCLDSW

!-- 0.3  Set up local constants

  zcaa = 0.0059
  zcab = 0.003102
  ztmelt = tmelt
  zfpi = 1.2

  IF (lamip2) THEN
    zasic = 0.89
  ELSE
    zasic = 0.91
  END IF

  zhey = 1.
  zrh2o = 1000.
  zhpbl = 800.
  nexp = 2
  zrex = 1./3.
  zref = 1.E6*(3.E-9/(4.*api*zrh2o))**zrex
  ! Measurement temperature cfc's (296 k)
  tlab = 296.

!-- 1. Set-up input quantities for radiation

  zcco2 = ccardi

  DO jl = 1, kdlon
    zflux(jl,1,kflev+1) = 0.
    zflux(jl,2,kflev+1) = 0.
    zfluxc(jl,1,kflev+1) = 0.
    zfluxc(jl,2,kflev+1) = 0.
    zpsol(jl) = psurp(jl)
    zpmb(jl,1) = psurp(jl)/100.
  END DO

!-- 1.1 Initialize various fields

  DO jl = 1, kdlon
    zalbsu(jl,1) = palbs(jl)
    zalbsu(jl,2) = palbs(jl)
    IF (laland(jl) .AND. ( .NOT. laglac(jl))) THEN
      IF (lamip2) THEN
        zn1(jl) = 100.
        zn2(jl) = 500.
      ELSE
        zn1(jl) = 50.
        zn2(jl) = 220.
      END IF
      zkap(jl) = 1.143
    ELSE
      zn1(jl) = 50.
      IF (lamip2) THEN
        zn2(jl) = 150.
      ELSE
        zn2(jl) = 100.
      END IF
      zkap(jl) = 1.077
    END IF
  END DO

  DO jk = 1, kflev
    jkl = kflev + 1 - jk
    jklp1 = jkl + 1
    DO jl = 1, kdlon
      zpmb(jl,jk+1) = paph(jl,jkl)/100.
      zflux(jl,1,jk) = 0.
      zflux(jl,2,jk) = 0.
      zfluxc(jl,1,jk) = 0.
      zfluxc(jl,2,jk) = 0.
    END DO
  END DO

  DO jk = 1, kflev
    jkl = kflev + 1 - jk
    jklp1 = jkl + 1
    DO jl = 1, kdlon
      ztl(jl,jk) = pth(jl,jklp1)
      ztave(jl,jk) = pt(jl,jkl)
    END DO
  END DO
  DO jl = 1, kdlon
    ztl(jl,kflev+1) = pth(jl,1)
  END DO

!-- 2. Cloud, aerosol and cfc parameters

  DO jk = 1, kflev
    ikl = kflev + 1 - jk

    ! --- Initialize optical properties to clear sky values
    DO jl = 1, kdlon
      zcldsw(jl,jk) = 0.
      ztau(jl,1,jk) = 0.
      ztau(jl,2,jk) = 0.
      zomega(jl,1,jk) = 1.
      zomega(jl,2,jk) = 1.
      zcg(jl,1,jk) = 0.
      zcg(jl,2,jk) = 0.
      zcldlw(jl,jk) = 0.
    END DO

    DO jl = 1, kdlon
      zpmbm = 0.5*(zpmb(jl,jk)+zpmb(jl,jk+1))
      pclfr(jl,ikl) = MAX(zepsc,MIN(pclfr(jl,ikl),1.-zepsc))
      lo2 = pclfr(jl,ikl) > (2.*zepclc)

      ! --- Liquid water content (g.m-3) and liquid water path (g.m-2)
      zlwgkg = pqw(jl,ikl)*1000.
      zlwgkg = MERGE(zlwgkg/pclfr(jl,ikl),0.,lo2)
      zflwp(jl) = zlwgkg*pdp(jl,ikl)/g
      zlwc = zlwgkg*100.*zpmbm/(rd*ztave(jl,jk))

      ! --- Effective radius for water and ice particles (micrometer)
      !     according to heymsfield (1977) and mc farlane et al (1992)

      zlogi = LOG10(MAX(1.E-4,zhey*zlwc))
      zlogp2 = zlogi*zlogi
      zlogp3 = zlogp2*zlogi
      IF (zpmbm<zhpbl) THEN
        zprat = (MIN(2.,zhpbl/zpmbm))**nexp
        zcdnc = zn1(jl) + (zn2(jl)-zn1(jl))*(EXP(1.-zprat))
      ELSE
        zcdnc = zn2(jl)
      END IF
      zmeanr = 0.001*(0.698+0.366*zlogi+0.122*zlogp2+0.0136*zlogp3)
      zeffir = zfpi*5640.*zmeanr**0.786
      zefflr = zref*zkap(jl)*(zlwc/zcdnc)**zrex
      zradip(jl) = MIN(80.,zeffir)
      zradlp(jl) = MAX(4.,MIN(24.,zefflr))

      ! --- Partitioning between liquid and ice water droplets

      ztc = ztave(jl,jk) - ztmelt
      IF (ztave(jl,jk)<ztmelt) THEN
        zfice(jl) = 1. - (zcaa+(1.-zcaa)*EXP(-zcab*ztc**2))
      ELSE
        zfice(jl) = 0.
      END IF
    END DO

    ! sw optical parameters
    ! liquid water and ice (Rockel, 1993)

    DO jl = 1, kdlon
      IF (ABS(zflwp(jl)) > 0.) THEN
        zicewp = zflwp(jl)*zfice(jl)
        zliqwp = zflwp(jl)*(1.-zfice(jl))
        zlp1 = LOG10(zradlp(jl))
        zlp2 = zlp1*zlp1
        zlp3 = zlp2*zlp1
        zto1l = zliqwp*1.8706058417*zradlp(jl)**(-1.0756364457)
        zto2l = zliqwp*1.9655460426*zradlp(jl)**(-1.0778999732)
        zg1l = 0.78756640717 + 0.10660598895*zlp1 - 0.031012468401*zlp2
        zg2l = 0.79208639802 - 0.044930076174*zlp1 + 0.18980672305*zlp2 - &
               0.082590933352*zlp3
        zo1l = 0.99999999
        zo2l = 0.9854369057 + 0.013584242533*zlp1 - 0.024856960461*zlp2 + &
               0.0055918314369*zlp3

        zip1 = LOG10(zradip(jl))
        zip2 = zip1*zip1
        zip3 = zip2*zip1
        zip4 = zip2*zip2
        zto1i = zicewp*1.9056067246*zradip(jl)**(-1.0318784654)
        zto2i = zicewp*2.1666771102*zradip(jl)**(-1.0634702711)
        zg1i = 0.7700034985 + 0.19598466851*zip1 - 0.11836420885*zip2 + &
               0.025209205131*zip3
        zg2i = 0.83631171237 - 0.19965998649*zip1 + 0.46130320487*zip2 - &
               0.29719270332*zip3 + 0.062554483594*zip4
        zg1i = zg1i*zasic
        zg2i = zg2i*zasic
        zo1i = 0.99999999
        zo2i = 0.98475089485 + 0.0053152066002*zip1 - 0.0061150583857*zip2 - &
               0.0032775655896*zip3

        ! - MIX of WATER and ICE CLOUDS
        ztaumx1 = zto1l + zto1i
        ztaumx2 = zto2l + zto2i
        zomgmx1 = zto1l*zo1l + zto1i*zo1i
        zomgmx2 = zto2l*zo2l + zto2i*zo2i
        zasymx1 = zto1l*zo1l*zg1l + zto1i*zo1i*zg1i
        zasymx2 = zto2l*zo2l*zg2l + zto2i*zo2i*zg2i

        zasymx1 = zasymx1/zomgmx1
        zasymx2 = zasymx2/zomgmx2
        zomgmx1 = zomgmx1/ztaumx1
        zomgmx2 = zomgmx2/ztaumx2

        ! --- sw final cloud optical parameters

        zcldsw(jl,jk) = pclfr(jl,ikl)
        ztau(jl,1,jk) = ztaumx1
        ztau(jl,2,jk) = ztaumx2
        zomega(jl,1,jk) = zomgmx1
        zomega(jl,2,jk) = zomgmx2
        zcg(jl,1,jk) = zasymx1
        zcg(jl,2,jk) = zasymx2

        ! lw optical parameters

        ! Liquid water and ice (Rockel, 1993)

        zmacl = 0.025520637 + 0.2854650784*EXP(-0.088968393014*zradlp(jl))
        zmaci = 0.020219423 + 0.2058619832*EXP(-0.067631070625*zradip(jl))
        zcldlw(jl,jk) = pclfr(jl,ikl)*(1.-EXP(-zmacl*zliqwp)*EXP(-zmaci*zicewp))
      END IF

    END DO
  END DO

  DO jl = 1, kdlon
    zpmb(jl,kflev+1) = 0.0
  END DO

  DO jm = 1, 4
    zcfcabs(jm) = 0.
  END DO

  ! cfc's

  IF (kcfc/=0) THEN
    zfcfc = rd*tlab*1.E+03/101325.
    DO jm = 1, 4
      DO jcfc = 1, ncfc
        zcfcabs(jm) = zcfcabs(jm) + zcfc(jcfc)*cfcwmo(jm,jcfc)
      END DO
      zcfcabs(jm) = zcfcabs(jm)*zfcfc
    END DO
  END IF

  zri0 = zsct

!-- 3. Call longwave radiation code

  CALL lw(kdlon,kflev,kewaer,kaerh,kaer,kcfc,zcco2,zcldlw,pdp,zpmb,pozon,ztl, &
          paer,ztave,pq,zflux,zfluxc,zcfcabs)

!-- 4. Call shortwave radiation code

  CALL sw(kdlon,kflev,kewaer,kaerh,zri0,zcco2,zpsol,zalbsu,pq,pqs,pmu0,zcg, &
          zcldsw,zomega,pozon,zpmb,ztau,ztave,paer,zfsdwn,zfsup)

!-- 5. Fill up the model net lw and sw radiative fluxes

  DO jkl = 1, kflev + 1
    jk = kflev + 1 + 1 - jkl
    DO jl = 1, kdlon
      pfls(jl,jkl) = zfsdwn(jl,jk) - zfsup(jl,jk)
      pflt(jl,jkl) = -zflux(jl,1,jk) - zflux(jl,2,jk)
    END DO
  END DO

  ! Long-wave clear sky fluxes

  DO jl = 1, kdlon
    pfltc(jl,1) = -zfluxc(jl,1,kflev+1) - zfluxc(jl,2,kflev+1)
    pfltc(jl,2) = -zfluxc(jl,1,1) - zfluxc(jl,2,1)
  END DO

!-- 6. Short wave clear sky fluxes

  IF (kmode==1) THEN
!    DO jl = 1, kdlon*kflev
!      zcldsw(jl,1) = 0.
!    END DO
    zcldsw = 0.

    CALL sw(kdlon,kflev,kewaer,kaerh,zri0,zcco2,zpsol,zalbsu,pq,pqs,pmu0,zcg, &
            zcldsw,zomega,pozon,zpmb,ztau,ztave,paer,zfsdwn,zfsup)

    DO jl = 1, kdlon
      pflsc(jl,1) = zfsdwn(jl,kflev+1) - zfsup(jl,kflev+1)
      pflsc(jl,2) = zfsdwn(jl,1) - zfsup(jl,1)
    END DO
  ELSE
    DO jl = 1, kdlon
      pflsc(jl,1) = 0.
      pflsc(jl,2) = 0.
    END DO
  END IF

  RETURN
END SUBROUTINE radlsw
