!+ shortwave radiation, 2nd spectral interval
!+ $Id: sw2s.f90,v 1.10 1999/12/10 10:11:02 m214003 Exp $

SUBROUTINE sw2s(kdlon,kflev,kewaer,kaerh,knu,ppmb,paer,paki,palbs,pcg,pcldsw, &
&      pdsig,pomega,prmu,psec,ptau,pud,pum,pwv,pqs,pfdown,pfup)

  ! Description:
  !
  ! This routine computes the shortwave radiation fluxes in the
  ! second spectral interval following fouquart and bonnel (1980).
  !
  ! Method:
  !
  ! *sw2s* is called from *sw*.
  !
  ! 1. Computes reflectivity/transmissivity corresponding to
  !    continuum scattering
  ! 2. Computes reflectivity/transmissivity corresponding for
  !    a grey molecular absorption
  ! 3. Laplace transform on the previous to get effective amounts
  !    of absorbers
  ! 4. Apply H2O and u.m.g. transmission functions
  ! 5. Multiply by ozone transmission function
  !
  ! Reference:
  ! See radiation's part of the ECMWF research department
  ! documentation, and fouquart and bonnel (1980)
  !
  ! Authors:
  !
  ! J.-J. Morcrette, ECMWF, July 1989, original source
  ! U. Schlese, DKRZ, July 1993, changed
  ! R. Van Dorland, KNMI, May 1995, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_shortwave
  USE mo_radiation

  IMPLICIT NONE

  !  Scalar arguments 
  INTEGER :: kdlon, kewaer, kflev, knu

  !  Array arguments 
  REAL :: paer(kdlon,kflev,5+kewaer), paki(kdlon,2), palbs(kdlon,2), &
&      pcg(kdlon,2,kflev), pcldsw(kdlon,kflev), pdsig(kdlon,kflev), &
&      pfdown(kdlon,kflev+1), pfup(kdlon,kflev+1), pomega(kdlon,2,kflev), &
&      ppmb(kdlon,kflev+1), pqs(kdlon,kflev), prmu(kdlon), psec(kdlon), &
&      ptau(kdlon,2,kflev), pud(kdlon,3,kflev+1), pum(kdlon,kflev+1), &
&      pwv(kdlon,kflev)
  INTEGER :: kaerh(kdlon,kflev)

  !  Local scalars: 
  REAL :: zaa, zbb, zcltest, zcneb, zr21, zr22, zre11, zrki, zrmum1, zwh2o
  INTEGER :: ikl, jabs, jaj, jajp, jk, jkki, jkkp4, jkl, jklp1, jkm1, jl, jn, &
&      jn2j, kref

  !  Local arrays: 
  REAL :: zcgaz(kdlon,kflev), zg(kdlon), zgg(kdlon), zpizaz(kdlon,kflev), &
&      zr1(kdlon), zr2(kdlon), zray1(kdlon,kflev+1), zray2(kdlon,kflev+1), &
&      zrayl(kdlon), zre1(kdlon), zre2(kdlon), zref(kdlon), &
&      zrefz(kdlon,2,kflev+1), zrj(kdlon,6,kflev+1), zrk(kdlon,6,kflev+1), &
&      zrl(kdlon,8), zrmue(kdlon,kflev+1), zrmuz(kdlon), zrneb(kdlon), &
&      zs(kdlon), ztauaz(kdlon,kflev), zto1(kdlon), ztr(kdlon,2,kflev+1), &
&      ztr1(kdlon), ztr2(kdlon), ztra1(kdlon,kflev+1), ztra2(kdlon,kflev+1), &
&      zw(kdlon), zw1(kdlon), zw2(kdlon)

  !  External subroutines 
  EXTERNAL swde, swr, swtt

  !  Intrinsic functions 
#ifdef ECLIB
  REAL :: EXPHF,ALOGHF
!DIR$ VFUNCTION EXPHF
!DIR$ VFUNCTION ALOGHF
#define EXP(x)  EXPHF(x)
#define LOG(x)  ALOGHF(x)
#else
  INTRINSIC EXP, LOG
#endif
  INTRINSIC MAX, MIN, SUM

!DIR$ NOBOUNDS

  !  Executable statements 

!-- 1. Second spectral interval (0.68-4.00 micron)

!-- 1.1 Optical thickness for rayleigh scattering

  DO jl = 1, kdlon
    zrmum1 = 1. - prmu(jl)
    zrayl(jl) = rayc(knu,1) + zrmum1*(rayc(knu,2)+zrmum1*(rayc(knu,3) &
&        +zrmum1*(rayc(knu,4)+zrmum1*(rayc(knu,5)+zrmum1*rayc(knu,6)))))
  END DO

!-- 2. Continuum scattering calculations

  CALL swr(kdlon,kflev,kewaer,kaerh,knu,ppmb,paer,palbs,pcg,pcldsw,pdsig, &
&      pomega,zrayl,psec,ptau,zcgaz,zpizaz,zray1,zray2,zrefz,zrj,zrk,zrmue, &
&      ztauaz,ztra1,ztra2)

!-- 3. Scattering calculations with grey molecular absorption

  jn = 2

  DO jabs = 1, 2

!-- 3.1 Surface conditions

    DO jl = 1, kdlon
      zrefz(jl,2,1) = palbs(jl,knu)
      zrefz(jl,1,1) = palbs(jl,knu)
    END DO

!-- 3.2 Introducing cloud effects

    DO jk = 2, kflev + 1
      jkm1 = jk - 1
      ikl = kflev + 1 - jkm1
      DO jl = 1, kdlon
        zrneb(jl) = pcldsw(jl,jkm1)
        IF (jabs==1 .AND. zrneb(jl)>2.*zeelog) THEN
          zwh2o = MAX(pwv(jl,ikl),zeelog)
          zcneb = MAX(zeelog,MIN(zrneb(jl),1.-zeelog))
          zbb = pud(jl,jabs,jkm1)*pqs(jl,ikl)/zwh2o
          zaa = MAX((pud(jl,jabs,jkm1)-zcneb*zbb)/(1.-zcneb),zeelog)
        ELSE
          zaa = pud(jl,jabs,jkm1)
          zbb = zaa
        END IF
        zrki = paki(jl,jabs)
        zs(jl) = EXP(-zrki*zaa*1.66)
        zg(jl) = EXP(-zrki*zaa/zrmue(jl,jk))
        ztr1(jl) = 0.
        zre1(jl) = 0.
        ztr2(jl) = 0.
        zre2(jl) = 0.

        zw(jl) = pomega(jl,knu,jkm1)
        zto1(jl) = ptau(jl,knu,jkm1)/zw(jl) + ztauaz(jl,jkm1)/zpizaz(jl,jkm1) &
&            + zbb*zrki
        zr21 = ptau(jl,knu,jkm1) + ztauaz(jl,jkm1)
        zr22 = ptau(jl,knu,jkm1)/zr21
        zgg(jl) = zr22*pcg(jl,knu,jkm1) + (1.-zr22)*zcgaz(jl,jkm1)
        zw(jl) = zr21/zto1(jl)
        zref(jl) = zrefz(jl,1,jkm1)
        zrmuz(jl) = zrmue(jl,jk)
      END DO

      zcltest = SUM(zrneb(1:kdlon))
      IF (zcltest > zepsc*1.E3) THEN

        CALL swde(kdlon,zgg,zref,zrmuz,zto1,zw,zre1,zre2,ztr1,ztr2)

      END IF

      DO jl = 1, kdlon

        zrefz(jl,2,jk) = (1.-zrneb(jl))*(zray1(jl,jkm1)+zrefz(jl,2,jkm1)* &
&            ztra1(jl,jkm1)*ztra2(jl,jkm1))*zg(jl)*zs(jl) + zrneb(jl)*zre1(jl)

        ztr(jl,2,jkm1) = zrneb(jl)*ztr1(jl) + (ztra1(jl,jkm1))*zg(jl)*(1.-zrneb(jl))

        zrefz(jl,1,jk) = (1.-zrneb(jl))*(zray1(jl,jkm1)+zrefz(jl,1,jkm1)* &
&            ztra1(jl,jkm1)*ztra2(jl,jkm1)/(1.-zray2(jl,jkm1)* &
&            zrefz(jl,1,jkm1)))*zg(jl)*zs(jl) + zrneb(jl)*zre2(jl)

        ztr(jl,1,jkm1) = zrneb(jl)*ztr2(jl) + (ztra1(jl,jkm1)/(1.-zray2(jl,jkm1) &
&            *zrefz(jl,1,jkm1)))*zg(jl)*(1.-zrneb(jl))

      END DO
    END DO

!-- 3.3 Reflect./transmissivity between surface and level

    DO kref = 1, 2

      jn = jn + 1

      DO jl = 1, kdlon
        zrj(jl,jn,kflev+1) = 1.
        zrk(jl,jn,kflev+1) = zrefz(jl,kref,kflev+1)
      END DO

      DO jk = 1, kflev
        jkl = kflev + 1 - jk
        jklp1 = jkl + 1
        DO jl = 1, kdlon
          zre11 = zrj(jl,jn,jklp1)*ztr(jl,kref,jkl)
          zrj(jl,jn,jkl) = zre11
          zrk(jl,jn,jkl) = zre11*zrefz(jl,kref,jkl)
        END DO
      END DO
    END DO
  END DO

!-- 4. Invert grey and continuum fluxes

!-- 4.1 Upward (zrk) and downward (zrj) pseudo-fluxes

  DO jk = 1, kflev + 1
    DO jaj = 1, 5, 2
      jajp = jaj + 1
      DO jl = 1, kdlon
        zrj(jl,jaj,jk) = zrj(jl,jaj,jk) - zrj(jl,jajp,jk)
        zrk(jl,jaj,jk) = zrk(jl,jaj,jk) - zrk(jl,jajp,jk)
      END DO
    END DO
  END DO

  zrj(:,:,:) = MAX(zrj(:,:,:), zeelog)
  zrk(:,:,:) = MAX(zrk(:,:,:), zeelog)

!-- 4.2 Effective absorber amounts by inverse Laplace

  DO jk = 1, kflev + 1
    jkki = 1
    DO jaj = 1, 2
      DO jn = 1, 2
        jn2j = jn + 2*jaj
        jkkp4 = jkki + 4

!-- 4.2.1 Effective absorber amounts

        DO jl = 1, kdlon
          zw1(jl) = LOG(zrj(jl,jn,jk)/zrj(jl,jn2j,jk))/paki(jl,jaj)
        END DO

!-- 4.2.2 Transmission function

        CALL swtt(kdlon,knu,jaj,zw1,zr1)

        DO jl = 1, kdlon
          zrl(jl,jkki) = zr1(jl)
          zw2(jl) = LOG(zrk(jl,jn,jk)/zrk(jl,jn2j,jk))/paki(jl,jaj)
        END DO

        CALL swtt(kdlon,knu,jaj,zw2,zr2)

        DO jl = 1, kdlon
          zrl(jl,jkkp4) = zr2(jl)
        END DO

        jkki = jkki + 1
      END DO
    END DO

!-- 4.3 Upward and downward fluxes with H2O and umg absorption

    DO jl = 1, kdlon
      pfdown(jl,jk) = zrj(jl,1,jk)*zrl(jl,1)*zrl(jl,3) + &
                      zrj(jl,2,jk)*zrl(jl,2)*zrl(jl,4)
      pfup(jl,jk)   = zrk(jl,1,jk)*zrl(jl,5)*zrl(jl,7) + &
                      zrk(jl,2,jk)*zrl(jl,6)*zrl(jl,8)
    END DO
  END DO

!-- 5. Introduction of ozone absorption

  jabs = 3
  DO jk = 1, kflev + 1
    DO jl = 1, kdlon
      zw1(jl) = pud(jl,jabs,jk)
    END DO

    CALL swtt(kdlon,knu,jabs,zw1,zr1)

    DO jl = 1, kdlon
      pfdown(jl,jk) = zr1(jl)*pfdown(jl,jk)*sun_data(knu)
      zw2(jl) = pum(jl,jk)
    END DO

    CALL swtt(kdlon,knu,jabs,zw2,zr2)

    DO jl = 1, kdlon
      pfup(jl,jk) = zr2(jl)*pfup(jl,jk)*sun_data(knu)
    END DO
  END DO

  RETURN
END SUBROUTINE sw2s
