!+ continuum scattering computations
!+ $Id: swr.f90,v 1.7 1999/07/20 14:21:44 m214003 Exp $

SUBROUTINE swr(kdlon,kflev,kewaer,kaerh,knu,ppmb,paer,palbs,pcg,pcldsw,pdsig, &
&      pomega,prayl,psec,ptau,pcgaz,ppizaz,pray1,pray2,prefz,prj,prk,prmue, &
&      ptauaz,ptra1,ptra2)

  ! Description:
  !
  ! Computes the reflectivity and transmissivity in case of
  ! continuum scattering
  !
  ! Method:
  !
  ! *swr* is called either from *sw1s* or from *sw2s*
  !
  ! 1. Computes continuum fluxes corresponding to aerosol
  !    or/and rayleigh scattering (no molecular gas absorption)
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

  USE mo_aerosols
  USE mo_radiation
  USE mo_radint
  USE mo_constants

  IMPLICIT NONE

  !  Scalar arguments 
  INTEGER :: kdlon, kewaer, kflev, knu

  !  Array arguments 
  REAL :: paer(kdlon,kflev,5+kewaer), palbs(kdlon,2), pcg(kdlon,2,kflev), &
&      pcgaz(kdlon,kflev), pcldsw(kdlon,kflev), pdsig(kdlon,kflev), &
&      pomega(kdlon,2,kflev), ppizaz(kdlon,kflev), ppmb(kdlon,kflev+1), &
&      pray1(kdlon,kflev+1), pray2(kdlon,kflev+1), prayl(kdlon), &
&      prefz(kdlon,2,kflev+1), prj(kdlon,6,kflev+1), prk(kdlon,6,kflev+1), &
&      prmue(kdlon,kflev+1), psec(kdlon), ptau(kdlon,2,kflev), &
&      ptauaz(kdlon,kflev), ptra1(kdlon,kflev+1), ptra2(kdlon,kflev+1)
  INTEGER :: kaerh(kdlon,kflev)

  !  Local scalars: 
  REAL :: zbmu0, zbmu1, zcltest, zcorae, zcorcd, zden, zden1, zeodae, zfacoa, &
&      zfacoc, zff, zgap, zgar, zmu1, zmue, zr21, zr22, zratio, zre11, zrhodz, &
&      zss1, zto, ztray, zww
  INTEGER :: iae, icae, ih, jae, jaer, jaj, jk, jkl, jklp1, jkm1, jl

  !  Local arrays: 
  REAL :: zc1i(kdlon,kflev+1), zgg(kdlon), zre1(kdlon), zre2(kdlon), &
&      zref(kdlon), zrmuz(kdlon), zrneb(kdlon), zto1(kdlon), &
&      ztr(kdlon,2,kflev+1), ztr1(kdlon), ztr2(kdlon), zw(kdlon)

  !  External subroutines 
  EXTERNAL swde

  !  Intrinsic functions 
#ifdef ECLIB
  REAL :: EXPHF
!DIR$ VFUNCTION EXPHF
#define EXP(x)  EXPHF(x)
#else
  INTRINSIC EXP
#endif
  INTRINSIC MAX, MIN, SUM


  !  Executable statements 

!-- 0.9 Threshold for zero aerosol optical depth

  zeodae = 1.E-15

!-- 1. Optical parameters for aerosols and rayleigh

  DO jk = 1, kflev + 1
    DO ja = 1, 6
      DO jl = 1, kdlon
        prj(jl,ja,jk) = 0.
        prk(jl,ja,jk) = 0.
      END DO
    END DO
  END DO

  DO jk = 1, kflev
    DO jl = 1, kdlon
      pcgaz(jl,jk) = 0.
      ppizaz(jl,jk) = 0.
      ptauaz(jl,jk) = 0.
    END DO
    DO jae = 1, 5
      DO jl = 1, kdlon
        ptauaz(jl,jk) = ptauaz(jl,jk) + paer(jl,jk,jae)*taua(knu,jae)
        ppizaz(jl,jk) = ppizaz(jl,jk) + paer(jl,jk,jae)*taua(knu,jae)*piza( &
&            knu,jae)
        pcgaz(jl,jk) = pcgaz(jl,jk) + paer(jl,jk,jae)*taua(knu,jae)*piza(knu, &
&            jae)*cga(knu,jae)
      END DO
    END DO

    ! Contribution gads aerosols
    DO jaer = 1, kewaer
      icae = ndfaer(jaer)
      iae = jaer + 5
      DO jl = 1, kdlon
        ih = kaerh(jl,kflev+1-jk)
        zrhodz = 10.*(ppmb(jl,jk)-ppmb(jl,jk+1))/g
        ptauaz(jl,jk) = ptauaz(jl,jk) + paer(jl,jk,iae)*tauan(ih,knu,icae)* &
&            fcvaer(icae)*zrhodz
        ppizaz(jl,jk) = ppizaz(jl,jk) + paer(jl,jk,iae)*tauan(ih,knu,icae)* &
&            pizan(ih,knu,icae)*fcvaer(icae)*zrhodz
        pcgaz(jl,jk) = pcgaz(jl,jk) + paer(jl,jk,iae)*tauan(ih,knu,icae)* &
&            pizan(ih,knu,icae)*cgan(ih,knu,icae)*fcvaer(icae)*zrhodz
      END DO
    END DO

    DO jl = 1, kdlon
      IF (ptauaz(jl,jk)>zeodae) THEN
        pcgaz(jl,jk) = pcgaz(jl,jk)/ppizaz(jl,jk)
        ppizaz(jl,jk) = ppizaz(jl,jk)/ptauaz(jl,jk)
      ELSE
        ptauaz(jl,jk) = zeodae
        ppizaz(jl,jk) = 0.5
        pcgaz(jl,jk) = 0.5
      END IF
      ztray = prayl(jl)*pdsig(jl,jk)
      zratio = ztray/(ztray+ptauaz(jl,jk))
      zgar = pcgaz(jl,jk)
      zff = zgar*zgar
      ptauaz(jl,jk) = ztray + ptauaz(jl,jk)*(1.-ppizaz(jl,jk)*zff)
      pcgaz(jl,jk) = zgar*(1.-zratio)/(1.+zgar)
      ppizaz(jl,jk) = zratio + (1.-zratio)*ppizaz(jl,jk)*(1.-zff)/(1.-ppizaz( &
&          jl,jk)*zff)

    END DO
  END DO

!-- 2. Total effective cloudiness above a given level

  DO jl = 1, kdlon
    zc1i(jl,kflev+1) = 0.
  END DO

  DO jk = 1, kflev
    jkl = kflev + 1 - jk
    jklp1 = jkl + 1
    DO jl = 1, kdlon
      zfacoa = 1. - ppizaz(jl,jkl)*pcgaz(jl,jkl)*pcgaz(jl,jkl)
      zfacoc = 1. - pomega(jl,knu,jkl)*pcg(jl,knu,jkl)*pcg(jl,knu,jkl)
      zcorae = zfacoa*ptauaz(jl,jkl)*psec(jl)
      zcorcd = zfacoc*ptau(jl,knu,jkl)*psec(jl)
      zr21 = EXP(-zcorae)
      zr22 = EXP(-zcorcd)
      zss1 = pcldsw(jl,jkl)*(1.0-zr21*zr22) + (1.0-pcldsw(jl,jkl))*(1.0-zr21)
      zc1i(jl,jkl) = 1.0 - (1.0-zss1)*(1.0-MAX(zss1,zc1i(jl, &
&          jklp1)))/(1.0-MIN(zc1i(jl,jklp1),1.-zepsec))

    END DO
  END DO

!-- 3. Reflectivity/transmissivity for pure scattering

  DO jl = 1, kdlon
    pray1(jl,kflev+1) = 0.
    pray2(jl,kflev+1) = 0.
    prefz(jl,2,1) = palbs(jl,knu)
    prefz(jl,1,1) = palbs(jl,knu)
    ptra1(jl,kflev+1) = 1.
    ptra2(jl,kflev+1) = 1.
  END DO

  DO jk = 2, kflev + 1
    jkm1 = jk - 1
    DO jl = 1, kdlon
      zrneb(jl) = pcldsw(jl,jkm1)
      zre1(jl) = 0.
      ztr1(jl) = 0.
      zre2(jl) = 0.
      ztr2(jl) = 0.

!-- 3.1 Equivalent zenith angle

      zmue = (1.-zc1i(jl,jk))*psec(jl) + zc1i(jl,jk)*1.66
      prmue(jl,jk) = 1./zmue

!-- 3.2 Reflect./transmissivity due to rayleigh and aerosols

      zgap = pcgaz(jl,jkm1)
      zbmu0 = 0.5 - 0.75*zgap/zmue
      zww = ppizaz(jl,jkm1)
      zto = ptauaz(jl,jkm1)
      zden = 1. + (1.-zww+zbmu0*zww)*zto*zmue + (1-zww)*(1.-zww+2.*zbmu0*zww) &
&          *zto*zto*zmue*zmue
      pray1(jl,jkm1) = zbmu0*zww*zto*zmue/zden
      ptra1(jl,jkm1) = 1./zden

      zmu1 = 0.5
      zbmu1 = 0.5 - 0.75*zgap*zmu1
      zden1 = 1. + (1.-zww+zbmu1*zww)*zto/zmu1 + (1-zww)*(1.-zww+2.*zbmu1*zww &
&          )*zto*zto/zmu1/zmu1
      pray2(jl,jkm1) = zbmu1*zww*zto/zmu1/zden1
      ptra2(jl,jkm1) = 1./zden1

!-- 3.3 Effect of cloud layer

      zw(jl) = pomega(jl,knu,jkm1)
      zto1(jl) = ptau(jl,knu,jkm1)/zw(jl) + ptauaz(jl,jkm1)/ppizaz(jl,jkm1)
      zr21 = ptau(jl,knu,jkm1) + ptauaz(jl,jkm1)
      zr22 = ptau(jl,knu,jkm1)/zr21
      zgg(jl) = zr22*pcg(jl,knu,jkm1) + (1.-zr22)*pcgaz(jl,jkm1)
      zw(jl) = zr21/zto1(jl)
      zref(jl) = prefz(jl,1,jkm1)
      zrmuz(jl) = prmue(jl,jk)
    END DO

    zcltest = SUM(zrneb(1:kdlon))
    IF (zcltest>zepsc*1.E3) THEN

      CALL swde(kdlon,zgg,zref,zrmuz,zto1,zw,zre1,zre2,ztr1,ztr2)

    END IF

    DO jl = 1, kdlon

      prefz(jl,1,jk) = (1.-zrneb(jl))*(pray1(jl,jkm1)+prefz(jl,1,jkm1)*ptra1( &
&          jl,jkm1)*ptra2(jl,jkm1)/(1.-pray2(jl,jkm1)*prefz(jl,1, &
&          jkm1))) + zrneb(jl)*zre2(jl)

      ztr(jl,1,jkm1) = zrneb(jl)*ztr2(jl) + (ptra1(jl,jkm1)/(1.-pray2(jl, &
&          jkm1)*prefz(jl,1,jkm1)))*(1.-zrneb(jl))

      prefz(jl,2,jk) = (1.-zrneb(jl))*(pray1(jl,jkm1)+prefz(jl,2,jkm1)*ptra1( &
&          jl,jkm1)*ptra2(jl,jkm1)) + zrneb(jl)*zre1(jl)

      ztr(jl,2,jkm1) = zrneb(jl)*ztr1(jl) + ptra1(jl,jkm1)*(1.-zrneb(jl))

    END DO
  END DO
  DO jl = 1, kdlon
    zmue = (1.-zc1i(jl,1))*psec(jl) + zc1i(jl,1)*1.66
    prmue(jl,1) = 1./zmue
  END DO

!-- 3.4 Reflect./transmissivity between surface and level

  IF (knu==1) THEN
    jaj = 2
    DO jl = 1, kdlon
      prj(jl,jaj,kflev+1) = 1.
      prk(jl,jaj,kflev+1) = prefz(jl,1,kflev+1)
    END DO

    DO jk = 1, kflev
      jkl = kflev + 1 - jk
      jklp1 = jkl + 1
      DO jl = 1, kdlon
        zre11 = prj(jl,jaj,jklp1)*ztr(jl,1,jkl)
        prj(jl,jaj,jkl) = zre11
        prk(jl,jaj,jkl) = zre11*prefz(jl,1,jkl)
      END DO
    END DO

  ELSE

    DO jaj = 1, 2
      DO jl = 1, kdlon
        prj(jl,jaj,kflev+1) = 1.
        prk(jl,jaj,kflev+1) = prefz(jl,jaj,kflev+1)
      END DO

      DO jk = 1, kflev
        jkl = kflev + 1 - jk
        jklp1 = jkl + 1
        DO jl = 1, kdlon
          zre11 = prj(jl,jaj,jklp1)*ztr(jl,jaj,jkl)
          prj(jl,jaj,jkl) = zre11
          prk(jl,jaj,jkl) = zre11*prefz(jl,jaj,jkl)
        END DO
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE swr
