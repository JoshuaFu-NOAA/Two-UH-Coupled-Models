!+ shortwave radiation, first spectral interval
!+ $Id: sw1s.f90,v 1.4 1998/10/28 12:34:07 m214003 Exp $

SUBROUTINE sw1s(kdlon,kflev,kewaer,kaerh,knu,ppmb,paer,palbs,pcg,pcldsw, &
&      pdsig,pomega,poz,prmu,psec,ptau,pud,pum,pfd,pfu)

  ! Description:
  !
  ! This routine computes the shortwave radiation fluxes in two
  ! spectral intervals following fouquart and bonnel (1980).
  !
  ! Method:
  !
  ! *sw1s* is called from *sw*.
  !
  ! 1. Computes upward and downward fluxes corresponding to
  !    continuum scattering
  ! 2. Multiply by ozone transmission function
  !
  ! Reference:
  ! See radiation's part of the ECMWF research department
  ! documentation, and fouquart and bonnel (1980)
  !
  ! Authors:
  !
  ! J.-J. Morcrette, ECMWF, July 1989, original source
  ! R. Van Dorland, KNMI, May 1995, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_shortwave
  USE mo_constants

  IMPLICIT NONE

  !  Scalar arguments 
  INTEGER :: kdlon, kewaer, kflev, knu

  !  Array arguments 
  REAL :: paer(kdlon,kflev,5+kewaer), palbs(kdlon,2), pcg(kdlon,2,kflev), &
&      pcldsw(kdlon,kflev), pdsig(kdlon,kflev), pfd(kdlon,kflev+1), &
&      pfu(kdlon,kflev+1), pomega(kdlon,2,kflev), poz(kdlon,kflev), &
&      ppmb(kdlon,kflev+1), prmu(kdlon), psec(kdlon), ptau(kdlon,2,kflev), &
&      pud(kdlon,3,kflev+1), pum(kdlon,kflev+1)
  INTEGER :: kaerh(kdlon,kflev)

  !  Local scalars: 
  REAL :: zozfac
  INTEGER :: ikl, ikm1, jaj, jk, jl

  !  Local arrays: 
  REAL :: zcgaz(kdlon,kflev), zpizaz(kdlon,kflev), zr1(kdlon), zr2(kdlon), &
&      zray1(kdlon,kflev+1), zray2(kdlon,kflev+1), zrayl(kdlon), &
&      zrefz(kdlon,2,kflev+1), zrj(kdlon,6,kflev+1), zrk(kdlon,6,kflev+1), &
&      zrmue(kdlon,kflev+1), ztauaz(kdlon,kflev), ztra1(kdlon,kflev+1), &
&      ztra2(kdlon,kflev+1), zw1(kdlon), zw2(kdlon)

  !  External subroutines 
  EXTERNAL swr, swtt


  !  Executable statements 

!-- 1. First spectral interval (0.25-0.68 micron)

!-- 1.1 Optical thickness for rayleigh scattering

  DO jl = 1, kdlon
    zrayl(jl) = rayc(knu,1) + prmu(jl)*(rayc(knu,2)+prmu(jl)*(rayc(knu,3) &
&      +prmu(jl)*(rayc(knu,4)+prmu(jl)*(rayc(knu,5)+prmu(jl)*rayc(knu,6)))))
  END DO

!-- 2. Continuum scattering calculations

  CALL swr(kdlon,kflev,kewaer,kaerh,knu,ppmb,paer,palbs,pcg,pcldsw,pdsig, &
&      pomega,zrayl,psec,ptau,zcgaz,zpizaz,zray1,zray2,zrefz,zrj,zrk,zrmue, &
&      ztauaz,ztra1,ztra2)

!-- 3. Ozone absorption

!-- 3.1 Downward fluxes

  jaj = 2
  zozfac = 46.6968/g

  DO jl = 1, kdlon
    zw1(jl) = 0.
    zw2(jl) = 0.
    zr1(jl) = 1.
    pfd(jl,kflev+1) = zrj(jl,jaj,kflev+1)*sun_data(knu)
    pud(jl,3,kflev+1) = 0.
  END DO
  DO jk = 1, kflev
    ikl = kflev + 1 - jk
    DO jl = 1, kdlon
      zw1(jl) = zw1(jl) + pud(jl,1,ikl)/zrmue(jl,ikl)
      zw2(jl) = zw2(jl) + poz(jl,jk)*zozfac/zrmue(jl,ikl)
    END DO

    CALL swtt(kdlon,knu,1,zw1,zr1)
    CALL swtt(kdlon,knu,3,zw2,zr2)

    DO jl = 1, kdlon
      pfd(jl,ikl) = zr1(jl)*zr2(jl)*zrj(jl,jaj,ikl)*sun_data(knu)
      pud(jl,3,ikl) = zw2(jl)
    END DO
  END DO

!-- 3.2 Upward fluxes

  DO jl = 1, kdlon
    pfu(jl,1) = palbs(jl,knu)*pfd(jl,1)
    pum(jl,1) = pud(jl,3,1)
  END DO

  DO jk = 2, kflev + 1
    ikm1 = jk - 1
    ikl = kflev + 2 - jk
    DO jl = 1, kdlon
      zw1(jl) = zw1(jl) + pud(jl,1,ikm1)*1.66
      zw2(jl) = zw2(jl) + poz(jl,ikl)*zozfac*1.66
    END DO

    CALL swtt(kdlon,knu,1,zw1,zr1)
    CALL swtt(kdlon,knu,3,zw2,zr2)

    DO jl = 1, kdlon
      pfu(jl,jk) = zr1(jl)*zr2(jl)*zrk(jl,jaj,jk)*sun_data(knu)
      pum(jl,jk) = zw2(jl)
    END DO
  END DO

  RETURN
END SUBROUTINE sw1s
