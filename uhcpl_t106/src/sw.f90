!+ computes the shortwave radiation fluxes.
!+ $Id: sw.f90,v 1.3 1998/10/28 12:34:05 m214003 Exp $

SUBROUTINE sw(kdlon,kflev,kewaer,kaerh,psct,pcardi,ppsol,palbs,pwv,pqs,prmu0, &
&      pcg,pcldsw,pomega,poz,ppmb,ptau,ptave,paer,pfdown,pfup)

  ! Description:
  !
  ! This routine computes the shortwave radiation fluxes in two
  ! spectral intervals following fouquart and bonnel (1980).
  !
  ! Method:
  !
  ! *sw* is called from *radite*
  !
  ! 1. Computes absorber amounts
  ! 2. Computes fluxes in 1st spectral interval
  ! 3. Computes fluxes in 2nd spectral interval
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

  IMPLICIT NONE

  !  Scalar arguments 
  REAL :: pcardi, psct
  INTEGER :: kdlon, kewaer, kflev

  !  Array arguments 
  REAL :: paer(kdlon,kflev,5+kewaer), palbs(kdlon,2), pcg(kdlon,2,kflev), &
&      pcldsw(kdlon,kflev), pfdown(kdlon,kflev+1), pfup(kdlon,kflev+1), &
&      pomega(kdlon,2,kflev), poz(kdlon,kflev), ppmb(kdlon,kflev+1), &
&      ppsol(kdlon), pqs(kdlon,kflev), prmu0(kdlon), ptau(kdlon,2,kflev), &
&      ptave(kdlon,kflev), pwv(kdlon,kflev)
  INTEGER :: kaerh(kdlon,kflev)

  !  Local scalars: 
  INTEGER :: inu, jk, jl

  !  Local arrays: 
  REAL :: zaki(kdlon,2), zdsig(kdlon,kflev), zfact(kdlon), &
&      zfd(kdlon,kflev+1), zfdown(kdlon,kflev+1), zfu(kdlon,kflev+1), &
&      zfup(kdlon,kflev+1), zrmu(kdlon), zsec(kdlon), zud(kdlon,3,kflev+1), &
&      zum(kdlon,kflev+1)

  !  External subroutines 
  EXTERNAL sw1s, sw2s, swu


  !  Executable statements 

!-- 1. Absorber amounts and other useful quantities

  CALL swu(kdlon,kflev,psct,pcardi,ppmb,ppsol,prmu0,ptave,pwv,zaki,zdsig, &
&      zfact,zrmu,zsec,zud)

!-- 2. First spectral interval (0.25-0.68 micron)

  inu = 1

  CALL sw1s(kdlon,kflev,kewaer,kaerh,inu,ppmb,paer,palbs,pcg,pcldsw,zdsig, &
&      pomega,poz,zrmu,zsec,ptau,zud,zum,zfd,zfu)

!-- 3. Second spectral interval (0.68-4.00 micron)

  inu = 2

  CALL sw2s(kdlon,kflev,kewaer,kaerh,inu,ppmb,paer,zaki,palbs,pcg,pcldsw, &
&      zdsig,pomega,zrmu,zsec,ptau,zud,zum,pwv,pqs,zfdown,zfup)

!-- 4. Fill the diagnostic arrays

  DO jk = 1, kflev + 1
    DO jl = 1, kdlon
      pfup(jl,jk) = (zfup(jl,jk)+zfu(jl,jk))*zfact(jl)
      pfdown(jl,jk) = (zfdown(jl,jk)+zfd(jl,jk))*zfact(jl)
    END DO
  END DO

  RETURN
END SUBROUTINE sw
