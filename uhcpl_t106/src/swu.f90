!+ shortwave radiation, absorber amounts
!+ $Id: swu.f90,v 1.10 1999/07/20 14:21:45 m214003 Exp $

SUBROUTINE swu(kdlon,kflev,psct,pcardi,ppmb,ppsol,prmu0,ptave,pwv,paki,pdsig, &
&      pfact,prmu,psec,pud)

  ! Description:
  !
  ! Computes the absorber amounts used in shortwave radiation
  ! calculations
  !
  ! Method:
  !
  ! *swu* is called by *sw*
  !
  ! 1. Computes absorber amounts with temperature and pressure
  !    scaling.
  !
  ! Reference:
  ! See radiation's part of the ECMWF research department
  ! documentation, and fouquart and bonnel (1980)
  !
  ! Authors:
  !
  ! J.-J. Morcrette, ECMWF, July 1989, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_shortwave
  USE mo_radiation

  IMPLICIT NONE

  !  Scalar arguments 
  REAL :: pcardi, psct
  INTEGER :: kdlon, kflev

  !  Array arguments 
  REAL :: paki(kdlon,2), pdsig(kdlon,kflev), pfact(kdlon), &
&      ppmb(kdlon,kflev+1), ppsol(kdlon), prmu(kdlon), prmu0(kdlon), &
&      psec(kdlon), ptave(kdlon,kflev), pud(kdlon,3,kflev+1), pwv(kdlon,kflev)

  !  Local scalars: 
  REAL :: zdsco2, zdsh2o, zn175, zn190, zrth, zrtu, zsign, zwh2o
  INTEGER :: jk, jkl, jkp1, jl

  !  Local arrays: 
  REAL :: zo175(kdlon), zo190(kdlon), zr1(kdlon), zr2(kdlon), zsigo(kdlon), &
&      zu1d(kdlon), zu2d(kdlon)

  !  External subroutines 
  EXTERNAL swtt

  !  Intrinsic functions 
#ifdef ECLIB
  REAL :: ALOGHF
!DIR$ VFUNCTION ALOGHF
#define LOG(x)  ALOGHF(x)
#else
  INTRINSIC LOG
#endif
  INTRINSIC MAX, SQRT


  !  Executable statements 

!-- 1. Computes amounts of absorbers

!-- 1.1 Initializes quantities

  DO jl = 1, kdlon
    pud(jl,1,kflev+1) = 0.
    pud(jl,2,kflev+1) = 0.
    pud(jl,3,kflev+1) = 0.
    pfact(jl) = prmu0(jl)*psct
    prmu(jl) = SQRT(1224.*prmu0(jl)*prmu0(jl)+1.)/35.
    psec(jl) = 1./prmu(jl)
  END DO

!-- 1.2 Amounts of absorbers

  DO jl = 1, kdlon
    zu1d(jl) = 0.
    zu2d(jl) = 0.
    zo175(jl) = ppsol(jl)**rpdu1
    zo190(jl) = ppsol(jl)**rpdh1
    zsigo(jl) = ppsol(jl)
  END DO

  DO jk = 1, kflev
    jkp1 = jk + 1
    jkl = kflev + 1 - jk
    DO jl = 1, kdlon
      zrth = (rth2o/ptave(jl,jk))**rtdh2o
      zrtu = (rtumg/ptave(jl,jk))**rtdumg
      zwh2o = MAX(pwv(jl,jkl),zepscq)
      zsign = 100.*ppmb(jl,jkp1)
      pdsig(jl,jk) = (zsigo(jl)-zsign)/ppsol(jl)
      zn175 = zsign**rpdu1
      zn190 = zsign**rpdh1
      zdsco2 = zo175(jl) - zn175
      zdsh2o = zo190(jl) - zn190
      pud(jl,1,jk) = rpnh*zdsh2o*zwh2o*zrth
      pud(jl,2,jk) = rpnu*zdsco2*pcardi*zrtu
      zu1d(jl) = zu1d(jl) + pud(jl,1,jk)
      zu2d(jl) = zu2d(jl) + pud(jl,2,jk)
      zsigo(jl) = zsign
      zo175(jl) = zn175
      zo190(jl) = zn190
    END DO
  END DO

!-- 1.3 Computes clear-sky grey absorption coefficients

  DO jl = 1, kdlon
    zu1d(jl) = zu1d(jl)*psec(jl)
    zu2d(jl) = zu2d(jl)*psec(jl)
  END DO

  CALL swtt(kdlon,2,1,zu1d,zr1)

  DO jl = 1, kdlon
    paki(jl,1) = -LOG(zr1(jl))/zu1d(jl)
  END DO

  CALL swtt(kdlon,2,2,zu2d,zr2)

  DO jl = 1, kdlon
    paki(jl,2) = -LOG(zr2(jl))/zu2d(jl)
  END DO

  RETURN
END SUBROUTINE swu
