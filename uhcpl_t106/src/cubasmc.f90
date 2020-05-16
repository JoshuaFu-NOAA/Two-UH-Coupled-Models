!+  Calculates cloud base values for midlevel convection
!+ $Id: cubasmc.f90,v 1.5 1999/09/17 11:18:36 m214089 Exp $

SUBROUTINE cubasmc(klp2,k2lp2,klon,klev,kk,pten,pqen,pqsen,puen,pven,ktrac, &
&      pxten,pxtu,pmfuxt,pverv,pgeo,pgeoh,ldcum,ktype,klab,pmfu,pmfub,pentr, &
&      kcbot,ptu,pqu,plu,puu,pvu,pmfus,pmfuq,pmful,pdmfup,pmfuu,pmfuv)

  ! Description:
  !
  ! This routine calculates cloud base values for midlevel convection
  !
  ! Method:
  !
  ! This routine is called from *cuasc*.
  ! Input are environmental values t,q etc.
  ! It returns cloudbase values for midlevel convection
  !
  ! Authors:
  !
  ! S. Tiedtke, ECMWF, in 1989, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_constants,    ONLY: cpd,     &! specific heat at constant pressure
                             rcpd,    &! rcpd=1./cpd
                             g         ! gravity acceleration
  USE mo_cumulus_flux, ONLY: cmfcmin, &! minimum massflux value (for safety)
                             cmfcmax, &! maximum massflux value allowed for
                             entrmid, &! entrainment rate for midlevel convect.
                             lmfdudv   ! true if cumulus friction is switched on

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: k2lp2, kk, klev, klon, klp2, ktrac

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: pgeo(klp2,klev), pgeoh(klp2,klev), pqen(klp2,klev),  &
&                      pqsen(klp2,klev), pten(klp2,klev), puen(k2lp2,klev), &
&                      pven(k2lp2,klev), pverv(klp2,klev), pxten(klon,klev,ktrac)
  LOGICAL, INTENT (IN) :: ldcum(klp2)

  !  Array arguments with intent(InOut):
  REAL, INTENT (INOUT) :: pdmfup(klp2,klev), pentr(klp2), plu(klp2,klev), &
&                         pmfu(klp2,klev), pmfub(klp2), pmful(klp2,klev), &
&                         pmfuq(klp2,klev), pmfus(klp2,klev), pmfuu(klp2), &
&                         pmfuv(klp2), pmfuxt(klon,klev,ktrac),  &
&                         pqu(klp2,klev), ptu(klp2,klev), puu(klp2,klev), &
&                         pvu(klp2,klev), pxtu(klon,klev,ktrac)
  INTEGER, INTENT (INOUT) :: kcbot(klp2), klab(klp2,klev), ktype(klp2)

  !  Local scalars: 
  REAL :: zzzmb
  INTEGER :: jl, jt

  !  Local arrays: 
  LOGICAL :: llo3(klon)

  !  Intrinsic functions 
  INTRINSIC MAX, MIN


  !  Executable statements 

!-- 1. CALCULATE ENTRAINMENT AND DETRAINMENT RATES

!DIR$ IVDEP
!OCL NOVREC
  DO jl = 1, klon
    llo3(jl) = .FALSE.
    IF ( .NOT. ldcum(jl) .AND. klab(jl,kk+1)==0 .AND. &
&          pqen(jl,kk)>0.90*pqsen(jl,kk)) THEN
      llo3(jl) = .TRUE.
      ptu(jl,kk+1) = (cpd*pten(jl,kk)+pgeo(jl,kk)-pgeoh(jl,kk+1))*rcpd
      pqu(jl,kk+1) = pqen(jl,kk)
      plu(jl,kk+1) = 0.
      zzzmb = MAX(cmfcmin,-pverv(jl,kk)/g)
      zzzmb = MIN(zzzmb,cmfcmax)
      pmfub(jl) = zzzmb
      pmfu(jl,kk+1) = pmfub(jl)
      pmfus(jl,kk+1) = pmfub(jl)*(cpd*ptu(jl,kk+1)+pgeoh(jl,kk+1))
      pmfuq(jl,kk+1) = pmfub(jl)*pqu(jl,kk+1)
      pmful(jl,kk+1) = 0.
      pdmfup(jl,kk+1) = 0.
      kcbot(jl) = kk
      klab(jl,kk+1) = 1
      ktype(jl) = 3
      pentr(jl) = entrmid
      IF (lmfdudv) THEN
        puu(jl,kk+1) = puen(jl,kk)
        pvu(jl,kk+1) = pven(jl,kk)
        pmfuu(jl) = pmfub(jl)*puu(jl,kk+1)
        pmfuv(jl) = pmfub(jl)*pvu(jl,kk+1)
      END IF
    END IF
  END DO

!DIR$ IVDEP
!OCL NOVREC
  DO jt = 1, ktrac
    DO jl = 1, klon
      IF (llo3(jl)) THEN
        pxtu(jl,kk+1,jt) = pxten(jl,kk,jt)
        pmfuxt(jl,kk+1,jt) = pmfub(jl)*pxtu(jl,kk+1,jt)
      END IF
    END DO
  END DO

  RETURN
END SUBROUTINE cubasmc
