!+ updates t and q tendencies, precipitation rates does global diagnostics
!+ $Id: cudtdq.f90,v 1.4 1998/10/28 12:28:37 m214003 Exp $

SUBROUTINE cudtdq(klp2,klon,klev,klevp1,ktopm2,paphp1,ldcum,pten,ptte,pqte, &
&      ktrac,pxtte,pmfuxt,pmfdxt,pxtec,pmfus,pmfds,pmfuq,pmfdq,pmful,pdmfup, &
&      pdmfdp,plude,pdpmel,prain,prfl,psfl,psrain,psevap,psheat,psmelt,prsfc, &
&      pssfc,paprc,paprs)

  ! Description:
  !
  ! Updates t and q tendencies, precipitation rates does global diagnostics
  !
  ! Method:
  !
  ! *cudtdq* is called from *cumastr*
  !
  ! Authors:
  !
  ! M. Tiedtke, ECMWF, in July 1986 and December 1989, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  !
  ! for more details see file AUTHORS
  !

  USE mo_control,   ONLY: twodt,  &! 2.*dtime
                          ldsst, lstratiform,  &!fu++
                          nrow     !current latit. line (one entry/task)

  USE mo_constants,  ONLY: rhoh2o, &! density of liquid water
                          tmelt,  &! temperature of fusion of ice
                          alv,    &! latent heat for vaporisation
                          g,      &! gravity acceleration
                          alf,    &! latent heat for fusion
                          als,    &! latent heat for sublimation
                          rcpd     ! rcpd=1./cpd
  USE mo_tracer,     ONLY: lxtconv
!fu++
  USE mo_dsst,       ONLY: prconv_diurnal
  USE mo_stratiform, ONLY: prconv_meso, prconv_dq, prconv_dt,&
                         meso_heat, dtmeso, dqmeso

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: klev, klevp1, klon, klp2, ktopm2, ktrac

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: paphp1(klp2,klevp1), pdmfdp(klp2,klev), &
&      pdmfup(klp2,klev), pdpmel(klp2,klev), plude(klp2,klev), &
&      pmfdq(klp2,klev), pmfds(klp2,klev), pmfdxt(klon,klev,ktrac), &
&      pmful(klp2,klev), pmfuq(klp2,klev), pmfus(klp2,klev), &
&      pmfuxt(klon,klev,ktrac), prain(klp2), prfl(klp2), &
&      psfl(klp2), pten(klp2,klev)
  LOGICAL, INTENT (IN) :: ldcum(klp2)

  !  Scalar arguments with intent(InOut):
  REAL, INTENT (INOUT) :: psevap, psheat, psmelt, psrain

  !  Array arguments with intent(InOut):
  REAL, INTENT (INOUT) :: paprc(klp2), paprs(klp2), pqte(klp2,klev), &
&                         prsfc(klp2), pssfc(klp2), ptte(klp2,klev), &
&                         pxtec(klp2,klev), pxtte(klon,klev,ktrac)

  !  Local scalars: 
  REAL :: zalv, zdiagt, zdiagw, zdqdt, zdtdt, zdxtdt
  INTEGER :: jk, jl, jt,jrow
  LOGICAL :: llo1

  !  Local arrays: 
  REAL :: zmelt(klp2), zsheat(klp2)

  !  Intrinsic functions 
  INTRINSIC MERGE


  !  Executable statements 

     jrow = nrow(2)  !fu++

!-- 1. Specify parameters

  zdiagt = 0.5*twodt
  zdiagw = zdiagt/rhoh2o

!-- 2. Incrementation of t and q tendencies

  DO jl = 1, klon
    zmelt(jl) = 0.
    zsheat(jl) = 0.
  END DO

  DO jk = ktopm2, klev

    IF (jk<klev) THEN
      DO jl = 1, klon
        IF (ldcum(jl)) THEN
          llo1 = (pten(jl,jk)-tmelt) > 0.
          zalv = MERGE(alv,als,llo1)
          zdtdt = (g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*rcpd*(pmfus(jl,jk+1)- &
&              pmfus(jl,jk)+pmfds(jl,jk+1)-pmfds(jl,jk)-alf*pdpmel(jl,jk)- &
&              zalv*(pmful(jl,jk+1)-pmful(jl,jk)-plude(jl, &
&              jk)-(pdmfup(jl,jk)+pdmfdp(jl,jk))))
          ptte(jl,jk) = ptte(jl,jk) + zdtdt
          zdqdt = (g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*(pmfuq(jl,jk+1)-pmfuq( &
&              jl,jk)+pmfdq(jl,jk+1)-pmfdq(jl,jk)+pmful(jl,jk+1)-pmful(jl,jk)- &
&              plude(jl,jk)-(pdmfup(jl,jk)+pdmfdp(jl,jk)))
          pqte(jl,jk) = pqte(jl,jk) + zdqdt
          pxtec(jl,jk) = (g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*plude(jl,jk)
          zsheat(jl) = zsheat(jl) + zalv*(pdmfup(jl,jk)+pdmfdp(jl,jk))
          zmelt(jl) = zmelt(jl) + pdpmel(jl,jk)

!fu++ pick up t/q tendencies for stratiform---------------------------------!
         IF (lstratiform) THEN
         prconv_dq(jl,jk,jrow)=zdqdt
         prconv_dt(jl,jk,jrow)=zdtdt
         END IF
!---------------------------------------------------------------------------!
        END IF
      END DO

      IF (lxtconv) THEN
        DO jt = 1, ktrac
          DO jl = 1, klon
            IF (ldcum(jl)) THEN
              zdxtdt = (g/(paphp1(jl,jk+1)-paphp1(jl, &
&                  jk)))*(pmfuxt(jl,jk+1,jt)-pmfuxt(jl,jk,jt)+ &
&                  pmfdxt(jl,jk+1,jt)-pmfdxt(jl,jk,jt))
              pxtte(jl,jk,jt) = pxtte(jl,jk,jt) + zdxtdt
            END IF
          END DO
        END DO
      END IF

    ELSE
      DO jl = 1, klon
        IF (ldcum(jl)) THEN
          llo1 = (pten(jl,jk)-tmelt) > 0.
          zalv = MERGE(alv,als,llo1)
          zdtdt = -(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*rcpd*(pmfus(jl,jk)+ &
&              pmfds(jl,jk)+alf*pdpmel(jl,jk)-zalv*(pmful(jl,jk)+pdmfup(jl, &
&              jk)+pdmfdp(jl,jk)+plude(jl,jk)))
          ptte(jl,jk) = ptte(jl,jk) + zdtdt
          zdqdt = -(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*(pmfuq(jl,jk)+pmfdq(jl &
&              ,jk)+plude(jl,jk)+(pmful(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)))
          pqte(jl,jk) = pqte(jl,jk) + zdqdt
          pxtec(jl,jk) = (g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*plude(jl,jk)
          zsheat(jl) = zsheat(jl) + zalv*(pdmfup(jl,jk)+pdmfdp(jl,jk))
          zmelt(jl) = zmelt(jl) + pdpmel(jl,jk)
!fu++ pick up t/q tendencies for stratiform---------------------------------!
         IF (lstratiform) THEN
         prconv_dq(jl,jk,jrow)=zdqdt
         prconv_dt(jl,jk,jrow)=zdtdt
         END IF
!---------------------------------------------------------------------------!
        END IF
      END DO

      IF (lxtconv) THEN
        DO jt = 1, ktrac
          DO jl = 1, klon
            IF (ldcum(jl)) THEN
              zdxtdt = -(g/(paphp1(jl,jk+1)-paphp1(jl, &
&                  jk)))*(pmfuxt(jl,jk,jt)+pmfdxt(jl,jk,jt))
              pxtte(jl,jk,jt) = pxtte(jl,jk,jt) + zdxtdt
            END IF
          END DO
        END DO
      END IF

    END IF
  END DO

!-- 3. Update surface fields and do global budgets

  DO jl = 1, klon
    prsfc(jl) = prfl(jl)
    pssfc(jl) = psfl(jl)
    paprc(jl) = paprc(jl) + zdiagw*(prfl(jl)+psfl(jl))
!fu++ SST diurnal cycle
    IF (ldsst) THEN
    prconv_diurnal(jl,jrow)=prconv_diurnal(jl,jrow)+zdiagw*(prfl(jl)+psfl(jl))
    END IF
!fu++ for stratiform parameterization 
    IF (lstratiform) THEN
    prconv_meso(jl,jrow)=prfl(jl)+psfl(jl)
    END IF

    paprs(jl) = paprs(jl) + zdiagw*psfl(jl)
    psheat = psheat + zsheat(jl)
    psrain = psrain + prain(jl)
    psevap = psevap - (prfl(jl)+psfl(jl))
    psmelt = psmelt + zmelt(jl)
  END DO

!========================================================================c
! fu++
!
  IF (lstratiform) THEN 
      call meso_heat
!
!fu++ add meso-heating/moistening here
! 
  DO jk = ktopm2, klev

    IF (jk<klev) THEN
      DO jl = 1, klon
        IF (ldcum(jl)) THEN
          ptte(jl,jk) = ptte(jl,jk)*0.7 + dtmeso(jl,jk,jrow)
          pqte(jl,jk) = pqte(jl,jk)*0.7 + dqmeso(jl,jk,jrow) 
        END IF
      END DO

    ELSE
      DO jl = 1, klon
        IF (ldcum(jl)) THEN
          ptte(jl,jk) = ptte(jl,jk)*0.7 + dtmeso(jl,jk,jrow)
          pqte(jl,jk) = pqte(jl,jk)*0.7 + dqmeso(jl,jk,jrow) 
        END IF
      END DO

    END IF
  END DO

  END IF
!fu++------------------------end------------------------------------------------!

  psevap = psevap + psrain

  RETURN
END SUBROUTINE cudtdq
