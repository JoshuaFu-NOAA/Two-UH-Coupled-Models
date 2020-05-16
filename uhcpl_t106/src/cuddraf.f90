!+ calculates cumulus downdraft descent
!+ $Id: cuddraf.f90,v 1.5 1998/11/03 14:32:01 m214003 Exp $

SUBROUTINE cuddraf(klp2,k2lp2,klon,klev,klevp1,ptenh,pqenh,puen,pven,ktrac, &
&      pxtenh,pxtd,pmfdxt,pgeoh,paphp1,prfl,ptd,pqd,pud,pvd,pmfd,pmfds,pmfdq, &
&      pdmfdp,lddraf)

  ! Description:
  !
  ! Calculates cumulus downdraft descent.
  !
  ! Method:
  !
  ! Produce the vertical profiles for cumulus downdrafts
  ! (i.e. t,q,u and v and fluxes)
  !
  ! This routine is called from *cumastr*.
  ! Input is t,q,p,phi,u,v at half levels.
  ! It returns fluxes of s,q and evaporation rate
  ! and u,v at levels where downdraft occurs
  !
  ! Calculate moist descent for entraining/detraining plume by
  !    a) moving air dry-adiabatically to next level below and
  !    b) correcting for evaporation to obtain saturated state.
  !
  ! *cuadjtq* for adjusting t and q due to evaporation in saturated descent
  !
  ! 1. calculate moist descent for cumulus downdraft by
  !    (a) calculating entrainment rates, assuming
  !        linear decrease of massflux in pbl
  !    (b) doing moist descent - evaporative cooling
  !        and moistening is calculated in *cuadjtq*
  !    (c) checking for negative buoyancy and
  !        specifying final t,q,u,v and downward fluxes
  !
  ! Authors:
  !
  ! M. Tiedtke, ECMWF, in December 1986 and December 1989, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_constants,    ONLY: rd,      &! gas constant for dry air
                             cpd,     &! specific heat at constant pressure
                             rcpd,    &! rcpd=1./cpd
                             vtmpc1,  &! vtmpc1=rv/rd-1
                             g         ! gravity acceleration
  USE mo_cumulus_flux, ONLY: entrdd,  &! entrainment rate for cumulus downdrafts
                             cmfcmin, &! minimum massflux value (for safety)
                             lmfdudv   ! true if cumulus friction is switched on

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: k2lp2, klev, klevp1, klon, klp2, ktrac

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: paphp1(klp2,klevp1), pgeoh(klp2,klev), pqenh(klp2,klev), &
&                      ptenh(klp2,klev), puen(k2lp2,klev), pven(k2lp2,klev), &
&                      pxtenh(klon,klev,ktrac)
  LOGICAL, INTENT (IN) :: lddraf(klp2)

  !  Array arguments with intent(InOut):
  REAL, INTENT (INOUT) :: pdmfdp(klp2,klev), pmfd(klp2,klev), pmfdq(klp2,klev), &
&                         pmfds(klp2,klev), pmfdxt(klon,klev,ktrac), &
&                         pqd(klp2,klev), prfl(klp2), ptd(klp2,klev), &
&                         pud(klp2,klev), pvd(klp2,klev), pxtd(klon,klev,ktrac)

  !  Local scalars: 
  REAL :: zbuo, zdmfdp, zentr, zmfdqk, zmfdsk, zmfduk, zmfdvk, zmfdxtk, &
&      zqdde, zqeen, zsdde, zseen, zxtdde, zxteen
  INTEGER :: icall, ik, is, itopde, jk, jl, jt
  LOGICAL :: llo1

  !  Local arrays: 
  REAL :: zcond(klp2), zdmfde(klp2), zdmfen(klp2), zph(klp2)
  LOGICAL :: llo2(klp2)

  !  External subroutines 
  EXTERNAL cuadjtq

  !  Intrinsic functions 
  INTRINSIC MAX, MERGE, MIN

  !  Executable statements
 
  DO jk = 3, klev
    is = 0
    DO jl = 1, klon
      zph(jl) = paphp1(jl,jk)
      llo2(jl) = lddraf(jl) .AND. pmfd(jl,jk-1) < 0.
      is = is + MERGE(1,0,llo2(jl))
    END DO
    IF (is==0) CYCLE
    DO jl = 1, klon
      IF (llo2(jl)) THEN
        zentr = entrdd*pmfd(jl,jk-1)*rd*ptenh(jl,jk-1)/(g*paphp1(jl,jk-1))* &
&            (paphp1(jl,jk)-paphp1(jl,jk-1))
        zdmfen(jl) = zentr
        zdmfde(jl) = zentr
      END IF
    END DO
    itopde = klev - 2
    IF (jk>itopde) THEN
      DO jl = 1, klon
        IF (llo2(jl)) THEN
          zdmfen(jl) = 0.
          zdmfde(jl) = pmfd(jl,itopde)*(paphp1(jl,jk)-paphp1(jl,jk-1))/ &
&              (paphp1(jl,klevp1)-paphp1(jl,itopde))
        END IF
      END DO
    END IF

    DO jl = 1, klon
      IF (llo2(jl)) THEN
        pmfd(jl,jk) = pmfd(jl,jk-1) + zdmfen(jl) - zdmfde(jl)
        zseen = (cpd*ptenh(jl,jk-1)+pgeoh(jl,jk-1))*zdmfen(jl)
        zqeen = pqenh(jl,jk-1)*zdmfen(jl)
        zsdde = (cpd*ptd(jl,jk-1)+pgeoh(jl,jk-1))*zdmfde(jl)
        zqdde = pqd(jl,jk-1)*zdmfde(jl)
        zmfdsk = pmfds(jl,jk-1) + zseen - zsdde
        zmfdqk = pmfdq(jl,jk-1) + zqeen - zqdde
        pqd(jl,jk) = zmfdqk*(1./MIN(-cmfcmin,pmfd(jl,jk)))
        ptd(jl,jk) = (zmfdsk*(1./MIN(-cmfcmin,pmfd(jl,jk)))-pgeoh(jl,jk))* &
&            rcpd
        ptd(jl,jk) = MIN(400.,ptd(jl,jk))
        ptd(jl,jk) = MAX(100.,ptd(jl,jk))
        zcond(jl) = pqd(jl,jk)
      END IF
    END DO

    DO jt = 1, ktrac
      DO jl = 1, klon
        IF (llo2(jl)) THEN
          zxteen = pxtenh(jl,jk-1,jt)*zdmfen(jl)
          zxtdde = pxtd(jl,jk-1,jt)*zdmfde(jl)
          zmfdxtk = pmfdxt(jl,jk-1,jt) + zxteen - zxtdde
          pxtd(jl,jk,jt) = zmfdxtk*(1./MIN(-cmfcmin,pmfd(jl,jk)))
        END IF
      END DO
    END DO

    ik = jk
    icall = 2
    CALL cuadjtq(klp2,klon,klev,ik,zph,ptd,pqd,llo2,icall)

    DO jl = 1, klon
      IF (llo2(jl)) THEN
        zcond(jl) = zcond(jl) - pqd(jl,jk)
        zbuo = ptd(jl,jk)*(1.+vtmpc1*pqd(jl,jk)) - ptenh(jl,jk)*(1.+vtmpc1* &
&            pqenh(jl,jk))
        llo1 = zbuo < 0. .AND. (prfl(jl)-pmfd(jl,jk)*zcond(jl)>0.)
        pmfd(jl,jk) = MERGE(pmfd(jl,jk),0.,llo1)
        pmfds(jl,jk) = (cpd*ptd(jl,jk)+pgeoh(jl,jk))*pmfd(jl,jk)
        pmfdq(jl,jk) = pqd(jl,jk)*pmfd(jl,jk)
        zdmfdp = -pmfd(jl,jk)*zcond(jl)
        pdmfdp(jl,jk-1) = zdmfdp
        prfl(jl) = prfl(jl) + zdmfdp
      END IF
    END DO

    DO jt = 1, ktrac
      DO jl = 1, klon
        IF (llo2(jl)) THEN
          pmfdxt(jl,jk,jt) = pxtd(jl,jk,jt)*pmfd(jl,jk)
        END IF
      END DO
    END DO

    IF (lmfdudv) THEN
      DO jl = 1, klon
        IF (llo2(jl) .AND. pmfd(jl,jk)<0.) THEN
          zmfduk = pmfd(jl,jk-1)*pud(jl,jk-1) + zdmfen(jl)*puen(jl,jk-1) - &
&              zdmfde(jl)*pud(jl,jk-1)
          zmfdvk = pmfd(jl,jk-1)*pvd(jl,jk-1) + zdmfen(jl)*pven(jl,jk-1) - &
&              zdmfde(jl)*pvd(jl,jk-1)
          pud(jl,jk) = zmfduk*(1./MIN(-cmfcmin,pmfd(jl,jk)))
          pvd(jl,jk) = zmfdvk*(1./MIN(-cmfcmin,pmfd(jl,jk)))
        END IF
      END DO
    END IF

  END DO

  RETURN
END SUBROUTINE cuddraf
