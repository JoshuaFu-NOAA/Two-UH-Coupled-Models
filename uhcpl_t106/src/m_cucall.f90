!+ master routine - provides interface for: cumastr (cumulus parameterization)

!#define DEBUG

!MODULE m_cucall
!
!CONTAINS

SUBROUTINE cucall(klp2,k2lp2,klon,klev,klevp1,klevm1,ilab,ktrac,pxtm1,pxtte,               &
                  ptm1,pqm1,pum1,pvm1,pxm1,ptte,pqte,pvom,pvol,pxte,pverv,pqhfl,pxtec,     &
                  papp1,paphp1,pgeo,prsfc,pssfc,paprc,paprs,ktype,ldland,ptopmax,ptopmaxm, &
                  psrain,psevap,psheat,psdiss,psmelt)

  ! ------------------> ldland included for AMIP2

  ! Description:
  !
  ! Master routine - provides interface for: cumastr (cumulus parameterization)
  !
  ! Method:
  !
  ! *cucall* - interface for *cumastr*:
  ! Provides input for cumastr.
  ! Receives updated tendencies, precipitation.
  !
  ! *cucall* is called from *physc*
  !
  ! Authors:
  !
  ! M. Tiedtke, ECMWF, December 1989, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! A. Rhodin, MPI, December 1998, lookup tables removed
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_start_dataset,  ONLY: nstep,   &! current time step
                               nstart    ! time step for start/restart
  USE mo_control,        ONLY: twodt     ! 2.*dtime
  USE mo_constants,      ONLY: vtmpc1    ! vtmpc1=rv/rd-1
#ifndef NOLOOKUP
  USE mo_convect_tables, ONLY: tlucua, jptlucu1, jptlucu2, &
                               lookuperror, lookupoverflow
#else
  USE mo_constants,      ONLY: c2es,    &! constants used for computation 
                               c3ies,   &!   of saturation mixing ratio
                               c3les,   &!
                               c4ies,   &!
                               c4les,   &!
                               tmelt     ! temperature of fusion of ice
#endif

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: k2lp2, klev, klevm1, klevp1, klon, klp2, ktrac

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: paphp1(klp2,klevp1), papp1(klp2,klev), pgeo(klp2,klev), &
                       pqhfl(klp2), pqm1(k2lp2,klev), ptm1(k2lp2,klev),        &
                       ptopmaxm(klp2), pum1(k2lp2,klev), pverv(klp2,klev),     &
                       pvm1(k2lp2,klev), pxm1(klp2,klev), pxte(klp2,klev),     &
                       pxtm1(klon,klev,ktrac)
  LOGICAL, INTENT (IN) :: ldland(klp2)    ! INCLUDED FOR AMIP2 !!

  !  Scalar arguments with intent(InOut):
  REAL, INTENT (INOUT) :: psdiss, psevap, psheat, psmelt, psrain

  !  Array arguments with intent(InOut):
  REAL, INTENT (INOUT) :: paprc(klp2), paprs(klp2), pqte(klp2,klev),           &
                          prsfc(klp2), pssfc(klp2), ptopmax(klp2),             &
                          ptte(klp2,klev), pvol(k2lp2,klev), pvom(k2lp2,klev), &
                          pxtec(klp2,klev), pxtte(klon,klev,ktrac)
  INTEGER, INTENT (INOUT) :: ilab(klp2,klev), ktype(klp2)

  !  Local scalars: 
  REAL    :: ztmst
  INTEGER :: ilevmin, jk, jl, jt
#ifndef NOLOOKUP
  INTEGER :: it
#else
  LOGICAL :: lo
#endif

  !  Local arrays: 
  REAL :: zlu(klp2,klev), zlude(klp2,klev), zmfd(klp2,klev), zmfu(klp2,klev), &
          zqp1(klp2,klev), zqsat(klp2,klev), zqu(klp2,klev), zrain(klp2), &
          ztp1(klp2,klev), ztu(klp2,klev), zup1(klp2,klev), zvp1(klp2,klev), &
          zxp1(klp2,klev), zxtp1(klon,klev,ktrac), zxtu(klon,klev,ktrac)
  INTEGER :: icbot(klp2), ictop(klp2), itopec2(klp2)
  LOGICAL :: locum(klp2)

  !  External subroutines 
  EXTERNAL cumastr

  !  Intrinsic functions 
  INTRINSIC MIN

  !  Executable statements 

#ifndef NOLOOKUP
  lookupoverflow = .FALSE.
#endif

!-- 1. Calculate t,q and qs at main levels

  ztmst = twodt
  IF (nstep==nstart) ztmst = 0.5*twodt

  DO jk = 1, klev
    DO jl = 1, klon
      ztp1(jl,jk) = ptm1(jl,jk) + ptte(jl,jk)*ztmst
      zqp1(jl,jk) = pqm1(jl,jk) + pqte(jl,jk)*ztmst
      zxp1(jl,jk) = pxm1(jl,jk) + pxte(jl,jk)*ztmst
      zup1(jl,jk) = pum1(jl,jk) + pvom(jl,jk)*ztmst
      zvp1(jl,jk) = pvm1(jl,jk) + pvol(jl,jk)*ztmst
#ifndef NOLOOKUP
      it = INT(ztp1(jl,jk)*1000.)
      IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
#ifdef DEBUG
      IF (it<jptlucu1 .OR. it>jptlucu2) THEN
        WRITE(*,*) 'cucall : 1 it=',it,jl,jk
      END IF
#endif
      it = MAX(MIN(it,jptlucu2),jptlucu1)
      zqsat(jl,jk) = tlucua(it)/papp1(jl,jk)
#else
      lo = ztp1(jl,jk) > tmelt
      zqsat(jl,jk) = c2es*EXP(MERGE(c3les,c3ies,lo)*(ztp1(jl,jk)-tmelt) &
               / (ztp1(jl,jk)-MERGE(c4les,c4ies,lo)))/papp1(jl,jk)
#endif
      zqsat(jl,jk) = MIN(0.5,zqsat(jl,jk))
      zqsat(jl,jk) = zqsat(jl,jk)/(1.-vtmpc1*zqsat(jl,jk))
    END DO
#ifndef NOLOOKUP
    IF (lookupoverflow) CALL lookuperror ('cucall      ')
#endif

    DO jt = 1, ktrac
      DO jl = 1, klon
        zxtp1(jl,jk,jt) = pxtm1(jl,jk,jt) + pxtte(jl,jk,jt)*ztmst
      END DO
    END DO

  END DO

  DO jl = 1, klon
    zrain(jl) = 0.
    locum(jl) = .FALSE.
  END DO

!-- 2. Call 'cumastr'(master-routine for cumulus parameterization)
!      -----------------------------------------------------------

  CALL cumastr(klp2,k2lp2,klon,klev,klevp1,klevm1,ilab,ztp1,zqp1,zxp1,zup1,          &
               zvp1,ktrac,ldland,zxtp1,zxtu,pxtte,pverv,zqsat,pqhfl,papp1,paphp1,pgeo,ptte, &
               pqte,pvom,pvol,prsfc,pssfc,paprc,paprs,pxtec,locum,ktype,icbot,ictop, &
               ztu,zqu,zlu,zlude,zmfu,zmfd,zrain,psrain,psevap,psheat,psdiss,psmelt)

!-- 3. Pressure altitude of convective cloud tops

  ilevmin = klev - 4

  DO jl = 1, klon
    itopec2(jl) = klevp1
  END DO

  DO jk = 1, ilevmin
    DO jl = 1, klon
      IF (ilab(jl,jk)==2 .AND. itopec2(jl)==klevp1) THEN
        itopec2(jl) = jk
      END IF
    END DO
  END DO

  DO jl = 1, klon
    IF (itopec2(jl)==1) THEN
      ptopmax(jl) = papp1(jl,1)
    ELSE IF (itopec2(jl)/=klevp1) THEN
      ptopmax(jl) = paphp1(jl,itopec2(jl))
    ELSE
      ptopmax(jl) = 99999.
    END IF
    ptopmax(jl) = MIN(ptopmax(jl),ptopmaxm(jl))
  END DO

  RETURN
END SUBROUTINE cucall

!END MODULE m_cucall
