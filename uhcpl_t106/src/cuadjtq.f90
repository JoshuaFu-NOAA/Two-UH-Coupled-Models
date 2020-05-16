!+ produce t,q and l values for cloud ascent
!+ $Id: cuadjtq.f90,v 1.14 2000/09/14 13:14:12 m214003 Exp $

!#define DEBUG

SUBROUTINE cuadjtq(klp2,klon,klev,kk,pp,pt,pq,ldflag,kcall)

  ! Description:
  !
  ! Produce t,q and l values for cloud ascent
  !
  ! Method:
  !
  ! This routine is called from subroutines:
  !    *cubase*   (t and q at condenstion level)
  !    *cuasc*    (t and q at cloud levels)
  !    *cuini*    (environmental t and qs values at half levels)
  ! input are unadjusted t and q values,
  ! it returns adjusted values of t and q
  ! note: input parameter kcall defines calculation as
  !    kcall=0    env. t and qs in*cuini*
  !    kcall=1    condensation in updrafts  (e.g. cubase, cuasc)
  !    kcall=2    evaporation in downdrafts (e.g. cudlfs,cuddraf)
  !
  ! Externals:
  ! 3 lookup tables ( tlucua, tlucub, tlucuc )
  ! for condensation calculations.
  ! the tables are initialised in *setphys*.
  !
  ! Authors:
  !
  ! M. Tiedtke, ECMWF, December 1989, original source
  ! D. Salmond, CRAY(UK), August 1991, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! A. Rhodin, MPI, December 1998, lookup tables removed
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_constants,      ONLY: vtmpc1
#ifndef NOLOOKUP
  USE mo_convect_tables, ONLY: tlucua,   & ! table a
                               tlucub,   & ! table b
                               tlucuc,   & ! table c
                               jptlucu1, jptlucu2, &
                               lookuperror, lookupoverflow
#else
  USE mo_constants,      ONLY: alsdcp, alvdcp, c2es, c3ies, c3les, c4ies, &
                               c4les, c5alscp, c5alvcp, tmelt
#endif

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: kcall, kk, klev, klon, klp2

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: pp(klp2)
  LOGICAL, INTENT (IN) :: ldflag(klp2)

  !  Array arguments with intent(InOut):
  REAL, INTENT (INOUT) :: pq(klp2,klev), pt(klp2,klev)

  !  Local scalars: 
  REAL    :: zcond1, zcor, zqsat
  INTEGER :: isum, jl
#ifndef NOLOOKUP
  INTEGER :: it
#else
  LOGICAL :: lo
#endif

  !  Local arrays: 
  REAL :: zcond(klp2), zqp(klp2)

  !  Intrinsic functions 
  INTRINSIC MAX, MIN
#ifndef NOLOOKUP
  INTRINSIC INT
#endif

  !  Executable statements 

  ! 1.  Calculate condensation and adjust t and q accordingly

#ifndef NOLOOKUP
  lookupoverflow = .FALSE.
#endif

  zcond(:) = 0.

  IF (kcall==1) THEN

    isum = 0
!DIR$ IVDEP
!OCL NOVREC
    DO jl = 1, klon
      IF (ldflag(jl)) THEN
        zqp(jl) = 1./pp(jl)
#ifndef NOLOOKUP
        it = INT(pt(jl,kk)*1000.)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
#ifdef DEBUG
        IF (it<jptlucu1 .OR. it>jptlucu2) THEN
          WRITE(*,*) 'cuadjtq: 1 it=',it,jl,kk
        END IF
#endif
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqsat = tlucua(it)*zqp(jl)
        zqsat = MIN(0.5,zqsat)
        zcor  = 1./(1.-vtmpc1*zqsat)
        zqsat = zqsat*zcor
        zcond(jl) = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor*tlucub(it))
        zcond(jl) = MAX(zcond(jl),0.)
        pt(jl,kk) = pt(jl,kk) + tlucuc(it)*zcond(jl)
#else
        lo = pt(jl,kk) > tmelt
        zqsat = c2es*EXP(MERGE(c3les,c3ies,lo)*(pt(jl,kk)-tmelt) &
            / (pt(jl,kk)-MERGE(c4les,c4ies,lo)))*zqp(jl)
        zqsat = MIN(0.5,zqsat)
        zcor = 1./(1.-vtmpc1*zqsat)
        zqsat = zqsat*zcor
        zcond(jl) = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor &
        * (MERGE(c5alvcp,c5alscp,lo) / (pt(jl,kk)-MERGE(c4les,c4ies,lo))**2))
        zcond(jl) = MAX(zcond(jl),0.)
        pt(jl,kk) = pt(jl,kk) + MERGE(alvdcp, alsdcp, lo)*zcond(jl)
#endif
        pq(jl,kk) = pq(jl,kk) - zcond(jl)
        IF (ABS(zcond(jl)) > 0.) isum = isum + 1
      END IF
    END DO
#ifndef NOLOOKUP
    IF (lookupoverflow) CALL lookuperror ('cuadjtq (1) ')
#endif

    IF (isum/=0) THEN
!DIR$ IVDEP
!OCL NOVREC
      DO jl = 1, klon
        IF (ldflag(jl) .AND. ABS(zcond(jl)) > 0.) THEN
#ifndef NOLOOKUP
          it = INT(pt(jl,kk)*1000.)
          IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
#ifdef DEBUG
          IF (it<jptlucu1 .OR. it>jptlucu2) THEN
            WRITE(*,*) 'cuadjtq: 2 it=',it,jl,kk
          END IF
#endif
          it = MAX(MIN(it,jptlucu2),jptlucu1)
          zqsat = tlucua(it)*zqp(jl)
          zqsat = MIN(0.5,zqsat)
          zcor = 1./(1.-vtmpc1*zqsat)
          zqsat = zqsat*zcor
          zcond1 = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor*tlucub(it))
          pt(jl,kk) = pt(jl,kk) + tlucuc(it)*zcond1
#else
          lo = pt(jl,kk) > tmelt
          zqsat = c2es*EXP(MERGE(c3les,c3ies,lo)*(pt(jl,kk)-tmelt) &
              / (pt(jl,kk)-MERGE(c4les,c4ies,lo)))*zqp(jl)
          zqsat = MIN(0.5,zqsat)
          zcor = 1./(1.-vtmpc1*zqsat)
          zqsat = zqsat*zcor
          zcond1 = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor &
          * (MERGE(c5alvcp,c5alscp,lo) / (pt(jl,kk)-MERGE(c4les,c4ies,lo))**2) )
          pt(jl,kk) = pt(jl,kk) + MERGE(alvdcp, alsdcp, lo) * zcond1
#endif
          pq(jl,kk) = pq(jl,kk) - zcond1
        END IF
      END DO
#ifndef NOLOOKUP
      IF (lookupoverflow) CALL lookuperror ('cuadjtq (2) ')
#endif
    END IF

  END IF

  IF (kcall==2) THEN

    isum = 0
!DIR$ IVDEP
!OCL NOVREC
    DO jl = 1, klon
      IF (ldflag(jl)) THEN
        zqp(jl) = 1./pp(jl)
#ifndef NOLOOKUP
        it = INT(pt(jl,kk)*1000.)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
#ifdef DEBUG
        IF (it<jptlucu1 .OR. it>jptlucu2) THEN
          WRITE(*,*) 'cuadjtq: 3 it=',it,jl,kk
        END IF
#endif
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqsat = tlucua(it)*zqp(jl)
        zqsat = MIN(0.5,zqsat)
        zcor = 1./(1.-vtmpc1*zqsat)
        zqsat = zqsat*zcor
        zcond(jl) = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor*tlucub(it))
        zcond(jl) = MIN(zcond(jl),0.)
        pt(jl,kk) = pt(jl,kk) + tlucuc(it)*zcond(jl)
#else
        lo = pt(jl,kk) > tmelt
        zqsat = c2es*EXP(MERGE(c3les,c3ies,lo)*(pt(jl,kk)-tmelt) &
            / (pt(jl,kk)-MERGE(c4les,c4ies,lo))) *zqp(jl)
        zqsat = MIN(0.5,zqsat)
        zcor = 1./(1.-vtmpc1*zqsat)
        zqsat = zqsat*zcor
        zcond(jl) = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor &
        * (MERGE(c5alvcp,c5alscp,lo) / (pt(jl,kk)-MERGE(c4les,c4ies,lo))**2))
        zcond(jl) = MIN(zcond(jl),0.)
        pt(jl,kk) = pt(jl,kk) + MERGE(alvdcp, alsdcp, lo) * zcond(jl)
#endif
        pq(jl,kk) = pq(jl,kk) - zcond(jl)
        IF (ABS(zcond(jl)) > 0.) isum = isum + 1
      END IF
    END DO
#ifndef NOLOOKUP
    IF (lookupoverflow) CALL lookuperror ('cuadjtq (3) ')
#endif

    IF (isum/=0) THEN
!DIR$ IVDEP
!OCL NOVREC
      DO jl = 1, klon
        IF (ldflag(jl) .AND. ABS(zcond(jl)) > 0.) THEN
#ifndef NOLOOKUP
          it = INT(pt(jl,kk)*1000.)
          IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
#ifdef DEBUG
          IF (it<jptlucu1 .OR. it>jptlucu2) THEN
            WRITE(*,*) 'cuadjtq: 4 it=',it,jl,kk
          END IF
#endif
          it = MAX(MIN(it,jptlucu2),jptlucu1)
          zqsat = tlucua(it)*zqp(jl)
          zqsat = MIN(0.5,zqsat)
          zcor = 1./(1.-vtmpc1*zqsat)
          zqsat = zqsat*zcor
          zcond1 = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor*tlucub(it))
          pt(jl,kk) = pt(jl,kk) + tlucuc(it)*zcond1
#else
          lo = pt(jl,kk) > tmelt
          zqsat = c2es*EXP(MERGE(c3les,c3ies,lo)*(pt(jl,kk)-tmelt) &
              / (pt(jl,kk)-MERGE(c4les,c4ies,lo))) *zqp(jl)
          zqsat = MIN(0.5,zqsat)
          zcor = 1./(1.-vtmpc1*zqsat)
          zqsat = zqsat*zcor
          zcond1 = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor &
          * (MERGE(c5alvcp,c5alscp,lo) / (pt(jl,kk)-MERGE(c4les,c4ies,lo))**2))
          pt(jl,kk) = pt(jl,kk) + MERGE(alvdcp, alsdcp, lo) *zcond1
#endif
          pq(jl,kk) = pq(jl,kk) - zcond1
        END IF
      END DO
#ifndef NOLOOKUP
      IF (lookupoverflow) CALL lookuperror ('cuadjtq (4) ')
#endif

    END IF

  END IF

  IF (kcall==0) THEN

    isum = 0
!DIR$ IVDEP
!OCL NOVREC
    DO jl = 1, klon
      zqp(jl) = 1./pp(jl)
#ifndef NOLOOKUP
      it = INT(pt(jl,kk)*1000.)
      IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
#ifdef DEBUG
      IF (it<jptlucu1 .OR. it>jptlucu2) THEN
        WRITE(*,*) 'cuadjtq: 5 it=',it,jl,kk
      END IF
#endif
      it = MAX(MIN(it,jptlucu2),jptlucu1)
      zqsat = tlucua(it)*zqp(jl)
      zqsat = MIN(0.5,zqsat)
      zcor = 1./(1.-vtmpc1*zqsat)
      zqsat = zqsat*zcor
      zcond(jl) = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor*tlucub(it))
      pt(jl,kk) = pt(jl,kk) + tlucuc(it)*zcond(jl)
#else
      lo = pt(jl,kk) > tmelt
      zqsat = c2es*EXP(MERGE(c3les,c3ies,lo)*(pt(jl,kk)-tmelt) &
         /      (pt(jl,kk)-MERGE(c4les,c4ies,lo))) *zqp(jl)
      zqsat = MIN(0.5,zqsat)
      zcor = 1./(1.-vtmpc1*zqsat)
      zqsat = zqsat*zcor
      zcond(jl) = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor &
      * (MERGE(c5alvcp,c5alscp,lo) / (pt(jl,kk)-MERGE(c4les,c4ies,lo))**2))
      pt(jl,kk) = pt(jl,kk) + MERGE(alvdcp, alsdcp, lo)*zcond(jl)
#endif
      pq(jl,kk) = pq(jl,kk) - zcond(jl)
      IF (ABS(zcond(jl)) > 0.) isum = isum + 1
    END DO
#ifndef NOLOOKUP
    IF (lookupoverflow) CALL lookuperror ('cuadjtq (5) ')
#endif

    IF (isum/=0) THEN
!DIR$ IVDEP
!OCL NOVREC
      DO jl = 1, klon
#ifndef NOLOOKUP
        it = INT(pt(jl,kk)*1000.)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
#ifdef DEBUG
        IF (it<jptlucu1 .OR. it>jptlucu2) THEN
          WRITE(*,*) 'cuadjtq: 6 it=',it,jl,kk
        END IF
#endif
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqsat = tlucua(it)*zqp(jl)
        zqsat = MIN(0.5,zqsat)
        zcor = 1./(1.-vtmpc1*zqsat)
        zqsat = zqsat*zcor
        zcond1 = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor*tlucub(it))
        pt(jl,kk) = pt(jl,kk) + tlucuc(it)*zcond1
#else
        lo = pt(jl,kk) > tmelt
        zqsat = c2es*EXP(MERGE(c3les,c3ies,lo)*(pt(jl,kk)-tmelt) &
         /      (pt(jl,kk)-MERGE(c4les,c4ies,lo)))*zqp(jl)
        zqsat = MIN(0.5,zqsat)
        zcor = 1./(1.-vtmpc1*zqsat)
        zqsat = zqsat*zcor
        zcond1 = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor &
         * (MERGE(c5alvcp,c5alscp,lo) / (pt(jl,kk)-MERGE(c4les,c4ies,lo))**2))
        pt(jl,kk) = pt(jl,kk) + MERGE(alvdcp, alsdcp, lo) * zcond1
#endif
        pq(jl,kk) = pq(jl,kk) - zcond1
      END DO
#ifndef NOLOOKUP
      IF (lookupoverflow) CALL lookuperror ('cuadjtq (6) ')
#endif

    END IF

  END IF

  IF (kcall==4) THEN

!DIR$ IVDEP
!OCL NOVREC
    DO jl = 1, klon
      zqp(jl) = 1./pp(jl)
#ifndef NOLOOKUP
      it = INT(pt(jl,kk)*1000.)
      IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
#ifdef DEBUG
      IF (it<jptlucu1 .OR. it>jptlucu2) THEN
        WRITE(*,*) 'cuadjtq: 7 it=',it,jl,kk
      END IF
#endif
      it = MAX(MIN(it,jptlucu2),jptlucu1)
      zqsat = tlucua(it)*zqp(jl)
      zqsat = MIN(0.5,zqsat)
      zcor = 1./(1.-vtmpc1*zqsat)
      zqsat = zqsat*zcor
      zcond(jl) = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor*tlucub(it))
      pt(jl,kk) = pt(jl,kk) + tlucuc(it)*zcond(jl)
#else
      lo = pt(jl,kk) > tmelt
      zqsat = c2es*EXP(MERGE(c3les,c3ies,lo)*(pt(jl,kk)-tmelt) &
         /      (pt(jl,kk)-MERGE(c4les,c4ies,lo)))*zqp(jl)
      zqsat = MIN(0.5,zqsat)
      zcor = 1./(1.-vtmpc1*zqsat)
      zqsat = zqsat*zcor
      zcond(jl) = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor &
       * (MERGE(c5alvcp,c5alscp,lo) / (pt(jl,kk)-MERGE(c4les,c4ies,lo))**2))
      pt(jl,kk) = pt(jl,kk) + MERGE(alvdcp, alsdcp, lo) * zcond(jl)
#endif
      pq(jl,kk) = pq(jl,kk) - zcond(jl)
    END DO
#ifndef NOLOOKUP
    IF (lookupoverflow) CALL lookuperror ('cuadjtq (7) ')
#endif

!DIR$ IVDEP
!OCL NOVREC
    DO jl = 1, klon
#ifndef NOLOOKUP
      it = INT(pt(jl,kk)*1000.)
      IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
#ifdef DEBUG
      IF (it<jptlucu1 .OR. it>jptlucu2) THEN
        WRITE(*,*) 'cuadjtq: 8 it=',it,jl,kk
      END IF
#endif
      it = MAX(MIN(it,jptlucu2),jptlucu1)
      zqsat = tlucua(it)*zqp(jl)
      zqsat = MIN(0.5,zqsat)
      zcor = 1./(1.-vtmpc1*zqsat)
      zqsat = zqsat*zcor
      zcond1 = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor*tlucub(it))
      pt(jl,kk) = pt(jl,kk) + tlucuc(it)*zcond1
#else
      lo = pt(jl,kk) > tmelt
      zqsat = c2es*EXP(MERGE(c3les,c3ies,lo)*(pt(jl,kk)-tmelt) &
         /      (pt(jl,kk)-MERGE(c4les,c4ies,lo)))*zqp(jl)
      zqsat = MIN(0.5,zqsat)
      zcor = 1./(1.-vtmpc1*zqsat)
      zqsat = zqsat*zcor
      zcond1 = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor &
       * (MERGE(c5alvcp,c5alscp,lo) / (pt(jl,kk)-MERGE(c4les,c4ies,lo))**2))
      pt(jl,kk) = pt(jl,kk) + MERGE(alvdcp, alsdcp, lo) * zcond1
#endif
      pq(jl,kk) = pq(jl,kk) - zcond1
    END DO

#ifndef NOLOOKUP
    IF (lookupoverflow) CALL lookuperror ('cuadjtq (8) ')
#endif

  END IF

  RETURN
END SUBROUTINE cuadjtq
