!+ interpolates large-scale fields of t,q etc.
!+ $Id: cuini.f90,v 1.4 1998/10/28 12:28:50 m214003 Exp $

SUBROUTINE cuini(klp2,k2lp2,klon,klev,klevp1,klevm1,pten,pqen,pqsen,pxen, &
&      puen,pven,ktrac,pxten,pxtenh,pxtu,pxtd,pmfuxt,pmfdxt,pverv,pgeo,paphp1, &
&      pgeoh,ptenh,pqenh,pqsenh,pxenh,klwmin,ptu,pqu,ptd,pqd,puu,pvu,pud,pvd, &
&      pmfu,pmfd,pmfus,pmfds,pmfuq,pmfdq,pdmfup,pdmfdp,pdpmel,plu,plude,klab)

  ! Description:
  !
  ! Interpolates large-scale fields of t,q etc.
  !
  ! Method:
  !
  ! This routine interpolates large-scale fields of t,q etc.
  ! to half levels (i.e. grid for massflux scheme),
  ! determines level of maximum vertical velocity
  ! and initializes values for updrafts and downdrafts
  !
  ! This routine is called from *cumastr*.
  !
  ! For extrapolation to half levels see tiedtke(1989)
  !
  ! *cuadjtq* to specify qs at half levels
  !
  ! Authors:
  !
  ! M. Tiedtke, ECMWF, December 1989, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_constants, only: cpd, &! specific heat at constant pressure
                          rcpd  ! rcpd=1./cpd

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: k2lp2, klev, klevm1, klevp1, klon, klp2, ktrac

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: paphp1(klp2,klevp1), pgeo(klp2,klev), pqen(klp2,klev),  &
&                      pqsen(klp2,klev), pten(klp2,klev), puen(k2lp2,klev), &
&                      pven(k2lp2,klev), pverv(klp2,klev), pxen(klp2,klev), &
&                      pxten(klon,klev,ktrac)

   !  Array arguments with intent(InOut):
  REAL, INTENT (INOUT) :: pdmfdp(klp2,klev), pdmfup(klp2,klev), &
&                         pdpmel(klp2,klev), pgeoh(klp2,klev), plu(klp2,klev), &
&                         plude(klp2,klev), pmfd(klp2,klev), pmfdq(klp2,klev), &
&                         pmfds(klp2,klev), pmfdxt(klon,klev,ktrac), &
&                         pmfu(klp2,klev), pmfuq(klp2,klev), pmfus(klp2,klev), &
&                         pmfuxt(klon,klev,ktrac), pqd(klp2,klev), &
&                         pqenh(klp2,klev), pqsenh(klp2,klev), pqu(klp2,klev), &
&                         ptd(klp2,klev), ptenh(klp2,klev), ptu(klp2,klev), &
&                         pud(klp2,klev), puu(klp2,klev), pvd(klp2,klev), &
&                         pvu(klp2,klev), pxenh(klp2,klev), pxtd(klon,klev,ktrac), &
&                         pxtenh(klon,klev,ktrac), pxtu(klon,klev,ktrac)
  INTEGER, INTENT (INOUT) :: klab(klp2,klev), klwmin(klp2)

  !  Local scalars: 
  REAL :: zdp, zzs
  INTEGER :: icall, ik, jk, jl, jt

  !  Local arrays: 
  REAL :: zph(klp2), zwmax(klp2)
  LOGICAL :: loflag(klp2)

  !  External subroutines 
  EXTERNAL cuadjtq

  !  Intrinsic functions 
  INTRINSIC MAX, MIN


  !  Executable statements 

!-- 1. Specify large scale parameters at half levels adjust temperature
!      fields if staticly unstable find level of maximum vertical velocity

  zdp = 0.5
  DO jk = 2, klev
    DO jl = 1, klon
      pgeoh(jl,jk) = pgeo(jl,jk) + (pgeo(jl,jk-1)-pgeo(jl,jk))*zdp
      ptenh(jl,jk) = (MAX(cpd*pten(jl,jk-1)+pgeo(jl,jk-1),cpd*pten(jl, &
&          jk)+pgeo(jl,jk))-pgeoh(jl,jk))*rcpd
      pqsenh(jl,jk) = pqsen(jl,jk-1)
      zph(jl) = paphp1(jl,jk)
      loflag(jl) = .TRUE.
    END DO

    DO jt = 1, ktrac
      DO jl = 1, klon
        pxtenh(jl,jk,jt) = (pxten(jl,jk,jt)+pxten(jl,jk-1,jt))*zdp
      END DO
    END DO

    ik = jk
    icall = 0
    CALL cuadjtq(klp2,klon,klev,ik,zph,ptenh,pqsenh,loflag,icall)

    DO jl = 1, klon
      pxenh(jl,jk) = (pxen(jl,jk)+pxen(jl,jk-1))*zdp
      pqenh(jl,jk) = MIN(pqen(jl,jk-1),pqsen(jl,jk-1)) + &
&          (pqsenh(jl,jk)-pqsen(jl,jk-1))
      pqenh(jl,jk) = MAX(pqenh(jl,jk),0.)
    END DO
  END DO

  DO jl = 1, klon
    ptenh(jl,klev) = (cpd*pten(jl,klev)+pgeo(jl,klev)-pgeoh(jl,klev))*rcpd
    pxenh(jl,klev) = pxen(jl,klev)
    pqenh(jl,klev) = pqen(jl,klev)
    ptenh(jl,1) = pten(jl,1)
    pxenh(jl,1) = pxen(jl,1)
    pqenh(jl,1) = pqen(jl,1)
    pgeoh(jl,1) = pgeo(jl,1)
    klwmin(jl) = klev
    zwmax(jl) = 0.
  END DO

  DO jt = 1, ktrac
    DO jl = 1, klon
      pxtenh(jl,klev,jt) = pxten(jl,klev,jt)
      pxtenh(jl,1,jt) = pxten(jl,1,jt)
    END DO
  END DO

  DO jk = klevm1, 2, -1
    DO jl = 1, klon
      zzs = MAX(cpd*ptenh(jl,jk)+pgeoh(jl,jk),cpd*ptenh(jl,jk+1)+pgeoh(jl,jk+ &
&          1))
      ptenh(jl,jk) = (zzs-pgeoh(jl,jk))*rcpd
    END DO
  END DO

  DO jk = klev, 3, -1
!DIR$ IVDEP
!OCL NOVREC
    DO jl = 1, klon
      IF (pverv(jl,jk)<zwmax(jl)) THEN
        zwmax(jl) = pverv(jl,jk)
        klwmin(jl) = jk
      END IF
    END DO
  END DO

!-- 2. Initialize values for updrafts and downdrafts

  DO jk = 1, klev
    ik = jk - 1
    IF (jk==1) ik = 1
    DO jl = 1, klon
      ptu(jl,jk) = ptenh(jl,jk)
      ptd(jl,jk) = ptenh(jl,jk)
      pqu(jl,jk) = pqenh(jl,jk)
      pqd(jl,jk) = pqenh(jl,jk)
      plu(jl,jk) = 0.
      puu(jl,jk) = puen(jl,ik)
      pud(jl,jk) = puen(jl,ik)
      pvu(jl,jk) = pven(jl,ik)
      pvd(jl,jk) = pven(jl,ik)
      pmfu(jl,jk) = 0.
      pmfd(jl,jk) = 0.
      pmfus(jl,jk) = 0.
      pmfds(jl,jk) = 0.
      pmfuq(jl,jk) = 0.
      pmfdq(jl,jk) = 0.
      pdmfup(jl,jk) = 0.
      pdmfdp(jl,jk) = 0.
      pdpmel(jl,jk) = 0.
      plude(jl,jk) = 0.
      klab(jl,jk) = 0
    END DO

    DO jt = 1, ktrac
      DO jl = 1, klon
        pxtu(jl,jk,jt) = pxtenh(jl,jk,jt)
        pxtd(jl,jk,jt) = pxtenh(jl,jk,jt)
        pmfuxt(jl,jk,jt) = 0.
        pmfdxt(jl,jk,jt) = 0.
      END DO

    END DO
  END DO

  RETURN
END SUBROUTINE cuini
