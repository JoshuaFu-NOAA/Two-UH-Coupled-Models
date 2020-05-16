MODULE mo_hyb

  USE mo_parameters

  IMPLICIT NONE

  ! ------------------------------------------------------------------
  !
  ! module *mo_hyb* - *loop indices and surface-pressure independent
  !        variables associated with the vertical finite-difference scheme.
  !
  ! a.j.simmons     e.c.m.w.f.     16/11/81.
  !
  ! ------------------------------------------------------------------

  INTEGER :: nlevm1         !  (number of levels)-1.
  INTEGER :: nplev          !  *number of pressure levels.
  INTEGER :: nplvp1         !  *nplev+1.*
  INTEGER :: nplvp2         !  *nplev+2.*
  INTEGER :: nplvpa         !  *nplvp1,* or 2 if *nplev=0.*
  INTEGER :: nlmsgl         !  *nlev* - (number of sigma levels).
  INTEGER :: nlmslp         !  *nlmsgl+1.*
  INTEGER :: nlmsla         !  *nlmslp,* or 2 if *nlmslp=1.*

  REAL :: apzero            !  *reference pressure for computation of the
  !   hybrid vertical levels.
  REAL :: rpr               !  *reciprocal of reference surface pressure.
  REAL :: rdtr              !   rd*(reference temperature).
  REAL :: apsurf            !   fixed global mean of surface pressure
  REAL :: t0icao            !  *surface temperatur of reference atmosphere
  REAL :: tsticao           !  *stratospheric temperature of reference atmosphere
  REAL :: rdtstic           !  *rd*tsticao
  REAL :: rdlnp0i           !  *rd*ln(surface pressure) of reference atmosphere
  REAL :: alrrdic           !  *lapse-rate parameter of reference atmosphere
  REAL :: rdt0ral           !  *rd*t0icao/alphaic
  REAL :: ptricao           !  *tropopause pressure of reference atmosphere
  REAL :: rdlnpti           !  *rd*ln(ptricao)
  REAL :: gsticao           !  *constant used in geopotential calculation
  REAL, ALLOCATABLE :: ralpha(:)    !   rd*alpha at pressure and sigma levels.
  REAL, ALLOCATABLE :: rlnpr(:)     !   rd*ln(p(k+.5)/p(k-.5)) at pressure and sigma levels.
  REAL, ALLOCATABLE :: dela(:)      !   a(k+.5)-a(k-.5).
  REAL, ALLOCATABLE :: delb(:)      !   b(k+.5)-b(k-.5).
  REAL, ALLOCATABLE :: rddelb(:)    !   rd*delb.
  REAL, ALLOCATABLE :: cpg(:)       !   a(k+.5)*b(k-.5)-b(k+.5)*a(k-.5).
  REAL, ALLOCATABLE :: delpr(:)     !   p(k+.5)-p(k-.5) for reference surface pressure.
  REAL, ALLOCATABLE :: rdelpr(:)    !  *reciprocal of *delpr.*
  REAL, ALLOCATABLE :: ralphr(:)    !  *constant array for use by pgrad.
  REAL, ALLOCATABLE :: alpham(:)    !  *constant array for use by dyn.
  REAL, ALLOCATABLE :: ardprc(:)    !  *constant array for use by dyn.
  REAL, ALLOCATABLE :: rlnmar(:)    !  *constant array for use by pgrad.
  REAL, ALLOCATABLE :: aktlrd(:)    !  *constant array for use by conteq.
  REAL, ALLOCATABLE :: altrcp(:)    !  *constant array for use by conteq.
  REAL, ALLOCATABLE :: ceta(:)      !  *full hybrid vertical levels.
  REAL, ALLOCATABLE :: cetah(:)     !  *half hybrid vertical levels.
  REAL, ALLOCATABLE :: bb(:,:)      !  *gravity wave matrix

CONTAINS

  SUBROUTINE inihyb

    ! Description:
    !
    ! Initializes constants for vertical coordinate calculations.
    !
    ! Method:
    !
    ! Compute loop indices and surface-pressure independent
    ! variables associated with the vertical finite-difference scheme.
    !
    ! Output is in module *mo_hyb*
    !
    ! Authors:
    !
    ! A. J. Simmons, ECMWF, November 1981, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! A. Rhodin, MPI, Jan 1999, subroutine inihyb -> module mo_hyb
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_constants, ONLY: g, rcpd, rd
    USE mo_control,   ONLY: nlev, nlevp1, nvclev, vct, lamip2

    IMPLICIT NONE

    !  Local scalars: 
    REAL :: za, zb, zetam, zetap, zp, zp0icao, zpp, zrd, zs, zsm
    INTEGER :: ilev, ilevp1, iplev, iplvp1, is, ism, ist, jk, jlev

    !  Intrinsic functions 
    INTRINSIC EXP, LOG


    !  Executable statements 

    !-- 1. Initialize variables

    apzero = 101325.
    zrd = rd
    ralpha(1) = zrd*LOG(2.)
    rlnpr(1)  = 2.*ralpha(1)
    ilev   = nlev
    ilevp1 = ilev + 1
    nlevp1 = ilevp1
    nlevm1 = ilev - 1
    iplev  = 0
    iplvp1 = 1
    is  = nvclev + ilevp1
    ism = is - 1
    zpp = vct(1)
    zsm = vct(is)

    IF (lamip2) THEN
      apsurf = 98610.
    ELSE
      apsurf = 98200.
    END IF

    t0icao  = 288.
    tsticao = 216.5
    zp0icao = 101320.
    rdlnp0i = rd*LOG(zp0icao)
    rdtstic = rd*tsticao
    alrrdic = 0.0065/g
    rdt0ral = t0icao/alrrdic
    rdlnpti = rdlnp0i + (LOG(tsticao/t0icao))/alrrdic
    ptricao = EXP(rdlnpti/rd)
    gsticao = tsticao*(rdlnpti-1./alrrdic)

    !-- 2. Calculate pressure-level values

10  CONTINUE

    zb = vct(nvclev+iplvp1+1)
    IF (zb>0.) THEN
       nplev  = iplev
       nplvp1 = iplvp1
       nplvp2 = iplvp1 + 1
       IF (iplev==0) THEN
          nplvpa = 2
       ELSE
          nplvpa = iplvp1
       END IF
       GO TO 20
    ELSE
       iplev  = iplvp1
       iplvp1 = iplev + 1
       IF (iplvp1==ilevp1) GO TO 40
       zp  = zpp
       zpp = vct(iplvp1)
       delpr(iplev) = zpp - zp
       rdelpr(iplev) = 1./delpr(iplev)
       IF (iplev>1) THEN
          rlnpr(iplev)  = zrd*LOG(zpp/zp)
          ralpha(iplev) = zrd - zp*rlnpr(iplev)/delpr(iplev)
       END IF
       alpham(iplev) = ralpha(iplev)*rcpd
       ardprc(iplev) = rlnpr(iplev)*rdelpr(iplev)*rcpd
       GO TO 10
    END IF

    !-- 3. Calculate sigma-level values

20  CONTINUE

    za = vct(ism-nvclev)
    IF (za>0.) THEN
       nlmsgl = ism - nvclev
       nlmslp = nlmsgl + 1
       nlmsla = nlmslp
       GO TO 30
    ELSE
       is  = ism
       ism = is - 1
       ist = is - nvclev
       zs  = zsm
       zsm = vct(is)
       IF (ist==1) THEN
          nlmsgl = 0
          nlmslp = 1
          nlmsla = 2
          GO TO 30
       ELSE
          rlnpr(ist)  = zrd*LOG(zs/zsm)
          ralpha(ist) = zrd - zsm*rlnpr(ist)/(zs-zsm)
       END IF
       GO TO 20
    END IF

    !-- 4. Calculate dela, delb, rddelb, cpg, and complete alphdb

30  CONTINUE

    DO jk = 1, nlev
       dela(jk)   = vct(jk+1) - vct(jk)
       delb(jk)   = vct(nvclev+jk+1) - vct(nvclev+jk)
       rddelb(jk) = rd*delb(jk)
       cpg(jk)    = vct(nvclev+jk)*vct(jk+1) - vct(nvclev+jk+1)*vct(jk)
    END DO

    DO jk = nlmslp, nlev
       alpham(jk) = ralpha(jk)*delb(jk)
    END DO

    !-- 5. Compute full level values of the hybrid coordinate

40  CONTINUE

    zetam = vct(1)/apzero + vct(nvclev+1)
    cetah(1) = zetam

    DO jlev = 1, nlev
       zetap = vct(jlev+1)/apzero + vct(nvclev+1+jlev)
       ceta(jlev) = (zetam+zetap)*.5
       cetah(jlev+1) = zetap
       zetam = zetap
    END DO

    RETURN
  END SUBROUTINE inihyb

END MODULE mo_hyb
