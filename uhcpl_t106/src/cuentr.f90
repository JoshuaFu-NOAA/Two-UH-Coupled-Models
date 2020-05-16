!+ Calculates entrainment/detrainment rates
!+ $Id: cuentr.f90,v 1.6 1998/10/28 12:28:44 m214003 Exp $

SUBROUTINE cuentr(klp2,klon,klev,klevp1,kk,ptenh,paphp1,papp1,klwmin,ldcum, &
&      ktype,kcbot,kctop0,ppbase,pmfu,pentr,podetr,khmin,pgeoh,pdmfen,pdmfde)

  ! Description:
  !
  ! Calculates entrainment/detrainment rates.
  !
  ! Method:
  !
  ! This routine calculates entrainment/detrainment rates
  ! for updrafts in cumulus parameterization
  !
  ! This routine is called from *cuasc*.
  ! Input are environmental values t,q etc and updraft values t,q etc.
  ! It returns entrainment/detrainment rates.
  !
  ! Authors:
  !
  ! S. Tiedtke, ECMWF, in 1989, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_constants, ONLY : g,  & ! gravity acceleration
                           rd    ! gas constant for dry air

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: kk, klev, klevp1, klon, klp2

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: paphp1(klp2,klevp1), papp1(klp2,klev), pentr(klp2), &
&                      pgeoh(klp2,klev), pmfu(klp2,klev), ptenh(klp2,klev)
  INTEGER, INTENT (IN) :: kcbot(klp2), kctop0(klp2), khmin(klp2),  &
&                         klwmin(klp2), ktype(klp2)
  LOGICAL, INTENT (IN) :: ldcum(klp2)

  !  Array arguments with intent(InOut):
  REAL, INTENT (INOUT) :: pdmfde(klp2), pdmfen(klp2), podetr(klp2,klev),  &
&                         ppbase(klp2)

  !  Local scalars: 
  REAL :: arg, zdprho, zentr, zorgde, zpmid, zrg, zrrho, ztmzk, zzmzk
  INTEGER :: ikb, ikh, iklwmin, ikt, jl
  LOGICAL :: llo1, llo2

  !  Intrinsic functions 
  INTRINSIC MAX, MERGE, MIN, TAN


  !  Executable statements 

!-- 1.  Calculate entrainment and detrainment rates

!-- 1.1 Specify entrainment rates for shallow clouds

!-- 1.2 Specify entrainment rates for deep clouds

  zrg = 1./g
  DO jl = 1, klon
    ppbase(jl) = paphp1(jl,kcbot(jl))
    zrrho = (rd*ptenh(jl,kk+1))/paphp1(jl,kk+1)
    zdprho = (paphp1(jl,kk+1)-paphp1(jl,kk))*zrg
    zpmid = 0.5*(ppbase(jl)+paphp1(jl,kctop0(jl)))
    zentr = pentr(jl)*pmfu(jl,kk+1)*zdprho*zrrho
    llo1 = kk < kcbot(jl) .AND. ldcum(jl)
    pdmfde(jl) = MERGE(zentr,0.,llo1)
    llo2 = llo1 .AND. ktype(jl) == 2 .AND. (ppbase(jl)-paphp1(jl,kk)<0.2E5 &
&        .OR. paphp1(jl,kk)>zpmid)
    pdmfen(jl) = MERGE(zentr,0.,llo2)
    iklwmin = MAX(klwmin(jl),kctop0(jl)+2)
    llo2 = llo1 .AND. ktype(jl) == 3 .AND. (kk>=iklwmin .OR. papp1(jl,kk)> &
&        zpmid)
    IF (llo2) pdmfen(jl) = zentr
    llo2 = llo1 .AND. ktype(jl) == 1
    ! Turbulent entrainment
    IF (llo2) pdmfen(jl) = zentr
    ! Organized detrainment, detrainment starts at khmin
    ikb = kcbot(jl)
    podetr(jl,kk) = 0.
    IF (llo2 .AND. kk<=khmin(jl) .AND. kk>=kctop0(jl)) THEN
      ikt = kctop0(jl)
      ikh = khmin(jl)
      IF (ikh>ikt) THEN
        zzmzk = -(pgeoh(jl,ikh)-pgeoh(jl,kk))*zrg
        ztmzk = -(pgeoh(jl,ikh)-pgeoh(jl,ikt))*zrg
        arg = 3.1415*(zzmzk/ztmzk)*0.5
        zorgde = TAN(arg)*3.1415*0.5/ztmzk
        zdprho = (paphp1(jl,kk+1)-paphp1(jl,kk))*(zrg*zrrho)
        podetr(jl,kk) = MIN(zorgde,1.E-3)*pmfu(jl,kk+1)*zdprho
      END IF
    END IF
  END DO

  RETURN
END SUBROUTINE cuentr
