!+ longwave radiation, vertical integration
!+ $Id: lwv.f90,v 1.5 1998/12/11 13:23:37 m214003 Exp $

SUBROUTINE lwv(kdlon,kflev,kuaer,ktraer,kaer,kcfc,pabcu,pbint,pbtop,pdbsl, &
&      pga,pgb,pgasur,pgbsur,pgatop,pgbtop,pgc,pgd,pgcsur,pgdsur,pgctop, &
&      pgdtop,pcntrb,pfluc)

  ! Description:
  !
  ! Carries out the vertical integration to give longwave
  ! fluxes or radiances.
  !
  ! Method:
  !
  ! Explicit arguments :
  ! ==== inputs ===
  ! kx.    (kdlon,.                 ; temperature indices
  ! pabcu : (kdlon,nua,3*kflev+1)   ; absorber amounts
  ! pbint  : (kdlon,0:kflev)        ; half-level planck function
  ! pbtop  : (kdlon,nsint)          ; t.o.a. spectral planck function
  ! pdbsl  : (kdlon,kflev*2)        ; sub-layer planck function gradient
  ! ptave  : (kdlon,kflev)          ; temperature
  ! ==== outputs ===
  ! pcntrb : (kdlon,kflev+1,kflev+1); clear-sky energy exchange matrix
  ! if kmode = 0, 1, 2
  ! pfluc(kdlon,2,kflev)            ; radiative fluxes clear-sky:
  !                      1  ==>  upward   flux total
  !
  ! 1. Performs the vertical integration distinguishing between
  !    contributions by -  the nearby layers
  !                     -  the distant layers
  !                     -  the boundary terms
  ! 2. Computes the clear-sky downward and upward emissivities.
  !
  ! Reference:
  ! See radiation's part of the model's documentation and
  ! ECMWF research department documentation of the "in core model"
  !
  ! Authors:
  !
  ! J.-J. Morcrette, ECMWF, July 1989, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_longwave

  IMPLICIT NONE

  !  Scalar arguments 
  INTEGER :: kaer, kcfc, kdlon, kflev, ktraer, kuaer

  !  Array arguments 
  REAL :: pabcu(kdlon,nua,3*kflev+1), pbint(kdlon,kflev+1), &
&      pbtop(kdlon,nsint), pcntrb(kdlon,kflev+1,kflev+1), &
&      pdbsl(kdlon,nsint,kflev*2), pfluc(kdlon,2,kflev+1), &
&      pga(kdlon,8,2,kflev), pgasur(kdlon,8,2), pgatop(kdlon,8,2), &
&      pgb(kdlon,8,2,kflev), pgbsur(kdlon,8,2), pgbtop(kdlon,8,2), &
&      pgc(kdlon,5,2,kflev), pgcsur(kdlon,5,2), pgctop(kdlon,5,2), &
&      pgd(kdlon,5,2,kflev), pgdsur(kdlon,5,2), pgdtop(kdlon,5,2)

  !  Local arrays: 
  REAL :: zadjd(kdlon,kflev+1), zadju(kdlon,kflev+1), &
&      zdbdt(kdlon,nsint,kflev), zdisd(kdlon,kflev+1), zdisu(kdlon,kflev+1), &
&      zdwfsu(kdlon,nsint)

  !  External subroutines 
  EXTERNAL lwvb, lwvd, lwvn


  !  Executable statements 

!DIR$ NOBOUNDS

!-- 1. Initialization

!-- 1.1 Initialize layer contributions

  zadjd(:,:) = 0.
  zadju(:,:) = 0.
  zdisd(:,:) = 0.
  zdisu(:,:) = 0.

  zdwfsu(:,:) = 0.

!-- 2. Vertical integration

!-- 2.1 Contribution from adjacent layers

  CALL lwvn(kdlon,kflev,kuaer,kaer,kcfc,pabcu,pdbsl,pga,pgb,pgc,pgd,zadjd, &
&      zadju,pcntrb,zdbdt,zdwfsu)

!-- 2.2 Contribution from distant layers

  CALL lwvd(kdlon,kflev,kuaer,ktraer,kaer,kcfc,pabcu,zdbdt,pga,pgb,pgc,pgd, &
&      pcntrb,zdisd,zdisu,zdwfsu)

!-- 2.3 Exchange with the boundaries

  CALL lwvb(kdlon,kflev,kuaer,kaer,kcfc,pabcu,zadjd,zadju,pbint,pbtop,zdisd, &
&      zdisu,pgasur,pgbsur,pgatop,pgbtop,pgcsur,pgdsur,pgctop,pgdtop,pfluc, &
&      zdwfsu)

  RETURN
END SUBROUTINE lwv
