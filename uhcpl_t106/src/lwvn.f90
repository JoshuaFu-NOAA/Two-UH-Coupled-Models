!+ l.w., vertical integration, nearby layers
!+ $Id: lwvn.f90,v 1.6 1999/07/19 19:07:10 m214030 Exp $

SUBROUTINE lwvn(kdlon,kflev,kuaer,kaer,kcfc,pabcu,pdbsl,pga,pgb,pgc,pgd, &
&      padjd,padju,pcntrb,pdbdt,pdwfsu)

  ! Description:
  !
  ! Carries out the vertical integration on nearby layers
  ! to give longwave fluxes or radiances.
  !
  ! Method:
  !
  ! Explicit arguments :
  ! ==== inputs ===
  ! pabcu  : (kdlon,nua,3*kflev+1)  ; absorber amounts
  ! pdbsl  : (kdlon,kflev*2)        ; sub-layer planck function gradient
  ! ==== outputs ===
  ! padj   : (kdlon,kflev+1)        ; contribution of adjacent layers
  ! pcntrb : (kdlon,kflev+1,kflev+1); clear-sky energy exchange matrix
  ! pdbdt  : (kdlon,nua,kflev)      ; layer planck function gradient
  ! pdwfsu : (kdlon,nsint)          ; downward band flux at surface
  !
  ! 1. Performs the vertical integration corresponding to the
  !    contributions of the adjacent layers using a gaussian quadrature
  !
  ! Reference
  ! See radiation's part of the model's documentation and
  ! ECMWF research department documentation of the "in core model"
  !
  ! Authors:
  !
  ! J.-J. Morcrette, ECMWF, July 1989, original source
  ! U. Schlese, DKRZ, July 1993, changed
  ! M. Giorgetta, MPI, July 1993, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_longwave

  IMPLICIT NONE

  !  Scalar arguments 
  INTEGER :: kaer, kcfc, kdlon, kflev, kuaer

  !  Array arguments 
  REAL :: pabcu(kdlon,nua,3*kflev+1), padjd(kdlon,kflev+1), &
&      padju(kdlon,kflev+1), pcntrb(kdlon,kflev+1,kflev+1), &
&      pdbdt(kdlon,nsint,kflev), pdbsl(kdlon,nsint,kflev*2), &
&      pdwfsu(kdlon,nsint), pga(kdlon,8,2,kflev), pgb(kdlon,8,2,kflev), &
&      pgc(kdlon,5,2,kflev), pgd(kdlon,5,2,kflev)

  !  Local scalars: 
  REAL :: zwtr
  INTEGER :: ig, jk, jk1, jk2, jl, kbs, kdd, km12, kmu, knd, knu, kxd, kxu

  !  Local arrays: 
  REAL :: zglayd(kdlon), zglayu(kdlon), ztt(kdlon,ntra), zttnc(kdlon,5), &
&      zuu(kdlon,nua)

  !  External subroutines 
  EXTERNAL lwtt

!DIR$ NOBOUNDS

  !  Executable statements 

!-- 1. Initialization

!-- 1.1 Initialize layer contributions

  padjd(:,:) = 0.
  padju(:,:) = 0.

!-- 2. Vertical integration

!-- 2.1 Contribution from adjacent layers

  DO jk = 1, kflev

!-- 2.1.1 Downward layers

    km12 = 2*(jk-1)
    knd = (jk-1)*ng1p1 + 1
    kxd = knd
    knu = jk*ng1p1 + 1
    kxu = knd

    DO jl = 1, kdlon
      zglayd(jl) = 0.
      zglayu(jl) = 0.
    END DO

    DO ig = 1, ng1
      kbs = km12 + ig
      kdd = kxd + ig
!      DO jl = 1, kdlon*kuaer
!        zuu(jl,1) = pabcu(jl,1,knd) - pabcu(jl,1,kdd)
!      END DO
      zuu = pabcu(:,:,knd) - pabcu(:,:,kdd)

      CALL lwtt(kdlon,kaer,kcfc,pga(1,1,1,jk),pgb(1,1,1,jk),pgc(1,1,1,jk), &
&          pgd(1,1,1,jk),zuu,ztt,zttnc)

      IF (jk==1) THEN
        DO jl = 1, kdlon
          pdwfsu(jl,1) = pdwfsu(jl,1) + wg1(ig)*pdbsl(jl,1,kbs)*ztt(jl,1)*ztt &
&              (jl,10)
          pdwfsu(jl,2) = pdwfsu(jl,2) + wg1(ig)*pdbsl(jl,2,kbs)*ztt(jl,2)*ztt &
&              (jl,7)*ztt(jl,11)
          pdwfsu(jl,3) = pdwfsu(jl,3) + wg1(ig)*pdbsl(jl,3,kbs)*ztt(jl,4)*ztt &
&              (jl,8)*ztt(jl,12)
          pdwfsu(jl,4) = pdwfsu(jl,4) + wg1(ig)*pdbsl(jl,4,kbs)*ztt(jl,5)*ztt &
&              (jl,9)*ztt(jl,13)
          pdwfsu(jl,5) = pdwfsu(jl,5) + wg1(ig)*pdbsl(jl,5,kbs)*ztt(jl,3)*ztt &
&              (jl,14)
          pdwfsu(jl,6) = pdwfsu(jl,6) + wg1(ig)*pdbsl(jl,6,kbs)*ztt(jl,6)*ztt &
&              (jl,15)
        END DO
      END IF

      DO jl = 1, kdlon
        zwtr = pdbsl(jl,1,kbs)*ztt(jl,1)*ztt(jl,10) + &
&            pdbsl(jl,2,kbs)*ztt(jl,2)*ztt(jl,7)*ztt(jl,11) + &
&            pdbsl(jl,3,kbs)*ztt(jl,4)*ztt(jl,8)*ztt(jl,12) + &
&            pdbsl(jl,4,kbs)*ztt(jl,5)*ztt(jl,9)*ztt(jl,13) + &
&            pdbsl(jl,5,kbs)*ztt(jl,3)*ztt(jl,14) + pdbsl(jl,6,kbs)*ztt(jl,6)* &
&            ztt(jl,15)
        zglayd(jl) = zglayd(jl) + zwtr*wg1(ig)
      END DO

!-- 2.1.2 Upward layers

      kmu = kxu + ig
!      DO jl = 1, kdlon*kuaer
!        zuu(jl,1) = pabcu(jl,1,kmu) - pabcu(jl,1,knu)
!      END DO
      zuu(:,:) = pabcu(:,:,kmu) - pabcu(:,:,knu)

      CALL lwtt(kdlon,kaer,kcfc,pga(1,1,1,jk),pgb(1,1,1,jk),pgc(1,1,1,jk), &
&          pgd(1,1,1,jk),zuu,ztt,zttnc)

      DO jl = 1, kdlon
        zwtr = pdbsl(jl,1,kbs)*ztt(jl,1)*ztt(jl,10) + &
&            pdbsl(jl,2,kbs)*ztt(jl,2)*ztt(jl,7)*ztt(jl,11) + &
&            pdbsl(jl,3,kbs)*ztt(jl,4)*ztt(jl,8)*ztt(jl,12) + &
&            pdbsl(jl,4,kbs)*ztt(jl,5)*ztt(jl,9)*ztt(jl,13) + &
&            pdbsl(jl,5,kbs)*ztt(jl,3)*ztt(jl,14) + pdbsl(jl,6,kbs)*ztt(jl,6)* &
&            ztt(jl,15)
        zglayu(jl) = zglayu(jl) + zwtr*wg1(ig)

      END DO
    END DO

    DO jl = 1, kdlon
      padjd(jl,jk) = zglayd(jl)
      pcntrb(jl,jk,jk+1) = zglayd(jl)
      padju(jl,jk+1) = zglayu(jl)
      pcntrb(jl,jk+1,jk) = zglayu(jl)
      pcntrb(jl,jk,jk) = 0.0

    END DO
  END DO

  DO jk = 1, kflev
    jk2 = 2*jk
    jk1 = jk2 - 1
!    DO jl = 1, kdlon*nsint
!      pdbdt(jl,1,jk) = pdbsl(jl,1,jk1) + pdbsl(jl,1,jk2)
!    END DO
    pdbdt(:,:nsint,jk) = pdbsl(:,:nsint,jk1) + pdbsl(:,:nsint,jk2)
  END DO

  RETURN
END SUBROUTINE lwvn
