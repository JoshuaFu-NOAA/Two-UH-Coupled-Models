!+ l.w., vertical integration, distant layers
!+ $Id: lwvd.f90,v 1.8 1999/07/19 19:07:10 m214030 Exp $

SUBROUTINE lwvd(kdlon,kflev,kuaer,ktraer,kaer,kcfc,pabcu,pdbdt,pga,pgb,pgc, &
&      pgd,pcntrb,pdisd,pdisu,pdwfsu)

  ! Description:
  !
  ! Carries out the vertical integration on the distant layers.
  !
  ! Method:
  !
  ! Explicit arguments :
  ! ==== inputs ===
  ! pabcu  : (kdlon,nua,3*kflev+1)  ; absorber amounts
  ! pdbdt  : (kdlon,kflev)          ; layer planck function gradient
  ! ==== outputs ===
  ! pdis   : (kdlon,kflev+1)        ; contribution by distant layers
  ! pcntrb : (kdlon,kflev+1,kflev+1); energy exchange matrix
  ! pdwfsu : (kdlon,nsint)          ; downward band flux at surface
  !
  ! 1. Performs the vertical integration corresponding to the
  !    contributions of the distant layers using trapezoidal rule
  !
  ! Reference:
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
  INTEGER :: kaer, kcfc, kdlon, kflev, ktraer, kuaer

  !  Array arguments 
  REAL :: pabcu(kdlon,nua,3*kflev+1), pcntrb(kdlon,kflev+1,kflev+1), &
&      pdbdt(kdlon,nsint,kflev), pdisd(kdlon,kflev+1), pdisu(kdlon,kflev+1), &
&      pdwfsu(kdlon,nsint), pga(kdlon,8,2,kflev), pgb(kdlon,8,2,kflev), &
&      pgc(kdlon,5,2,kflev), pgd(kdlon,5,2,kflev)

  !  Local scalars: 
  REAL :: zt1, zt10, zt11, zt12, zt13, zt14, zt15, zt2, zt3, zt4, zt5, zt6, &
&      zt7, zt8, zt9, zww
  INTEGER :: jk, jkj, jkl, jkp1, jl, jlk, kd1, kd2, kj, kjp1, km1, kn, &
&      ku1, ku2

  !  Local arrays: 
  REAL :: ztt(kdlon,ntra), ztt1(kdlon,ntra), ztt2(kdlon,ntra), &
&      zttnc(kdlon,5), zuu(kdlon,nua)

  !  External subroutines 
  EXTERNAL lwtt


!DIR$ NOBOUNDS

  !  Executable statements 

!-- 1. Initialization

!-- 1.1 Initialize layer contributions

    pdisd(:,:) = 0.
    pdisu(:,:) = 0.

!-- 2. Vertical integration

!-- 2.1 Contribution from distant layers

!-- 2.1.1 Distant and above layers

!-- 2.1.2 First upper level

  DO jk = 1, kflev - 1
    jkp1 = jk + 1
    kn = (jk-1)*ng1p1 + 1
    kd1 = jk*ng1p1 + 1

!    DO ja = 1, kdlon*kuaer
!      zuu(ja,1) = pabcu(ja,1,kn) - pabcu(ja,1,kd1)
!    END DO
    zuu(:,:) = pabcu(:,:,kn) - pabcu(:,:,kd1)

    CALL lwtt(kdlon,kaer,kcfc,pga(1,1,1,jk),pgb(1,1,1,jk),pgc(1,1,1,jk), &
&        pgd(1,1,1,jk),zuu,ztt1,zttnc)

!-- 2.1.3 Higher up

    DO jkj = jkp1, kflev
      kjp1 = jkj + 1
      kd2 = jkj*ng1p1 + 1

!      DO ja = 1, kdlon*kuaer
!        zuu(ja,1) = pabcu(ja,1,kn) - pabcu(ja,1,kd2)
!      END DO
      zuu(:,:) = pabcu(:,:,kn) - pabcu(:,:,kd2)

      CALL lwtt(kdlon,kaer,kcfc,pga(1,1,1,jkj),pgb(1,1,1,jkj),pgc(1,1,1,jkj), &
&          pgd(1,1,1,jkj),zuu,ztt2,zttnc)

!      DO ja = 1, kdlon*ktraer
!        ztt(ja,1) = (ztt1(ja,1)+ztt2(ja,1))*0.5
!      END DO
      ztt(:,:) = (ztt1(:,:)+ztt2(:,:))*0.5

      IF (jk==1) THEN
        DO jl = 1, kdlon
          pdwfsu(jl,1) = pdwfsu(jl,1) + pdbdt(jl,1,jkj)*ztt(jl,1)*ztt(jl,10)
          pdwfsu(jl,2) = pdwfsu(jl,2) + pdbdt(jl,2,jkj)*ztt(jl,2)*ztt(jl,7)* &
&              ztt(jl,11)
          pdwfsu(jl,3) = pdwfsu(jl,3) + pdbdt(jl,3,jkj)*ztt(jl,4)*ztt(jl,8)* &
&              ztt(jl,12)
          pdwfsu(jl,4) = pdwfsu(jl,4) + pdbdt(jl,4,jkj)*ztt(jl,5)*ztt(jl,9)* &
&              ztt(jl,13)
          pdwfsu(jl,5) = pdwfsu(jl,5) + pdbdt(jl,5,jkj)*ztt(jl,3)*ztt(jl,14)
          pdwfsu(jl,6) = pdwfsu(jl,6) + pdbdt(jl,6,jkj)*ztt(jl,6)*ztt(jl,15)
        END DO
      END IF

      DO jl = 1, kdlon
        zww = pdbdt(jl,1,jkj)*ztt(jl,1)*ztt(jl,10) + &
&             pdbdt(jl,2,jkj)*ztt(jl,2)*ztt(jl,7)*ztt(jl,11) + &
&             pdbdt(jl,3,jkj)*ztt(jl,4)*ztt(jl,8)*ztt(jl,12) + &
&             pdbdt(jl,4,jkj)*ztt(jl,5)*ztt(jl,9)*ztt(jl,13) + &
&             pdbdt(jl,5,jkj)*ztt(jl,3)*ztt(jl,14) + pdbdt(jl,6,jkj)*ztt(jl,6)* &
&             ztt(jl,15)

        pdisd(jl,jk) = pdisd(jl,jk) + zww
        pcntrb(jl,jk,kjp1) = zww
      END DO

!      DO ja = 1, kdlon*ktraer
!        ztt1(ja,1) = ztt2(ja,1)
!      END DO
      ztt1(:,:) = ztt2(:,:)
    END DO
  END DO

!-- 2.1.4 Distant and below layers

!-- 2.1.5 First lower level

  DO jk = 3, kflev + 1
    kn = (jk-1)*ng1p1 + 1
    km1 = jk - 1
    kj = jk - 2
    ku1 = kj*ng1p1 + 1

!    DO ja = 1, kdlon*kuaer
!      zuu(ja,1) = pabcu(ja,1,ku1) - pabcu(ja,1,kn)
!    END DO
    zuu(:,:) = pabcu(:,:,ku1) - pabcu(:,:,kn)

    CALL lwtt(kdlon,kaer,kcfc,pga(1,1,1,km1),pgb(1,1,1,km1),pgc(1,1,1,km1), &
&        pgd(1,1,1,km1),zuu,ztt1,zttnc)

!-- 2.1.6 Down below

    DO jlk = 1, kj
      jkl = km1 - jlk
      ku2 = (jkl-1)*ng1p1 + 1

!      DO ja = 1, kdlon*kuaer
!        zuu(ja,1) = pabcu(ja,1,ku2) - pabcu(ja,1,kn)
!      END DO
      zuu(:,:) = pabcu(:,:,ku2) - pabcu(:,:,kn)

      CALL lwtt(kdlon,kaer,kcfc,pga(1,1,1,jkl),pgb(1,1,1,jkl),pgc(1,1,1,jkl), &
&          pgd(1,1,1,jkl),zuu,ztt2,zttnc)

      DO jl = 1, kdlon
        zt1 = (ztt1(jl,1)+ztt2(jl,1))*0.5
        zt2 = (ztt1(jl,2)+ztt2(jl,2))*0.5
        zt3 = (ztt1(jl,3)+ztt2(jl,3))*0.5
        zt4 = (ztt1(jl,4)+ztt2(jl,4))*0.5
        zt5 = (ztt1(jl,5)+ztt2(jl,5))*0.5
        zt6 = (ztt1(jl,6)+ztt2(jl,6))*0.5
        zt7 = (ztt1(jl,7)+ztt2(jl,7))*0.5
        zt8 = (ztt1(jl,8)+ztt2(jl,8))*0.5
        zt9 = (ztt1(jl,9)+ztt2(jl,9))*0.5
        zt10 = (ztt1(jl,10)+ztt2(jl,10))*0.5
        zt11 = (ztt1(jl,11)+ztt2(jl,11))*0.5
        zt12 = (ztt1(jl,12)+ztt2(jl,12))*0.5
        zt13 = (ztt1(jl,13)+ztt2(jl,13))*0.5
        zt14 = (ztt1(jl,14)+ztt2(jl,14))*0.5
        zt15 = (ztt1(jl,15)+ztt2(jl,15))*0.5

        zww = pdbdt(jl,1,jkl)*zt1*zt10 + pdbdt(jl,2,jkl)*zt2*zt7*zt11 + &
&            pdbdt(jl,3,jkl)*zt4*zt8*zt12 + pdbdt(jl,4,jkl)*zt5*zt9*zt13 + &
&            pdbdt(jl,5,jkl)*zt3*zt14 + pdbdt(jl,6,jkl)*zt6*zt15
        pdisu(jl,jk) = pdisu(jl,jk) + zww
        pcntrb(jl,jk,jkl) = zww
      END DO

!      DO ja = 1, kdlon*ktraer
!        ztt1(ja,1) = ztt2(ja,1)
!      END DO
      ztt1(:,:) = ztt2(:,:)
    END DO
  END DO

  RETURN

END SUBROUTINE lwvd
