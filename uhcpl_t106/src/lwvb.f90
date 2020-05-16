!+ l.w., vertical integration, exchange with boundaries

SUBROUTINE lwvb(kdlon,kflev,kuaer,kaer,kcfc,pabcu,padjd,padju,pbint,pbtop,           &
                pdisd,pdisu,pgasur,pgbsur,pgatop,pgbtop,pgcsur,pgdsur,pgctop,pgdtop, &
                pfluc,pdwfsu)

  ! Description:
  !
  ! Introduces the effects of the boundaries in the vertical
  ! integration.
  !
  ! Method:
  !
  ! 1. Computes the energy exchange with top and surface of the
  !    atmosphere
  ! 2. Computes the cooling-to-space and heating-from-ground
  !    terms for the approximate cooling rate above 10 hpa
  ! 3. Adds up all contributions to get the clear-sky fluxes
  !
  ! For reference see radiation's part of the model's documentation and
  ! ECMWF research department documentation of the "in core model"
  !
  ! Explicit arguments :
  ! ==== inputs ===
  ! kx.    (kdlon,.               ; temperature indices
  ! pabcu  : (kdlon,nua,3*kflev+1); absorber amounts
  ! padj   : (kdlon,kflev+1)      ; contribution by adjacent layers
  ! pbint  : (kdlon,kflev+1)      ; half-level planck functions
  ! pbsui  : (kdlon)              ; surface planck function
  ! pbsur  : (kdlon,nsint)        ; surface spectral planck function
  ! pbtop  : (kdlon,nsint)        ; t.o.a. spectral planck function
  ! pdis   : (kdlon,kflev+1)      ; contribution by distant layers
  ! ==== outputs ===
  ! if kmode = 0, 1, 2
  ! pfluc(kdlon,2,kflev)          ; radiative fluxes clear-sky:
  !                      1  ==>  upward   flux total
  ! pdwfsu : (kdlon,nsint)        ; downward band flux at surface
  !
  ! Reference:
  ! See radiation's part of the model's documentation and
  ! ECMWF research department documentation of the "in core model"
  !
  ! Authors:
  !
  ! J.-J. Morcrette, ECMWF, July 1989, original source
  ! U. Schlese, DKRZ, unknown, changed
  ! M. Giorgetta, MPI, unknown, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_longwave
  USE mo_radint

  IMPLICIT NONE

  !  Scalar arguments 
  INTEGER :: kaer, kcfc, kdlon, kflev, kuaer

  !  Array arguments 
  REAL :: pabcu(kdlon,nua,3*kflev+1), padjd(kdlon,kflev+1),                &
          padju(kdlon,kflev+1), pbint(kdlon,kflev+1), pbtop(kdlon,nsint),  &
          pdisd(kdlon,kflev+1), pdisu(kdlon,kflev+1), pdwfsu(kdlon,nsint), &
          pfluc(kdlon,2,kflev+1), pgasur(kdlon,8,2), pgatop(kdlon,8,2),    &
          pgbsur(kdlon,8,2), pgbtop(kdlon,8,2), pgcsur(kdlon,5,2),         &
          pgctop(kdlon,5,2), pgdsur(kdlon,5,2), pgdtop(kdlon,5,2)

  !  Local scalars: 
  REAL :: zcnsol, zcntop
  INTEGER :: jk, jl, kn

  !  Local arrays: 
  REAL :: ztt(kdlon,ntra), zttnc(kdlon,5), zuu(kdlon,nua)

  !  External subroutines 
  EXTERNAL lwtt

!DIR$ NOBOUNDS

  !  Executable statements 

!-- 1. Vertical integration

!-- 1.1 Exchange with top of the atmosphere

  DO jk = 1, kflev
    kn = (jk-1)*ng1p1 + 1

!    DO jl = 1, kdlon*kuaer
!      zuu(jl,1) = pabcu(jl,1,kn)
!    END DO
    zuu(:,:) = pabcu(:,:,kn)

    CALL lwtt(kdlon,kaer,kcfc,pgatop(1,1,1),pgbtop(1,1,1),pgctop(1,1,1), &
              pgdtop(1,1,1),zuu,ztt,zttnc)

    IF (jk==1) THEN
      DO jl = 1, kdlon
        pdwfsu(jl,1) = pdwfsu(jl,1) - pbtop(jl,1)*ztt(jl,1)*ztt(jl,10)
        pdwfsu(jl,2) = pdwfsu(jl,2) - pbtop(jl,2)*ztt(jl,2)*ztt(jl,7)*ztt(jl,11)
        pdwfsu(jl,3) = pdwfsu(jl,3) - pbtop(jl,3)*ztt(jl,4)*ztt(jl,8)*ztt(jl,12)
        pdwfsu(jl,4) = pdwfsu(jl,4) - pbtop(jl,4)*ztt(jl,5)*ztt(jl,9)*ztt(jl,13)
        pdwfsu(jl,5) = pdwfsu(jl,5) - pbtop(jl,5)*ztt(jl,3)*ztt(jl,14)
        pdwfsu(jl,6) = pdwfsu(jl,6) - pbtop(jl,6)*ztt(jl,6)*ztt(jl,15)
      END DO
    END IF

    DO jl = 1, kdlon
      zcntop = pbtop(jl,1)*ztt(jl,1)*ztt(jl,10) + pbtop(jl,2)*ztt(jl,2)*ztt(jl,7) &
             *ztt(jl,11) + pbtop(jl,3)*ztt(jl,4)*ztt(jl,8)*ztt(jl,12) +           &
              pbtop(jl,4)*ztt(jl,5)*ztt(jl,9)*ztt(jl,13) +                        &
              pbtop(jl,5)*ztt(jl,3)*ztt(jl,14) + pbtop(jl,6)*ztt(jl,6)*ztt(jl,15)

      pfluc(jl,2,jk) = zcntop - pbint(jl,jk) - pdisd(jl,jk) - padjd(jl,jk)

    END DO
  END DO

  jk = kflev + 1
  kn = (jk-1)*ng1p1 + 1

  DO jl = 1, kdlon
    zcntop = pbtop(jl,1) + pbtop(jl,2) + pbtop(jl,3) + pbtop(jl,4) + &
             pbtop(jl,5) + pbtop(jl,6)
    pfluc(jl,2,jk) = zcntop - pbint(jl,jk) - pdisd(jl,jk) - padjd(jl,jk)
  END DO

!-- 2.2 Exchange with lower limit

  jk = 1
  kn = (jk-1)*ng1p1 + 1

  DO jl = 1, kdlon

    zcnsol = pdwfsu(jl,1) + pdwfsu(jl,2) + pdwfsu(jl,3) + pdwfsu(jl,4) + &
             pdwfsu(jl,5) + pdwfsu(jl,6)
    zcnsol = (1.-cemiss)*zcnsol

    pfluc(jl,1,jk) = zcnsol + pbint(jl,jk) - pdisu(jl,jk) - padju(jl,jk)
  END DO

  DO jk = 2, kflev + 1
    kn = (jk-1)*ng1p1 + 1

!    DO jl = 1, kdlon*kuaer
!      zuu(jl,1) = pabcu(jl,1,1) - pabcu(jl,1,kn)
!    END DO
    zuu(:,:) = pabcu(:,:,1) - pabcu(:,:,kn)

    CALL lwtt(kdlon,kaer,kcfc,pgasur(1,1,1),pgbsur(1,1,1),pgcsur(1,1,1), &
              pgdsur(1,1,1),zuu,ztt,zttnc)

    DO jl = 1, kdlon

      zcnsol = pdwfsu(jl,1)*ztt(jl,1)*ztt(jl,10) + pdwfsu(jl,2)*ztt(jl,2)*          &
               ztt(jl,7)*ztt(jl,11) + pdwfsu(jl,3)*ztt(jl,4)*ztt(jl,8)*ztt(jl,12) + &
               pdwfsu(jl,4)*ztt(jl,5)*ztt(jl,9)*ztt(jl,13) +                        &
               pdwfsu(jl,5)*ztt(jl,3)*ztt(jl,14) + pdwfsu(jl,6)*ztt(jl,6)*          &
               ztt(jl,15)
      zcnsol = (1.-cemiss)*zcnsol

      pfluc(jl,1,jk) = zcnsol + pbint(jl,jk) - pdisu(jl,jk) - padju(jl,jk)

    END DO
  END DO

  RETURN
END SUBROUTINE lwvb
