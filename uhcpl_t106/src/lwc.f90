!+ longwave radiation, cloud effects

SUBROUTINE lwc(kdlon,kflev,pbint,pbsuin,pcldlw,pcntrb,pfluc,pflux)

  ! Description:
  !
  ! Introduces cloud effects on longwave fluxes or radiances.
  !
  ! Method:
  !
  ! Explicit arguments :
  ! ==== inputs ===
  ! pbint  : (kdlon,0:kflev)        ; half level Planck function
  ! pbsuin : (kdlon)                ; surface Planck function
  ! pcldlw : (kdlon,kflev)          ; cloud fractional cover
  ! pcntrb : (kdlon,0:kflev,0:kflev); clear-sky energy exchange
  ! pcts   : (kdlon,kflev)          ; clear-sky layer cooling-to-space
  ! pfluc  : (kdlon,2,kflev+1)      ; clear-sky fluxes
  ! ==== outputs ===
  ! pflux(kdlon,2,kflev)            ; radiative fluxes :
  !                      1  ==>  upward   flux total
  !                      2  ==>  downward flux total
  !
  ! 1. Initializes all fluxes to clear-sky values
  ! 2. Effect of one overcast unity emissivity cloud layer
  ! 3. Effect of semi-transparent, partial or multi-layered clouds
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

  USE mo_radiation
  USE mo_radint

  IMPLICIT NONE

  !  Scalar arguments 
  INTEGER :: kdlon, kflev

  !  Array arguments 
  REAL :: pbint(kdlon,kflev+1), pbsuin(kdlon), pcldlw(kdlon,kflev), &
          pcntrb(kdlon,kflev+1,kflev+1), pfluc(kdlon,2,kflev+1),    &
          pflux(kdlon,2,kflev+1)

  !  Local scalars: 
  REAL :: zcfrac
  INTEGER :: imaxc, jcloud, jk, jk1, jk2, jkc, jkcp1, jkj, jkm1, jkp1, jl

  !  Local arrays: 
  REAL :: zclear(kdlon), zclm(kdlon,kflev+1,kflev+1), zcloud(kdlon), &
          zdnf(kdlon,kflev+1,kflev+1), zfd(kdlon), zfu(kdlon),       &
          zupf(kdlon,kflev+1,kflev+1)
  INTEGER :: imx(kdlon), imxp(kdlon)

  !  Intrinsic functions 
  INTRINSIC MAX, MIN


  !  Executable statements 

 !-- 1. Initialization

  imaxc = 0

  DO jl = 1, kdlon
    imx(jl) = 0
    imxp(jl) = 0
    zcloud(jl) = 0.
  END DO

!-- 1.1 Search the layer index of the highest cloud

  DO jk = 1, kflev
    DO jl = 1, kdlon
      IF (pcldlw(jl,jk) > zepsc) THEN
         imxp(jl) = jk
      ELSE
         imxp(jl) = imx(jl)
      END IF  
    END DO
!    imaxc = MAX(MAXVAL(imxp(:)),imaxc)
    imx(:) = imxp(:)
  END DO

  ! this change is necessary for keeping consistent results on different
  ! numbers of PEs and follows the changes in the parallel IFS code. 
  imaxc = kflev
 
  pflux(:,:,:) = pfluc(:,:,:)

!-- 2. Effect of cloudiness on longwave fluxes

  IF (imaxc > 0) THEN

!-- 2.0 Initialize to clear-sky fluxes

    DO jk1 = 1, kflev + 1
      DO jk2 = 1, kflev + 1
        DO jl = 1, kdlon
          zupf(jl,jk2,jk1) = pfluc(jl,1,jk1)
          zdnf(jl,jk2,jk1) = pfluc(jl,2,jk1)
        END DO
      END DO
    END DO

!-- 2.1 Fluxes for one overcast unity emissivity cloud

    DO jkc = 1, imaxc
      jcloud = jkc

!-- 2.1.1 Above the cloud

      jkcp1 = jcloud + 1

      DO jk = jkcp1, kflev + 1
        jkm1 = jk - 1
        DO jl = 1, kdlon
          zfu(jl) = 0.
        END DO
        IF (jk>jkcp1) THEN
          DO jkj = jkcp1, jkm1
            DO jl = 1, kdlon
              zfu(jl) = zfu(jl) + pcntrb(jl,jk,jkj)
            END DO
          END DO
        END IF

        DO jl = 1, kdlon
          zupf(jl,jkcp1,jk) = pbint(jl,jk) - zfu(jl)
        END DO
      END DO

!-- 2.1.2 Below the cloud

      DO jk = 1, jcloud
        jkp1 = jk + 1
        DO jl = 1, kdlon
          zfd(jl) = 0.
        END DO

        IF (jk<jcloud) THEN
          DO jkj = jkp1, jcloud
            DO jl = 1, kdlon
              zfd(jl) = zfd(jl) + pcntrb(jl,jk,jkj)
            END DO
          END DO
        END IF
        DO jl = 1, kdlon
          zdnf(jl,jkcp1,jk) = -pbint(jl,jk) - zfd(jl)
        END DO

      END DO
    END DO

!-- 2.2 Cloud cover matrix

!   zclm(jk1,jk2) is the obscuration factor by cloud layers between
!   half-levels jk1 and jk2 as seen from jk1

    DO jk1 = 1, kflev + 1
      DO jk2 = 1, kflev + 1
        DO jl = 1, kdlon
          zclm(jl,jk1,jk2) = 0.
        END DO
      END DO
    END DO

!-- 2.3 Cloud cover below the level of calculation

    DO jk1 = 2, kflev + 1
      DO jl = 1, kdlon
        zclear(jl) = 1.
        zcloud(jl) = 0.
      END DO
      DO jk = jk1 - 1, 1, -1
        DO jl = 1, kdlon
          zclear(jl) = zclear(jl)*(1.0-MAX(pcldlw(jl,jk),zcloud(jl))) &
                     / (1.0-MIN(zcloud(jl),1.-zepsec))
          zclm(jl,jk1,jk) = 1.0 - zclear(jl)
          zcloud(jl) = pcldlw(jl,jk)
        END DO
      END DO
    END DO

!-- 2.4 Cloud cover above the level of calculation

    DO jk1 = 1, kflev
      DO jl = 1, kdlon
        zclear(jl) = 1.
        zcloud(jl) = 0.
      END DO
      DO jk = jk1, kflev
        DO jl = 1, kdlon
          zclear(jl) = zclear(jl)*(1.0-MAX(pcldlw(jl,jk),zcloud(jl))) &
                     / (1.0-MIN(zcloud(jl),1.-zepsec))
          zclm(jl,jk1,jk) = 1.0 - zclear(jl)
          zcloud(jl) = pcldlw(jl,jk)
        END DO
      END DO
    END DO

!-- 3. Fluxes for partial/multiple layered cloudiness

!-- 3.1 Downward fluxes

    DO jl = 1, kdlon
      pflux(jl,2,kflev+1) = 0.
    END DO

    DO jk1 = kflev, 1, -1

      ! Contribution from clear-sky fraction

      DO jl = 1, kdlon
        zfd(jl) = (1.-zclm(jl,jk1,kflev))*zdnf(jl,1,jk1)
      END DO

      ! Contribution from adjacent cloud

      DO jl = 1, kdlon
        zfd(jl) = zfd(jl) + zclm(jl,jk1,jk1)*zdnf(jl,jk1+1,jk1)
      END DO

      ! Contribution from other cloudy fractions

      DO jk = kflev - 1, jk1, -1
        DO jl = 1, kdlon
          zcfrac = zclm(jl,jk1,jk+1) - zclm(jl,jk1,jk)
          zfd(jl) = zfd(jl) + zcfrac*zdnf(jl,jk+2,jk1)
        END DO
      END DO

      DO jl = 1, kdlon
        pflux(jl,2,jk1) = zfd(jl)
      END DO
    END DO

!-- 3.2 Upward flux at the surface

    DO jl = 1, kdlon
      pflux(jl,1,1) = cemiss*pbsuin(jl) - (1.-cemiss)*pflux(jl,2,1)
    END DO

!-- 3.3 Upward fluxes

    DO jk1 = 2, kflev + 1

      ! Contribution from clear-sky fraction

      DO jl = 1, kdlon
        zfu(jl) = (1.-zclm(jl,jk1,1))*zupf(jl,1,jk1)
      END DO

      ! Contribution from adjacent cloud

      DO jl = 1, kdlon
        zfu(jl) = zfu(jl) + zclm(jl,jk1,jk1-1)*zupf(jl,jk1,jk1)
      END DO

      ! Contribution from other cloudy fractions

      DO jk = 2, jk1 - 1
        DO jl = 1, kdlon
          zcfrac  = zclm(jl,jk1,jk-1) - zclm(jl,jk1,jk)
          zfu(jl) = zfu(jl) + zcfrac*zupf(jl,jk,jk1)
        END DO
      END DO

      DO jl = 1, kdlon
        pflux(jl,1,jk1) = zfu(jl)
      END DO
    END DO

  END IF

!-- 3.4 End of cloud effect computations

END SUBROUTINE lwc
