!+ compute diagnostics of physical variables.

SUBROUTINE statp

  ! Description:
  !
  ! This subroutine computes some diagnostics for the physics.
  !
  ! Method:
  !
  ! *statp* is called from *gpc*.
  !
  ! Results:
  ! The diagnostics are computed from information stored
  ! in *mo_diagnostics*.
  !
  ! Authors:
  !
  ! J. F. Geleyn, ECMWF, September 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_parameters,    ONLY: jpnlp2
  USE mo_control,       ONLY: nrow, lamip2
  USE mo_gaussgrid,     ONLY: budw
  USE mo_stat_zonal,    ONLY: gsnz, gtdz, gtsz, gwdz, gwsz
  USE mo_memory_g3a,    ONLY: slmm, td3m, td4m, wsm, snm, wlm
  USE mo_decomposition, ONLY: dc => local_decomposition
  USE mo_physc2,        ONLY: csncri

  IMPLICIT NONE

  !  Local array bounds
  INTEGER :: nglon
  !  Local scalars: 
  REAL :: zcsno, zcts, zdel1, zdel2, zdel3, zdel4, zdel5, zrosn, zsn, zsncri, &
          ztd, zts, zw, zwd, zws
  INTEGER :: irow, jrow, jl

  !  Local arrays: 
  REAL :: zcdtdm(jpnlp2), zcdtsm(jpnlp2)

  !  Intrinsic functions 
  INTRINSIC DOT_PRODUCT


  !  Executable statements 

  !  Local array bounds
  nglon = dc% nglon ! number of longitudes

  irow = nrow(1)
  jrow = nrow(2)

!-- 1. Accumulate diagnostics

  zw = budw(irow)
  zdel1 = 0.065
  zdel2 = 0.254
  zdel3 = 0.913
  zdel4 = 2.902
  zdel5 = 5.700

  zcts  = 0.
  zrosn = 1000./300.
  zcsno = 634500.

  IF (lamip2) THEN
    zsncri = csncri
  ELSE
    zsncri = 0.025
  END IF

  DO jl = 1, nglon
    zcdtsm(jl) = zdel1*zcts*td3m(jl,jrow)
    zcdtdm(jl) = zdel2*zcts*td4m(jl,jrow)
  END DO

  zts = zw*DOT_PRODUCT(slmm(1:nglon,jrow),zcdtsm(1:nglon))
  ztd = zw*DOT_PRODUCT(slmm(1:nglon,jrow),zcdtdm(1:nglon))
  zws = zw*DOT_PRODUCT(slmm(1:nglon,jrow),wsm(1:nglon,jrow)) +  &
        zw*DOT_PRODUCT(slmm(1:nglon,jrow),wlm(1:nglon,jrow))
  zwd = 0.
  zsn = zw*DOT_PRODUCT(slmm(1:nglon,jrow),snm(1:nglon,jrow))

  gtsz(irow) = zts
  gtdz(irow) = ztd
  gwsz(irow) = zws
  gwdz(irow) = zwd
  gsnz(irow) = zsn

  RETURN

END SUBROUTINE statp
