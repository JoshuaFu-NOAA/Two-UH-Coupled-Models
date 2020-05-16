!+ for solar zenith angle and relative daylength.
!+ $Id: m_solang.f90,v 1.4 1999/07/23 08:21:23 m214030 Exp $

!MODULE m_solang
!
!CONTAINS

SUBROUTINE solang(klon,klp2,klof,ptim1,ptim2,ptim3,pamu0,prdayl)

  ! Description:
  !
  ! For solar zenith angle and relative daylength.
  !
  ! Method:
  !
  ! This routine gives different results depending on a logical
  ! switch. If ldiur is true one obtains actual solar zenith angles
  ! and values of one or zero depending on the sign of the former. If
  ! ldiur is false one gets the same answers at all points, i.e. mean
  ! value of the daytime solar zenith angle and relative length of
  ! the day.
  !
  ! *solang* is called from *radmod* and from *radheat*.
  ! There are three dummy arguments: *ptim1*, *ptim2* and *ptim3*
  ! are latitude dependent parameter about the sun's position.
  ! The routine returns solar zenith angles and relative day
  ! lengths to the long term storage.
  !
  ! Staightforward in the case "on". For the case "off" the
  ! type of "on" computation is repeated  with 128 points and the
  ! relevant mean values are computed and stored.
  !
  ! Authors:
  !
  ! J. F. Geleyn, ECMWF, June 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_gaussgrid,    ONLY: coslon, sinlon
  USE mo_rad_switches, ONLY: ldiur
  USE mo_physc1,       ONLY: cosrad, sinrad

  IMPLICIT NONE

  !  Scalar arguments 
  REAL    ,INTENT(in) :: ptim1, ptim2, ptim3
  INTEGER ,INTENT(in) :: klon ! number of longitudes
  INTEGER ,INTENT(in) :: klp2 ! number of longitudes allocated
  INTEGER ,INTENT(in) :: klof ! longitude index offset to sin/coslon

  !  Array arguments 
  REAL ,INTENT(out) :: pamu0(klp2), prdayl(klp2)

  !  Local scalars: 
  REAL :: zs1, zs2, ztim1, ztim2, ztim3
  INTEGER :: jl
  LOGICAL :: lo

  !  Local arrays: 
  REAL :: zmu0(128), zrdayl(128)

  !  Intrinsic functions 
  INTRINSIC MERGE, SUM


  !  Executable statements 

  ztim1 = ptim1
  ztim2 = ptim2
  ztim3 = ptim3

!-- 1. Computations if diurnal cycle "on"

  IF (ldiur) THEN
    DO jl = 1, klon
      pamu0(jl) = ztim1 + ztim2*coslon(jl+klof) + ztim3*sinlon(jl+klof)
      lo = pamu0(jl) >= 0.
      pamu0(jl) = MERGE(pamu0(jl),0.,lo)
      prdayl(jl) = MERGE(1.,0.,lo)
    END DO

!-- 2. Computations if diurnal cycle "off"

  ELSE
    DO jl = 1, 128
      zmu0(jl) = ztim1 + ztim2*cosrad(jl) + ztim3*sinrad(jl)
      lo = zmu0(jl) >= 0.
      zmu0(jl) = MERGE(zmu0(jl),0.,lo)
      zrdayl(jl) = MERGE(1.,0.,lo)
    END DO
    zs1 = SUM(zmu0(1:128))
    zs2 = SUM(zrdayl(1:128))
    IF (ABS(zs2) > 0.) THEN
      zs1 = zs1/zs2
      zs2 = zs2/128.
    END IF
    DO jl = 1, klon
      pamu0(jl) = zs1
      prdayl(jl) = zs2
    END DO
  END IF

  RETURN
END SUBROUTINE solang

!END MODULE m_solang
