!+ contribution of current latitude to global integral
!+ $Id: qmassd.f90,v 1.9 1999/09/08 16:38:42 m214030 Exp $

SUBROUTINE qmassd(cwava,w,q1,q2,pdel,etamid,kftype,hwn,jlat)

  ! Description:
  !
  ! Compute contribution of current latitude to global integral of
  ! F*q2*|q2 - q1|**Beta.
  !
  ! Method:
  !
  ! Compute contribution of current latitude to global integral of
  ! F*q2*|q2 - q1|**Beta.
  ! This is a measure of the difference between the fields before and
  ! after the SLT "forecast".
  ! It is used in the "fixer" which enforces conservation in constituent
  ! fields transport via SLT.
  ! Two options are available:
  ! 1. kftype=1 : F=1.  and Beta=1.5
  ! 2. kftype=2 : F=eta and Beta=1.
  !
  ! Reference Rasch and Williamson, 1991
  !
  ! Authors:
  !
  ! Original version:  J. Olson
  ! Standardized:      J. Rosinski, June 1992
  ! Reviewed:          P. Rasch, D. Williamson, August 1992
  ! Modified:          U. Schlese, April 1995  (optional fixers)
  ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
  ! parallel version:  T. Diehl, DKRZ, July 1999
  !                    A. Rhodin, MPI, Sept 1999
  !
  ! for more details see file AUTHORS
  !

  USE mo_grid,          ONLY: plev, plond, pcnst, plon, i1
  USE mo_global_op,     ONLY: sum_zonal_sl

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  REAL, INTENT (IN) :: cwava, w
  INTEGER, INTENT (IN) :: jlat ! local latitude index N->S

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: etamid(plev), pdel(plond,plev), q1(plond,plev,pcnst), &
       &      q2(plond,plev,pcnst)
  INTEGER :: kftype(pcnst)

  !  Array arguments with intent(InOut):
  REAL, INTENT (INOUT) :: hwn(pcnst)

  ! cwava   l/(g*plon)
  ! w       Gaussian weight.
  ! q1      Untransported q-field.
  ! q2      Transported   q-field.
  ! pdel    array of pressure differences between layer interfaces
  !         (used for mass weighting ).
  ! hwn     Mass averaged constituent in units of kg/m**2.

  !  Local scalars: 
  REAL :: hwava
  INTEGER :: k, m

  !  Local arrays:
  REAL :: zhwava (1:plon,plev)

  !  Intrinsic functions 
  INTRINSIC ABS, SQRT

  !  Executable statements 

  ! accumulator

  DO m = 1, pcnst

     IF (kftype(m)==1) THEN
     hwava = sum_zonal_sl(&
       (q2(i1:plon+i1-1,:,m)*(SQRT(ABS(q1(i1:plon+i1-1,:,m)-&
        q2(i1:plon+i1-1,:,m))))**3)*pdel(i1:plon+i1-1,:), jlat)

!      more efficient but not identical to serial version:
!      zhwava = 0.
!      DO k = 1, plev
!        zhwava(:) = zhwava(:) + &
!          (q2(i1:plon+i1-1,k,m)*(SQRT(ABS(q1(i1:plon+i1-1,k,m)-&
!           q2(i1:plon+i1-1,k,m))))**3)*pdel(i1:plon+i1-1,k)
!      END DO
!      hwava = sum_zonal ( zhwava, jlat )

     ELSE IF (kftype(m)==2) THEN
       zhwava = 0.
       DO k = 1, plev
         zhwava(:,k)=(q2(i1:plon+i1-1,k,m)*etamid(k)*ABS(q1(i1:plon+i1-1,k,m)-&
                      q2(i1:plon+i1-1,k,m)))*pdel(i1:plon+i1-1,k)
       END DO
       hwava = sum_zonal_sl(zhwava, jlat )

!      more efficient but not identical to serial version:
!      zhwava = 0.
!      DO k = 1, plev
!        zhwava(:) = zhwava(:) + &
!          (q2(i1:plon+i1-1,k,m)*etamid(k)*ABS(q1(i1:plon+i1-1,k,m)-&
!           q2(i1:plon+i1-1,k,m)))*pdel(i1:plon+i1-1,k)
!      END DO
!      hwava = sum_zonal ( zhwava, jlat )

     END IF

     ! The 0.5 factor arises because gaussian weights sum to 2

     hwn(m) = hwn(m) + cwava*w*hwava*0.5
  END DO

  RETURN
END SUBROUTINE qmassd
