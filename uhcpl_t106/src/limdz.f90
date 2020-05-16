!+ apply SCMO limiter to vertical derivative estimates on a vertical slice
!+ $Id: limdz.f90,v 1.5 1999/08/27 17:00:28 m214030 Exp $

SUBROUTINE limdz(f,dsig,fst,fsb)

  ! Description:
  !
  ! Apply SCMO limiter to vertical derivative estimates on a vertical
  ! slice.
  !
  ! Method:
  !
  ! Authors:
  !
  ! Original version:  J. Olson
  ! Standardized:      J. Rosinski, June 1992
  ! Reviewed:          D. Williamson, P. Rasch, August 1992
  ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_grid, ONLY: plon, plev, pcnst, plevm1

  IMPLICIT NONE

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: f(plon,plev,pcnst), dsig(plev)

  !  Array arguments with intent(InOut):
  REAL, INTENT (INOUT) :: fst(plon,plev,pcnst), fsb(plon,plev,pcnst)

  ! f       Field values used to compute the discrete differences for
  !         each interval in the vertical grid.
  ! dsig    Increment in the sigma-coordinate value for each interval.
  ! fst     Limited derivative at the top of each interval.
  ! fsb     Limited derivative at the bottom of each interval.

  !  Local scalars: 
  REAL :: rdsig                ! 1./dsig
  INTEGER :: i                 ! longitude   index
  INTEGER :: k                 ! vertical    index
  INTEGER :: m                 ! constituent index

  !  Local arrays: 
  REAL :: deli(plon,plevm1)    ! simple linear derivative

  !  External subroutines 
  EXTERNAL scm0


  !  Executable statements 

  ! Loop over fields.

  DO m = 1, pcnst
    DO k = 1, plevm1
      rdsig = 1.0/dsig(k)
      DO i = 1, plon
        deli(i,k) = (f(i,k+1,m)-f(i,k,m))*rdsig
      END DO
    END DO
    CALL scm0(plon*plevm1,deli,fst(1,1,m),fsb(1,1,m))
  END DO

  RETURN
END SUBROUTINE limdz
