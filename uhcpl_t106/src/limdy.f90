!+ Limit the y-derivative estimates
!+ $Id: limdy.f90,v 1.5 1999/08/27 17:00:27 m214030 Exp $

SUBROUTINE limdy(pf,fint,dy,jdp,fyb,fyt)

  ! Description:
  !
  ! Limit the y-derivative estimates.
  !
  ! Method:
  !
  ! Limit the y-derivative estimates so they satisy the SCM0 for the
  ! x-interpolated data corresponding to the departure points of a single
  ! latitude slice in the global grid, that is, they are monotonic, but
  ! spline has only C0 continuity
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

  USE mo_grid, ONLY: plon, plev, platd
  USE mo_slt,  ONLY: ppdy

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: pf

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: fint(plon,plev,ppdy,pf), dy(platd)
  INTEGER, INTENT (IN) :: jdp(plon,plev)

  !  Array arguments with intent(InOut):
  REAL, INTENT (INOUT) :: fyb(plon,plev,pf), fyt(plon,plev,pf)

  ! pf      Number of fields being interpolated.
  ! fint    (fint(i,k,j,m),j=1,ppdy) contains the x interpolants at each
  !         latitude needed for the y derivative estimates at the
  !         endpoints of the interval that contains the departure point
  !         for grid point (i,k).  The last index of fint allows for
  !         interpolation of multiple fields.  fint is generated by a
  !         call to herxin.
  ! dy      Increment in the y-coordinate value for each interval in the
  !         extended array.
  ! jdp     jdp(i,k) is the index of the y-interval that contains the
  !         departure point corresponding to global grid point (i,k) in
  !         the latitude slice being forecasted.
  !         Suppose yb contains the y-coordinates of the extended array
  !         and ydp(i,k) is the y-coordinate of the departure point
  !         corresponding to grid point (i,k).  Then,
  !         yb(jdp(i,k)) .le. ydp(i,k) .lt. yb(jdp(i,k)+1) .
  ! fyb     fyb(i,k,.) is the limited derivative at the bot of the y
  !         interval that contains the departure point of global grid
  !         point (i,k).
  ! fyt     fyt(i,k,.) is the limited derivative at the top of the y
  !         interval that contains the departure point of global grid
  !         point (i,k).

  !  Local scalars: 
  INTEGER :: i, k, m          ! indices
  INTEGER :: jb               ! index for bottom of interval
  INTEGER :: jt               ! index for top    of interval

  !  Local arrays: 
  REAL :: rdy(plon,plev)      ! 1./dy
  REAL :: deli(plon,plev)     ! simple linear derivative

  !  External subroutines 
  EXTERNAL scm0


  !  Executable statements 

  jb = ppdy/2
  jt = jb + 1

  DO k = 1, plev
    DO i = 1, plon
      rdy(i,k) = 1./dy(jdp(i,k))
    END DO
  END DO

  ! Loop over fields.

  DO m = 1, pf
    DO k = 1, plev
      DO i = 1, plon
        deli(i,k) = (fint(i,k,jt,m)-fint(i,k,jb,m))*rdy(i,k)
      END DO
    END DO

    ! Limiter

    CALL scm0(plon*plev,deli,fyb(1,1,m),fyt(1,1,m))
  END DO

  RETURN
END SUBROUTINE limdy