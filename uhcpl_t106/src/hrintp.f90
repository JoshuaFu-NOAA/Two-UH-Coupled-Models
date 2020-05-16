!+ interpolate 2-d field to departure point
!+ $Id: hrintp.f90,v 1.6 1999/08/27 17:00:23 m214030 Exp $

SUBROUTINE hrintp(pf,fb,fxl,fxr,x,y,dy,wdy,xdp,ydp,idp,jdp,limitd, &
&      fint,fyb,fyt,fdp)

  ! Description:
  !
  ! Interpolate 2-d field to departure point
  !
  ! Method:
  !
  ! Interpolate 2-d field to departure point using tensor product
  ! Hermite cubic interpolation.
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

  USE mo_grid, ONLY: platd, plond, plev, plon
  USE mo_slt,  ONLY: ppdy

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: pf
  LOGICAL, INTENT (IN) :: limitd    ! flag for shape-preservation

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: dy(platd), fb(plond,plev,pf,platd),  &
&      fxl(plond,plev,pf,platd), fxr(plond,plev,pf,platd),  &
&      wdy(4,2,platd), x(plond), xdp(plon,plev), y(platd), ydp(plon,plev)
  INTEGER, INTENT (IN) :: idp(plon,plev), jdp(plon,plev)

  !  Array arguments with intent(Out):
  REAL, INTENT (OUT) :: fint(plon,plev,ppdy,pf), fyb(plon,plev,pf),  &
&                       fyt(plon,plev,pf), fdp(plon,plev,pf)

  ! pf      Number of fields being interpolated.
  ! fb      Extended array of data to be interpolated.
  ! fxl     x-derivatives at the left  edge of each interval containing
  !         the departure point.
  ! fxr     x-derivatives at the right edge of each interval containing
  !         the departure point.
  ! x       Equally spaced x grid values in extended arrays.
  ! y       y-coordinate (latitude) values in the extended array.
  ! dy      Increment in the y-coordinate value for each interval in the
  !         extended array.
  ! wdy     Weights for Lagrange cubic derivative estimates on the
  !         unequally spaced y-grid.  If grid interval j (in extended
  !         array is surrounded by a 4 point stencil, then the
  !         derivative at the "bottom" of the interval uses the weights
  !         wdy(1,1,j),wdy(2,1,j), wdy(3,1,j), and wdy(4,1,j).  The
  !         derivative at the "top" of the interval uses wdy(1,2,j),
  !         wdy(2,2,j), wdy(3,2,j) and wdy(4,2,j).
  ! xdp     xdp(i,k) is the x-coordinate of the departure point that
  !         corresponds to global grid point (i,k) in the latitude slice
  !         being forecasted.
  ! ydp     ydp(i,k) is the y-coordinate of the departure point that
  !         corresponds to global grid point (i,k) in the latitude slice
  !         being forecasted.
  ! idp     idp(i,k) is the index of the x-interval that contains the
  !         departure point corresponding to global grid point (i,k) in
  !         the latitude slice being forecasted.
  !         Note that
  !         x(idp(i,k)) .le. xdp(i,k) .lt. x(idp(i,k)+1) .
  ! jdp     jdp(i,k) is the index of the y-interval that contains the
  !         departure point corresponding to global grid point (i,k) in
  !         the latitude slice being forecasted.
  !         Suppose yb contains the y-coordinates of the extended array
  !         and ydp(i,k) is the y-coordinate of the departure point
  !         corresponding to grid point (i,k).  Then,
  !         yb(jdp(i,k)) .le. ydp(i,k) .lt. yb(jdp(i,k)+1) .
  !         limitd  Logical flag to specify whether or not the y-derivatives
  !         will be limited.
  ! fint    WORK ARRAY, results not used on return
  ! fyb     WORK ARRAY, results not used on return
  ! fyt     WORK ARRAY, results not used on return
  ! fdp     Value of field at the horizontal departure points.

  !  External aubroutines 
  EXTERNAL cubydr, herxin, heryin, limdy


  !  Executable atatements 

  ! Hermite cubic interpolation to the x-coordinate of each
  ! departure point at each y-coordinate required to compute the
  ! y-derivatives.

  CALL herxin(pf,fb,fxl,fxr,x,xdp,idp,jdp,fint)

  ! Compute y-derivatives.

  CALL cubydr(pf,fint,wdy,jdp,fyb,fyt)
  IF (limitd) THEN
    CALL limdy(pf,fint,dy,jdp,fyb,fyt)
  END IF

  ! Hermite cubic interpolation in the y-coordinate.

  CALL heryin(pf,fint,fyb,fyt,y,dy,ydp,jdp,fdp)

  RETURN
END SUBROUTINE hrintp
